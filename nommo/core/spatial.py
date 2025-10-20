"""
Spatial indexing for efficient neighbor finding in molecular dynamics.

Implements cell lists and Verlet neighbor lists for O(N) scaling
instead of naive O(N²) pairwise distance calculations.
"""

from dataclasses import dataclass, field

import numpy as np
from numba import njit, prange

from nommo.utils.logging import get_logger, get_performance_logger

logger = get_logger("spatial")
perf_logger = get_performance_logger()


@dataclass
class NeighborList:
    """Verlet neighbor list for tracking particle pairs."""

    pairs: list[tuple[int, int]] = field(default_factory=list)
    skin_distance: float = 0.3
    cutoff_distance: float = 1.0
    max_displacement: np.ndarray = field(default_factory=lambda: np.zeros(1))
    reference_positions: np.ndarray | None = None

    @property
    def rebuild_threshold(self) -> float:
        """Maximum displacement before rebuild needed."""
        return self.skin_distance / 2.0

    def needs_update(self, positions: np.ndarray) -> bool:
        """Check if neighbor list needs rebuilding."""
        if self.reference_positions is None:
            return True

        displacements = positions - self.reference_positions
        max_disp = np.max(np.sqrt(np.sum(displacements**2, axis=1)))
        return max_disp > self.rebuild_threshold

    def update_reference(self, positions: np.ndarray):
        """Update reference positions after rebuild."""
        self.reference_positions = positions.copy()
        self.max_displacement[:] = 0.0


class CellList:
    """
    Cell list spatial indexing for O(N) neighbor finding.

    Divides simulation box into cells and tracks particles in each cell.
    Only checks neighboring cells for interactions.
    """

    def __init__(self, box_size: np.ndarray, cutoff: float, skin: float = 0.3):
        """
        Initialize cell list.

        Args:
            box_size: Simulation box dimensions (nm)
            cutoff: Interaction cutoff distance (nm)
            skin: Verlet list skin distance (nm)
        """
        self.box_size = np.asarray(box_size)
        self.cutoff = cutoff
        self.skin = skin
        self.cell_size = cutoff + skin

        self.n_cells = np.maximum(1, (self.box_size / self.cell_size).astype(int))
        self.cell_size = self.box_size / self.n_cells

        self.cells: dict[tuple[int, int, int], list[int]] = {}
        self.particle_cells: dict[int, tuple[int, int, int]] = {}

        self._init_cells()

        logger.debug(
            f"Initialized cell list: {self.n_cells} cells, "
            f"cell_size={self.cell_size}, cutoff={cutoff}"
        )

    def _init_cells(self):
        """Initialize empty cell structure."""
        self.cells.clear()
        for i in range(self.n_cells[0]):
            for j in range(self.n_cells[1]):
                for k in range(self.n_cells[2]):
                    self.cells[(i, j, k)] = []

    def _get_cell_index(self, position: np.ndarray) -> tuple[int, int, int]:
        """Get cell indices for a position."""
        indices = (position / self.cell_size).astype(int)
        indices = np.clip(indices, 0, self.n_cells - 1)
        return tuple(indices)

    def update(self, positions: np.ndarray):
        """
        Update cell list with new particle positions.

        Args:
            positions: Particle positions (N x 3)
        """
        self._init_cells()
        self.particle_cells.clear()

        for idx, pos in enumerate(positions):
            cell_idx = self._get_cell_index(pos)
            self.cells[cell_idx].append(idx)
            self.particle_cells[idx] = cell_idx

    def get_neighbors_for_cell(self, cell_idx: tuple[int, int, int]) -> list[tuple[int, int, int]]:
        """Get all neighboring cells (including self)."""
        neighbors = []
        i, j, k = cell_idx

        for di in [-1, 0, 1]:
            for dj in [-1, 0, 1]:
                for dk in [-1, 0, 1]:
                    ni = (i + di) % self.n_cells[0]
                    nj = (j + dj) % self.n_cells[1]
                    nk = (k + dk) % self.n_cells[2]
                    neighbors.append((ni, nj, nk))

        return neighbors

    def get_potential_pairs(self) -> list[tuple[int, int]]:
        """Get all potential interacting particle pairs."""
        pairs = []

        for cell_idx, particles in self.cells.items():
            if not particles:
                continue

            for i in range(len(particles)):
                for j in range(i + 1, len(particles)):
                    pairs.append((particles[i], particles[j]))

            neighbor_cells = self.get_neighbors_for_cell(cell_idx)
            for neighbor_idx in neighbor_cells:
                if neighbor_idx <= cell_idx:
                    continue

                neighbor_particles = self.cells.get(neighbor_idx, [])
                for p1 in particles:
                    for p2 in neighbor_particles:
                        pairs.append((p1, p2))

        return pairs


@njit(parallel=True)
def build_neighbor_list_numba(
    positions: np.ndarray, box_size: np.ndarray, cutoff: float, skin: float
) -> list[tuple[int, int]]:
    """
    Build neighbor list using Numba for performance.

    Args:
        positions: Particle positions (N x 3)
        box_size: Box dimensions
        cutoff: Interaction cutoff
        skin: Verlet skin distance

    Returns:
        List of particle index pairs
    """
    n_particles = positions.shape[0]
    cutoff_sq = (cutoff + skin) ** 2
    pairs = []

    for i in prange(n_particles - 1):
        for j in range(i + 1, n_particles):
            dr = positions[j] - positions[i]

            for k in range(3):
                if dr[k] > box_size[k] / 2:
                    dr[k] -= box_size[k]
                elif dr[k] < -box_size[k] / 2:
                    dr[k] += box_size[k]

            dist_sq = np.sum(dr * dr)

            if dist_sq < cutoff_sq:
                pairs.append((i, j))

    return pairs


@njit
def minimum_image_distance(
    r1: np.ndarray, r2: np.ndarray, box_size: np.ndarray
) -> tuple[np.ndarray, float]:
    """
    Calculate minimum image distance with periodic boundaries.

    Args:
        r1: Position 1
        r2: Position 2
        box_size: Box dimensions

    Returns:
        Distance vector and scalar distance
    """
    dr = r2 - r1

    for i in range(3):
        if dr[i] > box_size[i] / 2:
            dr[i] -= box_size[i]
        elif dr[i] < -box_size[i] / 2:
            dr[i] += box_size[i]

    dist = np.sqrt(np.sum(dr * dr))
    return dr, dist


class SpatialIndex:
    """
    High-level spatial indexing combining cell lists and Verlet lists.
    """

    def __init__(
        self, box_size: np.ndarray, cutoff: float, skin: float = 0.3, use_cell_list: bool = True
    ):
        """
        Initialize spatial index.

        Args:
            box_size: Simulation box dimensions
            cutoff: Interaction cutoff distance
            skin: Verlet list skin distance
            use_cell_list: Whether to use cell lists
        """
        self.box_size = np.asarray(box_size)
        self.cutoff = cutoff
        self.skin = skin
        self.use_cell_list = use_cell_list

        if use_cell_list:
            self.cell_list = CellList(box_size, cutoff, skin)
        else:
            self.cell_list = None

        self.neighbor_list = NeighborList(skin_distance=skin, cutoff_distance=cutoff)

        self._update_count = 0
        self._rebuild_count = 0

    @perf_logger.timer("neighbor_list_update")
    def update(self, positions: np.ndarray, force_rebuild: bool = False):
        """
        Update spatial index with new positions.

        Args:
            positions: Particle positions
            force_rebuild: Force neighbor list rebuild
        """
        self._update_count += 1

        if force_rebuild or self.neighbor_list.needs_update(positions):
            self._rebuild_neighbor_list(positions)
            self._rebuild_count += 1

    def _rebuild_neighbor_list(self, positions: np.ndarray):
        """Rebuild the neighbor list."""
        if self.use_cell_list:
            self.cell_list.update(positions)
            pairs = self.cell_list.get_potential_pairs()

            filtered_pairs = []
            for i, j in pairs:
                dr, dist = minimum_image_distance(positions[i], positions[j], self.box_size)
                if dist < self.cutoff + self.skin:
                    filtered_pairs.append((i, j))

            self.neighbor_list.pairs = filtered_pairs
        else:
            self.neighbor_list.pairs = build_neighbor_list_numba(
                positions, self.box_size, self.cutoff, self.skin
            )

        self.neighbor_list.update_reference(positions)

        logger.debug(
            f"Rebuilt neighbor list: {len(self.neighbor_list.pairs)} pairs, "
            f"update {self._update_count}, rebuild {self._rebuild_count}"
        )

    def get_pairs(self) -> list[tuple[int, int]]:
        """Get current neighbor pairs."""
        return self.neighbor_list.pairs

    def get_neighbors(self, particle_idx: int) -> list[int]:
        """Get neighbors of a specific particle."""
        neighbors = []
        for i, j in self.neighbor_list.pairs:
            if i == particle_idx:
                neighbors.append(j)
            elif j == particle_idx:
                neighbors.append(i)
        return neighbors

    @property
    def efficiency_ratio(self) -> float:
        """Ratio of updates to rebuilds (higher is better)."""
        if self._rebuild_count == 0:
            return float("inf")
        return self._update_count / self._rebuild_count
