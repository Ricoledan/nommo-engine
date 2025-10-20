"""
Force calculations for molecular dynamics using Lennard-Jones potential.

Implements efficient force calculations with Numba acceleration and
proper handling of periodic boundary conditions.
"""


import numba
import numpy as np
from numba import njit, prange

from nommo.core.spatial import minimum_image_distance
from nommo.physics.constants import MAX_FORCE, MIN_DISTANCE
from nommo.utils.logging import get_logger, get_performance_logger

logger = get_logger("forces")
perf_logger = get_performance_logger()


@njit
def lennard_jones_force(r: float, epsilon: float, sigma: float) -> tuple[float, float]:
    """
    Calculate Lennard-Jones force and potential.

    V(r) = 4ε[(σ/r)¹² - (σ/r)⁶]
    F(r) = 24ε/r [(2(σ/r)¹² - (σ/r)⁶)]

    Args:
        r: Distance between particles (nm)
        epsilon: Potential well depth (kJ/mol)
        sigma: Zero-crossing distance (nm)

    Returns:
        Force magnitude (kJ/mol/nm) and potential energy (kJ/mol)
    """
    if r < MIN_DISTANCE:
        r = MIN_DISTANCE

    sigma_r = sigma / r
    sigma_r6 = sigma_r**6
    sigma_r12 = sigma_r6**2

    force_mag = 24 * epsilon / r * (2 * sigma_r12 - sigma_r6)
    potential = 4 * epsilon * (sigma_r12 - sigma_r6)

    force_mag = min(force_mag, MAX_FORCE)

    return force_mag, potential


@njit
def lennard_jones_force_shifted(
    r: float, epsilon: float, sigma: float, cutoff: float
) -> tuple[float, float]:
    """
    Lennard-Jones with force shifting for smooth cutoff.

    Ensures force goes smoothly to zero at cutoff distance.
    """
    if r >= cutoff:
        return 0.0, 0.0

    force_mag, potential = lennard_jones_force(r, epsilon, sigma)

    force_cut, pot_cut = lennard_jones_force(cutoff, epsilon, sigma)

    alpha = -force_cut / cutoff
    force_mag = force_mag - force_cut - alpha * (r - cutoff)
    potential = potential - pot_cut - force_cut * (r - cutoff) - 0.5 * alpha * (r - cutoff) ** 2

    return force_mag, potential


@njit(parallel=True)
def calculate_forces_vectorized(
    positions: np.ndarray,
    types: np.ndarray,
    box_size: np.ndarray,
    epsilon_matrix: np.ndarray,
    sigma_matrix: np.ndarray,
    cutoff: float,
    neighbor_pairs: numba.typed.List,
) -> tuple[np.ndarray, float]:
    """
    Vectorized force calculation using neighbor list.

    Args:
        positions: Particle positions (N x 3)
        types: Particle type indices (N,)
        box_size: Simulation box dimensions
        epsilon_matrix: LJ epsilon parameters by type pair
        sigma_matrix: LJ sigma parameters by type pair
        cutoff: Force cutoff distance
        neighbor_pairs: List of interacting particle pairs

    Returns:
        Forces (N x 3) and total potential energy
    """
    n_particles = positions.shape[0]
    forces = np.zeros((n_particles, 3), dtype=np.float64)
    total_potential = 0.0

    for pair_idx in prange(len(neighbor_pairs)):
        i, j = neighbor_pairs[pair_idx]

        dr, r = minimum_image_distance(positions[i], positions[j], box_size)

        if r < cutoff and r > MIN_DISTANCE:
            type_i = types[i]
            type_j = types[j]
            epsilon = epsilon_matrix[type_i, type_j]
            sigma = sigma_matrix[type_i, type_j]

            force_mag, potential = lennard_jones_force_shifted(r, epsilon, sigma, cutoff)

            force_vec = force_mag * dr / r

            forces[i] -= force_vec
            forces[j] += force_vec

            total_potential += potential

    return forces, total_potential


class LennardJonesForce:
    """
    Lennard-Jones force field implementation.

    Manages force calculations with proper mixing rules and
    efficient neighbor list usage.
    """

    def __init__(
        self, epsilon: float = 1.0, sigma: float = 0.3, cutoff: float = 1.0, shift: bool = True
    ):
        """
        Initialize Lennard-Jones force field.

        Args:
            epsilon: Default potential well depth (kJ/mol)
            sigma: Default zero-crossing distance (nm)
            cutoff: Force cutoff distance (nm)
            shift: Whether to use force shifting
        """
        self.epsilon_default = epsilon
        self.sigma_default = sigma
        self.cutoff = cutoff
        self.cutoff_squared = cutoff**2
        self.shift = shift

        self.epsilon_matrix = None
        self.sigma_matrix = None
        self.n_types = 0

        logger.info(
            f"Initialized LJ force: epsilon={epsilon}, sigma={sigma}, "
            f"cutoff={cutoff}, shift={shift}"
        )

    def setup_parameters(
        self,
        n_types: int,
        epsilon_values: np.ndarray | None = None,
        sigma_values: np.ndarray | None = None,
    ):
        """
        Setup interaction parameters for particle types.

        Args:
            n_types: Number of particle types
            epsilon_values: Epsilon values by type (optional)
            sigma_values: Sigma values by type (optional)
        """
        self.n_types = n_types

        if epsilon_values is None:
            epsilon_values = np.full(n_types, self.epsilon_default)
        if sigma_values is None:
            sigma_values = np.full(n_types, self.sigma_default)

        self.epsilon_matrix = np.zeros((n_types, n_types))
        self.sigma_matrix = np.zeros((n_types, n_types))

        for i in range(n_types):
            for j in range(n_types):
                self.epsilon_matrix[i, j] = np.sqrt(epsilon_values[i] * epsilon_values[j])
                self.sigma_matrix[i, j] = (sigma_values[i] + sigma_values[j]) / 2

    @perf_logger.timer("force_calculation")
    def calculate(
        self,
        positions: np.ndarray,
        types: np.ndarray,
        box_size: np.ndarray,
        neighbor_pairs: list[tuple[int, int]],
    ) -> tuple[np.ndarray, float]:
        """
        Calculate forces on all particles.

        Args:
            positions: Particle positions
            types: Particle type indices
            box_size: Box dimensions
            neighbor_pairs: Neighbor list pairs

        Returns:
            Forces and total potential energy
        """
        if self.epsilon_matrix is None:
            max_type = np.max(types) + 1
            self.setup_parameters(max_type)

        typed_pairs = numba.typed.List()
        for pair in neighbor_pairs:
            typed_pairs.append(pair)

        forces, potential = calculate_forces_vectorized(
            positions,
            types,
            box_size,
            self.epsilon_matrix,
            self.sigma_matrix,
            self.cutoff,
            typed_pairs,
        )

        return forces, potential

    def calculate_single_pair(
        self, r: float, type1: int = 0, type2: int = 0
    ) -> tuple[float, float]:
        """
        Calculate force for a single particle pair.

        Args:
            r: Distance between particles
            type1: Type of particle 1
            type2: Type of particle 2

        Returns:
            Force magnitude and potential energy
        """
        if self.epsilon_matrix is None:
            epsilon = self.epsilon_default
            sigma = self.sigma_default
        else:
            epsilon = self.epsilon_matrix[type1, type2]
            sigma = self.sigma_matrix[type1, type2]

        if self.shift:
            return lennard_jones_force_shifted(r, epsilon, sigma, self.cutoff)
        else:
            if r >= self.cutoff:
                return 0.0, 0.0
            return lennard_jones_force(r, epsilon, sigma)

    def get_potential_curve(
        self, r_values: np.ndarray, type1: int = 0, type2: int = 0
    ) -> tuple[np.ndarray, np.ndarray]:
        """
        Get potential and force curves for plotting.

        Args:
            r_values: Distance values
            type1: Type of particle 1
            type2: Type of particle 2

        Returns:
            Force and potential arrays
        """
        forces = np.zeros_like(r_values)
        potentials = np.zeros_like(r_values)

        for i, r in enumerate(r_values):
            if r > MIN_DISTANCE:
                forces[i], potentials[i] = self.calculate_single_pair(r, type1, type2)

        return forces, potentials


class CoulombicForce:
    """
    Electrostatic force calculations.

    Implements Coulomb's law with proper treatment of
    periodic boundaries using Ewald summation (simplified).
    """

    def __init__(
        self, coulomb_constant: float = 138.935, cutoff: float = 1.0, dielectric: float = 1.0
    ):
        """
        Initialize Coulombic force.

        Args:
            coulomb_constant: Coulomb constant (kJ·nm/mol/e²)
            cutoff: Cutoff distance (nm)
            dielectric: Dielectric constant
        """
        self.coulomb_constant = coulomb_constant / dielectric
        self.cutoff = cutoff
        self.cutoff_squared = cutoff**2

    def calculate(
        self,
        positions: np.ndarray,
        charges: np.ndarray,
        box_size: np.ndarray,
        neighbor_pairs: list[tuple[int, int]],
    ) -> tuple[np.ndarray, float]:
        """
        Calculate electrostatic forces.

        Args:
            positions: Particle positions
            charges: Particle charges
            box_size: Box dimensions
            neighbor_pairs: Neighbor list

        Returns:
            Forces and electrostatic potential energy
        """
        n_particles = positions.shape[0]
        forces = np.zeros((n_particles, 3))
        total_potential = 0.0

        for i, j in neighbor_pairs:
            dr, r = minimum_image_distance(positions[i], positions[j], box_size)

            if r < self.cutoff and r > MIN_DISTANCE:
                charge_product = charges[i] * charges[j]

                # Coulomb force magnitude (positive for repulsion, negative for attraction)
                force_mag = self.coulomb_constant * charge_product / (r**2)
                potential = self.coulomb_constant * charge_product / r

                # Force vector points from i to j for repulsion, j to i for attraction
                # dr points from i to j, so force on i is along -dr for repulsion
                force_vec = force_mag * dr / r

                # Apply equal and opposite forces
                forces[i] += force_vec
                forces[j] -= force_vec

                total_potential += potential

        return forces, total_potential


class CompositeForce:
    """
    Combines multiple force fields.
    """

    def __init__(self):
        """Initialize composite force."""
        self.force_fields = []

    def add_force_field(self, force_field):
        """Add a force field component."""
        self.force_fields.append(force_field)

    def calculate(
        self,
        positions: np.ndarray,
        properties: dict,
        box_size: np.ndarray,
        neighbor_pairs: list[tuple[int, int]],
    ) -> tuple[np.ndarray, float]:
        """
        Calculate total forces from all components.

        Args:
            positions: Particle positions
            properties: Particle properties (types, charges, etc.)
            box_size: Box dimensions
            neighbor_pairs: Neighbor list

        Returns:
            Total forces and potential energy
        """
        n_particles = positions.shape[0]
        total_forces = np.zeros((n_particles, 3))
        total_potential = 0.0

        for ff in self.force_fields:
            if isinstance(ff, LennardJonesForce):
                forces, potential = ff.calculate(
                    positions,
                    properties.get("types", np.zeros(n_particles, dtype=int)),
                    box_size,
                    neighbor_pairs,
                )
            elif isinstance(ff, CoulombicForce):
                forces, potential = ff.calculate(
                    positions,
                    properties.get("charges", np.zeros(n_particles)),
                    box_size,
                    neighbor_pairs,
                )
            else:
                continue

            total_forces += forces
            total_potential += potential

        return total_forces, total_potential
