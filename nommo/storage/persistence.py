"""
HDF5-based persistence for universe states and simulation data.

Provides efficient storage with compression and versioning.
"""

import json
from datetime import datetime
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from nommo.core.particle import Particle
from nommo.core.types import UniverseParameters
from nommo.core.universe import Universe, UniverseMetrics
from nommo.utils.logging import get_logger

logger = get_logger("storage")


class HDF5Storage:
    """
    HDF5-based storage for universe states.

    Structure:
    /universe_id/
        /metadata/
            parameters (JSON)
            creation_time
            last_modified
        /snapshots/
            /tick_XXXXX/
                particles/
                    positions (N x 3)
                    velocities (N x 3)
                    forces (N x 3)
                    properties (structured array)
                metrics (JSON)
                timestamp
        /trajectory/
            positions (T x N x 3)
            energies (T x 3)
            temperatures (T,)
    """

    def __init__(self, storage_dir: Path):
        """
        Initialize HDF5 storage.

        Args:
            storage_dir: Directory for storage files
        """
        self.storage_dir = Path(storage_dir)
        self.storage_dir.mkdir(parents=True, exist_ok=True)

    def save_universe(
        self, universe: Universe, filename: str | None = None, compress: bool = True
    ) -> Path:
        """
        Save universe state to HDF5 file.

        Args:
            universe: Universe to save
            filename: Optional filename (defaults to universe ID)
            compress: Whether to use compression

        Returns:
            Path to saved file
        """
        if filename is None:
            filename = f"{universe.params.name}_{universe.id[:8]}.h5"

        filepath = self.storage_dir / filename

        compression = "gzip" if compress else None
        compression_opts = 4 if compress else None

        with h5py.File(filepath, "w") as f:
            universe_group = f.create_group(universe.id)

            metadata = universe_group.create_group("metadata")
            metadata.attrs["parameters"] = json.dumps(universe.params.model_dump())
            metadata.attrs["creation_time"] = datetime.now().isoformat()
            metadata.attrs["tick"] = universe.tick
            metadata.attrs["time"] = universe.time

            snapshot = universe_group.create_group(f"snapshots/tick_{universe.tick:08d}")

            if universe.particles:
                positions = np.array([p.position for p in universe.particles])
                velocities = np.array([p.velocity for p in universe.particles])
                forces = np.array([p.forces for p in universe.particles])

                snapshot.create_dataset(
                    "particles/positions",
                    data=positions,
                    compression=compression,
                    compression_opts=compression_opts,
                )
                snapshot.create_dataset(
                    "particles/velocities",
                    data=velocities,
                    compression=compression,
                    compression_opts=compression_opts,
                )
                snapshot.create_dataset(
                    "particles/forces",
                    data=forces,
                    compression=compression,
                    compression_opts=compression_opts,
                )

                properties_dtype = np.dtype(
                    [
                        ("id", "S36"),
                        ("mass", "f8"),
                        ("charge", "f8"),
                        ("radius", "f8"),
                        ("type_index", "i4"),
                        ("age", "i4"),
                        ("generation", "i4"),
                    ]
                )

                properties = np.zeros(len(universe.particles), dtype=properties_dtype)
                for i, p in enumerate(universe.particles):
                    properties[i] = (
                        p.id.encode(),
                        p.mass,
                        p.charge,
                        p.radius,
                        p.type_index,
                        p.age,
                        p.generation,
                    )

                snapshot.create_dataset(
                    "particles/properties",
                    data=properties,
                    compression=compression,
                    compression_opts=compression_opts,
                )

                bond_data = []
                for p in universe.particles:
                    for bond_id in p.bonds:
                        bond_data.append((p.id, bond_id))

                if bond_data:
                    bond_dtype = np.dtype([("particle1", "S36"), ("particle2", "S36")])
                    bonds = np.array(bond_data, dtype=bond_dtype)
                    snapshot.create_dataset(
                        "bonds",
                        data=bonds,
                        compression=compression,
                        compression_opts=compression_opts,
                    )

            snapshot.attrs["metrics"] = json.dumps(universe.metrics.to_dict())
            snapshot.attrs["timestamp"] = datetime.now().isoformat()

        logger.info(f"Saved universe to {filepath}")
        return filepath

    def load_universe(self, filepath: Path, universe_id: str | None = None) -> Universe:
        """
        Load universe from HDF5 file.

        Args:
            filepath: Path to HDF5 file
            universe_id: Optional universe ID to load

        Returns:
            Loaded universe
        """
        with h5py.File(filepath, "r") as f:
            if universe_id is None:
                universe_id = list(f.keys())[0]

            universe_group = f[universe_id]

            params_json = universe_group["metadata"].attrs["parameters"]
            params = UniverseParameters(**json.loads(params_json))

            universe = Universe(params)
            universe.id = universe_id
            universe.tick = universe_group["metadata"].attrs["tick"]
            universe.time = universe_group["metadata"].attrs["time"]

            latest_snapshot = sorted(universe_group["snapshots"].keys())[-1]
            snapshot = universe_group[f"snapshots/{latest_snapshot}"]

            if "particles" in snapshot:
                positions = snapshot["particles/positions"][:]
                velocities = snapshot["particles/velocities"][:]
                forces = snapshot["particles/forces"][:]
                properties = snapshot["particles/properties"][:]

                universe.particles.clear()
                universe.particle_dict.clear()

                for i, prop in enumerate(properties):
                    particle = Particle(
                        id=prop["id"].decode(),
                        mass=prop["mass"],
                        charge=prop["charge"],
                        radius=prop["radius"],
                        position=positions[i],
                        velocity=velocities[i],
                        forces=forces[i],
                        type_index=prop["type_index"],
                        age=prop["age"],
                        generation=prop["generation"],
                    )
                    universe.particles.append(particle)
                    universe.particle_dict[particle.id] = particle

                if "bonds" in snapshot:
                    bonds = snapshot["bonds"][:]
                    for bond in bonds:
                        p1_id = bond["particle1"].decode()
                        p2_id = bond["particle2"].decode()
                        if p1_id in universe.particle_dict:
                            universe.particle_dict[p1_id].bonds.add(p2_id)

            if "metrics" in snapshot.attrs:
                metrics_dict = json.loads(snapshot.attrs["metrics"])
                universe.metrics = UniverseMetrics(**metrics_dict)

        logger.info(f"Loaded universe from {filepath}")
        return universe

    def save_trajectory(
        self, universe_id: str, history: list[UniverseMetrics], filename: str | None = None
    ) -> Path:
        """
        Save simulation trajectory data.

        Args:
            universe_id: Universe identifier
            history: List of metrics over time
            filename: Optional filename

        Returns:
            Path to saved file
        """
        if filename is None:
            filename = f"trajectory_{universe_id[:8]}.h5"

        filepath = self.storage_dir / filename

        with h5py.File(filepath, "w") as f:
            traj_group = f.create_group("trajectory")

            times = np.array([m.time for m in history])
            temperatures = np.array([m.current_temperature for m in history])
            kinetic = np.array([m.total_kinetic_energy for m in history])
            potential = np.array([m.total_potential_energy for m in history])
            total = np.array([m.total_energy for m in history])
            particles = np.array([m.total_particles for m in history])
            entropy = np.array([m.shannon_entropy for m in history])

            traj_group.create_dataset("times", data=times)
            traj_group.create_dataset("temperatures", data=temperatures)
            traj_group.create_dataset("kinetic_energy", data=kinetic)
            traj_group.create_dataset("potential_energy", data=potential)
            traj_group.create_dataset("total_energy", data=total)
            traj_group.create_dataset("particle_count", data=particles)
            traj_group.create_dataset("entropy", data=entropy)

            traj_group.attrs["universe_id"] = universe_id
            traj_group.attrs["n_frames"] = len(history)
            traj_group.attrs["total_time"] = times[-1] if len(times) > 0 else 0

        logger.info(f"Saved trajectory to {filepath}")
        return filepath

    def list_saves(self) -> list[dict[str, Any]]:
        """List all saved universes."""
        saves = []

        for filepath in self.storage_dir.glob("*.h5"):
            if filepath.stem.startswith("trajectory_"):
                continue

            try:
                with h5py.File(filepath, "r") as f:
                    for universe_id in f:
                        metadata = f[f"{universe_id}/metadata"]
                        saves.append(
                            {
                                "filename": filepath.name,
                                "universe_id": universe_id,
                                "tick": metadata.attrs.get("tick", 0),
                                "time": metadata.attrs.get("time", 0.0),
                                "creation_time": metadata.attrs.get("creation_time", ""),
                            }
                        )
            except Exception as e:
                logger.warning(f"Could not read {filepath}: {e}")

        return saves


class CheckpointManager:
    """
    Manages automatic checkpointing of simulations.
    """

    def __init__(self, storage: HDF5Storage, interval: int = 1000, keep_last: int = 5):
        """
        Initialize checkpoint manager.

        Args:
            storage: HDF5 storage instance
            interval: Ticks between checkpoints
            keep_last: Number of checkpoints to keep
        """
        self.storage = storage
        self.interval = interval
        self.keep_last = keep_last
        self.checkpoints: list[Path] = []

    def checkpoint(self, universe: Universe) -> Path | None:
        """
        Create a checkpoint if needed.

        Args:
            universe: Universe to checkpoint

        Returns:
            Path to checkpoint file if created
        """
        if universe.tick % self.interval != 0:
            return None

        filename = f"checkpoint_{universe.id[:8]}_{universe.tick:08d}.h5"
        filepath = self.storage.save_universe(universe, filename)

        self.checkpoints.append(filepath)

        if len(self.checkpoints) > self.keep_last:
            old_checkpoint = self.checkpoints.pop(0)
            try:
                old_checkpoint.unlink()
                logger.debug(f"Removed old checkpoint: {old_checkpoint}")
            except Exception as e:
                logger.warning(f"Could not remove checkpoint: {e}")

        return filepath

    def restore_latest(self, universe_id: str) -> Universe | None:
        """
        Restore from latest checkpoint.

        Args:
            universe_id: Universe ID to restore

        Returns:
            Restored universe or None
        """
        pattern = f"checkpoint_{universe_id[:8]}_*.h5"
        checkpoints = sorted(self.storage.storage_dir.glob(pattern))

        if not checkpoints:
            return None

        latest = checkpoints[-1]
        return self.storage.load_universe(latest, universe_id)
