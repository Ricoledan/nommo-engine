"""
Time integration algorithms for molecular dynamics.

Implements Velocity Verlet and other symplectic integrators that
conserve energy and are time-reversible.
"""

from abc import ABC, abstractmethod
from collections.abc import Callable

import numpy as np
from numba import njit

from nommo.core.types import BoundaryType
from nommo.physics.constants import MAX_FORCE
from nommo.utils.logging import get_logger, get_performance_logger

logger = get_logger("integrator")
perf_logger = get_performance_logger()


@njit
def velocity_verlet_step(
    positions: np.ndarray,
    velocities: np.ndarray,
    forces: np.ndarray,
    masses: np.ndarray,
    dt: float,
    new_forces: np.ndarray,
) -> None:
    """
    Velocity Verlet integration step (Numba optimized).

    Updates positions and velocities in place.

    Args:
        positions: Particle positions (N x 3)
        velocities: Particle velocities (N x 3)
        forces: Current forces (N x 3)
        masses: Particle masses (N,)
        dt: Timestep
        new_forces: New forces after position update (N x 3)
    """
    dt_sq = dt * dt

    for i in range(positions.shape[0]):
        acceleration = forces[i] / masses[i]

        positions[i] += velocities[i] * dt + 0.5 * acceleration * dt_sq

        new_acceleration = new_forces[i] / masses[i]
        velocities[i] += 0.5 * (acceleration + new_acceleration) * dt


@njit
def apply_periodic_boundaries(positions: np.ndarray, box_size: np.ndarray) -> None:
    """Apply periodic boundary conditions."""
    for i in range(positions.shape[0]):
        for j in range(3):
            if positions[i, j] < 0:
                positions[i, j] += box_size[j]
            elif positions[i, j] >= box_size[j]:
                positions[i, j] -= box_size[j]


@njit
def apply_reflective_boundaries(
    positions: np.ndarray, velocities: np.ndarray, box_size: np.ndarray
) -> None:
    """Apply reflective boundary conditions."""
    for i in range(positions.shape[0]):
        for j in range(3):
            if positions[i, j] < 0:
                positions[i, j] = -positions[i, j]
                velocities[i, j] = -velocities[i, j]
            elif positions[i, j] > box_size[j]:
                positions[i, j] = 2 * box_size[j] - positions[i, j]
                velocities[i, j] = -velocities[i, j]


class Integrator(ABC):
    """Base class for time integration algorithms."""

    @abstractmethod
    def step(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        forces: np.ndarray,
        masses: np.ndarray,
        force_calculator: Callable,
        dt: float,
    ) -> float:
        """
        Perform one integration step.

        Args:
            positions: Particle positions
            velocities: Particle velocities
            forces: Current forces
            masses: Particle masses
            force_calculator: Function to calculate new forces
            dt: Timestep

        Returns:
            Potential energy after step
        """
        pass

    @abstractmethod
    def get_kinetic_energy(self, velocities: np.ndarray, masses: np.ndarray) -> float:
        """Calculate total kinetic energy."""
        pass


class VelocityVerletIntegrator(Integrator):
    """
    Velocity Verlet integrator.

    Second-order accurate, symplectic, time-reversible.
    Best for constant energy (NVE) simulations.
    """

    def __init__(
        self,
        timestep: float,
        box_size: np.ndarray,
        boundary_type: BoundaryType = BoundaryType.PERIODIC,
    ):
        """
        Initialize Velocity Verlet integrator.

        Args:
            timestep: Integration timestep (ps)
            box_size: Simulation box dimensions
            boundary_type: Boundary condition type
        """
        self.dt = timestep
        self.dt_sq = timestep**2
        self.box_size = box_size
        self.boundary_type = boundary_type

        logger.info(f"Initialized Velocity Verlet: dt={timestep} ps")

    @perf_logger.timer("integration_step")
    def step(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        forces: np.ndarray,
        masses: np.ndarray,
        force_calculator: Callable,
        dt: float | None = None,
    ) -> float:
        """
        Perform one Velocity Verlet step.

        Algorithm:
        1. x(t+dt) = x(t) + v(t)dt + 0.5*a(t)dt²
        2. Calculate forces at new positions
        3. v(t+dt) = v(t) + 0.5*[a(t) + a(t+dt)]dt
        """
        if dt is None:
            dt = self.dt

        accelerations = forces / masses[:, np.newaxis]

        positions += velocities * dt + 0.5 * accelerations * dt * dt

        if self.boundary_type == BoundaryType.PERIODIC:
            apply_periodic_boundaries(positions, self.box_size)
        elif self.boundary_type == BoundaryType.REFLECTIVE:
            apply_reflective_boundaries(positions, velocities, self.box_size)

        new_forces, potential_energy = force_calculator(positions)

        new_forces = np.clip(new_forces, -MAX_FORCE, MAX_FORCE)

        new_accelerations = new_forces / masses[:, np.newaxis]
        velocities += 0.5 * (accelerations + new_accelerations) * dt

        forces[:] = new_forces

        return potential_energy

    def get_kinetic_energy(self, velocities: np.ndarray, masses: np.ndarray) -> float:
        """Calculate total kinetic energy."""
        v_squared = np.sum(velocities**2, axis=1)
        return 0.5 * np.sum(masses * v_squared) * 0.01


class AdaptiveTimestepIntegrator(VelocityVerletIntegrator):
    """
    Velocity Verlet with adaptive timestep.

    Adjusts timestep based on maximum force to maintain stability.
    """

    def __init__(
        self,
        base_timestep: float,
        box_size: np.ndarray,
        boundary_type: BoundaryType = BoundaryType.PERIODIC,
        min_timestep: float = 1e-5,
        max_timestep: float = 0.01,
        force_tolerance: float = 1000.0,
    ):
        """
        Initialize adaptive timestep integrator.

        Args:
            base_timestep: Base integration timestep
            box_size: Simulation box dimensions
            boundary_type: Boundary conditions
            min_timestep: Minimum allowed timestep
            max_timestep: Maximum allowed timestep
            force_tolerance: Force threshold for timestep adjustment
        """
        super().__init__(base_timestep, box_size, boundary_type)
        self.base_timestep = base_timestep
        self.min_timestep = min_timestep
        self.max_timestep = max_timestep
        self.force_tolerance = force_tolerance
        self.current_timestep = base_timestep

        logger.info(
            f"Initialized adaptive integrator: dt={base_timestep} "
            f"({min_timestep}-{max_timestep})"
        )

    def calculate_timestep(self, forces: np.ndarray, masses: np.ndarray) -> float:
        """
        Calculate appropriate timestep based on forces.

        Args:
            forces: Current forces
            masses: Particle masses

        Returns:
            Adjusted timestep
        """
        max_force = np.max(np.linalg.norm(forces, axis=1))
        min_mass = np.min(masses)

        if max_force > 0:
            characteristic_time = np.sqrt(min_mass / max_force)
            suggested_dt = 0.01 * characteristic_time
        else:
            suggested_dt = self.base_timestep

        dt = np.clip(suggested_dt, self.min_timestep, self.max_timestep)

        if abs(dt - self.current_timestep) / self.current_timestep > 0.1:
            logger.debug(f"Adjusted timestep: {self.current_timestep:.6f} -> {dt:.6f}")

        self.current_timestep = dt
        return dt

    def step(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        forces: np.ndarray,
        masses: np.ndarray,
        force_calculator: Callable,
        dt: float | None = None,
    ) -> float:
        """Perform integration with adaptive timestep."""
        if dt is None:
            dt = self.calculate_timestep(forces, masses)

        return super().step(positions, velocities, forces, masses, force_calculator, dt)


class LeapfrogIntegrator(Integrator):
    """
    Leapfrog integrator.

    Similar accuracy to Velocity Verlet but with different
    phase space trajectory.
    """

    def __init__(
        self,
        timestep: float,
        box_size: np.ndarray,
        boundary_type: BoundaryType = BoundaryType.PERIODIC,
    ):
        """Initialize Leapfrog integrator."""
        self.dt = timestep
        self.box_size = box_size
        self.boundary_type = boundary_type
        self.first_step = True

    def step(
        self,
        positions: np.ndarray,
        velocities: np.ndarray,
        forces: np.ndarray,
        masses: np.ndarray,
        force_calculator: Callable,
        dt: float | None = None,
    ) -> float:
        """
        Perform one leapfrog step.

        v(t+dt/2) = v(t-dt/2) + a(t)*dt
        x(t+dt) = x(t) + v(t+dt/2)*dt
        """
        if dt is None:
            dt = self.dt

        if self.first_step:
            velocities -= 0.5 * dt * forces / masses[:, np.newaxis]
            self.first_step = False

        velocities += dt * forces / masses[:, np.newaxis]

        positions += velocities * dt

        if self.boundary_type == BoundaryType.PERIODIC:
            apply_periodic_boundaries(positions, self.box_size)
        elif self.boundary_type == BoundaryType.REFLECTIVE:
            apply_reflective_boundaries(positions, velocities, self.box_size)

        new_forces, potential_energy = force_calculator(positions)
        forces[:] = new_forces

        return potential_energy

    def get_kinetic_energy(self, velocities: np.ndarray, masses: np.ndarray) -> float:
        """Calculate kinetic energy (accounting for staggered velocities)."""
        v_squared = np.sum(velocities**2, axis=1)
        return 0.5 * np.sum(masses * v_squared) * 0.01
