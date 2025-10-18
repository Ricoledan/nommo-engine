"""
Temperature control methods for molecular dynamics simulations.

Implements various thermostat algorithms to maintain constant temperature
(NVT ensemble) simulations.
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Optional
from numba import njit

from nommo.physics.constants import KB_KJMOL_K, ENERGY_CONVERSION
from nommo.utils.logging import get_logger

logger = get_logger("thermostat")


@njit
def calculate_temperature(
    velocities: np.ndarray,
    masses: np.ndarray
) -> float:
    """
    Calculate instantaneous temperature from velocities.
    
    T = 2*E_kinetic / (N_dof * k_B)
    
    Args:
        velocities: Particle velocities (N x 3)
        masses: Particle masses (N,)
        
    Returns:
        Temperature in Kelvin
    """
    n_particles = velocities.shape[0]
    n_dof = 3 * n_particles
    
    kinetic_energy = 0.0
    for i in range(n_particles):
        v_squared = np.sum(velocities[i] ** 2)
        kinetic_energy += 0.5 * masses[i] * v_squared * ENERGY_CONVERSION
    
    if n_dof > 0:
        temperature = 2 * kinetic_energy / (n_dof * KB_KJMOL_K)
    else:
        temperature = 0.0
        
    return temperature


@njit
def rescale_velocities(
    velocities: np.ndarray,
    current_temp: float,
    target_temp: float
) -> None:
    """
    Simple velocity rescaling.
    
    Args:
        velocities: Particle velocities (modified in place)
        current_temp: Current temperature
        target_temp: Target temperature
    """
    if current_temp > 1e-10:
        scale_factor = np.sqrt(target_temp / current_temp)
        velocities *= scale_factor


class Thermostat(ABC):
    """Base class for thermostat implementations."""
    
    @abstractmethod
    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        dt: float
    ) -> float:
        """
        Apply thermostat to velocities.
        
        Args:
            velocities: Particle velocities (modified in place)
            masses: Particle masses
            dt: Timestep
            
        Returns:
            Current temperature after thermostat
        """
        pass
    
    @property
    @abstractmethod
    def name(self) -> str:
        """Thermostat name."""
        pass


class BerendsenThermostat(Thermostat):
    """
    Berendsen thermostat (weak coupling).
    
    Rescales velocities to gradually approach target temperature.
    Good for equilibration but doesn't generate correct canonical ensemble.
    
    Reference: Berendsen et al. (1984) J. Chem. Phys. 81, 3684
    """
    
    def __init__(
        self,
        target_temperature: float,
        coupling_time: float = 0.1
    ):
        """
        Initialize Berendsen thermostat.
        
        Args:
            target_temperature: Target temperature (K)
            coupling_time: Coupling time constant (ps)
        """
        self.T_target = target_temperature
        self.tau = coupling_time
        
        logger.info(
            f"Initialized Berendsen thermostat: T={target_temperature}K, tau={coupling_time}ps"
        )
    
    @property
    def name(self) -> str:
        return "Berendsen"
    
    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        dt: float
    ) -> float:
        """
        Apply Berendsen thermostat.
        
        λ = sqrt(1 + (dt/τ)(T_target/T_current - 1))
        v_new = λ * v_old
        """
        T_current = calculate_temperature(velocities, masses)
        
        if T_current < 1e-6:
            T_current = self.T_target
            self._initialize_velocities(velocities, masses)
            return T_current
        
        lambda_factor = np.sqrt(
            1 + (dt / self.tau) * (self.T_target / T_current - 1)
        )
        
        lambda_factor = np.clip(lambda_factor, 0.8, 1.25)
        
        velocities *= lambda_factor
        
        return calculate_temperature(velocities, masses)
    
    def _initialize_velocities(
        self,
        velocities: np.ndarray,
        masses: np.ndarray
    ):
        """Initialize velocities with Maxwell-Boltzmann distribution."""
        n_particles = velocities.shape[0]
        
        for i in range(n_particles):
            sigma = np.sqrt(KB_KJMOL_K * self.T_target / (masses[i] * ENERGY_CONVERSION))
            velocities[i] = np.random.normal(0, sigma, 3)
        
        velocities -= np.mean(velocities, axis=0)


class VelocityRescaleThermostat(Thermostat):
    """
    Simple velocity rescaling thermostat.
    
    Instantly rescales velocities to exact target temperature.
    Useful for initial equilibration but disrupts dynamics.
    """
    
    def __init__(
        self,
        target_temperature: float,
        interval: int = 10
    ):
        """
        Initialize velocity rescale thermostat.
        
        Args:
            target_temperature: Target temperature (K)
            interval: Steps between rescaling
        """
        self.T_target = target_temperature
        self.interval = interval
        self.step_count = 0
        
    @property
    def name(self) -> str:
        return "VelocityRescale"
    
    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        dt: float
    ) -> float:
        """Apply velocity rescaling."""
        self.step_count += 1
        
        T_current = calculate_temperature(velocities, masses)
        
        if self.step_count % self.interval == 0:
            rescale_velocities(velocities, T_current, self.T_target)
            T_current = self.T_target
            
        return T_current


class AndersenThermostat(Thermostat):
    """
    Andersen thermostat (stochastic collisions).
    
    Randomly reassigns velocities from Maxwell-Boltzmann distribution.
    Generates correct canonical ensemble but disrupts dynamics.
    
    Reference: Andersen (1980) J. Chem. Phys. 72, 2384
    """
    
    def __init__(
        self,
        target_temperature: float,
        collision_frequency: float = 0.1
    ):
        """
        Initialize Andersen thermostat.
        
        Args:
            target_temperature: Target temperature (K)
            collision_frequency: Collision frequency (1/ps)
        """
        self.T_target = target_temperature
        self.nu = collision_frequency
        
        logger.info(
            f"Initialized Andersen thermostat: T={target_temperature}K, nu={collision_frequency}/ps"
        )
    
    @property
    def name(self) -> str:
        return "Andersen"
    
    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        dt: float
    ) -> float:
        """
        Apply Andersen thermostat.
        
        Each particle has probability P = 1 - exp(-nu*dt) of collision.
        """
        n_particles = velocities.shape[0]
        collision_prob = 1 - np.exp(-self.nu * dt)
        
        for i in range(n_particles):
            if np.random.random() < collision_prob:
                sigma = np.sqrt(KB_KJMOL_K * self.T_target / (masses[i] * ENERGY_CONVERSION))
                velocities[i] = np.random.normal(0, sigma, 3)
        
        return calculate_temperature(velocities, masses)


class NoseHooverThermostat(Thermostat):
    """
    Nosé-Hoover thermostat (extended system).
    
    Introduces additional degree of freedom for heat bath.
    Generates correct canonical ensemble with smooth dynamics.
    
    References:
    - Nosé (1984) Mol. Phys. 52, 255
    - Hoover (1985) Phys. Rev. A 31, 1695
    """
    
    def __init__(
        self,
        target_temperature: float,
        coupling_time: float = 0.1
    ):
        """
        Initialize Nosé-Hoover thermostat.
        
        Args:
            target_temperature: Target temperature (K)
            coupling_time: Coupling time constant (ps)
        """
        self.T_target = target_temperature
        self.tau = coupling_time
        self.Q = None
        self.xi = 0.0
        
        logger.info(
            f"Initialized Nosé-Hoover thermostat: T={target_temperature}K, tau={coupling_time}ps"
        )
    
    @property
    def name(self) -> str:
        return "NoseHoover"
    
    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        dt: float
    ) -> float:
        """
        Apply Nosé-Hoover thermostat.
        
        Extended system dynamics with friction coefficient xi.
        """
        n_particles = velocities.shape[0]
        n_dof = 3 * n_particles
        
        if self.Q is None:
            self.Q = n_dof * KB_KJMOL_K * self.T_target * self.tau ** 2
        
        T_current = calculate_temperature(velocities, masses)
        
        kinetic_energy = 0.5 * n_dof * KB_KJMOL_K * T_current
        
        self.xi += dt * (kinetic_energy - 0.5 * n_dof * KB_KJMOL_K * self.T_target) / self.Q
        
        velocities *= np.exp(-self.xi * dt)
        
        self.xi += dt * (kinetic_energy - 0.5 * n_dof * KB_KJMOL_K * self.T_target) / self.Q
        
        return calculate_temperature(velocities, masses)


class AdaptiveThermostat(Thermostat):
    """
    Adaptive thermostat that switches between aggressive and gentle modes.
    
    Uses aggressive rescaling during equilibration and switches to
    gentler method for production runs.
    """
    
    def __init__(
        self,
        target_temperature: float,
        equilibration_steps: int = 1000,
        tolerance: float = 10.0
    ):
        """
        Initialize adaptive thermostat.
        
        Args:
            target_temperature: Target temperature (K)
            equilibration_steps: Steps for initial equilibration
            tolerance: Temperature tolerance for switching (K)
        """
        self.T_target = target_temperature
        self.equilibration_steps = equilibration_steps
        self.tolerance = tolerance
        self.step_count = 0
        
        self.aggressive = VelocityRescaleThermostat(target_temperature, interval=1)
        self.gentle = BerendsenThermostat(target_temperature, coupling_time=0.5)
        
    @property
    def name(self) -> str:
        return "Adaptive"
    
    def apply(
        self,
        velocities: np.ndarray,
        masses: np.ndarray,
        dt: float
    ) -> float:
        """Apply adaptive thermostat."""
        self.step_count += 1
        T_current = calculate_temperature(velocities, masses)
        
        if (
            self.step_count < self.equilibration_steps or
            abs(T_current - self.T_target) > self.tolerance
        ):
            return self.aggressive.apply(velocities, masses, dt)
        else:
            return self.gentle.apply(velocities, masses, dt)


def create_thermostat(
    thermostat_type: str,
    target_temperature: float,
    **kwargs
) -> Optional[Thermostat]:
    """
    Factory function to create thermostat.
    
    Args:
        thermostat_type: Type of thermostat
        target_temperature: Target temperature
        **kwargs: Additional thermostat parameters
        
    Returns:
        Thermostat instance or None
    """
    thermostat_map = {
        "berendsen": BerendsenThermostat,
        "velocity_rescale": VelocityRescaleThermostat,
        "andersen": AndersenThermostat,
        "nose_hoover": NoseHooverThermostat,
        "adaptive": AdaptiveThermostat,
    }
    
    thermostat_class = thermostat_map.get(thermostat_type.lower())
    if thermostat_class:
        return thermostat_class(target_temperature, **kwargs)
    elif thermostat_type.lower() == "none":
        return None
    else:
        logger.warning(f"Unknown thermostat type: {thermostat_type}")
        return None