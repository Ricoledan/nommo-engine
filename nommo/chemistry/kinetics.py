"""
Arrhenius kinetics and collision theory implementation.

This module implements temperature-dependent reaction kinetics based on
the Arrhenius equation and collision theory principles.

References:
- Steinfeld et al. (1999) Chemical Kinetics and Dynamics
- Laidler (1987) Chemical Kinetics
- Atkins & de Paula (2017) Physical Chemistry
"""

from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from nommo.core.particle import Particle

from nommo.utils.logging import get_logger

logger = get_logger("kinetics")


class ArrheniusKinetics:
    """
    Temperature-dependent reaction kinetics using the Arrhenius equation.

    The Arrhenius equation describes how reaction rates depend on temperature:
    k(T) = A * exp(-E_a / RT)

    Where:
    - k(T): Rate constant at temperature T
    - A: Pre-exponential factor (frequency factor)
    - E_a: Activation energy (kJ/mol)
    - R: Gas constant (8.314 J/(mol·K) = 0.008314 kJ/(mol·K))
    - T: Temperature (K)

    This implementation also includes collision theory for determining
    whether individual particle collisions have sufficient energy to react.
    """

    # Physical constants
    R = 0.008314  # Gas constant in kJ/(mol·K)

    def __init__(
        self,
        activation_energy: float,
        pre_exponential: float = 1e13,
        steric_factor: float = 1.0
    ):
        """
        Initialize Arrhenius kinetics calculator.

        Args:
            activation_energy: Energy barrier for reaction (kJ/mol)
            pre_exponential: Attempt frequency factor (1/ps)
            steric_factor: Probability factor for correct orientation (0-1)
        """
        self.E_a = activation_energy
        self.A = pre_exponential
        self.steric_factor = steric_factor

        # Validate parameters
        if self.E_a < 0:
            raise ValueError("Activation energy must be non-negative")
        if self.A <= 0:
            raise ValueError("Pre-exponential factor must be positive")
        if not 0 <= self.steric_factor <= 1:
            raise ValueError("Steric factor must be between 0 and 1")

        logger.debug(f"Initialized Arrhenius kinetics: E_a={self.E_a:.2f} kJ/mol, "
                    f"A={self.A:.2e} 1/ps, steric={self.steric_factor:.3f}")

    def rate_constant(self, temperature: float) -> float:
        """
        Calculate the rate constant at given temperature.

        Args:
            temperature: Temperature in Kelvin

        Returns:
            Rate constant in 1/ps
        """
        if temperature <= 0:
            raise ValueError("Temperature must be positive")

        exponent = -self.E_a / (self.R * temperature)

        # Prevent numerical overflow for very large activation energies
        if exponent < -700:  # exp(-700) ≈ 0
            return 0.0

        return self.A * self.steric_factor * np.exp(exponent)

    def reaction_probability(self, temperature: float, timestep: float) -> float:
        """
        Calculate probability of reaction occurring in a timestep.

        Uses first-order kinetics approximation:
        P = 1 - exp(-k * dt)

        For small k*dt, this approximates to k*dt.

        Args:
            temperature: Temperature in Kelvin
            timestep: Time interval in picoseconds

        Returns:
            Probability of reaction (0-1)
        """
        k = self.rate_constant(temperature)
        rate_arg = k * timestep

        # For numerical stability, use linear approximation when rate_arg is small
        if rate_arg < 1e-6:
            return rate_arg
        else:
            return 1 - np.exp(-rate_arg)

    def has_sufficient_energy(
        self,
        particle1: "Particle",
        particle2: "Particle",
        temperature: float
    ) -> bool:
        """
        Determine if a collision has sufficient energy for reaction.

        Uses collision theory: compares the kinetic energy of relative motion
        to the activation energy, with Maxwell-Boltzmann thermal fluctuations.

        Args:
            particle1: First colliding particle
            particle2: Second colliding particle
            temperature: System temperature in Kelvin

        Returns:
            True if collision can overcome activation barrier
        """
        # Calculate relative velocity
        v_rel = particle1.velocity - particle2.velocity
        v_rel_magnitude = np.linalg.norm(v_rel)

        # Calculate reduced mass (effective mass for collision)
        m_reduced = (particle1.mass * particle2.mass) / (
            particle1.mass + particle2.mass
        )

        # Kinetic energy of relative motion in simulation units
        E_rel_sim = 0.5 * m_reduced * v_rel_magnitude**2

        # Convert to kJ/mol (conversion factor from simulation units)
        # Assuming mass in amu, velocity in nm/ps
        E_rel_kjmol = E_rel_sim * 0.01 * 6.022e23 / 1000

        # If collision energy exceeds activation energy, reaction can proceed
        if E_rel_kjmol >= self.E_a:
            return True

        # Otherwise, use Maxwell-Boltzmann probability for thermal activation
        energy_deficit = self.E_a - E_rel_kjmol

        # Boltzmann factor for overcoming energy deficit
        if energy_deficit > 100 * self.R * temperature:  # Avoid numerical underflow
            return False

        boltzmann_prob = np.exp(-energy_deficit / (self.R * temperature))
        return np.random.random() < boltzmann_prob

    def collision_frequency(
        self,
        particle1: "Particle",
        particle2: "Particle",
        temperature: float
    ) -> float:
        """
        Calculate collision frequency between two particles.

        Based on kinetic theory of gases:
        Z = σ * v_rel * n

        Where σ is collision cross-section and v_rel is relative velocity.

        Args:
            particle1: First particle
            particle2: Second particle
            temperature: System temperature in Kelvin

        Returns:
            Collision frequency in 1/ps
        """
        # Collision cross-section (geometric)
        sigma = np.pi * (particle1.radius + particle2.radius)**2

        # Relative velocity magnitude
        v_rel = np.linalg.norm(particle1.velocity - particle2.velocity)

        # For dilute gas approximation, frequency is proportional to v_rel
        # The exact prefactor depends on particle density, which we approximate as 1
        frequency = sigma * v_rel / (np.pi * 0.1)  # Rough normalization

        return frequency

    def effective_rate_constant(
        self,
        particle1: "Particle",
        particle2: "Particle",
        temperature: float
    ) -> float:
        """
        Calculate effective rate constant for specific particle pair.

        Combines Arrhenius kinetics with collision theory for the
        specific kinetic states of two particles.

        Args:
            particle1: First particle
            particle2: Second particle
            temperature: System temperature in Kelvin

        Returns:
            Effective rate constant in 1/ps
        """
        base_rate = self.rate_constant(temperature)

        # Collision frequency factor
        freq_factor = self.collision_frequency(particle1, particle2, temperature)

        # Energy availability factor
        if self.has_sufficient_energy(particle1, particle2, temperature):
            energy_factor = 1.0
        else:
            # Reduced rate for insufficient energy
            v_rel_mag = np.linalg.norm(particle1.velocity - particle2.velocity)
            thermal_velocity = np.sqrt(8 * self.R * temperature / np.pi /
                                     ((particle1.mass + particle2.mass) / 2))
            energy_factor = min(1.0, v_rel_mag / thermal_velocity)

        return base_rate * freq_factor * energy_factor

    def temperature_coefficient(self, temperature: float) -> float:
        """
        Calculate temperature coefficient Q10 at given temperature.

        Q10 = k(T+10) / k(T) shows how much rate increases per 10K.

        Args:
            temperature: Temperature in Kelvin

        Returns:
            Q10 coefficient
        """
        k_t = self.rate_constant(temperature)
        k_t_plus_10 = self.rate_constant(temperature + 10)

        if k_t == 0:
            return float('inf') if k_t_plus_10 > 0 else 1.0

        return k_t_plus_10 / k_t

    def arrhenius_plot_data(self, temp_range: tuple, num_points: int = 50):
        """
        Generate data for Arrhenius plot (ln(k) vs 1/T).

        Args:
            temp_range: (T_min, T_max) temperature range in Kelvin
            num_points: Number of data points to generate

        Returns:
            tuple: (1/T array, ln(k) array) for plotting
        """
        temperatures = np.linspace(temp_range[0], temp_range[1], num_points)
        inv_temps = 1.0 / temperatures

        ln_rates = []
        for T in temperatures:
            k = self.rate_constant(T)
            if k > 0:
                ln_rates.append(np.log(k))
            else:
                ln_rates.append(-np.inf)

        return inv_temps, np.array(ln_rates)

    def __repr__(self) -> str:
        return (f"ArrheniusKinetics(E_a={self.E_a:.2f} kJ/mol, "
                f"A={self.A:.2e} 1/ps, steric={self.steric_factor:.3f})")
