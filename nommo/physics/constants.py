"""
Physical constants and conversion factors for molecular dynamics.

All constants in SI-derived units appropriate for nanoscale simulations.
"""

import numpy as np

KB_KJMOL_K = 0.00831446
R_GAS = 8.314462618

NA = 6.02214076e23

AMU_TO_KG = 1.66053906660e-27
NM_TO_M = 1e-9
PS_TO_S = 1e-12

ENERGY_CONVERSION = 0.01

COULOMB_CONSTANT = 138.935

DEFAULT_TIMESTEP = 0.001
DEFAULT_TEMPERATURE = 300.0
DEFAULT_CUTOFF = 1.0

MAX_FORCE = 1e6
MIN_DISTANCE = 0.01

EPSILON_WATER = 1.0
SIGMA_WATER = 0.316


def temperature_to_energy(temperature: float) -> float:
    """Convert temperature (K) to thermal energy (kJ/mol)."""
    return KB_KJMOL_K * temperature


def energy_to_temperature(energy: float) -> float:
    """Convert thermal energy (kJ/mol) to temperature (K)."""
    return energy / KB_KJMOL_K


def velocity_from_temperature(temperature: float, mass: float) -> float:
    """
    Calculate RMS velocity from temperature.

    Args:
        temperature: Temperature in K
        mass: Mass in amu

    Returns:
        RMS velocity in nm/ps
    """
    return float(np.sqrt(3 * KB_KJMOL_K * temperature / (mass * ENERGY_CONVERSION)))
