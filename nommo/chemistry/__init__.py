"""
Chemistry module for molecular reactions and bond dynamics.

This module implements realistic chemical kinetics based on:
- Arrhenius equation for temperature-dependent reaction rates
- Collision theory for activation energy barriers
- Thermodynamic principles for bond stability

References:
- Steinfeld et al. (1999) Chemical Kinetics and Dynamics
- Atkins & de Paula (2017) Physical Chemistry
"""

from .bonds import Bond, BondManager
from .kinetics import ArrheniusKinetics
from .reactions import ReactionEngine, ReactionRule
from .thermodynamics import ThermodynamicsCalculator

__all__ = [
    "ArrheniusKinetics",
    "Bond",
    "BondManager",
    "ReactionRule",
    "ReactionEngine",
    "ThermodynamicsCalculator",
]
