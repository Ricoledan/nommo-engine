"""
Bond formation and breaking based on thermodynamics and kinetics.

This module implements realistic chemical bonding using:
- Thermodynamic favorability (Gibbs free energy)
- Kinetic barriers (activation energy)
- Distance-dependent bond strength
- Thermal bond breaking

References:
- Atkins & de Paula (2017) Physical Chemistry
- McQuarrie & Simon (1997) Molecular Thermodynamics
"""

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    from nommo.core.particle import Particle

from nommo.utils.logging import get_logger

from .kinetics import ArrheniusKinetics

logger = get_logger("bonds")


@dataclass
class Bond:
    """
    Represents a chemical bond between two particles.

    Bonds have thermodynamic and kinetic properties that determine
    their stability and reactivity.
    """

    particle1_id: str
    particle2_id: str
    bond_energy: float  # kJ/mol - positive for stable bonds
    bond_length: float  # nm - equilibrium distance
    age: int = 0  # simulation ticks since formation
    bond_type: str = "single"  # single, double, triple, etc.

    # Dynamic properties
    current_length: float | None = None  # current distance
    strain_energy: float = 0.0  # excess energy due to stretching/compression

    def __post_init__(self):
        """Validate bond parameters."""
        if self.bond_energy < 0:
            logger.warning(f"Bond with negative energy: {self.bond_energy} kJ/mol")
        if self.bond_length <= 0:
            raise ValueError("Bond length must be positive")

    def update_geometry(self, distance: float):
        """
        Update bond geometry and calculate strain energy.

        Uses harmonic potential for bond stretching:
        E_strain = (1/2) * k * (r - r0)²

        Args:
            distance: Current distance between bonded particles
        """
        self.current_length = distance

        # Spring constant scales with bond energy
        # Typical C-C bond: k ≈ 1000 kJ/(mol·nm²)
        spring_constant = self.bond_energy * 50  # Rough approximation

        length_deviation = distance - self.bond_length
        self.strain_energy = 0.5 * spring_constant * length_deviation**2

    def is_overstretched(self, breaking_threshold: float = 2.0) -> bool:
        """
        Check if bond is stretched beyond breaking point.

        Args:
            breaking_threshold: Multiple of equilibrium length for breaking

        Returns:
            True if bond should break due to overstretching
        """
        if self.current_length is None:
            return False
        return self.current_length > breaking_threshold * self.bond_length

    def total_energy(self) -> float:
        """Calculate total bond energy including strain."""
        return self.bond_energy + self.strain_energy

    def __str__(self) -> str:
        return (
            f"Bond({self.particle1_id[:8]}↔{self.particle2_id[:8]}, {self.bond_energy:.1f} kJ/mol)"
        )


class BondManager:
    """
    Manages bond formation and breaking in the simulation.

    This class implements realistic chemical bonding based on:
    1. Thermodynamic favorability (ΔG < 0)
    2. Kinetic accessibility (activation energy)
    3. Geometric constraints (distance, angles)
    4. Thermal stability (bond dissociation)

    Formation criteria:
    - Particles within bonding distance
    - Sufficient collision energy (> E_activation)
    - Thermodynamically favorable (ΔG < 0)
    - Available bonding sites (valency)

    Breaking criteria:
    - Thermal fluctuation (Boltzmann probability)
    - Energetic collision (> bond energy)
    - Mechanical overstretching
    - External energy input
    """

    def __init__(
        self,
        bond_energy: float = 20.0,  # kJ/mol
        activation_energy: float = 10.0,  # kJ/mol
        bonding_distance: float = 0.3,  # nm
        breaking_threshold: float = 2.0,  # multiple of bond length
        thermal_breaking_rate: float = 1e-6,  # base probability per timestep
    ):
        """
        Initialize bond manager.

        Args:
            bond_energy: Standard bond strength (kJ/mol)
            activation_energy: Energy barrier for bond formation (kJ/mol)
            bonding_distance: Maximum distance for bond formation (nm)
            breaking_threshold: Distance threshold for mechanical breaking
            thermal_breaking_rate: Base rate for thermal bond breaking
        """
        self.bond_energy = bond_energy
        self.activation_energy = activation_energy
        self.bonding_distance = bonding_distance
        self.breaking_threshold = breaking_threshold
        self.thermal_breaking_rate = thermal_breaking_rate

        # Physical constants
        self.k_B = 0.00831446  # Boltzmann constant in kJ/(mol·K)

        # Kinetics calculator for bond formation
        self.formation_kinetics = ArrheniusKinetics(
            activation_energy=activation_energy,
            pre_exponential=1e12,  # Typical for bond formation
            steric_factor=0.1,  # Account for orientation requirements
        )

        # Track all bonds in the system
        self.bonds: dict[str, Bond] = {}  # bond_id -> Bond
        self.particle_bonds: dict[str, set[str]] = {}  # particle_id -> set of bond_ids

        logger.debug(
            f"Initialized BondManager: E_bond={bond_energy} kJ/mol, "
            f"E_act={activation_energy} kJ/mol, d_bond={bonding_distance} nm"
        )

    def can_form_bond(self, p1: "Particle", p2: "Particle", temperature: float) -> bool:
        """
        Check if two particles can form a bond.

        Args:
            p1: First particle
            p2: Second particle
            temperature: System temperature (K)

        Returns:
            True if bond formation is possible
        """
        # Check basic requirements
        if not p1.can_bond_with(p2):
            return False

        # Check distance constraint
        distance = p1.distance_to(p2)
        if distance > self.bonding_distance:
            return False

        # Check if already bonded
        bond_id = self._get_bond_id(p1.id, p2.id)
        if bond_id in self.bonds:
            return False

        # Check activation energy requirement
        if not self.formation_kinetics.has_sufficient_energy(p1, p2, temperature):
            return False

        # Check thermodynamic favorability
        return self._is_thermodynamically_favorable(p1, p2, temperature)

    def attempt_bond_formation(
        self, p1: "Particle", p2: "Particle", temperature: float, timestep: float
    ) -> Bond | None:
        """
        Attempt to form a bond between two particles.

        Args:
            p1: First particle
            p2: Second particle
            temperature: System temperature (K)
            timestep: Time interval (ps)

        Returns:
            New Bond if successful, None otherwise
        """
        if not self.can_form_bond(p1, p2, temperature):
            return None

        # Calculate formation probability
        formation_prob = self.formation_kinetics.reaction_probability(temperature, timestep)

        # Check stochastic formation
        if np.random.random() > formation_prob:
            return None

        # Create bond
        bond = self._create_bond(p1, p2)

        # Update particle bond lists
        p1.add_bond(p2)
        p2.add_bond(p1)

        # Register in manager
        bond_id = self._get_bond_id(p1.id, p2.id)
        self.bonds[bond_id] = bond

        # Update particle bond tracking
        if p1.id not in self.particle_bonds:
            self.particle_bonds[p1.id] = set()
        if p2.id not in self.particle_bonds:
            self.particle_bonds[p2.id] = set()

        self.particle_bonds[p1.id].add(bond_id)
        self.particle_bonds[p2.id].add(bond_id)

        logger.debug(f"Bond formed: {bond}")
        return bond

    def should_break_bond(
        self, bond: Bond, p1: "Particle", p2: "Particle", temperature: float
    ) -> bool:
        """
        Determine if a bond should break.

        Args:
            bond: Bond to evaluate
            p1: First bonded particle
            p2: Second bonded particle
            temperature: System temperature (K)

        Returns:
            True if bond should break
        """
        # Update bond geometry
        distance = p1.distance_to(p2)
        bond.update_geometry(distance)

        # Check mechanical breaking (overstretching)
        if bond.is_overstretched(self.breaking_threshold):
            logger.debug(
                f"Bond breaking due to overstretching: {distance:.3f} nm > "
                f"{self.breaking_threshold * bond.bond_length:.3f} nm"
            )
            return True

        # Check energetic breaking (high-energy collision)
        if self._has_breaking_collision(bond, p1, p2, temperature):
            logger.debug("Bond breaking due to energetic collision")
            return True

        # Check thermal breaking (random dissociation)
        if self._thermal_breaking_probability(bond, temperature):
            logger.debug("Bond breaking due to thermal fluctuation")
            return True

        return False

    def break_bond(self, bond_id: str, p1: "Particle", p2: "Particle") -> bool:
        """
        Break a specific bond.

        Args:
            bond_id: ID of bond to break
            p1: First bonded particle
            p2: Second bonded particle

        Returns:
            True if bond was broken successfully
        """
        if bond_id not in self.bonds:
            return False

        # Remove from particles
        p1.remove_bond(p2.id)
        p2.remove_bond(p1.id)

        # Remove from manager
        bond = self.bonds.pop(bond_id)

        # Update tracking
        if p1.id in self.particle_bonds:
            self.particle_bonds[p1.id].discard(bond_id)
        if p2.id in self.particle_bonds:
            self.particle_bonds[p2.id].discard(bond_id)

        logger.debug(f"Bond broken: {bond}")
        return True

    def update_bonds(
        self, particles: list["Particle"], temperature: float, timestep: float
    ) -> dict[str, int]:
        """
        Update all bonds in the system.

        Args:
            particles: All particles in system
            temperature: System temperature (K)
            timestep: Time interval (ps)

        Returns:
            Dictionary with formation/breaking statistics
        """
        particle_dict = {p.id: p for p in particles}
        bonds_formed = 0
        bonds_broken = 0

        # Age existing bonds and check for breaking
        bonds_to_break = []
        for bond_id, bond in self.bonds.items():
            bond.age += 1

            p1 = particle_dict.get(bond.particle1_id)
            p2 = particle_dict.get(bond.particle2_id)

            if p1 is None or p2 is None:
                # Particle was removed from system
                bonds_to_break.append(bond_id)
                continue

            if self.should_break_bond(bond, p1, p2, temperature):
                bonds_to_break.append(bond_id)

        # Break bonds marked for breaking
        for bond_id in bonds_to_break:
            bond = self.bonds[bond_id]
            p1 = particle_dict.get(bond.particle1_id)
            p2 = particle_dict.get(bond.particle2_id)
            if p1 and p2:
                self.break_bond(bond_id, p1, p2)
                bonds_broken += 1

        return {
            "bonds_formed": bonds_formed,
            "bonds_broken": bonds_broken,
            "total_bonds": len(self.bonds),
        }

    def get_bond_network(self) -> dict[str, list[str]]:
        """
        Get adjacency list representation of bond network.

        Returns:
            Dictionary mapping particle_id -> list of bonded particle_ids
        """
        network = {}
        for bond in self.bonds.values():
            # Add edges in both directions
            if bond.particle1_id not in network:
                network[bond.particle1_id] = []
            if bond.particle2_id not in network:
                network[bond.particle2_id] = []

            network[bond.particle1_id].append(bond.particle2_id)
            network[bond.particle2_id].append(bond.particle1_id)

        return network

    def get_particle_bonds(self, particle_id: str) -> list[Bond]:
        """Get all bonds for a specific particle."""
        bonds = []
        if particle_id in self.particle_bonds:
            for bond_id in self.particle_bonds[particle_id]:
                if bond_id in self.bonds:
                    bonds.append(self.bonds[bond_id])
        return bonds

    def _get_bond_id(self, p1_id: str, p2_id: str) -> str:
        """Generate consistent bond ID from particle IDs."""
        return f"{min(p1_id, p2_id)}↔{max(p1_id, p2_id)}"

    def _create_bond(self, p1: "Particle", p2: "Particle") -> Bond:
        """Create a new bond between particles."""
        distance = p1.distance_to(p2)

        # Bond strength depends on particle types
        bond_energy = self._calculate_bond_energy(p1, p2)

        return Bond(
            particle1_id=p1.id,
            particle2_id=p2.id,
            bond_energy=bond_energy,
            bond_length=distance * 0.9,  # Slightly shorter than formation distance
            bond_type="single",
        )

    def _calculate_bond_energy(self, p1: "Particle", p2: "Particle") -> float:
        """Calculate bond energy based on particle types."""
        # Base bond energy
        energy = self.bond_energy

        # Modify based on particle types
        type_factors = {"monomer_A": 1.0, "monomer_B": 1.2, "catalyst": 0.8, "default": 1.0}

        factor1 = type_factors.get(p1.particle_type, 1.0)
        factor2 = type_factors.get(p2.particle_type, 1.0)

        # Geometric mean for bond strength
        return energy * np.sqrt(factor1 * factor2)

    def _is_thermodynamically_favorable(
        self, p1: "Particle", p2: "Particle", temperature: float
    ) -> bool:
        """
        Check if bond formation is thermodynamically favorable.

        Simple approximation: ΔG = ΔH - TΔS
        Assume bond formation always decreases entropy (ΔS < 0)
        """
        # Enthalpy change (negative for bond formation)
        delta_H = -self._calculate_bond_energy(p1, p2)

        # Entropy change (negative for bond formation - loss of freedom)
        delta_S = -0.1  # kJ/(mol·K) - rough estimate

        # Gibbs free energy change
        delta_G = delta_H - temperature * delta_S

        # Favorable if ΔG < 0
        return delta_G < 0

    def _has_breaking_collision(
        self, bond: Bond, p1: "Particle", p2: "Particle", temperature: float
    ) -> bool:
        """Check if collision has enough energy to break bond."""
        # Relative kinetic energy
        v_rel = p1.velocity - p2.velocity
        v_rel_mag = np.linalg.norm(v_rel)

        m_reduced = (p1.mass * p2.mass) / (p1.mass + p2.mass)
        collision_energy = 0.5 * m_reduced * v_rel_mag**2

        # Convert to kJ/mol
        collision_energy_kjmol = collision_energy * 0.01 * 6.022e23 / 1000

        # Add strain energy contribution
        total_energy_barrier = bond.total_energy()

        return collision_energy_kjmol > total_energy_barrier

    def _thermal_breaking_probability(self, bond: Bond, temperature: float) -> bool:
        """Calculate probability of thermal bond breaking."""
        # Base thermal breaking rate
        base_rate = self.thermal_breaking_rate

        # Arrhenius factor for bond dissociation
        dissociation_energy = bond.total_energy()
        arrhenius_factor = np.exp(-dissociation_energy / (self.k_B * temperature))

        # Total breaking probability
        breaking_prob = base_rate * arrhenius_factor

        return np.random.random() < breaking_prob

    def __len__(self) -> int:
        """Number of bonds in the system."""
        return len(self.bonds)

    def __repr__(self) -> str:
        return f"BondManager({len(self.bonds)} bonds, E_bond={self.bond_energy} kJ/mol)"
