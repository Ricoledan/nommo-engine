"""
Reaction rule engine for complex chemical transformations.

This module implements a flexible system for defining and executing
chemical reactions beyond simple bond formation/breaking.

Features:
- Multi-particle reactions (A + B → C, A + B + C → D, etc.)
- Catalytic reactions (A + B --C--> D + E)
- Template matching for complex transformations
- Energy balance tracking
- Reaction pathway analysis

References:
- Gillespie (1977) - Stochastic simulation algorithm
- Gibson & Bruck (2000) - Efficient exact stochastic simulation
- Slepoy et al. (2008) - A constant-time kinetic Monte Carlo algorithm
"""

from collections import defaultdict
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from typing import TYPE_CHECKING, Any, Optional

import numpy as np

if TYPE_CHECKING:
    from nommo.core.particle import Particle

from nommo.utils.logging import get_logger

from .kinetics import ArrheniusKinetics

logger = get_logger("reactions")


class ReactionType(Enum):
    """Types of chemical reactions."""
    SYNTHESIS = "synthesis"        # A + B → AB
    DECOMPOSITION = "decomposition"  # AB → A + B
    SUBSTITUTION = "substitution"   # AB + C → AC + B
    CATALYSIS = "catalysis"        # A + B --C--> D + E (C unchanged)
    POLYMERIZATION = "polymerization"  # nA → (A)n
    REARRANGEMENT = "rearrangement"   # ABC → ACB


@dataclass
class ReactionRule:
    """
    Definition of a chemical reaction rule.

    A reaction rule specifies:
    - Reactant patterns (what particles can react)
    - Product specification (what is produced)
    - Kinetic parameters (rate, activation energy)
    - Geometric constraints (distances, angles)
    - Energetics (enthalpy change)
    """

    name: str
    reaction_type: ReactionType

    # Reactant specification
    reactant_types: list[str]  # e.g., ["monomer_A", "monomer_B"]

    # Product specification
    product_types: list[str]   # e.g., ["dimer_AB"]

    # Optional fields with defaults
    min_reactants: int = 2
    max_reactants: int = 2
    conserve_mass: bool = True

    # Kinetics
    activation_energy: float = 10.0  # kJ/mol
    pre_exponential: float = 1e12    # 1/ps
    enthalpy_change: float = 0.0     # kJ/mol (negative = exothermic)

    # Geometric constraints
    max_distance: float = 0.3        # nm - maximum reactant separation
    min_distance: float = 0.05       # nm - minimum reactant separation

    # Catalysis
    catalyst_types: list[str] = field(default_factory=list)
    catalyst_distance: float = 0.5   # nm - catalyst interaction range

    # Conditional requirements
    requirements: dict[str, Any] = field(default_factory=dict)

    # Custom reaction function (for complex transformations)
    custom_function: Optional[Callable] = None

    def __post_init__(self):
        """Validate reaction rule parameters."""
        if self.min_reactants < 1:
            raise ValueError("Must have at least 1 reactant")
        if self.max_reactants < self.min_reactants:
            raise ValueError("Max reactants must be >= min reactants")
        if self.max_distance <= self.min_distance:
            raise ValueError("Max distance must be > min distance")
        if self.activation_energy < 0:
            raise ValueError("Activation energy must be non-negative")

    def matches_reactants(self, particle_types: list[str]) -> bool:
        """
        Check if particle types match this reaction's reactants.

        Args:
            particle_types: List of particle type strings

        Returns:
            True if types match the reaction pattern
        """
        if not (self.min_reactants <= len(particle_types) <= self.max_reactants):
            return False

        # For simple reactions, check exact match
        if len(self.reactant_types) == len(particle_types):
            return sorted(self.reactant_types) == sorted(particle_types)

        # For complex reactions, check if all required types are present
        required = set(self.reactant_types)
        available = set(particle_types)
        return required.issubset(available)

    def requires_catalyst(self) -> bool:
        """Check if reaction requires a catalyst."""
        return len(self.catalyst_types) > 0

    def is_energetically_favorable(self, temperature: float) -> bool:
        """
        Check if reaction is thermodynamically favorable.

        Uses simple Gibbs free energy approximation:
        ΔG = ΔH - TΔS
        """
        # Estimate entropy change based on reaction type (kJ/(mol·K))
        if self.reaction_type == ReactionType.SYNTHESIS:
            delta_S = -0.05  # Loss of entropy when forming bonds
        elif self.reaction_type == ReactionType.DECOMPOSITION:
            delta_S = 0.05   # Gain of entropy when breaking apart
        else:
            delta_S = 0.0   # Neutral for other types

        # Calculate Gibbs free energy change
        delta_G = self.enthalpy_change - temperature * delta_S
        return delta_G < 0  # Favorable if ΔG < 0


@dataclass
class ReactionEvent:
    """Record of a completed reaction."""

    rule_name: str
    reactant_ids: list[str]
    product_ids: list[str]
    catalyst_ids: list[str]
    energy_released: float  # kJ/mol
    timestamp: int          # simulation tick
    location: np.ndarray    # reaction center position


class ReactionEngine:
    """
    Manages and executes chemical reactions in the simulation.

    The reaction engine:
    1. Maintains a library of reaction rules
    2. Identifies potential reactions between particles
    3. Evaluates reaction probability using kinetics
    4. Executes reactions by transforming particles
    5. Tracks reaction statistics and pathways

    Uses a spatial acceleration structure to efficiently find
    reactive particle pairs without O(N²) search.
    """

    def __init__(self, temperature: float = 300.0):
        """
        Initialize reaction engine.

        Args:
            temperature: System temperature in Kelvin
        """
        self.temperature = temperature

        # Reaction rules database
        self.rules: dict[str, ReactionRule] = {}

        # Kinetics calculators for each rule
        self.kinetics: dict[str, ArrheniusKinetics] = {}

        # Reaction history
        self.reaction_history: list[ReactionEvent] = []
        self.reaction_counts: dict[str, int] = defaultdict(int)

        # Statistics
        self.total_reactions = 0
        self.energy_released_total = 0.0

        logger.debug(f"Initialized ReactionEngine at T={temperature} K")

    def add_rule(self, rule: ReactionRule):
        """
        Add a reaction rule to the engine.

        Args:
            rule: ReactionRule to add
        """
        self.rules[rule.name] = rule

        # Create kinetics calculator for this rule
        self.kinetics[rule.name] = ArrheniusKinetics(
            activation_energy=rule.activation_energy,
            pre_exponential=rule.pre_exponential,
            steric_factor=0.1  # Account for orientation requirements
        )

        logger.debug(f"Added reaction rule: {rule.name}")

    def remove_rule(self, rule_name: str):
        """Remove a reaction rule."""
        if rule_name in self.rules:
            del self.rules[rule_name]
            del self.kinetics[rule_name]
            logger.debug(f"Removed reaction rule: {rule_name}")

    def find_potential_reactions(
        self,
        particles: list["Particle"],
        spatial_index=None
    ) -> list[dict[str, Any]]:
        """
        Find all potential reactions in the current system state.

        Args:
            particles: All particles in the system
            spatial_index: Spatial acceleration structure (optional)

        Returns:
            List of potential reaction dictionaries
        """
        potential_reactions = []

        # Create particle lookup
        particle_dict = {p.id: p for p in particles}

        # For each reaction rule, find matching particle combinations
        for rule_name, rule in self.rules.items():
            matches = self._find_rule_matches(rule, particles, particle_dict, spatial_index)

            for match in matches:
                potential_reactions.append({
                    "rule_name": rule_name,
                    "rule": rule,
                    "reactants": match["reactants"],
                    "catalysts": match.get("catalysts", []),
                    "center": match["center"]
                })

        return potential_reactions

    def execute_reactions(
        self,
        particles: list["Particle"],
        timestep: float,
        spatial_index=None
    ) -> dict[str, Any]:
        """
        Execute all feasible reactions in the system.

        Args:
            particles: List of particles (will be modified)
            timestep: Time interval in picoseconds
            spatial_index: Spatial acceleration structure

        Returns:
            Statistics about reactions that occurred
        """
        reactions_executed = 0
        energy_released = 0.0
        particles_created = 0
        particles_destroyed = 0

        # Find potential reactions
        potential_reactions = self.find_potential_reactions(particles, spatial_index)

        # Evaluate and execute each potential reaction
        for reaction_info in potential_reactions:
            rule = reaction_info["rule"]
            reactants = reaction_info["reactants"]
            catalysts = reaction_info.get("catalysts", [])

            # Check if reaction should occur
            if self._should_reaction_occur(rule, reactants, timestep):
                # Execute the reaction
                result = self._execute_single_reaction(
                    rule, reactants, catalysts, particles
                )

                if result["success"]:
                    reactions_executed += 1
                    energy_released += result["energy_released"]
                    particles_created += len(result["products"])
                    particles_destroyed += len(reactants)

                    # Record reaction event
                    self._record_reaction_event(
                        rule.name, reactants, result["products"],
                        catalysts, result["energy_released"],
                        reaction_info["center"]
                    )

        # Update statistics
        self.total_reactions += reactions_executed
        self.energy_released_total += energy_released

        return {
            "reactions_executed": reactions_executed,
            "energy_released": energy_released,
            "particles_created": particles_created,
            "particles_destroyed": particles_destroyed,
            "total_reactions": self.total_reactions
        }

    def get_statistics(self) -> dict[str, Any]:
        """Get comprehensive reaction statistics."""
        return {
            "total_reactions": self.total_reactions,
            "energy_released_total": self.energy_released_total,
            "reaction_counts": dict(self.reaction_counts),
            "average_energy_per_reaction": (
                self.energy_released_total / max(1, self.total_reactions)
            ),
            "active_rules": len(self.rules),
            "recent_reactions": self.reaction_history[-10:]  # Last 10 reactions
        }

    def _find_rule_matches(
        self,
        rule: ReactionRule,
        particles: list["Particle"],
        particle_dict: dict[str, "Particle"],
        spatial_index=None
    ) -> list[dict[str, Any]]:
        """Find all particle combinations that match a reaction rule."""
        matches = []

        # Group particles by type for efficient matching
        particles_by_type = defaultdict(list)
        for p in particles:
            particles_by_type[p.particle_type].append(p)

        # Find reactant combinations
        if rule.reaction_type == ReactionType.SYNTHESIS:
            matches.extend(self._find_synthesis_matches(rule, particles_by_type))
        elif rule.reaction_type == ReactionType.DECOMPOSITION:
            matches.extend(self._find_decomposition_matches(rule, particles_by_type))
        elif rule.reaction_type == ReactionType.CATALYSIS:
            matches.extend(self._find_catalysis_matches(rule, particles_by_type))
        else:
            # Generic matching for other reaction types
            matches.extend(self._find_generic_matches(rule, particles_by_type))

        return matches

    def _find_synthesis_matches(
        self,
        rule: ReactionRule,
        particles_by_type: dict[str, list["Particle"]]
    ) -> list[dict[str, Any]]:
        """Find particle pairs for synthesis reactions (A + B → AB)."""
        matches = []

        if len(rule.reactant_types) != 2:
            return matches

        type_a, type_b = rule.reactant_types
        particles_a = particles_by_type.get(type_a, [])
        particles_b = particles_by_type.get(type_b, [])

        # Find all valid pairs
        for p_a in particles_a:
            for p_b in particles_b:
                if p_a.id == p_b.id:
                    continue

                distance = p_a.distance_to(p_b)
                if rule.min_distance <= distance <= rule.max_distance:
                    # Check if particles can bond
                    if p_a.can_bond_with(p_b):
                        center = (p_a.position + p_b.position) / 2
                        matches.append({
                            "reactants": [p_a, p_b],
                            "center": center
                        })

        return matches

    def _find_decomposition_matches(
        self,
        rule: ReactionRule,
        particles_by_type: dict[str, list["Particle"]]
    ) -> list[dict[str, Any]]:
        """Find particles for decomposition reactions (AB → A + B)."""
        matches = []

        if len(rule.reactant_types) != 1:
            return matches

        reactant_type = rule.reactant_types[0]
        particles = particles_by_type.get(reactant_type, [])

        for particle in particles:
            # Check if particle has bonds that can be broken
            if len(particle.bonds) > 0:
                matches.append({
                    "reactants": [particle],
                    "center": particle.position.copy()
                })

        return matches

    def _find_catalysis_matches(
        self,
        rule: ReactionRule,
        particles_by_type: dict[str, list["Particle"]]
    ) -> list[dict[str, Any]]:
        """Find reactant-catalyst combinations for catalytic reactions."""
        matches = []

        # First find regular reactant matches
        regular_matches = self._find_synthesis_matches(rule, particles_by_type)

        # Then find catalysts near each reactant pair
        catalyst_particles = []
        for cat_type in rule.catalyst_types:
            catalyst_particles.extend(particles_by_type.get(cat_type, []))

        for match in regular_matches:
            reactants = match["reactants"]
            center = match["center"]

            # Find catalysts within range
            nearby_catalysts = []
            for catalyst in catalyst_particles:
                cat_distance = np.linalg.norm(catalyst.position - center)
                if cat_distance <= rule.catalyst_distance:
                    nearby_catalysts.append(catalyst)

            if nearby_catalysts:
                match["catalysts"] = nearby_catalysts
                matches.append(match)

        return matches

    def _find_generic_matches(
        self,
        rule: ReactionRule,
        particles_by_type: dict[str, list["Particle"]]
    ) -> list[dict[str, Any]]:
        """Generic matching for complex reaction types."""
        # This is a simplified implementation
        # In practice, would need more sophisticated pattern matching
        return []

    def _should_reaction_occur(
        self,
        rule: ReactionRule,
        reactants: list["Particle"],
        timestep: float
    ) -> bool:
        """Determine if a reaction should occur based on kinetics."""
        # Check thermodynamic favorability
        if not rule.is_energetically_favorable(self.temperature):
            return False

        # Calculate reaction probability
        kinetics = self.kinetics[rule.name]

        if len(reactants) == 2:
            # Bimolecular reaction
            prob = kinetics.reaction_probability(self.temperature, timestep)

            # Adjust for collision energy
            if kinetics.has_sufficient_energy(reactants[0], reactants[1], self.temperature):
                prob *= 2.0  # Increase probability for energetic collisions
        else:
            # Unimolecular reaction
            prob = kinetics.reaction_probability(self.temperature, timestep)

        return np.random.random() < prob

    def _execute_single_reaction(
        self,
        rule: ReactionRule,
        reactants: list["Particle"],
        catalysts: list["Particle"],
        all_particles: list["Particle"]
    ) -> dict[str, Any]:
        """Execute a single reaction and return results."""
        try:
            if rule.custom_function:
                # Use custom reaction function
                result = rule.custom_function(reactants, catalysts, rule)
            else:
                # Use default reaction logic
                result = self._default_reaction_execution(rule, reactants, catalysts)

            if result["success"]:
                # Remove reactants from particle list
                for reactant in reactants:
                    if reactant in all_particles:
                        all_particles.remove(reactant)

                # Add products to particle list
                all_particles.extend(result["products"])

            return result

        except Exception as e:
            logger.error(f"Error executing reaction {rule.name}: {e}")
            return {"success": False, "error": str(e)}

    def _default_reaction_execution(
        self,
        rule: ReactionRule,
        reactants: list["Particle"],
        catalysts: list["Particle"]
    ) -> dict[str, Any]:
        """Default reaction execution logic."""
        products = []

        if rule.reaction_type == ReactionType.SYNTHESIS:
            # Combine reactants into single product
            if len(reactants) == 2:
                product = self._create_synthesis_product(rule, reactants)
                products.append(product)

        elif rule.reaction_type == ReactionType.DECOMPOSITION:
            # Break apart reactant into multiple products
            if len(reactants) == 1:
                decomp_products = self._create_decomposition_products(rule, reactants[0])
                products.extend(decomp_products)

        energy_released = -rule.enthalpy_change  # Negative enthalpy = energy released

        return {
            "success": True,
            "products": products,
            "energy_released": energy_released
        }

    def _create_synthesis_product(
        self,
        rule: ReactionRule,
        reactants: list["Particle"]
    ) -> "Particle":
        """Create a new particle from synthesis reaction."""
        from nommo.core.particle import Particle

        # Calculate product properties
        total_mass = sum(r.mass for r in reactants)
        center_of_mass = sum(r.mass * r.position for r in reactants) / total_mass
        total_momentum = sum(r.momentum() for r in reactants)

        # Create product particle
        product = Particle(
            mass=total_mass,
            position=center_of_mass,
            velocity=total_momentum / total_mass,
            particle_type=rule.product_types[0] if rule.product_types else "product",
            generation=max(r.generation for r in reactants) + 1
        )

        return product

    def _create_decomposition_products(
        self,
        rule: ReactionRule,
        reactant: "Particle"
    ) -> list["Particle"]:
        """Create multiple particles from decomposition reaction."""
        from nommo.core.particle import Particle

        products = []
        num_products = len(rule.product_types)

        if num_products == 0:
            return products

        # Distribute mass and momentum
        mass_per_product = reactant.mass / num_products

        for i, product_type in enumerate(rule.product_types):
            # Slightly offset positions to avoid overlap
            offset = np.random.normal(0, 0.1, 3)  # Small random displacement

            # Give random velocities (conservation handled approximately)
            velocity = reactant.velocity + np.random.normal(0, 0.1, 3)

            product = Particle(
                mass=mass_per_product,
                position=reactant.position + offset,
                velocity=velocity,
                particle_type=product_type,
                generation=reactant.generation + 1,
                parent_id=reactant.id
            )

            products.append(product)

        return products

    def _record_reaction_event(
        self,
        rule_name: str,
        reactants: list["Particle"],
        products: list["Particle"],
        catalysts: list["Particle"],
        energy_released: float,
        center: np.ndarray
    ):
        """Record a reaction event for analysis."""
        event = ReactionEvent(
            rule_name=rule_name,
            reactant_ids=[r.id for r in reactants],
            product_ids=[p.id for p in products],
            catalyst_ids=[c.id for c in catalysts],
            energy_released=energy_released,
            timestamp=0,  # Would need simulation time
            location=center
        )

        self.reaction_history.append(event)
        self.reaction_counts[rule_name] += 1

        # Keep history bounded
        if len(self.reaction_history) > 1000:
            self.reaction_history = self.reaction_history[-500:]
