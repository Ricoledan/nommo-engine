"""
Thermodynamics calculations for chemical systems.

This module implements thermodynamic principles for:
- Gibbs free energy calculations
- Energy flow tracking
- Bond network analysis
- Autocatalytic set detection

References:
- Atkins & de Paula (2017) Physical Chemistry
- Kauffman (1986) - Autocatalytic sets
- Hordijk & Steel (2017) - Detection algorithms
"""

from collections import defaultdict
from dataclasses import dataclass
from typing import TYPE_CHECKING

import networkx as nx
import numpy as np

if TYPE_CHECKING:
    from nommo.core.particle import Particle

    from .bonds import Bond

from nommo.utils.logging import get_logger

logger = get_logger("thermodynamics")


@dataclass
class EnergyFlow:
    """Track energy flow in the system."""

    kinetic_energy: float = 0.0  # kJ/mol
    potential_energy: float = 0.0  # kJ/mol
    bond_energy: float = 0.0  # kJ/mol
    reaction_energy: float = 0.0  # kJ/mol released from reactions
    input_energy: float = 0.0  # kJ/mol added externally

    @property
    def total_energy(self) -> float:
        """Total energy in the system."""
        return (
            self.kinetic_energy
            + self.potential_energy
            + self.bond_energy
            + self.reaction_energy
            + self.input_energy
        )

    def to_dict(self) -> dict[str, float]:
        """Convert to dictionary for analysis."""
        return {
            "kinetic": self.kinetic_energy,
            "potential": self.potential_energy,
            "bond": self.bond_energy,
            "reaction": self.reaction_energy,
            "input": self.input_energy,
            "total": self.total_energy,
        }


@dataclass
class NetworkCluster:
    """A connected component in the bond network."""

    particle_ids: set[str]
    bonds: list["Bond"]
    size: int
    complexity: float  # Measure of structural complexity
    is_cyclic: bool  # Contains cycles

    @property
    def bond_count(self) -> int:
        """Number of bonds in cluster."""
        return len(self.bonds)

    @property
    def density(self) -> float:
        """Bond density: edges / max_possible_edges."""
        n = len(self.particle_ids)
        if n <= 1:
            return 0.0
        max_edges = n * (n - 1) / 2
        return len(self.bonds) / max_edges


class ThermodynamicsCalculator:
    """
    Calculate thermodynamic properties and analyze molecular networks.

    This class provides:
    1. Energy balance calculations
    2. Bond network analysis
    3. Autocatalytic set detection
    4. Complexity metrics
    5. Emergence indicators
    """

    def __init__(self, temperature: float = 300.0):
        """
        Initialize thermodynamics calculator.

        Args:
            temperature: System temperature in Kelvin
        """
        self.temperature = temperature
        self.k_B = 0.00831446  # Boltzmann constant in kJ/(mol·K)

        # Energy tracking
        self.energy_history: list[EnergyFlow] = []
        self.reaction_energy_released = 0.0
        self.bond_formation_energy = 0.0

        logger.debug(f"Initialized ThermodynamicsCalculator at T={temperature} K")

    def calculate_energy_flow(
        self, particles: list["Particle"], bonds: list["Bond"], reaction_energy: float = 0.0
    ) -> EnergyFlow:
        """
        Calculate comprehensive energy flow in the system.

        Args:
            particles: All particles in system
            bonds: All bonds in system
            reaction_energy: Energy released from reactions this timestep

        Returns:
            EnergyFlow object with all energy components
        """
        flow = EnergyFlow()

        # Kinetic energy
        flow.kinetic_energy = sum(p.kinetic_energy() for p in particles)

        # Potential energy (from inter-particle forces)
        flow.potential_energy = self._calculate_potential_energy(particles)

        # Bond energy
        flow.bond_energy = sum(bond.total_energy() for bond in bonds)

        # Reaction energy
        flow.reaction_energy = reaction_energy
        self.reaction_energy_released += reaction_energy

        # Store in history
        self.energy_history.append(flow)

        # Keep history bounded
        if len(self.energy_history) > 1000:
            self.energy_history = self.energy_history[-500:]

        return flow

    def analyze_bond_network(
        self, particles: list["Particle"], bonds: list["Bond"]
    ) -> dict[str, any]:
        """
        Analyze the bond network structure.

        Args:
            particles: All particles in system
            bonds: All bonds in system

        Returns:
            Dictionary with network analysis results
        """
        # Build network graph
        graph = self._build_bond_graph(particles, bonds)

        # Find connected components (clusters)
        clusters = self._find_clusters(graph, bonds)

        # Calculate network metrics
        metrics = self._calculate_network_metrics(graph, clusters)

        # Detect special structures
        autocatalytic_sets = self._detect_autocatalytic_sets(clusters, bonds)
        replicators = self._detect_replicators(particles, clusters)

        return {
            "total_particles": len(particles),
            "total_bonds": len(bonds),
            "clusters": clusters,
            "largest_cluster_size": max(len(c.particle_ids) for c in clusters) if clusters else 0,
            "network_metrics": metrics,
            "autocatalytic_sets": autocatalytic_sets,
            "replicators": replicators,
            "complexity_score": self._calculate_complexity_score(clusters),
        }

    def detect_emergence_signatures(
        self, network_analysis: dict[str, any], energy_flow: EnergyFlow, history_window: int = 100
    ) -> dict[str, any]:
        """
        Detect signatures of emergent behavior.

        Args:
            network_analysis: Results from bond network analysis
            energy_flow: Current energy state
            history_window: Number of timesteps to analyze

        Returns:
            Dictionary with emergence indicators
        """
        signatures = {}

        # Energy flow signatures
        if len(self.energy_history) >= history_window:
            recent_history = self.energy_history[-history_window:]
            signatures["energy_trend"] = self._analyze_energy_trend(recent_history)
            signatures["energy_efficiency"] = self._calculate_energy_efficiency(recent_history)

        # Structural signatures
        signatures["structural_complexity"] = network_analysis["complexity_score"]
        signatures["largest_cluster_growth"] = self._track_cluster_growth(network_analysis)

        # Autocatalytic signatures
        signatures["autocatalytic_activity"] = len(network_analysis["autocatalytic_sets"]) > 0
        signatures["self_replication"] = len(network_analysis["replicators"]) > 0

        # Information-theoretic signatures
        signatures["shannon_entropy"] = self._calculate_shannon_entropy(network_analysis)
        signatures["mutual_information"] = self._calculate_mutual_information(network_analysis)

        return signatures

    def _calculate_potential_energy(self, particles: list["Particle"]) -> float:
        """Calculate total potential energy from inter-particle forces."""
        # This is a simplified calculation
        # In practice, would integrate with force calculations from physics module
        total_pe = 0.0

        for i, p1 in enumerate(particles):
            for p2 in particles[i + 1 :]:
                distance = p1.distance_to(p2)
                if distance < 1.0:  # Within interaction range
                    # Simple Lennard-Jones-like potential
                    sigma = 0.3  # nm
                    epsilon = 1.0  # kJ/mol

                    if distance > 0.01:  # Avoid division by zero
                        r6 = (sigma / distance) ** 6
                        pe = 4 * epsilon * (r6**2 - r6)
                        total_pe += pe

        return total_pe

    def _build_bond_graph(self, particles: list["Particle"], bonds: list["Bond"]) -> nx.Graph:
        """Build NetworkX graph from particles and bonds."""
        graph = nx.Graph()

        # Add nodes
        for particle in particles:
            graph.add_node(particle.id, particle=particle)

        # Add edges
        for bond in bonds:
            if bond.particle1_id in graph.nodes and bond.particle2_id in graph.nodes:
                graph.add_edge(
                    bond.particle1_id, bond.particle2_id, bond=bond, weight=bond.bond_energy
                )

        return graph

    def _find_clusters(self, graph: nx.Graph, bonds: list["Bond"]) -> list[NetworkCluster]:
        """Find connected components (molecular clusters)."""
        clusters = []

        for component in nx.connected_components(graph):
            if len(component) < 2:
                continue  # Skip isolated particles

            # Get bonds for this cluster
            cluster_bonds = []
            for bond in bonds:
                if bond.particle1_id in component and bond.particle2_id in component:
                    cluster_bonds.append(bond)

            # Calculate complexity metrics
            subgraph = graph.subgraph(component)
            complexity = self._calculate_cluster_complexity(subgraph)
            is_cyclic = len(list(nx.simple_cycles(subgraph))) > 0

            cluster = NetworkCluster(
                particle_ids=set(component),
                bonds=cluster_bonds,
                size=len(component),
                complexity=complexity,
                is_cyclic=is_cyclic,
            )

            clusters.append(cluster)

        return sorted(clusters, key=lambda c: c.size, reverse=True)

    def _calculate_network_metrics(
        self, graph: nx.Graph, clusters: list[NetworkCluster]
    ) -> dict[str, float]:
        """Calculate network-level metrics."""
        metrics = {}

        if len(graph) == 0:
            return metrics

        # Basic metrics
        metrics["density"] = nx.density(graph)
        metrics["average_clustering"] = nx.average_clustering(graph)
        metrics["number_of_clusters"] = len(clusters)

        # Connectivity metrics
        if nx.is_connected(graph):
            metrics["diameter"] = nx.diameter(graph)
            metrics["average_path_length"] = nx.average_shortest_path_length(graph)
        else:
            # For disconnected graphs, analyze largest component
            largest = max(nx.connected_components(graph), key=len)
            if len(largest) > 1:
                subgraph = graph.subgraph(largest)
                metrics["largest_component_diameter"] = nx.diameter(subgraph)
                metrics["largest_component_path_length"] = nx.average_shortest_path_length(subgraph)

        # Degree distribution
        degrees = [d for n, d in graph.degree()]
        if degrees:
            metrics["average_degree"] = np.mean(degrees)
            metrics["degree_variance"] = np.var(degrees)

        return metrics

    def _calculate_cluster_complexity(self, subgraph: nx.Graph) -> float:
        """Calculate structural complexity of a cluster."""
        if len(subgraph) <= 1:
            return 0.0

        # Combine multiple complexity measures
        n_nodes = len(subgraph)
        n_edges = len(subgraph.edges)

        # Topological complexity
        if n_nodes > 1:
            max_edges = n_nodes * (n_nodes - 1) / 2
            edge_density = n_edges / max_edges
        else:
            edge_density = 0.0

        # Clustering coefficient
        clustering = nx.average_clustering(subgraph)

        # Branching factor
        degrees = [d for n, d in subgraph.degree()]
        avg_degree = np.mean(degrees) if degrees else 0.0

        # Combine into single complexity score
        complexity = (
            0.4 * edge_density
            + 0.3 * clustering
            + 0.3 * min(1.0, avg_degree / 4.0)  # Normalize by typical max degree
        )

        return complexity

    def _detect_autocatalytic_sets(
        self, clusters: list[NetworkCluster], bonds: list["Bond"]
    ) -> list[dict[str, any]]:
        """
        Detect autocatalytic sets in the molecular network.

        An autocatalytic set is a collection of molecules where:
        - Each molecule can be produced by reactions within the set
        - The set collectively catalyzes its own formation
        """
        autocatalytic_sets = []

        for cluster in clusters:
            if cluster.size >= 3 and cluster.is_cyclic and cluster.complexity > 0.3:
                # Simple heuristic: cyclic structures with sufficient complexity
                # may represent autocatalytic behavior (threshold for complexity)
                autocatalytic_sets.append(
                    {
                        "particle_ids": list(cluster.particle_ids),
                        "size": cluster.size,
                        "complexity": cluster.complexity,
                        "bond_count": cluster.bond_count,
                        "type": "cyclic_complex",
                    }
                )

        return autocatalytic_sets

    def _detect_replicators(
        self, particles: list["Particle"], clusters: list[NetworkCluster]
    ) -> list[dict[str, any]]:
        """Detect self-replicating structures."""
        replicators = []

        # Group particles by generation and parent
        generations = defaultdict(list)
        for particle in particles:
            if particle.generation > 0:
                generations[particle.generation].append(particle)

        # Look for patterns indicating replication
        for gen, gen_particles in generations.items():
            if gen >= 2:  # Need at least second generation
                # Group by parent
                by_parent = defaultdict(list)
                for p in gen_particles:
                    if p.parent_id:
                        by_parent[p.parent_id].append(p)

                # Look for parents with multiple offspring
                for parent_id, offspring in by_parent.items():
                    if len(offspring) >= 2 and self._are_structurally_similar(offspring):
                        # Check if offspring are structurally similar
                        replicators.append(
                            {
                                "parent_id": parent_id,
                                "offspring_ids": [p.id for p in offspring],
                                "generation": gen,
                                "offspring_count": len(offspring),
                                "type": "generational_replicator",
                            }
                        )

        return replicators

    def _are_structurally_similar(self, particles: list["Particle"]) -> bool:
        """Check if particles have similar bonding patterns."""
        if len(particles) < 2:
            return True

        # Compare bond counts and types
        first = particles[0]
        for particle in particles[1:]:
            if len(particle.bonds) != len(first.bonds):
                return False
            if particle.particle_type != first.particle_type:
                return False

        return True

    def _calculate_complexity_score(self, clusters: list[NetworkCluster]) -> float:
        """Calculate overall system complexity score."""
        if not clusters:
            return 0.0

        # Weight by cluster size and individual complexity
        total_weighted_complexity = 0.0
        total_particles = 0

        for cluster in clusters:
            weight = cluster.size
            total_weighted_complexity += weight * cluster.complexity
            total_particles += cluster.size

        if total_particles == 0:
            return 0.0

        return total_weighted_complexity / total_particles

    def _analyze_energy_trend(self, history: list[EnergyFlow]) -> dict[str, float]:
        """Analyze energy trends over time."""
        if len(history) < 2:
            return {"trend": 0.0}

        # Calculate trends for different energy components
        total_energies = [flow.total_energy for flow in history]
        reaction_energies = [flow.reaction_energy for flow in history]

        # Simple linear trend
        time_points = np.arange(len(total_energies))

        total_trend = (
            np.polyfit(time_points, total_energies, 1)[0] if len(total_energies) > 1 else 0.0
        )
        reaction_trend = (
            np.polyfit(time_points, reaction_energies, 1)[0] if len(reaction_energies) > 1 else 0.0
        )

        return {
            "total_energy_trend": total_trend,
            "reaction_energy_trend": reaction_trend,
            "energy_volatility": np.std(total_energies) if total_energies else 0.0,
        }

    def _calculate_energy_efficiency(self, history: list[EnergyFlow]) -> float:
        """Calculate energy conversion efficiency."""
        if not history:
            return 0.0

        total_input = sum(flow.input_energy for flow in history)
        total_useful = sum(flow.bond_energy + flow.reaction_energy for flow in history)

        if total_input == 0:
            return 0.0

        return total_useful / total_input

    def _track_cluster_growth(self, network_analysis: dict[str, any]) -> float:
        """Track growth rate of largest cluster."""
        # This would need historical data to calculate growth rate
        # For now, return current largest cluster size as proxy
        return network_analysis.get("largest_cluster_size", 0.0)

    def _calculate_shannon_entropy(self, network_analysis: dict[str, any]) -> float:
        """Calculate Shannon entropy of particle type distribution."""
        clusters = network_analysis.get("clusters", [])
        if not clusters:
            return 0.0

        # Count particle types across all clusters
        type_counts = defaultdict(int)
        total_particles = 0

        for cluster in clusters:
            total_particles += cluster.size
            # Would need particle type information for proper calculation
            type_counts["default"] += cluster.size

        if total_particles == 0:
            return 0.0

        # Calculate Shannon entropy
        entropy = 0.0
        for count in type_counts.values():
            if count > 0:
                p = count / total_particles
                entropy -= p * np.log2(p)

        return entropy

    def _calculate_mutual_information(self, network_analysis: dict[str, any]) -> float:
        """Calculate mutual information between structural features."""
        # Simplified implementation
        # Would need more detailed analysis for proper mutual information
        clusters = network_analysis.get("clusters", [])

        if len(clusters) < 2:
            return 0.0

        # Use cluster size distribution as proxy
        sizes = [c.size for c in clusters]
        complexities = [c.complexity for c in clusters]

        # Simple correlation as proxy for mutual information
        if len(sizes) > 1:
            correlation = np.corrcoef(sizes, complexities)[0, 1]
            return abs(correlation) if not np.isnan(correlation) else 0.0

        return 0.0
