"""
Main simulation engine for the Nommo universe.

Orchestrates particles, physics, chemistry, and time evolution.
"""

from dataclasses import dataclass, field
from typing import List, Dict, Optional, Callable, Any, Tuple
from pathlib import Path
import uuid
import time
import numpy as np
from collections import defaultdict, deque

from nommo.core.particle import Particle, ParticlePool
from nommo.core.types import (
    UniverseParameters, BoundaryType, EnergySource,
    ThermostatType, IntegratorType
)
from nommo.core.spatial import SpatialIndex
from nommo.physics.forces import LennardJonesForce, CompositeForce
from nommo.physics.integrator import (
    VelocityVerletIntegrator, AdaptiveTimestepIntegrator
)
from nommo.physics.thermostat import create_thermostat, calculate_temperature
from nommo.physics.constants import KB_KJMOL_K, velocity_from_temperature
from nommo.utils.logging import get_logger, get_simulation_logger, get_performance_logger

logger = get_logger("universe")
sim_logger = get_simulation_logger()
perf_logger = get_performance_logger()


@dataclass
class UniverseMetrics:
    """Real-time metrics for tracking emergence."""
    
    tick: int = 0
    time: float = 0.0
    
    total_particles: int = 0
    particles_by_type: Dict[str, int] = field(default_factory=dict)
    
    total_bonds: int = 0
    avg_bonds_per_particle: float = 0.0
    max_cluster_size: int = 0
    shannon_entropy: float = 0.0
    
    total_kinetic_energy: float = 0.0
    total_potential_energy: float = 0.0
    total_energy: float = 0.0
    current_temperature: float = 0.0
    
    reactions_this_tick: int = 0
    bonds_formed_this_tick: int = 0
    bonds_broken_this_tick: int = 0
    
    replicator_count: int = 0
    generation_max: int = 0
    
    is_active: bool = True
    has_replicators: bool = False
    complexity_trend: str = "stable"
    
    computation_time: float = 0.0
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "tick": self.tick,
            "time": self.time,
            "total_particles": self.total_particles,
            "particles_by_type": self.particles_by_type,
            "total_bonds": self.total_bonds,
            "avg_bonds_per_particle": self.avg_bonds_per_particle,
            "max_cluster_size": self.max_cluster_size,
            "shannon_entropy": self.shannon_entropy,
            "total_kinetic_energy": self.total_kinetic_energy,
            "total_potential_energy": self.total_potential_energy,
            "total_energy": self.total_energy,
            "current_temperature": self.current_temperature,
            "reactions_this_tick": self.reactions_this_tick,
            "bonds_formed_this_tick": self.bonds_formed_this_tick,
            "bonds_broken_this_tick": self.bonds_broken_this_tick,
            "replicator_count": self.replicator_count,
            "generation_max": self.generation_max,
            "is_active": self.is_active,
            "has_replicators": self.has_replicators,
            "complexity_trend": self.complexity_trend,
            "computation_time": self.computation_time
        }


class Universe:
    """
    Main simulation engine managing the entire universe.
    
    Coordinates particles, physics, chemistry, and evolution.
    """
    
    def __init__(self, parameters: UniverseParameters):
        """
        Initialize universe with given parameters.
        
        Args:
            parameters: Complete universe configuration
        """
        self.id = str(uuid.uuid4())
        self.params = parameters
        self.tick = 0
        self.time = 0.0
        
        self.particles: List[Particle] = []
        self.particle_dict: Dict[str, Particle] = {}
        self.particle_pool = ParticlePool(initial_size=parameters.max_particles // 10)
        
        self.box_size = np.array(parameters.box.dimensions)
        
        self._setup_physics()
        self._setup_chemistry()
        self._initialize_particles()
        
        self.metrics = UniverseMetrics()
        self.history: deque[UniverseMetrics] = deque(
            maxlen=parameters.history_buffer_size
        )
        
        self._energy_accumulator = 0.0
        self._last_checkpoint = 0
        
        logger.info(
            f"Created universe '{parameters.name}' with {len(self.particles)} particles"
        )
    
    def _setup_physics(self):
        """Initialize physics components."""
        params = self.params.physics
        
        self.spatial_index = SpatialIndex(
            self.box_size,
            params.cutoff_distance,
            params.neighbor_skin,
            use_cell_list=True
        )
        
        self.force_field = LennardJonesForce(
            epsilon=params.lj_epsilon,
            sigma=params.lj_sigma,
            cutoff=params.cutoff_distance,
            shift=True
        )
        
        if self.params.particle_types:
            n_types = len(self.params.particle_types)
            epsilon_values = np.zeros(n_types)
            sigma_values = np.zeros(n_types)
            
            for i, ptype in enumerate(self.params.particle_types.values()):
                epsilon_values[i] = ptype.lj_epsilon or params.lj_epsilon
                sigma_values[i] = ptype.lj_sigma or params.lj_sigma
                
            self.force_field.setup_parameters(n_types, epsilon_values, sigma_values)
        
        if params.adaptive_timestep:
            self.integrator = AdaptiveTimestepIntegrator(
                params.timestep,
                self.box_size,
                self.params.box.boundary_type
            )
        else:
            self.integrator = VelocityVerletIntegrator(
                params.timestep,
                self.box_size,
                self.params.box.boundary_type
            )
        
        thermo_params = self.params.thermodynamics
        self.thermostat = create_thermostat(
            thermo_params.thermostat.value,
            thermo_params.temperature,
            coupling_time=thermo_params.thermostat_coupling
        )
    
    def _setup_chemistry(self):
        """Initialize chemistry components."""
        self.reaction_rules = self.params.chemistry.reaction_rules
        self.bond_energy = self.params.chemistry.bond_energy
        self.bonding_distance = self.params.chemistry.bonding_distance
    
    def _initialize_particles(self):
        """Create initial particle population."""
        particle_types = list(self.params.particle_types.keys())
        type_configs = self.params.particle_types
        
        particle_index = 0
        for ptype, count in self.params.initial_composition.items():
            if ptype not in particle_types:
                logger.warning(f"Unknown particle type: {ptype}")
                continue
                
            config = type_configs[ptype]
            type_idx = particle_types.index(ptype)
            
            for _ in range(count):
                position = np.random.random(3) * self.box_size
                
                sigma_v = velocity_from_temperature(
                    self.params.thermodynamics.temperature,
                    config.mass
                )
                velocity = np.random.normal(0, sigma_v, 3)
                
                particle = self.particle_pool.acquire(
                    mass=config.mass,
                    charge=config.charge,
                    radius=config.radius,
                    position=position,
                    velocity=velocity,
                    particle_type=ptype,
                    type_index=type_idx,
                    max_bonds=config.max_bonds
                )
                
                self.particles.append(particle)
                self.particle_dict[particle.id] = particle
                particle_index += 1
        
        velocities = np.array([p.velocity for p in self.particles])
        velocities -= np.mean(velocities, axis=0)
        for i, p in enumerate(self.particles):
            p.update_velocity(velocities[i])
    
    @perf_logger.timer("simulation_step")
    def step(self) -> UniverseMetrics:
        """
        Advance simulation by one timestep.
        
        Returns:
            Current metrics after step
        """
        step_start = time.perf_counter()
        
        positions = np.array([p.position for p in self.particles])
        velocities = np.array([p.velocity for p in self.particles])
        forces = np.array([p.forces for p in self.particles])
        masses = np.array([p.mass for p in self.particles])
        types = np.array([p.type_index for p in self.particles])
        
        self.spatial_index.update(positions, force_rebuild=(self.tick % 10 == 0))
        neighbor_pairs = self.spatial_index.get_pairs()
        
        def force_calculator(pos):
            return self.force_field.calculate(pos, types, self.box_size, neighbor_pairs)
        
        potential_energy = self.integrator.step(
            positions, velocities, forces, masses,
            force_calculator, self.params.physics.timestep
        )
        
        for i, p in enumerate(self.particles):
            p.position = positions[i]
            p.velocity = velocities[i]
            p.forces = forces[i]
            p.age += 1
        
        if self.thermostat:
            current_temp = self.thermostat.apply(
                velocities, masses, self.params.physics.timestep
            )
            for i, p in enumerate(self.particles):
                p.velocity = velocities[i]
        else:
            current_temp = calculate_temperature(velocities, masses)
        
        if self.params.chemistry.enable_reactions:
            self._process_reactions(current_temp)
        
        if self.params.chemistry.enable_bonding:
            self._process_bonds(current_temp)
        
        self._apply_energy_source()
        
        self.tick += 1
        self.time += self.params.physics.timestep
        
        metrics = self._calculate_metrics(potential_energy, current_temp)
        metrics.computation_time = time.perf_counter() - step_start
        
        self.metrics = metrics
        self.history.append(metrics)
        
        if self.tick % 100 == 0:
            self._log_metrics(metrics)
        
        if self.tick - self._last_checkpoint >= self.params.checkpoint_interval:
            self._checkpoint()
        
        return metrics
    
    def _process_reactions(self, temperature: float):
        """Process chemical reactions."""
        pass
    
    def _process_bonds(self, temperature: float):
        """Process bond formation and breaking."""
        bonds_formed = 0
        bonds_broken = 0
        
        for pair in self.spatial_index.get_pairs():
            p1, p2 = self.particles[pair[0]], self.particles[pair[1]]
            distance = p1.distance_to(p2)
            
            if distance < self.bonding_distance:
                if p1.can_bond_with(p2):
                    if np.random.random() < self.params.chemistry.reaction_probability:
                        if p1.add_bond(p2):
                            bonds_formed += 1
            
            if p2.id in p1.bonds:
                break_prob = np.exp(-self.bond_energy / (KB_KJMOL_K * temperature))
                if np.random.random() < break_prob * self.params.physics.timestep:
                    p1.remove_bond(p2.id)
                    p2.remove_bond(p1.id)
                    bonds_broken += 1
        
        self.metrics.bonds_formed_this_tick = bonds_formed
        self.metrics.bonds_broken_this_tick = bonds_broken
    
    def _apply_energy_source(self):
        """Apply external energy input."""
        if self.params.thermodynamics.energy_amount <= 0:
            return
            
        energy_type = self.params.thermodynamics.energy_source
        
        if energy_type == EnergySource.CONSTANT:
            pass
        elif energy_type == EnergySource.PULSES:
            if self.tick % self.params.thermodynamics.energy_interval == 0:
                energy_per_particle = (
                    self.params.thermodynamics.energy_amount / len(self.particles)
                )
                for p in self.particles:
                    speed = np.linalg.norm(p.velocity)
                    if speed > 0:
                        p.velocity *= 1 + energy_per_particle / (p.mass * speed**2)
                        p.invalidate_cache()
    
    def _calculate_metrics(
        self,
        potential_energy: float,
        temperature: float
    ) -> UniverseMetrics:
        """Calculate current universe metrics."""
        metrics = UniverseMetrics()
        
        metrics.tick = self.tick
        metrics.time = self.time
        metrics.total_particles = len(self.particles)
        
        type_counts = defaultdict(int)
        for p in self.particles:
            type_counts[p.particle_type] += 1
        metrics.particles_by_type = dict(type_counts)
        
        total_bonds = sum(len(p.bonds) for p in self.particles)
        metrics.total_bonds = total_bonds // 2
        metrics.avg_bonds_per_particle = (
            total_bonds / len(self.particles) if self.particles else 0
        )
        
        kinetic_energy = sum(p.kinetic_energy() for p in self.particles)
        metrics.total_kinetic_energy = kinetic_energy
        metrics.total_potential_energy = potential_energy
        metrics.total_energy = kinetic_energy + potential_energy
        metrics.current_temperature = temperature
        
        if len(type_counts) > 0:
            total = metrics.total_particles
            entropy = 0.0
            for count in type_counts.values():
                p = count / total
                if p > 0:
                    entropy -= p * np.log2(p)
            metrics.shannon_entropy = entropy
        
        generations = [p.generation for p in self.particles]
        metrics.generation_max = max(generations) if generations else 0
        
        if len(self.history) >= 10:
            recent_complexity = [h.shannon_entropy for h in list(self.history)[-10:]]
            if recent_complexity[-1] > recent_complexity[0] * 1.1:
                metrics.complexity_trend = "increasing"
            elif recent_complexity[-1] < recent_complexity[0] * 0.9:
                metrics.complexity_trend = "decreasing"
        
        return metrics
    
    def _log_metrics(self, metrics: UniverseMetrics):
        """Log metrics periodically."""
        sim_logger.log_metrics(
            self.id,
            self.tick,
            {
                "particles": metrics.total_particles,
                "temperature": metrics.current_temperature,
                "energy": metrics.total_energy,
                "entropy": metrics.shannon_entropy,
                "bonds": metrics.total_bonds
            }
        )
    
    def _checkpoint(self):
        """Save checkpoint."""
        self._last_checkpoint = self.tick
    
    def run(
        self,
        ticks: int,
        callback: Optional[Callable[[UniverseMetrics], None]] = None,
        stop_condition: Optional[Callable[[UniverseMetrics], bool]] = None
    ) -> List[UniverseMetrics]:
        """
        Run simulation for multiple ticks.
        
        Args:
            ticks: Number of ticks to simulate
            callback: Optional callback after each tick
            stop_condition: Optional early stop condition
            
        Returns:
            List of metrics for each tick
        """
        metrics_list = []
        
        logger.info(f"Starting simulation for {ticks} ticks")
        start_time = time.perf_counter()
        
        for i in range(ticks):
            metrics = self.step()
            metrics_list.append(metrics)
            
            if callback:
                callback(metrics)
            
            if stop_condition and stop_condition(metrics):
                logger.info(f"Stop condition met at tick {self.tick}")
                break
            
            if (i + 1) % 1000 == 0:
                elapsed = time.perf_counter() - start_time
                ticks_per_second = (i + 1) / elapsed
                logger.info(
                    f"Progress: {i+1}/{ticks} ticks, "
                    f"{ticks_per_second:.1f} ticks/s"
                )
        
        elapsed = time.perf_counter() - start_time
        logger.info(
            f"Simulation complete: {len(metrics_list)} ticks in {elapsed:.2f}s"
        )
        
        return metrics_list
    
    def add_particle(self, **kwargs) -> Particle:
        """Add a new particle to the universe."""
        if len(self.particles) >= self.params.max_particles:
            logger.warning("Maximum particle count reached")
            return None
            
        particle = self.particle_pool.acquire(**kwargs)
        self.particles.append(particle)
        self.particle_dict[particle.id] = particle
        return particle
    
    def remove_particle(self, particle_id: str) -> bool:
        """Remove a particle from the universe."""
        if particle_id not in self.particle_dict:
            return False
            
        particle = self.particle_dict.pop(particle_id)
        self.particles.remove(particle)
        
        for other_id in list(particle.bonds):
            if other_id in self.particle_dict:
                self.particle_dict[other_id].remove_bond(particle_id)
        
        self.particle_pool.release(particle)
        return True
    
    def get_state(self) -> Dict[str, Any]:
        """Get current universe state."""
        return {
            "id": self.id,
            "tick": self.tick,
            "time": self.time,
            "parameters": self.params.model_dump(),
            "particles": [p.to_dict() for p in self.particles],
            "metrics": self.metrics.to_dict()
        }
    
    def load_state(self, state: Dict[str, Any]):
        """Load universe state."""
        self.id = state["id"]
        self.tick = state["tick"]
        self.time = state["time"]
        
        self.particles.clear()
        self.particle_dict.clear()
        
        for p_data in state["particles"]:
            particle = Particle.from_dict(p_data)
            self.particles.append(particle)
            self.particle_dict[particle.id] = particle