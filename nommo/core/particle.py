"""
Particle class representing molecular units with real physical properties.

Each particle represents a molecular unit (~10-100 atoms) with properties
based on real molecular dynamics principles.
"""

from dataclasses import dataclass, field
from typing import List, Optional, Set, Dict, Any
import uuid
import numpy as np

from nommo.utils.logging import get_logger

logger = get_logger("particle")


@dataclass
class Particle:
    """
    A particle represents a molecular unit in the simulation.
    
    Physical units:
    - mass: atomic mass units (amu)
    - position: nanometers (nm)
    - velocity: nm/picosecond (nm/ps)
    - charge: elementary charge units (e)
    - energy: kJ/mol
    - forces: kJ/(mol·nm)
    """
    
    id: str = field(default_factory=lambda: str(uuid.uuid4()))
    
    mass: float = 1.0
    charge: float = 0.0
    radius: float = 0.15
    
    position: np.ndarray = field(default_factory=lambda: np.zeros(3))
    velocity: np.ndarray = field(default_factory=lambda: np.zeros(3))
    forces: np.ndarray = field(default_factory=lambda: np.zeros(3))
    
    particle_type: str = "default"
    type_index: int = 0
    
    bonds: Set[str] = field(default_factory=set)
    max_bonds: int = 4
    
    age: int = 0
    generation: int = 0
    parent_id: Optional[str] = None
    
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    _kinetic_energy_cache: Optional[float] = field(default=None, init=False, repr=False)
    _cache_valid: bool = field(default=False, init=False, repr=False)
    
    def __post_init__(self):
        """Ensure numpy arrays and validate."""
        self.position = np.asarray(self.position, dtype=np.float64)
        self.velocity = np.asarray(self.velocity, dtype=np.float64)
        self.forces = np.asarray(self.forces, dtype=np.float64)
        
        if self.position.shape != (3,):
            raise ValueError(f"Position must be 3D vector, got shape {self.position.shape}")
        if self.velocity.shape != (3,):
            raise ValueError(f"Velocity must be 3D vector, got shape {self.velocity.shape}")
        if self.forces.shape != (3,):
            raise ValueError(f"Forces must be 3D vector, got shape {self.forces.shape}")
    
    def kinetic_energy(self) -> float:
        """
        Calculate kinetic energy: E_k = (1/2)mv²
        
        Returns:
            Kinetic energy in kJ/mol
        """
        if not self._cache_valid:
            v_squared = np.dot(self.velocity, self.velocity)
            self._kinetic_energy_cache = 0.5 * self.mass * v_squared * 0.01
            self._cache_valid = True
        return self._kinetic_energy_cache
    
    def momentum(self) -> np.ndarray:
        """
        Calculate momentum: p = mv
        
        Returns:
            Momentum vector in amu·nm/ps
        """
        return self.mass * self.velocity
    
    def distance_to(self, other: "Particle") -> float:
        """
        Calculate Euclidean distance to another particle.
        
        Args:
            other: Target particle
            
        Returns:
            Distance in nm
        """
        return np.linalg.norm(self.position - other.position)
    
    def distance_squared_to(self, other: "Particle") -> float:
        """
        Calculate squared distance (avoids sqrt for performance).
        
        Args:
            other: Target particle
            
        Returns:
            Squared distance in nm²
        """
        delta = self.position - other.position
        return np.dot(delta, delta)
    
    def distance_vector_to(self, other: "Particle") -> np.ndarray:
        """
        Calculate vector from this particle to another.
        
        Args:
            other: Target particle
            
        Returns:
            Vector pointing from self to other
        """
        return other.position - self.position
    
    def can_bond_with(self, other: "Particle") -> bool:
        """
        Check if this particle can form a bond with another.
        
        Args:
            other: Potential bonding partner
            
        Returns:
            True if bonding is possible
        """
        return (
            len(self.bonds) < self.max_bonds and
            len(other.bonds) < other.max_bonds and
            other.id not in self.bonds and
            self.id != other.id
        )
    
    def add_bond(self, other: "Particle") -> bool:
        """
        Form a bond with another particle.
        
        Args:
            other: Particle to bond with
            
        Returns:
            True if bond was formed, False otherwise
        """
        if not self.can_bond_with(other):
            return False
        
        self.bonds.add(other.id)
        other.bonds.add(self.id)
        return True
    
    def remove_bond(self, other_id: str) -> bool:
        """
        Break a bond with another particle.
        
        Args:
            other_id: ID of bonded particle
            
        Returns:
            True if bond was broken, False if no bond existed
        """
        if other_id in self.bonds:
            self.bonds.remove(other_id)
            return True
        return False
    
    def clear_bonds(self) -> int:
        """
        Remove all bonds.
        
        Returns:
            Number of bonds removed
        """
        count = len(self.bonds)
        self.bonds.clear()
        return count
    
    def invalidate_cache(self):
        """Invalidate cached values when velocity changes."""
        self._cache_valid = False
    
    def update_velocity(self, new_velocity: np.ndarray):
        """Update velocity and invalidate cache."""
        self.velocity = np.asarray(new_velocity, dtype=np.float64)
        self.invalidate_cache()
    
    def update_position(self, new_position: np.ndarray):
        """Update position."""
        self.position = np.asarray(new_position, dtype=np.float64)
    
    def apply_periodic_boundary(self, box_size: np.ndarray):
        """
        Apply periodic boundary conditions.
        
        Args:
            box_size: Simulation box dimensions
        """
        self.position = self.position % box_size
    
    def apply_reflective_boundary(self, box_size: np.ndarray):
        """
        Apply reflective boundary conditions.
        
        Args:
            box_size: Simulation box dimensions
        """
        for i in range(3):
            if self.position[i] < 0:
                self.position[i] = -self.position[i]
                self.velocity[i] = -self.velocity[i]
                self.invalidate_cache()
            elif self.position[i] > box_size[i]:
                self.position[i] = 2 * box_size[i] - self.position[i]
                self.velocity[i] = -self.velocity[i]
                self.invalidate_cache()
    
    def copy(self, new_id: bool = True) -> "Particle":
        """
        Create a copy of this particle.
        
        Args:
            new_id: Whether to generate a new ID
            
        Returns:
            Copied particle
        """
        p = Particle(
            id=str(uuid.uuid4()) if new_id else self.id,
            mass=self.mass,
            charge=self.charge,
            radius=self.radius,
            position=self.position.copy(),
            velocity=self.velocity.copy(),
            forces=self.forces.copy(),
            particle_type=self.particle_type,
            type_index=self.type_index,
            bonds=self.bonds.copy(),
            max_bonds=self.max_bonds,
            age=0,
            generation=self.generation + 1,
            parent_id=self.id,
            metadata=self.metadata.copy()
        )
        return p
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary for serialization."""
        return {
            "id": self.id,
            "mass": self.mass,
            "charge": self.charge,
            "radius": self.radius,
            "position": self.position.tolist(),
            "velocity": self.velocity.tolist(),
            "forces": self.forces.tolist(),
            "particle_type": self.particle_type,
            "type_index": self.type_index,
            "bonds": list(self.bonds),
            "max_bonds": self.max_bonds,
            "age": self.age,
            "generation": self.generation,
            "parent_id": self.parent_id,
            "metadata": self.metadata
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "Particle":
        """Create particle from dictionary."""
        return cls(
            id=data["id"],
            mass=data["mass"],
            charge=data["charge"],
            radius=data["radius"],
            position=np.array(data["position"]),
            velocity=np.array(data["velocity"]),
            forces=np.array(data["forces"]),
            particle_type=data["particle_type"],
            type_index=data.get("type_index", 0),
            bonds=set(data["bonds"]),
            max_bonds=data["max_bonds"],
            age=data["age"],
            generation=data.get("generation", 0),
            parent_id=data.get("parent_id"),
            metadata=data.get("metadata", {})
        )


class ParticlePool:
    """
    Object pool for efficient particle memory management.
    
    Reuses particle objects to reduce allocation overhead.
    """
    
    def __init__(self, initial_size: int = 100):
        """
        Initialize particle pool.
        
        Args:
            initial_size: Initial pool size
        """
        self._pool: List[Particle] = []
        self._active: Set[str] = set()
        
        for _ in range(initial_size):
            self._pool.append(Particle())
            
        logger.debug(f"Initialized particle pool with {initial_size} particles")
    
    def acquire(self, **kwargs) -> Particle:
        """
        Get a particle from the pool.
        
        Args:
            **kwargs: Particle initialization parameters
            
        Returns:
            Configured particle
        """
        if self._pool:
            particle = self._pool.pop()
        else:
            particle = Particle()
            
        particle.id = kwargs.get("id", str(uuid.uuid4()))
        particle.mass = kwargs.get("mass", 1.0)
        particle.charge = kwargs.get("charge", 0.0)
        particle.radius = kwargs.get("radius", 0.15)
        particle.position = np.asarray(kwargs.get("position", np.zeros(3)))
        particle.velocity = np.asarray(kwargs.get("velocity", np.zeros(3)))
        particle.forces = np.zeros(3)
        particle.particle_type = kwargs.get("particle_type", "default")
        particle.type_index = kwargs.get("type_index", 0)
        particle.bonds = set()
        particle.max_bonds = kwargs.get("max_bonds", 4)
        particle.age = 0
        particle.generation = kwargs.get("generation", 0)
        particle.parent_id = kwargs.get("parent_id")
        particle.metadata = kwargs.get("metadata", {})
        particle.invalidate_cache()
        
        self._active.add(particle.id)
        return particle
    
    def release(self, particle: Particle):
        """
        Return a particle to the pool.
        
        Args:
            particle: Particle to release
        """
        if particle.id in self._active:
            self._active.remove(particle.id)
            particle.bonds.clear()
            particle.metadata.clear()
            self._pool.append(particle)
    
    def release_all(self, particles: List[Particle]):
        """Release multiple particles."""
        for p in particles:
            self.release(p)
    
    @property
    def pool_size(self) -> int:
        """Current pool size."""
        return len(self._pool)
    
    @property
    def active_count(self) -> int:
        """Number of active particles."""
        return len(self._active)