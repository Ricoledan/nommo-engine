# **Nommo Engine - Complete Specification for Claude Code**

## **Project Overview**

**Name:** Nommo Engine  
**Type:** Scientific life emergence simulator with CLI interface  
**Language:** Python 3.11+  
**Purpose:** Explore how self-replicating systems emerge from non-living chemistry using real physics and thermodynamics

**Scientific Foundation:**
- Molecular dynamics (Lennard-Jones potentials, Velocity Verlet integration)
- Chemical kinetics (Arrhenius equation, collision theory)
- Thermodynamics (Gibbs free energy, Boltzmann distributions)
- Complexity theory (autocatalytic sets, emergence detection)

**Key References:**
- Kauffman (1986) - Autocatalytic sets
- Prigogine & Nicolis (1977) - Self-organization in nonequilibrium systems
- Frenkel & Smit (2001) - Molecular simulation methods
- Ruiz-Mirazo et al. (2014) - Prebiotic systems chemistry

---

## **1. Project Structure**

```
nommo-engine/
├── pyproject.toml              # Poetry configuration
├── README.md
├── docs/
│   ├── SCIENTIFIC_BASIS.md     # Theory and equations
│   ├── CLI_GUIDE.md            # User documentation
│   └── DEVELOPER.md            # Code documentation
├── nommo/
│   ├── __init__.py
│   ├── cli/
│   │   ├── __init__.py
│   │   ├── main.py             # Click CLI entry point
│   │   ├── universe.py         # Universe commands
│   │   ├── run.py              # Simulation commands
│   │   ├── analyze.py          # Analysis commands
│   │   └── utils.py            # CLI utilities (prompts, tables)
│   ├── core/
│   │   ├── __init__.py
│   │   ├── particle.py         # Particle class with physics
│   │   ├── universe.py         # Universe simulation container
│   │   ├── spatial.py          # Spatial indexing (quadtree/grid)
│   │   └── types.py            # Type definitions (enums, dataclasses)
│   ├── physics/
│   │   ├── __init__.py
│   │   ├── forces.py           # Lennard-Jones, Coulombic forces
│   │   ├── integrator.py       # Velocity Verlet integrator
│   │   ├── thermostat.py       # Temperature control (Berendsen)
│   │   └── constants.py        # Physical constants
│   ├── chemistry/
│   │   ├── __init__.py
│   │   ├── reactions.py        # Reaction rules and detection
│   │   ├── kinetics.py         # Arrhenius, activation energy
│   │   ├── bonds.py            # Bond formation/breaking
│   │   └── thermodynamics.py   # Gibbs free energy calculations
│   ├── analysis/
│   │   ├── __init__.py
│   │   ├── metrics.py          # Complexity, diversity, entropy
│   │   ├── detection.py        # Detect autocatalysis, replicators
│   │   └── visualization.py    # Matplotlib plotting
│   ├── storage/
│   │   ├── __init__.py
│   │   ├── config.py           # Configuration management
│   │   ├── persistence.py      # Save/load universe states
│   │   └── export.py           # Export data (CSV, JSON)
│   └── presets/
│       ├── __init__.py
│       ├── earth_like.py       # Earth-like parameters
│       ├── high_energy.py      # High energy environment
│       └── minimal.py          # Minimal complexity test
├── tests/
│   ├── __init__.py
│   ├── test_particle.py
│   ├── test_forces.py
│   ├── test_reactions.py
│   └── test_universe.py
└── data/                       # Created at runtime
    ├── universes/              # Saved universe states
    ├── exports/                # Exported data/plots
    └── config.json             # User configuration
```

---

## **2. Core Data Structures**

### **2.1 Particle Class**

```python
"""
nommo/core/particle.py

Represents a molecular particle with real physical properties.
"""

from dataclasses import dataclass
from typing import List, Optional
import numpy as np

@dataclass
class Particle:
    """
    A particle represents a molecular unit (~10-100 atoms).
    
    Physical units:
    - mass: atomic mass units (amu)
    - position: nanometers (nm)
    - velocity: nm/picosecond (nm/ps)
    - charge: elementary charge units (e)
    - energy: kJ/mol
    """
    
    # Unique identifier
    id: str
    
    # Physical properties
    mass: float              # amu
    charge: float            # elementary charges
    radius: float            # nm (for collision detection)
    
    # State
    position: np.ndarray     # [x, y, z] in nm
    velocity: np.ndarray     # [vx, vy, vz] in nm/ps
    forces: np.ndarray       # [fx, fy, fz] in kJ/(mol·nm)
    
    # Chemistry
    particle_type: str       # e.g., "monomer_A", "monomer_B"
    bonds: List[str]         # List of particle IDs bonded to
    max_bonds: int           # Valency
    
    # Bookkeeping
    age: int = 0             # Simulation ticks since creation
    parent_id: Optional[str] = None  # For tracking lineage
    
    def kinetic_energy(self) -> float:
        """Calculate kinetic energy: E_k = (1/2)mv²"""
        v_squared = np.dot(self.velocity, self.velocity)
        # Convert to kJ/mol: (amu * (nm/ps)²) * conversion_factor
        return 0.5 * self.mass * v_squared * 0.01  # Conversion factor
    
    def distance_to(self, other: 'Particle') -> float:
        """Calculate distance to another particle"""
        return np.linalg.norm(self.position - other.position)
    
    def can_bond_with(self, other: 'Particle') -> bool:
        """Check if this particle can form a bond"""
        return (
            len(self.bonds) < self.max_bonds and
            len(other.bonds) < other.max_bonds and
            other.id not in self.bonds
        )
```

### **2.2 Universe Parameters**

```python
"""
nommo/core/types.py

Type definitions and parameter structures.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, Tuple

class BoundaryType(Enum):
    """Boundary condition types"""
    PERIODIC = "periodic"     # Wrap around (infinite-like)
    REFLECTIVE = "reflective" # Bounce off walls
    ABSORBING = "absorbing"   # Particles disappear at edges

class EnergySource(Enum):
    """Energy input types"""
    CONSTANT = "constant"     # Steady heat bath
    PULSES = "pulses"         # Periodic bursts
    GRADIENT = "gradient"     # Spatial variation (hot/cold regions)
    CHAOTIC = "chaotic"       # Random spikes

@dataclass
class UniverseParameters:
    """
    Complete parameter set for a universe.
    Based on real physical chemistry principles.
    """
    
    # Identification
    name: str
    description: str
    
    # Physics
    temperature: float              # Kelvin (100-1000K typical)
    mass_density: float             # particles/nm³
    dimensions: Tuple[float, float, float]  # (Lx, Ly, Lz) in nm
    boundary_type: BoundaryType
    
    # Time integration
    timestep: float                 # Picoseconds (0.001-0.01 typical)
    
    # Force field (Lennard-Jones parameters)
    lj_epsilon: float               # Potential well depth (kJ/mol)
    lj_sigma: float                 # Zero-crossing distance (nm)
    cutoff_distance: float          # Force cutoff (nm)
    
    # Chemistry
    activation_energy: float        # kJ/mol (barrier to reactions)
    bond_energy: float              # kJ/mol (strength of bonds)
    reaction_probability: float     # Base probability (0-1)
    
    # Energy input
    energy_source: EnergySource
    energy_amount: float            # kJ/mol added per pulse
    energy_interval: int            # Ticks between pulses (if applicable)
    
    # Particle types
    particle_types: Dict[str, Dict] # Type definitions
    initial_composition: Dict[str, int]  # Initial particle counts by type
    
    # Limits
    max_particles: int              # Prevent runaway growth
    
    # Thermostat
    use_thermostat: bool
    thermostat_coupling: float      # Time constant (ps)
```

### **2.3 Universe State**

```python
"""
nommo/core/universe.py

Main simulation container.
"""

from dataclasses import dataclass
from typing import List, Dict
import numpy as np
from .particle import Particle
from .types import UniverseParameters

@dataclass
class UniverseMetrics:
    """Real-time metrics for tracking emergence"""
    
    # Population
    total_particles: int
    particles_by_type: Dict[str, int]
    
    # Complexity
    total_bonds: int
    avg_bonds_per_particle: float
    max_cluster_size: int           # Largest connected component
    shannon_entropy: float          # Diversity measure
    
    # Thermodynamics
    total_kinetic_energy: float     # kJ/mol
    total_potential_energy: float   # kJ/mol
    current_temperature: float      # K
    
    # Chemistry
    reactions_this_tick: int
    bonds_formed_this_tick: int
    bonds_broken_this_tick: int
    
    # Evolution markers
    replicator_count: int
    generation_max: int
    
    # State assessment
    is_active: bool                 # Recent chemical activity
    has_replicators: bool
    complexity_trend: str           # "increasing", "stable", "decreasing"

class Universe:
    """
    Main simulation engine.
    Manages particles, physics, chemistry, and time evolution.
    """
    
    def __init__(self, parameters: UniverseParameters):
        self.params = parameters
        self.particles: List[Particle] = []
        self.tick: int = 0
        self.metrics: UniverseMetrics = None
        self.history: List[UniverseMetrics] = []
        
        # Spatial indexing for performance
        self.spatial_index = None  # Will be Grid or QuadTree
        
        # Initialize universe
        self._initialize_particles()
        self._initialize_spatial_index()
    
    def step(self) -> UniverseMetrics:
        """
        Advance simulation by one timestep.
        
        1. Calculate forces (physics)
        2. Integrate motion (Verlet)
        3. Apply thermostat (if enabled)
        4. Detect/apply reactions (chemistry)
        5. Handle bonds (form/break)
        6. Update metrics
        7. Check for replication
        """
        pass
    
    def run(self, ticks: int, callback=None) -> List[UniverseMetrics]:
        """Run for N ticks, optionally calling callback each step"""
        pass
```

---

## **3. Physics Engine**

### **3.1 Force Calculations**

```python
"""
nommo/physics/forces.py

Real molecular forces: Lennard-Jones potential.
"""

import numpy as np
from typing import Tuple

class LennardJonesForce:
    """
    Lennard-Jones 12-6 potential:
    V(r) = 4ε[(σ/r)¹² - (σ/r)⁶]
    
    F(r) = -dV/dr = 24ε/r [(2(σ/r)¹² - (σ/r)⁶)]
    
    Parameters from Frenkel & Smit (2001)
    """
    
    def __init__(self, epsilon: float, sigma: float, cutoff: float):
        """
        Args:
            epsilon: Potential well depth (kJ/mol)
            sigma: Zero-crossing distance (nm)
            cutoff: Force cutoff distance (nm) - typically 2.5*sigma
        """
        self.epsilon = epsilon
        self.sigma = sigma
        self.cutoff = cutoff
        self.cutoff_squared = cutoff * cutoff
    
    def calculate(
        self, 
        r_vec: np.ndarray, 
        r: float
    ) -> Tuple[np.ndarray, float]:
        """
        Calculate force and potential energy.
        
        Args:
            r_vec: Vector from particle 1 to particle 2
            r: Distance between particles
            
        Returns:
            (force_vector, potential_energy)
        """
        if r > self.cutoff or r < 1e-10:
            return np.zeros(3), 0.0
        
        # Calculate terms
        sigma_r = self.sigma / r
        sigma_r6 = sigma_r**6
        sigma_r12 = sigma_r6**2
        
        # Force magnitude
        force_mag = 24 * self.epsilon / r * (2 * sigma_r12 - sigma_r6)
        
        # Force vector (pointing from 1 to 2)
        force_vec = force_mag * (r_vec / r)
        
        # Potential energy
        potential = 4 * self.epsilon * (sigma_r12 - sigma_r6)
        
        return force_vec, potential

def calculate_all_forces(
    particles: List[Particle],
    force_calculator: LennardJonesForce,
    spatial_index
) -> Tuple[np.ndarray, float]:
    """
    Calculate forces on all particles efficiently using spatial indexing.
    
    Returns:
        (force_array, total_potential_energy)
    """
    pass
```

### **3.2 Time Integration**

```python
"""
nommo/physics/integrator.py

Velocity Verlet integration - symplectic, time-reversible, energy-conserving.
From Frenkel & Smit (2001), Chapter 4.
"""

import numpy as np
from typing import List
from nommo.core.particle import Particle

class VelocityVerletIntegrator:
    """
    Velocity Verlet algorithm:
    
    1. x(t+Δt) = x(t) + v(t)Δt + ½a(t)Δt²
    2. Calculate forces at new positions
    3. v(t+Δt) = v(t) + ½[a(t) + a(t+Δt)]Δt
    
    This is more accurate than simple Euler and conserves energy better.
    """
    
    def __init__(self, timestep: float):
        """
        Args:
            timestep: Integration timestep in picoseconds
        """
        self.dt = timestep
        self.dt_sq = timestep * timestep
    
    def step(
        self, 
        particles: List[Particle],
        force_calculator,
        boundary_handler
    ):
        """
        Advance all particles by one timestep.
        
        Args:
            particles: List of particles to integrate
            force_calculator: Function to calculate forces
            boundary_handler: Handle periodic/reflective boundaries
        """
        # Step 1: Update positions
        for p in particles:
            acceleration = p.forces / p.mass
            p.position += p.velocity * self.dt + 0.5 * acceleration * self.dt_sq
            boundary_handler.apply(p)
        
        # Step 2: Calculate new forces
        new_forces, potential_energy = force_calculator(particles)
        
        # Step 3: Update velocities
        for i, p in enumerate(particles):
            old_acceleration = p.forces / p.mass
            new_acceleration = new_forces[i] / p.mass
            p.velocity += 0.5 * (old_acceleration + new_acceleration) * self.dt
            p.forces = new_forces[i]
        
        return potential_energy
```

### **3.3 Temperature Control**

```python
"""
nommo/physics/thermostat.py

Berendsen thermostat for temperature control.
Reference: Berendsen et al. (1984) J. Chem. Phys. 81, 3684
"""

import numpy as np
from typing import List
from nommo.core.particle import Particle

class BerendsenThermostat:
    """
    Rescale velocities to maintain target temperature.
    
    λ = sqrt(1 + (Δt/τ)(T_target/T_current - 1))
    v_new = λ * v_old
    
    Where τ is coupling time constant.
    """
    
    def __init__(self, target_temperature: float, coupling_time: float):
        """
        Args:
            target_temperature: Target temperature (K)
            coupling_time: Time constant for coupling (ps)
        """
        self.T_target = target_temperature
        self.tau = coupling_time
        self.k_B = 0.00831446  # Boltzmann constant in kJ/(mol·K)
    
    def get_temperature(self, particles: List[Particle]) -> float:
        """
        Calculate instantaneous temperature from kinetic energy.
        T = 2<E_k> / (N_df * k_B)
        """
        total_ke = sum(p.kinetic_energy() for p in particles)
        n_degrees = 3 * len(particles)  # 3 per particle in 3D
        return (2 * total_ke) / (n_degrees * self.k_B)
    
    def apply(self, particles: List[Particle], dt: float):
        """Rescale velocities to approach target temperature"""
        T_current = self.get_temperature(particles)
        
        if T_current < 1e-6:  # Avoid division by zero
            return
        
        lambda_factor = np.sqrt(
            1 + (dt / self.tau) * (self.T_target / T_current - 1)
        )
        
        for p in particles:
            p.velocity *= lambda_factor
```

---

## **4. Chemistry Engine**

### **4.1 Reaction Kinetics**

```python
"""
nommo/chemistry/kinetics.py

Arrhenius kinetics and collision theory.
Reference: Steinfeld et al. (1999) Chemical Kinetics and Dynamics
"""

import numpy as np

class ArrheniusKinetics:
    """
    Reaction rate: k(T) = A * exp(-E_a / RT)
    
    Where:
    - A: Pre-exponential factor (frequency factor)
    - E_a: Activation energy (kJ/mol)
    - R: Gas constant (8.314 J/(mol·K))
    - T: Temperature (K)
    """
    
    def __init__(self, activation_energy: float, pre_exponential: float = 1e13):
        """
        Args:
            activation_energy: Barrier to reaction (kJ/mol)
            pre_exponential: Attempt frequency (1/ps)
        """
        self.E_a = activation_energy
        self.A = pre_exponential
        self.R = 0.008314  # kJ/(mol·K)
    
    def rate_constant(self, temperature: float) -> float:
        """Calculate rate constant at temperature T"""
        return self.A * np.exp(-self.E_a / (self.R * temperature))
    
    def reaction_probability(self, temperature: float, dt: float) -> float:
        """Probability of reaction in timestep dt"""
        k = self.rate_constant(temperature)
        return 1 - np.exp(-k * dt)
    
    def has_sufficient_energy(
        self,
        particle1,
        particle2,
        temperature: float
    ) -> bool:
        """
        Collision theory: Check if collision has E > E_a
        
        Uses relative kinetic energy of collision.
        """
        # Relative velocity
        v_rel = particle1.velocity - particle2.velocity
        v_rel_mag = np.linalg.norm(v_rel)
        
        # Reduced mass
        m_reduced = (particle1.mass * particle2.mass) / (
            particle1.mass + particle2.mass
        )
        
        # Collision energy (convert to kJ/mol)
        E_collision = 0.5 * m_reduced * v_rel_mag**2
        E_collision_kjmol = E_collision * 0.01 * 6.022e23 / 1000
        
        # Boltzmann factor
        if E_collision_kjmol < self.E_a:
            prob = np.exp(-(self.E_a - E_collision_kjmol) / (self.R * temperature))
            return np.random.random() < prob
        
        return True
```

### **4.2 Bond Management**

```python
"""
nommo/chemistry/bonds.py

Bond formation and breaking based on thermodynamics.
"""

from dataclasses import dataclass
from typing import Optional
import numpy as np

@dataclass
class Bond:
    """Represents a chemical bond between two particles"""
    particle1_id: str
    particle2_id: str
    energy: float      # kJ/mol
    age: int = 0       # Ticks since formation

class BondManager:
    """
    Manages bond formation and breaking.
    
    Formation criteria:
    1. Particles within bonding distance
    2. Sufficient collision energy (E > E_activation)
    3. ΔG < 0 (thermodynamically favorable)
    4. Available bonding sites
    
    Breaking criteria:
    1. Thermal fluctuation (Boltzmann probability)
    2. Energetic collision
    3. External energy input
    """
    
    def __init__(
        self,
        bond_energy: float,
        activation_energy: float,
        bonding_distance: float
    ):
        self.bond_energy = bond_energy
        self.activation_energy = activation_energy
        self.bonding_distance = bonding_distance
        self.k_B = 0.00831446  # kJ/(mol·K)
    
    def can_form_bond(
        self,
        p1,
        p2,
        temperature: float
    ) -> bool:
        """Check if two particles can form a bond"""
        # Must be close enough
        distance = p1.distance_to(p2)
        if distance > self.bonding_distance:
            return False
        
        # Must have available bonding sites
        if not p1.can_bond_with(p2):
            return False
        
        # Must have sufficient energy (activation)
        # Use collision theory
        return self._has_activation_energy(p1, p2, temperature)
    
    def break_bond_probability(self, bond: Bond, temperature: float) -> float:
        """
        Probability of bond breaking due to thermal fluctuation.
        P_break = exp(-E_bond / k_B*T)
        """
        return np.exp(-bond.energy / (self.k_B * temperature))
```

---

## **5. Analysis & Detection**

### **5.1 Complexity Metrics**

```python
"""
nommo/analysis/metrics.py

Quantify complexity, diversity, and emergence.
References:
- Adami (2002) "What is complexity?"
- Shannon (1948) Information theory
"""

import numpy as np
from collections import Counter
from typing import List, Dict
from nommo.core.particle import Particle

class ComplexityAnalyzer:
    """Calculate various complexity measures"""
    
    @staticmethod
    def shannon_entropy(particles: List[Particle]) -> float:
        """
        Shannon entropy of particle type distribution.
        H = -Σ p_i * log(p_i)
        
        High entropy = diverse population
        Low entropy = dominated by few types
        """
        types = [p.particle_type for p in particles]
        counts = Counter(types)
        total = len(particles)
        
        entropy = 0.0
        for count in counts.values():
            p_i = count / total
            if p_i > 0:
                entropy -= p_i * np.log2(p_i)
        
        return entropy
    
    @staticmethod
    def average_bonds_per_particle(particles: List[Particle]) -> float:
        """Mean number of bonds per particle"""
        if not particles:
            return 0.0
        return sum(len(p.bonds) for p in particles) / len(particles)
    
    @staticmethod
    def max_cluster_size(particles: List[Particle]) -> int:
        """
        Size of largest connected component.
        Indicator of molecular complexity.
        """
        # Build adjacency graph from bonds
        # Find connected components (DFS/BFS)
        # Return size of largest
        pass
    
    @staticmethod
    def lineage_depth(particles: List[Particle]) -> int:
        """
        Maximum generation number.
        Indicates evolutionary time.
        """
        max_gen = 0
        for p in particles:
            gen = 0
            current = p
            while current.parent_id is not None:
                gen += 1
                # Would need to track parent references
            max_gen = max(max_gen, gen)
        return max_gen
```

### **5.2 Emergence Detection**

```python
"""
nommo/analysis/detection.py

Detect autocatalytic sets and self-replicators.
References:
- Kauffman (1986) - Autocatalytic sets
- Hordijk & Steel (2017) - Detection algorithms
"""

from typing import List, Set, Dict
from nommo.core.particle import Particle

class EmergenceDetector:
    """Detect life-like phenomena"""
    
    @staticmethod
    def detect_replicators(
        particles: List[Particle],
        history: List[List[Particle]]
    ) -> List[Particle]:
        """
        Find particles that have successfully copied themselves.
        
        Criteria:
        1. Parent-offspring relationship
        2. Structural similarity (similar bonding pattern)
        3. Multiple generations
        """
        replicators = []
        for p in particles:
            if p.parent_id is not None:
                # Has a parent - check if it's a copy
                if EmergenceDetector._is_similar_structure(p, p.parent_id, particles):
                    replicators.append(p)
        return replicators
    
    @staticmethod
    def detect_autocatalytic_set(
        reactions: List[Dict],
        particles: List[Particle]
    ) -> bool:
        """
        Detect if a self-sustaining reaction network exists.
        
        An autocatalytic set is a collection of molecules where:
        - Each molecule can be produced by reactions within the set
        - The set collectively catalyzes its own formation
        
        Algorithm from Hordijk & Steel (2017)
        """
        # Build reaction network graph
        # Check for catalytic closure
        # Return True if found
        pass
    
    @staticmethod
    def exponential_growth_detected(
        population_history: List[int],
        window: int = 100
    ) -> bool:
        """
        Check if population shows exponential growth.
        Signature of replication.
        """
        if len(population_history) < window:
            return False
        
        recent = population_history[-window:]
        # Fit exponential: N(t) = N0 * exp(r*t)
        # Check R² > threshold
        pass
```

---

## **6. CLI Interface**

### **6.1 Main CLI Structure**

```python
"""
nommo/cli/main.py

Main CLI entry point using Click framework.
"""

import click
from rich.console import Console
from rich.table import Table

console = Console()

@click.group()
@click.version_option(version='0.1.0')
def cli():
    """
    Nommo Engine - A scientific life emergence simulator.
    
    Explore how self-replicating systems emerge from non-living
    chemistry using real physics and thermodynamics.
    """
    pass

# Add command groups
from nommo.cli import universe, run, analyze

cli.add_command(universe.universe)
cli.add_command(run.run_cmd)
cli.add_command(analyze.analyze)

if __name__ == '__main__':
    cli()
```

### **6.2 Universe Commands**

```python
"""
nommo/cli/universe.py

Commands for managing universes.
"""

import click
from rich.console import Console
from rich.prompt import Prompt, Confirm
from rich.table import Table

console = Console()

@click.group()
def universe():
    """Create and manage simulation universes"""
    pass

@universe.command()
@click.argument('name', required=False)
@click.option('--preset', help='Use a preset configuration')
@click.option('--interactive', '-i', is_flag=True, help='Interactive setup')
def create(name, preset, interactive):
    """
    Create a new universe.
    
    Examples:
        nommo universe create earth-1 --preset earth-like
        nommo universe create -i
    """
    if interactive:
        console.print("[bold cyan]Creating new universe[/bold cyan]")
        name = Prompt.ask("Universe name")
        preset = Prompt.ask(
            "Preset",
            choices=["earth-like", "high-energy", "minimal", "custom"],
            default="earth-like"
        )
    
    # Create universe with parameters
    # Save to data/universes/
    console.print(f"✓ Universe '{name}' created", style="bold green")

@universe.command()
@click.option('--verbose', '-v', is_flag=True)
def list(verbose):
    """
    List all universes.
    
    Examples:
        nommo universe list
        nommo universe list -v
    """
    # Load universes from data/universes/
    table = Table(title="Active Universes")
    table.add_column("ID", style="cyan")
    table.add_column("Name", style="magenta")
    table.add_column("Status", style="green")
    table.add_column("Particles")
    table.add_column("Tick")
    
    # Add rows from loaded data
    console.print(table)

@universe.command()
@click.argument('universe_id')
def show(universe_id):
    """
    Show universe details.
    
    Examples:
        nommo universe show earth-1
    """
    # Load universe
    # Display parameters, metrics, status
    pass

@universe.command()
@click.argument('universe_id')
@click.option('--force', '-f', is_flag=True)
def delete(universe_id, force):
    """Delete a universe"""
    if not force:
        if not Confirm.ask(f"Delete universe '{universe_id}'?"):
            return
    
    # Delete from storage
    console.print(f"✓ Universe '{universe_id}' deleted", style="bold green")
```

### **6.3 Run Commands**

```python
"""
nommo/cli/run.py

Commands for running simulations.
"""

import click
from rich.console import Console
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn
from rich.live import Live
from rich.panel import Panel
from rich.layout import Layout

console = Console()

@click.group(name='run')
def run_cmd():
    """Execute simulations"""
    pass

@run_cmd.command()
@click.argument('universe_id')
@click.option('--ticks', type=int, default=1000, help='Number of steps')
@click.option('--watch', is_flag=True, help='Live progress display')
@click.option('--speed', type=float, default=1.0, help='Simulation speed multiplier')
def start(universe_id, ticks, watch, speed):
    """
    Run a simulation.
    
    Examples:
        nommo run start earth-1 --ticks 10000 --watch
        nommo run start earth-1 --ticks 5000 --speed 2.0
    """
    # Load universe
    # Run simulation with progress bar
    
    if watch:
        with Progress(
            SpinnerColumn(),
            TextColumn("[progress.description]{task.description}"),
            BarColumn(),
            TextColumn("[progress.percentage]{task.percentage:>3.0f}%"),
        ) as progress:
            task = progress.add_task(f"Running {universe_id}", total=ticks)
            
            for i in range(ticks):
                # universe.step()
                # Update metrics display
                progress.update(task, advance=1)
                
                if i % 100 == 0:
                    # Display current stats
                    pass
    else:
        # Run without display
        pass
    
    console.print("✓ Simulation complete", style="bold green")
```

### **6.4 Analysis Commands**

```python
"""
nommo/cli/analyze.py

Commands for analyzing simulation results.
"""

import click
from rich.console import Console
from rich.table import Table
import matplotlib.pyplot as plt

console = Console()

@click.group()
def analyze():
    """Analyze simulation results"""
    pass

@analyze.command()
@click.argument('universe_id')
def stats(universe_id):
    """
    Show statistics.
    
    Examples:
        nommo analyze stats earth-1
    """
    # Load universe and history
    # Calculate metrics
    
    table = Table(title=f"Statistics for {universe_id}")
    table.add_column("Metric", style="cyan")
    table.add_column("Value", style="magenta")
    
    # Add metrics
    table.add_row("Total Particles", "1,247")
    table.add_row("Complexity", "3.2")
    table.add_row("Replicators", "6")
    
    console.print(table)

@analyze.command()
@click.argument('universe_id')
@click.argument('metric')
@click.option('--output', help='Save plot to file')
def plot(universe_id, metric, output):
    """
    Plot metric over time.
    
    Examples:
        nommo analyze plot earth-1 complexity
        nommo analyze plot earth-1 particles --output plot.png
    """
    # Load history
    # Create plot with matplotlib
    # Save or display
    pass

@analyze.command()
@click.argument('universe_ids', nargs=-1, required=True)
def compare(universe_ids):
    """
    Compare multiple universes.
    
    Examples:
        nommo analyze compare earth-1 earth-2 earth-3
    """
    # Load all universes
    # Create comparison table/plots
    pass
```

---

## **7. Dependencies**

```toml
# pyproject.toml

[tool.poetry]
name = "nommo-engine"
version = "0.1.0"
description = "Scientific life emergence simulator"
authors = ["Your Name <you@example.com>"]
readme = "README.md"
packages = [{include = "nommo"}]

[tool.poetry.dependencies]
python = "^3.11"
numpy = "^1.26.0"
scipy = "^1.11.0"
matplotlib = "^3.8.0"
click = "^8.1.0"
rich = "^13.7.0"
pydantic = "^2.5.0"
pyyaml = "^6.0.0"
h5py = "^3.10.0"  # For efficient data storage

[tool.poetry.group.dev.dependencies]
pytest = "^7.4.0"
pytest-cov = "^4.1.0"
black = "^23.12.0"
ruff = "^0.1.0"
mypy = "^1.7.0"

[tool.poetry.scripts]
nommo = "nommo.cli.main:cli"

[build-system]
requires = ["poetry-core"]
build-backend = "poetry.core.masonry.api"
```

---

## **8. Implementation Priorities**

### **Phase 1: Core Simulation (Week 1-2)**

```markdown
Priority 1: Basic Physics
- [ ] Particle class with position, velocity, forces
- [ ] Lennard-Jones force calculation
- [ ] Velocity Verlet integrator
- [ ] Periodic boundary conditions
- [ ] Test: 100 particles bouncing around

Priority 2: Basic Chemistry
- [ ] Bond formation (distance-based)
- [ ] Bond breaking (energy-based)
- [ ] Simple reaction rules (A + B → AB)
- [ ] Test: Molecules forming and breaking

Priority 3: Temperature Control
- [ ] Berendsen thermostat
- [ ] Temperature calculation
- [ ] Test: Maintain constant temperature

Deliverable: Working molecular dynamics simulation
```

### **Phase 2: Chemistry & Emergence (Week 3-4)**

```markdown
Priority 4: Advanced Chemistry
- [ ] Arrhenius kinetics
- [ ] Activation energy barriers
- [ ] Multiple particle types
- [ ] Test: Temperature-dependent reactions

Priority 5: Complexity Metrics
- [ ] Shannon entropy
- [ ] Bond count tracking
- [ ] Cluster size analysis
- [ ] Test: Complexity increases over time

Priority 6: Basic CLI
- [ ] Click framework setup
- [ ] Universe create/list/show commands
- [ ] Run command with progress bar
- [ ] Test: Create and run universe via CLI

Deliverable: Chemical reactions with metrics tracking
```

### **Phase 3: Analysis & Detection (Week 5-6)**

```markdown
Priority 7: Emergence Detection
- [ ] Replicator detection
- [ ] Exponential growth detection
- [ ] Lineage tracking
- [ ] Test: Detect when replication emerges

Priority 8: Analysis Tools
- [ ] Plotting with matplotlib
- [ ] Statistics calculation
- [ ] Data export (CSV, JSON)
- [ ] Test: Generate analysis plots

Priority 9: Full CLI
- [ ] All universe commands
- [ ] All run commands
- [ ] All analyze commands
- [ ] Interactive prompts with Rich

Deliverable: Complete analysis pipeline
```

### **Phase 4: Optimization & Polish (Week 7-8)**

```markdown
Priority 10: Performance
- [ ] Spatial indexing (cell lists or quadtree)
- [ ] Neighbor lists
- [ ] Profiling and optimization
- [ ] Test: Handle 10,000+ particles

Priority 11: Presets & Documentation
- [ ] Earth-like preset
- [ ] High-energy preset
- [ ] Minimal preset
- [ ] Complete documentation

Priority 12: Testing
- [ ] Unit tests for all modules
- [ ] Integration tests
- [ ] Scientific validation (reproduce known results)

Deliverable: Production-ready simulator
```

---

## **9. Testing Strategy**

### **Unit Tests**

```python
# tests/test_forces.py

import pytest
import numpy as np
from nommo.physics.forces import LennardJonesForce

def test_lj_force_at_sigma():
    """At r = σ, force should be zero"""
    lj = LennardJonesForce(epsilon=1.0, sigma=1.0, cutoff=2.5)
    r_vec = np.array([1.0, 0.0, 0.0])
    force, potential = lj.calculate(r_vec, 1.0)
    
    assert np.allclose(force, 0.0, atol=1e-10)
    assert np.isclose(potential, 0.0, atol=1e-10)

def test_lj_force_attractive():
    """At r > σ, force should be attractive"""
    lj = LennardJonesForce(epsilon=1.0, sigma=1.0, cutoff=2.5)
    r_vec = np.array([1.5, 0.0, 0.0])
    force, potential = lj.calculate(r_vec, 1.5)
    
    # Force should point toward origin (negative x)
    assert force[0] < 0
    assert potential < 0  # Attractive well

def test_lj_force_repulsive():
    """At r < σ, force should be repulsive"""
    lj = LennardJonesForce(epsilon=1.0, sigma=1.0, cutoff=2.5)
    r_vec = np.array([0.5, 0.0, 0.0])
    force, potential = lj.calculate(r_vec, 0.5)
    
    # Force should point away from origin (positive x)
    assert force[0] > 0
    assert potential > 0  # Repulsive wall
```

### **Integration Tests**

```python
# tests/test_universe.py

def test_energy_conservation():
    """Without thermostat, total energy should be conserved"""
    params = UniverseParameters(
        name="test",
        temperature=300.0,
        use_thermostat=False,
        # ... other params
    )
    universe = Universe(params)
    
    initial_energy = universe.total_energy()
    
    for _ in range(1000):
        universe.step()
    
    final_energy = universe.total_energy()
    
    assert np.isclose(initial_energy, final_energy, rtol=0.01)

def test_temperature_control():
    """Thermostat should maintain target temperature"""
    params = UniverseParameters(
        name="test",
        temperature=300.0,
        use_thermostat=True,
        # ... other params
    )
    universe = Universe(params)
    
    for _ in range(1000):
        universe.step()
    
    current_temp = universe.metrics.current_temperature
    
    assert np.isclose(current_temp, 300.0, rtol=0.05)
```

---

## **10. Scientific Validation**

```markdown
# Validation Experiments

## Test 1: Reproduce Ideal Gas Law
Simulate non-interacting particles (ε = 0)
Measure: PV = NkT relationship
Expected: Should match within 1%

## Test 2: Phase Transition
Vary temperature from 100K to 1000K
Measure: Particle clustering
Expected: Gas → Liquid → Solid behavior

## Test 3: Reaction Rate vs Temperature
Measure reaction rate at different T
Plot: ln(k) vs 1/T (Arrhenius plot)
Expected: Linear relationship

## Test 4: Emergence Reproducibility
Run 100 identical universes
Measure: % that develop replicators
Expected: Consistent emergence rate

## Test 5: Parameter Sensitivity
Vary one parameter at a time
Measure: Impact on emergence time
Expected: Identify critical parameters
```

---

## **11. Example Usage**

```bash
# Create a new universe
nommo universe create earth-1 --preset earth-like

# Run simulation with live display
nommo run start earth-1 --ticks 10000 --watch

# Analyze results
nommo analyze stats earth-1
nommo analyze plot earth-1 complexity --output complexity.png

# Compare multiple universes
nommo universe create hot --preset high-energy
nommo run start hot --ticks 10000
nommo analyze compare earth-1 hot

# Interactive creation
nommo universe create -i
```

---

## **12. Success Criteria**

```markdown
✅ Minimal Success (MVP):
- Particles interact via realistic forces
- Bonds form and break thermodynamically
- Temperature control works
- Basic CLI functional
- Can run 1000 particle simulation

✅ Full Success (v1.0):
- Self-replicating systems emerge
- Multiple universe types supported
- Complete analysis tools
- Scientific validation passed
- Performance: 10,000 particles at 10+ fps
- Documentation complete

✅ Excellence (Beyond):
- Published insights from simulations
- Used by researchers/educators
- Beautiful visualizations
- Community presets
- Reproducible science
```

---

This specification provides Claude Code with:
1. **Complete architecture** - Every module defined
2. **Real science** - Actual equations and algorithms
3. **Clear implementation path** - Phased approach
4. **Testable components** - Unit and integration tests
5. **User interface** - Complete CLI specification
6. **Validation strategy** - How to verify correctness

**Ready to provide to Claude Code!**
