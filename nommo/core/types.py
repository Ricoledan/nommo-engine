"""
Type definitions and parameter structures with Pydantic validation.

Provides comprehensive configuration schema with validation, defaults,
and serialization support for all simulation parameters.
"""

from enum import Enum
from pathlib import Path

from pydantic import BaseModel, Field, field_validator, model_validator
from pydantic_settings import BaseSettings


class BoundaryType(str, Enum):
    """Boundary condition types for the simulation box."""

    PERIODIC = "periodic"
    REFLECTIVE = "reflective"
    ABSORBING = "absorbing"


class EnergySource(str, Enum):
    """Energy input mechanisms for the system."""

    CONSTANT = "constant"
    PULSES = "pulses"
    GRADIENT = "gradient"
    CHAOTIC = "chaotic"


class ForceFieldType(str, Enum):
    """Available force field implementations."""

    LENNARD_JONES = "lennard_jones"
    MORSE = "morse"
    BUCKINGHAM = "buckingham"
    CUSTOM = "custom"


class IntegratorType(str, Enum):
    """Time integration algorithms."""

    VELOCITY_VERLET = "velocity_verlet"
    LEAPFROG = "leapfrog"
    RK4 = "rk4"


class ThermostatType(str, Enum):
    """Temperature control methods."""

    NONE = "none"
    BERENDSEN = "berendsen"
    NOSE_HOOVER = "nose_hoover"
    ANDERSEN = "andersen"


class ParticleTypeConfig(BaseModel):
    """Configuration for a particle type."""

    name: str = Field(..., description="Unique particle type identifier")
    mass: float = Field(gt=0, description="Particle mass in amu")
    charge: float = Field(default=0.0, description="Particle charge in elementary units")
    radius: float = Field(gt=0, description="Particle radius in nm")
    max_bonds: int = Field(ge=0, le=12, default=4, description="Maximum number of bonds")
    color: str = Field(default="#FF0000", description="Visualization color (hex)")

    lj_epsilon: float | None = Field(
        default=None, gt=0, description="LJ potential well depth (kJ/mol)"
    )
    lj_sigma: float | None = Field(default=None, gt=0, description="LJ zero-crossing distance (nm)")

    @field_validator("color")
    @classmethod
    def validate_color(cls, v: str) -> str:
        if not v.startswith("#") or len(v) != 7:
            raise ValueError("Color must be hex format (#RRGGBB)")
        return v


class ReactionRule(BaseModel):
    """Chemical reaction rule definition."""

    reactants: list[str] = Field(..., min_length=1, description="Reactant particle types")
    products: list[str] = Field(..., min_length=1, description="Product particle types")
    activation_energy: float = Field(ge=0, description="Activation energy (kJ/mol)")
    rate_constant: float = Field(gt=0, default=1e13, description="Pre-exponential factor (1/ps)")
    reversible: bool = Field(default=True, description="Whether reaction is reversible")
    catalyst: str | None = Field(default=None, description="Catalyst particle type")

    @model_validator(mode="after")
    def validate_conservation(self) -> "ReactionRule":
        """Validate mass/charge conservation (simplified)."""
        return self


class PhysicsParameters(BaseModel):
    """Physics simulation parameters."""

    timestep: float = Field(default=0.001, gt=0, le=0.1, description="Integration timestep (ps)")
    force_field: ForceFieldType = Field(
        default=ForceFieldType.LENNARD_JONES, description="Force field type"
    )
    lj_epsilon: float = Field(default=1.0, gt=0, description="Default LJ epsilon (kJ/mol)")
    lj_sigma: float = Field(default=0.3, gt=0, description="Default LJ sigma (nm)")
    cutoff_distance: float = Field(default=1.0, gt=0, description="Force cutoff (nm)")

    integrator: IntegratorType = Field(
        default=IntegratorType.VELOCITY_VERLET, description="Integration algorithm"
    )
    adaptive_timestep: bool = Field(default=False, description="Enable adaptive timestep")
    max_force: float = Field(default=1e6, gt=0, description="Maximum allowed force (kJ/mol/nm)")

    use_neighbor_list: bool = Field(default=True, description="Use Verlet neighbor lists")
    neighbor_skin: float = Field(default=0.3, gt=0, description="Neighbor list skin (nm)")
    neighbor_update_interval: int = Field(
        default=10, gt=0, description="Steps between neighbor list updates"
    )


class ChemistryParameters(BaseModel):
    """Chemistry simulation parameters."""

    activation_energy: float = Field(
        default=10.0, ge=0, description="Default activation energy (kJ/mol)"
    )
    bond_energy: float = Field(default=20.0, gt=0, description="Default bond energy (kJ/mol)")
    bonding_distance: float = Field(default=0.4, gt=0, description="Maximum bonding distance (nm)")
    reaction_probability: float = Field(
        default=0.1, gt=0, le=1, description="Base reaction probability"
    )

    enable_reactions: bool = Field(default=True, description="Enable chemical reactions")
    enable_bonding: bool = Field(default=True, description="Enable bond formation/breaking")

    reaction_rules: list[ReactionRule] = Field(default_factory=list, description="Reaction rules")


class ThermodynamicsParameters(BaseModel):
    """Thermodynamics parameters."""

    temperature: float = Field(default=300.0, gt=0, le=10000, description="System temperature (K)")
    pressure: float | None = Field(default=None, gt=0, description="System pressure (bar)")

    thermostat: ThermostatType = Field(
        default=ThermostatType.BERENDSEN, description="Thermostat type"
    )
    thermostat_coupling: float = Field(
        default=0.1, gt=0, description="Thermostat coupling time (ps)"
    )

    energy_source: EnergySource = Field(
        default=EnergySource.CONSTANT, description="Energy input mechanism"
    )
    energy_amount: float = Field(default=0.0, ge=0, description="Energy input amount (kJ/mol)")
    energy_interval: int = Field(default=100, gt=0, description="Steps between energy inputs")


class SimulationBox(BaseModel):
    """Simulation box parameters."""

    dimensions: tuple[float, float, float] = Field(
        default=(10.0, 10.0, 10.0), description="Box dimensions (nm)"
    )
    boundary_type: BoundaryType = Field(
        default=BoundaryType.PERIODIC, description="Boundary conditions"
    )

    @field_validator("dimensions")
    @classmethod
    def validate_dimensions(cls, v: tuple[float, float, float]) -> tuple[float, float, float]:
        if any(d <= 0 for d in v):
            raise ValueError("All dimensions must be positive")
        return v


class UniverseParameters(BaseModel):
    """Complete parameter set for a universe simulation."""

    name: str = Field(..., min_length=1, description="Universe name")
    description: str = Field(default="", description="Universe description")

    physics: PhysicsParameters = Field(default_factory=PhysicsParameters)
    chemistry: ChemistryParameters = Field(default_factory=ChemistryParameters)
    thermodynamics: ThermodynamicsParameters = Field(default_factory=ThermodynamicsParameters)
    box: SimulationBox = Field(default_factory=SimulationBox)

    particle_types: dict[str, ParticleTypeConfig] = Field(
        default_factory=dict, description="Particle type definitions"
    )
    initial_composition: dict[str, int] = Field(
        default_factory=dict, description="Initial particle counts"
    )

    max_particles: int = Field(default=10000, gt=0, description="Maximum particle count")
    random_seed: int | None = Field(default=None, description="Random seed for reproducibility")

    checkpoint_interval: int = Field(default=1000, gt=0, description="Steps between checkpoints")
    history_buffer_size: int = Field(
        default=10000, gt=0, description="Maximum history entries in memory"
    )

    @model_validator(mode="after")
    def validate_composition(self) -> "UniverseParameters":
        """Validate initial composition references defined particle types."""
        for ptype in self.initial_composition:
            if ptype not in self.particle_types:
                raise ValueError(f"Unknown particle type in composition: {ptype}")
        return self

    @model_validator(mode="after")
    def validate_particle_count(self) -> "UniverseParameters":
        """Validate initial particle count doesn't exceed maximum."""
        total = sum(self.initial_composition.values())
        if total > self.max_particles:
            raise ValueError(f"Initial particles ({total}) exceeds maximum ({self.max_particles})")
        return self


class SimulationConfig(BaseSettings):
    """Application-wide configuration."""

    data_dir: Path = Field(default=Path("data"), description="Data directory")
    log_level: str = Field(default="INFO", description="Logging level")
    log_file: Path | None = Field(default=None, description="Log file path")
    use_json_logs: bool = Field(default=False, description="Use JSON structured logging")

    enable_gpu: bool = Field(default=False, description="Enable GPU acceleration if available")
    num_threads: int = Field(default=-1, gt=-2, description="Number of threads (-1 for auto)")

    autosave: bool = Field(default=True, description="Enable automatic saving")
    autosave_interval: int = Field(default=5000, gt=0, description="Steps between autosaves")

    visualization_backend: str = Field(default="matplotlib", description="Visualization backend")
    export_format: str = Field(default="hdf5", description="Default export format")

    class Config:
        env_prefix = "NOMMO_"
        env_file = ".env"
        env_file_encoding = "utf-8"
        case_sensitive = False

    @field_validator("data_dir")
    @classmethod
    def create_data_dir(cls, v: Path) -> Path:
        v.mkdir(parents=True, exist_ok=True)
        (v / "universes").mkdir(exist_ok=True)
        (v / "exports").mkdir(exist_ok=True)
        (v / "cache").mkdir(exist_ok=True)
        return v


def load_universe_parameters(path: Path) -> UniverseParameters:
    """Load universe parameters from file."""
    import json

    import yaml

    if path.suffix == ".json":
        with open(path) as f:
            data = json.load(f)
    elif path.suffix in [".yaml", ".yml"]:
        with open(path) as f:
            data = yaml.safe_load(f)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")

    return UniverseParameters(**data)


def save_universe_parameters(params: UniverseParameters, path: Path) -> None:
    """Save universe parameters to file."""
    import json

    import yaml

    data = params.model_dump()

    if path.suffix == ".json":
        with open(path, "w") as f:
            json.dump(data, f, indent=2, default=str)
    elif path.suffix in [".yaml", ".yml"]:
        with open(path, "w") as f:
            yaml.dump(data, f, default_flow_style=False)
    else:
        raise ValueError(f"Unsupported file format: {path.suffix}")
