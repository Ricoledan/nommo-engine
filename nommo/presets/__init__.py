"""Preset universe configurations."""

from nommo.core.types import UniverseParameters, ParticleTypeConfig


def get_preset(name: str) -> UniverseParameters:
    """Get a preset universe configuration."""
    presets = {
        "earth-like": earth_like_preset(),
        "high-energy": high_energy_preset(),
        "minimal": minimal_preset(),
    }
    return presets.get(name, minimal_preset())


def earth_like_preset() -> UniverseParameters:
    """Earth-like conditions preset."""
    return UniverseParameters(
        name="earth-like",
        description="Earth-like temperature and pressure conditions",
        particle_types={
            "monomer_A": ParticleTypeConfig(
                name="monomer_A",
                mass=50.0,
                radius=0.2,
                max_bonds=2,
                color="#FF6B6B"
            ),
            "monomer_B": ParticleTypeConfig(
                name="monomer_B",
                mass=75.0,
                radius=0.25,
                max_bonds=3,
                color="#4ECDC4"
            ),
            "catalyst": ParticleTypeConfig(
                name="catalyst",
                mass=100.0,
                radius=0.3,
                max_bonds=4,
                color="#95E77E"
            ),
        },
        initial_composition={
            "monomer_A": 200,
            "monomer_B": 200,
            "catalyst": 50,
        }
    )


def high_energy_preset() -> UniverseParameters:
    """High energy environment preset."""
    params = earth_like_preset()
    params.name = "high-energy"
    params.description = "High temperature, energetic environment"
    params.thermodynamics.temperature = 500.0
    params.thermodynamics.energy_amount = 10.0
    params.chemistry.activation_energy = 5.0
    return params


def minimal_preset() -> UniverseParameters:
    """Minimal test preset."""
    return UniverseParameters(
        name="minimal",
        description="Minimal configuration for testing",
        particle_types={
            "particle": ParticleTypeConfig(
                name="particle",
                mass=50.0,
                radius=0.2,
                max_bonds=2
            ),
        },
        initial_composition={
            "particle": 100,
        }
    )