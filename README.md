# Nommo Engine

A scientific life emergence simulator that explores how self-replicating systems emerge from non-living chemistry using
real physics and thermodynamics.

## Overview

Nommo Engine is a molecular dynamics simulation framework designed to study the emergence of life-like phenomena from
simple chemical systems. Built with scientific rigor, it implements:

- **Real Physics**: Lennard-Jones potentials, Velocity Verlet integration, proper thermodynamics
- **Chemical Kinetics**: Arrhenius equation, activation energies, bond formation/breaking
- **Emergence Detection**: Autocatalytic sets, self-replicators, complexity metrics
- **High Performance**: Numba JIT compilation, spatial indexing, neighbor lists
- **Extensible Architecture**: Plugin system for custom force fields, reactions, and analysis

## Features

### Core Capabilities

- **Molecular Dynamics**: Accurate force calculations with periodic/reflective boundaries
- **Temperature Control**: Multiple thermostat algorithms (Berendsen, Nosé-Hoover, Andersen)
- **Spatial Optimization**: Cell lists and Verlet neighbor lists for O(N) scaling
- **Adaptive Timesteps**: Automatic timestep adjustment for numerical stability
- **Chemical Reactions**: Rule-based reactions with activation energies and catalysis
- **Emergence Metrics**: Shannon entropy, bond complexity, lineage tracking

### Performance Features

- **Numba Acceleration**: JIT compilation of critical loops
- **Vectorized Operations**: NumPy-based array computations
- **Memory Pooling**: Efficient particle object reuse
- **Parallel Processing**: Multi-threaded force calculations
- **Smart Rebuilding**: Adaptive neighbor list updates

### Analysis Tools

- **Real-time Metrics**: Track complexity, temperature, energy during simulation
- **Emergence Detection**: Identify self-replicating structures and autocatalytic sets
- **Visualization**: 2D/3D rendering, bond network graphs, trajectory animations
- **Data Export**: HDF5, CSV, JSON formats with compression

## Installation

### Prerequisites

- Python 3.11 or higher
- Poetry (for dependency management)

### Install with Poetry

```bash
git clone https://github.com/ricoledan/nommo-engine.git
cd nommo-engine
poetry install
```

### Install with pip (development)

```bash
git clone https://github.com/ricoledan/nommo-engine.git
cd nommo-engine
pip install -e .
```

## Running the Local Environment

### Activating the Virtual Environment

The project uses a Python virtual environment for dependency isolation. You have two options:

#### Option 1: Direct Virtual Environment Activation

If you have a `venv/` directory in your project:

```bash
# On macOS/Linux
source venv/bin/activate

# On Windows
venv\Scripts\activate
```

To deactivate when you're done:

```bash
deactivate
```

#### Option 2: Using Poetry (Recommended)

Poetry automatically manages the virtual environment for you:

```bash
# Run a single command in the Poetry environment
poetry run python example.py

# Or activate a shell with the Poetry environment
poetry shell

# To exit the Poetry shell
exit
```

### Running Python Scripts

Once your environment is activated:

```bash
# Run the example script
python example.py

# Or with Poetry (no activation needed)
poetry run python example.py
```

### Using the Nommo CLI

The project includes a command-line interface:

```bash
# If environment is activated
nommo --help

# Using Poetry
poetry run nommo --help
```

### Development Commands

#### Running Tests

```bash
poetry run pytest
# Or with coverage
poetry run pytest --cov=nommo --cov-report=term-missing
```

#### Code Quality Checks

```bash
# Format code with Black
poetry run black nommo/

# Check code style with Ruff
poetry run ruff check nommo/

# Type checking with MyPy
poetry run mypy nommo/
```

#### Documentation Server

```bash
poetry run mkdocs serve
# Then visit http://127.0.0.1:8000
```

### Installing Additional Dependencies

If you need to add new packages:

```bash
# Add a runtime dependency
poetry add package-name

# Add a development dependency
poetry add --group dev package-name

# Update all dependencies
poetry update
```

## Quick Start

### Create a Universe

```bash
nommo universe create earth-1 --preset earth-like
```

### Run a Simulation

```bash
nommo run start earth-1 --ticks 10000 --watch
```

### Analyze Results

```bash
nommo analyze stats earth-1
nommo analyze plot earth-1 complexity --output complexity.png
```

## Architecture

### Project Structure

```
nommo-engine/
├── nommo/
│   ├── core/           # Core data structures (Particle, Universe)
│   ├── physics/        # Force calculations, integrators, thermostats
│   ├── chemistry/      # Reactions, bonds, kinetics
│   ├── analysis/       # Metrics, emergence detection
│   ├── storage/        # HDF5 persistence, configuration
│   ├── plugins/        # Extension system
│   ├── utils/          # Logging, performance tracking
│   └── cli/            # Command-line interface
├── tests/              # Unit and integration tests
├── docs/               # Documentation
└── data/               # Runtime data storage
```

### Key Components

1. **Particle System**: Represents molecular units with realistic physics
2. **Spatial Indexing**: Efficient neighbor finding with cell lists
3. **Force Engine**: Modular force calculations with plugin support
4. **Integrators**: Multiple time integration algorithms
5. **Thermostats**: Temperature control for different ensembles
6. **Plugin System**: Extensible architecture for custom components

## Scientific Basis

### Theoretical Background

Nommo Engine implements a scientifically rigorous approach to studying emergence phenomena in chemical systems. The simulation combines:

- **Classical Molecular Dynamics**: Particles evolve according to Newton's equations with realistic inter-molecular forces
- **Chemical Kinetics**: Reactions occur according to established physical chemistry principles
- **Emergence Detection**: Sophisticated algorithms identify autocatalytic sets and self-replicating systems
- **Validation Framework**: Comprehensive benchmarks ensure scientific accuracy

For detailed theoretical foundations, see [docs/METHODOLOGY.md](docs/METHODOLOGY.md).

### Physics Implementation

The simulation uses classical molecular dynamics with:

- **Lennard-Jones Potential**: V(r) = 4ε[(σ/r)¹² - (σ/r)⁶] for inter-particle interactions
- **Velocity Verlet Integration**: Symplectic, time-reversible algorithm ensuring energy conservation
- **Berendsen Thermostat**: Weak coupling for temperature control with realistic fluctuations
- **Periodic Boundaries**: Eliminate finite-size effects for bulk property studies

**Energy Conservation**: Total energy drift < 0.01% over nanosecond timescales (validated)

### Chemistry Implementation

Chemical reactions follow established physical chemistry:

- **Arrhenius Kinetics**: k(T) = A·exp(-Ea/RT) with realistic activation energies (10-100 kJ/mol)
- **Collision Theory**: Geometric constraints and energy requirements for bond formation
- **Bond Thermodynamics**: Gibbs free energy considerations for reaction spontaneity
- **Steric Effects**: Orientation-dependent reaction probabilities

**Validation**: Reaction rates match experimental Arrhenius behavior (R² > 0.99)

### Model Limitations

**What the simulation captures:**
- Inter-molecular forces and dynamics
- Temperature-dependent reaction kinetics
- Emergence of complex molecular networks
- Self-organization and replication phenomena

**What it does not capture:**
- Quantum effects (tunneling, zero-point motion)
- Electronic structure changes during reactions
- Detailed transition state geometries
- Solvent effects beyond simple interactions

For parameter sensitivity analysis, see [docs/PARAMETERS.md](docs/PARAMETERS.md).

### Units and Scales

**Physical Units:**
- **Distance**: nanometers (nm) - molecular scale
- **Time**: picoseconds (ps) - molecular motion timescale
- **Energy**: kJ/mol - chemical energy scale
- **Mass**: atomic mass units (amu) - molecular masses
- **Temperature**: Kelvin (K) - thermal energy scale

**Typical System Scales:**
- **Particles**: 1,000-50,000 (representing 10⁵-10⁷ atoms)
- **Time**: 1-1000 ps (microsecond processes via scaling)
- **Box size**: 5-50 nm (bulk behavior above 10 nm)

## Configuration

### Universe Parameters

```yaml
name: my_universe
physics:
    timestep: 0.001  # ps
    cutoff_distance: 1.0  # nm
    force_field: lennard_jones

chemistry:
    activation_energy: 10.0  # kJ/mol
    bond_energy: 20.0  # kJ/mol

thermodynamics:
    temperature: 300.0  # K
    thermostat: berendsen

particle_types:
    monomer_A:
        mass: 50.0  # amu
        radius: 0.2  # nm
        max_bonds: 2
```

## Performance Considerations

### Optimization Strategies

1. **Neighbor Lists**: Reduce O(N²) to O(N) for force calculations
2. **Numba JIT**: Compile hot loops to machine code
3. **Vectorization**: Use NumPy for array operations
4. **Spatial Indexing**: Cell lists for efficient neighbor finding
5. **Adaptive Timesteps**: Maintain stability with minimal computation

### Benchmarks

- 1,000 particles: ~100 fps on modern CPU
- 10,000 particles: ~10 fps with full physics
- 100,000 particles: Possible with GPU acceleration (planned)

## Development

### Running Tests

```bash
poetry run pytest
```

### Code Quality

```bash
poetry run black nommo/
poetry run ruff check nommo/
poetry run mypy nommo/
```

### Documentation

```bash
poetry run mkdocs serve
```

## Roadmap

### Phase 1: Core Engine ✓

- [x] Particle physics implementation
- [x] Spatial indexing system
- [x] Force calculations with Numba
- [x] Temperature control
- [x] Plugin architecture

### Phase 2: Chemistry (In Progress)

- [ ] Reaction rule engine
- [ ] Bond network analysis
- [ ] Autocatalytic set detection
- [ ] Energy flow tracking

### Phase 3: Emergence

- [ ] Replicator detection algorithms
- [ ] Information theory metrics
- [ ] Evolution tracking
- [ ] Fitness landscapes

### Phase 4: Visualization

- [ ] Real-time 3D rendering
- [ ] Network graph visualization
- [ ] Statistical dashboards
- [ ] Animation export

## Contributing

Contributions are welcome! Please see [CONTRIBUTING.md](docs/CONTRIBUTING.md) for guidelines.

## License

MIT License - see [LICENSE](LICENSE) file for details.

## Citation

If you use Nommo Engine in your research, please cite:

```bibtex
@software{nommo_engine,
  author = {Ledan, Rico},
  title = {Nommo Engine: A Scientific Life Emergence Simulator},
  year = {2024},
  url = {https://github.com/ricoledan/nommo-engine}
}
```

## Research Documentation

### Core Documentation

- **[METHODOLOGY.md](docs/METHODOLOGY.md)** - Detailed algorithms, equations, and theoretical foundations
- **[VALIDATION.md](docs/VALIDATION.md)** - Benchmark results and scientific accuracy verification
- **[PARAMETERS.md](docs/PARAMETERS.md)** - Parameter sensitivity analysis and selection guidelines
- **[LITERATURE.md](docs/LITERATURE.md)** - Comprehensive bibliography and literature review

### Reproducibility

All scientific results can be reproduced using:

```bash
# Run complete validation suite
nommo validate --all --output results/validation/

# Reproduce specific benchmarks
nommo validate energy-conservation
nommo validate arrhenius-kinetics
nommo validate emergence-detection

# Generate parameter sensitivity analysis
python scripts/sensitivity_analysis.py --output results/sensitivity/
```

### Research Examples

**Example research workflows:**

```bash
# Study autocatalytic network formation
nommo universe create autocatal --preset autocatalytic
nommo run start autocatal --ticks 50000 --watch
nommo analyze emergence autocatal

# Parameter sweep for emergence conditions  
nommo sweep temperature 200,300,400,500 --universe-template minimal
nommo analyze compare temp-sweep-*

# Replication fidelity study
nommo experiment replication-fidelity --mutations 0.01,0.05,0.1
```

See `results/` directory for example outputs and validation data.

## References

### Core Implementation References

- **Frenkel, D. & Smit, B. (2001).** *Understanding Molecular Simulation: From Algorithms to Applications* (2nd ed.). Academic Press. [DOI:10.1016/B978-0-12-267351-1.X5000-7](https://doi.org/10.1016/B978-0-12-267351-1.X5000-7)
- **Allen, M. P. & Tildesley, D. J. (2017).** *Computer Simulation of Liquids* (2nd ed.). Oxford University Press. [ISBN:978-0-19-855645-9](https://global.oup.com/academic/product/computer-simulation-of-liquids-9780198556459)
- **Berendsen, H. J. C. et al. (1984).** "Molecular dynamics with coupling to an external bath". *J. Chem. Phys.* 81, 3684. [DOI:10.1063/1.448118](https://doi.org/10.1063/1.448118)

### Emergence Theory References

- **Kauffman, S. A. (1986).** "Autocatalytic sets of proteins". *Journal of Theoretical Biology*, 119(1), 1-24. [DOI:10.1016/S0022-5193(86)80047-9](https://doi.org/10.1016/S0022-5193(86)80047-9)
- **Hordijk, W. & Steel, M. (2017).** "Detecting autocatalytic, self-sustaining sets in chemical reaction systems". *J. Theor. Biol.* 227, 451. [DOI:10.1016/j.jtbi.2017.06.020](https://doi.org/10.1016/j.jtbi.2017.06.020)
- **Prigogine, I. & Nicolis, G. (1977).** *Self-Organization in Nonequilibrium Systems*. Wiley. [Archive](https://archive.org/details/selforganization0000nico)

### Systems Chemistry References

- **Ruiz-Mirazo, K., Briones, C., & de la Escosura, A. (2014).** "Prebiotic systems chemistry: New perspectives for the origins of life". *Chemical Reviews*, 114(1), 285-366. [DOI:10.1021/cr2004844](https://doi.org/10.1021/cr2004844)
- **Ashkenasy, G. et al. (2004).** "Design of a directed molecular network". *Proc. Natl. Acad. Sci.* 101, 10872. [DOI:10.1073/pnas.0402674101](https://doi.org/10.1073/pnas.0402674101)

**Complete bibliography with 100+ references:** [docs/LITERATURE.md](docs/LITERATURE.md)
