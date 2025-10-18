# Results Directory

This directory contains scientific results, validation data, and example outputs from Nommo Engine simulations.

## Directory Structure

### `/validation/`
Validation test results demonstrating scientific accuracy:
- `energy_conservation.png` - Energy drift over time
- `temperature_control.png` - Thermostat performance
- `equation_of_state.png` - PV = NkT validation
- `phase_transition.png` - Gas-liquid transition
- `arrhenius_kinetics.png` - Reaction rate vs temperature
- `equilibrium_constants.png` - Thermodynamic equilibrium
- `diffusion_coefficient.png` - Transport properties
- `viscosity.png` - Fluid properties
- `autocatalytic_sets.png` - Emergence validation
- `replicator_evolution.png` - Self-replication dynamics
- `performance_scaling.png` - Computational scaling

### `/benchmarks/`
Performance benchmarks and scaling studies:
- `cpu_performance.csv` - Timestep performance vs system size
- `memory_usage.csv` - Memory scaling measurements
- `algorithm_comparison.csv` - Different algorithm performance
- `optimization_impact.csv` - Effect of various optimizations

### `/examples/`
Example simulation outputs demonstrating key features:
- `simple_gas/` - Basic molecular dynamics
- `liquid_formation/` - Phase transition example
- `chemical_reactions/` - Bond formation/breaking
- `autocatalysis/` - Self-catalyzing networks
- `replication/` - Template-based copying
- `emergence_timeline/` - Complexity evolution

### `/sensitivity/`
Parameter sensitivity analysis results:
- `temperature_sweep/` - Temperature vs emergence time
- `density_effects/` - Particle density impact
- `activation_energy/` - Barrier height effects
- `bond_strength/` - Stability vs dynamics
- `system_size/` - Finite size effects
- `timestep_convergence/` - Integration accuracy

### `/emergence/`
Detailed emergence phenomenon studies:
- `autocatalytic_networks/` - Network formation dynamics
- `replicator_competition/` - Evolution and selection
- `information_flow/` - Information propagation
- `complexity_metrics/` - Various complexity measures
- `phase_transitions/` - Emergence phase diagrams

## File Formats

### Data Files
- **CSV format** for tabular data (time series, parameter sweeps)
- **HDF5 format** for large datasets (trajectories, snapshots)
- **JSON format** for metadata and configuration summaries

### Plots
- **PNG format** for figures (300 DPI for publications)
- **SVG format** for vector graphics when needed
- **PDF format** for multi-panel figures

### Analysis Scripts
- **Python scripts** for generating plots and analysis
- **Jupyter notebooks** for interactive exploration
- **Configuration files** for reproducing results

## Usage Examples

### Validation Results
```bash
# View all validation plots
ls results/validation/*.png

# Check energy conservation data
cat results/validation/energy_conservation.csv

# Reproduce validation plots
python scripts/generate_validation_plots.py
```

### Parameter Sensitivity
```bash
# Temperature sweep results
ls results/sensitivity/temperature_sweep/

# Plot emergence time vs temperature
python scripts/plot_temperature_sensitivity.py
```

### Example Simulations
```bash
# Run simple gas example
nommo universe create simple-gas --preset minimal
nommo run start simple-gas --ticks 10000
nommo analyze plot simple-gas temperature

# Run autocatalysis example
nommo universe create autocatal --preset autocatalytic
nommo run start autocatal --ticks 50000 --watch
nommo analyze emergence autocatal
```

## Reproducing Results

All results can be reproduced using the following workflow:

1. **Install dependencies**:
   ```bash
   poetry install
   ```

2. **Run validation suite**:
   ```bash
   nommo validate --all --output results/validation/
   ```

3. **Generate benchmarks**:
   ```bash
   python scripts/run_benchmarks.py --output results/benchmarks/
   ```

4. **Create example outputs**:
   ```bash
   python scripts/generate_examples.py --output results/examples/
   ```

5. **Perform sensitivity analysis**:
   ```bash
   python scripts/sensitivity_analysis.py --output results/sensitivity/
   ```

## Publication-Ready Figures

Selected figures suitable for publication:

- `validation/energy_conservation.png` - Demonstrates simulation stability
- `validation/arrhenius_kinetics.png` - Validates chemical kinetics
- `emergence/autocatalytic_networks.png` - Shows network formation
- `emergence/replicator_evolution.png` - Demonstrates self-replication
- `benchmarks/performance_scaling.png` - Computational efficiency

## Data Provenance

Each result includes metadata documenting:
- Software version used
- Parameter values
- Random seed for reproducibility
- Compute environment details
- Generation timestamp

Example metadata (embedded in HDF5 files):
```
{
  "software": "nommo-engine",
  "version": "0.1.0",
  "git_commit": "a1b2c3d4",
  "parameters": {...},
  "random_seed": 12345,
  "timestamp": "2024-01-15T10:30:00Z",
  "environment": {
    "platform": "Linux-5.4.0",
    "python": "3.11.2",
    "cpu": "Intel Core i7-10700K"
  }
}
```

## Contributing Results

When adding new results:

1. **Use descriptive names** that indicate the phenomenon studied
2. **Include metadata** for reproducibility
3. **Add documentation** explaining the significance
4. **Use consistent formats** matching existing files
5. **Update this README** with new categories as needed

## Archive Policy

- **Current results** (< 6 months): Keep in main directory
- **Historical results**: Archive to `archive/` subdirectory
- **Large datasets** (> 1GB): Store compressed with documentation
- **Failed experiments**: Keep in `failed/` with notes on issues

This directory serves as both a demonstration of Nommo Engine capabilities and a scientific record of validation studies.