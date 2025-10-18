# Parameter Guide and Sensitivity Analysis

## Overview

This document provides comprehensive guidance on parameter selection for Nommo Engine simulations, including valid ranges, sensitivity analysis, and recommendations for different research scenarios.

## Core Physics Parameters

### 1. Time Integration

#### 1.1 Timestep (dt)

**Parameter:** `timestep` in UniverseParameters
**Units:** picoseconds (ps)
**Valid Range:** 0.0001 - 0.01 ps

**Guidelines:**
```
dt < 0.1 * sqrt(m_min / k_max)
```

Where:
- `m_min` is the lightest particle mass
- `k_max` is the steepest force constant

**Recommended Values:**
- **Conservative:** 0.001 ps (most stable)
- **Standard:** 0.002 ps (good balance)
- **Aggressive:** 0.005 ps (faster, check energy drift)

**Sensitivity Analysis:**

| Timestep (ps) | Energy Drift (%/ns) | Performance | Stability |
|---------------|---------------------|-------------|-----------|
| 0.0005        | < 0.001             | 0.5x        | Excellent |
| 0.001         | < 0.01              | 1.0x        | Very Good |
| 0.002         | < 0.05              | 2.0x        | Good      |
| 0.005         | 0.1-0.5             | 5.0x        | Marginal  |
| 0.01          | > 1.0               | 10.0x       | Unstable  |

**Impact on Results:**
- **Too small:** Unnecessary computational cost
- **Too large:** Energy drift, artificial heating, simulation instability

---

#### 1.2 Integration Algorithm

**Parameter:** Currently Velocity Verlet (fixed)
**Alternatives:** Leap-frog, Runge-Kutta (future implementations)

**Velocity Verlet Advantages:**
- Symplectic (conserves energy)
- Time-reversible
- Second-order accuracy
- Industry standard for MD

---

### 2. Force Field Parameters

#### 2.1 Lennard-Jones Epsilon (ε)

**Parameter:** `lj_epsilon` in UniverseParameters
**Units:** kJ/mol
**Valid Range:** 0.1 - 50 kJ/mol

**Physical Meaning:** Depth of potential well (interaction strength)

**Typical Values:**
- **Noble gases:** 0.1-1 kJ/mol
- **Small molecules:** 1-5 kJ/mol
- **Organic compounds:** 5-15 kJ/mol
- **Strong interactions:** 15-50 kJ/mol

**Sensitivity Analysis:**

| ε (kJ/mol) | Phase Behavior | Reaction Rate | Complexity |
|------------|----------------|---------------|------------|
| 0.5        | Gas-like       | Fast          | Low        |
| 2.0        | Liquid-like    | Moderate      | Medium     |
| 8.0        | Dense liquid   | Slow          | High       |
| 20.0       | Solid-like     | Very slow     | Very high  |

**Emergence Dependence:**
- **Low ε (< 2):** Fast dynamics, low stability
- **Medium ε (2-8):** Optimal for emergence studies
- **High ε (> 15):** Slow but stable structures

---

#### 2.2 Lennard-Jones Sigma (σ)

**Parameter:** `lj_sigma` in UniverseParameters
**Units:** nanometers (nm)
**Valid Range:** 0.2 - 1.0 nm

**Physical Meaning:** Effective particle diameter

**Typical Values:**
- **Small atoms:** 0.2-0.3 nm
- **Small molecules:** 0.3-0.4 nm
- **Medium molecules:** 0.4-0.6 nm
- **Large molecules:** 0.6-1.0 nm

**Sensitivity Analysis:**

| σ (nm) | Density | Diffusion | Packing |
|--------|---------|-----------|---------|
| 0.25   | High    | Fast      | Tight   |
| 0.35   | Medium  | Medium    | Normal  |
| 0.50   | Low     | Slow      | Loose   |
| 0.70   | Very low| Very slow | Sparse  |

**Impact on Chemistry:**
- **Small σ:** Higher collision frequency, faster reactions
- **Large σ:** Lower collision frequency, slower reactions
- **Ratio σ/box_size:** Should be < 0.1 for realistic behavior

---

#### 2.3 Cutoff Distance

**Parameter:** `cutoff_distance` in UniverseParameters
**Units:** nanometers (nm)
**Valid Range:** 2.0σ - 5.0σ

**Recommended:** 2.5σ (standard choice)

**Trade-offs:**

| Cutoff/σ | Accuracy | Performance | Energy Drift |
|----------|----------|-------------|--------------|
| 2.0      | Poor     | Fast        | High         |
| 2.5      | Good     | Medium      | Low          |
| 3.0      | Better   | Slow        | Very low     |
| 4.0      | Excellent| Very slow   | Minimal      |

**Guidelines:**
- Minimum: 2.0σ (acceptable for exploratory studies)
- Standard: 2.5σ (good balance)
- High accuracy: 3.0σ (for quantitative studies)

---

### 3. Thermodynamic Parameters

#### 3.1 Temperature

**Parameter:** `temperature` in UniverseParameters
**Units:** Kelvin (K)
**Valid Range:** 50 - 2000 K

**Physical Ranges:**
- **Cryogenic:** 50-100 K (very slow dynamics)
- **Room temperature:** 250-350 K (moderate dynamics)
- **Elevated:** 400-600 K (fast dynamics)
- **High temperature:** 600-1000 K (very fast, gas-like)
- **Extreme:** 1000+ K (plasma-like, use with caution)

**Sensitivity Analysis:**

| Temperature (K) | Reaction Rate | Stability | Emergence Time |
|-----------------|---------------|-----------|----------------|
| 150             | Very slow     | High      | Very long      |
| 300             | Moderate      | Good      | Long           |
| 500             | Fast          | Medium    | Medium         |
| 800             | Very fast     | Low       | Short          |
| 1200            | Extreme       | Very low  | Very short     |

**Temperature Effects:**
```
Rate increase ∝ exp(-Ea/RT)
```

**Rule of thumb:** Every 100K increase approximately doubles reaction rates.

---

#### 3.2 Thermostat Parameters

**Parameter:** `thermostat_coupling` in UniverseParameters
**Units:** picoseconds (ps)
**Valid Range:** 0.01 - 10 ps

**Physical Meaning:** Time constant for temperature coupling

**Recommended Values:**
- **Tight coupling:** 0.1 ps (quick equilibration)
- **Standard:** 0.5 ps (good balance)
- **Loose coupling:** 2.0 ps (more natural fluctuations)
- **Very loose:** 5.0+ ps (minimal intervention)

**Sensitivity Analysis:**

| τ (ps) | Equilibration | Fluctuations | Artifacts |
|--------|---------------|--------------|-----------|
| 0.05   | Very fast     | Suppressed   | High      |
| 0.1    | Fast          | Reduced      | Medium    |
| 0.5    | Medium        | Natural      | Low       |
| 2.0    | Slow          | Enhanced     | Minimal   |
| 10.0   | Very slow     | Large        | None      |

---

## Chemical Parameters

### 4. Reaction Kinetics

#### 4.1 Activation Energy (Ea)

**Parameter:** `activation_energy` in UniverseParameters
**Units:** kJ/mol
**Valid Range:** 0 - 200 kJ/mol

**Physical Ranges:**
- **Barrierless:** 0-5 kJ/mol (diffusion-limited)
- **Low barrier:** 5-20 kJ/mol (fast at room temp)
- **Medium barrier:** 20-50 kJ/mol (moderate rates)
- **High barrier:** 50-100 kJ/mol (slow rates)
- **Very high:** 100+ kJ/mol (very slow, high temp needed)

**Temperature Dependence:**
```
Rate ∝ exp(-Ea/RT)
```

**Sensitivity Analysis (at 300K):**

| Ea (kJ/mol) | Relative Rate | Half-life | Practical |
|-------------|---------------|-----------|-----------|
| 0           | 1.0           | Instant   | Diffusion |
| 10          | 0.018         | ~1 ps     | Very fast |
| 25          | 4.5×10⁻⁵      | ~100 ps   | Fast      |
| 50          | 5.3×10⁻⁹      | ~1 ns     | Moderate  |
| 75          | 6.1×10⁻¹³     | ~1 μs     | Slow      |
| 100         | 7.0×10⁻¹⁷     | ~1 ms     | Very slow |

---

#### 4.2 Bond Energy

**Parameter:** `bond_energy` in UniverseParameters
**Units:** kJ/mol
**Valid Range:** 5 - 500 kJ/mol

**Physical Ranges:**
- **Weak bonds:** 5-20 kJ/mol (hydrogen bonds)
- **Medium bonds:** 20-100 kJ/mol (van der Waals)
- **Strong bonds:** 100-300 kJ/mol (covalent single)
- **Very strong:** 300+ kJ/mol (multiple bonds)

**Stability Analysis:**

| Bond Energy (kJ/mol) | Lifetime at 300K | Stability | Use Case |
|----------------------|------------------|-----------|----------|
| 10                   | 0.1 ps           | Very low  | Dynamic  |
| 25                   | 100 ps           | Low       | Temporary|
| 50                   | 1 μs             | Medium    | Moderate |
| 100                  | 1 s              | High      | Stable   |
| 200                  | Years            | Very high | Permanent|

---

#### 4.3 Reaction Probability

**Parameter:** `reaction_probability` in UniverseParameters
**Units:** dimensionless (0-1)
**Valid Range:** 0.001 - 1.0

**Physical Meaning:** Steric factor (geometry/orientation effects)

**Typical Values:**
- **Gas phase:** 0.1-1.0 (unhindered)
- **Solution:** 0.01-0.1 (some hindrance)
- **Solid surface:** 0.001-0.01 (high hindrance)
- **Enzyme active site:** 0.1-1.0 (optimized geometry)

**Impact Analysis:**

| Probability | Effective Rate | Character | Use Case |
|-------------|----------------|-----------|----------|
| 1.0         | Maximum        | Ideal     | Gas phase|
| 0.1         | 10% of max     | Realistic | Solution |
| 0.01        | 1% of max      | Hindered  | Crowded  |
| 0.001       | 0.1% of max    | Rare      | Specific |

---

## System Parameters

### 5. Box Geometry

#### 5.1 Box Dimensions

**Parameter:** `dimensions` in UniverseParameters
**Units:** nanometers (nm)
**Valid Range:** 5-1000 nm per dimension

**Size Guidelines:**
```
L > 2 * cutoff_distance (minimum)
L > 10 * σ (recommended minimum)
```

**System Size Effects:**

| Box Size (nm) | Particles | Regime     | Finite Size Effects |
|---------------|-----------|------------|-------------------|
| 5×5×5         | 100-500   | Nano       | Large             |
| 10×10×10      | 500-2000  | Small      | Medium            |
| 20×20×20      | 2000-8000 | Medium     | Small             |
| 50×50×50      | 50000+    | Large      | Negligible        |

**Aspect Ratio Considerations:**
- **Cubic:** Isotropic, no preferred direction
- **Rectangular:** Can study anisotropic effects
- **Slab geometry:** For surface/interface studies

---

#### 5.2 Boundary Conditions

**Parameter:** `boundary_type` in UniverseParameters
**Options:** periodic, reflective, absorbing

**Periodic Boundaries:**
- **Pros:** No wall effects, infinite-like system
- **Cons:** Artificial correlations across box
- **Use:** Bulk properties, homogeneous systems

**Reflective Boundaries:**
- **Pros:** Confined system, surface effects
- **Cons:** Wall artifacts, energy buildup
- **Use:** Droplets, clusters, confined systems

**Absorbing Boundaries:**
- **Pros:** Open system, no back-scattering
- **Cons:** Particle loss, non-equilibrium
- **Use:** Flow systems, evaporation studies

---

### 6. Particle Composition

#### 6.1 Particle Density

**Parameter:** `initial_composition` counts / box volume
**Units:** particles/nm³
**Valid Range:** 0.01 - 2.0 particles/nm³

**Density Regimes:**

| Density (particles/nm³) | Phase    | Behavior           | Emergence |
|-------------------------|----------|--------------------|-----------|
| 0.01-0.1               | Gas      | Rare collisions    | Unlikely  |
| 0.1-0.5                | Dilute   | Frequent contact   | Possible  |
| 0.5-1.0                | Liquid   | Continuous contact | Likely    |
| 1.0-2.0                | Dense    | Crowded           | Complex   |

**Optimal Ranges:**
- **Emergence studies:** 0.3-0.8 particles/nm³
- **Kinetics studies:** 0.1-0.5 particles/nm³
- **Structure studies:** 0.5-1.0 particles/nm³

---

#### 6.2 Type Ratios

**Parameter:** Relative counts in `initial_composition`

**Strategy Guidelines:**

**Equal Ratios (1:1:1...):**
- **Pros:** Unbiased evolution, maximum diversity
- **Cons:** May not reach critical concentrations
- **Use:** Exploratory studies, evolution experiments

**Stoichiometric Ratios:**
- **Pros:** Matched to specific chemistry
- **Cons:** Predetermined outcomes
- **Use:** Targeted reaction studies

**Biased Ratios (e.g., 10:1:1):**
- **Pros:** Excess reactants, drives reactions
- **Cons:** May overwhelm rare species
- **Use:** Catalysis studies, autocatalysis

---

## Emergence-Specific Parameters

### 7. Detection Thresholds

#### 7.1 Replicator Detection

**Template Similarity Threshold:** 0.7-0.95
- **Low (0.7):** Detects loose copies (high false positives)
- **Medium (0.8):** Balanced detection
- **High (0.95):** Only exact copies (high false negatives)

**Minimum Lineage Depth:** 3-10 generations
- **Low (3):** Early detection, more noise
- **High (10):** Late but confident detection

#### 7.2 Autocatalytic Set Criteria

**Closure Threshold:** 0.8-1.0
- **0.8:** Partial closure acceptable
- **1.0:** Perfect closure required

**Minimum Set Size:** 2-5 species
- **Small (2):** Simple autocatalysis
- **Large (5+):** Complex networks

---

## Parameter Sensitivity Rankings

### High Sensitivity (Major Impact)
1. **Temperature** - Affects all rates exponentially
2. **Activation Energy** - Controls reaction feasibility
3. **Timestep** - Affects stability and accuracy
4. **Epsilon (ε)** - Controls phase behavior

### Medium Sensitivity (Moderate Impact)
5. **Bond Energy** - Affects structure stability
6. **Density** - Controls collision frequency
7. **Box Size** - Finite size effects
8. **Thermostat Coupling** - Temperature fluctuations

### Low Sensitivity (Minor Impact)
9. **Sigma (σ)** - Within reasonable range
10. **Cutoff Distance** - Above 2.5σ
11. **Reaction Probability** - Linear scaling
12. **Boundary Type** - For bulk properties

---

## Recommended Parameter Sets

### 1. Default Exploration Set
```yaml
temperature: 300.0          # Room temperature
timestep: 0.002            # Safe timestep
lj_epsilon: 2.0            # Moderate interactions
lj_sigma: 0.35             # Typical molecular size
cutoff_distance: 0.875     # 2.5 * sigma
activation_energy: 25.0    # Moderate barrier
bond_energy: 50.0          # Stable but breakable
density: 0.5               # Liquid-like
box_size: [10, 10, 10]     # Medium system
```

### 2. Fast Exploration Set
```yaml
temperature: 500.0          # Higher temperature
timestep: 0.005            # Larger timestep
lj_epsilon: 1.0            # Weaker interactions
activation_energy: 15.0    # Lower barrier
bond_energy: 30.0          # Weaker bonds
density: 0.3               # Lower density
```

### 3. Detailed Study Set
```yaml
temperature: 300.0          # Standard
timestep: 0.001            # Conservative
lj_epsilon: 3.0            # Strong interactions
cutoff_distance: 1.05      # 3.0 * sigma
activation_energy: 30.0    # Realistic barrier
bond_energy: 75.0          # Stable structures
density: 0.7               # Dense liquid
box_size: [15, 15, 15]     # Larger system
```

### 4. High-Performance Set
```yaml
temperature: 400.0          # Faster dynamics
timestep: 0.003            # Balanced
lj_epsilon: 1.5            # Moderate
cutoff_distance: 0.875     # Standard
activation_energy: 20.0    # Lower barrier
density: 0.4               # Medium density
box_size: [12, 12, 12]     # Reasonable size
```

---

## Parameter Validation Protocol

### 1. Stability Checks
- Run 10,000 timesteps
- Monitor energy drift < 0.1%
- Check temperature stability < 5%
- Verify no particle overlaps

### 2. Physical Reasonableness
- Diffusion coefficients in realistic range
- Phase behavior matches expectations
- Reaction rates temperature-dependent
- Equilibrium constants reasonable

### 3. Emergence Feasibility
- Sufficient particle density
- Reasonable reaction rates
- Stable but dynamic structures
- Observable timescales

### 4. Computational Efficiency
- Timestep stability margin
- Performance benchmarks
- Memory usage acceptable
- Scalability to target size

---

## Common Parameter Problems

### Problem: Simulation Explodes
**Symptoms:** Huge forces, temperature spikes, particles overlap
**Causes:** Timestep too large, particles too close initially
**Solutions:** Reduce timestep, add soft potential, equilibrate slowly

### Problem: Nothing Happens
**Symptoms:** No reactions, static structure, low energy
**Causes:** Temperature too low, activation energy too high, density too low
**Solutions:** Increase temperature, reduce barriers, increase density

### Problem: Too Fast Evolution
**Symptoms:** Instant reactions, immediate complexity, unrealistic rates
**Causes:** No activation barriers, temperature too high, density too high
**Solutions:** Add realistic barriers, reduce temperature, lower density

### Problem: Energy Drift
**Symptoms:** Total energy increasing/decreasing over time
**Causes:** Timestep too large, cutoff too small, integration errors
**Solutions:** Reduce timestep, increase cutoff, check implementation

### Problem: Temperature Fluctuations
**Symptoms:** Large temperature swings, poor control
**Causes:** Thermostat coupling wrong, system too small, poor statistics
**Solutions:** Adjust coupling time, increase system size, longer averaging

---

## References

1. Allen, M. P. & Tildesley, D. J. (2017). *Computer Simulation of Liquids*. Oxford University Press.
2. Frenkel, D. & Smit, B. (2001). *Understanding Molecular Simulation*. Academic Press.
3. Leach, A. R. (2001). *Molecular Modelling: Principles and Applications*. Pearson.
4. Tuckerman, M. E. (2010). *Statistical Mechanics: Theory and Molecular Simulation*. Oxford University Press.
5. Rapaport, D. C. (2004). *The Art of Molecular Dynamics Simulation*. Cambridge University Press.