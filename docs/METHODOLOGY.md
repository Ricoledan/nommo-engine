# Methodology

## Overview

Nommo Engine implements a scientifically rigorous molecular dynamics simulation framework specifically designed to study emergence phenomena in chemical systems. This document details the theoretical foundations, algorithms, and implementation choices that underpin the simulation.

## Theoretical Framework

### 1. Classical Molecular Dynamics

The simulation is based on classical mechanics, where particles represent molecular units (~10-100 atoms) evolving according to Newton's equations of motion:

```
F_i = m_i * a_i = m_i * d²r_i/dt²
```

Where:
- `F_i` is the force on particle i
- `m_i` is the mass of particle i  
- `r_i` is the position vector of particle i
- `t` is time

### 2. Force Field: Lennard-Jones Potential

Inter-particle interactions are modeled using the Lennard-Jones 12-6 potential:

```
V(r) = 4ε[(σ/r)¹² - (σ/r)⁶]
```

The corresponding force is:

```
F(r) = -dV/dr = 24ε/r[(2(σ/r)¹² - (σ/r)⁶)]
```

**Parameters:**
- `ε` (epsilon): Potential well depth (kJ/mol) - controls interaction strength
- `σ` (sigma): Zero-crossing distance (nm) - effective particle size
- `r_cutoff`: Force cutoff distance (typically 2.5σ) - computational efficiency

**Physical Justification:**
The Lennard-Jones potential captures the essential physics of inter-molecular interactions:
- Repulsive r⁻¹² term: Pauli exclusion principle at short range
- Attractive r⁻⁶ term: van der Waals forces at long range

### 3. Chemical Reactions

#### 3.1 Arrhenius Kinetics

Reaction rates follow the Arrhenius equation:

```
k(T) = A * exp(-E_a / RT)
```

Where:
- `k(T)`: Rate constant at temperature T
- `A`: Pre-exponential factor (frequency factor) ~ 10¹³ s⁻¹
- `E_a`: Activation energy barrier (kJ/mol)
- `R`: Gas constant (8.314 J/(mol·K))
- `T`: Temperature (K)

#### 3.2 Collision Theory

For a reaction to occur between particles i and j:

1. **Proximity criterion**: Distance < bonding_threshold
2. **Energy criterion**: Collision energy > activation energy
3. **Steric criterion**: Available bonding sites
4. **Thermodynamic criterion**: ΔG < 0 for spontaneous reactions

**Collision Energy Calculation:**
```
E_collision = ½μv_rel²
```

Where:
- `μ = (m_i * m_j)/(m_i + m_j)` is the reduced mass
- `v_rel = |v_i - v_j|` is the relative velocity magnitude

#### 3.3 Bond Thermodynamics

Bond formation/breaking follows Gibbs free energy:

```
ΔG = ΔH - TΔS
```

Where:
- `ΔH`: Enthalpy change (bond energy)
- `ΔS`: Entropy change
- `T`: Temperature

**Bond Breaking Probability:**
```
P_break = exp(-E_bond / k_B*T)
```

## Numerical Methods

### 1. Time Integration: Velocity Verlet Algorithm

The Velocity Verlet algorithm provides superior energy conservation compared to simple Euler methods:

**Algorithm:**
1. `r(t+Δt) = r(t) + v(t)Δt + ½a(t)Δt²`
2. Calculate new forces: `F(t+Δt)`
3. `v(t+Δt) = v(t) + ½[a(t) + a(t+Δt)]Δt`

**Advantages:**
- Symplectic (preserves phase space volume)
- Time-reversible
- Second-order accuracy in position and velocity
- Superior energy conservation

**Stability Criterion:**
Maximum stable timestep is limited by:
```
Δt_max ≈ 0.1 * sqrt(m/k_max)
```

Where `k_max` is the maximum force constant in the system.

### 2. Temperature Control: Berendsen Thermostat

Temperature is maintained using the Berendsen weak-coupling thermostat:

```
λ = sqrt(1 + (Δt/τ)(T_target/T_current - 1))
v_new = λ * v_old
```

**Parameters:**
- `τ`: Coupling time constant (typically 0.1-1.0 ps)
- `T_target`: Target temperature
- `T_current`: Instantaneous temperature from kinetic energy

**Temperature Calculation:**
```
T = 2⟨E_kinetic⟩ / (N_df * k_B)
```

Where `N_df = 3N` is the number of degrees of freedom for N particles in 3D.

### 3. Spatial Optimization

#### 3.1 Neighbor Lists

To reduce computational complexity from O(N²) to O(N), the simulation uses Verlet neighbor lists:

1. For each particle, maintain a list of neighbors within `r_cutoff + r_skin`
2. Update lists only when any particle moves more than `r_skin/2`
3. Force calculations only consider particles in neighbor lists

**Typical Parameters:**
- `r_skin = 0.3σ` (skin distance)
- Update frequency: every 10-20 timesteps

#### 3.2 Cell Lists

For very large systems (N > 10,000), cell-based spatial decomposition:

1. Divide simulation box into cells of size ≥ `r_cutoff`
2. Each particle belongs to one cell
3. Force calculations only between particles in neighboring cells

### 4. Boundary Conditions

#### 4.1 Periodic Boundaries

Particles leaving one side of the box re-enter from the opposite side:

```
if x > L_x: x = x - L_x
if x < 0:   x = x + L_x
```

**Minimum Image Convention:**
For distance calculations, use the shortest distance considering periodic images:

```
dx = x_j - x_i
if dx > L_x/2:  dx = dx - L_x
if dx < -L_x/2: dx = dx + L_x
```

#### 4.2 Reflective Boundaries

Particles bounce elastically off walls:
```
if x > L_x: 
    x = 2*L_x - x
    v_x = -v_x
```

## Emergence Detection Algorithms

### 1. Autocatalytic Set Detection

Based on Hordijk & Steel (2017) algorithm:

1. **Build reaction network graph**
   - Nodes: molecular species
   - Edges: reactions producing each species

2. **Find catalytic closure**
   - Iteratively add reactions that are catalyzed by existing species
   - Continue until no new reactions can be added

3. **Check self-sustaining criterion**
   - Every species in the set can be produced by reactions within the set
   - The set collectively catalyzes its own formation

### 2. Replicator Detection

**Criteria for identifying replicators:**

1. **Template matching**: Offspring structurally similar to parent
2. **Heredity**: Copying fidelity > threshold (typically 0.8)
3. **Multiplication**: Population growth over time
4. **Variation**: Mutations introduce diversity

**Algorithm:**
```
for each particle p:
    if p.parent_id exists:
        similarity = structural_similarity(p, parent)
        if similarity > threshold:
            mark as replicator
```

### 3. Complexity Metrics

#### 3.1 Shannon Entropy

Measures diversity of particle types:

```
H = -Σ p_i * log₂(p_i)
```

Where `p_i` is the fraction of particles of type i.

#### 3.2 Network Complexity

Based on bond network topology:

```
C_network = (E * log(E)) / (N * log(N))
```

Where:
- `E`: Number of bonds (edges)
- `N`: Number of particles (nodes)

#### 3.3 Cluster Size Distribution

Power-law distribution indicates criticality:

```
P(s) ∝ s^(-τ)
```

Where `s` is cluster size and `τ` is the critical exponent.

## Validation Procedures

### 1. Energy Conservation

Without thermostat, total energy should be conserved to within numerical precision:

```
|E(t) - E(0)| / E(0) < 10⁻⁶
```

### 2. Temperature Equilibration

With thermostat, temperature should converge to target:

```
|T(t) - T_target| / T_target < 0.01
```

### 3. Equation of State

For ideal gas (ε = 0), verify:

```
PV = Nk_BT
```

Where pressure is calculated from virial theorem:

```
P = ρk_BT + (1/3V)⟨Σr_ij · F_ij⟩
```

### 4. Diffusion Coefficient

For simple liquids, Einstein relation:

```
D = lim(t→∞) ⟨|r(t) - r(0)|²⟩ / (6t)
```

Should match experimental/theoretical values for comparable systems.

## Computational Complexity

### Time Complexity

- **Force calculation**: O(N) with neighbor lists, O(N²) without
- **Integration**: O(N)
- **Thermostat**: O(N)
- **Bond detection**: O(N_neighbors)
- **Overall per timestep**: O(N)

### Space Complexity

- **Particle storage**: O(N)
- **Neighbor lists**: O(N * N_neighbors)
- **Spatial indexing**: O(N)
- **Overall**: O(N)

### Scaling Performance

**Target performance (single CPU core):**
- 1,000 particles: ~100 timesteps/second
- 10,000 particles: ~10 timesteps/second
- 100,000 particles: ~1 timestep/second

## Implementation Considerations

### 1. Numerical Precision

- Use double precision (64-bit) for positions and forces
- Single precision (32-bit) acceptable for velocities in some cases
- Monitor energy drift as quality indicator

### 2. Parameter Selection

**Typical values for molecular systems:**
- Timestep: 0.001-0.01 ps
- Temperature: 100-1000 K
- Density: 0.1-1.0 particles/nm³
- LJ epsilon: 0.1-10 kJ/mol
- LJ sigma: 0.2-0.5 nm

### 3. Stability Considerations

- Monitor maximum force magnitude
- Adjust timestep if forces become too large
- Use velocity rescaling to prevent "explosion"
- Check for overlapping particles (r < 0.1σ)

### 4. Memory Management

- Pre-allocate arrays for better performance
- Use memory pools for particle creation/destruction
- Consider structure-of-arrays vs array-of-structures

## References

1. Frenkel, D. & Smit, B. (2001). *Understanding Molecular Simulation*. Academic Press.
2. Allen, M. P. & Tildesley, D. J. (2017). *Computer Simulation of Liquids*. Oxford University Press.
3. Berendsen, H. J. C. et al. (1984). "Molecular dynamics with coupling to an external bath". *J. Chem. Phys.* 81, 3684.
4. Verlet, L. (1967). "Computer experiments on classical fluids". *Phys. Rev.* 159, 98.
5. Hordijk, W. & Steel, M. (2017). "Detecting autocatalytic, self-sustaining sets in chemical reaction systems". *J. Theor. Biol.* 227, 451.
6. Kauffman, S. A. (1986). "Autocatalytic sets of proteins". *J. Theor. Biol.* 119, 1.