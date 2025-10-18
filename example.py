#!/usr/bin/env python
"""
Example script demonstrating Nommo Engine capabilities.

This script creates a simple universe, runs a simulation,
and demonstrates the key features of the engine.
"""

import numpy as np
from pathlib import Path

from nommo.core.types import UniverseParameters, ParticleTypeConfig
from nommo.core.universe import Universe
from nommo.storage.persistence import HDF5Storage


def create_test_universe():
    """Create a simple test universe."""
    params = UniverseParameters(
        name="test-universe",
        description="Example simulation demonstrating Nommo Engine",
        particle_types={
            "type_A": ParticleTypeConfig(
                name="type_A",
                mass=50.0,
                radius=0.2,
                max_bonds=2,
                lj_epsilon=1.0,
                lj_sigma=0.3
            ),
            "type_B": ParticleTypeConfig(
                name="type_B",
                mass=75.0,
                radius=0.25,
                max_bonds=3,
                lj_epsilon=1.2,
                lj_sigma=0.35
            ),
        },
        initial_composition={
            "type_A": 50,
            "type_B": 50,
        }
    )
    
    # Adjust some physics parameters
    params.physics.timestep = 0.001
    params.physics.cutoff_distance = 1.0
    params.thermodynamics.temperature = 300.0
    params.box.dimensions = (10.0, 10.0, 10.0)
    
    return Universe(params)


def main():
    """Run example simulation."""
    print("=" * 60)
    print("Nommo Engine - Example Simulation")
    print("=" * 60)
    
    # Create universe
    print("\nCreating universe...")
    universe = create_test_universe()
    print(f"  Universe ID: {universe.id}")
    print(f"  Particles: {len(universe.particles)}")
    print(f"  Box size: {universe.box_size}")
    print(f"  Temperature: {universe.params.thermodynamics.temperature} K")
    
    # Run simulation
    print("\nRunning simulation...")
    n_steps = 100
    
    for i in range(n_steps):
        metrics = universe.step()
        
        if i % 20 == 0:
            print(f"  Step {i:4d}: "
                  f"T={metrics.current_temperature:.1f}K, "
                  f"E={metrics.total_energy:.2f} kJ/mol, "
                  f"Particles={metrics.total_particles}")
    
    print(f"\nSimulation complete!")
    print(f"  Final time: {universe.time:.3f} ps")
    print(f"  Final temperature: {metrics.current_temperature:.1f} K")
    print(f"  Total energy: {metrics.total_energy:.2f} kJ/mol")
    print(f"  Kinetic energy: {metrics.total_kinetic_energy:.2f} kJ/mol")
    print(f"  Potential energy: {metrics.total_potential_energy:.2f} kJ/mol")
    print(f"  Shannon entropy: {metrics.shannon_entropy:.3f}")
    
    # Save universe
    print("\nSaving universe...")
    storage = HDF5Storage(Path("data/universes"))
    filepath = storage.save_universe(universe, "example_universe.h5")
    print(f"  Saved to: {filepath}")
    
    # Save trajectory
    if universe.history:
        print("\nSaving trajectory...")
        trajectory_path = storage.save_trajectory(
            universe.id, 
            list(universe.history),
            "example_trajectory.h5"
        )
        print(f"  Trajectory saved to: {trajectory_path}")
    
    print("\n" + "=" * 60)
    print("Example complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()