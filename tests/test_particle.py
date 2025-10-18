"""
Tests for Particle class and related functionality.
"""

import pytest
import numpy as np
from hypothesis import given, strategies as st

from nommo.core.particle import Particle, ParticlePool


class TestParticle:
    """Test Particle class."""
    
    def test_particle_creation(self):
        """Test basic particle creation."""
        p = Particle(
            mass=50.0,
            charge=1.0,
            radius=0.2,
            position=np.array([1.0, 2.0, 3.0]),
            velocity=np.array([0.1, 0.2, 0.3])
        )
        
        assert p.mass == 50.0
        assert p.charge == 1.0
        assert p.radius == 0.2
        assert np.allclose(p.position, [1.0, 2.0, 3.0])
        assert np.allclose(p.velocity, [0.1, 0.2, 0.3])
        assert len(p.id) > 0
    
    def test_kinetic_energy(self):
        """Test kinetic energy calculation."""
        p = Particle(
            mass=50.0,
            velocity=np.array([1.0, 0.0, 0.0])
        )
        
        expected_ke = 0.5 * 50.0 * 1.0**2 * 0.01
        assert np.isclose(p.kinetic_energy(), expected_ke)
    
    def test_kinetic_energy_caching(self):
        """Test that kinetic energy is cached properly."""
        p = Particle(
            mass=50.0,
            velocity=np.array([1.0, 0.0, 0.0])
        )
        
        ke1 = p.kinetic_energy()
        assert p._cache_valid
        
        ke2 = p.kinetic_energy()
        assert ke1 == ke2
        
        p.update_velocity(np.array([2.0, 0.0, 0.0]))
        assert not p._cache_valid
        
        ke3 = p.kinetic_energy()
        assert ke3 != ke1
        assert p._cache_valid
    
    def test_distance_calculations(self):
        """Test distance calculation methods."""
        p1 = Particle(position=np.array([0.0, 0.0, 0.0]))
        p2 = Particle(position=np.array([3.0, 4.0, 0.0]))
        
        assert p1.distance_to(p2) == 5.0
        assert p1.distance_squared_to(p2) == 25.0
        assert np.allclose(p1.distance_vector_to(p2), [3.0, 4.0, 0.0])
    
    def test_bonding(self):
        """Test bond formation and breaking."""
        p1 = Particle(max_bonds=2)
        p2 = Particle(max_bonds=2)
        p3 = Particle(max_bonds=1)
        
        assert p1.can_bond_with(p2)
        assert p1.add_bond(p2)
        assert p2.id in p1.bonds
        assert p1.id in p2.bonds
        
        assert not p1.add_bond(p2)
        
        assert p1.can_bond_with(p3)
        assert p1.add_bond(p3)
        
        p4 = Particle(max_bonds=2)
        assert not p1.can_bond_with(p4)
        
        assert p1.remove_bond(p2.id)
        assert p2.id not in p1.bonds
        
        count = p1.clear_bonds()
        assert count == 1
        assert len(p1.bonds) == 0
    
    def test_boundary_conditions(self):
        """Test periodic and reflective boundary conditions."""
        box_size = np.array([10.0, 10.0, 10.0])
        
        p = Particle(
            position=np.array([11.0, -1.0, 5.0]),
            velocity=np.array([1.0, -1.0, 0.0])
        )
        
        p.apply_periodic_boundary(box_size)
        assert np.allclose(p.position, [1.0, 9.0, 5.0])
        
        p2 = Particle(
            position=np.array([11.0, -1.0, 5.0]),
            velocity=np.array([1.0, -1.0, 0.0])
        )
        
        p2.apply_reflective_boundary(box_size)
        assert p2.position[0] < box_size[0]
        assert p2.position[1] >= 0
        assert p2.velocity[0] < 0
        assert p2.velocity[1] > 0
    
    def test_serialization(self):
        """Test particle serialization."""
        p1 = Particle(
            mass=50.0,
            position=np.array([1.0, 2.0, 3.0]),
            velocity=np.array([0.1, 0.2, 0.3]),
            particle_type="test",
            generation=5
        )
        p1.bonds.add("bond1")
        p1.bonds.add("bond2")
        
        data = p1.to_dict()
        assert data["mass"] == 50.0
        assert data["position"] == [1.0, 2.0, 3.0]
        assert "bond1" in data["bonds"]
        assert data["generation"] == 5
        
        p2 = Particle.from_dict(data)
        assert p2.mass == p1.mass
        assert np.allclose(p2.position, p1.position)
        assert p2.bonds == p1.bonds
        assert p2.generation == p1.generation
    
    def test_particle_copy(self):
        """Test particle copying."""
        p1 = Particle(
            mass=50.0,
            position=np.array([1.0, 2.0, 3.0]),
            generation=2
        )
        p1.bonds.add("bond1")
        
        p2 = p1.copy(new_id=True)
        assert p2.id != p1.id
        assert p2.mass == p1.mass
        assert np.allclose(p2.position, p1.position)
        assert p2.generation == 3
        assert p2.parent_id == p1.id
        assert "bond1" in p2.bonds
        
        p1.position[0] = 5.0
        assert p2.position[0] == 1.0
    
    @given(
        mass=st.floats(min_value=0.1, max_value=1000.0),
        velocity=st.tuples(
            st.floats(min_value=-10, max_value=10),
            st.floats(min_value=-10, max_value=10),
            st.floats(min_value=-10, max_value=10)
        )
    )
    def test_kinetic_energy_property(self, mass, velocity):
        """Property test: kinetic energy is always non-negative."""
        p = Particle(
            mass=mass,
            velocity=np.array(velocity)
        )
        assert p.kinetic_energy() >= 0


class TestParticlePool:
    """Test ParticlePool memory management."""
    
    def test_pool_creation(self):
        """Test pool initialization."""
        pool = ParticlePool(initial_size=10)
        assert pool.pool_size == 10
        assert pool.active_count == 0
    
    def test_acquire_release(self):
        """Test particle acquisition and release."""
        pool = ParticlePool(initial_size=5)
        
        p1 = pool.acquire(mass=50.0, charge=1.0)
        assert p1.mass == 50.0
        assert p1.charge == 1.0
        assert pool.pool_size == 4
        assert pool.active_count == 1
        
        p2 = pool.acquire()
        assert pool.pool_size == 3
        assert pool.active_count == 2
        
        pool.release(p1)
        assert pool.pool_size == 4
        assert pool.active_count == 1
        
        p3 = pool.acquire(mass=75.0)
        assert p3.mass == 75.0
    
    def test_pool_exhaustion(self):
        """Test behavior when pool is exhausted."""
        pool = ParticlePool(initial_size=2)
        
        particles = []
        for i in range(5):
            particles.append(pool.acquire())
        
        assert pool.pool_size == 0
        assert pool.active_count == 5
        
        pool.release_all(particles)
        assert pool.pool_size == 5
        assert pool.active_count == 0