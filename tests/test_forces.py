"""
Tests for force calculations and physics.
"""

import pytest
import numpy as np
from hypothesis import given, strategies as st

from nommo.physics.forces import (
    lennard_jones_force, LennardJonesForce,
    CoulombicForce, CompositeForce
)
from nommo.physics.constants import MIN_DISTANCE


class TestLennardJonesForce:
    """Test Lennard-Jones force calculations."""
    
    def test_lj_force_at_sigma(self):
        """At r = σ, force should be zero and potential should be zero."""
        force_mag, potential = lennard_jones_force(
            r=1.0, epsilon=1.0, sigma=1.0
        )
        
        assert np.isclose(potential, 0.0, atol=1e-10)
    
    def test_lj_force_at_minimum(self):
        """At r = 2^(1/6) * σ, potential should be minimum."""
        r_min = 2**(1/6)
        force_mag, potential = lennard_jones_force(
            r=r_min, epsilon=1.0, sigma=1.0
        )
        
        assert np.isclose(force_mag, 0.0, atol=1e-10)
        assert np.isclose(potential, -1.0, atol=1e-10)
    
    def test_lj_force_attractive(self):
        """At r > σ, force should be attractive (negative)."""
        force_mag, potential = lennard_jones_force(
            r=1.5, epsilon=1.0, sigma=1.0
        )
        
        assert potential < 0
    
    def test_lj_force_repulsive(self):
        """At r < σ, force should be repulsive (positive)."""
        force_mag, potential = lennard_jones_force(
            r=0.5, epsilon=1.0, sigma=1.0
        )
        
        assert force_mag > 0
        assert potential > 0
    
    def test_lj_force_minimum_distance(self):
        """Force should handle very small distances gracefully."""
        force_mag, potential = lennard_jones_force(
            r=1e-12, epsilon=1.0, sigma=1.0
        )
        
        assert not np.isnan(force_mag)
        assert not np.isinf(force_mag)
    
    def test_lj_force_field_class(self):
        """Test LennardJonesForce class."""
        lj = LennardJonesForce(epsilon=1.0, sigma=0.3, cutoff=1.0)
        
        positions = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [2.0, 0.0, 0.0]
        ])
        types = np.array([0, 0, 0])
        box_size = np.array([10.0, 10.0, 10.0])
        neighbor_pairs = [(0, 1), (0, 2), (1, 2)]
        
        forces, potential = lj.calculate(
            positions, types, box_size, neighbor_pairs
        )
        
        assert forces.shape == (3, 3)
        assert not np.any(np.isnan(forces))
        assert potential != 0
        
        assert np.allclose(forces[0], -forces[1], rtol=1e-5)
    
    def test_potential_curve(self):
        """Test potential curve generation."""
        lj = LennardJonesForce(epsilon=1.0, sigma=1.0, cutoff=2.5)
        
        r_values = np.linspace(0.5, 3.0, 100)
        forces, potentials = lj.get_potential_curve(r_values)
        
        assert len(forces) == len(r_values)
        assert len(potentials) == len(r_values)
        
        min_idx = np.argmin(potentials)
        r_min = r_values[min_idx]
        assert np.isclose(r_min, 2**(1/6), rtol=0.01)
    
    def test_force_shifting(self):
        """Test force shifting at cutoff."""
        lj = LennardJonesForce(epsilon=1.0, sigma=1.0, cutoff=2.5, shift=True)
        
        force_cut, pot_cut = lj.calculate_single_pair(2.5)
        assert force_cut == 0.0
        assert pot_cut == 0.0
        
        force_near, pot_near = lj.calculate_single_pair(2.49)
        assert abs(force_near) < 0.1
    
    @given(
        epsilon=st.floats(min_value=0.1, max_value=10.0),
        sigma=st.floats(min_value=0.1, max_value=2.0),
        r=st.floats(min_value=0.5, max_value=5.0)
    )
    def test_lj_force_properties(self, epsilon, sigma, r):
        """Property tests for LJ force."""
        force_mag, potential = lennard_jones_force(r, epsilon, sigma)
        
        assert not np.isnan(force_mag)
        assert not np.isnan(potential)
        
        if r < 2**(1/6) * sigma:
            assert force_mag > 0
        else:
            assert force_mag <= 0


class TestCoulombicForce:
    """Test electrostatic force calculations."""
    
    def test_coulomb_force(self):
        """Test basic Coulomb force calculation."""
        coulomb = CoulombicForce(coulomb_constant=138.935, cutoff=2.0)
        
        positions = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0]
        ])
        charges = np.array([1.0, -1.0])
        box_size = np.array([10.0, 10.0, 10.0])
        neighbor_pairs = [(0, 1)]
        
        forces, potential = coulomb.calculate(
            positions, charges, box_size, neighbor_pairs
        )
        
        assert forces.shape == (2, 3)
        
        assert forces[0, 0] < 0
        assert forces[1, 0] > 0
        assert np.allclose(forces[0], -forces[1])
        
        assert potential < 0
    
    def test_coulomb_like_charges(self):
        """Like charges should repel."""
        coulomb = CoulombicForce(cutoff=2.0)  # Use larger cutoff
        
        positions = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0]
        ])
        charges = np.array([1.0, 1.0])
        box_size = np.array([10.0, 10.0, 10.0])
        neighbor_pairs = [(0, 1)]
        
        forces, potential = coulomb.calculate(
            positions, charges, box_size, neighbor_pairs
        )
        
        assert forces[0, 0] > 0
        assert forces[1, 0] < 0
        assert potential > 0


class TestCompositeForce:
    """Test composite force field."""
    
    def test_composite_force(self):
        """Test combining multiple force fields."""
        composite = CompositeForce()
        composite.add_force_field(
            LennardJonesForce(epsilon=1.0, sigma=0.3, cutoff=1.0)
        )
        composite.add_force_field(
            CoulombicForce(coulomb_constant=138.935, cutoff=1.0)
        )
        
        positions = np.array([
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0]
        ])
        properties = {
            "types": np.array([0, 0]),
            "charges": np.array([1.0, -1.0])
        }
        box_size = np.array([10.0, 10.0, 10.0])
        neighbor_pairs = [(0, 1)]
        
        forces, potential = composite.calculate(
            positions, properties, box_size, neighbor_pairs
        )
        
        assert forces.shape == (2, 3)
        assert not np.any(np.isnan(forces))
        assert potential != 0