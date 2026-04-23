"""
Tests for the public Python API of pyalkcalc.
"""

# ==============================================================================
# Imports
# ==============================================================================

import pytest
import numpy as np
from pyalkcalc import (
        energy,
        state,
        radial_matrix_element,
        clebsch_gordan_coefficient,
        spinor_uncoupled_basis,
        spinor_coupled_basis,
        oscillator_strength,
        lifetime,
        State,
        CG,
        Spinor,
)


# ==============================================================================
# Helpers and Configuration
# ==============================================================================

### The following species is used for testing atomic data functions.
### To test a different species (e.g., "88SR+"), change this value.
TEST_SPECIES = "40CA+"

def has_test_data(species: str = TEST_SPECIES):
    """
    Check if data for the test species is present.
    """
    try:
        ### Try to get energy for a common state
        energy(species, 4, 0, .5)
        return True
    except (RuntimeError, SystemExit, Exception):
        ### AlkCalc often exits or errors if data is missing
        return False

DATA_REASON = f"Data for {TEST_SPECIES} not present in AlkCalc installation"


# ==============================================================================
# Tests: Atomic Data
# ==============================================================================

@pytest.mark.skipif(not has_test_data(), reason=DATA_REASON)
def test_energy_data():
    """
    Test the energy function.
    """
    en = energy(TEST_SPECIES, 4, 0, .5)
    assert isinstance(en, float)
    ### Physical states should have negative energy
    assert en < 0

@pytest.mark.skipif(not has_test_data(), reason=DATA_REASON)
def test_state_data():
    """
    Test the state function.
    """
    s = state(TEST_SPECIES, 4, 0, .5)
    assert isinstance(s, State)
    assert s.n == 4
    assert s.l == 0
    assert s.j == .5
    assert s.t.size > 0
    assert s.fnlsj.size > 0

@pytest.mark.skipif(not has_test_data(), reason=DATA_REASON)
def test_radial_matrix_element_data():
    """
    Test radial matrix element.
    """
    ### Overlap of state with itself
    overlap = radial_matrix_element(
            TEST_SPECIES, 4, 0, .5, .0, 4, 0, .5
    )
    assert pytest.approx(overlap, abs=1e-2) == 1.

@pytest.mark.skipif(not has_test_data(), reason=DATA_REASON)
def test_oscillator_strength_data():
    """
    Test oscillator strength (4S -> 4P).
    """
    f = oscillator_strength(TEST_SPECIES, 4, 0, .5, 4, 1, 1.5)
    assert isinstance(f, float)
    assert f > 0

@pytest.mark.skipif(not has_test_data(), reason=DATA_REASON)
def test_lifetime_data():
    """
    Test lifetime for 4P state.
    """
    tau = lifetime(.0, TEST_SPECIES, 4, 0, 1, 1.5)
    assert isinstance(tau, float)
    assert tau > 0


# ==============================================================================
# Tests: Dataclasses
# ==============================================================================

def test_state_dataclass():
    """
    Test State dataclass instantiation.
    """
    s = State(
            N=10,
            dim=8,
            n=1,
            l=0,
            j=.5,
            t=np.zeros(10),
            h=np.zeros(9),
            fnlsj=np.zeros(8),
    )
    assert s.N == 10
    assert s.dim == 8
    assert s.n == 1

def test_cg_dataclass():
    """
    Test CG dataclass instantiation.
    """
    cg = CG(sign=1, numerator=1, denominator=2)
    assert cg.sign == 1
    assert cg.numerator == 1
    assert cg.denominator == 2

def test_spinor_dataclass():
    """
    Test Spinor dataclass instantiation.
    """
    sp = Spinor(u=1.+0j, d=.0+1j)
    assert sp.u == 1.+0j
    assert sp.d == .0+1j


# ==============================================================================
# Tests: Angular Momentum
# ==============================================================================

def test_clebsch_gordan_analytic():
    """
    Test analytic Clebsch-Gordan coefficient.
    """

    ### C(1, 1, 0, 0, 1, 1) = 1
    cg = clebsch_gordan_coefficient(1., 1., .0, .0, 1., 1.)

    assert isinstance(cg, CG)
    assert cg.sign == 1
    assert cg.numerator == 1
    assert cg.denominator == 1

def test_clebsch_gordan_half_integer():
    """
    Test CG with half-integer values.
    """

    ### C(.5, .5, .5, -.5, 1., .0) = 1/sqrt(2)
    cg = clebsch_gordan_coefficient(.5, .5, .5, -.5, 1., .0)

    assert cg.sign == 1
    assert cg.numerator == 1
    assert cg.denominator == 2

def test_clebsch_gordan_numeric():
    """
    Test numeric Clebsch-Gordan coefficient.
    """

    ### C(1, 1, 0, 0, 1, 1) = 1
    val = clebsch_gordan_coefficient(
            1., 1., .0, .0, 1., 1., result="numeric"
    )

    assert isinstance(val, float)
    assert pytest.approx(val) == 1.

def test_clebsch_gordan_error():
    """
    Test Clebsch-Gordan error handling.
    """
    with pytest.raises(ValueError):
        clebsch_gordan_coefficient(
                1., 1., .0, .0, 1., 1., result="invalid"
        )

def test_spinor_uncoupled():
    """
    Test spinor in uncoupled basis.
    """

    ### ms = .5 (up)
    sp = spinor_uncoupled_basis(0, 0, .5, .0, .0)

    assert isinstance(sp, Spinor)
    assert isinstance(sp.u, complex)
    assert isinstance(sp.d, complex)

def test_spinor_coupled():
    """
    Test spinor in coupled basis.
    """
    sp = spinor_coupled_basis(0, .5, .5, .0, .0)
    assert isinstance(sp, Spinor)
    assert isinstance(sp.u, complex)
    assert isinstance(sp.d, complex)
