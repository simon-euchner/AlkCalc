"""
Public Python API for AlkCalc.
"""

# ==============================================================================
# Imports
# ==============================================================================

from ._core import _enlsj_core, _StateCore, _tau_core
from dataclasses import dataclass


# ==============================================================================
# Dataclasses
# ==============================================================================

@dataclass(slots=True, frozen=True)
class State:
    """
    Container for state data.

    This dataclass represents the metadata associated with a computed
    radial eigenstate for a given set of quantum numbers. It provides a
    lightweight, immutable snapshot of the state returned by the underlying C
    implementation.

    Attributes
    ----------
    N : int
        Total number of discretisation points.
    dim : int
        Dimension of the basis used in the calculation.
    n : int
        Principal quantum number.
    l : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1.
    j : float
        Total angular momentum quantum number j = |l - 1/2| or |l + 1/2|.

    Notes
    -----
    The total wave-function is defined as the product of the angular spinor and
    the radial eigenstate Rnlsj. This dataclass represents Rnlsj in terms of
    fnlsj(t) = sqrt(aB) * r * Rnlsj(r), where t = r / aB with Bohr's radius aB.
    Note that fnlsj is dimensionless, because the dimension of Rnlsj is
    1 / sqrt(aB)**3.
    """
    N: int
    dim: int
    n: int
    l: int
    j: float


# ==============================================================================
# Public Functions
# ==============================================================================

def energy(species: str, n: int, l: int, j: float) -> float:
    """
    Eigenenergy of fine-structure state.

    Parameters
    ----------
    species : str
        String to specify atom/ion species, e.g., 1H for Hydrogen, or 88SR+ for
        the 88Sr+ ion.
    n : int
        Principal quantum number.
    l : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1.
    j : float
        Total angular momentum quantum number j = |l - 1/2|, |l + 1/2|.

    Returns
    -------
    float
        Eigenenergy corresponding to the quantum numbers n, l, s = 1/2, and j in
        units of Hartree.
    """
    species_ascii = species.encode("ascii")
    return _enlsj_core(species_ascii, n, l, j)

def state(species: str, n: int, l: int, j: float, result: str = "f") -> State:
    """
    Radial eigenstate.

    The total wave-function is defined as the product of the angular spinor and
    the radial eigenstate Rnlsj. This function returns the state
    fnlsj(t) = sqrt(aB) * r * Rnlsj(r) , where t = r / aB with Bohr's radius aB.
    Note that fnlsj is dimensionless, because the dimension of Rnlsj is
    1/aB**(3/2).

    Parameters
    ----------
    species : str
        String to specify atom/ion species, e.g., 1H for Hydrogen, or 88SR+ for
        the 88Sr+ ion.
    n : int
        Principal quantum number.
    l : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1.
    j : float
        Total angular momentum quantum number j = |l - 1/2|, |l + 1/2|.
    result : str
        Character for deciding if discretisation data is returned. If it is
        `"f"` (full), the full result is returned. If result is `"p"` (partial),
        the discretisation data is not returned. In the latter case, the arrays
        `t` and `h` are empty (see documentation of the dataclass: State).

    Returns
    -------
    State
        Radial eigenstate.
    """
    species_ascii = species.encode("ascii")
    result_ascii = result.encode("ascii")
    core = _StateCore(species_ascii, n, l, j, result_ascii)
    return State(
            N=core.N,
            dim=core.dim,
            n=core.n,
            l=core.l,
            j=core.j,
            )

def lifetime(T: float, species: str, n: int, dn: int, l: int,
             j: float) -> float:
    """
    Lifetime of fine-structure state.

    Parameters
    ----------
    T : float
        Temperature of black-body excitation spectrum in Kelvin.
    species : str
        String to specify atom/ion species, e.g., 1H for Hydrogen, or 88SR+ for
        the 88Sr+ ion.
    n : int
        Principal quantum number.
    dn : int
        At finite temperature, decay also takes place to states with higher
        energy. The number `dn` specifies how many states, with energy above
        the requested one, are are taken into account for calculating the
        lifetime. For instance, if `dn` is 5, decay up to (and including) states
        with n+5 is considered.
    l : int
        Orbital angular momentum quantum number l = 0, 1, ..., n-1.
    j : float
        Total angular momentum quantum number j = |l - 1/2|, |l + 1/2|.

    Returns
    -------
    float
        Lifetime of the fine-structure state with quantum numbers n, l, s = 1/2,
        and j in nanoseconds.
    """
    species_ascii = species.encode("ascii")
    return _tau_core(T, species_ascii, n, dn, l, j)
