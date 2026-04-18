"""
Public Python API for AlkCalc.
"""

# ==============================================================================
# Imports
# ==============================================================================

from ._core import _Enlsj_core


# ==============================================================================
# Public Functions
# ==============================================================================

def Enlsj(species: str, n: int, l: int, j: float) -> float:
    """
    Eigenenergies in units of Hartree.

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
        Eigenenergy corresponding to the quantum numbers n, l, s = 1/2, and j.
    """
    return _Enlsj_core(species.encode("ascii"), n, l, j)


# ==============================================================================
# Helper functions
# ==============================================================================
