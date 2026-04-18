################################################################################
### Python wrapper for the C library AlkCalc                                 ###
###                                                                          ###
### Author of this file: Simon Euchner                                       ###
################################################################################

"""
Python wrapper for the C library AlkCalc

This is a Python interfeace for the precompiled shared library 'libalkcalc.so'.
"""

# ------------------------------------------------------------------------------
# Comiler directives
# ------------------------------------------------------------------------------

# cython: language_level=3
# distutils: language = c

# ---------------------------------------------------------------------------- #
# Imports
# ---------------------------------------------------------------------------- #

from PyAlkCalc._alkcalc cimport alkcalc_Enlsj
from libc.stdint cimport int32_t

# ------------------------------------------------------------------------------
# Functions to isolate C calls from Python logic
# ------------------------------------------------------------------------------

cdef double _Enlsj(const char *species, int32_t n, int32_t l, double j):
    return alkcalc_Enlsj(species, n, l, j)

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

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

    ### Convert j to exact half-integer float
    #j = _convert(j)

    return _Enlsj(species, n, l, j)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------
