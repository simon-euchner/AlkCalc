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

cimport alkcalc
from libc.stdint cimport int32_t
import numpy as np
cimport numpy as cnp

# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

cpdef double alkcalc_Enlsj(object species, int32_t n, int32_t l, double j):
    """
    Eigenenergies in units of Hartree.

    Parameters
    ----------
    species : str or bytes
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
    cdef bytes species_buf
    return alkcalc.alkcalc_Enlsj(_to_c_string(species, &species_buf), n, l, j)

# ------------------------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------------------------

cdef inline const char * _to_c_string(object s, bytes *buffer) except *:
    """
    Convert Python string to charactar pointer.

    Parameters
    ----------
    s : str
        Input Python string.
    buffer : bytes *
        Python object holding the UTF-8 encoded string.

    Returns
    -------
    const char *
        Pointer to UTF-8 encoded C string. This pointer is valid as long as
        the buffer is valid.
    """
    buffer[0] = s.encode("utf-8")
    return buffer[0]
