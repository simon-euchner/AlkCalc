"""
Python wrapper for the C library AlkCalc
"""


# ==============================================================================
# Imports
# ==============================================================================

# cython: language_level=3
# distutils: language = c

from libc.stdint cimport int32_t


# ==============================================================================
# C Library Wrapper
# ==============================================================================

cpdef double _Enlsj_core(bytes species, int32_t n, int32_t l, double j):
    cdef const char *species_c = species
    return alkcalc_Enlsj(species_c, n, l, j)
