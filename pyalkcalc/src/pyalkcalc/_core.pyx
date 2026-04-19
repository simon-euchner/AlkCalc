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
# C Library Wrapper (functions)
# ==============================================================================

cpdef double _enlsj_core(bytes species, int32_t n, int32_t l, double j):
    cdef const char *species_c = <const char *>species
    return alkcalc_Enlsj(species_c, n, l, j)


# ==============================================================================
# C Library Wrapper (classes)
# ==============================================================================

cdef class _StateCore:
    cdef alkcalc_state *state_ptr

    def __cinit__(self, bytes species, int32_t n, int32_t l, double j,
                  bytes result):
        cdef const char *species_c = <const char *>species
        cdef const char *result_c = <const char *>result
        cdef char result_char = result_c[0]
        self.state_ptr = alkcalc_fnlsj(result_char, species_c, n, l, j)

    def __dealloc__(self):
        if self.state_ptr != NULL:
            alkcalc_state_free(self.state_ptr)
            self.state_ptr = NULL

    @property
    def N(self):
        return self.state_ptr.N

    @property
    def dim(self):
        return self.state_ptr.dim

    @property
    def n(self):
        return self.state_ptr.n

    @property
    def l(self):
        return self.state_ptr.l

    @property
    def j(self):
        return self.state_ptr.j
