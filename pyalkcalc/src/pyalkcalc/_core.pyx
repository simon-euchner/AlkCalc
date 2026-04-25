"""
Python wrapper for the C library AlkCalc
"""


# ==============================================================================
# Imports
# ==============================================================================

# cython: language_level=3
# distutils: language = c

from libc.stdint cimport int32_t
from libc.string cimport memcpy

import numpy as np
cimport numpy as np
from numpy import empty, float64, ndarray


# ==============================================================================
# C Library Wrapper (functions)
# ==============================================================================

cpdef double _enlsj_core(bytes species, int32_t n, int32_t l, double j):
    cdef const char *species_c = <const char *>species
    return alkcalc_Enlsj(species_c, n, l, j)

cpdef double _tau_core(
        double T,
        bytes species,
        int32_t n,
        int32_t dn,
        int32_t l,
        double j,
):
    cdef const char *species_c = <const char *>species
    return alkcalc_tau(T, species_c, n, dn, l, j)

cpdef double _rp_core(
        bytes species,
        int32_t nb,
        int32_t lb,
        double jb,
        double p,
        int32_t nk,
        int32_t lk,
        double jk,
):
    cdef const char *species_c = <const char *>species
    return alkcalc_rp(species_c, nb, lb, jb, p, nk, lk, jk)

cpdef double _fitof_core(
        bytes species,
        int32_t ni,
        int32_t li,
        double ji,
        int32_t nf,
        int32_t lf,
        double jf,
):
    cdef const char *species_c = <const char *>species
    return alkcalc_fitof(species_c, ni, li, ji, nf, lf, jf)


# ==============================================================================
# C Library Wrapper (classes)
# ==============================================================================

cdef class _StateCore:
    cdef alkcalc_state *state_ptr

    def __cinit__(
            self,
            bytes species,
            int32_t n,
            int32_t l,
            double j,
            bytes result,
    ):
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

    @property
    def t(self):
        if self.state_ptr.t == NULL: return empty(0, dtype=float64)
        cdef int32_t ln = self.state_ptr.N
        cdef np.ndarray[double, ndim=1, mode="c"] out = empty(ln, dtype=float64)
        memcpy(&out[0], self.state_ptr.t, ln * sizeof(double))
        return out

    @property
    def h(self):
        if self.state_ptr.h == NULL: return empty(0, dtype=float64)
        cdef int32_t ln = self.state_ptr.N - 1
        cdef np.ndarray[double, ndim=1, mode="c"] out = empty(ln, dtype=float64)
        memcpy(&out[0], self.state_ptr.h, ln * sizeof(double))
        return out

    @property
    def fnlsj(self):
        cdef int32_t ln = self.state_ptr.dim
        cdef np.ndarray[double, ndim=1, mode="c"] out = empty(ln, dtype=float64)
        memcpy(&out[0], self.state_ptr.fnlsj, ln * sizeof(double))
        return out

cdef class _CGCore:
    cdef alkcalc_cg cg

    def __cinit__(
            self,
            double j1,
            double m1,
            double j2,
            double m2,
            double j,
            double mj,
    ):
        self.cg = alkcalc_cj1m1j2m2jmj(j1, m1, j2, m2, j, mj)

    @property
    def sign(self):
        return self.cg.sign

    @property
    def numerator(self):
        return self.cg.numerator

    @property
    def denominator(self):
        return self.cg.denominator

cdef class _SpinorUncoupledBasis:
    cdef alkcalc_spinor spinor

    def __cinit__(
            self,
            int32_t l,
            int32_t ml,
            double ms,
            double theta,
            double phi,
    ):
        self.spinor = alkcalc_YlmlXsms(l, ml, ms, theta, phi)

    @property
    def u(self):
        return self.spinor.u

    @property
    def d(self):
        return self.spinor.d

cdef class _SpinorCoupledBasis:
    cdef alkcalc_spinor spinor

    def __cinit__(
            self,
            int32_t l,
            double j,
            double mj,
            double theta,
            double phi,
    ):
        self.spinor = alkcalc_Philsjmj(l, j, mj, theta, phi)

    @property
    def u(self):
        return self.spinor.u

    @property
    def d(self):
        return self.spinor.d
