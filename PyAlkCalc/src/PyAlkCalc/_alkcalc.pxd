################################################################################
### Declaration of the C library AlkCalc to Cython                           ###
###                                                                          ###
### Author of this file: Simon Euchner                                       ###
################################################################################

# ------------------------------------------------------------------------------
# Imports
# ------------------------------------------------------------------------------

from libc.stdint cimport int64_t, int32_t, int8_t
from libc.complex cimport complex

# ------------------------------------------------------------------------------
# C library objects
# ------------------------------------------------------------------------------
cdef extern from "../../../interface/alkcalc.h":

    # --------------------------------------------------------------------------
    # Types
    # --------------------------------------------------------------------------

    cdef struct alkcalc_state:
        int32_t N
        int32_t dim
        int32_t n
        int32_t l
        double j
        double *t
        double *h
        double *fnlsj

    cdef struct alkcalc_cg:
        int8_t s
        int64_t numerator
        int64_t denominator

    cdef struct alkcalc_spinor:
        double complex u
        double complex d

    # --------------------------------------------------------------------------
    # Functions
    # --------------------------------------------------------------------------

    double alkcalc_Enlsj(
            const char *, int32_t, int32_t, double
            );

    alkcalc_state *alkcalc_fnlsj(
            char, const char *, int32_t, int32_t, double
            );

    void alkcalc_state_free(
            alkcalc_state *
            );

    double alkcalc_rp(
            const char *, int32_t, int32_t, double, double, int32_t, int32_t,
            double
            );

    alkcalc_cg alkcalc_cj1m1j2m2jmj(
            double, double, double, double, double, double
            );

    alkcalc_spinor alkcalc_YlmlXsms(
            int32_t, int32_t, double, double, double
            );

    alkcalc_spinor alkcalc_Philsjmj(
            int32_t, double, double, double, double
            );

    double alkcalc_fitof(
            const char *, int32_t, int32_t, double, int32_t, int32_t, double
            );

    double alkcalc_tau(
            double, const char *, int32_t, int32_t, int32_t, double
            );
