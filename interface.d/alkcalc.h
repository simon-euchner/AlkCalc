/* -------------------------------------------------------------------------- *
 * Main header file for AlkCalc                                               *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#ifndef ALCKALC_H
#define ALCKALC_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <float.h>

#define PI 3.141592653589793238462643383279502884 /* Pi */

/* -------------------------------------------------------------------------- *
 * Absolute path to 'AlkCalc' (Always put '/' at the end!)                    *
 * -------------------------------------------------------------------------- */
#define PATH_TO_ALKCALC "/home/simon/Files/GitHub/AlkCalc/"
#define LEN_PATH_TO_ALKCALC 33

/* -------------------------------------------------------------------------- *
 * Absolute path to location where radial eigenstates are stored (Always put  *
 * '/' at the end!)                                                           *
 * -------------------------------------------------------------------------- */
#define PATH_TO_STATES "/home/simon/Files/GitHub/AlkCalc/data.d/"
#define LEN_PATH_TO_STATES 40

/* -------------------------------------------------------------------------- *
 * Macros for numerical applications and error handling                       *
 * -------------------------------------------------------------------------- */
#define CONVERT(X) (int32_t)floor(2.*(X)+.5)
#define INTEGER_ABS(X) (((X) > 0) ? (X): -1*(X))
#define MAX(X, Y) ((X) > (Y)) ? (X): (Y)
#define MIN(X, Y) ((X) > (Y)) ? (Y): (X)
#define COMPLEX(X, Y) ((X)+(Y)*I)

/* -------------------------------------------------------------------------- *
 * Data type for radial eigenstates                                           *
 * -------------------------------------------------------------------------- */
typedef struct _AlkcalcState {
    int32_t N; /* Number of discretisation points */
    int32_t dim; /* Dimension of 'fnlsj', dim = N-2 */
    int32_t n; /* Prinzipal quantum number */
    int32_t l; /* Orbital angular momentum quantum number */
    double j; /* Total angular momentum quantum number */
    double *t; /* Discretisation points, t[0]=0, ..., t[N-1]=tmax */
    double *h; /* Stepsizes, h[0]=t1-t0, ..., h[N-2]=tN-1-tN-2 */
    double *fnlsj; /* Radial eigenstate in basis of elements, dimension 'dim' */
} alkcalc_state;

/* -------------------------------------------------------------------------- *
 * Data type for Clebsch-Gordan coefficients, exactly representing the        *
 * following real number:                                                     *
 *     sign * sqrt( numerator / denominator )                                 *
 * -------------------------------------------------------------------------- */
typedef struct _AlkcalcCg {
    int8_t sign;
    int64_t numerator;
    int64_t denominator;
} alkcalc_cg;

/* -------------------------------------------------------------------------- *
 * Data type for two-component spinors                                        *
 * -------------------------------------------------------------------------- */
typedef struct _AlkcalcSpinor {
    double complex u; /* First (spin-up) spinor component */
    double complex d; /* Second (spin-down= spinor component */
} alkcalc_spinor;

/* -------------------------------------------------------------------------- *
 * Eigenenergy in units of Hartree (27.211386245981(30) eV, Ref. [5])         *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * n       : Principal quantum number n = 1, 2, 3, ...                        *
 * l       : Orbital angular momentum l = 0, 1, ..., n-1                      *
 * s       : Spin (Not an argument, since we always have s = 1/2!)            *
 * j       : Total angular momentum quantum number j = |l-1/2|, l+1/2         *
 * -------------------------------------------------------------------------- */
double alkcalc_Enlsj(char *species, int32_t n, int32_t l, double j);

/* -------------------------------------------------------------------------- *
 * Radial eigenstate times radius: Rnlsj(r) = fnlsj(r/aB) / (aB^(1/3) r/aB)   *
 * Result owned by caller, destroy with 'alkcalc_state_free' after usage      *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * result  : 'f': full result; 'p': partial result (only 'fnlsj' not NULL)    *
 * species : String specifying atom/ion species                               *
 * n       : Principal quantum number n = 1, 2, 3, ...                        *
 * l       : Orbital angular momentum l = 0, 1, ..., n-1                      *
 * s       : Spin (Not an argument, since we always have s = 1/2!)            *
 * j       : Total angular momentum quantum number j = |l-1/2|, l+1/2         *
 * -------------------------------------------------------------------------- */
alkcalc_state *alkcalc_fnlsj(char result, char *species, int32_t n, int32_t l,
                             double j);

/* -------------------------------------------------------------------------- *
 * Free for data type 'alkcalc_state'                                         *
 * -------------------------------------------------------------------------- */
void alkcalc_state_free(alkcalc_state *state);

/* -------------------------------------------------------------------------- *
 * Radial matrix element <n,l,s,j|r^p|n',l',s',j'> (s=s'=1/2)                 *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * nb      : Principal quantum number of bra                                  *
 * lb      : Orbital angular momentum l = 0, 1, ..., n-1 of bra               *
 * sb      : Spin of bra (Not an argument, since we always have s = 1/2!)     *
 * jb      : Total angular momentum quantum number j = |l-1/2|, l+1/2 of bra  *
 * p       : Power of radius operator in matrix element                       *
 * nk      : Principal quantum number of ket                                  *
 * lk      : Orbital angular momentum l = 0, 1, ..., n-1, of ket              *
 * sk      : Spin of ket (Not an argument, since we always have s = 1/2!)     *
 * jk      : Total angular momentum quantum number j = |l-1/2|, l+1/2, of ket *
 * -------------------------------------------------------------------------- */
double alkcalc_rp(char *species, int32_t nb, int32_t ln, double jb, double p,
                  int32_t nk, int32_t lk, double jk);

/* -------------------------------------------------------------------------- *
 * Clebsch-Gordan coefficients (see 'theory.d/theory.pdf', section 'Manual')  *
 * -------------------------------------------------------------------------- */
alkcalc_cg alkcalc_cj1m1j2m2jmj(double j1, double m1, double j2, double m2,
                                double j, double mj);

/* -------------------------------------------------------------------------- *
 * Angular eigenstate in uncoupled basis (dimensionless)                      *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * l       : Orbital angular momentum l = 0, 1, ..., n-1                      *
 * ml      : Magnetic quantum number, ml = -l, ..., l                         *
 * s       : Spin (Not an argument, since we always have s = 1/2!)            *
 * ms      : Spin-projection quantum number ms = -1/2, 1/2                    *
 * theta   : Polar angle (Zenitwinkel), theta in [0, pi]                      *
 * phi     : Azimuthal angle (Azimut), phi in [0, 2pi]                        *
 * -------------------------------------------------------------------------- */
alkcalc_spinor alkcalc_YlmlXsms(int32_t l, int32_t ml, double ms, double theta,
                                double phi);

/* -------------------------------------------------------------------------- *
 * Angular eigenstate in coupled basis (dimensionless)                        *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * l       : Orbital angular momentum l = 0, 1, ..., n-1                      *
 * s       : Spin (Not an argument, since we always have s = 1/2!)            *
 * j       : Total angular momentum quantum number j = |l-1/2|, l+1/2         *
 * mj      : Total magnetic quantum number, mj = -j, ..., j                   *
 * theta   : Polar angle (Zenitwinkel), theta in [0, pi]                      *
 * phi     : Azimuthal angle (Azimut), phi in [0, 2pi]                        *
 * -------------------------------------------------------------------------- */
alkcalc_spinor alkcalc_Philsjmj(int32_t l, double j, double mj, double theta,
                                double phi);

/* -------------------------------------------------------------------------- *
 * Oscillator strength between fine-structure states (dimensionless)          *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * ni      : Principal quantum number of initial state (i)                    *
 * li      : Orbital angular momentum l = 0, 1, ..., n-1 of (i)               *
 * si      : Spin of (i) (Not an argument, since we always have s = 1/2!)     *
 * ji      : Total angular momentum quantum number j = |l-1/2|, l+1/2  of (i) *
 * nf      : Principal quantum number of final state (f)                      *
 * lf      : Orbital angular momentum l = 0, 1, ..., n-1 of (f)               *
 * sf      : Spin of (f) (Not an argument, since we always have s = 1/2!)     *
 * jf      : Total angular momentum quantum number j = |l-1/2|, l+1/2, of (f) *
 * -------------------------------------------------------------------------- */
double alkcalc_fitof(char *species, int32_t ni, int32_t li, double ji,
                     int32_t nf, int32_t lf, double jf);

/* -------------------------------------------------------------------------- *
 * Lifetime of fine-structure state (nanoseconds)                             *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * T       : Temperature of black-body excitation spectrum in Kelvin (K)      *
 * species : String specifying atom/ion species                               *
 * n       : Principal quantum number n = 1, 2, 3, ...                        *
 * dn      : Consider up to (including) n+dn for absorption                   *
 * l       : Orbital angular momentum l = 0, 1, ..., n-1                      *
 * s       : Spin (Not an argument, since we always have s = 1/2!)            *
 * j       : Total angular momentum quantum number j = |l-1/2|, l+1/2         *
 * -------------------------------------------------------------------------- */
double alkcalc_tau(double T, char *species, int32_t n, int32_t dn, int32_t l,
                   double j);

#endif
