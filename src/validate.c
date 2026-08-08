/* -------------------------------------------------------------------------- *
 * Validation of input parameters                                             *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../inc/eigensolver.h"

/* Validate input parameters                                                  *
 *                                                                            *
 * In AlkCalc the user specifies parameters in interface/settings.c. However, *
 * not every possible choice of parameters is valid: there are different      *
 * constraints, physical ones such as nmax > l, but also technical            *
 * constraints, for instance on the minimal possible order of the B-splines.  *
 * This function validates the parameters set by the user and returns an      *
 * error in case a constraint is not met.                                     */
void validate_settings(int32_t nl) {

    int32_t k, N, nmax, l, J, Jlower, Jupper;
    double j, rmax;

    /* Extract parameters from settings */
    k = settings.k;
    N = settings.N;
    nmax = settings.nmax;
    l = settings.l;
    j = settings.j;
    rmax = settings.rmax;

    /* B-spline order (k)                                                     *
     *                                                                        *
     * Constraint(s) : k >= 3                                                 *
     * Information   : This makes sure that the second derivatives of         *
     *                 B-splines are still polynomials with degree larger     *
     *                 than or equal to zero. As a result, the components of  *
     *                 the stiffness matrix K in theory/theory.pdf are        *
     *                 well-defined.                                          */
    if (k < 3) {
        ERROR("B-SPLINE ORDER K = %" PRId32 " IS TOO SMALL", k);
    }

    /* Number of unique knots (N)                                             *
     *                                                                        *
     * Constraint(s) : N >= 2                                                 *
     * Information   : The number N specifies how many knots there are,       *
     *                 without counting multiplicities. The constraint        *
     *                 ensures that there is a first and a final knot. This   *
     *                 ensures that the dimension of the generalised          *
     *                 eigenvalue problem in theory/theory.pdf is larger than *
     *                 or equal to one.                                       */
    if (N < 2) {
        ERROR("TOO FEW KNOTS, N MUST SATISFY: N(%" PRId32 ") >= 2", N);
    }

    /* Maximal principal quantum number (nmax)                                *
     *                                                                        *
     * Constraint(s) : nmax >= nl, where nl is the minimal possible principal *
     *                 quantum number for fixed l                             *
     * Information   : For fixed l, a physical constraint on n is that        *
     *                 l <= n - 1. Therefore, nl = l + 1 is the minimal       *
     *                 principal quantum number that can host l. When         *
     *                 nmax < nl there does not exist an eigenstate.          */
    if (nmax < nl) {
        ERROR("INVALID NMAX: NMAX(%" PRId32 ") < NL(%" PRId32 ")", nmax, nl);
    }

    /* Orbital angular momentum quantum number (l)                            *
     *                                                                        *
     * Constraint(s) : 0 <= l < nmax                                          *
     * Information   : For fixed n, a physical constraint is l < n.           *
     *                 Therefore, l must be less than nmax.                   */
    if (l < 0 || l >= nmax) {
        ERROR("INVALID L: 0 <= L(%" PRId32 ") < NMAX(%" PRId32 ")", l, nmax);
    }

    /* Total angular momentum quantum number (j)                              *
     *                                                                        *
     * Constraint(s) : j = |l - 1 / 2|, l + 1 / 2                             *
     * Information   : By coupling l with s = 1 / 2, the only possible total  *
     *                 angular momentum quantum numbers are j = |l - s| and   *
     *                 j = l + s.                                             */
    J = (int32_t)floor(2. * j + .5);
    Jlower = (int32_t)floor(2. * fabs(l - .5) + .5);
    Jupper = (int32_t)floor(2. * (l + .5) + .5);
    if (J != Jlower && J != Jupper) {
        ERROR("INVALID J: J IS NEITHER |L - S| NOR L + S");
    }

    /* Maximally considered radius (rmax)                                     *
     *                                                                        *
     * Constraint(s) : rmax > 0                                               *
     * Information   : In theory/theory.pdf it is assumed that rmax > 0. This *
     *                 is a technical constraint which is applied to simplify *
     *                 the theory. Of course, it also does not make sense to  *
     *                 choose rmax <= 0.                                      */
    if (!(rmax > 0)) {
        ERROR("MAXIMAL RADIUS RMAX MUST BE LARGER THAN ZERO");
    }
}
