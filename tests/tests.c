/* -------------------------------------------------------------------------- *
 * Tests for the library functions of AlkCalc                                 *
 *                                                                            *
 * Compile: gcc -L../lib/ tests.c -lalkcalc -lm -Wl,-rpath,../lib/            *
 *          (Or, simply run 'make test' in this directory instead.)           *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * IMPORTANT: THE MASS CORRECTION MUST NOT BE INCLUDED!                       *
 *                                                                            *
 *     All tests collected in this file are performed EXCLUSIVELY for the     *
 *     species 1H, i.e., for the Hydrogen atom, because for 1H the reference  *
 *     values are known in closed form. Every reference value used below is   *
 *     the EXACT analytical resul for an INFINITELY HEAVY nucleus, that is,   *
 *     for a reduced mass which is equal to the electron's mass, me.          *
 *                                                                            *
 *     Concretely, the mass correction C, i.e., rpar[7] in src/potential.c,   *
 *     must be equal to unity, which is the DEFAULT behaviour of AlkCalc. In  *
 *     the file src/potential.c the two relevant lines must therefore read    *
 *     (the first line commented out, the second line active):
 *                                                                            *
 *         rpar[7] = 1. / (1. + ME / rpar[6]);   <- MUST STAY COMMENTED OUT   *
 *         rpar[7] = 1.;                         <- MUST BE THE ACTIVE LINE   *
 *                                                                            *
 * Data required by the tests.                                                *
 *                                                                            *
 *     Before the tests can be run, the eigenenergies and radial eigenstates  *
 *     of 1H must be generated (see README.txt, section 'Generating           *
 *     eigenenergies and radial eigenstates') for the pairs                   *
 *                                                                            *
 *         (l, j) = (0, 1/2), (1, 1/2), (1, 3/2), (2, 3/2), (2, 5/2),         *
 *                  (3, 7/2)                                                  *
 *                                                                            *
 *     keeping the parameters species, k, N, nmax, and rmax in                *
 *     interface/settings.c fixed. For 1H the minimal principal quantum       *
 *     number obeys the Hydrogenic law, nl = l + 1, so offset = -nl + 1 = -l  *
 *     be set; that is, offset = 0 for the S series, offset = -1 for the P    *
 *     series, offset = -2 for the D series, and offset = -3 for the F        *
 *     series. The maximal principal quantum number must satisfy nmax >= 10,  *
 *     and rmax must be large enough to support the state with n = nmax. The  *
 *     tolerances used below were fixed with data generated using k = 8,      *
 *     N = 2000, nmax = 20, and rmax = 10000 (default settings in             *
 *     interface/settings.c).                                                 *
 *                                                                            *
 * Tolerances.                                                                *
 *                                                                            *
 *     A test passes if the deviation from the reference value, relative for  *
 *     reference values of magnitude larger than unity and absolute           *
 *     otherwise, is smaller than the tolerance stated below. The tolerances  *
 *     are not universal: they reflect the quality of the eigensolver         *
 *     settings (k, N, rmax), the number of digits with which eigenenergies   *
 *     are stored on disk, and the fact that AlkCalc includes the spin-orbit  *
 *     coupling term VR (see src/potential.c), whereas the non-relativistic   *
 *     reference values do not. The latter is the reason why states with      *
 *     l > 0 are tested with a looser tolerance than states with l = 0.       *
 * -------------------------------------------------------------------------- */

#include <stdio.h>
#include <math.h>
#include "../interface/alkcalc.h"

#define SPECIES "1H" /* All tests are performed for the Hydrogen atom only */

#define PI 3.141592653589793238462643383279502884 /* Pi */

/* Conversion factor 2 * alpha**3 * EH / hbar in units of 1 / ns, as used by  *
 * AlkCalc in src/alkcalc.c. It is used here as well, so that the lifetime    *
 * tests probe the numerics of AlkCalc and not the values of the fundamental  *
 * constants.                                                                 */
#define TONS 32.1300103

/* Tolerances (see header)                                                    *
 *                                                                            *
 * States with l > 0 are shifted by the spin-orbit coupling term VR (see      *
 * src/potential.c), which the non-relativistic reference values do not       *
 * contain. Quantities involving such states are therefore tested with the    *
 * looser tolerances TOL_EFS, TOL_RFS, and TOL_FFS.                           */
#define TOL_EXACT 1e-13 /* Results that are exact up to rounding errors */
#define TOL_E     1e-7  /* Eigenenergies, l = 0 */
#define TOL_EFS   1e-5  /* Eigenenergies, l > 0 */
#define TOL_R     1e-7  /* Radial matrix elements, l = 0 only */
#define TOL_RFS   1e-4  /* Radial matrix elements, l > 0 involved */
#define TOL_F     1e-8  /* Radial eigenfunctions, l = 0 */
#define TOL_FFS   1e-4  /* Radial eigenfunctions, l > 0 */
#define TOL_OS    1e-4  /* Oscillator strengths */
#define TOL_TAU   1e-4  /* Lifetimes */

/* Number of tests performed and number of tests failed */
static int32_t ntests = 0, nfails = 0;

static void section(const char *);
static void check(const char *, double, double, double);
static void checkz(const char *, double complex, double complex, double);
static double cgtof(alkcalc_cg);
static double Eex(int32_t);
static double rpex(int32_t, int32_t, int32_t);
static double f1s(double);
static double f2s(double);
static double f2p(double);
static double complex Yex(int32_t, int32_t, double, double);
static void test_eigenenergies(void);
static void test_radial_matrix_elements(void);
static void test_states(void);
static void test_oscillator_strengths(void);
static void test_lifetimes(void);
static void test_clebsch_gordan(void);
static void test_uncoupled_basis(void);
static void test_coupled_basis(void);

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    printf("%s\n", "--- Tests for AlkCalc (species: " SPECIES ")");

    test_eigenenergies();
    test_radial_matrix_elements();
    test_states();
    test_oscillator_strengths();
    test_lifetimes();
    test_clebsch_gordan();
    test_uncoupled_basis();
    test_coupled_basis();

    printf("\n%s\n\n", "Summary");
    printf("TESTS PERFORMED: %" PRId32 "\n", ntests);
    printf("TESTS PASSED   : %" PRId32 "\n", ntests - nfails);
    printf("TESTS FAILED   : %" PRId32 "\n", nfails);
    printf("\n%s\n", nfails ? "--- End (FAILURE)" : "--- End (SUCCESS)");

    return nfails ? 1 : 0;
}

/* -------------------------------------------------------------------------- *
 * Tests                                                                      *
 * -------------------------------------------------------------------------- */

/* Eigenenergies in units of Hartree                                          *
 *                                                                            *
 * For 1H the exact eigenenergies are E(n) = -1 / (2 * n**2), independent of  *
 * the quantum numbers l and j. The degeneracy in l and j is lifted by the    *
 * spin-orbit coupling term VR (see src/potential.c) only, which is why       *
 * states with l > 0 are tested with the looser tolerance TOL_EFS.            */
static void test_eigenenergies(void) {

    int32_t n;
    char name[81];

    section("Eigenenergies [Hartree]");

    /* S states */
    for (n = 1; n <= 10; n++) {
        (void)sprintf(name, "E(n=%2" PRId32 ",l=0,j=1/2)", n);
        check(name, alkcalc_Enlsj(SPECIES, n, 0, .5), Eex(n), TOL_E);
    }

    /* P states */
    for (n = 2; n <= 6; n++) {
        (void)sprintf(name, "E(n=%2" PRId32 ",l=1,j=1/2)", n);
        check(name, alkcalc_Enlsj(SPECIES, n, 1, .5), Eex(n), TOL_EFS);
        (void)sprintf(name, "E(n=%2" PRId32 ",l=1,j=3/2)", n);
        check(name, alkcalc_Enlsj(SPECIES, n, 1, 1.5), Eex(n), TOL_EFS);
    }

    /* D states */
    for (n = 3; n <= 6; n++) {
        (void)sprintf(name, "E(n=%2" PRId32 ",l=2,j=3/2)", n);
        check(name, alkcalc_Enlsj(SPECIES, n, 2, 1.5), Eex(n), TOL_EFS);
        (void)sprintf(name, "E(n=%2" PRId32 ",l=2,j=5/2)", n);
        check(name, alkcalc_Enlsj(SPECIES, n, 2, 2.5), Eex(n), TOL_EFS);
    }
}

/* Radial matrix elements <n,l,s,j|r**p|n',l',s,j'> in units of Bohr's radius *
 *                                                                            *
 * The diagonal elements are compared to the exact expectation values of 1H   *
 * (see the function rpex below). For p = 0 the diagonal element is the norm  *
 * of the radial eigenstate and must equal unity.                             *
 *                                                                            *
 * An eigenvector is fixed by the eigensolver up to a global factor of -1     *
 * only. Therefore only the modulus of an off-diagonal radial matrix element  *
 * is a well-defined quantity, and only the modulus is tested here.           */
static void test_radial_matrix_elements(void) {

    section("Radial matrix elements [Bohr's radius]");

    /* Norm of the radial eigenstates */
    check("<1s|r^0|1s>", alkcalc_rp(SPECIES, 1, 0, .5, 0., 1, 0, .5),
          rpex(1, 0, 0), TOL_R);
    check("<2s|r^0|2s>", alkcalc_rp(SPECIES, 2, 0, .5, 0., 2, 0, .5),
          rpex(2, 0, 0), TOL_R);
    check("<2p|r^0|2p>", alkcalc_rp(SPECIES, 2, 1, 1.5, 0., 2, 1, 1.5),
          rpex(2, 1, 0), TOL_RFS);
    check("<3d|r^0|3d>", alkcalc_rp(SPECIES, 3, 2, 2.5, 0., 3, 2, 2.5),
          rpex(3, 2, 0), TOL_RFS);

    /* Expectation value of the radius */
    check("<1s|r^1|1s>", alkcalc_rp(SPECIES, 1, 0, .5, 1., 1, 0, .5),
          rpex(1, 0, 1), TOL_R);
    check("<3s|r^1|3s>", alkcalc_rp(SPECIES, 3, 0, .5, 1., 3, 0, .5),
          rpex(3, 0, 1), TOL_R);
    check("<2p|r^1|2p>", alkcalc_rp(SPECIES, 2, 1, 1.5, 1., 2, 1, 1.5),
          rpex(2, 1, 1), TOL_RFS);
    check("<3d|r^1|3d>", alkcalc_rp(SPECIES, 3, 2, 2.5, 1., 3, 2, 2.5),
          rpex(3, 2, 1), TOL_RFS);

    /* Expectation value of the squared radius */
    check("<1s|r^2|1s>", alkcalc_rp(SPECIES, 1, 0, .5, 2., 1, 0, .5),
          rpex(1, 0, 2), TOL_R);
    check("<2s|r^2|2s>", alkcalc_rp(SPECIES, 2, 0, .5, 2., 2, 0, .5),
          rpex(2, 0, 2), TOL_R);
    check("<2p|r^2|2p>", alkcalc_rp(SPECIES, 2, 1, 1.5, 2., 2, 1, 1.5),
          rpex(2, 1, 2), TOL_RFS);

    /* Negative powers of the radius */
    check("<1s|r^-1|1s>", alkcalc_rp(SPECIES, 1, 0, .5, -1., 1, 0, .5),
          rpex(1, 0, -1), TOL_R);
    check("<2p|r^-1|2p>", alkcalc_rp(SPECIES, 2, 1, 1.5, -1., 2, 1, 1.5),
          rpex(2, 1, -1), TOL_RFS);
    check("<2p|r^-2|2p>", alkcalc_rp(SPECIES, 2, 1, 1.5, -2., 2, 1, 1.5),
          rpex(2, 1, -2), TOL_RFS);
    check("<3d|r^-2|3d>", alkcalc_rp(SPECIES, 3, 2, 2.5, -2., 3, 2, 2.5),
          rpex(3, 2, -2), TOL_RFS);
    check("<2p|r^-3|2p>", alkcalc_rp(SPECIES, 2, 1, 1.5, -3., 2, 1, 1.5),
          rpex(2, 1, -3), TOL_RFS);
    check("<3d|r^-3|3d>", alkcalc_rp(SPECIES, 3, 2, 2.5, -3., 3, 2, 2.5),
          rpex(3, 2, -3), TOL_RFS);

    /* Off-diagonal elements (moduli only, see comment above) */
    check("|<1s|r^1|2p>|",
          fabs(alkcalc_rp(SPECIES, 1, 0, .5, 1., 2, 1, 1.5)),
          128. * sqrt(6.) / 243., TOL_RFS);
    check("|<2s|r^1|2p>|",
          fabs(alkcalc_rp(SPECIES, 2, 0, .5, 1., 2, 1, 1.5)),
          3. * sqrt(3.), TOL_RFS);
    check("|<2p|r^1|3d>|",
          fabs(alkcalc_rp(SPECIES, 2, 1, 1.5, 1., 3, 2, 2.5)),
          165888. * sqrt(5.) / 78125., TOL_RFS);
}

/* Radial eigenfunctions fnlsj(t) = t * Rnl(t), t in units of Bohr's radius   *
 *                                                                            *
 * As stated above, an eigenvector is determined up to a global factor of -1  *
 * only. Therefore the modulus of fnlsj is tested.                            */
static void test_states(void) {

    int32_t i, ltevals;
    char name[81];
    double tevals[6], t[] = {.1, .5, 1., 2., 5., 10.};

    ltevals = 6;

    section("Radial eigenfunctions |fnlsj(t)| [1 / sqrt(Bohr's radius)]");

    /* Ground state, 1S(j = 1 / 2) */
    for (i = 0; i < ltevals; i++) { tevals[i] = t[i]; }
    alkcalc_fnlsj_eval(SPECIES, 1, 0, .5, tevals, ltevals);
    for (i = 0; i < ltevals; i++) {
        (void)sprintf(name, "|f(n=1,l=0,j=1/2;t=%4.1f)|", t[i]);
        check(name, fabs(tevals[i]), fabs(f1s(t[i])), TOL_F);
    }

    /* First excited S state, 2S(j = 1 / 2) */
    for (i = 0; i < ltevals; i++) { tevals[i] = t[i]; }
    alkcalc_fnlsj_eval(SPECIES, 2, 0, .5, tevals, ltevals);
    for (i = 0; i < ltevals; i++) {
        (void)sprintf(name, "|f(n=2,l=0,j=1/2;t=%4.1f)|", t[i]);
        check(name, fabs(tevals[i]), fabs(f2s(t[i])), TOL_F);
    }

    /* Lowest P state, 2P(j = 3 / 2) */
    for (i = 0; i < ltevals; i++) { tevals[i] = t[i]; }
    alkcalc_fnlsj_eval(SPECIES, 2, 1, 1.5, tevals, ltevals);
    for (i = 0; i < ltevals; i++) {
        (void)sprintf(name, "|f(n=2,l=1,j=3/2;t=%4.1f)|", t[i]);
        check(name, fabs(tevals[i]), fabs(f2p(t[i])), TOL_FFS);
    }
}

/* Oscillator strengths (dimensionless)                                       *
 *                                                                            *
 * Summed over the total angular momentum quantum number of the final state,  *
 * the oscillator strength reduces to the result of LS coupling, which is     *
 * known in closed form for 1H. The individual fine-structure components are  *
 * obtained from the LS-coupled result by multiplication with the             *
 * corresponding angular factors (see theory/theory.pdf).                     */
static void test_oscillator_strengths(void) {

    double dE, r, fsum, fh, fl;

    section("Oscillator strengths [dimensionless]");

    /* 1S(j = 1 / 2) -> 2P(j' = 1 / 2, 3 / 2); the sum of the two components  *
     * is the well-known value f(1s -> 2p) = 2**13 / 3**9 = 0.4162            */
    fl = alkcalc_fitof(SPECIES, 1, 0, .5, 2, 1, .5);
    fh = alkcalc_fitof(SPECIES, 1, 0, .5, 2, 1, 1.5);
    fsum = 8192. / 19683.;
    check("f(1s,1/2 -> 2p,1/2)", fl, fsum / 3., TOL_OS);
    check("f(1s,1/2 -> 2p,3/2)", fh, 2. * fsum / 3., TOL_OS);
    check("f(1s,1/2 -> 2p), sum over j'", fl + fh, fsum, TOL_OS);

    /* 2P(j = 3 / 2) -> 1S(j' = 1 / 2), emission, hence negative */
    check("f(2p,3/2 -> 1s,1/2)", alkcalc_fitof(SPECIES, 2, 1, 1.5, 1, 0, .5),
          -8192. / 59049., TOL_OS);

    /* 2P(j = 1 / 2, 3 / 2) -> 3D(j' = 3 / 2, 5 / 2) */
    dE = Eex(3) - Eex(2); r = 165888. * sqrt(5.) / 78125.;
    fsum = 2. / 3. * dE * r * r * 2. / 3.;
    fl = alkcalc_fitof(SPECIES, 2, 1, .5, 3, 2, 1.5);
    fh = alkcalc_fitof(SPECIES, 2, 1, 1.5, 3, 2, 1.5)
       + alkcalc_fitof(SPECIES, 2, 1, 1.5, 3, 2, 2.5);
    check("f(2p,1/2 -> 3d), sum over j'", fl, fsum, TOL_OS);
    check("f(2p,3/2 -> 3d), sum over j'", fh, fsum, TOL_OS);

    /* Selection rules, all of the following transitions are forbidden */
    check("f(1s,1/2 -> 2s,1/2)", alkcalc_fitof(SPECIES, 1, 0, .5, 2, 0, .5),
          0., TOL_EXACT);
    check("f(1s,1/2 -> 3d,3/2)", alkcalc_fitof(SPECIES, 1, 0, .5, 3, 2, 1.5),
          0., TOL_EXACT);
    check("f(2p,1/2 -> 3d,5/2)", alkcalc_fitof(SPECIES, 2, 1, .5, 3, 2, 2.5),
          0., TOL_EXACT);
}

/* Lifetimes in units of nanoseconds                                          *
 *                                                                            *
 * The temperature of the black-body excitation spectrum is set to T = 0 K,   *
 * so that only spontaneous emission contributes. The reference values are    *
 * assembled from the exact eigenenergies and radial dipole matrix elements   *
 * of 1H, together with the angular factors of the decay channels involved    *
 * (see theory/theory.pdf). For orientation, the reference lifetimes are      *
 * 1.595 ns for 2P, 5.268 ns for 3P, and 15.459 ns for 3D.                    */
static void test_lifetimes(void) {

    double dE1, dE2, r1, r2, Gamma;

    section("Lifetimes at T = 0 K [ns]");

    /* 2P decays into 1S only, the angular factor of the channel is 1 / 3 */
    dE1 = Eex(2) - Eex(1); r1 = 128. * sqrt(6.) / 243.;
    Gamma = dE1 * dE1 * (2. / 3. * dE1 * r1 * r1 / 3.);
    check("tau(n=2,l=1,j=1/2)", alkcalc_tau(0., SPECIES, 2, 0, 1, .5),
          1. / (TONS * Gamma), TOL_TAU);
    check("tau(n=2,l=1,j=3/2)", alkcalc_tau(0., SPECIES, 2, 0, 1, 1.5),
          1. / (TONS * Gamma), TOL_TAU);

    /* 3P decays into 1S and 2S, both channels have angular factor 1 / 3 */
    dE1 = Eex(3) - Eex(1); r1 = 27. * sqrt(6.) / 128.;
    dE2 = Eex(3) - Eex(2); r2 = 27648. * sqrt(3.) / 15625.;
    Gamma = dE1 * dE1 * (2. / 3. * dE1 * r1 * r1 / 3.)
          + dE2 * dE2 * (2. / 3. * dE2 * r2 * r2 / 3.);
    check("tau(n=3,l=1,j=1/2)", alkcalc_tau(0., SPECIES, 3, 0, 1, .5),
          1. / (TONS * Gamma), TOL_TAU);
    check("tau(n=3,l=1,j=3/2)", alkcalc_tau(0., SPECIES, 3, 0, 1, 1.5),
          1. / (TONS * Gamma), TOL_TAU);

    /* 3D(j = 5 / 2) decays into 2P(j' = 3 / 2) only, since the channel       *
     * 3D(j = 5 / 2) -> 2P(j' = 1 / 2) is forbidden. The angular factor of    *
     * the allowed channel is 2 / 5.                                          */
    dE1 = Eex(3) - Eex(2); r1 = 165888. * sqrt(5.) / 78125.;
    Gamma = dE1 * dE1 * (2. / 3. * dE1 * r1 * r1 * 2. / 5.);
    check("tau(n=3,l=2,j=5/2)", alkcalc_tau(0., SPECIES, 3, 0, 2, 2.5),
          1. / (TONS * Gamma), TOL_TAU);
}

/* Clebsch-Gordan coefficients (dimensionless, exact)                         *
 *                                                                            *
 * These coefficients are computed in exact integer arithmetic and are        *
 * independent of the atom/ion species.                                       */
static void test_clebsch_gordan(void) {

    double cu, cd;

    section("Clebsch-Gordan coefficients [dimensionless]");

    /* Some coefficients with known values */
    check("c(2,-1;1,0|1,-1)",
          cgtof(alkcalc_cj1m1j2m2jmj(2., -1., 1., 0., 1., -1.)),
          -sqrt(3. / 10.), TOL_EXACT);
    check("c(1/2,-1/2;1/2,1/2|0,0)",
          cgtof(alkcalc_cj1m1j2m2jmj(.5, -.5, .5, .5, 0., 0.)),
          -sqrt(.5), TOL_EXACT);
    check("c(1/2,1/2;1/2,-1/2|0,0)",
          cgtof(alkcalc_cj1m1j2m2jmj(.5, .5, .5, -.5, 0., 0.)),
          sqrt(.5), TOL_EXACT);
    check("c(3/2,-1/2;1,1|3/2,1/2)",
          cgtof(alkcalc_cj1m1j2m2jmj(1.5, -.5, 1., 1., 1.5, .5)),
          -sqrt(8. / 15.), TOL_EXACT);
    check("c(1,1;1,-1|0,0)",
          cgtof(alkcalc_cj1m1j2m2jmj(1., 1., 1., -1., 0., 0.)),
          sqrt(1. / 3.), TOL_EXACT);
    check("c(1,0;1,0|2,0)",
          cgtof(alkcalc_cj1m1j2m2jmj(1., 0., 1., 0., 2., 0.)),
          sqrt(2. / 3.), TOL_EXACT);
    check("c(1,0;1,0|1,0)",
          cgtof(alkcalc_cj1m1j2m2jmj(1., 0., 1., 0., 1., 0.)),
          0., TOL_EXACT);

    /* Coefficients of the coupled basis states with l = 1 */
    cu = cgtof(alkcalc_cj1m1j2m2jmj(1., 0., .5, .5, 1.5, .5));
    cd = cgtof(alkcalc_cj1m1j2m2jmj(1., 1., .5, -.5, 1.5, .5));
    check("c(1,0;1/2,1/2|3/2,1/2)", cu, sqrt(2. / 3.), TOL_EXACT);
    check("c(1,1;1/2,-1/2|3/2,1/2)", cd, sqrt(1. / 3.), TOL_EXACT);
    check("Norm of |l=1,j=3/2,mj=1/2>", cu * cu + cd * cd, 1., TOL_EXACT);
    cu = cgtof(alkcalc_cj1m1j2m2jmj(1., 0., .5, .5, .5, .5));
    cd = cgtof(alkcalc_cj1m1j2m2jmj(1., 1., .5, -.5, .5, .5));
    check("c(1,0;1/2,1/2|1/2,1/2)", cu, -sqrt(1. / 3.), TOL_EXACT);
    check("c(1,1;1/2,-1/2|1/2,1/2)", cd, sqrt(2. / 3.), TOL_EXACT);
    check("Norm of |l=1,j=1/2,mj=1/2>", cu * cu + cd * cd, 1., TOL_EXACT);

    /* Coefficients that vanish, since mj is not equal to m1 + m2 */
    check("c(1,0;1/2,1/2|3/2,3/2)",
          cgtof(alkcalc_cj1m1j2m2jmj(1., 0., .5, .5, 1.5, 1.5)),
          0., TOL_EXACT);
    check("c(1,1;1/2,1/2|1/2,1/2)",
          cgtof(alkcalc_cj1m1j2m2jmj(1., 1., .5, .5, .5, .5)),
          0., TOL_EXACT);
}

/* Angular eigenstates in the uncoupled basis (dimensionless)                 *
 *                                                                            *
 * The spherical harmonics are tested against their closed-form expressions   *
 * in the Condon-Shortley convention, and the addition theorem                *
 * sum_ml |Y(l,ml)|**2 = (2 * l + 1) / (4 * pi) is verified. These tests are  *
 * independent of the atom/ion species.                                       */
static void test_uncoupled_basis(void) {

    int32_t l, ml;
    char name[81];
    double theta, phi, sum;
    alkcalc_spinor yu, yd;

    theta = .736; phi = 3.57;

    section("Uncoupled basis states [dimensionless]");

    for (l = 0; l <= 2; l++) {
        sum = 0.;
        for (ml = -l; ml <= l; ml++) {
            yu = alkcalc_YlmlXsms(l, ml, .5, theta, phi);
            yd = alkcalc_YlmlXsms(l, ml, -.5, theta, phi);

            /* Upper component for ms = 1 / 2, lower one for ms = -1 / 2 */
            (void)sprintf(name, "Y(%" PRId32 ",%+" PRId32 ")X(1/2,+1/2)_u",
                          l, ml);
            checkz(name, yu.u, Yex(l, ml, theta, phi), TOL_EXACT);
            (void)sprintf(name, "Y(%" PRId32 ",%+" PRId32 ")X(1/2,-1/2)_d",
                          l, ml);
            checkz(name, yd.d, Yex(l, ml, theta, phi), TOL_EXACT);

            /* Vanishing components */
            (void)sprintf(name, "Y(%" PRId32 ",%+" PRId32 ")X(1/2,+1/2)_d",
                          l, ml);
            checkz(name, yu.d, 0., TOL_EXACT);
            (void)sprintf(name, "Y(%" PRId32 ",%+" PRId32 ")X(1/2,-1/2)_u",
                          l, ml);
            checkz(name, yd.u, 0., TOL_EXACT);

            sum += creal(yu.u * conj(yu.u));
        }

        /* Addition theorem */
        (void)sprintf(name, "sum_ml |Y(l=%" PRId32 ",ml)|^2", l);
        check(name, sum, (2. * l + 1.) / (4. * PI), TOL_EXACT);
    }
}

/* Angular eigenstates in the coupled basis (dimensionless)                   *
 *                                                                            *
 * The coupled basis states are tested against the closed-form expressions    *
 *                                                                            *
 *     |l,1/2,l+1/2,mj> =  sqrt((l+mj+1/2) / (2l+1)) |ml=mj-1/2,ms=+1/2>      *
 *                       + sqrt((l-mj+1/2) / (2l+1)) |ml=mj+1/2,ms=-1/2>      *
 *     |l,1/2,l-1/2,mj> = -sqrt((l-mj+1/2) / (2l+1)) |ml=mj-1/2,ms=+1/2>      *
 *                       + sqrt((l+mj+1/2) / (2l+1)) |ml=mj+1/2,ms=-1/2>      *
 *                                                                            *
 * which are orthonormal for fixed mj, as it must be.                         *
 *                                                                            *
 * and the sum rule sum_mj |Phi(l,s,j,mj)|**2 = (2 * j + 1) / (4 * pi) is     *
 * verified. These tests are independent of the atom/ion species.             */
static void test_coupled_basis(void) {

    int32_t l, MJ;
    char name[81];
    double theta, phi, j, mj, sum;
    alkcalc_spinor s;

    theta = .736; phi = 3.57;

    section("Coupled basis states [dimensionless]");

    /* Phi(l = 1, s = 1 / 2, j = 3 / 2, mj = -1 / 2) */
    s = alkcalc_Philsjmj(1, 1.5, -.5, theta, phi);
    checkz("Phi(l=1,j=3/2,mj=-1/2)_u", s.u,
           sqrt(1. / 3.) * Yex(1, -1, theta, phi), TOL_EXACT);
    checkz("Phi(l=1,j=3/2,mj=-1/2)_d", s.d,
           sqrt(2. / 3.) * Yex(1, 0, theta, phi), TOL_EXACT);

    /* Phi(l = 1, s = 1 / 2, j = 1 / 2, mj = 1 / 2) */
    s = alkcalc_Philsjmj(1, .5, .5, theta, phi);
    checkz("Phi(l=1,j=1/2,mj=+1/2)_u", s.u,
           -sqrt(1. / 3.) * Yex(1, 0, theta, phi), TOL_EXACT);
    checkz("Phi(l=1,j=1/2,mj=+1/2)_d", s.d,
           sqrt(2. / 3.) * Yex(1, 1, theta, phi), TOL_EXACT);

    /* Phi(l = 2, s = 1 / 2, j = 5 / 2, mj = 5 / 2) (stretched state) */
    s = alkcalc_Philsjmj(2, 2.5, 2.5, theta, phi);
    checkz("Phi(l=2,j=5/2,mj=+5/2)_u", s.u, Yex(2, 2, theta, phi), TOL_EXACT);
    checkz("Phi(l=2,j=5/2,mj=+5/2)_d", s.d, 0., TOL_EXACT);

    /* Sum rule */
    for (l = 1; l <= 2; l++) {
        for (j = l - .5; j < l + 1.; j += 1.) {
            sum = 0.;
            for (MJ = -(int32_t)(2. * j); MJ <= (int32_t)(2. * j); MJ += 2) {
                mj = .5 * MJ;
                s = alkcalc_Philsjmj(l, j, mj, theta, phi);
                sum += creal(s.u * conj(s.u)) + creal(s.d * conj(s.d));
            }
            (void)sprintf(name,
                          "sum_mj |Phi(l=%" PRId32 ",j=%3.1f,mj)|^2", l, j);
            check(name, sum, (2. * j + 1.) / (4. * PI), TOL_EXACT);
        }
    }
}

/* -------------------------------------------------------------------------- *
 * Helper functions                                                           *
 * -------------------------------------------------------------------------- */

/* Print the headline of a section of tests                                   */
static void section(const char *title) {
    printf("\n%s\n\n", title);
    printf("%-4s %-32s %-14s  %-14s  %s\n", "", "QUANTITY", "IS", "SHOULD BE",
           "ERROR");
}

/* Compare a real result to its reference value and report the outcome        *
 *                                                                            *
 * The error is the deviation relative to the reference value, unless the     *
 * magnitude of the reference value is smaller than unity, in which case the  *
 * absolute deviation is used.                                                */
static void check(const char *name, double is, double should, double tol) {

    double scale, err;

    ntests++;
    scale = (fabs(should) > 1.) ? fabs(should) : 1.;
    err = fabs(is - should) / scale;
    if (err > tol) { nfails++; }
    printf("%-4s %-32s %+1.7E  %+1.7E  %1.1E\n", (err > tol) ? "FAIL" : "PASS",
           name, is, should, err);
}

/* Compare a complex result to its reference value and report the outcome     */
static void checkz(const char *name, double complex is, double complex should,
                   double tol) {

    double scale, err;

    ntests++;
    scale = (cabs(should) > 1.) ? cabs(should) : 1.;
    err = cabs(is - should) / scale;
    if (err > tol) { nfails++; }
    printf("%-4s %-32s (%+1.3E,%+1.3E)   %1.1E\n",
           (err > tol) ? "FAIL" : "PASS", name, creal(is), cimag(is), err);
}

/* Convert a symbolic Clebsch-Gordan coefficient into a floating-point number */
static double cgtof(alkcalc_cg c) {
    return c.sign * sqrt(c.numerator / (double)c.denominator);
}

/* Exact eigenenergy of 1H in units of Hartree (infinitely heavy nucleus)     */
static double Eex(int32_t n) {
    return -.5 / (double)(n * n);
}

/* Exact expectation value <n,l|r**p|n,l> of 1H in units of Bohr's radius     */
static double rpex(int32_t n, int32_t l, int32_t p) {

    double dn, dl, result;

    dn = (double)n; dl = (double)l;
    switch (p) {
        case 0:
            result = 1.;
            break;
        case 1:
            result = .5 * (3. * dn * dn - dl * (dl + 1.));
            break;
        case 2:
            result = .5 * dn * dn * (5. * dn * dn + 1. - 3. * dl * (dl + 1.));
            break;
        case -1:
            result = 1. / (dn * dn);
            break;
        case -2:
            result = 1. / (dn * dn * dn * (dl + .5));
            break;
        case -3:
            result = 1. / (dn * dn * dn * dl * (dl + .5) * (dl + 1.));
            break;
        default:
            result = 0.;
            break;
    }

    return result;
}

/* Exact radial eigenfunction t * R(n = 1, l = 0) of 1H                       */
static double f1s(double t) {
    return 2. * t * exp(-t);
}

/* Exact radial eigenfunction t * R(n = 2, l = 0) of 1H                       */
static double f2s(double t) {
    return t * (2. - t) * exp(-.5 * t) / (2. * sqrt(2.));
}

/* Exact radial eigenfunction t * R(n = 2, l = 1) of 1H                       */
static double f2p(double t) {
    return t * t * exp(-.5 * t) / (2. * sqrt(6.));
}

/* Exact spherical harmonics for l = 0, 1, 2 (Condon-Shortley convention)     */
static double complex Yex(int32_t l, int32_t ml, double theta, double phi) {

    double st, ct, pf;
    double complex ac, result;

    st = sin(theta); ct = cos(theta);
    ac = cos(ml * phi) + sin(ml * phi) * I;
    pf = (ml & 1) ? -1. : 1.; /* Condon-Shortley phase, ml > 0 only */
    if (ml < 0) { pf = 1.; }

    switch (10 * l + (ml < 0 ? -ml : ml)) {
        case 0: /* l = 0, ml = 0 */
            result = .5 / sqrt(PI);
            break;
        case 10: /* l = 1, ml = 0 */
            result = sqrt(.75 / PI) * ct;
            break;
        case 11: /* l = 1, |ml| = 1 */
            result = pf * sqrt(.375 / PI) * st;
            break;
        case 20: /* l = 2, ml = 0 */
            result = sqrt(5. / (16. * PI)) * (3. * ct * ct - 1.);
            break;
        case 21: /* l = 2, |ml| = 1 */
            result = pf * sqrt(15. / (8. * PI)) * st * ct;
            break;
        case 22: /* l = 2, |ml| = 2 */
            result = sqrt(15. / (32. * PI)) * st * st;
            break;
        default:
            result = 0.;
            break;
    }

    return result * ac;
}
