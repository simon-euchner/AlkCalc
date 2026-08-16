/* -------------------------------------------------------------------------- *
 * Main source file for AlkCalc                                               *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <cblas.h>
#include "../interface/alkcalc.h"
#include "../GAUSSQ/inc/gaussq.h"
#include "../BSPLINES/inc/bsplines.h"
#include "../MVMBLAS/inc/mvmblas.h"

#define PI 3.141592653589793238462643383279502884 /* Pi */

#define CONVERT(X) (int32_t)floor(2. * (X) + .5)
#define INTEGER_ABS(X) (((X) > 0) ? (X): -(X))
#define MAX(X, Y) (((X) > (Y)) ? (X): (Y)) /* Careful with X++ and alike! */
#define MIN(X, Y) (((X) > (Y)) ? (Y): (X))
#define COMPLEX(X, Y) ((X) + (Y) * I)
#define ERROR(...) do { \
    fprintf(stderr, "ERROR (%s:%d): ", __FILE__, __LINE__); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    exit(EXIT_FAILURE); \
} while (0)

static void move(FILE *, int32_t);
static inline double parse(const char *, int32_t);
static alkcalc_cg w3jm(int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
static int64_t s64imul(int64_t, int64_t);
static int64_t ns64imul(int32_t, const int64_t *);
static int64_t s64iadd(int64_t, int64_t);
static int64_t fac(int64_t);
static int64_t euclid(int64_t, int64_t);
static double complex Ylml(int32_t, int32_t, double, double);
static double cgtofloat(alkcalc_cg);
static double thermal_photon_occupation(double, double);
static void nextrm(const char *, int32_t *, int32_t *, int32_t, double);

/* -------------------------------------------------------------------------- *
 * Eigenenergy in units of Hartree (27.211386245981(30) eV Ref. [5])          *
 * (see theory/theory.pdf, section Manual)                                    *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * n       : Principal quantum number n = 1, 2, 3, ...                        *
 * l       : Orbital angular momentum l = 0, 1, ..., n - 1                    *
 * s       : Spin (Not an argument, since s = 1 / 2!)                         *
 * j       : Total angular momentum quantum number j = |l - 1 / 2|, l + 1 / 2 *
 * -------------------------------------------------------------------------- */
double alkcalc_Enlsj(const char *species, int32_t n, int32_t l, double j) {

    char file[LEN_PATH_TO_ALKCALC + 101], filename[101];
    int32_t J, nl, nmax, dummy;
    double E;
    FILE *fd;

    /* Open file for reading */
    J = CONVERT(j);
    (void)sprintf(filename,
                  "data/energies-%s-%03" PRId32 "-%03" PRId32 ".dat", species,
                  l, J);
    (void)strcpy(file, PATH_TO_ALKCALC);
    (void)strcat(file, filename);
    if (!(fd = fopen(file, "r"))) {
        ERROR("REQUESTED EIGENENERGY NOT AVAILABLE");
    }

    /* Extract energy */
    move(fd, 3);
    (void)fscanf(fd, "MINIMAL PRINCIPAL QUANTUM NUMBER: %" SCNd32 " ",
                 &nl);
    if (n < nl) {
        ERROR("REQUESTED EIGENENERGY DOES NOT EXIST");
    }
    (void)fscanf(fd, "MAXIMAL PRINCIPAL QUANTUM NUMBER: %" SCNd32 " ",
                 &nmax);
    if (nmax < n) {
        ERROR("REQUESTED EIGENENERGY NOT AVAILABLE");
    }
    move(fd, n + 13);
    (void)fscanf(fd, "%" SCNd32 "     %lf ", &dummy, &E);

    /* Clean up */
    fclose(fd); fd = NULL;

    return E;
}

/* -------------------------------------------------------------------------- *
 * Radial eigenstate times radius                                             *
 * Result owned by caller, destroy with alkcalc_state_free after usage        *
 * (see theory/theory.pdf, section Manual)                                    *
 *                                                                            *
 * result  : 'f': full result; 'p': partial result (only fnlsj is not NULL)   *
 * species : String specifying atom/ion species                               *
 * n       : Principal quantum number n = 1, 2, 3, ...                        *
 * l       : Orbital angular momentum l = 0, 1, ..., n - 1                    *
 * s       : Spin (Not an argument, since s = 1 / 2!)                         *
 * j       : Total angular momentum quantum number j = |l - 1 / 2|, l + 1 / 2 *
 * -------------------------------------------------------------------------- */
alkcalc_state *alkcalc_fnlsj(char result, const char *species, int32_t n,
                             int32_t l, double j) {

    char file[LEN_PATH_TO_STATES + 101], filename[101], *buffer, *bfr;
    int32_t J, k, Nks, N, Nbs, ndi, ndf, a, b, c, d, i;
    double *t, *h, *fnlsj, tmax;
    FILE *fd;
    alkcalc_state *state;

    /* Open file with requested state */
    J = CONVERT(j);
    (void)sprintf(filename,
                  "state-%s-%03" PRId32 "-%03" PRId32 "-%03" PRId32 ".dat",
                  species, n, l, J);
    (void)strcpy(file, PATH_TO_STATES);
    (void)strcat(file, filename);
    if (!(fd = fopen(file, "r"))) {
        ERROR("REQUESTED RADIAL EIGENSTATE NOT FOUND");
    }

    /* Read in metadata and move file pointer to data */
    move(fd, 5);
    (void)fscanf(fd, "ORDER OF B-SPLINES (K): %" SCNd32 " ", &k);
    (void)fscanf(fd, "TOTAL NUMBER OF KNOTS (NKS): %" SCNd32 " ", &Nks);
    (void)fscanf(fd, "NUMBER OF B-SPLINES (NBS): %" SCNd32 " ", &Nbs);
    move(fd, 6);

    /* Allocate memory for result */
    N = Nks - 2 * (k - 1);
    state = (alkcalc_state *)malloc(sizeof(alkcalc_state));
    switch (result) {
        case 'f':
            state->t = t = (double *)calloc(Nks, sizeof(double));
            state->h = h = (double *)malloc((N - 1) * sizeof(double));
            break;
        case 'p':
            state->t = t = NULL;
            state->h = h = NULL;
            break;
        default:
            ERROR("INVALID VALUE TO FIRST ARGUMENT OF ALKCALC_FNLSJ");
            break;
    }
    state->fnlsj = fnlsj = (double *)malloc(Nbs * sizeof(double));

    /* Number of digits                                                       *
     *                                                                        *
     * The number of digits used per integer (ndi) and floating-point value   *
     * (ndf) in the data files states-... and knotdata-.... For instance, if  *
     * the floating-point values are of the form +1.23E+45, ndf = 3. If the   *
     * integers, numbering the discretisation points, are of the form 0123,   *
     * ndi = 4. The integers a, b, c, and d repeatedly appear in the code.    *
     * They depend only on ndi and ndf.                                       */
    ndi = 8; ndf = 15;
    a = ndf + 7;
    b = Nbs * a - 1;
    c = ndi + 1;
    d = (N - 1) * (ndi + 1 + 2 * a) - 1;

    /* Read in radial eigenstate */
    (void)fread(buffer = (char *)malloc(b), 1, b, fd);
    for (i = 0; i < Nbs; i++) { state->fnlsj[i] = parse(buffer + a * i, ndf); }
    free(buffer); buffer = NULL;

    /* Close file */
    fclose(fd); fd = NULL;

    /* Add quantum numbers to result */
    state->n = n; state->l = l; state->j = J / 2.;

    /* Add information on B-spline basis to result */
    state->k = k; state->Nks = Nks; state->N = N; state->Nbs = Nbs;

    /* Read in knotdata (if requested, i.e., if result = 'f') */
    if (result == 'p') { return state; }
    file[0] = filename[0] = '\0';
    (void)sprintf(filename, "data/knotdata-%s.dat", species);
    (void)strcpy(file, PATH_TO_ALKCALC);
    (void)strcat(file, filename);
    if (!(fd = fopen(file, "r"))) {
        ERROR("REQUESTED KNOTDATA DOES NOT EXIST");
    }
    move(fd, 12 + k);
    (void)fread(bfr = buffer = (char *)malloc(d), 1, d, fd);
    for (i = 0; i < N - 2; i++) {
        state->t[k + i] = parse(bfr += c, ndf);
        state->h[i] = parse(bfr += a, ndf);
        bfr += a;
    }
    tmax = state->t[k + N - 2] = parse(bfr += c, ndf);
    state->h[N - 2] = parse(bfr += a, ndf);
    for (i = Nks - (k - 1); i < Nks; state->t[i++] = tmax);
    free(buffer); buffer = NULL; bfr = NULL;

    /* Close file */
    fclose(fd); fd = NULL;

    return state;
}

/* -------------------------------------------------------------------------- *
 * Free for data type alkcalc_state                                           *
 * -------------------------------------------------------------------------- */
void alkcalc_state_free(alkcalc_state *state) {
    free(state->t); state->t = NULL;
    free(state->h); state->h = NULL;
    free(state->fnlsj); state->fnlsj = NULL;
    free(state); state = NULL;
}

/* -------------------------------------------------------------------------- *
 * Evaluate radial eigenfunction fnlsj                                        *
 * (see theory/theory.pdf, section Manual)                                    *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * n       : Principal quantum number n = 1, 2, 3, ...                        *
 * l       : Orbital angular momentum l = 0, 1, ..., n - 1                    *
 * s       : Spin (Not an argument, since s = 1 / 2!)                         *
 * j       : Total angular momentum quantum number j = |l - 1 / 2|, l + 1 / 2 *
 * tevals  : Array containing points where to evaulate fnlsj                  *
 * ltevals : Length of array tevals                                           */
void alkcalc_fnlsj_eval(char *species, int32_t n, int32_t l, double j,
                        double *tevals, int32_t ltevals) {

    int32_t k, Nks, nderiv, i, ilo, ileft, mflag, a;
    double *ts, *fnlsj, *vnikx, *work, teval;
    alkcalc_state *state;

    /* Load requested radial eigenfunction */
    state = alkcalc_fnlsj('f', species, n, l, j);

    /* Extract data */
    k = state->k; Nks = state->Nks; ts = state->t; fnlsj = state->fnlsj;

    /* Allocate memory */
    nderiv = 1; /* No derivatives of B-splines are needed (see DBSPVD) */
    vnikx = (double *)malloc(k * nderiv * sizeof(double));
    work = (double *)malloc(((k + 1) * (k + 2)) / 2 * sizeof(double));

    /* Overwrite entries of tevals with fnlsj(tevals) */
    ilo = 1;
    for (i = 0; i < ltevals; i++) {

        /* Point at which to evaluate fnlsj */
        teval = tevals[i];

        /* Find largest integer satisfying ts[ileft] <= teval */
        dintrv_c(ts, &Nks, &teval, &ilo, &ileft, &mflag);

        /* Evaluate B-splines at teval */
        dbspvd_c(ts, &k, &nderiv, &teval, &ileft, vnikx, work);

        /* Compute fnlsj(teval) and store the result in tevals[i]             *
         *                                                                    *
         * On the interval [ts[ileft], ts[ileft + 1]] only the B-splines with *
         * indices a = ileft - d, ..., ileft are non-zero, where d = k - 1 is *
         * the polynomial degree of the B-splines.                            */
        tevals[i] = 0.;
        for (a = 0; a < k; a++) {
            tevals[i] += fnlsj[ileft - 1 - (k - 1) + a] * vnikx[a];
        }
    }

    /* Clean up */
    alkcalc_state_free(state); state = NULL;
    free(vnikx); vnikx = NULL;
    free(work); work = NULL;
}

/* -------------------------------------------------------------------------- *
 * Radial matrix element <n,l,s,j|r^p|n',l',s',j'> (s = s' = 1 / 2)           *
 * (see theory/theory.pdf, section Manual)                                    *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * nb      : Principal quantum number of bra                                  *
 * lb      : Orbital angular momentum l = 0, 1, ..., n - 1 of bra             *
 * sb      : Spin of bra (Not an argument, since s = 1 / 2!)                  *
 * jb      : Total angular momentum quantum number j = |l - 1 / 2|, l + 1 / 2 *
 *           of bra                                                           *
 * p       : Power of radius operator in matrix element                       *
 * nk      : Principal quantum number of ket                                  *
 * lk      : Orbital angular momentum l = 0, 1, ..., n - 1 of ket             *
 * sk      : Spin of ket (Not an argument, since s = 1 / 2!)                  *
 * jk      : Total angular momentum quantum number j = |l - 1 / 2|, l + 1 / 2 *
 *           of ket                                                           *
 * -------------------------------------------------------------------------- */
double alkcalc_rp(const char *species, int32_t nb, int32_t lb, double jb,
                  double p, int32_t nk, int32_t lk, double jk) {

    int32_t k, N, Nbs, dim, nRp, ip, i, imin, j, ileft, nderiv, a, ia, b, ib,
            iarr, d, ldRp;
    double *ts, *trs, *hs, *brafnlsj, *ketfnlsj, *ketfnlsj_cpy, *Rp, *wRp, *xRp,
           *vnikx, *work, w, t, alpha, beta, rp;
    alkcalc_state *bra, *ket;

    /* Load states */
    bra = alkcalc_fnlsj('f', species, nb, lb, jb);
    ket = alkcalc_fnlsj('p', species, nk, lk, jk);

    /* Extract data (trs: knot vector without mulitplicities) */
    k = bra->k; N = bra->N; Nbs = bra->Nbs; dim = Nbs - 2;
    ts = bra->t; trs = bra->t + bra->k - 1; hs = bra->h;
    brafnlsj = bra->fnlsj; ketfnlsj = ket->fnlsj;

    /* Order nRp of Gauss-Legendre quadrature for the matrix Rp               *
     *                                                                        *
     * When the power p is (almost; see code below for meaning of 'almost') a *
     * non-negative integer, the matrix component can be computed exactly up  *
     * to machine precision. In this case, nRp = k - 1 + p / 2 yields the     *
     * correct result for even p and nRp = k - 1 + (p - 1) / 2 yields the     *
     * correct result for odd p. In all other cases, it is best to            *
     * precision exact result for odd p. In all other cases the integral is   *
     * approximated. For this is it best to choose a high quadrature order    *
     * nRp. This is the point in the code where this order is hard-coded. It  *
     * be adjusted by the user, if necessary.                                 */
    if ((ip = (int32_t)round(p)) >= 0. && fabs(p - ip) < 1e-11) { /* Integer */
        if (ip & 1) { /* Odd */
            nRp = k - 1 + (ip - 1) / 2;
        } else { /* Even */
            nRp = k - 1 + ip / 2;
        }
    } else { /* General */
        nRp = 1000;
    }

    /* Allocate memory */
    nderiv = 1; /* No derivatives of B-splines are needed (see DBSPVD) */
    ketfnlsj_cpy = (double *)malloc(dim * sizeof(double));
    Rp = (double *)calloc(k * dim, sizeof(double));
    wRp = (double *)malloc(nRp * sizeof(double));
    xRp = (double *)malloc(nRp * sizeof(double));
    vnikx = (double *)malloc(k * nderiv * sizeof(double));
    work = (double *)malloc(((k + 1) * (k + 2)) / 2 * sizeof(double));

    /* Compute weights and points for Gauss-Legendre quadrature rule          *
     *                                                                        *
     * The weights, wRp, and the points, xRp, are computed. Both are arrays   *
     * of length nRp. The weights and points are computed for integrals over  *
     * the interval [-1, 1]. For more information, see theory/theory.pdf.     */
    gaussq_c(&nRp, xRp, wRp); /* Call to GAUSSQ (Quadrature order nRp) */

    /* Construct matrix Rp (see theory/theory.pdf) */
    for (i = 1; i < N; i++) { /* Loop over intervals [ts[i - 1], ts[i]] */

        /* The interval [trs[i - 1], trs[i]] correpsonds to the interval      *
         * [ts[i + k - 2], ts[i + k - 1]] in terms of the full knot           *
         * vector ts. On this interval only the B-splines with indices        *
         * imin = i - 1, ..., i + k - 2 are non-zero.                         */
        imin = i - 1;

        for (j = 0; j < nRp; j++) { /* Loop over quadrature points */

            /* Compute to interval adjusted weight, w, and point, t */
            w = .5 * hs[i - 1] * wRp[j];
            t = .5 * (hs[i - 1] * xRp[j] + trs[i] + trs[i - 1]);

            /* Largest integer satisfying ts[ileft] <= t */
            ileft = k - 1 + i - 1; /* ts[ileft] = trs[i - 1] */

            /* Evaluate B-splines at quadrature point t */
            ileft += 1; /* Add one, because in FORTRAN counting starts at ONE */
            dbspvd_c(ts, &k, &nderiv, &t, &ileft, vnikx, work);

            /* Accumulate matrix components of K and M (order: column-major) */
            for (a = 0; a < k; a++) {
                ia = imin + a; /* Index of relevant B-spline */
                if (ia == 0 || ia == Nbs - 1) { continue; }
                for (b = 0; b <= a; b++) {
                    ib = imin + b; /* Index of relevant B-spline */
                    if (ib == 0 || ib == Nbs - 1) { continue; }

                    /* Array index for column-major upper (U) (see EIGLAPACK) */
                    iarr = (ia - 1) * k + k - (a - b) - 1;

                    /* Accumulate matrix components */
                    Rp[iarr] += w * vnikx[a] * pow(t, p) * vnikx[b];
                }
            }
        }
    }

    /* Compute action of Rp on bra; dsbmv: y -> y = alpha * A * x + beta * y */
    d = k - 1; alpha = 1.; ldRp = k; beta = 0.;
    for (i = 0; i < dim; i++) { ketfnlsj_cpy[i] = ketfnlsj[1 + i]; }
    dsbmv_c(&dim, &d, &alpha, Rp, &ldRp, ketfnlsj_cpy, &beta, ketfnlsj + 1);

    /* Compute radial matrix element */
    rp = 0.;
    for (i = 1; i < dim + 1; i++) {
        rp += brafnlsj[i] * ketfnlsj[i];
    }

    /* Clean up */
    free(ketfnlsj_cpy); ketfnlsj_cpy = NULL;
    free(Rp); Rp = NULL;
    free(wRp); wRp = NULL;
    free(xRp); xRp = NULL;
    free(vnikx); vnikx = NULL;
    free(work); work = NULL;
    alkcalc_state_free(bra); bra = NULL;
    alkcalc_state_free(ket); ket = NULL;

    return rp;
}

/* -------------------------------------------------------------------------- *
 * Clebsch-Gordan coefficients (see theory/theory.pdf, section Manual)        *
 * -------------------------------------------------------------------------- */
alkcalc_cg alkcalc_cj1m1j2m2jmj(double j1, double m1, double j2, double m2,
                                double j, double mj) {

    int32_t J1, M1, J2, M2, J, MJ, gcd;
    alkcalc_cg result;

    /* Convert to half-integers and multiply by 2 */
    J1 = CONVERT(j1); M1 = CONVERT(m1);
    J2 = CONVERT(j2); M2 = CONVERT(m2);
    J = CONVERT(j); MJ = CONVERT(mj);

    /* Wigner's 3jm symbol (no '/ 2' necessary, see source for 'w3jm') */
    result = w3jm(J1, M1, J2, M2, J, -MJ);

    /* Add phase and scaling factor */
    result.sign *= (((-J1 + J2 - MJ) / 2) & 1) ? -1 : 1;
    result.numerator = s64imul(result.numerator, J + 1);

    /* Clean up result */
    gcd = euclid(result.numerator, result.denominator);
    result.numerator /= gcd; result.denominator /= gcd;

    return result;
}

/* -------------------------------------------------------------------------- *
 * Angular eigenstate in uncoupled basis (dimensionless)                      *
 * (see theory/theory.pdf, section Manual)                                    *
 *                                                                            *
 * l       : Orbital angular momentum l = 0, 1, ..., n - 1                    *
 * ml      : Magnetic quantum number, ml = -l, ..., l                         *
 * s       : Spin (Not an argument, since s = 1 / 2!)                         *
 * ms      : Spin-projection quantum number ms = -1 / 2, 1 / 2                *
 * theta   : Polar angle (Zenitwinkel), theta in [0, pi]                      *
 * phi     : Azimuthal angle (Azimut), phi in [0, 2pi]                        *
 * -------------------------------------------------------------------------- */
alkcalc_spinor alkcalc_YlmlXsms(int32_t l, int32_t ml, double ms, double theta,
                                double phi) {

    double complex y;
    alkcalc_spinor spinor;

    /* Check input validity */
    if (l < 0 || INTEGER_ABS(ml) > l) {
        ERROR("INVALID ORBITAL ANGULAR MOMENTUM");
    }
    if (theta < 0 || PI < theta) {
        ERROR("INVALID POLAR ANGLE");
    }
    if (phi < 0 || 2.*PI < phi) {
        ERROR("INVALID AZIMUTHAL ANGLE");
    }

    /* Compute value of spherical harmonic */
    y = Ylml(l, ml, theta, phi);

    /* Spinor structure (only sign of spin-projection quantum number matters) */
    if (ms < 0) {
        spinor.u = COMPLEX(0., 0.); spinor.d = y;
    } else {
        spinor.u = y; spinor.d = COMPLEX(0., 0.);
    }

    return spinor;
}

/* -------------------------------------------------------------------------- *
 * Angular eigenstate in coupled basis (dimensionless)                        *
 * (see theory/theory.pdf, section Manual)                                    *
 *                                                                            *
 * l       : Orbital angular momentum l = 0, 1, ..., n - 1                    *
 * s       : Spin (Not an argument, since s = 1 / 2!)                         *
 * j       : Total angular momentum quantum number j = |l - 1 / 2|, l + 1 / 2 *
 * mj      : Total magnetic quantum number, mj = -j, ..., j                   *
 * theta   : Polar angle (Zenitwinkel), theta in [0, pi]                      *
 * phi     : Azimuthal angle (Azimut), phi in [0, 2pi]                        *
 * -------------------------------------------------------------------------- */
alkcalc_spinor alkcalc_Philsjmj(int32_t l, double j, double mj, double theta,
                                double phi) {

    int32_t ll, J, MJ;
    double cgu, cgd;
    double complex yu, yd;
    alkcalc_spinor spinor;

    /* Check input validity */
    ll = 2 * l; J = CONVERT(j); MJ = CONVERT(mj);
    if (J < INTEGER_ABS(MJ)) { /* Ensure mj <= j */
        ERROR("J MUST BE LARGER EQUAL ABSOLUTE VALUE OF MJ");
    }
    if ((J & 1) != (MJ & 1)) { /* Check J, MJ both (half-) integer */
        ERROR("HALF-INTEGER J/MJ BUT INTEGER MJ/J");
    }
    if (J != ll - 1 && J != ll + 1) { /* Check |l - 1 / 2| <= j <= l + 1 / 2 */
        ERROR("J MUST BE |L - 1 / 2| OR L + 1 / 2");
    }
    if (theta < 0 || PI < theta) {
        ERROR("INVALID POLAR ANGLE");
    }
    if (phi < 0 || 2.*PI < phi) {
        ERROR("INVALID AZIMUTHAL ANGLE");
    }

    /* Compute Clebsch-Gordan coefficients (spin up [u] and down [d]) */
    cgu = cgtofloat(alkcalc_cj1m1j2m2jmj(l, mj - .5, .5, .5, j, mj));
    cgd = cgtofloat(alkcalc_cj1m1j2m2jmj(l, mj + .5, .5, -.5, j, mj));

    /* Compute spherical harmonics */
    yu = Ylml(l, (MJ - 1) / 2, theta, phi);
    yd = Ylml(l, (MJ + 1) / 2, theta, phi);

    /* Assemble result */
    spinor.u = cgu*yu; spinor.d = cgd*yd;

    return spinor;
}

/* -------------------------------------------------------------------------- *
 * Oscillator strength between fine-structure states (dimensionless)          *
 * (see theory/theory.pdf, section Manual)                                    *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * ni      : Principal quantum number of initial state (i)                    *
 * li      : Orbital angular momentum l = 0, 1, ..., n - 1, of (i)            *
 * si      : Spin of (i) (Not an argument, since s = 1 / 2!)                  *
 * ji      : Total angular momentum quantum number j = |l - 1 / 2|, l + 1 / 2 *
 *           of (i)                                                           *
 * nf      : Principal quantum number of final state (f)                      *
 * lf      : Orbital angular momentum l = 0, 1, ..., n - 1 of (f)             *
 * sf      : Spin of (f) (Not an argument, since s = 1 / 2!)                  *
 * jf      : Total angular momentum quantum number j = |l - 1 / 2|, l + 1 / 2 *
 *           of (f)                                                           *
 * -------------------------------------------------------------------------- */
double alkcalc_fitof(const char *species, int32_t ni, int32_t li, double ji,
                     int32_t nf, int32_t lf, double jf) {

    int32_t JI, JF, llp1, llm1;
    double Efi, r, al, fitof;

    /* Apply selection rules */
    JI = CONVERT(ji); JF = CONVERT(jf);
    if (JI < 0 || JF < 0 || li < 0 || lf < 0) { return 0; }
    if (INTEGER_ABS(JF-JI) > 2 || INTEGER_ABS(lf-li) != 1) { return 0.; }

    /* Energy difference between initial (i) and final (f) state in Hartree */
    Efi = alkcalc_Enlsj(species, nf, lf, jf)
        - alkcalc_Enlsj(species, ni, li, ji);

    /* Radial dipole-transition matrix element between (i) and (f) */
    r = alkcalc_rp(species, ni, li, ji, 1., nf, lf, jf);

    /* Angular factor */
    llp1 = 2 * li + 1; llm1 = 2 * li - 1;
    if (lf == li + 1) { /* lf - li = 1 */
        if (JI == llp1 && JF == llp1) {
            al = 1. / ((llp1 + 2.) * llp1);
        } else
        if (JI == llm1 && JF == llp1) {
            al = (li + 1.) / llp1;
        } else
        if (JI == llp1 && JF == llp1 + 2) {
            al = (li + 2.) / (llp1 + 2);
        } else {
            return 0;
        }
    } else { /* lf - li = -1 */
        if (JI == llm1 && JF == llm1) {
            al = 1. / (llp1 * llm1);
        } else
        if (JI == llp1 && JF == llm1) {
            al = (double)li / llp1;
        } else
        if (JI == llm1 && JF == llm1 - 2) {
            al = (li - 1.) / llm1;
        } else {
            return 0;
        }
    }

    /* Assemble result */
    fitof = 2. / 3. * Efi * r * r * al;

    return fitof;
}

/* -------------------------------------------------------------------------- *
 * Lifetime of fine-structure state (nanoseconds)                             *
 * (see theory/theory.pdf, section Manual)                                    *
 *                                                                            *
 * T       : Temperature of black-body excitation spectrum in Kelvin (K)      *
 * species : String specifying atom/ion species                               *
 * n       : Principal quantum number n = 1, 2, 3, ...                        *
 * dn      : Consider up to (including) n+dn for absorption                   *
 * l       : Orbital angular momentum l = 0, 1, ..., n - 1                    *
 * s       : Spin (Not an argument, since s = 1 / 2!)                         *
 * j       : Total angular momentum quantum number j = |l - 1 / 2|, l + 1 / 2 *
 * -------------------------------------------------------------------------- */
double alkcalc_tau(double T, const char *species, int32_t n, int32_t dn,
                   int32_t l, double j) {

    int32_t lp, lm, nmnlp1, nmxlp1, nlp1, nmnlm1, nmxlm1, nlm1, J, k;
    double En, Gamma, jp, jm, hnu, fftoi, nocc, tau;

    /* Get lowest n' such that E(n,l,s,l+s) < E(n',l',s,l'+s) is still true */
    lp = l + 1; lm = l - 1;
    En = alkcalc_Enlsj(species, n, l, j);
    nextrm(species, &nmnlp1, &nmxlp1, lp, lp + .5); /* l' = l + 1 */
    nlp1 = (n < nmnlp1) ? nmnlp1 : n;
    if (alkcalc_Enlsj(species, nlp1, lp, lp + .5) > En) {
        while (    nmnlp1 < --nlp1
                && alkcalc_Enlsj(species, nlp1, lp, lp + .5) > En );
        nlp1++;
    } else {
        while (alkcalc_Enlsj(species, ++nlp1, lp, lp + .5) < En);
    }
    if (!l) { nlm1 = -1; nmnlm1 = 0; goto SkipedSState; } /* l' = l - 1 */
    nextrm(species, &nmnlm1, &nmxlm1, lm, lm + .5);
    nlm1 = (n < nmnlm1) ? nmnlm1 : n;
    if (alkcalc_Enlsj(species, nlm1, lm, lm + .5) > En) {
        while (    nmnlm1 < --nlm1
                && alkcalc_Enlsj(species, nlm1, lm, lm + .5) > En );
        nlm1++;
    } else {
        while (alkcalc_Enlsj(species, ++nlm1, lm, lm + .5) < En);
    }
SkipedSState:

    /* Compute decay rate Gamma */
    Gamma = 0.;
    J = CONVERT(j); jp = l + .5; jm = l - .5;
    if (J == 2 * l + 1) { /* j = l + s */

        /* Emission: l' = l + 1 */
        for (k = nlp1 - 1; k >= nmnlp1; k--) {

            /* j' = l + s */
            hnu = En - alkcalc_Enlsj(species, k, lp, jp);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lp, jp);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu * hnu * fftoi * (1. + nocc);

            /* j' = l + 3s */
            hnu = En - alkcalc_Enlsj(species, k, lp, jp + 1.);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lp, jp + 1.);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu * hnu * fftoi * (1. + nocc);
        }

        /* Emission: l' = l - 1 */
        for (k = nlm1 - 1; k >= nmnlm1; k--) {

            /* j' = l - s */
            hnu = En - alkcalc_Enlsj(species, k, lm, jm);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lm, jm);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu * hnu * fftoi * (1. + nocc);
        }

        /* Absorption: l' = l + 1 */
        for (k = nlp1; k < nlp1 + dn; k++) {

            /* j' = l + s */
            hnu = alkcalc_Enlsj(species, k, lp, jp) - En;
            fftoi = alkcalc_fitof(species, n, l, j, k, lp, jp);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu * hnu * fftoi * nocc;

            /* j' = l + 3 * s */
            hnu = alkcalc_Enlsj(species, k, lp, jp + 1.) - En;
            fftoi = alkcalc_fitof(species, n, l, j, k, lp, jp + 1.);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu * hnu * fftoi * nocc;
        }

        /* Absorption: l' = l - 1 */
        if (l) {
            for (k = nlm1; k < nlm1 + dn; k++) {

                /* j' = l - s */
                hnu = alkcalc_Enlsj(species, k, lm, jm) - En;
                fftoi = alkcalc_fitof(species, n, l, j, k, lm, jm);
                nocc = thermal_photon_occupation(hnu, T);
                Gamma += hnu * hnu * fftoi * nocc;
            }
        }
    } else { /* j = l - s */

        /* Emission: l' = l + 1 */
        for (k = nlp1 - 1; k >= nmnlp1; k--) {

            /* j' = l + s */
            hnu = En - alkcalc_Enlsj(species, k, lp, jp);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lp, jp);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu * hnu * fftoi * (1. + nocc);
        }

        /* Emission: l = l - 1 */
        for (k = nlm1 - 1; k >= nmnlm1; k--) {

            /* j' = l - s */
            hnu = En - alkcalc_Enlsj(species, k, lm, jm);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lm, jm);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu * hnu * fftoi * (1. + nocc);

            /* j' = l - 3s */
            if (l > 1) { /* P(j = 1 / 2) -> S(j' = -1 / 2) is not possible */
                hnu = En - alkcalc_Enlsj(species, k, lm, jm - 1.);
                fftoi = -alkcalc_fitof(species, n, l, j, k, lm, jm - 1.);
                nocc = thermal_photon_occupation(hnu, T);
                Gamma += hnu * hnu * fftoi * (1. + nocc);
            }
        }

        /* Absorption: l' = l + 1 */
        for (k = nlp1; k < nlp1 + dn; k++) {

            /* j' = l + s */
            hnu = alkcalc_Enlsj(species, k, lp, jp) - En;
            fftoi = alkcalc_fitof(species, n, l, j, k, lp, jp);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu * hnu * fftoi * nocc;
        }

        /* Absorption: l = l - 1 */
        if (l) {
            for (k = nlm1; k < nlm1 + dn; k++) {

                /* j' = l - s */
                hnu = alkcalc_Enlsj(species, k, lm, jm) - En;
                fftoi = alkcalc_fitof(species, n, l, j, k, lm, jm);
                nocc = thermal_photon_occupation(hnu, T);
                Gamma += hnu * hnu * fftoi * nocc;

                /* j' = l - 3 * s */
                if (l > 1) { /* P(j = 1 / 2) -> S(j' = -1 / 2) not possible */
                    hnu = alkcalc_Enlsj(species, k, lm, jm - 1.) - En;
                    fftoi = alkcalc_fitof(species, n, l, j, k, lm, jm - 1.);
                    nocc = thermal_photon_occupation(hnu, T);
                    Gamma += hnu * hnu * fftoi * nocc;
                }
            }
        }
    }

    /* Compute lifetime in units of nanoseconds                               *
     *                                                                        *
     * The conversion factor used below is 2 * alpha**3 * EH / hbar, where    *
     * alpha is the fine-structure constant, EH is the Hartree, and hbar is   *
     * the reduced Planck constant; for their values, see Ref. [5].           */
    tau = 1. / (32.1300103 * Gamma);

    return tau;
}

/* -------------------------------------------------------------------------- *
 * Helper functions                                                           *
 * -------------------------------------------------------------------------- */

/* Move file descriptor down by nlines many lines                             */
static void move(FILE *fd, int32_t nlines) {
    int c;
    int32_t k;
    for (k = 0; k < nlines; k++) {
        while ((c = fgetc(fd)) != '\n' && c != EOF);
        if (c == EOF) { break; };
    }
}

/* Optimised parser for reading data files quickly                            */
static inline double parse(const char *str, int32_t nd) {

    int8_t s, i;
    int32_t k;
    uint64_t p10, dec;
    double r;

    /* Get sign */
    s = (str[0] == '-') ? -1 : 1;

    /* Leading integer */
    i = str[1] - '0';

    /* Get decimal places */
    p10 = 1; dec = str[3] - '0'; /* First digit */
    for (k = 4; k < nd + 2; k++) { p10 *= 10; dec = 10 * dec + str[k] - '0'; }

    /* Result without power */
    r = s * (i + (double)dec * .1 / p10);

    /* Get exponent */
    s = (str[nd + 3] == '-') ? -1 : 1;
    i = s*((str[nd + 4] - '0') * 10 + (str[nd + 5] - '0'));

    /* Add power to result */
    while (i-- > 0) r *= 10.;
    while (++i < 0) r *= .1;

    return r;
}

/* Wigner's 3jm symbols (arguments must be TWICE the desired arguments)       */
static alkcalc_cg w3jm(int32_t j1, int32_t m1, int32_t j2, int32_t m2,
                       int32_t j3, int32_t m3) {

    /* IMPORTANT: The arguments j1, m1, j2, m2, j3, and m3 must be TWICE the  *
     * desired argument, i.e., the following equality between Wigner's 3jm    *
     * symbols and the function w3jm holds:                                   *
     *                                                                        *
     *   / j1 j2 j3  \                                                        *
     *  |            | = w3jm(2 * j1, 2 * m1, 2 * j2, 2 * m2, 2 * j3, 2 * m3) *
     *  \  m1 m2 m3 /                                                         *
     *                                                                        *
     * This allows to work with integers throughout this function.            */

    int8_t phase, s;
    int64_t f[7], fs, K, N, g[6], *A, k, l, m, SN, SD;

    /* Prepare result */
    alkcalc_cg result;
    result.sign = 1; result.numerator = 0; result.denominator = 1;

    /* Check if ji and mi (i = 1, 2, 3) are compatible */
    if (    (j1 & 1 && !(m1 & 1)) || (m1 & 1 && !(j1 & 1))
         || (j2 & 1 && !(m2 & 1)) || (m2 & 1 && !(j2 & 1))
         || (j3 & 1 && !(m3 & 1)) || (m3 & 1 && !(j3 & 1)) ) { return result; }

    /* Kronecker delta */
    if (m1 + m2 + m3) { return result; }

    /* Phase */
    phase = (((j1 - j2 - m3) / 2) & 1) ? -1 : 1;

    /* Compute factorials in prefactors of sum */
    f[0] = fac((j1 + j2 - j3) / 2);
    f[1] = fac((j1 - j2 + j3) / 2);
    f[2] = fac((-j1 + j2 + j3) / 2);
    f[3] = s64imul(fac((j1 - m1) / 2), fac((j1 + m1) / 2));
    f[4] = s64imul(fac((j2 - m2) / 2), fac((j2 + m2) / 2));
    f[5] = s64imul(fac((j3 - m3) / 2), fac((j3 + m3) / 2));
    f[6] = fac((j1 + j2 + j3) / 2 + 1); /* Used later */
    if (!(fs = ns64imul(6, f))) { return result; }

    /* Bounds for summation */
    K = MAX(0, MAX((j2 - j3 - m1) / 2, (j1 - j3 + m2) / 2));
    N = MIN((j1 + j2 - j3) / 2, MIN((j1 - m1) / 2, (j2 + m2) / 2));

    /* Summation */
    if (N < K) { return result; }
    A = (int64_t *)malloc((N - K + 1) * sizeof(int64_t));
    for (k = K; k < N + 1; k++) {
        s = (k & 1) ? -1 : 1;
        g[0] = fac(k);
        g[1] = fac((j1 + j2 - j3) / 2 - k);
        g[2] = fac((j1 - m1) / 2 - k);
        g[3] = fac((j2 + m2) / 2 - k);
        g[4] = fac((j3 - j2 + m1) / 2 + k);
        g[5] = fac((j3 - j1 - m2) / 2 + k);
        A[k - K] = s64imul(s, ns64imul(6, g)); /* Denominators of summands */
    }
    SN = 0; /* Numerator of sum */
    for (k = 0; k < N - K + 1; k++) {
        m = 1;
        for (l = 0; l < N - K + 1; l++) {
            if (l == k) { continue; }
            m = s64imul(m, A[l]);
        }
        SN = s64iadd(SN, m);
    }
    SD = 1; /* Denominator of sum */
    for (k = 0; k < N - K + 1; k++) { SD = s64imul(SD, A[k]); }
    free(A);

    /* Assemble result */
    phase *= (SN < 0) ? -1 : 1; phase *= (SD < 0) ? -1 : 1;
    result.sign = phase;
    result.numerator = s64imul(fs, s64imul(SN, SN));
    result.denominator = s64imul(f[6], s64imul(SD, SD));

    return result;
}

/* Secure 64-bit integer multiplication                                       */
static int64_t s64imul(int64_t a, int64_t b) {
    if (!a || !b) { return 0; }
    if (a == INT64_MIN) {
        if (b == 1) { return a; } else { goto s64imulOverflow; }
    }
    if (b == INT64_MIN) {
        if (a == 1) { return b; } else { goto s64imulOverflow; }
    }
    if (    (a > 0 && b > 0 && a <= INT64_MAX / b)
         || (a < 0 && b < 0 && -a <= INT64_MAX/ (-b))
         || (a > 0 && b < 0 && -a >= INT64_MIN/ (-b))
         || (a < 0 && b > 0 && -b >= INT64_MIN/ (-a)) ) {
        return a*b;
    }
s64imulOverflow:
    ERROR("OVERFLOW IN INTEGER MULTIPLICATION");
}

/* Secure 64-bit integer multiplication (n times)                             */
static int64_t ns64imul(int32_t n, const int64_t *a) {
    int32_t k; int64_t result = 1;
    for (k = 0; k < n; result = s64imul(result, a[k++]));
    return result;
}

/* Secure 64-bit integer addition                                             */
static int64_t s64iadd(int64_t a, int64_t b) {

    /* It might look dangerous to do INT64_MIN - a when a can be equal to     *
     * INT64_MIN. However, in the particular order the subtraction is         *
     * performed, the C99 standard guarantees that the expression evaluates   *
     * to zero (see Sec. 6.5.5 in Ref. [11]).                                 */

    if (a >= 0 && b <= INT64_MAX - a) { return a + b; }
    if (a < 0 && a >= INT64_MIN && b >= INT64_MIN - a) { return a + b; }
    ERROR("OVERFLOW IN INTEGER ADDITION");
}

/* Integer factorial                                                          */
static int64_t fac(int64_t n) {
    if (n < 0) { return 0; }
    int64_t l, m = 1;
    for (l = 0; l < n-1; l++) { m = s64imul(m, n - l); }
    return m;
}

/* Euclidean algorithm                                                        */
static int64_t euclid(int64_t a, int64_t b) {
    int64_t c;
    while (b != 0) { c = a % b; a = b; b = c; }
    return a;
}

/* Spherical harmonics (see definition in theory/theory.pdf)                  */
static double complex Ylml(int32_t l, int32_t ml, double theta, double phi) {

    int8_t sml, phase;
    int32_t k;
    double pf, x, Pk, Pkm1, Pkm2;
    double complex ac, y;

    /* Azimuthal contribution (here, the sign of ml is still need) */
    ac = COMPLEX(cos(ml * phi), sin(ml * phi));

    /* Check input regime and react accordingly */
    sml = (ml < 0) ? -1 : 1; ml *= sml;
    if (l < ml) { return COMPLEX(0., 0.); }

    /* Prefactor (the Condon-Shortley phase is in Legendre polynomials) */
    pf = sqrt((2 * l + 1) * fac(l - ml) / (4. * PI * fac(l + ml)));

    /* Polar contribution                                                     *
     *                                                                        *
     * The associated Legendre polynomial's value is computed for the         *
     * absolute value of ml. The conversion formula, described in             *
     * theory/theory.pdf, allows to obtain the value for negative ml, and the *
     * conversion factor is not just a phase. Please note that this factor is *
     * not missing here, but included in the prefactor, pf, already, because  *
     * at the point where pf is computed, ml is already rendered              *
     * non-negative. This means automatically everything (up to a phase       *
     * included later) is correct. It is probably vital to view               *
     * theory/theory.pdf to understand this part.                             */
    phase = (ml & 1) ? -1 : 1; x = cos(theta);
    if (l == ml) {
        Pk = phase * fac(2 * ml) / fac(ml) * pow(.25 * (1. - x * x), .5 * ml);
    } else
    if (l == ml + 1) {
        Pkm1 = phase * fac(2 * ml) / fac(ml) * pow(.25 * (1. - x * x), .5 * ml);
        Pk = (2 * (ml + 1) - 1) * x * Pkm1;
    } else {
        Pkm2 = phase * fac(2 * ml) / fac(ml) * pow(.25 * (1. - x * x), .5 * ml);
        Pkm1 = (2 * ml + 1) * x * Pkm2; k = ml + 2;
        do {
            Pk = ((2 * k - 1) * x * Pkm1 - (k + ml - 1) * Pkm2) / (k - ml);
            Pkm2 = Pkm1; Pkm1 = Pk;
        } while (k++ < l);
    }

    /* Assemble result (here, the phase is included, as noted above) */
    y = (sml < 0) ? phase : 1.; /* Phase from conversion: ml -> -ml */
    y *= pf * Pk * ac;

    return y;
}

/* Convert symbolic Clebsch-Gordan coefficient into a floating-point number */
static double cgtofloat(alkcalc_cg c) {
    return c.sign * sqrt(c.numerator / (double)c.denominator);
}

/* Thermal photon-occupation number at energy hnu and temperature T */
static double thermal_photon_occupation(double hnu, double T) {

    double r, x;

    /* Zero T case */
    if (T <= 0.) { return 0.; }

    /* Ratio of photon and thermal energy                                     *
     *                                                                        *
     * - Photon energy h x nu (hnu) in units of Hartree (EH)                  *
     * - Temperature T in units of Kelvin (K)                                 */
    r = hnu / (3.166811e-6 * T); /* For Boltzmann's constant see Ref. [5] */

    /* Compute photon occupation number according to Planck's law */
    if (r < cbrt(720. * DBL_EPSILON)) { /* High T */
        return 1. / r - .5 + 1. / 12. * r;
    }
    if (r > -1. / 3. * log(DBL_EPSILON)) { /* Low T */
        x = exp(-r); return x + x * x;
    }
    return 1. / expm1(r); /* Intermediate regime */
}

/* Fetch extremal principal quantum numbers */
static void nextrm(const char *species, int32_t *nmin, int32_t *nmax, int32_t l,
                   double j) {

    char file[LEN_PATH_TO_ALKCALC + 101], filename[101];
    int32_t J;
    FILE *fd;

    /* Open file for reading */
    J = CONVERT(j);
    (void)sprintf(filename,
                  "data/energies-%s-%03" PRId32 "-%03" PRId32 ".dat", species,
                  l, J);
    (void)strcpy(file, PATH_TO_ALKCALC);
    (void)strcat(file, filename);
    if (!(fd = fopen(file, "r"))) {
        ERROR("REQUESTED ENERGY SERIES IS NOT AVAILABLE");
    }

    /* Extract information */
    move(fd, 3);
    (void)fscanf(fd, "MINIMAL PRINCIPAL QUANTUM NUMBER: %" SCNd32 " ", nmin);
    (void)fscanf(fd, "MAXIMAL PRINCIPAL QUANTUM NUMBER: %" SCNd32 " ",
                 nmax);

    /* Clean up */
    fclose(fd); fd = NULL;
}
