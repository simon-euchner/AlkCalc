/* -------------------------------------------------------------------------- *
 * Main program for computing eigenenergies and radial eigenstates            *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * For more information please see theory/theory.pdf.                         *
 * -------------------------------------------------------------------------- */

#include "../inc/eigensolver.h"

static double step(int32_t);
//static void save_energies(eigensolver_data *, const double *);
//static void save_states(eigensolver_data *, const double *);
//static void save_discretisation();
//static void fmt_2d_exp(char *, int32_t, double);

/* --- MAIN ----------------------------------------------------------------- */
int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    /* Initialise generalised eigenvalue problem */
    eigensolver_data *data = eigensolver_data_init();

    /* Solve eigenproblem and save result */
    solve(data);

    /* Clean up */
    //eigensolver_data_free(data);

    return 0;
}
/* -------------------------------------------------------------------------- */

/* Initialise generalised eigenvalue problem (result owned by caller)         */
eigensolver_data *eigensolver_data_init() {

    int32_t k, N, Nks, Nbs, dim, nderivKM, nderivW, nderivm, nKM, nW, *ipar, i,
            im1, imin, j, ileft, a, ia, b, ib, iarr, l;
    double rmax, *ts, *trs, *hs, *vnikx, *work, *K, *M, *W, *H, *wKM, *xKM, *wW,
           *xW, *rpar, w, t, C, vt;
    clock_t tstart, tend;
    eigensolver_data *data;

    /* Start measuremet of execution time */
    tstart = clock();

    /* Allocate memeory for result */
    data = (eigensolver_data *)malloc(sizeof(eigensolver_data));

    /* Constants                                                              *
     *                                                                        *
     * k        : Order of B-splines                                          *
     * N        : Number of knots without counting multiplicities             *
     * rmax     : Maximal radius in units of Bohr's radius                    *
     * nW       : See code below.                                             *
     * Nks      : Number of knots including mutiplicities                     *
     * Nbs      : Number of B-splines                                         *
     * dim      : Dimension of generalised eigenvalue problem                 *
     * nderivKM : Control parameter for derivatives (see DBSPVD) for K, and M *
     * nderivW  : Control parameter for derivatives (see DBSPVD) for W        *
     * nderivm  : Maximum of nderivm = max(nderivKM, nderivW)                 *
     * nKM      : Quadrature order for computing components of K and M        */
    k = settings.k;
    N = settings.N;
    rmax = settings.rmax;
    data->Nks = Nks = N - 2 + 2 * k;
    data->Nbs = Nbs = N + k - 2;
    data->dim = dim = Nbs - 2;
    nderivKM = 1 + 1;
    nderivW = 1;
    nderivm = (nderivKM > nderivW) ? nderivKM: nderivW;
    nKM = k;

    /* Order nW of Gauss-Legendre quadrature for the potential matrix W       *
     *                                                                        *
     * The effective potential in theory/theory.pdf which determines the      *
     * components of the potential matrix W is not a polynomial. Therefore,   *
     * to compute the integral it is best to choose a high quadrature order   *
     * nW. This is the point in the code where this order is hard-coded. It   *
     * be adjusted by the user, if necessary.                                 */
    nW = 1000;

    /* Allocate memory                                                        *
     *                                                                        *
     * ts    : Knots including multiplicities (k: knot vector)                *
     * trs   : Knots excluding mutliplicities (r: reduced knot vector)        *
     * hs    : Steps h[i] = trs[i] - trs[i - 1], excluding mutliplicities     *
     * vnikx : Array to store values of B-splines and their derivatives       *
     * work  : Working space for DBSPVD                                       *
     * K     : Array to hold relevant components of the stiffness matrix K    *
     * M     : Array to hold relevant components of the stiffness matrix M    *
     * W     : Array to hold relevant components of the potenial matrix W     *
     * H     : Array to hold relevant components of the matrix H = K + W      *
     * wKM   : Quadrature weights for computing components of K and M         *
     * xKM   : Quadrature points for computing components of K and M          *
     * wW    : Quadrature weights for computing components of W               *
     * xW    : Quadrature points for computing components of W                */
    ts = (double *)calloc(Nks, sizeof(double));
    trs = (double *)malloc(N * sizeof(double));
    hs = (double *)malloc((N  - 1) * sizeof(double));
    vnikx = (double *)malloc(k * nderivm * sizeof(double));
    work = (double *)malloc(((k + 1) * (k + 2)) / 2 * sizeof(double));
    K = (double *)calloc(k * dim, sizeof(double));
    data->M = M = (double *)calloc(k * dim, sizeof(double));
    W = (double *)calloc(k * dim, sizeof(double));
    data->H = H = (double *)calloc(k * dim, sizeof(double));
    wKM = (double *)malloc(nKM * sizeof(double));
    xKM = (double *)malloc(nKM * sizeof(double));
    wW = (double *)malloc(nW * sizeof(double));
    xW = (double *)malloc(nW * sizeof(double));

    /* Initialise parametric model potential V (see src/potential.c) */
    potential_initpar(ipar = data->ipar, rpar = data->rpar);

    /* Validate settings (see src/validate.c) */
    validate_settings(ipar[3]);

    /* Procedure                                                              *
     *                                                                        *
     * The goal is to construct the matrices K, M, W, and H defined in        *
     * theory/theory.pdf. All of these matrices are real-symmetric banded     *
     * matrices. The number of superdiagonals is given by the degree, d, of   *
     * the B-splines. Because of these properties, the information carried by *
     * each of the matrices K, M, W, and H is fully encoded in the d          *
     * superdiagonals and the one diagonal. Therefore, each of the matrices   *
     * can be represented by arrays of dimension (d + 1) * dim = k * dim,     *
     * where k = d + 1 is the order of the B-splines and dim is the dimension *
     * of the generalised eigenvalue problem                                  *
     *                                                                        *
     *     H fbar = lambda M f                                                *
     *                                                                        *
     * derived in theory/theory.pdf.                                          *
     *                                                                        *
     * To construct the matrices Gauss-Legendre quadratures are employed (see *
     * theory/theory.pdf). For the mass matrix M, a quadrature rule of order  *
     * n = d + 1 yields the exact result (to machine precision) for the       *
     * components. In constrast, the stiffness matrix only requires a         *
     * quadrature rule of order d. However, here, both cases are treated with *
     * order d + 1. This is numerically more efficient, because this way,     *
     * both matrices K and M can be constructed with one loop over quadrature *
     * points. This reduced the number of calls to de Boor's algorithm.       *
     *                                                                        *
     * Finally, the matrix W, which represents the effective potential, is    *
     * constructed. This is done separately, because the computation of its   *
     * components requires high quadrature orders n, as the potential is only *
     * approximately a polynomial.                                            *
     *                                                                        *
     * In the final step, the matrix H = K + W is constructed. With this, all *
     * matrices to formulate the generalised eigenvalue problem are           *
     * constructed. The next step is then to numerically solve the            *
     * generalised eigenvalue problem.                                        */

    /* Knots and step sizes */
    trs[0] = 0.; /* First knot is zero */
    for (i = 1; i < N; i++) {
        im1 = i - 1;
        hs[im1] = step(i);
        trs[i] = trs[im1] + hs[im1];
        ts[k - 1 + i] = trs[i];
    }
    for (i = N - 2 + k; i < Nks; ts[i++] = rmax); /* Add multiplicity */

    /* Compute weights and points for Gauss-Legendre quadrature rule          *
     *                                                                        *
     * The weights, wKM, and the points, xKM, are computed. Both are arrays   *
     * of length nKM. The weights and points are computed for integrals over  *
     * the interval [-1, 1]. For more information, see theory/theory.pdf.     */
    gaussq_c(&nKM, xKM, wKM); /* Call to GAUSSQ (Quadrature order n = k) */

    /* Construct stiffness matrix K and mass matrix M */
    for (i = 1; i < N; i++) { /* Loop over intervals [ts[i - 1], ts[i]] */

        /* The interval [trs[i - 1], trs[i]] correpsonds to the interval      *
         * [ts[i + k - 2], ts[i + k - 1]] in terms of the full knot           *
         * vector ts. On this interval only the B-splines with indices        *
         * imin = i - 1, ..., i + k - 2 are non-zero.                         */
        imin = i - 1;

        for (j = 0; j < nKM; j++) { /* Loop over quadrature points */

            /* Compute to interval adjusted weight, w, and point, t */
            w = .5 * hs[i - 1] * wKM[j];
            t = .5 * (hs[i - 1] * xKM[j] + trs[i] + trs[i - 1]);

            /* Largest integer satisfying ts[ileft] <= t */
            ileft = k - 1 + i - 1; /* ts[ileft] = trs[i - 1] */

            /* Evaluate B-splines and derivatives at quadrature point t */
            ileft += 1; /* Add one, because in FORTRAN counting starts at ONE */
            dbspvd_c(ts, &k, &nderivKM, &t, &ileft, vnikx, work);

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
                    M[iarr] += w * vnikx[a] * vnikx[b];
                    K[iarr] += .5 * w * vnikx[k + a] * vnikx[k + b];
                }
            }
        }
    }

    /* Compute weights and points for Gauss-Legendre quadrature rule          */
    gaussq_c(&nW, xW, wW); /* Call to GAUSSQ (Quadrature order n = nW) */

    /* Construct potential matrix W */
    C = rpar[7]; l = settings.l;
    for (i = 1; i < N; i++) { /* Loop over intervals [ts[i - 1], ts[i]] */

        imin = i - 1;

        for (j = 0; j < nW; j++) { /* Loop over quadrature points */

            /* Compute to interval adjusted weight, w, and point, t */
            w = .5 * hs[i - 1] * wW[j];
            t = .5 * (hs[i - 1] * xW[j] + trs[i] + trs[i - 1]);

            /* Largest integer satisfying ts[ileft] <= t */
            ileft = k - 1 + i - 1; /* ts[ileft] = trs[i - 1] */

            /* Evaluate B-splines and derivatives at quadrature point t */
            ileft += 1; /* Add one, because in FORTRAN counting starts at ONE */
            dbspvd_c(ts, &k, &nderivW, &t, &ileft, vnikx, work);

            /* Evaluate the effective potential (see theory/theory.pdf) at t */
            vt = C * V(t, ipar, rpar) + (l * (l + 1)) / (2. * t * t);

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
                    W[iarr] += w * vnikx[a] * vt * vnikx[b];
                }
            }
        }
    }

    /* Construct matrix H = K + M */
    for (i = 0; i < k * dim; i++) { H[i] = K[i] + W[i]; }

    /* Clean up */
    free(trs); trs = NULL;
    free(vnikx); vnikx = NULL;
    free(work); work = NULL;
    free(xKM); xKM = NULL;
    free(K); K = NULL;
    free(W); W = NULL;
    free(wKM); wKM = NULL;
    free(xKM); xKM = NULL;
    free(wW); wW = NULL;
    free(xW); xW = NULL;

    /* Compute and save partial execution time */
    tend = clock(); data->runtime = (tend - tstart)/(double)CLOCKS_PER_SEC;

    return data;
}

/* Free data of type eigensolver_data                                         */
void eigensolver_data_free(eigensolver_data *data) {
    free(data->ts); data->ts = NULL;
    free(data->hs); data->hs = NULL;
    free(data->M); data->M = NULL;
    free(data->H); data->H = NULL;
    free(data); data = NULL;
}

/* Solve generalised eigenvalue problem (result owned by caller)              */
void solve(eigensolver_data *data) {

    clock_t tstart, tend;

    /* Start measuremet of execution time */
    tstart = clock();

    /* Compute and save total execution time */
    tend = clock(); data->runtime = (tend - tstart) / (double)CLOCKS_PER_SEC;

    /* Print information */
    printf("ALGORITHM FINISHED SUCCESSFULLY (RUNTIME: %.3f S)\n\n"
           "SAVING DATA\n\n", data->runtime);

    /* Prepare and save eigenenergies */
    //iC = 1. / data->rpar[7];
    //for (k = 0; k < nev; k++) { d[k] = iC * (d[k] - offset); }
    //save_energies(data, d);

    /* Save radial eigenstates */
    //save_states(data, z);

    /* Clean up */

    /* Save discretisation points and step sizes to file */
    //save_discretisation();

    /* Print information */
    printf("DATA SAVED SUCCESSFULLY\n\n"
           "EXITING\n\n");
}

/* -------------------------------------------------------------------------- *
 * Helper functions                                                           *
 * -------------------------------------------------------------------------- */

/* Function returning step sizes (i = 1, ..., N - 1)                          */
static double step(int32_t i) {

    /* Information                                                            *
     *                                                                        *
     * Here, the step sizes hi = ti - tim1, i = 1, ..., N - 1 (im1 means      *
     * i - 1) are defined, where ti is the i-th knot WITHOUT considering      *
     * multiplicities. Note that the index i for labelling was shifted here   *
     * compared to theory/theory.pdf, e.g., ti here, in theory/theory.pdf     *
     * would be tj with j = d + i = k + i - 1, where k is the order of the    *
     * B-splines and d their degree.                                          *
     *                                                                        *
     * Visualisation.                                                         *
     *                                                                        *
     * In theory/theory.pdf:                                                  *
     *     t[0], ..., t[d], t[k], ..., t[k + N - 3], tmax, ..., tmax          *
     *     --------------                            ---------------          *
     *        k knots                                   k knots               *
     *                                                                        *
     * Here:                                                                  *
     *            0 = t[0], t[1], ..., t[N - 2]    , t[N - 1] = tmax          *
     *                                                                        *
     * The step sizes hi must be chosen such that their sum is equal to the   *
     * maximal radius rmax (see interface/settings.c). The reason the step    *
     * sizes are set, rather than the actual knots ti, is that                *
     * algorithmically the step sizes are needed. However, computing them     *
     * from the knots ti can be numerically unstable due to catastrophic      *
     * cancellation in computing the differences ti - tim1. Therefore, it is  *
     * numerically more stable to define the step sizes instead of the knots. *
     *                                                                        *
     * Note: All multiplication and addition should be done in double         *
     * precision rather than 32-bit signed integer. This is to avoid integer  *
     * overflow for large N.                                                  */

    int32_t N;
    double rmax, hi;

    /* Extract N and rmax from settings */
    N = settings.N; rmax = settings.rmax;

    /* Definition of step sizes */
    hi = rmax * (2. * i - 1.) / ((double)(N - 1) * (N - 1));

    return hi;
}

//      /* Save computed eigenenergies to file                                        */
//      static void save_energies(eigensolver_data *data, const double *energies) {
//      
//          char file[71], filename[51], buffer[101];
//          int32_t *ipar, nl, lo, jj, runtime, n;
//          double EGS, dti, dtf;
//          FILE *fd;
//      
//          /* Open file for writing */
//          nl = (ipar = data->ipar)[3]; lo = ipar[2]; jj = 2 * (int32_t)j + 1;
//          EGS = data->rpar[9]; runtime = (int32_t)data->runtime;
//          dti = step(1); dtf = step(N - 1);
//          (void)sprintf(filename, "energies-%s-%03" PRId32 "-%03" PRId32 ".dat",
//                        species, lo, jj);
//          (void)strcpy(file, "./data/");
//          (void)strcat(file, filename);
//          if (!(fd = fopen(file, "w"))) {
//              ERROR("COULD NOT OPEN FILE '%s' FOR WRITING", filename);
//          }
//      
//          /* Save metadata */
//          (void)fprintf(fd,
//                        "EIGENENERGIES FOR %s [HARTREE]\n\n"
//                        "CPU TIME TO GENERATE DATA SET [S]: %" PRId32 "\n"
//                        "GROUND STATE ENERGY [HARTREE]: %1.8lf\n"
//                        "ORBITAL ANGULAR MOMENTUM [HBAR]: %" PRId32 "\n"
//                        "TOTAL ANGULAR MOMENTUM [HBAR]: %" PRId32 "/2\n"
//                        "RMAX [BOHR'S RADIUS]: %1.3E\n"
//                        "NUMBER OF DISCRETISATION POINTS: %" PRId32 "\n"
//                        "FIRST, FINAL STEP SIZE: %1.3E, %1.3E\n"
//                        "MINIMAL PRINCIPAL QUANTUM NUMBER: %" PRId32 "\n"
//                        "MAXIMAL PRINCIPAL QUANTUM NUMBER (N): %" PRId32 "\n\n\n\n"
//                        "N   ENERGY\n\n",
//                        species, runtime, EGS, lo, jj, rmax, N, dti, dtf, nl, nmax);
//      
//          /* Save eigenenergies */
//          n = 0;
//          while (++n < nl) { (void)fprintf(fd, "%03" PRId32 "\n", n); }
//          while (n++ < nmax + 1) {
//              fmt_2d_exp(buffer, 9, energies[n - (nl - 1) - 2]);
//              (void)fprintf(fd, "%03" PRId32 " %s\n", n - 1, buffer);
//          }
//      
//          /* Close file */
//          fclose(fd); fd = NULL;
//      }
//      
//      /* Save computed radial eigenstates to file                                   */
//      static void save_states(eigensolver_data *data, const double *z) {
//      
//          char file[LEN_PATH_TO_STATES + 101], filename[101], buffer[101];
//          int32_t *ipar, nl, lo, jj, dim, n, k;
//          FILE *fd;
//      
//          /* Open file for writing */
//          nl = (ipar = data->ipar)[3]; lo = ipar[2]; jj = 2 * (int32_t)j + 1;
//          dim = data->dim;
//          for (n = nl; n < nmax + 1; n++) {
//      
//              /* Open file for writing */
//              file[0] = filename[0] = '\0';
//              (void)sprintf(filename,
//                            "state-%s-%03" PRId32 "-%03" PRId32 "-%03" PRId32 ".dat",
//                            species, n, lo, jj);
//              (void)strcpy(file, PATH_TO_STATES);
//              (void)strcat(file, filename);
//      
//              /* Save metadata */
//              if (!(fd = fopen(file, "w"))) {
//                  ERROR("COULD NOT OPEN FILE '%s' FOR WRITING", filename);
//              }
//              (void)fprintf(fd,
//                            "RADIAL EIGENSTATE FOR '%s' [DIMENSIONLESS]\n\n"
//                            "COEFFICIENTS FK (K = 1, ..., N-2)\n"
//                            "PRINCIPAL QUANTUM NUMBER (N): %" PRId32 "\n"
//                            "ORBITAL ANGULAR MOMENTUM [HBAR]: %" PRId32 "\n"
//                            "TOTAL ANGULAR MOMENTUM [HBAR]: %" PRId32 "/2\n"
//                            "RMAX [BOHR'S RADIUS]: %1.3E\n"
//                            "NUMBER OF DISCRETISATION POINTS: %" PRId32 "\n\n\n\n"
//                            "FK\n\n",
//                            species, n, lo, jj, rmax, N);
//      
//              /* Save radial eigenstate */
//              for (k = 0; k < dim; k++) {
//                  fmt_2d_exp(buffer, 15, z[dim * (n - nl) + k]);
//                  (void)fprintf(fd, "%s\n", buffer);
//              }
//      
//              /* Close file */
//              fclose(fd); fd = NULL;
//          }
//      }
//      
//      /* Save discretisation points and step sizes to file                          */
//      static void save_discretisation() {
//      
//          char file[71], filename[51], buffer_tk[101], buffer_hk[101];
//          int32_t k;
//          double tk, hk;
//          FILE *fd;
//      
//          /* Check if file already exists */
//          (void)sprintf(filename, "discretisation-%s.dat", species);
//          (void)strcpy(file, "./data/");
//          (void)strcat(file, filename);
//          if ((fd = fopen(file, "r"))) { fclose(fd); fd = NULL; return; }
//      
//          /* Save metadata */
//          if (!(fd = fopen(file, "w"))) {
//              ERROR("COULD NOT WRITE DICRETISATION DATA");
//          }
//          (void)fprintf(fd,
//                        "DISCRETISATION DATA FOR SPECIES %s\n\n"
//                        "NUMBER OF DISCRETISATION POINTS: %" PRId32 "\n\n\n\n"
//                        "K        T                     H\n\n",
//                        species, N);
//      
//          /* Save discretisation data */
//          fmt_2d_exp(buffer_tk, 15, 0.);
//          fprintf(fd, "%08" PRId32 " %s\n", 0, buffer_tk); tk = 0.;
//          for (k = 1; k < N; k++) {
//              tk += (hk = step(k));
//              fmt_2d_exp(buffer_tk, 15, tk); fmt_2d_exp(buffer_hk, 15, hk);
//              fprintf(fd, "%08" PRId32 " %s %s\n", k, buffer_tk, buffer_hk);
//          }
//      
//          /* Close file */
//          fclose(fd); fd = NULL;
//      }
//      
//      /* Formatter to ensure two-digit exponent (Number of digits: 1.234 -> nd = 4) */
//      static void fmt_2d_exp(char *buffer, int32_t nd, double x) {
//      
//          char *d;
//          int32_t len;
//          double y;
//      
//          /* IMPORTANT: The C99 standard specifies (Sec. 7.19.6.1 and               *
//           * Sec. 7.19.6.6 in Ref. [11]):                                           *
//           *                                                                        *
//           *     The sprintf function is equivalent to fprintf, ...                 *
//           *                                                                        *
//           *     ... The exponent always contains at least two digits, and only as  *
//           *     many more digits as necessary to represent the exponent.           *
//           *                                                                        *
//           * Because of this, the checks below allow one to assume that the         *
//           * exponent is printed with exactly two digits on systems that strictly   *
//           * follow the C99 standard. However, Windows does not always do this,     *
//           * which is the reason for the shift logic below. It trims a three-digit  *
//           * exponent, typically employed by Windows, to a two-digit one. Strictly  *
//           * speaking, this is not necessary; it is a nicety offered to Windows     *
//           * users.                                                                 */
//      
//          /* Check if a two-digit exponent is able to capture the number */
//          y = (x < 0) ? -x: x;
//          if (y > 1e98) {
//              ERROR("IMPOSSIBLE NUMBER DETECTED: EXPONENT OUT OF BOUNDS");
//          }
//          if ( y < 1e-98) { x = 0.; }
//      
//          /* Get the total length of the string representing x */
//          len = (int32_t)sprintf(buffer, "%+1.*E", nd - 1, x);
//      
//          /* Trim leading zero in a three-digit exponent */
//          if (len > nd + 6) {
//              d = buffer + nd + 4;
//              d[0] = d[1]; d[1] = d[2]; d[2] = '\0';
//          }
//      }
