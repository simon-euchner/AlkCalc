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
    //solve(data);

    /* Clean up */
    //eigensolver_data_free(data);

    return 0;
}
/* -------------------------------------------------------------------------- */

/* Initialise generalised eigenvalue problem (result owned by caller)         */
eigensolver_data *eigensolver_data_init() {

    int32_t nW, i, im1, k, N, Nkns, Nbs, dim, *ipar, nKM;
    double *ts, *tkns, *hs, *K, *W, *M, *H, *wKM, *xKM, *wW, *xW, *rpar;
    clock_t tstart, tend;
    eigensolver_data *data;

    /* Set order nW of Gauss-Legendre quadrature for the potential matrix W   *
     *                                                                        *
     * The effective potential in theory/theory.pdf which determines the      *
     * components of the potential matrix W is not a polynomial. Therefore,   *
     * to compute the integral it is best to choose a high quadrature order   *
     * nW. This is the point in the code where this order is hard-coded. It   *
     * be adjusted by the user, if necessary.                                 */
    nW = 1000;

    /* Allocate memory */
    data = (eigensolver_data *)malloc(sizeof(eigensolver_data));
    k = settings.k; N = settings.N;
    data->Nkns = Nkns = N - 2 + 2 * k; /* Number of knots with multiplicities */
    data->Nbs = Nbs = N + k - 2; /* Number of B-splines */
    data->dim = dim = Nbs - 2; /* Dimension of generalised eigenvalue problem */
    ts = (double *)malloc(N * sizeof(double)); /* Knots, no multipl. */
    tkns = (double *)calloc(Nkns, sizeof(double)); /* Full knot vector */
    hs = (double *)malloc((N  - 1) * sizeof(double)); /* Steps, no multipl. */
    K = (double *)calloc(k * dim, sizeof(double)); /* Stiffness matrix */
    W = (double *)calloc(k * dim, sizeof(double)); /* Potential matrix */
    data->M = M = (double *)calloc(k * dim, sizeof(double)); /* Mass matrix */
    data->H = H = (double *)calloc(k * dim, sizeof(double)); /* H = K + W */
    wKM = (double *)malloc(k * sizeof(double)); /* Weights for K, M integrals */
    xKM = (double *)malloc(k * sizeof(double)); /* Points for K, M integrals */
    wW = (double *)malloc(nW * sizeof(double)); /* Weights for W integrals */
    xW = (double *)malloc(nW * sizeof(double)); /* Points for W integrals */

    /* Start measuremet of execution time */
    tstart = clock();

    /* Initialise parametric model potential V */
    potential_initpar(ipar = data->ipar, rpar = data->rpar);

    /* Validate settings */
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
    ts[0] = 0.; /* First knot is zero */
    for (i = 1; i < N; i++) {
        im1 = i - 1;
        hs[im1] = step(i);
        ts[i] = ts[im1] + hs[im1];
        tkns[k - 1 + i] = ts[i];
    }
    for (i = N - 2 + k; i < Nkns; tkns[i++] = settings.rmax); /* Multiplicity */

    for (i = 0; i < Nkns; i++) {
        printf("%1.3E\n", tkns[i]);
    }


    /* Compute weights and points for Gauss-Legendre quadrature rule          *
     *                                                                        *
     * The weights, wKM, and the points, xKM, are computed. Both are arrays   *
     * of length k. The weights and points are computed for integrals over    *
     * the interval [-1, 1]. For more information, see theory/theory.pdf.     */
    nKM = k;
    gaussq_c(&nKM, xKM, wKM); /* Call to GAUSSQ */

















    /* Compute weights and points for Gauss-Legendre quadrature rule          *
     *                                                                        *
     * The weights, wW, and the points, xW, are computed. Both are arrays of  *
     * length nW. The weights and points are computed for integrals over the  *
     * interval [-1, 1]. For more information, see theory/theory.pdf.         */
    gaussq_c(&nW, xW, wW); /* Call to GAUSSQ */

















    /* Clean up */
    free(hs); hs = NULL;
    free(ts); ts = NULL;
    free(wKM); wKM = NULL;
    free(xKM); xKM = NULL;
    free(wW); wW = NULL;
    free(xW); xW = NULL;
    free(K); K = NULL;
    free(W); W = NULL;

    /* Compute and save partial execution time */
    tend = clock(); data->runtime = (tend - tstart)/(double)CLOCKS_PER_SEC;

    return data;
}

/* Free data of type eigensolver_data                                         */
void eigensolver_data_free(eigensolver_data *data) {
    free(data->M); data->M = NULL;
    free(data->H); data->H = NULL;
    free(data); data = NULL;
}


//      
//      /* Solve eigenproblem                                                         */
//      void solve(eigensolver_data *data) {
//      
//          int ido, nerr, info;
//          int32_t dim, n, nl, nev, ncv, ldv, ldz, lworkl, *iparam, *ipntr, k;
//          double tol, *resid, *v, *workd, *workl, *select, *d, *z, sigma, runtime,
//                 *dummy, iC;
//          clock_t tstart, tend;
//      
//          /* Initialise variables for Lanczos algorithm */
//          dim = data->dim;
//          ido = 0; n = dim; nl = data->ipar[3]; nev = nmax - nl + 1;
//          if ((ncv = 2 * nev + 1) < 20) { ncv = 20; }
//          if (ncv > dim) { ncv = dim; }
//          ldv = dim; ldz = dim; lworkl = ncv * (8 + ncv); info = 0; tol = 1e-12;
//          iparam = (int32_t *)calloc(11, sizeof(int32_t));
//          ipntr = (int32_t *)calloc(11, sizeof(int32_t));
//          resid = (double *)malloc(dim * sizeof(double));
//          v = (double *)malloc(dim * ncv * sizeof(double));
//          workd = (double *)malloc(3 * dim * sizeof(double));
//          workl = (double *)malloc(lworkl * sizeof(double));
//          select = (double *)calloc(ncv, sizeof(double)); /* Ritz value ordering */
//          d = (double *)malloc(nev * sizeof(double));
//          z = (double *)malloc(nev * dim * sizeof(double));
//          sigma = shift;
//          iparam[0] = 1;
//          iparam[2] = 1000000000; /* Large enough to avoid becoming a problem */
//          iparam[3] = 1;
//          iparam[6] = 3; /* Shift-invert mode */
//      
//          /* Iterative calls to 'DSAUPD' */
//          tstart = clock(); nerr = 0;
//          do {
//      
//              /* Call 'DSAUPD' */
//              dsaupd_c(&ido, &n, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr,
//                       workd, workl, &lworkl, &info); data->info = info;
//      
//              /* Check if call was successful */
//              if (ido != 1 && ido != -1 && ido != 2 && ido != 99) {
//                  ERROR("ERROR DURING ITERATION: IDO = %d", ido);
//                  nerr++;
//              }
//              if (info != 0 && info != 1) {
//                  ERROR("ERROR DURING ITERATION: INFO = %d", info);
//                  nerr++;
//              }
//              if (info == 1) {
//                  ERROR("REACHED MAXIMAL NUMBER OF ITERATIONS");
//                  nerr++;
//              }
//              if (nerr) { ERROR("DSAUPD ENDED WITH NERR = %d", nerr); }
//      
//              /* React to instructions from 'DSAUPD' */
//              if (ido == 1) { /* Compute action of shift-inverted Hamiltonian */
//                  for (k = 0; k < dim; k++) {
//                      workd[ipntr[1] - 1 + k] = workd[ipntr[2] - 1 + k];
//                  }
//                  shift_invert_f(data, &workd[ipntr[1] - 1]); /* Result in argument */
//              } else
//              if (ido == 2) { /* Compute action of mass matrix */
//                  mass_matrix_f(data, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
//              } else { /* Initialisation step */
//                  mass_matrix_f(data, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
//                  shift_invert_f(data, &workd[ipntr[1] - 1]); /* Result in argument */
//              }
//      
//          } while (ido == 1 || ido == 2 || ido == -1);
//      
//          /* Call 'DSEUPD' to extract results */
//          dseupd_c(select, d, z, &ldz, &n, &nev, &tol, resid, &ncv, v, &ldv, &sigma,
//                   iparam, ipntr, workd, workl, &lworkl, &info);
//          tend = clock();
//          runtime = (tend - tstart) / (double)CLOCKS_PER_SEC;
//          data->runtime += runtime;
//          printf("ALGORITHM FINISHED SUCCESSFULLY (RUNTIME: %.3f S)\n\n", runtime);
//      
//          /* Dummy data to give 'Destroy_Dense_Matrix' something to free */
//          dummy = (double *)calloc(1, sizeof(double));
//          ((DNformat *)(data->B.Store))->nzval = dummy;
//      
//          /* Prepare and save eigenenergies */
//          iC = 1. / data->rpar[7];
//          for (k = 0; k < nev; k++) { d[k] = iC * (d[k] - offset); }
//          save_energies(data, d);
//      
//          /* Save radial eigenstates */
//          save_states(data, z);
//      
//          /* Clean up */
//          free(iparam); iparam = NULL;
//          free(ipntr); ipntr = NULL;
//          free(resid); resid = NULL;
//          free(v); v = NULL;
//          free(workd); workd = NULL;
//          free(workl); workl = NULL;
//          free(select); select = NULL;
//          free(d); d = NULL;
//          free(z); z = NULL;
//      
//          /* Save discretisation points and step sizes to file */
//          save_discretisation();
//      }

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
     * would be tj with j = d + i = k + (i - 1), where k is the order of the  *
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
     *                t[0], t[1], ..., t[N - 2]    , t[N - 1] = tmax          *
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

//      
//      /* Compute action of mass matrix (result stored in y)                         */
//      static void mass_matrix_f(const eigensolver_data *data, const double *x,
//                                double *y) {
//      
//          int32_t k, k0, dim = data->dim;
//          double *Mdata = data->Mdata;
//      
//          y[0] = Mdata[0] * x[0] + Mdata[1] * x[1];
//          for (k = 1; k < dim - 1; k++) {
//              k0 = 2 + 3 * (k - 1);
//              y[k] = Mdata[k0] * x[k - 1]
//                   + Mdata[k0 + 1] * x[k]
//                   + Mdata[k0 + 2] * x[k + 1];
//          }
//          k0 = 2 + 3 * (dim - 2);
//          y[dim - 1] = Mdata[k0] * x[dim - 2] + Mdata[k0 + 1] * x[dim - 1];
//      }
//      
//      /* Compute action of shift-inverted Hamiltonian (result stored in x)          */
//      static void shift_invert_f(eigensolver_data *data, double *x) {
//      
//          /* Prepare input */
//          ((DNformat *)(data->B.Store))->nzval = x;
//      
//          /* Solve the system (H-sigma*M) * vout = vin */
//          dgstrs(NOTRANS, &data->L, &data->U, data->perm_c, data->perm_r, &data->B,
//                 &data->stat, &data->info);
//      }
//      
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
//                        "EIGENENERGIES FOR '%s' [HARTREE]\n\n"
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
//                            "COEFFICIENTS 'FK' (K = 1, ..., N-2)\n"
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
//                        "DISCRETISATION DATA FOR SPECIES '%s'\n\n"
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
