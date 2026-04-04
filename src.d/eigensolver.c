/* -------------------------------------------------------------------------- *
 * Main program for computing eigenenergies and radial eigenstates            *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * All eigenenergies are in units of Hartree and the radial eigenstates are   *
 * represented by a vector in the basis of the finite elements; see           *
 * 'theory.d/theory.pdf'.                                                     *
 * -------------------------------------------------------------------------- */

#include "../inc.d/eigensolver.h"

/* --- MAIN ----------------------------------------------------------------- */
int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    /* Initialise eigenproblem */
    eigensolver_data *data = eigensolver_data_init();

    /* Compute eigenenergies */
    solve(data);

    /* Save discretisation points and step sizes to file */
    save_discretisation();

    /* Clean up */
    eigensolver_data_free(data);

    return 0;
}
/* -------------------------------------------------------------------------- */

/* Initialise eigenproblem (result owned by caller)                           */
eigensolver_data *eigensolver_data_init() {

    int *info;
    int32_t dim, *Hrs, *Hcs, *perm_r, *perm_c, *ipar, lo, k, k0;
    double *vs, *hs, *Mdata, *Hdata, *b, *rpar, C, tk, runtime;
    clock_t tstart, tend;
    SuperMatrix H, *B;
    superlu_options_t options;
    SuperLUStat_t *stat;
    eigensolver_data *data;

    /* Allocate memory */
    data = (eigensolver_data *)malloc(sizeof(eigensolver_data));
    data->dim = dim = N - 2;
    vs = (double *)malloc(dim * sizeof(double));
    hs = (double *)malloc((N - 1) * sizeof(double));
    data->Mdata = Mdata = (double *)malloc((3 * dim - 2) * sizeof(double));
    Hdata = (double *)malloc((3 * dim - 2) * sizeof(double));
    Hrs = (int32_t *)malloc((3 * dim - 2) * sizeof(int32_t));
    Hcs = (int32_t *)malloc((dim + 1) * sizeof(int32_t));
    b = (double *)calloc(dim, sizeof(double));
    data->perm_r = perm_r = (int32_t *)malloc(dim * sizeof(int32_t));
    data->perm_c = perm_c = (int32_t *)malloc(dim * sizeof(int32_t));

    /* Initialise potential (see '../interface.d/settings.c') */
    tstart = clock();
    vint_initpar(rpar = data->rpar, ipar = data->ipar);

    /* Sanity check for specified maximal principal quantum number */
    if (nmax < ipar[3]) {
        ERROR("TOO SMALL MAXIMAL PRINCIPAL QUANTUM NUMBER");
    }

    /* Sanity check for number of discretisation points */
    if (N < 3) {
        ERROR("TOO FEW DISCRETISATION POINTS");
    }

    /* Compute potential vector and step sizes */
    lo = ipar[2]; /* Orbital angular momentum */
    C = rpar[7]; /* Mass correction */
    tk = 0.; /* t0 = 0 */
    for (k = 1; k < N - 1; k++) {
        tk += (hs[k - 1] = step(k));
        vs[k - 1] = lo * (lo + 1.) / (2. * tk * tk)
                  + C * vint(tk, rpar, ipar) + offset-shift;
    }
    hs[N - 2] = step(N - 1);

    /* Mass matrix (scalar product for generalised eigenproblem [row-major])  */
    Mdata[0] = (hs[0] + hs[1]) / 3.;
    Mdata[1] = hs[1] / 6.;
    for (k = 1; k < dim - 1; k++) {
        k0 = 2 + 3 * (k - 1);
        Mdata[k0] = hs[k] / 6.;
        Mdata[k0 + 1] = (hs[k] + hs[k + 1]) / 3.;
        Mdata[k0 + 2] = hs[k + 1] / 6.;
    }
    k0 = 2 + 3 * (dim - 2);
    Mdata[k0] = hs[dim - 1] / 6.;
    Mdata[k0 + 1] = (hs[dim - 1] + hs[dim]) / 3.;

    /* Hamiltonian in CSC format */
    Hcs[0] = 0;
    Hdata[0] = .5 / hs[0] + .5 / hs[1] + vs[0] * Mdata[0];
    Hdata[1] = -.5 / hs[1] + .5 * (vs[0] + vs[1]) * Mdata[1];
    Hrs[0] = 0; Hrs[1] = 1;
    Hcs[1] = 2;
    for (k = 1; k < dim - 1; k++) {
        k0 = 2 + 3 * (k - 1);
        Hdata[k0] = -.5 / hs[k] + .5 * (vs[k - 1] + vs[k]) * Mdata[k0];
        Hdata[k0 + 1] = .5 / hs[k] + .5 / hs[k + 1] + vs[k] * Mdata[k0 + 1];
        Hdata[k0 + 2] = -.5 / hs[k + 1]
                      + .5* (vs[k] + vs[k + 1]) * Mdata[k0 + 2];
        Hrs[k0] = k - 1; Hrs[k0 + 1] = k; Hrs[k0 + 2] = k + 1;
        Hcs[k + 1] = Hcs[k] + 3;
    }
    k0 = 2 + 3 * (dim - 2);
    Hdata[k0] = -.5 / hs[N - 3] + .5* (vs[N - 4] + vs[N - 3]) * Mdata[k0];
    Hdata[k0 + 1] = .5 / hs[N - 3] + .5 / hs[N - 2] + vs[N - 3] * Mdata[k0 + 1];
    Hrs[k0] = dim - 2; Hrs[k0 + 1] = dim - 1;
    Hcs[dim] += Hcs[dim - 1] + 2;
    dCreate_CompCol_Matrix(&H, dim, dim, 3 * dim - 2, Hdata, Hrs, Hcs, SLU_NC,
                           SLU_D, SLU_GE);

    /* Initialise solver */
    set_default_options(&options); options.ColPerm = NATURAL;
    StatInit(stat=&data->stat);

    /* Dummy right-hand side (b is already initialised to hold zeros) */
    dCreate_Dense_Matrix(B=&data->B, dim, 1, b, dim, SLU_DN, SLU_D, SLU_GE);

    /* Perform LU decomposition */
    tend = clock();
    runtime = (tend-tstart)/(double)CLOCKS_PER_SEC; data->runtime = runtime;
    tstart = clock(); info = &data->info;
    dgssv(&options, &H, perm_c, perm_r, &data->L, &data->U, B, stat, info);
    tend = clock(); runtime = (tend-tstart)/(double)CLOCKS_PER_SEC;
    if (*info) {
        ERROR(("LU-DECOMPOSITION FAILED WITH 'INFO = %d'"), *info);
    } else {
        printf("LU-DECOMPOSITION SUCCESSFUL (RUNTIME: %1.3f S)\n\n", runtime);
    }
    data->runtime += runtime;

    /* Clean up */
    free(vs); vs = NULL;
    free(hs); hs = NULL;
    free(b); b = NULL;
    Destroy_CompCol_Matrix(&H);

    return data;
}

/* Free data of type 'eigensolver_data'                                       */
void eigensolver_data_free(eigensolver_data *data) {
    SUPERLU_FREE(data->perm_r); data->perm_r = NULL;
    SUPERLU_FREE(data->perm_c); data->perm_c = NULL;
    free(data->Mdata); data->Mdata = NULL;
    Destroy_Dense_Matrix(&data->B);
    Destroy_SuperNode_Matrix(&data->L);
    Destroy_CompCol_Matrix(&data->U);
    StatFree(&data->stat);
    free(data); data = NULL;
}

/* Function returning the k-th (k = 1, ..., N - 1) step size                  */
double step(int32_t k) {

    /* Here the step sizes hk = tk - tkm1 (km1 means 'k - 1') are defined.    *
     * The step sizes must be such that their sum is rmax (see                *
     * 'interface.d/settings.c'). Note here the step sizes are set rather     *
     * than the discretisation mesh itself. This is to avoid so-called        *
     * 'catastrophic cancellation' when computing the step sizes, which       *
     * ensures numerical stability.                                           *
     *                                                                        *
     * Note: One must be careful with overflow in integer multiplication and  *
     * addition here, if N is very large.                                     */

    double hk;

    hk = rmax*(2 * k - 1)/((double)(N - 1) * (N - 1));

    return hk;
}

/* Solve eigenproblem                                                         */
void solve(eigensolver_data *data) {

    int ido, nerr, info;
    int32_t dim, n, nl, nev, ncv, ldv, ldz, lworkl, *iparam, *ipntr, k;
    double tol, *resid, *v, *workd, *workl, *select, *d, *z, sigma, runtime,
           *dummy, iC;
    clock_t tstart, tend;

    /* Initialise variables for Lanczos algorithm */
    dim = data->dim;
    ido = 0; n = dim; nl = data->ipar[3]; nev = nmax-nl+1;
    if ((ncv = 2 * nev + 1) < 20) { ncv = 20; }
    if (ncv > dim) { ncv = dim; }
    ldv = dim; ldz = dim; lworkl = ncv * (8 + ncv); info = 0; tol = 1e-12;
    iparam = (int32_t *)calloc(11, sizeof(int32_t));
    ipntr = (int32_t *)calloc(11, sizeof(int32_t));
    resid = (double *)malloc(dim * sizeof(double));
    v = (double *)malloc(dim * ncv * sizeof(double));
    workd = (double *)malloc(3 * dim * sizeof(double));
    workl = (double *)malloc(lworkl * sizeof(double));
    select = (double *)calloc(ncv, sizeof(double)); /* Ritz value ordering */
    d = (double *)malloc(nev * sizeof(double));
    z = (double *)malloc(nev * dim * sizeof(double));
    sigma = shift;
    iparam[0] = 1;
    iparam[2] = 1000000000; /* Large enough to avoid becoming a problem */
    iparam[3] = 1;
    iparam[6] = 3; /* Shift-invert mode */

    /* Iterative calls to 'DSAUPD' */
    tstart = clock(); nerr = 0;
    do {

        /* Call 'DSAUPD' */
        dsaupd_c(&ido, &n, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr,
                 workd, workl, &lworkl, &info); data->info = info;

        /* Check if call was successful */
        if (ido != 1 && ido != -1 && ido != 2 && ido != 99) {
            ERROR("ERROR DURING ITERATION: IDO = %d", ido);
            nerr++;
        }
        if (info != 0 && info != 1) {
            ERROR("ERROR DURING ITERATION: INFO = %d", info);
            nerr++;
        }
        if (info == 1) {
            ERROR("REACHED MAXIMAL NUMBER OF ITERATIONS");
            nerr++;
        }
        if (nerr) { ERROR("DSAUPD ENDED WITH NERR = %d", nerr); }

        /* React to instructions from 'DSAUPD' */
        if (ido == 1) { /* Compute action of shift-inverted Hamiltonian */
            for (k = 0; k < dim; k++) {
                workd[ipntr[1] - 1 + k] = workd[ipntr[2] - 1 + k];
            }
            shift_invert_f(data, &workd[ipntr[1] - 1]); /* Result in argument */
        } else
        if (ido == 2) { /* Compute action of mass matrix */
            mass_matrix_f(data, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
        } else { /* Initialisation step */
            mass_matrix_f(data, &workd[ipntr[0] - 1], &workd[ipntr[1] - 1]);
            shift_invert_f(data, &workd[ipntr[1] - 1]); /* Result in argument */
        }

    } while (ido == 1 || ido == 2 || ido == -1);

    /* Call 'DSEUPD' to extract results */
    dseupd_c(select, d, z, &ldz, &n, &nev, &tol, resid, &ncv, v, &ldv, &sigma,
             iparam, ipntr, workd, workl, &lworkl, &info);
    tend = clock();
    runtime = (tend - tstart) / (double)CLOCKS_PER_SEC;
    data->runtime += runtime;
    printf("ALGORITHM FINISHED SUCCESSFULLY (RUNTIME: %.3f S)\n\n", runtime);

    /* Dummy data to give 'Destroy_Dense_Matrix' something to free */
    dummy = (double *)calloc(1, sizeof(double));
    ((DNformat *)(data->B.Store))->nzval = dummy;

    /* Prepare and save eigenenergies */
    iC = 1. / data->rpar[7];
    for (k = 0; k < nev; k++) { d[k] = iC * (d[k] - offset); }
    save_energies(data, d);

    /* Save radial eigenstates */
    save_states(data, z);

    /* Clean up */
    free(iparam); iparam = NULL;
    free(ipntr); ipntr = NULL;
    free(resid); resid = NULL;
    free(v); v = NULL;
    free(workd); workd = NULL;
    free(workl); workl = NULL;
    free(select); select = NULL;
    free(d); d = NULL;
    free(z); z = NULL;
}

/* Compute action of mass matrix (result stored in y)                         */
void mass_matrix_f(eigensolver_data *data, const double *x, double *y) {

    int32_t k, k0, dim = data->dim;
    double *Mdata = data->Mdata;

    y[0] = Mdata[0] * x[0] + Mdata[1] * x[1];
    for (k = 1; k < dim - 1; k++) {
        k0 = 2 + 3 * (k - 1);
        y[k] = Mdata[k0] * x[k - 1]
             + Mdata[k0 + 1] * x[k]
             + Mdata[k0 + 2] * x[k + 1];
    }
    k0 = 2 + 3 * (dim - 2);
    y[dim - 1] = Mdata[k0] * x[dim - 2] + Mdata[k0 + 1] * x[dim - 1];
}

/* Compute action of shift-inverted Hamiltonian (result stored in x)          */
void shift_invert_f(eigensolver_data *data, double *x) {

    /* Prepare input */
    ((DNformat *)(data->B.Store))->nzval = x;

    /* Solve the system (H-sigma*M) * vout = vin */
    dgstrs(NOTRANS, &data->L, &data->U, data->perm_c, data->perm_r, &data->B,
           &data->stat, &data->info);
}

/* Save computed eigenenergies to file                                        */
void save_energies(eigensolver_data *data, double *energies) {

    char file[71], filename[51];
    int32_t *ipar, nl, lo, jj, runtime, n;
    double EGS, dti, dtf;
    FILE *fd;

    /* Open file for writing */
    nl = (ipar = data->ipar)[3]; lo = ipar[2]; jj = 2 * (int32_t)j + 1;
    EGS = data->rpar[9]; runtime = (int32_t)data->runtime;
    dti = step(1); dtf = step(N - 1);
    (void)sprintf(filename, "energies-%s-%02" PRId32 "-%02" PRId32 ".dat",
                  species, lo, jj);
    (void)strcpy(file, "./data.d/");
    (void)strcat(file, filename);
    if (!(fd = fopen(file, "w"))) {
        ERROR("COULD NOT OPEN FILE '%s' FOR WRITING", filename);
    }

    /* Save metadata */
    (void)fprintf(fd,
                  "EIGENENERGIES FOR '%s' [HARTREE]\n\n"
                  "CPU TIME TO GENERATE DATA SET [S]: %" PRId32 "\n"
                  "GROUND STATE ENERGY [HARTREE]: %1.8lf\n"
                  "ORBITAL ANGULAR MOMENTUM [HBAR]: %c\n"
                  "TOTAL ANGULAR MOMENTUM [HBAR]: %" PRId32 "/2\n"
                  "RMAX [BOHR'S RADIUS]: %1.3E\n"
                  "NUMBER OF DISCRETISATION POINTS: %" PRId32 "\n"
                  "FIRST, FINAL STEP SIZE: %1.3E, %1.3E\n"
                  "MINIMAL PRINCIPAL QUANTUM NUMBER: %" PRId32 "\n"
                  "MAXIMAL PRINCIPAL QUANTUM NUMBER (N): %" PRId32 "\n\n\n\n"
                  "N   ENERGY\n\n",
                  species, runtime, EGS, l, jj, rmax, N, dti, dtf, nl, nmax);

    /* Save eigenenergies */
    n = 0;
    while (++n < nl) { (void)fprintf(fd, "%03" PRId32 "\n", n); }
    while (n++ < nmax + 1) {
        (void)fprintf(fd, "%03" PRId32 " %+1.8E\n", n - 1,
                      energies[n - (nl - 1) - 2]);
    }

    /* Close file */
    fclose(fd); fd = NULL;
}

/* Save computed radial eigenstates to file                                   */
void save_states(eigensolver_data *data, double *z) {

    char file[LEN_PATH_TO_STATES + 101], filename[101];
    int32_t *ipar, nl, lo, jj, dim, n, k;
    FILE *fd;

    /* Open file for writing */
    nl = (ipar = data->ipar)[3]; lo = ipar[2]; jj = 2 * (int32_t)j + 1;
    dim = data->dim;
    for (n = nl; n < nmax + 1; n++) {

        /* Open file for writing */
        file[0] = filename[0] = '\0';
        (void)sprintf(filename,
                      "state-%s-%03" PRId32 "-%02" PRId32 "-%02" PRId32 ".dat",
                      species, n, lo, jj);
        (void)strcpy(file, PATH_TO_STATES);
        (void)strcat(file, filename);

        /* Save metadata */
        if (!(fd = fopen(file, "w"))) {
            ERROR("COULD NOT OPEN FILE '%s' FOR WRITING", filename);
        }
        (void)fprintf(fd,
                      "RADIAL EIGENSTATE FOR '%s' [DIMENSIONLESS]\n\n"
                      "COEFFICIENTS 'FK' (K = 1, ..., N-2)\n"
                      "PRINCIPAL QUANTUM NUMBER (N): %" PRId32 "\n"
                      "ORBITAL ANGULAR MOMENTUM [HBAR]: %c\n"
                      "TOTAL ANGULAR MOMENTUM [HBAR]: %" PRId32 "/2\n"
                      "RMAX [BOHR'S RADIUS]: %1.3E\n"
                      "NUMBER OF DISCRETISATION POINTS: %" PRId32 "\n\n\n\n"
                      "FK\n\n",
                      species, n, l, jj, rmax, N);

        /* Save radial eigenstates */
        for (k = 0; k < dim; k++) {
            (void)fprintf(fd, "%+1.14E\n", z[dim * (n - nl) + k]);
        }

        /* Close file */
        fclose(fd); fd = NULL;
    }
}

/* Save discretisation points and step sizes to file                          */
void save_discretisation() {

    char file[71], filename[51];
    int32_t k;
    double tk, hk;
    FILE *fd;

    /* Check if file already exists */
    (void)sprintf(filename, "discretisation-%s.dat", species);
    (void)strcpy(file, "./data.d/");
    (void)strcat(file, filename);
    if ((fd = fopen(file, "r"))) { return; }

    /* Save metadata */
    (void)fprintf(fd = fopen(file, "w"),
                  "DISCRETISATION DATA FOR SPECIES '%s'\n\n"
                  "NUMBER OF DISCRETISATION POINTS: %" PRId32 "\n\n\n\n"
                  "K        T                     H\n\n",
                  species, N);

    /* Save discretisation data */
    fprintf(fd, "%08" PRId32 " %+1.14E\n", 0, 0.); tk = 0.;
    for (k = 1; k < N; k++) {
        tk += (hk = step(k));
        fprintf(fd, "%08" PRId32 " %+1.14E %+1.14E\n", k, tk, hk);
    }

    /* Close file */
    fclose(fd); fd = NULL;
}
