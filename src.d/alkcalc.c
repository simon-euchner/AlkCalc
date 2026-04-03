/* -------------------------------------------------------------------------- *
 * Main source file for AlkCalc                                               *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface.d/alkcalc.h"

static void move(FILE *, int32_t);
static inline double parse(const char *, int32_t);
static alkcalc_cg w3jm(int32_t, int32_t, int32_t, int32_t, int32_t, int32_t);
static int64_t s64imul(int64_t, int64_t);
static int64_t ns64imul(int32_t, int64_t *);
static int64_t s64iadd(int64_t, int64_t);
static int64_t fac(int64_t);
static int64_t euclid(int64_t, int64_t);
static double complex Ylml(int32_t, int32_t, double, double);
static double cgtofloat(alkcalc_cg);
static double thermal_photon_occupation(double, double);
static void nextrm(char *, int32_t *, int32_t *, int32_t, double);

/* -------------------------------------------------------------------------- *
 * Eigenenergy in units of Hartree (27.211386245981(30) eV Ref. [5])          *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * n       : Principal quantum number n = 1, 2, 3, ...                        *
 * l       : Orbital angular momentum l = 0, 1, ..., n-1                      *
 * s       : Spin (Not an argument, since we always have s = 1/2!)            *
 * j       : Total angular momentum quantum number j = |l-1/2|, l+1/2         *
 * -------------------------------------------------------------------------- */
double alkcalc_Enlsj(char *species, int32_t n, int32_t l, double j) {

    char file[LEN_PATH_TO_ALKCALC+101], filename[101];
    int32_t J, nl, nmax, dummy;
    double E;
    FILE *fd;

    /* Open file for reading */
    J = CONVERT(j);
    (void)sprintf(filename, "data.d/energies-%s-%02d-%02d.dat", species, l, J);
    (void)strcpy(file, PATH_TO_ALKCALC);
    (void)strcat(file, filename);
    if (!(fd = fopen(file, "r"))) {
        printf("%s\n", "ERROR: REQUESTED EIGENENERGY NOT AVAILABLE");
        exit(1);
    }

    /* Extract energy */
    move(fd, 9);
    (void)fscanf(fd, "MINIMAL PRINCIPAL QUANTUM NUMBER: %d\n", &nl);
    if (n < nl) {
        printf("%s\n", "ERROR: REQUESTED EIGENENERGY DOES NOT EXIST");
        exit(1);
    }
    (void)fscanf(fd, "MAXIMAL PRINCIPAL QUANTUM NUMBER (N): %d\n", &nmax);
    if (nmax < n) {
        printf("%s\n", "ERROR: REQUESTED EIGENENERGY NOT AVAILABLE");
        exit(1);
    }
    move(fd, n+1);
    (void)fscanf(fd, "%d     %lf\n", &dummy, &E);

    /* Clean up */
    fclose(fd); fd = NULL;

    return E;
}

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
                             double j) {

    char file[LEN_PATH_TO_STATES+101], filename[101], *buffer, *bfr;
    int32_t J, N, dim, dummy, ndi, ndf, a, b, c, d, k;
    double *t, *h, *fnlsj;
    FILE *fd;
    alkcalc_state *state;

    /* Open file with requested state */
    J = CONVERT(j);
    (void)sprintf(filename, "state-%s-%03d-%02d-%02d.dat", species, n, l, J);
    (void)strcpy(file, PATH_TO_STATES);
    (void)strcat(file, filename);
    if (!(fd = fopen(file, "r"))) {
        printf("%s\n", "ERROR: REQUESTED RADIAL EIGENSTATE NOT FOUND");
        exit(1);
    }

    /* Read in metadata and move file pointer to data */
    move(fd, 7);
    (void)fscanf(fd, "NUMBER OF DISCRETISATION POINTS: %d\n", &N);
    move(fd, 2);

    /* Allocate memory for result */
    state = (alkcalc_state *)malloc(sizeof(alkcalc_state));
    state->dim = dim = (state->N = N)-2;
    switch (result) {
        case 'f':
            state->t = t = (double *)malloc(N*sizeof(double));
            state->h = h = (double *)malloc((N-1)*sizeof(double));
            break;
        case 'p':
            state->t = t = NULL;
            state->h = h = NULL;
            break;
        default:
            printf("%s\n", "ERROR: INVALID ARGUMENT FOR 'RESULT'");
            exit(1);
            break;
    }
    state->fnlsj = fnlsj = (double *)malloc(dim*sizeof(double));

    /* Number of digits                                                       *
     *                                                                        *
     * The number of digits used per integer (ndi) and floating-point value   *
     * (ndf) in the data files 'states-...' and 'discretisation-...'. For     *
     * instance, if the floating-point values are of the form '+1.23E+45',    *
     * ndf=3. If the integers, numbering the discretisation points, are of    *
     * the form '0123', ndi=4. The integers a, b, c, and d repeatedly appear  *
     * in the code. They only dependent on ndi and ndf.                       */
    ndi = 8; ndf = 15;
    a = ndf+7; b = dim*a-1; c = ndi+1; d = (N-1)*(ndi+1+2*a)-1;

    /* Read in radial eigenstate */
    (void)fread(buffer = (char *)malloc(b), 1, b, fd);
    for (k=0; k<dim; k++) state->fnlsj[k] = parse(buffer+a*k, ndf);
    free(buffer); buffer = NULL;

    /* Close file */
    fclose(fd); fd = NULL;

    /* Add quantum numbers to result */
    state->n = n; state->l = l; state->j = J/2.;

    /* Read in discretisation data (if requested, i.e., if result = 'f') */
    if (result == 'p') return state;
    file[0] = filename[0] = '\0';
    (void)sprintf(filename, "data.d/discretisation-%s.dat", species);
    (void)strcpy(file, PATH_TO_ALKCALC);
    (void)strcat(file, filename);
    if (!(fd = fopen(file, "r"))) {
        printf("%s\n", "ERROR: REQUESTED DISCRETISATION DATA DOES NOT EXIST");
        exit(1);
    }
    move(fd, 8);
    (void)fscanf(fd, "%d %lf\n", &dummy, t);
    (void)fread(bfr = buffer = (char *)malloc(d), 1, d, fd);
    for (k=0; k<N-2; k++) {
        state->t[k+1] = parse(bfr += c, ndf);
        state->h[k] = parse(bfr += a, ndf);
        bfr += a;
    }
    state->t[N-2+1] = parse(bfr += c, ndf);
    state->h[N-2] = parse(bfr += a, ndf);
    free(buffer); buffer = NULL;

    /* Close file */
    fclose(fd); fd = NULL;

    return state;
}

/* -------------------------------------------------------------------------- *
 * Free for data type 'alkcalc_state'                                         *
 * -------------------------------------------------------------------------- */
void alkcalc_state_free(alkcalc_state *state) {
    free(state->t); state->t = NULL;
    free(state->h); state->h = NULL;
    free(state->fnlsj); state->fnlsj = NULL;
    free(state); state = NULL;
}

/* -------------------------------------------------------------------------- *
 * Radial matrix element <n,l,s,j|r^p|n',l',s',j'> (s=s'=1/2)                 *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * nb      : Principal quantum number of bra                                  *
 * lb      : Orbital angular momentum l = 0, 1, ..., n-1 of bra               *
 * sb      : Spin of bra (Not an argument, since we always have s = 1/2!)     *
 * jb      : Total angular momentum quantum number j = |l-1/2|, l+1/2, of bra *
 * p       : Power of radius operator in matrix element                       *
 * nk      : Principal quantum number of ket                                  *
 * lk      : Orbital angular momentum l = 0, 1, ..., n-1 of ket               *
 * sk      : Spin of ket (Not an argument, since we always have s = 1/2!)     *
 * jk      : Total angular momentum quantum number j = |l-1/2|, l+1/2, of ket *
 * -------------------------------------------------------------------------- */
double alkcalc_rp(char *species, int32_t nb, int32_t lb, double jb, double p,
                  int32_t nk, int32_t lk, double jk) {

    int32_t dim, k;
    double *t, *h, *Rpd, *Rpo, isr3, tkm1, tk, tkp1, hk, hkp1, tm, tp, pm, pp,
           I1, I2, I3, *f, *g, rp;
    alkcalc_state *bra, *ket;

    /* Load states */
    bra = alkcalc_fnlsj('f', species, nb, lb, jb);
    ket = alkcalc_fnlsj('p', species, nk, lk, jk);

    /* Extract discretisation data */
    t = bra->t; h = bra->h;

    /* Allocate memory for matrix Rp */
    Rpd = (double *)malloc((dim = bra->dim)*sizeof(double)); /* Diagonal */
    Rpo = (double *)malloc((dim-1)*sizeof(double)); /* Off-diagonal */

    /* Compute components of matrix Rp (diagonal [d] and off-diagonal [o]) */
    isr3 = 1./sqrt(3.);
    for (k=1; k<dim; k++) {

        /* Discretisation data */
        tkm1 = t[k-1]; tk = t[k]; tkp1 = t[k+1]; hk = h[k-1]; hkp1 = h[k];

        /* Integral for (Rp)k,kp1 */
        tm = .5*(tk+tkp1-isr3*hkp1); tp = .5*(tk+tkp1+isr3*hkp1);
        pm = pow(tm, p); pp = pow(tp, p);
        I1 = .5*((tkp1-tm)*(tm-tk)*pm+(tkp1-tp)*(tp-tk)*pp)/hkp1;
        Rpo[k-1] = I1;

        /* Integrals for (Rp)k,k */
        I2 = .5*((tkp1-tm)*(tkp1-tm)*pm+(tkp1-tp)*(tkp1-tp)*pp)/hkp1;
        tm = .5*(tkm1+tk-isr3*hk); tp = .5*(tkm1+tk+isr3*hk);
        pm = pow(tm, p); pp = pow(tp, p);
        I3 = .5*((tm-tkm1)*(tm-tkm1)*pm+(tp-tkm1)*(tp-tkm1)*pp)/hk;
        Rpd[k-1] = I2 + I3;
    }
    tkm1 = t[dim-1]; tk = t[dim]; tkp1 = t[dim+1]; hk = h[dim-1]; hkp1 = h[dim];
    tm = .5*(tk+tkp1-isr3*hkp1); tp = .5*(tk+tkp1+isr3*hkp1);
    pm = pow(tm, p); pp = pow(tp, p);
    I2 = .5*((tkp1-tm)*(tkp1-tm)*pm+(tkp1-tp)*(tkp1-tp)*pp)/hkp1;
    tm = .5*(tkm1+tk-isr3*hk); tp = .5*(tkm1+tk+isr3*hk);
    pm = pow(tm, p); pp = pow(tp, p);
    I3 = .5*((tm-tkm1)*(tm-tkm1)*pm+(tp-tkm1)*(tp-tkm1)*pp)/hk;
    Rpd[dim-1] = I2 + I3;

    /* Avoid warning by compiler because Rpo[0] is not initialised if dim is  *
     * less than 2, i.e., if N < 4. Practically, this of course not a         *
     * problem, because N must be large, e.g. 100000, to yield sensible       *
     * results. This code just defines the behavior for dim = 1 unambigously. *
     * The case dim < 1 are excluded by the condition N > 2, which is         *
     * presented to the user in 'interface.d/settings.c'.                     */
    f = bra->fnlsj; g = ket->fnlsj;
    if (dim < 2) {
        Rpo[0] = 0.; rp = f[0]*Rpd[0]*g[0];
    } else {
        rp = f[0]*(Rpd[0]*g[0]+Rpo[0]*g[1]);
    }

    /* Compute matrix element */
    for (k=1; k<dim-1; k++)
        rp += f[k]*(Rpo[k-1]*g[k-1]+Rpd[k]*g[k]+Rpo[k]*g[k+1]);
    rp += f[dim-1]*(Rpo[dim-2]*g[dim-2]+Rpd[dim-1]*g[dim-1]);

    /* Clean up */
    alkcalc_state_free(bra); bra = NULL;
    alkcalc_state_free(ket); ket = NULL;
    free(Rpd); Rpd = NULL;
    free(Rpo); Rpo = NULL;

    return rp;
}

/* -------------------------------------------------------------------------- *
 * Clebsch-Gordan coefficients (see 'theory.d/theory.pdf', section 'Manual')  *
 * -------------------------------------------------------------------------- */
alkcalc_cg alkcalc_cj1m1j2m2jmj(double j1, double m1, double j2, double m2,
                                double j, double mj) {

    int32_t J1, M1, J2, M2, J, MJ, gcd;
    alkcalc_cg result;

    /* Convert to half-integers and multiply by 2 */
    J1 = CONVERT(j1); M1 = CONVERT(m1);
    J2 = CONVERT(j2); M2 = CONVERT(m2);
    J = CONVERT(j); MJ = CONVERT(mj);

    /* Wigner's 3jm symbol (no '/2' necessary, see source for 'w3jm' */
    result = w3jm(J1, M1, J2, M2, J, -MJ);

    /* Add phase and scaling factor */
    result.sign *= (((-J1+J2-MJ)/2)%2) ? -1: 1;
    result.numerator = s64imul(result.numerator, J+1);

    /* Clean up result */
    gcd = euclid(result.numerator, result.denominator);
    result.numerator /= gcd; result.denominator /= gcd;

    return result;
}

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
                                double phi) {

    double complex y;
    alkcalc_spinor spinor;

    /* Check input validity */
    if (l < 0 || INTEGER_ABS(ml) > l) {
        printf("%s\n", "ERROR: INVALID ORBITAL ANGULAR MOMENTUM\n");
        exit(1);
    }
    if (theta < 0 || PI < theta) {
        printf("%s\n", "ERROR: INVALID POLAR ANGLE\n");
        exit(1);
    }
    if (phi < 0 || 2.*PI < phi) {
        printf("%s\n", "ERROR: INVALID AZIMUTHAL ANGLE\n");
        exit(1);
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
                                double phi) {

    int32_t ll, J, MJ;
    double cgu, cgd;
    double complex yu, yd;
    alkcalc_spinor spinor;

    /* Check input validity */
    ll = 2*l; J = CONVERT(j); MJ = CONVERT(mj);
    if (J < INTEGER_ABS(MJ)) { /* Ensure mj <= j */
        printf("%s\n", "ERROR: J MUST BE LARGER EQUAL ABSOLUTE VALUE OF MJ\n");
        exit(1);
    }
    if (J%2 != MJ%2) { /* Ensure (half-)int. */
        printf("%s\n", "ERROR: HALF-INTEGER J(MJ) BUT INTEGER MJ(J)\n");
        exit(1);
    }
    if (J != ll-1 && J != ll+1) { /* Ensure |l-1/2| <= j <= l+1/2 */
        printf("%s\n", "ERROR: J MUST BE |L-1/2| OR L+1/2\n");
        exit(1);
    }
    if (theta < 0 || PI < theta) {
        printf("%s\n", "ERROR: INVALID POLAR ANGLE\n");
        exit(1);
    }
    if (phi < 0 || 2.*PI < phi) {
        printf("%s\n", "ERROR: INVALID AZIMUTHAL ANGLE\n");
        exit(1);
    }

    /* Compute Clebsch-Gordan coefficients (spin up [u] and down [d]) */
    cgu = cgtofloat(alkcalc_cj1m1j2m2jmj(l, mj-.5, .5, .5, j, mj));
    cgd = cgtofloat(alkcalc_cj1m1j2m2jmj(l, mj+.5, .5, -.5, j, mj));

    /* Compute spherical harmonics */
    yu = Ylml(l, (MJ-1)/2, theta, phi);
    yd = Ylml(l, (MJ+1)/2, theta, phi);

    /* Assemble result */
    spinor.u = cgu*yu; spinor.d = cgd*yd;

    return spinor;
}

/* -------------------------------------------------------------------------- *
 * Oscillator strength between fine-structure states (dimensionless)          *
 * (see 'theory.d/theory.pdf', section 'Manual')                              *
 *                                                                            *
 * species : String specifying atom/ion species                               *
 * ni      : Principal quantum number of initial state (i)                    *
 * li      : Orbital angular momentum l = 0, 1, ..., n-1, of (i)              *
 * si      : Spin of (i) (Not an argument, since we always have s = 1/2!)     *
 * ji      : Total angular momentum quantum number j = |l-1/2|, l+1/2 of (i)  *
 * nf      : Principal quantum number of final state (f)                      *
 * lf      : Orbital angular momentum l = 0, 1, ..., n-1 of (f)               *
 * sf      : Spin of (f) (Not an argument, since we always have s = 1/2!)     *
 * jf      : Total angular momentum quantum number j = |l-1/2|, l+1/2, of (f) *
 * -------------------------------------------------------------------------- */
double alkcalc_fitof(char *species, int32_t ni, int32_t li, double ji,
                     int32_t nf, int32_t lf, double jf) {

    int32_t JI, JF, llp1, llm1;
    double Efi, r, al, fitof;

    /* Apply selection rules */
    JI = CONVERT(ji); JF = CONVERT(jf);
    if (JI < 0 || JF < 0 || li < 0 || lf < 0) return 0;
    if (INTEGER_ABS(JF-JI) > 2 || INTEGER_ABS(lf-li) != 1) return 0.;

    /* Energy difference between initial (i) and final (f) state in Hartree */
    Efi = alkcalc_Enlsj(species, nf, lf, jf)-alkcalc_Enlsj(species, ni, li, ji);

    /* Radial dipole-transition matrix element between (i) and (f) */
    r = alkcalc_rp(species, ni, li, ji, 1., nf, lf, jf);

    /* Angular factor */
    llp1 = 2*li+1; llm1 = 2*li-1;
    if (lf == li+1) { /* lf-li = 1 */
        if (JI == llp1 && JF == llp1) {
            al = 1./((llp1+2.)*llp1);
        } else
        if (JI == llm1 && JF == llp1) {
            al = (li+1.)/llp1;
        } else
        if (JI == llp1 && JF == llp1+2) {
            al = (li+2.)/(llp1+2);
        } else {
            return 0;
        }
    } else { /* lf-li = -1 */
        if (JI == llm1 && JF == llm1) {
            al = 1./(llp1*llm1);
        } else
        if (JI == llp1 && JF == llm1) {
            al = (double)li/llp1;
        } else
        if (JI == llm1 && JF == llm1-2) {
            al = (li-1.)/llm1;
        } else {
            return 0;
        }
    }

    /* Assemble result */
    fitof = 2./3. * Efi * r*r * al;

    return fitof;
}

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
                   double j) {

    int32_t lp, lm, nmnlp1, nmxlp1, nlp1, nmnlm1, nmxlm1, nlm1, J, k;
    double En, Gamma, jp, jm, hnu, fftoi, nocc, tau;

    /* Get lowest n' such that E(n,l,s,l+s) < E(n',l',s,l'+s) is still true */
    lp = l+1; lm = l-1;
    En = alkcalc_Enlsj(species, n, l, j);
    nextrm(species, &nmnlp1, &nmxlp1, lp, lp+.5); /* l'=l+1 */
    nlp1 = (n < nmnlp1) ? nmnlp1: n;
    if (alkcalc_Enlsj(species, nlp1, lp, lp+.5) > En) {
        while (nmnlp1 < --nlp1 && alkcalc_Enlsj(species, nlp1, lp, lp+.5) > En);
        nlp1++;
    } else {
        while (alkcalc_Enlsj(species, ++nlp1, lp, lp+.5) < En);
    }
    if (!l) { nlm1 = -1; nmnlm1 = 0; goto SkipedSState; } /* l'=l-1 */
    nextrm(species, &nmnlm1, &nmxlm1, lm, lm+.5);
    nlm1 = (n < nmnlm1) ? nmnlm1: n;
    if (alkcalc_Enlsj(species, nlm1, lm, lm+.5) > En) {
        while (nmnlm1 < --nlm1 && alkcalc_Enlsj(species, nlm1, lm, lm+.5) > En);
        nlm1++;
    } else {
        while (alkcalc_Enlsj(species, ++nlm1, lm, lm+.5) < En);
    }
SkipedSState:

    /* Compute decay rate Gamma */
    Gamma = 0.;
    J = CONVERT(j); jp = l+.5; jm = l-.5;
    if (J == 2*l+1) { /* j=l+s */

        /* Emission: l'=l+1 */
        for (k=nlp1-1; k>=nmnlp1; k--) {

            /* j'=l+s */
            hnu = En-alkcalc_Enlsj(species, k, lp, jp);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lp, jp);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu*hnu*fftoi*(1.+nocc);

            /* j'=l+3s */
            hnu = En-alkcalc_Enlsj(species, k, lp, jp+1.);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lp, jp+1.);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu*hnu*fftoi*(1.+nocc);
        }

        /* Emission: l'=l-1 */
        for (k=nlm1-1; k>=nmnlm1; k--) {

            /* j'=l-s */
            hnu = En-alkcalc_Enlsj(species, k, lm, jm);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lm, jm);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu*hnu*fftoi*(1.+nocc);
        }

        /* Absorption: l'=l+1 */
        for (k=nlp1; k<nlp1+dn; k++) {

            /* j'=l+s */
            hnu = alkcalc_Enlsj(species, k, lp, jp)-En;
            fftoi = alkcalc_fitof(species, n, l, j, k, lp, jp);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu*hnu*fftoi*nocc;

            /* j'=l+3s */
            hnu = alkcalc_Enlsj(species, k, lp, jp+1.)-En;
            fftoi = alkcalc_fitof(species, n, l, j, k, lp, jp+1.);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu*hnu*fftoi*nocc;
        }

        /* Absorption: l'=l-1 */
        if (l) {
            for (k=nlm1; k<nlm1+dn; k++) {

                /* j'=l-s */
                hnu = alkcalc_Enlsj(species, k, lm, jm)-En;
                fftoi = alkcalc_fitof(species, n, l, j, k, lm, jm);
                nocc = thermal_photon_occupation(hnu, T);
                Gamma += hnu*hnu*fftoi*nocc;
            }
        }
    } else { /* j=l-s */

        /* Emission: l'=l+1 */
        for (k=nlp1-1; k>=nmnlp1; k--) {

            /* j'=l+s */
            hnu = En-alkcalc_Enlsj(species, k, lp, jp);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lp, jp);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu*hnu*fftoi*(1.+nocc);
        }

        /* Emission: l=l-1 */
        for (k=nlm1-1; k>=nmnlm1; k--) {

            /* j'=l-s */
            hnu = En-alkcalc_Enlsj(species, k, lm, jm);
            fftoi = -alkcalc_fitof(species, n, l, j, k, lm, jm);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu*hnu*fftoi*(1.+nocc);

            /* j'=l-3s */
            if (l > 1) { /* P(j=1/2) -> S(j'=-1/2) is not possible */
                hnu = En-alkcalc_Enlsj(species, k, lm, jm-1.);
                fftoi = -alkcalc_fitof(species, n, l, j, k, lm, jm-1.);
                nocc = thermal_photon_occupation(hnu, T);
                Gamma += hnu*hnu*fftoi*(1.+nocc);
            }
        }

        /* Absorption: l'=l+1 */
        for (k=nlp1; k<nmnlp1+dn; k++) {

            /* j'=l+s */
            hnu = alkcalc_Enlsj(species, k, lp, jp)-En;
            fftoi = alkcalc_fitof(species, n, l, j, k, lp, jp);
            nocc = thermal_photon_occupation(hnu, T);
            Gamma += hnu*hnu*fftoi*nocc;
        }

        /* Absorption: l=l-1 */
        if (l) {
            for (k=nlm1; k<nmnlm1+dn; k++) {

                /* j'=l-s */
                hnu = alkcalc_Enlsj(species, k, lm, jm)-En;
                fftoi = alkcalc_fitof(species, n, l, j, k, lm, jm);
                nocc = thermal_photon_occupation(hnu, T);
                Gamma += hnu*hnu*fftoi*nocc;

                /* j'=l-3s */
                if (l > 1) { /* P(j=1/2) -> S(j'=-1/2) is not possible */
                    hnu = alkcalc_Enlsj(species, k, lm, jm-1.)-En;
                    fftoi = alkcalc_fitof(species, n, l, j, k, lm, jm-1.);
                    nocc = thermal_photon_occupation(hnu, T);
                    Gamma += hnu*hnu*fftoi*nocc;
                }
            }
        }
    }

    /* Compute lifetime in units of nanoseconds                               *
     *                                                                        *
     * The conversion factor used below is 2 x alpha**3 x EH / hbar, where    *
     * alpha is the fine-structure constant, EH is the Hartree, and hbar is   *
     * is the reduced Planck constant; for their values, see Ref. [5].        */
    tau = 1./(32.1300103*Gamma);

    return tau;
}

/* -------------------------------------------------------------------------- *
 * Helper functions                                                           *
 * -------------------------------------------------------------------------- */

/* Move file descriptor down by 'nlines' lines                                */
static void move(FILE *fd, int32_t nlines) {
    int32_t k, c;
    for (k=0; k<nlines; k++) {
        while ((c = fgetc(fd)) != '\n' && c != EOF);
        if (c == EOF) break;
    }
}

/* Optimised parser for reading data files quickly                            */
static inline double parse(const char *str, int32_t nd) {

    int8_t s, i;
    int32_t k;
    uint64_t p10, dec;
    double r;

    /* Get sign */
    s = (str[0] == '-') ? -1: 1;

    /* Leading integer */
    i = str[1]-'0';

    /* Get decimal places */
    p10 = 1; dec = str[3]-'0'; /* First digit */
    for (k=4; k<nd+2; k++) { p10 *= 10; dec = 10*dec + str[k]-'0'; }

    /* Result without power */
    r = s*(i+(double)dec*.1/p10);

    /* Get exponent */
    s = (str[nd+3] == '-') ? -1: 1;
    i = s*((str[nd+4]-'0')*10+(str[nd+5]-'0'));

    /* Add power to result */
    while (i-- > 0) r *= 10.;
    while (++i < 0) r *= .1;

    return r;
}

/* Wigner's 3jm symbols (arguments must be TWICE the desired arguments)       */
static alkcalc_cg w3jm(int32_t j1, int32_t m1, int32_t j2, int32_t m2,
                       int32_t j3, int32_t m3) {

    /* IMPORTANT: The arguments j1, m1, j2, m2, j3, and m3 must be TWICE the  *
     * desired argument, i.e., the following equatility between Wigner's 3jm  *
     * symbols and the function 'w3jm' holds:                                 *
     *                                                                        *
     *        / j1 j2 j3  \                                                   *
     *       |            | = w3jm(2*j1, 2*m1, 2*j2, 2*m2, 2*j3, 2*m3)        *
     *       \  m1 m2 m3 /                                                    *
     *                                                                        *
     * This allows to work with integers throughout this function.            */

    int8_t phase, s;
    int64_t f[7], fs, K, N, g[6], *A, k, l, m, SN, SD;

    /* Prepare result */
    alkcalc_cg result;
    result.sign = 1; result.numerator = 0; result.denominator = 1;

    /* Check if ji and mi (i=1,2,3) are compatible */
    if (   (j1%2 && !(m1%2)) || (m1%2 && !(j1%2))
        || (j2%2 && !(m2%2)) || (m2%2 && !(j2%2))
        || (j3%2 && !(m3%2)) || (m3%2 && !(j3%2))) return result;

    /* Kronecker delta */
    if (m1+m2+m3) return result;

    /* Phase */
    phase = (((j1-j2-m3)/2)%2) ? -1: 1;

    /* Compute factorials in prefactors of sum */
    f[0] = fac((j1+j2-j3)/2);
    f[1] = fac((j1-j2+j3)/2);
    f[2] = fac((-j1+j2+j3)/2);
    f[3] = s64imul(fac((j1-m1)/2), fac((j1+m1)/2));
    f[4] = s64imul(fac((j2-m2)/2), fac((j2+m2)/2));
    f[5] = s64imul(fac((j3-m3)/2), fac((j3+m3)/2));
    f[6] = fac((j1+j2+j3)/2+1); /* Used later */
    if (!(fs = ns64imul(6, f))) return result;

    /* Bounds for summation */
    K = MAX(0, MAX((j2-j3-m1)/2, (j1-j3+m2)/2));
    N = MIN((j1+j2-j3)/2, MIN((j1-m1)/2, (j2+m2)/2));

    /* Summation */
    if (N < K) return result;
    A = (int64_t *)malloc((N-K+1)*sizeof(int64_t));
    for (k=K; k<N+1; k++) {
        s = (k%2) ? -1: 1;
        g[0] = fac(k);
        g[1] = fac((j1+j2-j3)/2-k);
        g[2] = fac((j1-m1)/2-k);
        g[3] = fac((j2+m2)/2-k);
        g[4] = fac((j3-j2+m1)/2+k);
        g[5] = fac((j3-j1-m2)/2+k);
        A[k-K] = s64imul(s, ns64imul(6, g)); /* Denominators of summands */
    }
    SN = 0; /* Numerator of sum */
    for (k=0; k<N-K+1; k++) {
        m = 1;
        for (l=0; l<N-K+1; l++) {
            if (l == k) continue;
            m = s64imul(m, A[l]);
        }
        SN = s64iadd(SN, m);
    }
    SD = 1; /* Denominator of sum */
    for (k=0; k<N-K+1; k++) SD = s64imul(SD, A[k]);
    free(A);

    /* Assemble result */
    phase *= (SN < 0) ? -1: 1; phase *= (SD < 0) ? -1: 1;
    result.sign = phase;
    result.numerator = s64imul(fs, s64imul(SN, SN));
    result.denominator = s64imul(f[6], s64imul(SD, SD));

    return result;
}

/* Secure 64-bit integer multiplication                                       */
static int64_t s64imul(int64_t a, int64_t b) {
    if (!a || !b) return 0;
    if (a == INT64_MIN) {
        if (b == 1) { return a; } else { goto s64imulOverflow; }
    }
    if (b == INT64_MIN) {
        if (a == 1) { return b; } else { goto s64imulOverflow; }
    }
    if (   (a > 0 && b > 0 && a <= INT64_MAX/b)
        || (a < 0 && b < 0 && -a <= INT64_MAX/(-b))
        || (a > 0 && b < 0 && -a >= INT64_MIN/(-b))
        || (a < 0 && b > 0 && -b >= INT64_MIN/(-a)))
        return a*b;
s64imulOverflow:
    printf("%s\n", "ERROR: OVERFLOW IN INTEGER MULTIPLICATION");
    exit(1);
}

/* Secure 64-bit integer multiplication (n times)                             */
static int64_t ns64imul(int32_t n, int64_t *a) {
    int32_t k; int64_t result = 1;
    for (k=0; k<n; result = s64imul(result, a[k++]));
    return result;
}

/* Secure 64-bit integer addition                                             */
static int64_t s64iadd(int64_t a, int64_t b) {

    /* It might look dangerous to do 'INT64_MIN-a' when 'a' can be equal to   *
     * INT64_MIN. However, in the particular order we perform the             *
     * subtraction, the C99 standard guarantees that the expression evaluates *
     * to zero (see Sec. 6.5.5 in Ref. [11]).                                 */

    if (a >= 0 && b <= INT64_MAX-a) return a+b;
    if (a < 0 && a >= INT64_MIN && b >= INT64_MIN-a) return a+b;
    printf("%s\n", "ERROR: OVERFLOW IN INTEGER ADDITION");
    exit(1);
}

/* Integer factorial                                                          */
static int64_t fac(int64_t n) {
    if (n < 0) return 0;
    int64_t l, m = 1;
    for (l=0; l<n-1; l++) m = s64imul(m, n-l);
    return m;
}

/* Euclidean algorithm                                                        */
static int64_t euclid(int64_t a, int64_t b) {
    int64_t c;
    while (b != 0) { c = a%b; a = b; b = c; }
    return a;
}

/* Spherical harmonics (see definition in 'theory.d/theory.pdf')              */
static double complex Ylml(int32_t l, int32_t ml, double theta, double phi) {

    int8_t sml, phase;
    int32_t k;
    double pf, x, Pk, Pkm1, Pkm2;
    double complex ac, y;

    /* Azimuthal contribution (here we still need the sign of ml) */
    ac = COMPLEX(cos(ml*phi), sin(ml*phi));

    /* Check input regime and react accordingly */
    sml = (ml < 0) ? -1: 1; ml *= sml;
    if (l < ml) return COMPLEX(0., 0.);

    /* Prefactor (the Condon-Shortley phase is in Legendre polynomials) */
    pf = sqrt((2*l+1)*fac(l-ml)/(4.*PI*fac(l+ml)));

    /* Polar contribution                                                     *
     *                                                                        *
     * The associated Legendre polynomial's value is computed for the         *
     * absolute value of ml. The conversion formula, described in             *
     * 'theory.d/theory.pdf', allows to obtain the value for negative ml, and *
     * the conversion factor is not just a phase. Please note that this       *
     * factor is not missing here, but included in the prefactor, pf,         *
     * already, because at the point where pf is computed, ml is already      *
     * rendered non-negative. This means automatically everything (up to a    *
     * phase included later) is correct. It is probably vital to view         *
     * 'theory.d/theory.pdf' to understand this part.                         */
    phase = (ml%2) ? -1: 1; x = cos(theta);
    if (l == ml) {
        Pk = phase*fac(2*ml)/fac(ml)*pow(.25*(1.-x*x), .5*ml);
    } else
    if (l == ml+1) {
        Pkm1 = phase*fac(2*ml)/fac(ml)*pow(.25*(1.-x*x), .5*ml);
        Pk = (2*(ml+1)-1)*x*Pkm1;
    } else {
        Pkm2 = phase*fac(2*ml)/fac(ml)*pow(.25*(1.-x*x), .5*ml);
        Pkm1 = (2*ml+1)*x*Pkm2; k = ml+2;
        do {
            Pk = ((2*k-1)*x*Pkm1-(k+ml-1)*Pkm2)/(k-ml);
            Pkm2 = Pkm1; Pkm1 = Pk;
        } while (k++ < l);
    }

    /* Assemble result (here we also include the phase, as noted above) */
    y = (sml < 0) ? phase: 1.; /* Phase from 'ml -> -ml' conversion */
    y *= pf * Pk * ac;

    return y;
}

/* Convert symbolic Clebsch-Gordan coefficient into a floating-point number */
static double cgtofloat(alkcalc_cg c) {
    return c.sign*sqrt(c.numerator/(double)c.denominator);
}

/* Thermal photon-occupation number at energy hnu and temperature T */
static double thermal_photon_occupation(double hnu, double T) {

    double r, x;

    /* Ratio of photon and thermal energy                                     *
     *                                                                        *
     * - Photon energy h x nu (hnu) in units of Hartree (EH)                  *
     * - Temperature T in units of Kelvin (K)                                 */
    r = hnu/(3.166811e-6*T); /* For Boltzmann's constant see Ref. [5] */

    /* Compute photon occupation number according to Planck's law */
    if (T <= 0.) return 0.; /* Zero T case */
    if (r < cbrt(720.*DBL_EPSILON)) return 1./r-.5+1./12.*r; /* High T */
    if (r > -1./3.*log(DBL_EPSILON)) { x = exp(-r); return x+x*x; } /* Low T */
    return 1./expm1(r); /* Intermediate regime */
}

/* Fetch extremal principle quantum numbers */
static void nextrm(char *species, int32_t *nmin, int32_t *nmax, int32_t l,
                   double j) {

    char file[LEN_PATH_TO_ALKCALC+101], filename[101];
    int32_t J;
    FILE *fd;

    /* Open file for reading */
    J = CONVERT(j);
    (void)sprintf(filename, "data.d/energies-%s-%02d-%02d.dat", species, l, J);
    (void)strcpy(file, PATH_TO_ALKCALC);
    (void)strcat(file, filename);
    if (!(fd = fopen(file, "r"))) {
        printf("%s\n", "ERROR: REQUESTED ENERGY SERIES IS NOT AVAILABLE");
        exit(1);
    }

    /* Extract information */
    move(fd, 9);
    (void)fscanf(fd, "MINIMAL PRINCIPAL QUANTUM NUMBER: %d\n", nmin);
    (void)fscanf(fd, "MAXIMAL PRINCIPAL QUANTUM NUMBER (N): %d\n", nmax);

    /* Clean up */
    fclose(fd); fd = NULL;
}
