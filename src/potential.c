/* -------------------------------------------------------------------------- *
 * Parametric model potential V                                               *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * The parametric model potential consists of three terms:                    *
 *                                                                            *
 *     V = VC + VP + VR                                                       *
 *                                                                            *
 * Here, VC is a modified Coulomb potential, VP accounts for the              *
 * polarisability of the effective nucleus, and VR describes relativistic     *
 * spin-orbit coupling. For more information refer to Refs. [4,6,9,10].       *
 *                                                                            *
 * Units.                                                                     *
 *                                                                            *
 *     Charge : e > 0 (elementary charge)                                     *
 *     Mass   : me (electron's mass)                                          *
 *     Length : aB = hbar / (m * c * alpha) (Bohr's radius)                   *
 *     Energy : e**2 / (4 * pi * aB * varepsilon_0) = 27.211386245988(53) eV  *
 *              (Hartree), see Ref. [5].                                      *
 *                                                                            *
 * Each atom/ion species comes with a set of parameters (k1, k2, ...) that    *
 * fix the parametric model potential. In rpar the double precision and in    *
 * ipar the integer parameters are stored.                                    *
 *                                                                            *
 * ipar : integer parameters                                                  *
 * rpar : (real) double precision parameters                                  *
 *                                                                            *
 * For a definition of the parameters see Refs. [6,9].                        *
 * -------------------------------------------------------------------------- */

#include "../inc/eigensolver.h"

/* -------------------------------------------------------------------------- *
 * Parameters for potential                                                   *
 *                                                                            *
 * ipar: integer parameters                                                   *
 * rpar: (real) double precision parameters                                   *
 *                                                                            *
 * Information on parameters                                                  *
 *                                                                            *
 *     ipar[0], Z      : Nuclear charge, units of e > 0                       *
 *     ipar[1], Zc     : Charge of effective nucleus, units of e > 0          *
 *     ipar[2], l      : Orbital angular momentum quantum number (This is the *
 *                       l-value from the global variable settings.)          *
 *     ipar[3], nl     : Minimal principal quantum number for series l (This  *
 *                       is the nl-value from the global variable settings.)  *
 *     rpar[0], k1     : Fitting parameter, see Ref. [1], units of 1 / aB     *
 *     rpar[1], k2     : Fitting parameter, see Ref. [1], units of 1 / aB     *
 *     rpar[2], k3     : Fitting parameter, see Ref. [1], units of 1 / aB     *
 *     rpar[3], k4     : Fitting parameter, see Ref. [1], units of 1 / aB**2  *
 *     rpar[4], rc     : Cut-off radius, see Ref. [1], units of aB            *
 *     rpar[5], alphaD : Polarisability, see Ref. [1], units of               *
 *                       me * e**2 * aB**4 / hbar**2                          *
 *     rpar[6], M      : Total mass of atom/ion, units of me                  *
 *     rpar[7], C      : Mass correction, see theory/theory.pdf               *
 *     rpar[8], j      : Total angular momentum quantum number (This is the   *
 *                       j-value from the global variable settings, but       *
 *                       cleaned up to be an exact half integer.)             *
 *     rpar[9], EGS    : Ground state energy, see Ref. [5], units of Hartree  *
 * -------------------------------------------------------------------------- */
typedef struct potential_data_s {
    int32_t *ipar;
    double *rpar;
} potential_data;

static void move(FILE *, char *);
static double VC(double, const potential_data *);
static double VP(double, const potential_data *);
static double VR(double, const potential_data *, double, double);

/* Initialise parameters (ipar and rpar) for specified atom/ion species       */
void potential_initpar(int32_t *ipar, double *rpar) {

    char *species, id[101];
    int l, c, lread;
    double j, dummy;
    FILE *fd;

    /* Open data file */
    if (!(fd = fopen(SPECIES_DATA, "r"))) {
        ERROR("COULD NOT OPEN INTERFACE/SPECIES.DAT FOR READING");
    }

    /* Search for identifier corresponding to species */
    species = settings.species; id[0] = '*';
    while (strcmp(species, id)) {
        move(fd, id);
        if (!id[0]) {
            ERROR("REQUESTED SPECIES %s IS NOT KNOWN", species);
        }
    }

    /* Get orbital and total angular momentum quantum numbers from settings */
    l = settings.l; j = settings.j;

    /* Read data for atom/ion species */
    if ((c = fgetc(fd)) == '=') {
        ERROR("THERE MUST BE DATA FOR AT LEAST ONE ANGULAR MOMENTUM L");
    }
    do {
        (void)fscanf(fd, "%d %lf %lf %lf %lf %lf %" SCNd32 " ", &lread, rpar,
                     rpar + 1, rpar + 2, rpar + 3, rpar + 4, ipar + 3);

    } while ((c = fgetc(fd)) != '=' && lread != l);
    if (l < lread) { /* Case where it is unclear which parameters to use */
        ERROR("NO DATA FOUND FOR ANGULAR MOMENTUM L = %" PRId32, l);
    }
    if (l > lread) { /* Case where parameters for largest supplied l are used */
        ipar[3] = l + 1;
    }
    while ((c = fgetc(fd)) != '=');
    while ((c = fgetc(fd)) != '\n');
    (void)fscanf(fd, "Z %" SCNd32 " ", ipar);
    (void)fscanf(fd, "ZC %" SCNd32 " ", ipar + 1);
    (void)fscanf(fd, "ALPHAD %lf" " ", rpar + 5);
    (void)fscanf(fd, "M %lf(%lf) ", rpar + 6, &dummy);
    /* IMPORTANT: Here is the position in the code where the mass correction, *
     * i.e., the fact that the reduced mass is NOT the electron's mass, can   *
     * be accounted for. However, here the mass correction is not actually    *
     * included because the employed model parameters (see Refs. [6,8]) are   *
     * computed WITHOUT this correction. This is concluded from the fact that *
     * the computed ground-state energies better fit the ideal ionisation     *
     * energies when the mass correction is omitted. If one employs model     *
     * parameters that include the mass correction, the currently             *
     * commented-out version of the constant C, i.e., the value of rpar[7],   *
     * should be employed. Everything else will be handled automatically.     */
    /* rpar[7] = 1. / (1. + ME / rpar[6]); */
    rpar[7] = 1.;
    ipar[2] = l;
    rpar[8] = .5 * CONVERT(j);
    (void)fscanf(fd, "EGS %lf ", rpar + 9);

    /* Close file */
    fclose(fd); fd = NULL;
}

/* Parametric model potential V(r) (r in units of Bohr's radius)              *
 *                                                                            *
 * The potential V depends on the atom/ion species. To call V for some        *
 * species, e.g., 85RB, one first initialises the associated parameters       *
 * stored in ipar and rpar for that species using the function                *
 * potential_initpar. With ipar and rpar initialised for the desired species, *
 * the potential V at distance r is computed by calling v(r, ipar, rpar).     */
double V(double r, int32_t *ipar, double *rpar) {

    double vc, vp, vr, result;
    potential_data data;

    data.ipar = ipar; data.rpar = rpar;

    vc = VC(r, &data);
    vp = VP(r, &data);
    vr = VR(r, &data, vc, vp);

    result = vc + vp + vr;

    return result;
}

/* -------------------------------------------------------------------------- *
 * Helper functions                                                           *
 * -------------------------------------------------------------------------- */

/* Move filepointer to next entry and get identifier                          */
static void move(FILE *fd, char *id) {

    int c;

    while ((c = fgetc(fd)) != EOF && c != '$');
    if (c != EOF) {
        (void)fgetc(fd);
        (void)fscanf(fd, "ID %s ", id);
        while ((c = fgetc(fd)) != '\n');
        while ((c = fgetc(fd)) != '\n');
        while ((c = fgetc(fd)) != '\n');
    } else {
        id[0] = '\0';
    }
}

/* Modified Coulomb potential in units of Hartree                             */
static double VC(double r, const potential_data *data) {

    int32_t *ipar, Z, Zc;
    double *rpar, k1, k2, k3, k4, Zn, result;

    ipar = data->ipar;
    rpar = data->rpar;

    k1 = rpar[0];
    k2 = rpar[1];
    k3 = rpar[2];
    k4 = rpar[3];

    Z  = ipar[0];
    Zc = ipar[1];

    Zn = Zc + (Z - Zc) * exp(-k1 * r) + (k2 - k4 * r) * r * exp(-k3 * r);

    result = -Zn / r;

    return result;
}

/* Polarisation term in units of Hartree                                      */
static double VP(double r, const potential_data *data) {

    double *rpar, rc, alphaD, result;

    rpar = data->rpar;

    rc      = rpar[4];
    alphaD  = rpar[5];

    result = -.5 * alphaD / pow(r, 4.) * (1. - exp(-pow(r / rc, 6.)));

    return result;
}

/* Relativistic spin-orbit coupling in units of Hartree                       */
static double VR(double r, const potential_data *data, double vc, double vp) {

    int32_t *ipar, Z, Zc, l;
    double *rpar, k1, k2, k3, k4, rc, alphaD, C, j, alphaM, xpnt6, VNR, N,
           aBdZndr, aBdVCdr, aBdVPdr, aBdVNRdr, result;

    ipar = data->ipar;
    rpar = data->rpar;

    k1     = rpar[0];
    k2     = rpar[1];
    k3     = rpar[2];
    k4     = rpar[3];
    rc     = rpar[4];
    alphaD = rpar[5];
    C      = rpar[7];
    j      = rpar[8]; /* Total angular momentum */

    Z      = ipar[0];
    Zc     = ipar[1];
    l      = ipar[2]; /* Orbital angular momentum */

    alphaM   = FC / C;
    xpnt6    = exp(-pow(r / rc, 6.));
    VNR      = vc + vp;
    N        = 1. - .5 * FC * alphaM * VNR; N = N * N;
    aBdZndr  = -(Z - Zc) * k1 * exp(-k1 * r)
             + ((1. - k3 * r) * k2 - (2. - k3 * r) * r * k4) * exp(-k3 * r);
    aBdVCdr  = -(aBdZndr + vc) / r;
    aBdVPdr  = -4. / r * vp - 3. * alphaD / pow(rc, 6.) * r * xpnt6;
    aBdVNRdr = aBdVCdr + aBdVPdr;

    result = .25 * alphaM * alphaM * aBdVNRdr / (r * N)
           * (j * (j + 1.) - l * (l + 1.) - .75);

    return result;
}
