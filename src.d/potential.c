/* -------------------------------------------------------------------------- *
 * Model for interaction potential for a field free atom/ion in relative      *
 * co-ordinates                                                               *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * The model potential consists of three terms:                               *
 *                                                                            *
 *     Vint = VC + VP + VR                                                    *
 *                                                                            *
 * Here, VC is a modified Coulomb potential, VP accounts for polarisation of  *
 * the screened ionic core due to the valence electron, and VR is the         *
 * relativistic spin-orbit coupling. For more information refer to            *
 * Refs. [4,6,9,10] and references therein.                                   *
 *                                                                            *
 * Units.                                                                     *
 *                                                                            *
 *     Charge : e > 0 (elementary charge)                                     *
 *     Mass   : me (electron's mass)                                          *
 *     Length : aB = hbar / (m * c * alpha) (Bohr's radius)                   *
 *     Energy : e**2 / (4pi * aB * varepsilon_0) = 27.211386245988(53) eV     *
 *              (Hartree), see Ref. [5].                                      *
 *                                                                            *
 * Each atom/ion species comes with a set of parameters. In 'rpar' the double *
 * precision and in 'ipar' the integer ones are stored.                       *
 *                                                                            *
 * rpar : (real) double precision parameters                                  *
 * ipar : integer parameters                                                  *
 *                                                                            *
 * For a definition of the parameters see Refs. [6,9].                        *
 * -------------------------------------------------------------------------- */

#include "../inc.d/eigensolver.h"

/* -------------------------------------------------------------------------- *
 * Parameters for potential                                                   *
 *                                                                            *
 * rpar: double precision parameters                                          *
 * ipar: integer parameters                                                   *
 *                                                                            *
 * Information on parameters                                                  *
 *                                                                            *
 *     rpar[0], k1     : Fitting parameter, see Ref. [1], units of 1 / aB     *
 *     rpar[1], k2     : Fitting parameter, see Ref. [1], units of 1 / aB     *
 *     rpar[2], k3     : Fitting parameter, see Ref. [1], units of 1 / aB     *
 *     rpar[3], k4     : Fitting parameter, see Ref. [1], units of 1 / aB**2  *
 *     rpar[4], rc     : Cut-off radius, see Ref. [1], units of aB            *
 *     rpar[5], alphaD : Polarisability, see Ref. [1], units of               *
 *                       me * e**2 * aB**4 / hbar**2                          *
 *     rpar[6], M      : Total mass of atom/ion, units of me                  *
 *     rpar[7], C      : Mass correction, see 'theory.d/theory.pdf'           *
 *     rpar[8], j      : Total angular momentum quantum number (version of    *
 *                       globally defined j [see 'interface.d/settings.c']    *
 *                       but corrected to be exactly half integer).           *
 *     rpar[9], EGS    : Ground state energy, see Ref. [5], units of Hartree  *
 *                                                                            *
 *     ipar[0], Z      : Nuclear charge, units of e > 0                       *
 *     ipar[1], Zc     : Charge of screened core, units of e > 0              *
 *     ipar[2], l      : Orbital angular momentum quantum number              *
 *     ipar[3], nl     : Minimal principal quantum number for series 'l'      *
 * -------------------------------------------------------------------------- */
typedef struct _VintData {
    double *rpar;
    int32_t *ipar;
} vint_data;

static void move(FILE *, char *);
static double VC(double, vint_data *);
static double VP(double, vint_data *);
static double VR(double, vint_data *, double, double);

/* Initialise parameters (rpar, ipar), depending on atom/ion species          */
void vint_initpar(double *rpar, int32_t *ipar) {

    char id[101];
    int c;
    double dummy;
    FILE *fd;

    /* Open data file */
    if (!(fd = fopen(SPECIES_DATA, "r"))) {
        ERROR("COULD NOT OPEN 'SPECIES.DAT' FOR READING");
    }

    /* Search for identifier corresponding to 'species' */
    id[0] = '*';
    while (strcmp(species, id)) {
        move(fd, id);
        if (!id[0]) {
            ERROR("REQUESTED SPECIES '%s' IS NOT KNOWN", species);
        }
    }

    /* Read data for atom/ion species */
    while (l != (c = fgetc(fd)) && c != 'Z')
        while ((c = fgetc(fd)) != '\n');
    if (c == 'Z') {
        ERROR("REQUESTED QUANTUM NUMBER 'L = %c' IS NOT KNOWN", l);
    } else {
        (void)fgetc(fd);
        (void)fscanf(fd, "%lf %lf %lf %lf %lf %" SCNd32 " ", rpar, rpar + 1,
                     rpar + 2, rpar + 3, rpar + 4, ipar + 3);
    }
    while ((c = fgetc(fd)) != 'Z');
    (void)fscanf(fd, "%" SCNd32 " ", ipar);
    (void)fscanf(fd, "ZC %" SCNd32 " ", ipar + 1);
    (void)fscanf(fd, "ALPHAD %lf" " ", rpar + 5);
    (void)fscanf(fd, "M %lf(%lf) ", rpar + 6, &dummy);
    /* IMPORTANT: Here is the position in the code where the mass correction, *
     * i.e., the fact that the reduced mass is NOT the electron's mass, can   *
     * be accounted for. However, we do not actually include the mass         *
     * correction here, because the model parameters that we use (see         *
     * Refs. [6,8]) are computed WITHOUT this correction. This we conclude    *
     * from the fact that the computed ground state energies better fit the   *
     * ideal ionisation energy, if we omit the mass correction. If one        *
     * employs model parameters that include the mass correction, the         *
     * currently commented-out version for the constant C, i.e., the value of *
     * rpar[7], should be employed. Everything else will be handled           *
     * automatically.                                                         */
    /* rpar[7] = 1. / (1. + ME / rpar[6]); */
    rpar[7] = 1.;
    switch (l) {
        case 'S': ipar[2] =  0; break;
        case 'P': ipar[2] =  1; break;
        case 'D': ipar[2] =  2; break;
        case 'F': ipar[2] =  3; break;
        case 'G': ipar[2] =  4; break;
        case 'H': ipar[2] =  5; break;
        default: break;
    }
    rpar[8] = .5 * (2 * (int32_t)j + 1);
    (void)fscanf(fd, "EGS %lf ", rpar + 9);

    /* Close file */
    fclose(fd); fd = NULL;
}

/* Move filepointer to next dollar sign and get identifier                    */
static void move(FILE *fd, char *id) {

    int c;

    while ((c = fgetc(fd)) != EOF && c != '$');
    if (c != EOF) {
        (void)fgetc(fd);
        (void)fscanf(fd, "ID %s ", id);
        while ((c = fgetc(fd)) != '\n');
        while ((c = fgetc(fd)) != '\n');
    } else {
        id[0] = '\0';
    }
}

/* Interaction potential                                                      *
 * To call this function, first select an atom/ion species by initialising    *
 * 'rpar' and 'ipar' with 'vint_initpar'; argument in units of Bohr's radius  */
double vint(double r, double *rpar, int32_t *ipar) {

    double vc, vp, result;
    vint_data data;

    data.rpar = rpar; data.ipar = ipar;

    vc = VC(r, &data);
    vp = VP(r, &data);

    result = vc + vp + VR(r, &data, vc, vp);

    return result;
}

/* Modified Coulomb's potential in units of Hartree                           */
static double VC(double r, vint_data *data) {

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
static double VP(double r, vint_data *data) {

    double *rpar, rc, alphaD, result;

    rpar = data->rpar;

    rc      = rpar[4];
    alphaD  = rpar[5];

    result = -.5 * alphaD / pow(r, 4.) * (1. - exp(-pow(r / rc, 6.)));

    return result;
}

/* Relativistic spin-orbit coupling in units of Hartree                       */
static double VR(double r, vint_data *data, double vc, double vp) {

    /* In 'theory.d/theory.pdf' the variable K is called N. Here, it is named *
     * K to avoid clash with the global variables in 'interface.d/settings.h'.*/

    int32_t *ipar, Z, Zc, lo;
    double *rpar, k1, k2, k3, k4, rc, alphaD, C, jt, alphaM, xpnt6, VNR, K,
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
    jt     = rpar[8]; /* Total angular momentum */

    Z      = ipar[0];
    Zc     = ipar[1];
    lo     = ipar[2]; /* Orbital angular momentum */

    alphaM   = FC / C;
    xpnt6    = exp(-pow(r / rc, 6.));
    VNR      = vc + vp;
    K        = 1. - .5 * FC * alphaM * VNR; K = K * K;
    aBdZndr  = -(Z - Zc) * k1 * exp(-k1 * r)
             + ((1. - k3 * r) * k2 - (2. - k3 * r) * r * k4) * exp(-k3 * r);
    aBdVCdr  = -(aBdZndr + vc) / r;
    aBdVPdr  = -4. / r * vp - 3. * alphaD / pow(rc, 6.) * r * xpnt6;
    aBdVNRdr = aBdVCdr + aBdVPdr;

    result = .25 * alphaM * alphaM * aBdVNRdr / (r*K)
           * (jt * (jt + 1.) - lo * (lo + 1.) - .75);

    return result;
}
