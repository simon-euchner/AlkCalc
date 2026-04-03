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

#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <time.h>
#include "../interface.d/settings.h"
#include "../LUFac.d/inc.d/slu_ddefs.h"
#include "../LANCZOS.D/inc.d/lanczos.h"

#define FC 0.0072973525643 /* Fine-structure, 0.0072973525643(11) Ref. [5] */
#define ME 0.0005485799090441 /* me, 0.0005485799090441(97) u Ref. [5] */

#define SPECIES_DATA "./interface.d/species.dat"

/* Macro for error handling */
#define ERROR(...) do { \
    fprintf(stderr, "ERROR (%s:%d): ", __FILE__, __LINE__); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    exit(EXIT_FAILURE); \
} while (0)

/* Data type to store data for eigensolver */
typedef struct eigensolver_data_s {
    int32_t dim, *perm_r, *perm_c, ipar[4], info;
    double *Mdata, rpar[10], runtime;
    SuperMatrix B, L, U;
    SuperLUStat_t stat;
} eigensolver_data;

void vint_initpar(double *, int32_t *);
double vint(double, double *, int32_t *);
eigensolver_data *eigensolver_data_init();
void eigensolver_data_free(eigensolver_data *);
double step(int32_t);
void solve(eigensolver_data *);
void mass_matrix_f(eigensolver_data *, const double *, double *);
void shift_invert_f(eigensolver_data *, double *);
void save_energies(eigensolver_data *, double *);
void save_states(eigensolver_data *, double *);
void save_discretisation();

#endif
