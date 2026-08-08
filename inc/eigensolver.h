/* -------------------------------------------------------------------------- *
 * Main program for computing eigenenergies and radial eigenstates            *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * For more information please see theory/theory.pdf.                         *
 * -------------------------------------------------------------------------- */

#ifndef EIGENSOLVER_H
#define EIGENSOLVER_H

#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "../interface/settings.h"
#include "../GAUSSQ/inc/gaussq.h"
#include "../BSPLINES/inc/bsplines.h"
#include "../EIGLAPACK/inc/eiglapack.h"

#define FC 0.0072973525643 /* Fine-strct. cnst., 0.0072973525643(11) Ref. [5] */
#define ME 0.0005485799090441 /* me, 0.0005485799090441(97) u Ref. [5] */

#define SPECIES_DATA "./interface/species.dat"

/* Macro for error handling */
#define ERROR(...) do { \
    fprintf(stderr, "ERROR (%s:%d): ", __FILE__, __LINE__); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    exit(EXIT_FAILURE); \
} while (0)

/* Macro for converting a half-integer quantum number X to the integer 2 * X */
#define CONVERT(X) (int32_t)floor(2. * (X) + .5)

/* Data type to store data for eigensolver */
typedef struct eigensolver_data_s {
    int32_t Nbs, dim, ipar[4];
    double *M, *H, rpar[10], runtime;
} eigensolver_data;

void potential_initpar(int32_t *, double *);
void validate_settings(int32_t);
double V(double, int32_t *, double *);
eigensolver_data *eigensolver_data_init();
void eigensolver_data_free(eigensolver_data *);
void solve(eigensolver_data *);

#endif
