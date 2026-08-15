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

/* Macro for error handling */
#define ERROR(...) do { \
    fprintf(stderr, "ERROR (%s:%d): ", __FILE__, __LINE__); \
    fprintf(stderr, __VA_ARGS__); \
    fprintf(stderr, "\n"); \
    exit(EXIT_FAILURE); \
} while (0)

/* Macro for converting a half-integer quantum number X to the integer 2 * X */
#define CONVERT(X) (int32_t)floor(2. * (X) + .5)

void validate_settings(int32_t);
void potential_initpar(int32_t *, double *);
double V(double, int32_t *, double *);

#endif
