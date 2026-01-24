/* -------------------------------------------------------------------------- *
 * C header file for interface to ARPACK's Lanczos algorithm [3]              *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#ifndef LANCZOS_H
#define LANCZOS_H

#include <stdint.h>
#include <stdbool.h>

/* Note: I pass all arguments by reference, because, opposed to C, this is    *
 *       the default behaviour in FORTRAN.                                    */

/* DSAUPD */
void dsaupd_c(int32_t *, int32_t *, int32_t *, double *, double *, int32_t *,
              double *, int32_t *, int32_t *, int32_t *, double *, double *,
              int32_t *, int32_t *);

/* DSEUPD */
void dseupd_c(double *, double *, double *, int32_t *, int32_t *, int32_t *,
              double *, double *, int32_t *, double *, int32_t *, double *,
              int32_t *, int32_t *, double *, double *, int32_t *, int32_t *);

#endif
