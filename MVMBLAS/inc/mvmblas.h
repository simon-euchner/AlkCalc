/* -------------------------------------------------------------------------- *
 * C header file for banded matrix-vector multiplication                      *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#ifndef MVMBLAS_H
#define MVMBLAS_H

#include <stdint.h>

/* Note: Pass all arguments by reference, because, opposed to C, this is      *
 * default in FORTRAN.                                                        */

/* DSBMV */
void dsbmv_c(const int32_t *, const int32_t *, const double *, const double *,
             const int32_t *, const double *, const double *, double *);

#endif
