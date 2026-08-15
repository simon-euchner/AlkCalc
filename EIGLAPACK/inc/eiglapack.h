/* -------------------------------------------------------------------------- *
 * C header file for the eigensolver (LAPACK)                                 *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#ifndef EIGLAPACK_H
#define EIGLAPACK_H

#include <stdint.h>

/* Note: Pass all arguments by reference, because, opposed to C, this is      *
 * default in FORTRAN.                                                        */

/* DSGBVX */
void dsbgvx_c(const int32_t *, const int32_t *, const int32_t *, double *,
              const int32_t *, double *, const int32_t *, double *,
              const int32_t *, const int32_t *, const int32_t *, int32_t *,
              double *, double *, const int32_t *, double *, int32_t *,
              int32_t *, int32_t *);

#endif
