/* -------------------------------------------------------------------------- *
 * C header file for the interface to de Boor's B-splines                     *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#ifndef SPLINES_H
#define SPLINES_H

#include <stdint.h>

/* Note: Pass all arguments by reference, because, opposed to C, this is      *
 * default in FORTRAN.                                                        */

/* DBSPVD */
void dbspvd_c(const double *, const int32_t *, const int32_t *, const double *,
              const int32_t *, double *, double *);

/* DINTRV */
void dintrv_c(const double *, const int32_t *, const double *, int32_t *,
              int32_t *, int32_t *);

#endif
