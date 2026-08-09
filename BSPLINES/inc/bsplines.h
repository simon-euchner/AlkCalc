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

/* DINTRV */
void dintrv_c(double *, int32_t *, double *, int32_t *, int32_t *, int32_t *);

/* DBSPVD */
void dbspvd_c(double *, int32_t *, int32_t *, double *, int32_t *, double *,
              double *);

#endif
