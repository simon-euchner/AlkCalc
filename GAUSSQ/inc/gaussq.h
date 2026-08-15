/* -------------------------------------------------------------------------- *
 * C header file for the interface to Gaussian quadratures                    *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#ifndef GAUSSQ_H
#define GAUSSQ_H

#include <stdint.h>

/* Note: Pass all arguments by reference, because, opposed to C, this is      *
 * default in FORTRAN.                                                        */

/* GAUSSQ */
void gaussq_c(const int32_t *, double *, double *);

#endif
