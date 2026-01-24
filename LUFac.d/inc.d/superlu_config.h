
#ifndef SUPERLU_CONFIG_H
#define SUPERLU_CONFIG_H

/* Enable metis */
//#cmakedefine HAVE_METIS @HAVE_METIS@

/* Enable colamd */
//#cmakedefine HAVE_COLAMD @HAVE_COLAMD@

/* enable 64bit index mode */
//#cmakedefine XSDK_INDEX_SIZE @XSDK_INDEX_SIZE@

/* I put this such that it uses the vendor blas and doe not complai that the  *
 * functons 'dmatvec' and alike are missing (see 'dgstrs.c').                 */
#define USE_VENDOR_BLAS

/* Integer type for indexing sparse matrix meta structure */
#if defined(XSDK_INDEX_SIZE) && (XSDK_INDEX_SIZE == 64)
#include <stdint.h>
#define _LONGINT 1
typedef int64_t int_t;
#else
typedef int int_t; /* default */
#endif

#endif /* SUPERLU_CONFIG_H */

