/* -------------------------------------------------------------------------- *
 * Settings                                                                   *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * The purpose of this file is to provide an interface for the user to        *
 * specify settings for computing eigenenergies and radial eigenstates.       *
 *                                                                            *
 * Note that by running the main 'Makefile' with the argument 'solve', i.e.,  *
 * running the command 'make solve', this file is read, the settings applied, *
 * all specified eigenenergies and radial eigenstates computed, and the       *
 * results saved.                                                             *
 *                                                                            *
 * Note further that this file has no effect on the functions in 'alkcalc.h'. *
 * These merely use the already generated data.                               *
 * -------------------------------------------------------------------------- */

#include "./settings.h"

/* -------------------------------------------------------------------------- *
 * Settings                                                                   *
 *                                                                            *
 * species   : Atom or ion species                                            *
 * N         : Number of discretisation points (N > 2)                        *
 * nmax      : Maximal principle quantum number                               *
 * l         : Orbital angular momentum quantum number                        *
 * j         : Total angular momentum quantum number                          *
 * rmax      : Maximal radius in units of Bohr's radius, aB                   *
 * offset    : Energy offset to avoid relevant eigenenergies close to zero    *
 * shift     : Shift for shift-invert mode Ref. [3]                           *
 * -------------------------------------------------------------------------- */
const char      *species     =         "1H"                                    ;
const int32_t    N           =       100000                                    ;
const int32_t    nmax        =           70                                    ;
const char       l           =          'S'                                    ;
const double     j           =             .5                                  ;
const double     rmax        =        10000.                                   ;
const double     offset      =            1.6                                  ;
const double     shift       =            1.                                   ;
