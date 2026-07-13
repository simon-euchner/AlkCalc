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
 *                                                                            *
 * For convenience, spectroscopic notation (SN) can optionally be used for    *
 * specifying the orbital angular momentum quantum number l (see below)       *
 * within the range 0 <= l <= 20. For l > 20 one must use integers. To use    *
 * SN, write l = SN_S for l = 0, l = SN_P for l = 1, ..., and l = SN_Z for    *
 * l = 20. For the definition of 'SN_...' see 'interface/settings.h'.         *
 * -------------------------------------------------------------------------- */

#include "./settings.h"

/* -------------------------------------------------------------------------- *
 * Settings                                                                   *
 *                                                                            *
 * species   : Atom or ion species                                            *
 * N         : Number of discretisation points (N > 3)                        *
 * nmax      : Maximal principal quantum number                               *
 * l         : Orbital angular momentum quantum number                        *
 * j         : Total angular momentum quantum number                          *
 * rmax      : Maximal radius in units of Bohr's radius, aB                   *
 * offset    : Energy offset to avoid relevant eigenenergies close to zero    *
 * shift     : Shift for shift-invert mode Ref. [3]                           *
 * -------------------------------------------------------------------------- */
const char      *species     =         "1H"                                    ;
const int32_t    N           =      1000000                                    ;
const int32_t    nmax        =            9                                    ;
const int32_t    l           =            0                                    ;
const double     j           =             .5                                  ;
const double     rmax        =        20000.                                   ;
const double     offset      =            1.3                                  ;
const double     shift       =            1.                                   ;
