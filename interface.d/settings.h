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

#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdint.h>

/* -------------------------------------------------------------------------- *
 * Absolute path to location where radial eigenstates shall be stored (Always *
 * put '/' at the end!)                                                       *
 * -------------------------------------------------------------------------- */
#define PATH_TO_STATES "/home/simon/Files/GitHub/AlkCalc/data.d/"
#define LEN_PATH_TO_STATES 40

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
extern const char *species;
extern const int32_t N;
extern const int32_t nmax;
extern const char l;
extern const double j;
extern const double rmax;
extern const double offset;
extern const double shift;

#endif
