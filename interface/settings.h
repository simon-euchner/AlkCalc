/* -------------------------------------------------------------------------- *
 * Settings                                                                   *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- *
 *                                                                            *
 * The purpose of this file is to provide an interface for the user to        *
 * specify settings for computing eigenenergies and radial eigenstates.       *
 *                                                                            *
 * Note that by running the main Makefile with the argument solve, i.e.,      *
 * running the command 'make solve', this file is read, the settings applied, *
 * all specified eigenenergies and radial eigenstates computed, and the       *
 * results saved.                                                             *
 *                                                                            *
 * Note further that this file has no effect on the functions in alkcalc.h.   *
 * These merely use the already generated data.                               *
 * -------------------------------------------------------------------------- */

#ifndef SETTINGS_H
#define SETTINGS_H

#include <stdint.h>

/* -------------------------------------------------------------------------- *
 * Absolute path to location where radial eigenstates shall be stored (Always *
 * put '/' at the end!)                                                       *
 * -------------------------------------------------------------------------- */
#define PATH_TO_STATES "/home/simon/Files/GitHub/AlkCalc/data/"
#define LEN_PATH_TO_STATES 38

/* -------------------------------------------------------------------------- *
 * Definitions for spectroscopic notation (SN)                                *
 * -------------------------------------------------------------------------- */
#define SN_S  0
#define SN_P  1
#define SN_D  2
#define SN_F  3
#define SN_G  4
#define SN_H  5
#define SN_I  6
#define SN_K  7
#define SN_L  8
#define SN_M  9
#define SN_N 10
#define SN_O 11
#define SN_Q 12
#define SN_R 13
#define SN_T 14
#define SN_U 15
#define SN_V 16
#define SN_W 17
#define SN_X 18
#define SN_Y 19
#define SN_Z 20

/* Data type for settings */
typedef struct alkcalc_settings_s {
    char *species;
    int32_t k, N, nmax, l;
    double j, rmax;
} alkcalc_settings;

/* -------------------------------------------------------------------------- *
 * Settings                                                                   *
 *                                                                            *
 * species   : Atom or ion species                                            *
 * k         : Order of B-splines (k >= 3)                                    *
 * N         : Number of discretisation points (N >= k + 1)                   *
 * nmax      : Maximal principle quantum number                               *
 * l         : Orbital angular momentum quantum number                        *
 * j         : Total angular momentum quantum number                          *
 * rmax      : Maximal radius in units of Bohr's radius, aB                   *
 * -------------------------------------------------------------------------- */
extern const alkcalc_settings settings;

#endif
