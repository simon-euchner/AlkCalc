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
 * -------------------------------------------------------------------------- */
const alkcalc_settings settings;
settings.species             =                 "1H"                            ;
settings.N                   =                  100                            ;
settings.nmax                =                    9                            ;
settings.l                   =                    0                            ;
settings.j                   =                     .5                          ;
settings.rmax                =                20000.                           ;
