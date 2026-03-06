/* -------------------------------------------------------------------------- *
 * Example: Oscillator strengths                                              *
 *                                                                            *
 * Compile: gcc -L../lib.d/ oscillator_strengths.c -lalkcalc                  *
 *          -Wl,-rpath,../lib.d/                                              *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface.d/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    printf("%s\n", "--- Example: Oscillator strengths\n");

    /* Define initial and final state */
    int ni, li, nf, lf;
    double ji, mji, jf, mjf;
    ni = 3; li = 0; ji = .5; mji = .5;
    nf = nf; li = 1; jf = .5; mjf = .5;
    printf("fitof = %1.8lf\n\n", alkcalc_fitof);

    printf("%s\n", "--- End");

    return 0;
}
