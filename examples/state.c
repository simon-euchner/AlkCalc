/* -------------------------------------------------------------------------- *
 * Example: Evaluate radial eigenfunction fnslj                               *
 *                                                                            *
 * Compile: gcc -L../lib/ state.c -lalkcalc -Wl,-rpath,../lib/                *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include <stdio.h>
#include "../interface/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    printf("%s\n", "--- Example: Evaluate radial eigenfunction fnslj\n");

    /* Compute radial eigenfunction fnlsj at t = 0, 1, 2, 3, 4, 5 */
    int i, ltevals, n, l;
    double j, tevals_cpy[6], tevals[] = {0., 1., 2., 3., 4., 5.};

    /* Copy array tevals */
    for (i = 0; i < 6; i++) { tevals_cpy[i] = tevals[i]; }

    /* Evaluate fnlsj at t in tevals */
    n = 4, l = 0, j = .5; ltevals = 6;
    alkcalc_fnslj_eval("40CA+", n, l, j, tevals, ltevals);

    /* Print result */
    printf("%s\n", "t    fnlsj(t)");
    for (i = 0; i < 6; i++) {
        printf("%1.2f %+1.5E\n", tevals_cpy[i], tevals[i]);
    }
    printf("%c", '\n');

    printf("%s\n", "--- End");

    return 0;
}
