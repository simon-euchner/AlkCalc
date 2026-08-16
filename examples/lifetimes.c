/* -------------------------------------------------------------------------- *
 * Example: Lifetimes                                                         *
 *                                                                            *
 * Compile: gcc -L../lib/ lifetimes.c -lalkcalc -Wl,-rpath,../lib/            *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include <stdio.h>
#include "../interface/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    printf("%s\n", "--- Example: Lifetimes\n");

    /* Define initial and final state */
    int n, dn, l;
    double j, T;
    n = 30; dn = 0; l = 0; j = .5; T = 0.;
    printf("tau = %1.8lf\n\n", alkcalc_tau(T, "40CA+", n, dn, l, j));

    printf("%s\n", "--- End");

    return 0;
}
