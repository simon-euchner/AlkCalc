/* -------------------------------------------------------------------------- *
 * Example: Lifetimes                                                         *
 *                                                                            *
 * Compile: gcc -L../lib.d/ lifetimes.c -lalkcalc -Wl,-rpath,../lib.d/        *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface.d/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    printf("%s\n", "--- Example: Lifetimes\n");

    /* Define initial and final state */
    int n, dn, l;
    double j;
    n = 20; dn = 0; l = 0; j = .5;
    printf("tau = %1.8lf\n\n", alkcalc_tau(0., "88SR+", n, dn, l, j));

    printf("%s\n", "--- End");

    return 0;
}
