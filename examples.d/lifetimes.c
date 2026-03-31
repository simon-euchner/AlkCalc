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
    char *species = "88SR+";
    int n, dn, l;
    double j, T;
    n = 6; dn = 0; l = 1; j = .5; T = 0.;
    printf("tau = %1.8lf\n\n", alkcalc_tau(T, species, n, dn, l, j));

    printf("%s\n", "--- End");

    return 0;
}
