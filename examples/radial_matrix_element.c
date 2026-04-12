/* -------------------------------------------------------------------------- *
 * Example: Radial matrix elements                                            *
 *                                                                            *
 * Compile: gcc -L../lib/ radial_matrix_element.c -lalkcalc                   *
 *          -Wl,-rpath,../lib/                                                *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    printf("%s\n", "--- Example: Radial matrix elements\n");

    /* Compute expectation value of radius in S-states (no mass correction) */
    int nb, lb, nk, lk;
    double jb, jk, p;
    nb = 3; lb = 0; jb = .5;
    nk = nb; lk = lb; jk = .5;
    p = 1.;
    printf("Should be: %1.8lf\n", 1.5 * nb * nk);
    printf("Is: %1.8lf\n\n", alkcalc_rp("1H", nb, lb, jb, p, nk, lk, jk));

    printf("%s\n", "--- End");

    return 0;
}
