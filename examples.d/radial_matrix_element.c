/* -------------------------------------------------------------------------- *
 * Example: Radial matrix element                                             *
 *                                                                            *
 * Compile: gcc -L../lib.d/ radial_matrix_element.c -lalkcalc                 *
 *          -Wl,-rpath,../lib.d/                                              *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface.d/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    printf("%s\n", "--- Example: Radial matrix element\n");

    /* Compute expectation value of radius in S-states (no mass correction) */
    int nb, lb, nk, lk;
    double jb, jk, p;
    nb = 3; lb = 0; jb = .5;
    nk = nb; lk = lb; jk = .5;
    p = 1.;
    printf("Should be: %1.8lf\n", alkcalc_rp("1H", nb, lb, jb, p, nk, lk, jk));
    printf("Is: %1.8lf\n\n", 1.5*nb*nk);

    printf("%s\n", "--- End");

    return 0;
}
