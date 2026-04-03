/* -------------------------------------------------------------------------- *
 * Example: Eigenenergies                                                     *
 *                                                                            *
 * Compile: gcc -L../lib.d/ eigenenergies.c -lalkcalc -Wl,-rpath,../lib.d/    *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface.d/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    printf("%s\n\n", "--- Example: Eigenenergies");

    /* Print lowest five eigenenergies of Hydrogen (no mass correction) in    *
     * units of Hartree                                                       */
    printf("%s\n\n", "n       Energy             Energy (exact)");
    double E;
    for (int n = 1; n <= 5; n++) {
        E = alkcalc_Enlsj("1H", n, 0, .5);
        printf("%d       %lf       %lf\n", n, E, -.5 / (n * n));
    }

    printf("\n%s\n", "--- End");

    return 0;
}
