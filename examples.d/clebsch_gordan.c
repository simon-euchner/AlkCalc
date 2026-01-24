/* -------------------------------------------------------------------------- *
 * Example: Clebsch-Gordan coefficients                                       *
 *                                                                            *
 * Compile: gcc -L../lib.d/ clebsch_gordan.c -lalkcalc -Wl,-rpath,../lib.d/   *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface.d/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    /* Type that holds information to construct Clebsch-Gordon coefficient */
    alkcalc_cg c;

    printf("%s\n\n", "--- Example: Clebsch-Gordan coefficients");

    /* Clebsch-Gordan coefficient for j1=2, m1=-1, j2=1, m2=0, j=1, mj=-1 */
    c = alkcalc_cj1m1j2m2jmj(2., -1., 1., 0., 1., -1);
    printf("Should be: -1*sqrt(3/10)\n"
           "Is: %d*sqrt(%ld/%ld)\n", c.sign, c.numerator, c.denominator);

    /* Clebsch-Gordan coefficient for j1=j2=1/2, m1=-1/2, m2=1/2, j=0, mj=0 */
    c = alkcalc_cj1m1j2m2jmj(.5, -.5, .5, .5, 0., 0.);
    printf("Should be: -1*sqrt(1/2)\n"
           "Is: %d*sqrt(%ld/%ld)\n", c.sign, c.numerator, c.denominator);

    /* Clebsch-Gordan coefficient for j1=j2=1/2, m1=1/2, m2=-1/2, j=0, mj=0 */
    c = alkcalc_cj1m1j2m2jmj(.5, .5, .5, -.5, 0., 0.);
    printf("Should be: 1*sqrt(1/2)\n"
           "Is: %d*sqrt(%ld/%ld)\n", c.sign, c.numerator, c.denominator);

    /* Clebsch-Gordan coefficient for j1=j2=1/2, m1=1/2, m2=-1/2, j=0, mj=0 */
    c = alkcalc_cj1m1j2m2jmj(1.5, -.5, 1., 1., 1.5, .5);
    printf("Should be: -1*sqrt(8/15)\n"
           "Is: %d*sqrt(%ld/%ld)\n", c.sign, c.numerator, c.denominator);

    printf("\n%s\n", "--- End");

    return 0;
}
