/* -------------------------------------------------------------------------- *
 * Example: Coupled basis states                                              *
 *                                                                            *
 * Compile: gcc -L../lib/ coupled_basis.c -lalkcalc -Wl,-rpath,../lib/        *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    /* Orbital angular momentum */
    int l;
    double j, mj;
    l = 1;
    j = 1.5;
    mj = -.5;

    /* Fix some angles */
    double theta, phi;
    theta = .736;
    phi = 3.57;

    printf("%s\n", "--- Example: Coupled basis states\n\n(a, b) = a + i*b\n");

    /* Compute spinor components */
    alkcalc_spinor s;
    s = alkcalc_Philsjmj(l, j, mj, theta, phi);
    printf("Phi(l=P,s=1/2,j=3/2,ml=-1/2)_0 = (%+lf, %+lf)\n"
           "Phi(l=P,s=1/2,j=3/2,ml=-1/2)_1 = (%+lf, %+lf)\n",
           creal(s.u), cimag(s.u), creal(s.d), cimag(s.d));

    printf("\n%s\n", "--- End");

    return 0;
}
