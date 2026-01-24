/* -------------------------------------------------------------------------- *
 * Example: Uncoupled basis states                                            *
 *                                                                            *
 * Compile: gcc -L../lib.d/ uncoupled_basis.c -lalkcalc -Wl,-rpath,../lib.d/  *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface.d/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    /* Orbital angular momentum */
    int l, ml;
    l = 1;

    /* Fix some angles */
    double theta, phi;
    theta = .736;
    phi = 3.57;

    printf("%s\n", "--- Example: Uncoupled basis states\n\n(a, b) = a + i*b\n");

    /* Compute and print spinors for all possible ml and ms */
    alkcalc_spinor yp, ym;
    for (ml=-l; ml<=l; ml++) {
        yp = alkcalc_YlmlXsms(l, ml, .5, theta, phi);
        ym = alkcalc_YlmlXsms(l, ml, -.5, theta, phi);
        printf("Y(l=%+d,ml=%+d)*X_0(s=1/2,ms=+1/2) = (%+lf, %+lf)\n"
               "Y(l=%+d,ml=%+d)*X_1(s=1/2,ms=+1/2) = (%+lf, %+lf)\n\n"
               "Y(l=%+d,ml=%+d)*X_0(s=1/2,ms=-1/2) = (%+lf, %+lf)\n"
               "Y(l=%+d,ml=%+d)*X_1(s=1/2,ms=-1/2) = (%+lf, %+lf)\n\n",
               l, ml, creal(yp.u), cimag(yp.u), l, ml, 0., 0., l, ml, 0., 0., l,
               ml, creal(ym.d), cimag(ym.d));
    }

    printf("%s\n", "--- End");

    return 0;
}
