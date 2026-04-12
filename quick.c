/* -------------------------------------------------------------------------- *
 * C program to quickly get some matrix element and energy difference         *
 *                                                                            *
 * Compile: gcc -L./lib/ quick.c -lalkcalc -lm -Wl,-rpath,./lib/              *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "./interface/alkcalc.h"

#define CONVERT(X) (int32_t)floor(2.*(X)+.5)

void quantities(char *, int, int, double, int, int, int, double);

int main(int argc, char **argv)
{
    (void)argc;

    char *species;
    int i, nb, lb, p, nk, lk;
    double jb, jk;


    /* Species */
    species = argv[1];

    /* Bra */
    (void)sscanf(argv[2], "%d", &nb);
    (void)sscanf(argv[3], "%d", &lb);
    (void)sscanf(argv[4], "%lf", &jb);

    /* Power */
    (void)sscanf(argv[5], "%d", &p);

    /* Ket */
    (void)sscanf(argv[6], "%d", &nk);
    (void)sscanf(argv[7], "%d", &lk);
    (void)sscanf(argv[8], "%lf", &jk);

    quantities(species, nb, lb, jb, p, nk, lk, jk);

    return 0;
}

/* Computes <nb,lb,s,jb|r^p|nk,lk,s,jk> and E(nb,lb,jb)-E(nk,lk,jk)           */
void quantities(char *species, int nb, int lb, double jb, int p, int nk, int lk,
                double jk) {

    char sign;
    double dE, rp;

    /* Energy difference */
    dE = alkcalc_Enlsj(species, nb, lb, jb)-alkcalc_Enlsj(species, nk, lk, jk);
    sign = (dE < 0.) ? '-': '+';
    dE = 6.579683920499956e6*fabs(dE);

    /* Matrix element */
    rp = fabs(alkcalc_rp(species, nb, lb, jb, (double)p, nk, lk, jk));

    /* Information */
    printf("\n========================================================\n\n"
           "                         RESULT\n\n"
           " Species %s\n"
           "\n"
           " <nb,lb,sb,jb| = <%d,%d,1/2,%d/2|\n"
           " |nk,lk,sk,jk> = |%d,%d,1/2,%d/2>\n"
           "\n"
           " Energy difference\n"
           " E(nb,lb,jb)-E(nk,lk,jk) = %c ħ × 2π × %1.2f GHz\n"
           "\n"
           " Matrix element (a₀ is Bohr's radius)\n"
           " |<nb,lb,sb,jb|rᵖ|nk,lk,sk,jk>| = %1.6f × a₀ᵖ\n"
           "\n"
           "========================================================\n\n",
           species, nb, lb, 2*(int)jb+1, nk, lk, 2*(int)jk+1, sign, dE, rp);
}
