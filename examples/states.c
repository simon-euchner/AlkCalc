/* -------------------------------------------------------------------------- *
 * Example: Eigenenergies                                                     *
 *                                                                            *
 * Compile: gcc -L../lib/ eigenenergies.c -lalkcalc -Wl,-rpath,../lib/        *
 *                                                                            *
 * Author of this file: Simon Euchner                                         *
 * -------------------------------------------------------------------------- */

#include "../interface/alkcalc.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    char result, *species;
    int32_t n, l;
    double j;
    alkcalc_state *state;

    result = 'f';
    species = "88SR+";
    n = 10;
    l = 0;
    j = .5;

    state = alkcalc_fnlsj(result, species, n, l, j);

    alkcalc_state_free(state);

    return 0;
}
