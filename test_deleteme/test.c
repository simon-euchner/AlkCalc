#include "../inc/eigensolver.h"

int main(int argc, char **argv)
{
    (void)argc; (void)argv;

    int32_t ipar[4]; double rpar[10];

    potential_initpar(ipar, rpar);

    validate_settings(ipar[3]);

    double r;
    for (int i = 0; i < 100; i++) {
        r = 1e-6 + 1e-4 * i * i;
        printf("%+1.3E %+1.3E\n", r, V(r, ipar, rpar));
    }

    return 0;
}
