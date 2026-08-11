#!/bin/bash
#gcc -L../lib/ radial_matrix_element.c -lalkcalc -Wl,-rpath,../lib/
gcc -L../lib/ eigenenergies.c -lalkcalc -Wl,-rpath,../lib/
./a.out
rm a.out
