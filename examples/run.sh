#!/bin/bash
#gcc -L../lib/ -L../GAUSSQ/lib/ -L../BSPLINES/lib/ radial_matrix_element.c -lalkcalc -lbsplines -lgaussq -Wl,-rpath,../lib/ -Wl,-rpath,../GAUSSQ/lib/ -Wl,-rpath,../BSPLINES/lib/
gcc -L../lib/ -L../GAUSSQ/lib/ -L../BSPLINES/lib/ lifetimes.c -lalkcalc -lbsplines -lgaussq -Wl,-rpath,../lib/ -Wl,-rpath,../GAUSSQ/lib/ -Wl,-rpath,../BSPLINES/lib/
./a.out
rm a.out
