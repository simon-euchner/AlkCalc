#!/bin/sh
gcc -L./lib.d/ quick.c -lalkcalc -lm -Wl,-rpath,./lib.d/ -o quick
