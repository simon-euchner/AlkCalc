#!/bin/sh
gcc -L./lib/ quick.c -lalkcalc -lm -Wl,-rpath,./lib/ -o quick
