#!/bin/bash

gcc -L../lib/ states.c -lalkcalc -Wl,-rpath,../lib/
./a.out
rm a.out
