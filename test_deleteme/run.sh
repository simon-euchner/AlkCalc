#!/bin/bash
gcc --std=c99 -pedantic -Wall -Wextra -c test.c -o test.o
gcc potential.o test.o settings.o -o main -lm
./main
rm main
