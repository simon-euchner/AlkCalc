# ---------------------------------------------------------------------------- #
# Makefile for testing AlkCalc                                                 #
#                                                                              #
# Author of this file: Simon Euchner                                           #
# ---------------------------------------------------------------------------- #


### Variables

# Compiler and Linker (GNU C compiler and linker)
CC = gcc
LD = gcc
FLAGS = -Wall -pedantic -Wextra -O2 -std=c99

# Paths
OBJ = ../obj
LIB = ../lib
TMP = ../tmp

# File names
F0 = tests


### Fallback
all:
	@echo -e "\nVALID ARGUMENTS: test\n"


### Testing of AlkCalc (species 1H)
test: ${TMP}/${F0}
	@echo -e "\nRUNNING TESTS\n"
	- @${TMP}/${F0}
	@rm -f ${TMP}/${F0}
${TMP}/${F0}${TMSTMP}: ${OBJ}/${F0}.o
	@${LD} -o ${TMP}/${F0} -L${LIB}/ ${OBJ}/${F0}.o -lm -lalkcalc \
	-Wl,-rpath,$(abspath ${LIB}/)


### Compile

# tests.c
${OBJ}/${F0}.o: ${F0}.c
	${CC} ${FLAGS} -o ${OBJ}/${F0}.o -c ${F0}.c


### Cleanup

clean:
	- rm -f ${OBJ}/${F0}.o
	- rm -f ${TMP}/${F0}

.PHONY: clean
