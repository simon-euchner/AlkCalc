# ---------------------------------------------------------------------------- #
# Makefile for the tests of AlkCalc                                            #
#                                                                              #
# Author of this file: Simon Euchner                                           #
# ---------------------------------------------------------------------------- #


### Variables

# Timestamp
TMSTMP := $(shell date +"%d%m%Y%H%M%S%9N")

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


### Tests of the library functions (species 1H)
test: ${TMP}/tst${TMSTMP}
	@echo -e "\nRUNNING ALL TESTS FOR SPECIES 1H\n"
	- @${TMP}/tst${TMSTMP}
	@rm -f ${TMP}/tst${TMSTMP}
${TMP}/tst${TMSTMP}: ${OBJ}/${F0}.o
	@${LD} -o ${TMP}/tst${TMSTMP} -L${LIB}/ ${OBJ}/${F0}.o -lm -lalkcalc \
	-Wl,-rpath,$(abspath ${LIB}/)


### Compile

# tests.c
${OBJ}/${F0}.o: ${F0}.c
	${CC} ${FLAGS} -o ${OBJ}/${F0}.o -c ${F0}.c


### Cleanup

clean:
	- rm -f ${OBJ}/${F0}.o
	- rm -f ${TMP}/tst*

.PHONY: clean
