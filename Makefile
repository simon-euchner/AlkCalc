# ---------------------------------------------------------------------------- #
# Makefile for AlkCalc                                                         #
#                                                                              #
# Author of this file: Simon Euchner                                           #
# ---------------------------------------------------------------------------- #


### Variables

# Timestamp
TMSTMP := $(shell date +"%d%m%Y%H%M%S%9N")

# Compiler and Linker (GNU C compiler and linker)
CC = gcc
LD = gcc
FLAGS = -Wall -pedantic -Wextra -Ofast -std=c99

# Paths
SRC = ./src
OBJ = ./obj
LIB = ./lib
INF = ./interface
TMP = ./tmp
LUF = ./LUFac/lib
LCZ = ./LANCZOS/lib

# File names
F0 = potential
F1 = settings
F2 = eigensolver
F3 = alkcalc


### Fallback
all:
	@echo -e "\nVALID ARGUMENTS: 'lib', 'solve'\n"


### Library
lib: ${LIB}/libalkcalc.so
	@echo -e "\nBUILDING LIBRARY 'ALKCALC'\n"
${LIB}/libalkcalc.so: ${OBJ}/${F3}.o
	${LD} -shared -o ${LIB}/libalkcalc.so ${OBJ}/${F3}.o -lm


### Eigenenergies and radial eigenstates
solve: ${TMP}/e${TMSTMP}
	@echo -e "\nCOMPUTING EIGENENERGIES AND RADIAL EIGENSTATES\n"
	- @${TMP}/solve${TMSTMP}
	@rm -f ${TMP}/solve${TMSTMP}
${TMP}/e${TMSTMP}: ${OBJ}/${F0}.o ${OBJ}/${F1}.o ${OBJ}/${F2}.o
	@${LD} -o ${TMP}/solve${TMSTMP} -L${LUF}/ -L${LCZ}/ ${OBJ}/${F0}.o \
	${OBJ}/${F1}.o ${OBJ}/${F2}.o -lm -lblas -llapack -llufac -llanczos \
	-Wl,-rpath,{${LUF}/,${LCZ}/}


### Compile

# potential.c
${OBJ}/${F0}.o: ${SRC}/${F0}.c
	${CC} ${FLAGS} -o ${OBJ}/${F0}.o -c ${SRC}/${F0}.c

# settings.c
${OBJ}/${F1}.o: ${INF}/${F1}.c
	${CC} ${FLAGS} -o ${OBJ}/${F1}.o -c ${INF}/${F1}.c

# eigensolver.c
${OBJ}/${F2}.o: ${SRC}/${F2}.c
	${CC} ${FLAGS} -o ${OBJ}/${F2}.o -c ${SRC}/${F2}.c

# alkcalc.c
${OBJ}/${F3}.o: ${SRC}/${F3}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F3}.o -c ${SRC}/${F3}.c


### Cleanup

clean:
	- rm -f ${OBJ}/*.o
	- rm -f ${LIB}/*.so
	- rm -f ./tmp/solve*

.PHONY: clean
