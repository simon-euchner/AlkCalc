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
FLAGS = -Wall -pedantic -Wextra -O2 -std=c99

# Paths
SRC = ./src
OBJ = ./obj
LIB = ./lib
INF = ./interface
TMP = ./tmp
GAQ = ./GAUSSQ/lib
BSP = ./BSPLINES/lib
ELA = ./EIGLAPACK/lib

# File names
F0 = settings
F1 = potential
F2 = validate
F3 = eigensolver
F4 = alkcalc


### Fallback
all:
	@echo -e "\nVALID ARGUMENTS: lib, solve\n"


### Library
lib: ${LIB}/libalkcalc.so
	@echo -e "\nBUILDING LIBRARY ALKCALC\n"
${LIB}/libalkcalc.so: ${OBJ}/${F4}.o
	${LD} -shared -o ${LIB}/libalkcalc.so -L${GAQ}/ -L${BSP}/ ${OBJ}/${F4}.o \
	-lm -lcblas -lgaussq -lbsplines -Wl,-rpath,$(abspath ${GAQ}/) \
	-Wl,-rpath,$(abspath ${BSP})

### Eigenenergies and radial eigenstates
solve: ${TMP}/slv${TMSTMP}
	@echo -e "\nCOMPUTING EIGENENERGIES AND RADIAL EIGENSTATES\n"
	- @${TMP}/slv${TMSTMP}
	@rm -f ${TMP}/slv${TMSTMP}
${TMP}/slv${TMSTMP}: ${OBJ}/${F0}.o ${OBJ}/${F1}.o ${OBJ}/${F2}.o ${OBJ}/${F3}.o
	@${LD} -o ${TMP}/slv${TMSTMP} -L${GAQ}/ -L${BSP}/ -L${ELA}/ \
	${OBJ}/${F0}.o ${OBJ}/${F1}.o ${OBJ}/${F2}.o ${OBJ}/${F3}.o -lm -lgaussq \
	-lbsplines -leiglapack -Wl,-rpath,{${GAQ}/,${BSP}/,${ELA}/}


### Compile

# settings.c
${OBJ}/${F0}.o: ${INF}/${F0}.c
	${CC} ${FLAGS} -o ${OBJ}/${F0}.o -c ${INF}/${F0}.c

# potential.c
${OBJ}/${F1}.o: ${SRC}/${F1}.c
	${CC} ${FLAGS} -o ${OBJ}/${F1}.o -c ${SRC}/${F1}.c

# validate.c
${OBJ}/${F2}.o: ${SRC}/${F2}.c
	${CC} ${FLAGS} -o ${OBJ}/${F2}.o -c ${SRC}/${F2}.c

# eigensolver.c
${OBJ}/${F3}.o: ${SRC}/${F3}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F3}.o -c ${SRC}/${F3}.c

# alkcalc.c
${OBJ}/${F4}.o: ${SRC}/${F4}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F4}.o -c ${SRC}/${F4}.c


### Cleanup

clean:
	- rm -f ${OBJ}/*.o
	- rm -f ${LIB}/*.so
	- rm -f ./tmp/slv*

.PHONY: clean
