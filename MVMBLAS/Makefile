# ---------------------------------------------------------------------------- #
# Makefile for banded matrix-vector multiplication                             #
#                                                                              #
# Author of this file: Simon Euchner                                           #
# ---------------------------------------------------------------------------- #


### Variables

# Compiler and Linker (GNU FORTRAN compiler and linker)
GF = gfortran
LD = gfortran
FLAGS = -Wall -pedantic -Wextra -fcheck=all -ffixed-form -std=legacy -O2

# Paths
SRC = ./SRC
OBJ = ./OBJ
LIB = ./lib

# Library
TARGET = ${LIB}/libmvmblas.so

# File names
F0 = ICB
F1 = DSBMV
F2 = LSAME
F3 = XERBLA

# Object code
OBJ_FILENAMES = ${F0}.O ${F1}.O ${F2}.O ${F3}.O
OBJ_FILES = ${foreach file, ${OBJ_FILENAMES}, ${OBJ}/${file}}


### Print info
all: ${TARGET}
	@echo -e "\nBUILDING PROCESS COMPLETED SUCCESSFULLY\n"


### Link
${TARGET}: ${OBJ_FILES}
	${LD} -shared -o ${TARGET} ${wildcard ${OBJ}/*.O} ${LIBS}


### Compile

# ICB.F
${OBJ}/${F0}.O: ${SRC}/${F0}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F0}.O -c ${SRC}/${F0}.F

# DSBMV.F
${OBJ}/${F1}.O: ${SRC}/${F1}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F1}.O -c ${SRC}/${F1}.F

# LSAME.F
${OBJ}/${F2}.O: ${SRC}/${F2}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F2}.O -c ${SRC}/${F2}.F

# XERBLA.F
${OBJ}/${F3}.O: ${SRC}/${F3}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F3}.O -c ${SRC}/${F3}.F


### Cleanup

clean:
	rm -f ${OBJ}/*.O
	rm -f ${LIB}/*.so

.PHONY: clean
