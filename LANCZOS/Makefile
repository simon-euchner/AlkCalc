# ---------------------------------------------------------------------------- #
# Makefile for the C language interface to ARPACK's [3] Lanczos algorithm      #
#                                                                              #
# Author of this file: Simon Euchner                                           #
# ---------------------------------------------------------------------------- #


### Variables

# Compiler and Linker (GNU FORTRAN compiler and linker)
GF = gfortran
LD = gfortran
FLAGS = -Wall -pedantic -Wextra -fcheck=all -ffixed-form -std=legacy -Ofast

# Paths
SRC = ./SRC
OBJ = ./OBJ
LIB = ./lib

# Library
TARGET = ${LIB}/liblanczos.so

# File names
F0  = ICB
F1  = DSAUPD
F2  = DSAUP2
F3  = DSEUPD
F4  = DSAPPS
F5  = DSESRT
F6  = DSEIGT
F7  = DSGETS
F8  = DSCONV
F9  = DGETV0
FA  = DSAITR
FB  = DSTQRB
FC  = DSORTR
FD  = SECOND
FE  = IVOUT
FF  = DVOUT
FG  = DSTATS

# Object code
OBJ_FILENAMES = ${F0}.O ${F1}.O ${F2}.O ${F3}.O ${F4}.O ${F5}.O ${F6}.O ${F7}.O\
                ${F8}.O ${F9}.O ${FA}.O ${FB}.O ${FC}.O ${FD}.O ${FE}.O ${FF}.O\
				${FG}.O
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

# DSAUPD.F
${OBJ}/${F1}.O: ${SRC}/${F1}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F1}.O -c ${SRC}/${F1}.F

# DSAUP2.F
${OBJ}/${F2}.O: ${SRC}/${F2}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F2}.O -c ${SRC}/${F2}.F

# DSEUPD.F
${OBJ}/${F3}.O: ${SRC}/${F3}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F3}.O -c ${SRC}/${F3}.F

# DSAPPS.F
${OBJ}/${F4}.O: ${SRC}/${F4}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F4}.O -c ${SRC}/${F4}.F

# DSESRT.F
${OBJ}/${F5}.O: ${SRC}/${F5}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F5}.O -c ${SRC}/${F5}.F

# DSEIGT.F
${OBJ}/${F6}.O: ${SRC}/${F6}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F6}.O -c ${SRC}/${F6}.F

# DSGETS.F
${OBJ}/${F7}.O: ${SRC}/${F7}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F7}.O -c ${SRC}/${F7}.F

# DSCONV.F
${OBJ}/${F8}.O: ${SRC}/${F8}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F8}.O -c ${SRC}/${F8}.F

# DGETV0.F
${OBJ}/${F9}.O: ${SRC}/${F9}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F9}.O -c ${SRC}/${F9}.F

# DSAITR.F
${OBJ}/${FA}.O: ${SRC}/${FA}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FA}.O -c ${SRC}/${FA}.F

# DSTQRB.F
${OBJ}/${FB}.O: ${SRC}/${FB}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FB}.O -c ${SRC}/${FB}.F

# DSORTR.F
${OBJ}/${FC}.O: ${SRC}/${FC}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FC}.O -c ${SRC}/${FC}.F

# SECOND.F
${OBJ}/${FD}.O: ${SRC}/${FD}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FD}.O -c ${SRC}/${FD}.F

# DSTATS.F
${OBJ}/${FE}.O: ${SRC}/${FE}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FE}.O -c ${SRC}/${FE}.F

# DVOUT.F
${OBJ}/${FF}.O: ${SRC}/${FF}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FF}.O -c ${SRC}/${FF}.F

# IVOUT.F
${OBJ}/${FG}.O: ${SRC}/${FG}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FG}.O -c ${SRC}/${FG}.F


### Cleanup

clean:
	rm -f ${OBJ}/*.O
	rm -f ${LIB}/*.so

.PHONY: clean
