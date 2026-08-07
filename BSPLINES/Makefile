# ---------------------------------------------------------------------------- #
# Makefile for de Boor's B-splines                                             #
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
TARGET = ${LIB}/libbsplines.so

# File names
F0 = ICB
F1 = DINTRV
F2 = DBSPVD
F3 = DBSPVN
F4 = XERCNT
F5 = XERHLT
F6 = XERMSG
F7 = XERPRN
F8 = XERSVE
F9 = XGETUA
FA = J4SAVE
FB = I1MACH
FC = FDUMP

# Object code
OBJ_FILENAMES = ${F0}.O ${F1}.O ${F2}.O ${F3}.O ${F4}.O ${F5}.O ${F6}.O ${F7}.O\
                ${F8}.O ${F9}.O ${FA}.O ${FB}.O ${FC}.O
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

# DINTRV.F
${OBJ}/${F1}.O: ${SRC}/${F1}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F1}.O -c ${SRC}/${F1}.F

# DBSPVD.F
${OBJ}/${F2}.O: ${SRC}/${F2}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F2}.O -c ${SRC}/${F2}.F

# DBSPVN.F
${OBJ}/${F3}.O: ${SRC}/${F3}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F3}.O -c ${SRC}/${F3}.F

# XERCNT.F
${OBJ}/${F4}.O: ${SRC}/${F4}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F4}.O -c ${SRC}/${F4}.F

# XERHNT.F
${OBJ}/${F5}.O: ${SRC}/${F5}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F5}.O -c ${SRC}/${F5}.F

# XERMSG.F
${OBJ}/${F6}.O: ${SRC}/${F6}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F6}.O -c ${SRC}/${F6}.F

# XERPRN.F
${OBJ}/${F7}.O: ${SRC}/${F7}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F7}.O -c ${SRC}/${F7}.F

# XERSVE.F
${OBJ}/${F8}.O: ${SRC}/${F8}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F8}.O -c ${SRC}/${F8}.F

# XGETUA.F
${OBJ}/${F9}.O: ${SRC}/${F9}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F9}.O -c ${SRC}/${F9}.F

# J4SAVE.F
${OBJ}/${FA}.O: ${SRC}/${FA}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FA}.O -c ${SRC}/${FA}.F

# I1MACH.F
${OBJ}/${FB}.O: ${SRC}/${FB}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FB}.O -c ${SRC}/${FB}.F

# FDUMP.F
${OBJ}/${FC}.O: ${SRC}/${FC}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${FC}.O -c ${SRC}/${FC}.F

### Cleanup

clean:
	rm -f ${OBJ}/*.O
	rm -f ${LIB}/*.so

.PHONY: clean
