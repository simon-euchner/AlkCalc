# ---------------------------------------------------------------------------- #
# Makefile for eigensolver (LAPACK)                                            #
#                                                                              #
# Author of this file: Simon Euchner                                           #
# ---------------------------------------------------------------------------- #


### Variables

# Compiler and Linker (GNU FORTRAN compiler and linker)
GF = gfortran
LD = gfortran
FLAGS = -Wall -pedantic -Wextra -fcheck=all -ffixed-form -std=legacy -O2
FLAGSF90 = -Wall -pedantic -Wextra -fcheck=all -std=legacy -O2

# Paths
SRC = ./SRC
OBJ = ./OBJ
LIB = ./lib

# Library
TARGET = ${LIB}/libeiglapack.so

# File names
F00 = ICB
F01 = DASUM
F02 = DAXPY
F03 = DCOPY
F04 = DDOT
F05 = DGEMV
F06 = DGER
F07 = DISNAN
F08 = DLACPY
F09 = DLAE2
F10 = DLAEBZ
F11 = DLAEV2
F12 = DLAGTF
F13 = DLAGTS
F14 = DLAISNAN
F15 = DLAMCH
F16 = DLANST
F17 = DLAPY2
F18 = DLAR2V
F19 = DLARGV
F20 = DLARNV
F21 = DLARTG
F22 = DLARTV
F23 = DLARUV
F24 = DLASCL
F25 = DLASET
F26 = DLASR
F27 = DLASRT
F28 = DLASSQ
F29 = dnrm2
F30 = DPBSTF
F31 = DROT
F32 = DSBGST
F33 = DSBGVX
F34 = DSBTRD
F35 = DSCAL
F36 = DSTEBZ
F37 = DSTEIN
F38 = DSTEQR
F39 = DSTERF
F40 = DSWAP
F41 = DSYR
F42 = IDAMAX
F43 = IEEECK
F44 = ILAENV
F45 = IPARMQ
F46 = LSAME
F47 = XERBLA

# Object code
OBJ_FILENAMES = ${F00}.O ${F01}.O ${F02}.O ${F03}.O ${F04}.O ${F05}.O ${F06}.O\
				${F07}.O ${F08}.O ${F09}.O ${F10}.O ${F11}.O ${F12}.O ${F13}.O\
				${F14}.O ${F15}.O ${F16}.O ${F17}.O ${F18}.O ${F19}.O ${F20}.O\
				${F21}.O ${F22}.O ${F23}.O ${F24}.O ${F25}.O ${F26}.O ${F27}.O\
				${F28}.O ${F29}.o ${F30}.O ${F31}.O ${F32}.O ${F33}.O ${F34}.O\
				${F35}.O ${F36}.O ${F37}.O ${F38}.O ${F39}.O ${F40}.O ${F41}.O\
				${F42}.O ${F43}.O ${F44}.O ${F45}.O ${F46}.O ${F47}.O
OBJ_FILES = ${foreach file, ${OBJ_FILENAMES}, ${OBJ}/${file}}


### Print info
all: ${TARGET}
	@echo -e "\nBUILDING PROCESS COMPLETED SUCCESSFULLY\n"


### Link
${TARGET}: ${OBJ_FILES}
	${LD} -shared -o ${TARGET} ${wildcard ${OBJ}/*.O} ${OBJ}/dnrm2.o ${LIBS}


### Compile

# ICB.F
${OBJ}/${F00}.O: ${SRC}/${F00}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F00}.O -c ${SRC}/${F00}.F

# DASUM.F
${OBJ}/${F01}.O: ${SRC}/${F01}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F01}.O -c ${SRC}/${F01}.F

# DAXPY.F
${OBJ}/${F02}.O: ${SRC}/${F02}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F02}.O -c ${SRC}/${F02}.F

# DCOPY.F
${OBJ}/${F03}.O: ${SRC}/${F03}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F03}.O -c ${SRC}/${F03}.F

# DDOT.F
${OBJ}/${F04}.O: ${SRC}/${F04}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F04}.O -c ${SRC}/${F04}.F

# DGEMV.F
${OBJ}/${F05}.O: ${SRC}/${F05}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F05}.O -c ${SRC}/${F05}.F

# DGER.F
${OBJ}/${F06}.O: ${SRC}/${F06}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F06}.O -c ${SRC}/${F06}.F

# DISNAN.F
${OBJ}/${F07}.O: ${SRC}/${F07}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F07}.O -c ${SRC}/${F07}.F

# DLACPY.F
${OBJ}/${F08}.O: ${SRC}/${F08}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F08}.O -c ${SRC}/${F08}.F

# DLAE2.F
${OBJ}/${F09}.O: ${SRC}/${F09}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F09}.O -c ${SRC}/${F09}.F

# DLAEBZ.F
${OBJ}/${F10}.O: ${SRC}/${F10}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F10}.O -c ${SRC}/${F10}.F

# DLAEV2.F
${OBJ}/${F11}.O: ${SRC}/${F11}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F11}.O -c ${SRC}/${F11}.F

# DLAGTF.F
${OBJ}/${F12}.O: ${SRC}/${F12}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F12}.O -c ${SRC}/${F12}.F

# DLAGTS.F
${OBJ}/${F13}.O: ${SRC}/${F13}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F13}.O -c ${SRC}/${F13}.F

# DLAISNAN.F
${OBJ}/${F14}.O: ${SRC}/${F14}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F14}.O -c ${SRC}/${F14}.F

# DLAMCH.F
${OBJ}/${F15}.O: ${SRC}/${F15}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F15}.O -c ${SRC}/${F15}.F

# DLANST.F
${OBJ}/${F16}.O: ${SRC}/${F16}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F16}.O -c ${SRC}/${F16}.F

# DLAPY2.F
${OBJ}/${F17}.O: ${SRC}/${F17}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F17}.O -c ${SRC}/${F17}.F

# DLAR2V.F
${OBJ}/${F18}.O: ${SRC}/${F18}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F18}.O -c ${SRC}/${F18}.F

# DLARGV.F
${OBJ}/${F19}.O: ${SRC}/${F19}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F19}.O -c ${SRC}/${F19}.F

# DLARNV.F
${OBJ}/${F20}.O: ${SRC}/${F20}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F20}.O -c ${SRC}/${F20}.F

# DLARTG.F
${OBJ}/${F21}.O: ${SRC}/${F21}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F21}.O -c ${SRC}/${F21}.F

# DLARTV.F
${OBJ}/${F22}.O: ${SRC}/${F22}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F22}.O -c ${SRC}/${F22}.F

# DLARUV.F
${OBJ}/${F23}.O: ${SRC}/${F23}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F23}.O -c ${SRC}/${F23}.F

# DLASCL.F
${OBJ}/${F24}.O: ${SRC}/${F24}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F24}.O -c ${SRC}/${F24}.F

# DLASET.F
${OBJ}/${F25}.O: ${SRC}/${F25}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F25}.O -c ${SRC}/${F25}.F

# DLASR.F
${OBJ}/${F26}.O: ${SRC}/${F26}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F26}.O -c ${SRC}/${F26}.F

# DLASRT.F
${OBJ}/${F27}.O: ${SRC}/${F27}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F27}.O -c ${SRC}/${F27}.F

# DLASSQ.F
${OBJ}/${F28}.O: ${SRC}/${F28}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F28}.O -c ${SRC}/${F28}.F

# dnrm2.f90
${OBJ}/${F29}.o: ${SRC}/${F29}.f90
	${GF} ${FLAGSF90} -fPIC -o ${OBJ}/${F29}.o -c ${SRC}/${F29}.f90

# DPBSTF.F
${OBJ}/${F30}.O: ${SRC}/${F30}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F30}.O -c ${SRC}/${F30}.F

# DROT.F
${OBJ}/${F31}.O: ${SRC}/${F31}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F31}.O -c ${SRC}/${F31}.F

# DSBGST.F
${OBJ}/${F32}.O: ${SRC}/${F32}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F32}.O -c ${SRC}/${F32}.F

# DSBGVX.F
${OBJ}/${F33}.O: ${SRC}/${F33}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F33}.O -c ${SRC}/${F33}.F

# DSBTRD.F
${OBJ}/${F34}.O: ${SRC}/${F34}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F34}.O -c ${SRC}/${F34}.F

# DSCAL.F
${OBJ}/${F35}.O: ${SRC}/${F35}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F35}.O -c ${SRC}/${F35}.F

# DSTEBZ.F
${OBJ}/${F36}.O: ${SRC}/${F36}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F36}.O -c ${SRC}/${F36}.F

# DSTEIN.F
${OBJ}/${F37}.O: ${SRC}/${F37}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F37}.O -c ${SRC}/${F37}.F

# DSTEQR.F
${OBJ}/${F38}.O: ${SRC}/${F38}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F38}.O -c ${SRC}/${F38}.F

# DSTERF.F
${OBJ}/${F39}.O: ${SRC}/${F39}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F39}.O -c ${SRC}/${F39}.F

# DSWAP.F
${OBJ}/${F40}.O: ${SRC}/${F40}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F40}.O -c ${SRC}/${F40}.F

# DSYR.F
${OBJ}/${F41}.O: ${SRC}/${F41}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F41}.O -c ${SRC}/${F41}.F

# IDAMAX.F
${OBJ}/${F42}.O: ${SRC}/${F42}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F42}.O -c ${SRC}/${F42}.F

# IEEECK.F
${OBJ}/${F43}.O: ${SRC}/${F43}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F43}.O -c ${SRC}/${F43}.F

# ILAENV.F
${OBJ}/${F44}.O: ${SRC}/${F44}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F44}.O -c ${SRC}/${F44}.F

# IPARMQ.F
${OBJ}/${F45}.O: ${SRC}/${F45}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F45}.O -c ${SRC}/${F45}.F

# LSAME.F
${OBJ}/${F46}.O: ${SRC}/${F46}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F46}.O -c ${SRC}/${F46}.F

# XERBLA.F
${OBJ}/${F47}.O: ${SRC}/${F47}.F
	${GF} ${FLAGS} -fPIC -o ${OBJ}/${F47}.O -c ${SRC}/${F47}.F

### Cleanup

clean:
	rm -f ${OBJ}/*.O
	rm -f ${OBJ}/*.o
	rm -f ${LIB}/*.so

.PHONY: clean
