# ---------------------------------------------------------------------------- #
# Makefile for the part of SuperLU (see Ref. [7]) that implements sparse       #
# LU-decomposition for the type double 										   #
#                                                                              #
# Author of this file: Simon Euchner                                           #
# ---------------------------------------------------------------------------- #


### Variables

# Compiler and Linker (GNU C compiler and linker)
CC = gcc
LD = gcc
FLAGS = -Wall -pedantic -Wextra -std=c99 -Ofast

# Paths
SRC = ./src.d
OBJ = ./obj.d
LIB = ./lib.d

# Library
TARGET = ${LIB}/liblufac.so

# File names
F00 = dgssv
F01 = dmemory
F02 = memory
F03 = superlu_timer
F04 = input_error
F05 = sp_preorder
F06 = sp_ienv
F07 = dgstrf
F08 = dgstrs
F09 = get_perm_c
F10 = util
F11 = dutil
F12 = dpivotL
F13 = dpruneL
F14 = heap_relax_snode
F15 = sp_coletree
F16 = relax_snode
F17 = dsnode_dfs
F18 = dsnode_bmod
F19 = dpanel_bmod
F20 = dpanel_dfs
F21 = dcolumn_dfs
F22 = dcolumn_bmod
F23 = dcopy_to_ucol
F24 = colamd
F25 = dsp_blas2
F26 = dsp_blas3
F27 = mmd

# Object code
OBJ_FILENAMES = ${F00}.o ${F01}.o ${F02}.o ${F03}.o ${F04}.o ${F05}.o ${F06}.o \
			    ${F07}.o ${F08}.o ${F09}.o ${F10}.o ${F11}.o ${F12}.o ${F13}.o \
			    ${F14}.o ${F15}.o ${F16}.o ${F17}.o ${F18}.o ${F19}.o ${F20}.o \
			    ${F21}.o ${F22}.o ${F23}.o ${F24}.o ${F25}.o ${F26}.o ${F27}.o

OBJ_FILES = ${foreach file, ${OBJ_FILENAMES}, ${OBJ}/${file}}


### Print info
all: ${TARGET}
	@echo -e "\nBUILDING PROCESS COMPLETED SUCCESSFULLY\n"


### Link
${TARGET}: ${OBJ_FILES}
	${LD} -shared -o ${TARGET} ${wildcard ${OBJ}/*.o} ${LIBS}


### Compile

# dgssv.c
${OBJ}/${F00}.o: ${SRC}/${F00}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F00}.o -c ${SRC}/${F00}.c

# dmemory.c
${OBJ}/${F01}.o: ${SRC}/${F01}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F01}.o -c ${SRC}/${F01}.c

# memory.c
${OBJ}/${F02}.o: ${SRC}/${F02}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F02}.o -c ${SRC}/${F02}.c

# superlu_timer.c
${OBJ}/${F03}.o: ${SRC}/${F03}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F03}.o -c ${SRC}/${F03}.c

# input_error.c
${OBJ}/${F04}.o: ${SRC}/${F04}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F04}.o -c ${SRC}/${F04}.c

# sp_preorder.c
${OBJ}/${F05}.o: ${SRC}/${F05}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F05}.o -c ${SRC}/${F05}.c

# sp_ienv.c
${OBJ}/${F06}.o: ${SRC}/${F06}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F06}.o -c ${SRC}/${F06}.c

# dgstrf.c
${OBJ}/${F07}.o: ${SRC}/${F07}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F07}.o -c ${SRC}/${F07}.c

# dgstrs.c
${OBJ}/${F08}.o: ${SRC}/${F08}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F08}.o -c ${SRC}/${F08}.c

# get_perm_c.c
${OBJ}/${F09}.o: ${SRC}/${F09}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F09}.o -c ${SRC}/${F09}.c

# util.c
${OBJ}/${F10}.o: ${SRC}/${F10}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F10}.o -c ${SRC}/${F10}.c

# dutil.c
${OBJ}/${F11}.o: ${SRC}/${F11}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F11}.o -c ${SRC}/${F11}.c

# dpivotL.c
${OBJ}/${F12}.o: ${SRC}/${F12}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F12}.o -c ${SRC}/${F12}.c

# dpruneL.c
${OBJ}/${F13}.o: ${SRC}/${F13}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F13}.o -c ${SRC}/${F13}.c

# heap_relax_snod.c
${OBJ}/${F14}.o: ${SRC}/${F14}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F14}.o -c ${SRC}/${F14}.c

# sp_coletree.c
${OBJ}/${F15}.o: ${SRC}/${F15}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F15}.o -c ${SRC}/${F15}.c

# relax_snode.c
${OBJ}/${F16}.o: ${SRC}/${F16}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F16}.o -c ${SRC}/${F16}.c

# dsnode_dfs.c
${OBJ}/${F17}.o: ${SRC}/${F17}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F17}.o -c ${SRC}/${F17}.c

# dsnode_bmod.c
${OBJ}/${F18}.o: ${SRC}/${F18}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F18}.o -c ${SRC}/${F18}.c

# dpanel_bmod.c
${OBJ}/${F19}.o: ${SRC}/${F19}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F19}.o -c ${SRC}/${F19}.c

# dpanel_dfs.c
${OBJ}/${F20}.o: ${SRC}/${F20}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F20}.o -c ${SRC}/${F20}.c

# dcolumn_dfs.c
${OBJ}/${F21}.o: ${SRC}/${F21}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F21}.o -c ${SRC}/${F21}.c

# dcolumn_bmod.c
${OBJ}/${F22}.o: ${SRC}/${F22}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F22}.o -c ${SRC}/${F22}.c

# dcopy_to_ucol.c
${OBJ}/${F23}.o: ${SRC}/${F23}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F23}.o -c ${SRC}/${F23}.c

# colamd.c
${OBJ}/${F24}.o: ${SRC}/${F24}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F24}.o -c ${SRC}/${F24}.c

# dsp_blas2.c
${OBJ}/${F25}.o: ${SRC}/${F25}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F25}.o -c ${SRC}/${F25}.c

# dsp_blas3.c
${OBJ}/${F26}.o: ${SRC}/${F26}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F26}.o -c ${SRC}/${F26}.c

# mmd.c
${OBJ}/${F27}.o: ${SRC}/${F27}.c
	${CC} ${FLAGS} -fPIC -o ${OBJ}/${F27}.o -c ${SRC}/${F27}.c


### Cleanup

clean:
	rm -f ${OBJ}/*.o
	rm -f ${LIB}/*.so

.PHONY: clean
