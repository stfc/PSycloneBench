# This Makefile picks up the compiler to use, the location of the
# OpenCL libraries and any flags from environment variables. e.g. to
# use AMD's dispatch library:
# 
# export F90=gfortran
# export CC=gcc
# export OPENCL_LIBS="-L/opt/AMDAPPSDK-3.0/lib/x86_64 -lOpenCL"
# export OPENCL_INCLUDE=-I/opt/AMDAPPSDK-3.0/include
# export CFLAGS="-g -O0 -Wall"
# export LDFLAGS=-lm

LDFLAGS += -lOpenCL
LDFLAGS += ${OMPFLAGS}

.PHONY: all nemolite2d.exe
.DEFAULT_GOAL := nemolite2d.exe
include ../common.mk



# FortCL (provides OpenCL functionality in Fortran)
FCL_DIR = ${SHARED_DIR}/dl_esm_inf/external/FortCL/
FCL_INC = ${FCL_DIR}/src
FCL_LIB = ${FCL_INC}/libFortCL.a

COMMON_MODULES = time_step_mod.o nemolite2d.o

all: nemolite2d.exe

nemolite2d.exe: ${COMMON_LIB} inf_lib timer_lib fcl_lib ${COMMON_MODULES}
	${F90} -o $@ ${COMMON_MODULES} ${OPENCL_LIBS} ${COMMON_LIB} \
        ${FCL_LIB} ${INF_LIB} ${TIMER_LIB} $(LDFLAGS)

fcl_lib:
	${MAKE} -C ${FCL_DIR}

# We only need to compile clfortran in order to get the module interfaces.
# We *do not* link it into the final executable otherwise it collides
# with the real OpenCL run-time routines and nothing works!
opencl_utils_mod.o: clfortran.o


clean:
	rm -f *.o *~ gnu_opt_report.txt

allclean: clean
	rm -f *.exe *.mod
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean
	${MAKE} -C ${COMMON_DIR} allclean
	${MAKE} -C ${FCL_DIR} allclean
