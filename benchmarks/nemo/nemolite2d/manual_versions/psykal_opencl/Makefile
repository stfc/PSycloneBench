# Makefile for manual OpenCL version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler to use any flags from
# environment variables.

LDFLAGS += ${OMPFLAGS}
FFLAGS += ${OPENCL_INCLUDE}

.PHONY: all nemolite2d.exe device_binary
.DEFAULT_GOAL := nemolite2d.exe
include ../common.mk

COMMON_MODULES = time_step_mod.o nemolite2d.o

all: nemolite2d.exe

nemolite2d.exe: ${COMMON_LIB} inf_lib timer_lib fcl_lib ${COMMON_MODULES}
	${F90} -o $@ ${COMMON_MODULES} ${OPENCL_LIBS} ${COMMON_LIB} \
        ${FCL_LIB} ${INF_LIB} ${TIMER_LIB} -lOpenCL $(LDFLAGS)

device_binary: %.cl
	${OCL_COMPILER} $< -o $<.xo
	${OCL_LINKER} $<.xo -o $@

run:
	PSYCLONE_KERNELS_FILE=allkernels.cl ./nemolite2d.exe

clean:
	rm -f *.o *.mod *~ gnu_opt_report.txt *.optrpt *.dat

allclean: clean
	rm -f *.exe
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean
	${MAKE} -C ${COMMON_DIR} allclean
	${MAKE} -C ${FCL_DIR} allclean
