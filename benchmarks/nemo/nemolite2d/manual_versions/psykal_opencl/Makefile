# Makefile for manual OpenCL version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler and flags to use from environment
# variables. Beside the standard compiler and optimization flags, if JIT
# OpenCL is used the OPENCL_INCLUDE variable is required, if an offline OpenCL
# compiler is used the OCL_COMPILER and OCL_LINKER are also required.

LDFLAGS += ${OMPFLAGS}
FFLAGS += ${OPENCL_INCLUDE}

.PHONY: all device_binary device_binary_nohup device_binary_tasks \
	    device_binary_tasks_nohup
.DEFAULT_GOAL := nemolite2d.exe
include ../common.mk

COMMON_MODULES = time_step_mod.o nemolite2d.o

all: nemolite2d.exe

nemolite2d.exe: ${COMMON_LIB} inf_lib timer_lib fcl_lib ${COMMON_MODULES}
	${F90} -o $@ ${COMMON_MODULES} ${OPENCL_LIBS} ${COMMON_LIB} \
        ${FCL_LIB} ${INF_LIB} ${TIMER_LIB} -lOpenCL $(LDFLAGS)

# If OpenCL JIT is not used, this PHONY target can be invoked to
# offline-compile the kernels into a device binary object. This is
# the common workflow for FPGA devices.
# The generated device file can be used with:
# > FORTCL_KERNELS_FILE=device_binary ./nemolite2d.exe
device_binary: allkernels.cl
	${OCL_COMPILER} $< -o $<.xo
	${OCL_LINKER} $<.xo -o $@.xclbin

device_binary_nohup:
	nohup make device_binary > nohup.out &

device_binary_tasks: allkernels_tasks.cl
	${OCL_COMPILER} $< -o $<.xo
	${OCL_LINKER} $<.xo -o $@.xclbin

device_binary_tasks_nohup:
	nohup make device_binary_tasks > nohup.out &

run_jit:
	FORTCL_KERNELS_FILE=allkernels.cl ./nemolite2d.exe

clean:
	rm -f *.o *.mod *~ gnu_opt_report.txt *.optrpt *.dat *.log
	rm -rf _x device_binary.* allkernels.cl.* nohup.out .ipcache

allclean: clean
	rm -f *.exe
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean
	${MAKE} -C ${COMMON_DIR} allclean
	${MAKE} -C ${FCL_DIR} allclean

summary:
	@echo -n "${CC}\t\tundef\t"
	FORTCL_KERNELS_FILE=allkernels.cl ./nemolite2d.exe | awk \
		'BEGIN {ORS="\t"} \
		{if ($$2 == "checksum") {print $$4} } \
		{if ($$1 == "Time-stepping") {print $$5} } \
		END {ORS="\n"; print ""}'
