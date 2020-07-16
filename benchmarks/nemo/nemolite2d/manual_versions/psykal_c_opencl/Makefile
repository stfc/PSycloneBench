# Makefile for manual C OpenCL version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler to use any flags from
# environment variables. e.g.:
# export F90=gfortran
# export CXX=g++
# export CFLAGS="-O3"
# export F90FLAGS="-O3"

# Additionally the environment needs the 'gcc' compiler regardless of the compiler
# selected, this is to perform some preprocessor-only commands when preparing the
# single-file OpenCL kernels.

.PHONY: nemolite2d xilinxoffline alteraoffline
.DEFAULT_GOAL := nemolite2d

include ../common.mk

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}
COMMON_MODULES = time_step_mod.o ${COMMON_LIB} ${KERNELS} ${INF_LIB}

# Host targets
nemolite2d: ${COMMON_LIB} timer_lib inf_lib opencl_utils.o timing.o time_step.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} opencl_utils.o timing.o time_step.o" \
		nemolite2d.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order
time_step_mod.o: ${COMMON_LIB} inf_lib momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o

time_step.o: opencl_utils.o

# Linking the executable. It need the -lOpenCL library
%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) \
		$(LDFLAGS) ${OMPFLAGS} -lOpenCL

# C rules (since kernels are included directly any modification in the kernel
# files must trigger a recompilation)
%.o: %.c ../../kernels/c_family/*.c
	$(CC) $(CFLAGS) $(OMPFLAGS) -c $< -o $@

# Target for Offline compilation of Xilinx FPGA device kernels
PLATFORM ?= xilinx_u200_xdma_201830_2
EXECUTION_TARGET ?= hw # Options: sw_emu | hw_emu | hw
xilinxoffline:
	echo "// OpenCL Kernels" > allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/continuity_kern.c >> allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/momentum_u_kern.c >> allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/momentum_v_kern.c >> allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/boundary_conditions_kern.c >> allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/time_update_kern.c >> allkernels.cl
	v++ --compile --jobs 4 -O3 --platform ${PLATFORM} --target ${EXECUTION_TARGET} \
		allkernels.cl -o allkernels_xilinx.xo
	v++ --link --jobs 4 -O3 --platform ${PLATFORM} --target ${EXECUTION_TARGET} \
		allkernels_xilinx.xo -o allkernels_xilinx.xclbin

# Target for Offline compilation of Intel Altera FPGA device kernels
AOC_FLAGS = -v -board=p385a_sch_ax115 -D__OPENCL_VERSION__=100
#AOC_FLAGS += -march=emulator -emulator-channel-depth-model=strict
alteraoffline:
	$(error "Not supported yet")
	echo "// OpenCL Kernels" > allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/continuity_kern.c >> allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/momentum_u_kern.c >> allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/momentum_v_kern.c >> allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/boundary_conditions_kern.c >> allkernels.cl
	gcc -E -P -D__OPENCL_VERSION__ ../../kernels/c_family/time_update_kern.c >> allkernels.cl
	aoc ${AOC_FLAGS} allkernels.cl -o allkernels_altera.aocx

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ *.dat
	rm -f gnu_opt_report.txt *.optrpt *.log

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe
	rm -rf *_wpl_dir/ allkernels* _x/ .Xil/ .run/
	${MAKE} -C ${COMMON_DIR} allclean
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean

summary:
	@echo -n "${CC}\t\tJIT\t"
	@./nemolite2d.exe | awk 'BEGIN {ORS="\t"} \
		{if ($$2 == "checksum") {print $$4} } \
		{if ($$1 == "Time-stepping") {print $$5} } \
		END {ORS="\n"; print ""}'
