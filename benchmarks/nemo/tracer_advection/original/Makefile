# Makefile for the original form of the tra_adv benchmark. Contains
# targets for building the original (serial) version and three
# different PSyclone-processed versions:
#
#  1. Serial;
#  2. OpenACC kernels for GPU;
#  3. OpenACC loops for GPU;
#  4. OpenMP k-loops parallelism for CPU;
#  5. OpenMP all-loops parallelism for CPU;
#  6. OpenMP offload for GPU.
#
# Compiler and associated flags must be set through the following
# environment variables:
#
#   export F90=gfortran
#   export F90FLAGS=-O3
#   export OMPFLAGS=-fopenmp
#
# See the scripts in PSycloneBench/compiler_setup for examples.

VERSIONS := tra_adv_serial tra_adv_acc_kernels tra_adv_acc_loops tra_adv_omp \
		   	tra_adv_omp_gpu tra_adv_omp_gpu
.PHONY: all ${VERSIONS}

# The default is to build the original, serial CPU version
.DEFAULT_GOAL := tra_adv_serial
all: ${VERSIONS}

# Location of the dl_timer code
DL_TIMER_DIR = ../../../../shared/dl_timer
DL_TIMER_NAME = libdl_timer_omp.a

# Shorthand for invoking PSyclone.
PSYCLONE = psyclone -api nemo -l output
PSYCLONE_PROFILE ?=
#PSYCLONE_PROFILE ?= -profile invokes


tra_adv_serial: dl_timer ./tra_adv.F90
	mkdir -p $@
	cp tra_adv.F90 $@/.
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} -I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ${OMPFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

# The OpenACC-accelerated version has the existing dl_timer calls but also
# has profiling instrumentation added that uses NVTX.  We profile 'invokes'
# rather than 'kernels' since the whole outer loop gets put into a single
# ACC KERNELS region.
tra_adv_acc_kernels: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} --profile invokes -s ../scripts/acc_kernels_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} -I${PSYCLONE_NVIDIA_LIB_DIR} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME} ${PSYCLONE_NVIDIA_LIB_DIR}/libnvtx_prof.a" -C $@

tra_adv_acc_loops: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} --profile invokes -s ../scripts/acc_loops_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} -I${PSYCLONE_NVIDIA_LIB_DIR} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME} ${PSYCLONE_NVIDIA_LIB_DIR}/libnvtx_prof.a" -C $@

tra_adv_omp_levels: ../scripts/omp_trans.py dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/omp_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${OMPFLAGS} -I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_omp_cpu: ../scripts/omp_cpu_trans.py dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/omp_cpu_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${OMPFLAGS} -I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_omp_gpu: ../scripts/omp_gpu_trans.py dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/omp_gpu_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
            FORT_FLAGS="${F90FLAGS} ${OMPTARGETFLAGS} ${UMEMFLAGS} -I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ${OMPTARGETFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

dl_timer:
	make -C ${DL_TIMER_DIR} sm_lib

clean:
ifneq (,$(wildcard ./tra_adv_serial))
	${MAKE} -C tra_adv_serial $@
endif
ifneq (,$(wildcard ./tra_adv_acc_kernels))
	${MAKE} -C tra_adv_acc_kernels $@
endif
ifneq (,$(wildcard ./tra_adv_acc_loops))
	${MAKE} -C tra_adv_acc_loops $@
endif
ifneq (,$(wildcard ./tra_adv_omp_levels))
	${MAKE} -C tra_adv_omp_levels $@
endif
ifneq (,$(wildcard ./tra_adv_omp_cpu))
	${MAKE} -C tra_adv_omp_cpu $@
endif
ifneq (,$(wildcard ./tra_adv_omp_gpu))
	${MAKE} -C tra_adv_omp_gpu $@
endif

allclean:
	make -C ${DL_TIMER_DIR} allclean
	rm -rf tra_adv_serial
	rm -rf tra_adv_acc_kernels
	rm -rf tra_adv_acc_loops
	rm -rf tra_adv_omp_levels
	rm -rf tra_adv_omp_gpu
	rm -rf tra_adv_omp_cpu
