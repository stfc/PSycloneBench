# Makefile for the original form of the tra_adv benchmark. Contains
# targets for building the original (serial) version and three
# different PSyclone-processed versions:
#
#  1. OpenACC;
#  2. OpenMP for CPU;
#  3. OpenMP offload for GPU.
#
# Compiler and associated flags must be set through the following
# environment variables:
#
#   export F90=gfortran
#   export F90FLAGS=-O3
#   export OMPFLAGS=-fopenmp
#
# See the scripts in PSycloneBench/compiler_setup for examples.

.PHONY: all tra_adv_serial tra_adv_acc tra_adv_omp tra_adv_omp_gpu

# The default is to build the original, serial CPU version
all: tra_adv_serial

# Location of the dl_timer code
DL_TIMER_DIR = ../../../../shared/dl_timer
DL_TIMER_NAME = libdl_timer_omp.a

# Shorthand for invoking PSyclone.
PSYCLONE = psyclone -api nemo -l output

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
tra_adv_acc: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} --profile invokes -s ../scripts/kernels_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} -I${PSYCLONE_NVIDIA_LIB_DIR} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME} ${PSYCLONE_NVIDIA_LIB_DIR}/libnvtx_prof.a" -C $@

tra_adv_omp: ../scripts/omp_trans.py dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/omp_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${OMPFLAGS} -I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_omp_gpu: ../scripts/omp_gpu_trans.py dl_timer ./tra_adv.F90
ifndef UMEMFLAGS
	$(error The OMP offload target requires OpenMP unified memory but the UMEMFLAGS environment variable is not set)
endif
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
ifneq (,$(wildcard ./tra_adv_acc))
	${MAKE} -C tra_adv_acc $@
endif
ifneq (,$(wildcard ./tra_adv_omp))
	${MAKE} -C tra_adv_omp $@
endif
ifneq (,$(wildcard ./tra_adv_omp_gpu))
	${MAKE} -C tra_adv_omp_gpu $@
endif

allclean:
	make -C ${DL_TIMER_DIR} allclean
	rm -rf tra_adv_serial
	rm -rf tra_adv_acc
	rm -rf tra_adv_omp
	rm -rf tra_adv_omp_gpu
