# Makefile for the original form of the tra_adv benchmark. Contains
# targets for building the original (serial) version and multiple
# different PSyclone-processed versions:
#
#  1. Serial;
#  2. OpenMP k-loops parallelism for CPU;
#  3. OpenMP loop parallelism for CPU;
#  4. OpenACC kernels for GPU with Unified Memory;
#  5. OpenACC kernels for GPU with Explicit Data Movement;
#  6. OpenACC explicit loop parallelism for GPU with Unified Memory;
#  7. OpenACC explicit loop parallelism for GPU with Explicit Data Movement;
#  8. OpenACC mixed (kernels and explicit loops) parallelism for GPU with Unified Memory;
#  9. OpenACC mixed (kernels and explicit loops) parallelism for GPU with Explicit Data Movement;
#  10. OpenMP loop parallelism for GPU with Unified Memory.
#
# Compiler and associated flags must be set through the following
# environment variables:
#
#   export F90=gfortran
#   export F90FLAGS=-O3
#   export OMPFLAGS=-fopenmp
#
# See the scripts in PSycloneBench/compiler_setup for examples.
#
# The optional ENABLE_NVIDIA_PROFILE=yes flag can be set to enable
# GPU offloading versions to generate instrumented code. e.g.
#
#   make tra_adv_acc_kernels ENABLE_NVIDIA_PROFILE=yes
#
# this will require the PSYCLONE_NVIDIA_LIB_DIR to be properly set up.

VERSIONS := tra_adv_serial tra_adv_omp_cpu tra_adv_omp_cpu_levels \
			tra_adv_acc_kernels_unified_memory \
			tra_adv_acc_kernels_explicit_data_movement	\
			tra_adv_acc_loops_unified_memory \
			tra_adv_acc_loops_explicit_data_movement \
			tra_adv_acc_mixed_unified_memory \
			tra_adv_acc_mixed_explicit_data_movement \
			tra_adv_omp_gpu_unified_memory

.PHONY: all clean allclean compile-all run-all ${VERSIONS}

# The default is to build the original, serial CPU version
.DEFAULT_GOAL := tra_adv_serial
all: ${VERSIONS}

# Location of the dl_timer code
DL_TIMER_DIR = ../../../../shared/dl_timer
DL_TIMER_NAME = libdl_timer_omp.a

# Shorthand for invoking PSyclone.
PSYCLONE = psyclone -api nemo -l output ${PSYCLONE_PROFILE}

# Add necessary flags for Nvidia nvtx instrumentation
ifeq ($(ENABLE_NVIDIA_PROFILE),yes)
PSYCLONE_PROFILE := -p invokes
F90FLAGS += -I${PSYCLONE_NVIDIA_LIB_DIR}
LDFLAGS += -L${NVHPC_ROOT}/cuda/lib64 ${PSYCLONE_NVIDIA_LIB_DIR}/libnvtx_prof.a -lnvToolsExt -Mcuda
else
PSYCLONE_PROFILE :=
endif

# Instructions to compile CPU versions

tra_adv_serial: dl_timer ./tra_adv.F90
	mkdir -p $@
	cp tra_adv.F90 $@/.
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} -I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ${OMPFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_omp_cpu_levels: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/omp_cpu_levels_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${OMPFLAGS} -I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_omp_cpu: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/omp_cpu_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${OMPFLAGS} -I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

# Instructions to compile GPU versions

tra_adv_acc_kernels_unified_memory: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/acc_kernels_unified_memory_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${ACCFLAGS} ${UMEMFLAGS} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ${ACCFLAGS} ${UMEMFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_acc_kernels_explicit_data_movement: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/acc_kernels_explicit_data_movement_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${ACCFLAGS} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ${ACCFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_acc_loops_unified_memory: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/acc_loops_unified_memory_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${ACCFLAGS} ${UMEMFLAGS} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ${ACCFLAGS} ${UMEMFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_acc_loops_explicit_data_movement: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/acc_loops_explicit_data_movement_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${ACCFLAGS} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ${ACCFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_acc_mixed_unified_memory: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/acc_mixed_unified_memory_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${ACCFLAGS} ${UMEMFLAGS} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ${ACCFLAGS} ${UMEMFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_acc_mixed_explicit_data_movement: dl_timer ./tra_adv.F90
	mkdir -p $@
	${PSYCLONE} -s ../scripts/acc_mixed_explicit_data_movement_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} FORT_FLAGS="${F90FLAGS} ${ACCFLAGS} -I../${DL_TIMER_DIR}/src" \
           LDFLAGS="${LDFLAGS} ${ACCFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

tra_adv_omp_gpu_unified_memory: dl_timer ./tra_adv.F90
ifndef UMEMFLAGS
	$(error The OMP offload target requires OpenMP unified memory but the UMEMFLAGS environment variable is not set)
endif
	mkdir -p $@
	${PSYCLONE} -s ../scripts/omp_gpu_trans.py -opsy $@/tra_adv.f90 ./tra_adv.F90
	cp Makefile_gen $@/Makefile
	${MAKE} PROF_LIB_INC="-I../${DL_TIMER_DIR}/src" \
            FORT_FLAGS="${F90FLAGS} ${OMPTARGETFLAGS} ${UMEMFLAGS} -I../${DL_TIMER_DIR}/src" \
            LDFLAGS="${LDFLAGS} ${OMPTARGETFLAGS} ${UMEMFLAGS} ../${DL_TIMER_DIR}/${DL_TIMER_NAME}" -C $@

dl_timer:
	make -C ${DL_TIMER_DIR} sm_lib

clean:
	rm -rf ${VERSIONS}

allclean: clean
	make -C ${DL_TIMER_DIR} allclean

compile-all:
	for version in ${VERSIONS} ; do \
        echo "Compiling $$version" ; \
		make $$version ; \
    done

run-all:
	for version in ${VERSIONS} ; do \
        echo "Running $$version" ; \
		JPI=128 JPJ=128 JPK=75 IT=100 ./$$version/tra_adv.exe ; \
    done
