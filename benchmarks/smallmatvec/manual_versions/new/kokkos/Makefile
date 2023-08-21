
KOKKOS_PATH := ../../../../../shared/kokkos
# Import needed kokkos environment variables
include ${KOKKOS_PATH}/Makefile.kokkos
# Careful, KOKKOS_DEBUG can add 10x performance penalty in kernels.
KOKKOS_DEBUG := no

psy_layer.o: psy_layer.f90
	${F90} ${FFLAGS} -c $<

kokkos_cpu: KOKKOS_DEVICES="OpenMP"
kokkos_cpu: KOKKOS_OPTIONS="aggresive_vectorization"
kokkos_cpu: psy_layer.cpp psy_layer.o $(KOKKOS_CPP_DEPENDS)
	${CXX} $(KOKKOS_CXXFLAGS) ${CXXFLAGS} -DEXECUTION_SPACE=${KOKKOS_DEVICES} -c psy_layer.cpp -o $@_psy_layer.o

kokkos_gpu: psy_layer.cpp psy_layer.o $(KOKKOS_CPP_DEPENDS)
	KOKKOS_ARCH = Volta70
	# nvcc doesn't support all compiler flags, overwrite CFLAGS to just O3
	CFLAGS := -O3
	# Still use the selected compiler but using the Kokkos nvcc_wrapper
	NVCC_WRAPPER_DEFAULT_COMPILER = $(CXX)
	CXX := $(KOKKOS_PATH)/bin/nvcc_wrapper
	# The enable lambda option is necessary for the nvcc compiler to recognise
	# as CUDA kernels the lambda-inlined functions.
	KOKKOS_CUDA_OPTIONS = "enable_lambda"
	# If CUDA_LIB is not provided, infer path from the nvcc compiler location.
	CUDA_LIB ?= $(shell echo $(shell which nvcc) | sed 's/bin\/nvcc/lib64/g')
	${CXX} ${CXXFLAGS} -c psy_layer.cpp -o $@_psy_layer.o
	
clean: kokkos-clean
	rm -rf *.o *.mod
