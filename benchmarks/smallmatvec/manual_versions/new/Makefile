.DEFAULT_GOAL := fortran
.PHONY: clean fortran openmp_cpu openmp_gpu cpp openmp_cpu_cpp openmp_gpu_cpp kokkos_cpu kokkos_gpu raja

bin:
	mkdir bin

# Fortran targets
fortran: fortran/psy_layer.F90 kdriver.f90 bin fortran/Makefile
	$(MAKE) -C fortran
	${F90} ${F90FLAGS} -Ifortran kdriver.f90 fortran/psy_layer.o -o bin/$@

openmp_cpu: fortran/psy_layer.F90 kdriver.f90 bin fortran/Makefile
	$(MAKE) openmp_cpu -C fortran
	${F90} ${F90FLAGS} ${OMPFLAGS} -Ifortran kdriver.f90 fortran/openmp_cpu_psy_layer.o -o bin/$@
	
openmp_gpu: fortran/psy_layer.F90 kdriver.f90 bin fortran/Makefile
	$(MAKE) openmp_gpu -C fortran
	${F90} ${F90FLAGS} ${OMPTARGETFLAGS} -Ifortran kdriver.f90 fortran/openmp_gpu_psy_layer.o -o bin/$@

# C++ targets
cpp: cpp/psy_layer.f90 cpp/psy_layer.cpp kdriver.f90 bin cpp/Makefile
	$(MAKE) -C cpp
	${F90} ${F90FLAGS} -Icpp kdriver.f90 cpp/psy_layer.o cpp/cpp_psy_layer.o -lstdc++ -o bin/$@


openmp_cpu_cpp: cpp/psy_layer.f90 cpp/psy_layer.cpp kdriver.f90 bin cpp/Makefile
	$(MAKE) openmp_cpu_cpp -C cpp
	${F90} ${F90FLAGS} ${OMPFLAGS} -Icpp kdriver.f90  cpp/psy_layer.o cpp/openmp_cpu_psy_layer.o -lstdc++ -o bin/$@

openmp_gpu_cpp: cpp/psy_layer.f90 cpp/psy_layer.cpp kdriver.f90 bin cpp/Makefile
	$(MAKE) openmp_gpu_cpp -C cpp
	${F90} ${F90FLAGS} ${OMPTARGETFLAGS} ${UMEMFLAGS} -Icpp kdriver.f90 cpp/psy_layer.o cpp/openmp_gpu_psy_layer.o -lstdc++ -o bin/$@

# RAJA targets
RAJA_LIB_DIR=$(shell spack location -i raja -cuda)/lib
RAJA_LIB_DIR_CUDA=$(shell spack location -i raja +cuda)/lib
raja: raja/psy_layer.f90 raja/psy_layer.cpp kdriver.f90 bin
	$(MAKE) -C raja psy_layer
	${F90} ${F90FLAGS} ${OMPFLAGS} -Iraja kdriver.f90 raja/psy_layer.o raja/cpp_psy_layer.o -L${RAJA_LIB_DIR} -lRAJA -lstdc++ -o bin/$@

raja_omp: raja/psy_layer.f90 raja/psy_layer.cpp kdriver.f90 bin
	$(MAKE) -C raja psy_layer_omp
	${F90} ${F90FLAGS} ${OMPFLAGS} -Iraja kdriver.f90 raja/psy_layer.o raja/cpp_psy_layer_omp.o -L${RAJA_LIB_DIR} -lRAJA -lstdc++ -o bin/$@

raja_cuda: raja/psy_layer.f90 raja/psy_layer.cpp kdriver.f90 bin
	$(MAKE) -C raja psy_layer_cuda
	${F90} ${F90FLAGS} ${OMPFLAGS} ${UMEMFLAGS} -Iraja kdriver.f90 raja/psy_layer.o raja/cpp_psy_layer_cuda.o \
	    -L${RAJA_LIB_DIR_CUDA} -lRAJA -L$(shell spack location -i cuda@11)/lib64 -lcudart -lstdc++ -o bin/$@

# Kokkos targets
kokkos_cpu: kokkos/psy_layer.f90 kokkos/psy_layer.cpp kdriver.f90 bin kokkos/Makefile
	$(MAKE) kokkos_cpu -C kokkos
	${F90} ${F90FLAGS} -Ikokkos -Lkokkos kdriver.f90 kokkos/psy_layer.o kokkos/kokkos_cpu_psy_layer.o \
	    -lstdc++ ${OMPFLAGS} -lkokkos -ldl -latomic -lpthread -o bin/$@

# If CUDA_LIB is not provided, infer path from the nvcc compiler location.
CUDA_LIB ?= $(shell echo $(shell which nvcc) | sed 's/bin\/nvcc/lib64/g')
kokkos_gpu: kokkos/psy_layer.f90 kokkos/psy_layer.cpp kdriver.f90 bin kokkos/Makefile
	$(MAKE) kokkos_gpu -C kokkos
	${F90} ${F90FLAGS} -lstdc++ -Ikokkos kdriver.f90 kokkos/psy_layer.o kokkos/kokkos_gpu_psy_layer.o -o bin/$@

	
clean:
	$(MAKE) clean -C fortran
	$(MAKE) clean -C cpp
	$(MAKE) clean -C kokkos
	$(MAKE) clean -C raja
	rm -rf bin *.o *.mod
