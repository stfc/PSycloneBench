INC_FLAGS=-I$(shell spack location -i raja -cuda)/include -I$(shell spack location -i umpire -cuda)/include -I$(shell spack location -i camp -cuda +openmp)/include
INC_FLAGS_CUDA=-I$(shell spack location -i raja +cuda)/include -I$(shell spack location -i umpire +cuda)/include -I$(shell spack location -i camp +cuda)/include

psy_layer: psy_layer.f90 cpp_psy_layer.o
	${F90} ${FFLAGS} -c $<

psy_layer_omp: psy_layer.f90 cpp_psy_layer_omp.o
	${F90} ${FFLAGS} -c $<

psy_layer_cuda: psy_layer.f90 cpp_psy_layer_cuda.o
	        ${F90} ${FFLAGS} -c $<

cpp_psy_layer.o: psy_layer.cpp
	${CXX} ${CXXFLAGS} ${INC_FLAGS} ${OMPFLAGS} -fPIE -c $< -o $@

cpp_psy_layer_omp.o: psy_layer.cpp
	${CXX} ${CXXFLAGS} ${INC_FLAGS} ${OMPFLAGS} -fPIE -DRAJA_ENABLE_OPENMP=On -c $< -o $@

cpp_psy_layer_cuda.o: psy_layer.cpp
	 $(shell spack location -i cuda@11)/bin/nvcc -forward-unknown-to-host-compiler -ccbin=g++ -DCAMP_HAVE_CUDA ${INC_FLAGS_CUDA} \
	--extended-lambda -O3 -DNDEBUG --generate-code=arch=compute_70,code=[compute_70,sm_70] -Xcompiler=-fPIE -std=c++14 -fopenmp -x cu -c $< -o $@

clean:
	rm -rf *.o *.mod
