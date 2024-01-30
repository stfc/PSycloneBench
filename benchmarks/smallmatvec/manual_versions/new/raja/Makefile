INC_DIR =/apps/spack/spack-repo/opt/spack/linux-centos8-skylake_avx512/gcc-8.5.0/raja-2022.10.4-hjdblqvnrcjmczfjkclegnngdfg7f7md/include

psy_layer: psy_layer.f90 cpp_psy_layer.o
	${F90} ${FFLAGS} -c $<

psy_layer_omp: psy_layer.f90 cpp_psy_layer_omp.o
	${F90} ${FFLAGS} -c $<

psy_layer_cuda: psy_layer.f90 cpp_psy_layer_cuda.o
	        ${F90} ${FFLAGS} -c $<

cpp_psy_layer.o: psy_layer.cpp
	${CXX} ${CXXFLAGS} -I${INC_DIR} -fPIE -c $< -o $@

cpp_psy_layer_omp.o: psy_layer.cpp
	${CXX} ${CXXFLAGS} -I${INC_DIR} ${OMPFLAGS} -fPIE -DRAJA_ENABLE_OPENMP=On -c $< -o $@

cpp_psy_layer_cuda.o: psy_layer.cpp
	        ${CXX} ${CXXFLAGS} ${OMPFLAGS} -I${INC_DIR} -fPIE -DENABLE_CUDA=On -DENABLE_OPENMP=Off -DCMAKE_CUDA_ARCHITECTURES=70 -DCMAKE_CUDA_COMPILER=nvcc -DCUDA_TOOLKIT_ROOT_DIR=${CUDA_HOME} -DBLT_CXX_STD=c++14 -DCMAKE_BUILD_TYPE=Release -DRAJA_ENABLE_EXERCISES=On -DRAJA_ENABLE_OPENMP=Off -DCMAKE_CUDA_FLAGS=--extended-lambda -c $< -o $@

clean:
	rm -rf *.o *.mod
