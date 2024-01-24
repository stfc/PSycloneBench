

psy_layer: psy_layer.f90 cpp_psy_layer.o
	${F90} ${F90FLAGS} -c $<

cpp_psy_layer.o: psy_layer.cpp
	${CXX} ${CXXFLAGS} -fPIE -c $< -o $@

openmp_cpu_cpp: psy_layer.cpp psy_layer
	        ${CXX} ${CXXFLAGS} ${OMPFLAGS} -fPIE -c $< -o openmp_cpu_psy_layer.o

openmp_gpu_cpp: psy_layer.cpp psy_layer
	        ${CXX} ${CXXFLAGS} -DTARGET_GPU ${OMPTARGETFLAGS} ${UMEMFLAGS} -c $< -o openmp_gpu_psy_layer.o

clean:
	rm -rf *.o *.mod
