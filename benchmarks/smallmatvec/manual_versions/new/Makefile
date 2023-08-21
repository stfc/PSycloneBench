.DEFAULT_GOAL := fortran
.PHONY: clean fortran openmp_cpu openmp_gpu cpp kokkos_cpu kokkos_gpu

bin:
	mkdir bin

fortran: fortran/psy_layer.f90 kdriver.f90 bin fortran/Makefile
	$(MAKE) -C fortran
	${F90} ${F90FLAGS} -Ifortran kdriver.f90 fortran/psy_layer.o -o bin/$@

openmp_cpu: fortran/psy_layer.f90 kdriver.f90 bin fortran/Makefile
	$(MAKE) openmp_cpu -C fortran
	${F90} ${F90FLAGS} ${OMPFLAGS} -Ifortran kdriver.f90 fortran/openmp_cpu_psy_layer.o -o bin/$@
	
openmp_gpu: fortran/psy_layer.f90 kdriver.f90 bin fortran/Makefile
	$(MAKE) openmp_gpu -C fortran
	${F90} ${F90FLAGS} ${OMPTARGETFLAGS} -Ifortran kdriver.f90 fortran/openmp_gpu_psy_layer.o -o bin/$@
	
cpp: cpp/psy_layer.f90 cpp/psy_layer.cpp kdriver.f90 bin
	$(MAKE) -C cpp
	${F90} ${F90FLAGS} -lstdc++ -Icpp kdriver.f90 cpp/psy_layer.o cpp/cpp_psy_layer.o -o bin/$@

kokkos_cpu: kokkos/psy_layer.f90 kokkos/psy_layer.cpp kdriver.f90 bin kokkos/Makefile
	$(MAKE) kokkos_cpu -C kokkos
	${F90} ${F90FLAGS} -Ikokkos kdriver.f90 kokkos/psy_layer.o kokkos/kokkos_cpu_psy_layer.o -lstdc++ -o bin/$@

kokkos_gpu: kokkos/psy_layer.f90 kokkos/psy_layer.cpp kdriver.f90 bin kokkos/Makefile
	$(MAKE) kokkos_gpu -C kokkos
	${F90} ${F90FLAGS} -lstdc++ -Ikokkos kdriver.f90 kokkos/psy_layer.o kokkos/kokkos_gpu_psy_layer.o -o bin/$@

	
clean:
	$(MAKE) clean -C fortran
	$(MAKE) clean -C cpp
	$(MAKE) clean -C kokkos
	rm -rf bin *.o *.mod
