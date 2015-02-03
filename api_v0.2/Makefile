MODULES = argument_mod.o \
          kind_params_mod.o \
          global_parameters_mod.o \
          gocean_mod.o \
          grid_mod.o \
          kernel_mod.o \
          halo_mod.o \
          field_mod.o \
          omp_tiling_mod.o \
          region_mod.o \
          timing_mod.o \
          intel_timer_mod.o

all: ${API_LIB}

${API_LIB}: ${MODULES}
	${AR} ${ARFLAGS} ${API_LIB} ${MODULES}

clean:
	rm -f *.o *.mod *.a

# Inter-module dependencies, alphabetical order
argument_mod.o: global_parameters_mod.o

field_mod.o: kind_params_mod.o grid_mod.o \
             halo_mod.o omp_tiling_mod.o region_mod.o gocean_mod.o

grid_mod.o: kind_params_mod.o region_mod.o gocean_mod.o

halo_mod.o: region_mod.o

timing_mod.o: intel_timer_mod.o

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -c $<
