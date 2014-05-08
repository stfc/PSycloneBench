include ../Makefile.include

F90FLAGS += -I${NETCDF_INC}
# For timing module in new version of code
F90FLAGS += -I../

# add -mpi=mpilog to generate mpe logging
LDFLAGS = ../timing.o ../intel_timer_mod.o -L${NETCDF_LIB} -lnetcdf -lnetcdff

#EXECS = join_files shallow_base shallow_parallel
EXECS = shallow_base

all: $(EXECS)

shallow_base.o: ../timing.o

shallow_parallel.o: parallel_mod.o

shallow_parallel: shallow_parallel.o
	$(F90) $(F90FLAGS) -o $@ $^ parallel_mod.o $(LDFLAGS)

%: %.o 
	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm -f *.o *.mod *.MOD $(EXECS)
