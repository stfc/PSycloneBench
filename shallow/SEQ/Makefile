# rpc258
# 26 Jan 2010
# 20 Feb 2012
F90 = gfortran

F90FLAGS = -Wall -Wsurprising -fcheck=all -g -O0
#F90FLAGS = -O3
#F90FLAGS = -O3 -fopenmp
F90FLAGS += -I/usr/include

# add -mpi=mpilog to generate mpe logging
LDFLAGS = -L/usr/lib64 -lnetcdf -lnetcdff

#EXECS = join_files shallow_base shallow_parallel
EXECS = shallow_base

MODULES = model.o shallow_IO.o manual_invoke_initialise.o \
          timing.o time_smooth.o intel_timer_mod.o mesh.o

all: $(EXECS)

shallow_base.o: $(MODULES)

timing.o: intel_timer_mod.o

model.o: timing.o shallow_IO.o mesh.o time_smooth.o manual_invoke_initialise.o

%: %.o 
	$(F90) $(F90FLAGS) -o $@ $< $(MODULES) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm -f *.o *.mod *.MOD $(EXECS)

docs:
	doxygen shallow.doxy.config
