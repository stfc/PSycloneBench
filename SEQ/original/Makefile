# rpc258
# 26 Jan 2010
# 20 Feb 2012
F90 = gfortran

F90FLAGS = -Wall -Wsurprising -fcheck=all -g -O0
#F90FLAGS = -O3
#F90FLAGS = -O3 -fopenmp
F90FLAGS += -I/usr/include
# For timing module in new version of code
F90FLAGS += -I../

# add -mpi=mpilog to generate mpe logging
LDFLAGS = ../timing.o -L/usr/lib64 -lnetcdf -lnetcdff

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
