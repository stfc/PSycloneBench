# rpc258
# 26 Jan 2010
# 20 Feb 2012
F90 = mpif90
F90FLAGS = -O3
#F90FLAGS = -O3 -fopenmp
F90FLAGS += -I/home/grahamr/software/netcdf/include \
	    -I/home/cnc/MPI/mpich2-1.4.1-install/include
#	    -I/home/grahamr/software/mpich2/install/include
# add -mpi=mpilog to generate mpe logging
LDFLAGS = -L/home/grahamr/software/netcdf/lib/ -lnetcdf \
	   -L/home/cnc/MPI/mpich2-1.4.1-install/lib
#	   -L/home/grahamr/software/mpich2/install/lib
#EXECS = join_files shallow_base shallow_parallel
EXECS = shallow_base

all: $(EXECS)

shallow_parallel.o: parallel_mod.o

shallow_parallel: shallow_parallel.o
	$(F90) $(F90FLAGS) -o $@ $^ parallel_mod.o $(LDFLAGS)

%: %.o 
	$(F90) $(F90FLAGS) -o $@ $^ $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm *.o *.mod *.MOD $(EXECS)
