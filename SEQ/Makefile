# rpc258
# 26 Jan 2010
# 20 Feb 2012

UTILS_DIR = ../../../GungHo/src/utils

F90 = gfortran

F90FLAGS = -Wall -Wsurprising -fcheck=all -g -O0
#F90FLAGS = -O3
#F90FLAGS = -O3 -fopenmp
F90FLAGS += -I/usr/include -I${UTILS_DIR}

# add -mpi=mpilog to generate mpe logging
LDFLAGS = -L/usr/lib64 -lnetcdf -lnetcdff

#EXECS = join_files shallow_base shallow_parallel
EXECS = shallow_base

MODULES = kind_params.o physical_params.o field.o model.o shallow_IO.o \
          initial_conditions.o compute_z.o compute_h.o \
          compute_unew.o compute_vnew.o compute_pnew.o \
          timing.o time_smooth.o intel_timer_mod.o mesh.o          

GUNGHO_MODULES = ${UTILS_DIR}/global_parameters.o \
                 ${UTILS_DIR}/utils.o \
                 ${UTILS_DIR}/topology_module.o \
                 ${UTILS_DIR}/lfric.o \
                 ${UTILS_DIR}/debug.o \
                 ${UTILS_DIR}/debug_c.o \
                 argument.o kernel.o

all: $(EXECS)

shallow_base.o: $(MODULES)

timing.o: intel_timer_mod.o

model.o: field.o timing.o shallow_IO.o mesh.o time_smooth.o

time_smooth.o: $(GUNGHO_MODULES)

initial_conditions.o: model.o

%: %.o 
	$(F90) $(F90FLAGS) -o $@ $< $(GUNGHO_MODULES) $(MODULES) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm -f *.o *.mod *.MOD $(EXECS)

docs:
	doxygen shallow.doxy.config
