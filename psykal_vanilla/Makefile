# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file.
include ../common/Makefile.include

# For timing module
TIMING_DIR = ../common/

F90FLAGS += -I${TIMING_DIR}

# The target that this Makefile supports
EXECS = shallow_base

# The modules that are common to both targets
MODULES = kind_params_mod.o physical_params_mod.o field_mod.o model_mod.o \
          shallow_io_mod.o \
          initial_conditions_mod.o \
          compute_cu_mod.o compute_cv_mod.o \
          compute_z_mod.o compute_h_mod.o \
          compute_unew_mod.o compute_vnew_mod.o compute_pnew_mod.o \
          apply_bcs_mod.o time_step_mod.o \
          time_smooth_mod.o \
          ${TIMING_DIR}timing_mod.o ${TIMING_DIR}intel_timer_mod.o \
          mesh_mod.o topology_mod.o

# My local copies of GungHo modules. Avoids attempting to build
# the GungHo code which stretches some compilers...
GUNGHO_MODULES = global_parameters_mod.o argument_mod.o kernel_mod.o

COMMON_MODULES = $(GUNGHO_MODULES) $(MODULES)

all: $(EXECS)

# Normal targets
shallow_base: timing
	${MAKE} MODULE_LIST="shallow_base.o ${COMMON_MODULES}" shallow_base.exe

shallow_base.o: $(COMMON_MODULES)

timing:
	${MAKE} -C ${TIMING_DIR} F90="${F90}" F90FLAGS="${F90FLAGS}"

# Interdependencies between modules, alphabetical order

apply_bcs_mod.o: kind_params_mod.o topology_mod.o

compute_cu_mod.o: topology_mod.o

compute_fluxes_mod.o: compute_cu_mod.o \
                      compute_cv_mod.o \
                      compute_z_mod.o \
                      compute_h_mod.o \
                      topology_mod.o

field_mod.o: kind_params_mod.o topology_mod.o

initial_conditions_mod.o: physical_params_mod.o topology_mod.o model_mod.o

time_smooth_mod.o: kind_params_mod.o topology_mod.o $(GUNGHO_MODULES)

time_step_mod.o: kind_params_mod.o \
                 compute_cu_mod.o compute_cv_mod.o \
                 compute_z_mod.o compute_h_mod.o \
                 compute_unew_mod.o compute_vnew_mod.o compute_pnew_mod.o

mesh_mod.o: kind_params_mod.o

model_mod.o: field_mod.o shallow_io_mod.o mesh_mod.o \
             time_smooth_mod.o topology_mod.o

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -c $<

clean: 
	${MAKE} -C ${TIMING_DIR} clean
	rm -f *.o *.mod *.MOD *~ psy.f90 shallow_gen.f90

allclean: clean
	rm -f *.exe fparser.log
	rm -rf shallow_prog_lib

