# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file.
include Makefile.include

# Location of the psyclone source files (for doing the code
# generation)
PSYCLONE_DIR = ../../../PSyclone/src

# For linking with netcdf libraries
LDFLAGS += -lnetcdff -lnetcdf

# The two targets that this Makefile supports
#  - shallow_base is the version of the code with manual invokes
#  - shallow_gen uses PSyClone to generate the invokes
EXECS = shallow_base shallow_gen

# The modules that are common to both targets
MODULES = kind_params_mod.o physical_params_mod.o field_mod.o model_mod.o \
          shallow_io_mod.o \
          initial_conditions_mod.o \
          compute_fluxes_mod.o compute_cu_mod.o compute_cv_mod.o \
          compute_z_mod.o compute_h_mod.o \
          compute_new_fields_mod.o \
          compute_unew_mod.o compute_vnew_mod.o compute_pnew_mod.o \
          apply_bcs_mod.o \
          timing_mod.o time_smooth_mod.o intel_timer_mod.o \
          mesh_mod.o topology_mod.o

GENERATED_MODULES = psy.o shallow_gen.o

# My local copies of GungHo modules. Avoids attempting to build
# the GungHo code which stretches some compilers...
GUNGHO_MODULES = global_parameters_mod.o argument_mod.o kernel_mod.o

COMMON_MODULES = $(GUNGHO_MODULES) $(MODULES)

all: $(EXECS)

# Targets involving the code-generation framework
shallow_gen.f90: shallow_gocean.f90
	python ${PSYCLONE_DIR}/generator.py -api gocean shallow_gocean.f90 -oalg shallow_gen.f90 -opsy psy.f90

# psy.f90 is generated at the same time as shallow_gen.f90
psy.f90: shallow_gen.f90

# The generated code depends on the generated Psy middle-layer
shallow_gen:
	${MAKE} MODULE_LIST="${COMMON_MODULES} ${GENERATED_MODULES}" shallow_gen.exe

# Normal targets
shallow_base:
	${MAKE} MODULE_LIST="shallow_base.o ${COMMON_MODULES}" shallow_base.exe

shallow_base.o: $(COMMON_MODULES)

shallow_gen.o: $(COMMON_MODULES) ${GENERATED_MODULES}

# Interdependencies between modules, alphabetical order

apply_bcs_mod.o: kind_params_mod.o topology_mod.o

compute_cu_mod.o: topology_mod.o

field_mod.o: kind_params_mod.o topology_mod.o

initial_conditions_mod.o: physical_params_mod.o model_mod.o

time_smooth_mod.o: kind_params_mod.o topology_mod.o $(GUNGHO_MODULES)

timing_mod.o: intel_timer_mod.o

mesh_mod.o: kind_params_mod.o

model_mod.o: field_mod.o timing_mod.o shallow_io_mod.o mesh_mod.o \
             time_smooth_mod.o topology_mod.o

compute_fluxes_mod.o: compute_cu_mod.o \
                      compute_cv_mod.o \
                      compute_z_mod.o \
                      compute_h_mod.o \
                      topology_mod.o

compute_new_fields_mod.o: compute_unew_mod.o \
                          compute_vnew_mod.o \
                          compute_pnew_mod.o \
                          topology_mod.o

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $? $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm -f *.o *.mod *.MOD *~ psy.f90 shallow_gen.f90

allclean: clean
	rm -f *.exe fparser.log

docs:
	doxygen shallow.doxy.config
