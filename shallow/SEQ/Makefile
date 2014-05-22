# Makefile.include is a symbolic link to the desired
# Makefile.include.<compiler> file.
include Makefile.include

# Location of the psyclone source files (for doing the code
# generation)
PSYCLONE_DIR = ../../../PSyclone/src

# add -mpi=mpilog to generate mpe logging
LDFLAGS += -lnetcdff -lnetcdf

EXECS = shallow_base shallow_gen

MODULES = kind_params.o physical_params.o field.o model.o shallow_IO.o \
          initial_conditions.o compute_cu_mod.o compute_cv_mod.o \
          compute_z.o compute_h.o manual_invoke_compute_new_fields_mod.o \
          compute_unew.o compute_vnew.o compute_pnew.o \
          manual_invoke_apply_bcs_mod.o \
          apply_bcs_cf_mod.o apply_bcs_ct_mod.o \
          apply_bcs_cu_mod.o apply_bcs_cv_mod.o \
          timing.o time_smooth.o intel_timer_mod.o mesh.o topology_mod.o

# My local copies of GungHo modules. Avoids attempting to build
# the GungHo code which stretches some compilers...
GUNGHO_MODULES = global_parameters.o argument.o kernel.o

all: $(EXECS)

# Targets involving the code-generation framework
shallow_gen.f90: shallow_gocean.f90
	python ${PSYCLONE_DIR}/generator.py -api gocean shallow_gocean.f90 -oalg shallow_gen.f90 -opsy psy.f90

# psy.f90 is generated at the same time as shallow_gen.f90
psy.f90: shallow_gen.f90

# The generated code depends on the generated Psy middle-layer
shallow_gen.o: $(MODULES) psy.o

# Normal targets
shallow_base.o: $(MODULES)

timing.o: intel_timer_mod.o

model.o: field.o timing.o shallow_IO.o mesh.o time_smooth.o topology_mod.o

time_smooth.o: $(GUNGHO_MODULES)

initial_conditions.o: model.o

manual_invoke_compute_new_fields_mod.o: compute_unew.o compute_vnew.o compute_pnew.o
manual_invoke_apply_bcs_mod.o: apply_bcs_cf_mod.o apply_bcs_ct_mod.o \
          apply_bcs_cu_mod.o apply_bcs_cv_mod.o

%: %.o 
	$(F90) -o $@ $< $(GUNGHO_MODULES) $(MODULES) $(LDFLAGS)

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -c $<

clean: 
	rm -f *.o *.mod *.MOD *~ psy.f90 shallow_gen.f90

allclean: clean
	rm -f $(EXECS)

docs:
	doxygen shallow.doxy.config
