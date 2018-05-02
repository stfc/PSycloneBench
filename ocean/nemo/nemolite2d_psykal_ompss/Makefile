F90 = gfortran
F90FLAGS = -ffree-line-length-none

# Location of the code defining the GOcean API
API_VERSION=1.0
API_DIR = ../../api_v${API_VERSION}
API_LIB = ${API_DIR}/gocean_api.a

# Location of the dl_timer librar
TIMER_INC = /home/rupert/proj/dl_timer/src
TIMER_LIB = /home/rupert/proj/dl_timer/dl_timer_lib.a

RDIR=../nemolite2d_psykal

KERNELS = continuity_mod.o \
          momentum_mod.o \
	  boundary_conditions_mod.o \
          time_update_mod.o

SUPPORT = initialisation_mod.o \
	  model_mod.o \
	  physical_params_mod.o \
          gocean2d_io_mod.o \
          infrastructure_mod.o

ALGORITHM = nemolite2d_gen.o

PSY = psy.o

FOBJ = ${KERNELS} ${SUPPORT} ${PSY} ${ALGORITHM}

nemolite2d_ompss: ${FOBJ}
	${F90} -o $@ ${FOBJ} -fopenmp ${TIMER_LIB} ${API_LIB} ${LDFLAGS}

psy.o: psy.f90
	${F90} ${F90FLAGS} -I${API_DIR} -I${TIMER_INC} -c $<
clean:
	rm -f *.o *.mod nemolite2d_ompss

nemolite2d_gen.o: ${KERNELS} ${SUPPORT} ${PSY}

# Interdependencies between modules, alphabetical order

boundary_conditions_mod.o: physical_params_mod.o ${API_LIB} model_mod.o
continuity_mod.o: model_mod.o ${API_LIB}
gocean2d_io_mod.o: ${API_LIB}
infrastructure_mod.o: ${API_LIB}
model_mod.o: ${API_LIB} gocean2d_io_mod.o
momentum_mod.o: model_mod.o physical_params_mod.o ${API_DIR}/kind_params_mod.o
time_step_mod.o: ${API_LIB} momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o
time_update_mod.o: model_mod.o ${API_LIB}

# Generic rules

%.o: ${RDIR}/%.f90
	${F90} ${F90FLAGS} -I${API_DIR} -I${TIMER_INC} -c $<

