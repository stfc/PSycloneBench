# Makefile for manual MPI version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler to use any flags from
# environment variables. e.g. to use gfortran:
# 
# export F90=mpif90
# export F90FLAGS="-O3"
# export AR=ar

# Location of the dl_timer and infrastucture code
SHARED_DIR = ../../../../../shared

TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/dl_esm_inf_fd.a

COMMON_DIR = ../../common
COMMON_LIB = ${COMMON_DIR}/nemolite2d_common.a

KERNELS = boundary_conditions_mod.o \
          continuity_mod.o \
          momentum_mod.o \
          time_update_mod.o

${KERNELS}: ${COMMON_LIB}

# The modules that are common to both targets
MODULES = ${KERNELS} time_step_mod.o

# API lib is an archive that must come at the end of the list of objects
# passed to the linker
COMMON_MODULES = $(MODULES) ${INF_LIB}

.PHONY: nemolite2d nemolite2d-serial timer_lib inf_lib timer_lib_serial inf_lib_serial

all: nemolite2d

timer_lib:
	${MAKE} -C ${TIMER_DIR} F90=${F90} MPIF90=${F90} dm_lib

inf_lib:
	${MAKE} -C ${INF_DIR} F90=${F90} dm_fd_lib

timer_lib_serial:
	${MAKE} -C ${TIMER_DIR} F90=${F90} MPIF90=${F90} sm_lib

inf_lib_serial:
	${MAKE} -C ${INF_DIR} F90=${F90} fd_lib

${COMMON_LIB}:
	${MAKE} -C ${COMMON_DIR}

# Normal targets
nemolite2d: timer_lib inf_lib
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

nemolite2d-serial: timer_lib_serial inf_lib_serial
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" nemolite2d.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order

time_step_mod.o: ${COMMON_LIB} inf_lib

# Generic rules

%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) ${COMMON_LIB}  $(TIMER_LIB) $(LDFLAGS) ${OMPFLAGS}

%.o: %.f90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<
%.o: %.F90
	$(F90) $(F90FLAGS) ${OMPFLAGS} -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -c $<

# If we need a .f90 file that doesn't exist then it must be a kernel.
# Create a link to the required file...
%.f90:
	ln -sf ../../kernels/fortran/$@ .

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	make -C ${TIMER_DIR} allclean
	make -C ${INF_DIR} distclean
	make -C ${COMMON_DIR} allclean
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config
