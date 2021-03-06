# Makefile for the NEMOLite2D routines that are common to all of the
# various versions of the benchmark.

MPI ?= no

SHARED_DIR = ../../../../shared
TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src

ifeq ($(MPI), yes)
        INF_LIB_TARGET=dm_fd_lib
        TIMER_TARGET=dm_lib
else
        INF_LIB_TARGET=fd_lib
        TIMER_TARGET=sm_lib
endif

MODULES = physical_params_mod.o \
          gocean2d_io_mod.o \
          model_mod.o \
          initialisation_mod.o

nemolite2d_common.a: inf_lib timer_lib ${MODULES}
	${AR} ${ARFLAGS} $@ ${MODULES}

.PHONY: inf_lib timer_lib

timer_lib:
	${MAKE} -C ${TIMER_DIR} ${TIMER_TARGET}

inf_lib:
	${MAKE} -C ${INF_DIR} ${INF_LIB_TARGET}

%.o: %.f90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<
%.o: %.F90
	$(F90) $(F90FLAGS) -I${INF_INC} -I${TIMER_INC} -c $<

clean:
	rm -f *.o *~ gnu_opt_report.txt *.mod

allclean: clean
	rm -f nemolite2d_common.a
