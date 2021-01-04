
# Location of the dl_timer and infrastucture code
SHARED_DIR = ../../../../../shared

# dl_timer
TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/libdl_timer_omp.a
TIMER_LIB_DM = ${TIMER_DIR}/libdl_timer_dm.a
# dl_esm_inf
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/lib_fd.a
# Common (parts of NemoLite2D that are common to all implementations)
COMMON_DIR = ../../common
COMMON_LIB = ${COMMON_DIR}/nemolite2d_common.a
# FortCL (provides OpenCL functionality in Fortran)
FCL_DIR = ${SHARED_DIR}/FortCL/
FCL_INC = ${FCL_DIR}/src
FCL_LIB = ${FCL_INC}/libFortCL.a


# The kernels used by this application and their location
KERNEL_DIR = ../kernels/fortran
KERNELS = boundary_conditions_mod.o \
          continuity_mod.o \
          momentum_mod.o \
          time_update_mod.o \
          infrastructure_mod.o

# Generic rules

%.o: %.f90
	$(F90) $(F90FLAGS) -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -I${FCL_INC} -c $<

%.o: %.F90
	$(F90) $(F90FLAGS) -I${COMMON_DIR} -I${INF_INC} -I${TIMER_INC} -I${FCL_INC} -c $<

%.o: %.mod

# If we need a .f90 file that doesn't exist then it must be a kernel.
# Create a link to the required file.
%.f90:
	ln -sf ../../kernels/fortran/$@ .

# Common rules for compiling external libraries
.PHONY: timer_lib timer_lib timer_lib_parallel inf_lib_parallel fcl_lib ${COMMON_LIB}

timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

inf_lib:
	${MAKE} -C ${INF_DIR}

timer_lib_parallel:
	${MAKE} -C ${TIMER_DIR} MPIF90="${F90}" dm_lib

inf_lib_parallel:
	${MAKE} -C ${INF_DIR} F90="${F90}" dm_fd_lib

fcl_lib:
	${MAKE} -C ${FCL_DIR}

${COMMON_LIB}:
	${MAKE} -C ${COMMON_DIR} MPI=${MPI}
