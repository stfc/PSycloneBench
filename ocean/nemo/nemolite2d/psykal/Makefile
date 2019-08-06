# Makefile for nemolite2d using PSyclone to generate the
# PSy layer.
#
# This Makefile expects the following env. vars to be set:
#
# F90       - How to invoke the Fortran compiler
# F90FLAGS  - Flags to pass to the Fortran compiler
# OMPFLAGS  - Flags for compiling with OpenMP
# AR        - Command to use when creating an archive (.a)
#
# 'psyclone' must also be on your PATH.

# Location of the dl_timer, infrastucture and OpenCL-wrapper code
# (all submodules in the git repo).
SHARED_DIR = ../../../../shared

TIMER_DIR = ${SHARED_DIR}/dl_timer
TIMER_INC = ${TIMER_DIR}/src
TIMER_LIB = ${TIMER_DIR}/dl_timer_lib.a
INF_DIR = ${SHARED_DIR}/dl_esm_inf/finite_difference
INF_INC = ${INF_DIR}/src
INF_LIB = ${INF_DIR}/src/lib_fd.a
FORTCL_INC = ${SHARED_DIR}/dl_esm_inf/external/FortCL/src/
FORTCL_LIB = ${SHARED_DIR}/dl_esm_inf/external/FortCL/src/libFortCL.a
COMMON_DIR = ../common
COMMON_LIB = ${COMMON_DIR}/nemolite2d_common.a


# The targets that this Makefile supports
.PHONY: all inf_lib timer_lib

# We don't include the OpenACC or OpenCL targets here as they
# have specialist compiler requirements so we require the user
# to explicitly request that they be built.
all: nemolite2d_serial

# The kernels used by this application and their location
KERNEL_DIR = ../kernels/fortran

# Shorthand for invoking PSyclone
PSYCLONE = psyclone -api gocean1.0 -d ${KERNEL_DIR} -l

# Serial PSyclone code-generation of NemoLite2D
nemolite2d_serial: timer_lib serial_inf_lib nemolite2d_alg.f90
	@echo "Generating NemoLite2D Serial"; rm -rf $@; mkdir -p $@
	${PSYCLONE} -nodm -oalg $@/alg.f90 -opsy $@/psy.f90 nemolite2d_alg.f90
	cp ${INF_LIB} $@/. # Each generated folder needs its purposely built inf_lib
	cp namelist $@/.
	cp Makefile_gen $@/Makefile
	make -C $@

# OpenMP PSyclone code-generation of NemoLite2D
nemolite2d_omp: timer_lib serial_inf_lib nemolite2d_alg.f90 omp_transform.py
	@echo "Generating NemoLite2D with OpenMP using omp_transform.py"; rm -rf $@; mkdir -p $@
	${PSYCLONE} -nodm -s ./omp_transform.py -oalg $@/alg.f90 -opsy $@/psy.f90 nemolite2d_alg.f90
	cp ${INF_LIB} $@/. # Each generated folder needs its purposely built inf_lib
	cp namelist $@/.
	cp Makefile_gen $@/Makefile
	make -C $@


# MPI PSyclone code-generation of NemoLite2D
nemolite2d_mpi: timer_lib mpi_inf_lib nemolite2d_alg.f90
	@echo "Generating NemoLite2D with MPI"; rm -rf $@; mkdir -p $@
	${PSYCLONE} -dm -oalg $@/alg.f90 -opsy $@/psy.f90 nemolite2d_alg.f90
	cp ${INF_LIB} $@/. # Each generated folder needs its purposely built inf_lib
	cp namelist $@/.
	cp Makefile_gen $@/Makefile
	make -C $@

# Hybrid PSyclone code-generation of NemoLite2D
nemolite2d_hybrid: timer_lib mpi_inf_lib nemolite2d_alg.f90 omp_transform.py
	@echo "Generating NemoLite2D with hybrid MPI+OpenMP using omp_transform.py"
	rm -rf $@; mkdir -p $@
	${PSYCLONE} -dm -s ./omp_transform.py -oalg $@/alg.f90 -opsy $@/psy.f90 nemolite2d_alg.f90
	cp ${INF_LIB} $@/. # Each generated folder needs its purposely built inf_lib
	cp namelist $@/.
	cp Makefile_gen $@/Makefile
	make -C $@


# OpenACC PSyclone code-generation of NemoLite2D
nemolite2d_acc: timer_lib serial_inf_lib nemolite2d_alg.f90 acc_transform.py
	@echo "Generating NemoLite2D with OpenACC"; rm -rf $@; mkdir -p $@
	${PSYCLONE} -nodm -s ./acc_transform.py -oalg $@/alg.f90 -opsy $@/psy.f90 nemolite2d_alg.f90
	cp ${INF_LIB} $@/. # Each generated folder needs its purposely built inf_lib
	cp namelist $@/.
	cp Makefile_gen $@/Makefile
	make -C $@

# OpenCL PSyclone code-generation of NemoLite2D
nemolite2d_ocl: timer_lib ocl_inf_lib nemolite2d_alg.f90 ocl_transform.py
	@echo "Generating NemoLite2D with OpenCL"; rm -rf $@; mkdir -p $@
	${PSYCLONE} -nodm -s ./ocl_transform.py -okern $@ -oalg $@/alg.f90 -opsy $@/psy.f90 nemolite2d_alg.f90
	cp ${INF_LIB} $@/. # Each generated folder needs its purposely built inf_lib
	cp namelist $@/.
	cp Makefile_gen $@/Makefile
	# Currently we need all kernels and the global parameters in a single OpenCL file
	cat kernel_parameters.h > $@/allkernels
	cat $@/*.cl >> $@/allkernels
	mv $@/allkernels $@/allkernels.cl 
	make -C $@ nemolite2docl


# Build the infrastructure libraries
timer_lib:
	${MAKE} -C ${TIMER_DIR} sm_lib

serial_inf_lib:
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${INF_DIR} fd_lib

mpi_inf_lib:
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${INF_DIR} dm_fd_lib

ocl_inf_lib:
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${INF_DIR} cl_fd_lib

$(COMMON_LIB):
	${MAKE} -C ${COMMON_DIR}


clean: 
	rm -rf nemolite2d_serial nemolite2d_omp nemolite2d_mpi
	rm -rf nemolite2d_hybrid nemolite2d_acc nemolite2d_ocl

libclean:
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean
	${MAKE} -C ${COMMON_DIR} allclean

allclean: clean libclean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe fparser.log
	rm -f *.lst *.cg *.opt *.optrpt *.s
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config
