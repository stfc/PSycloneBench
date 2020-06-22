# Makefile for manual C++ kokkos version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler to use plus any flags from
# environment variables. e.g. to use g++:
# 
# export F90=gfortran
# export CXX=g++
# export CFLAGS="-O3"
# export F90FLAGS="-O3"
# export OMPFLAGS=" -fopenmp"
# export AR=ar

KOKKOS_PATH ?= $(HOME)/tools/kokkos
KOKKOS_DEBUG ?= yes

# By default us the OpenMP device
ifndef $(KOKKOS_DEVICES)
KOKKOS_DEVICES=OpenMP
endif

# Configuration specific to the selected Kokkos device 
ifeq ($(KOKKOS_DEVICES),OpenMP)
$(info "Using OpenMP device")
EXEC_SPACE = OpenMP
else ifeq ($(KOKKOS_DEVICES),Cuda)
$(info "Using CUDA device")
KOKKOS_ARCH = Pascal61
EXEC_SPACE = Cuda
# nvcc fails with some complex flags, overwirte CFLAGS to just O3
CFLAGS := "-O3"
# Still use the selected compiler but using the Kokkos nvcc_wrapper
NVCC_WRAPPER_DEFAULT_COMPILER = $(CXX)
CXX := $(KOKKOS_PATH)/bin/nvcc_wrapper
KOKKOS_CUDA_OPTIONS = "enable_lambda"
else
$(error "Unrecognised KOKKOS_DEVICES value: $(KOKKOS_DEVICES)")
endif

.PHONY: nemolite2d_kokkos
.DEFAULT_GOAL := nemolite2d_kokkos
include ../common.mk

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}
COMMON_MODULES = time_step_mod.o ${COMMON_LIB} ${KERNELS} ${INF_LIB}

# Normal targets
nemolite2d_kokkos: ${COMMON_LIB} timer_lib inf_lib time_step_kokkos.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_kokkos.o" $@.exe

nemolite2d_kokkos_rawpointers: ${COMMON_LIB} timer_lib inf_lib time_step_kokkos_rawpointers.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_kokkos_rawpointers.o" $@.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order
time_step_mod.o: ${COMMON_LIB} inf_lib momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o

# Import needed kokkos environment variables
include $(KOKKOS_PATH)/Makefile.kokkos

# Linking the executable needs the -lstdc++ to link also the std lib functions
# called from C++ code
%.exe: $(MODULE_LIST) $(KOKKOS_LINK_DEPENDS) 
	echo  $(KOKKOS_LINK_DEPENDS)
	echo  $(KOKKOS_LDFLAGS)
	$(F90) -o $@ $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) \
		$(LDFLAGS) ${OMPFLAGS} $(KOKKOS_LINK_DEPENDS) $(KOKKOS_LIBS) \
		-L/opt/nvidia/cuda-10/lib64 -L. -lstdc++

# Kokkos files rules
time_step_kokkos.o: time_step_kokkos.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) -DEXEC_SPACE=$(EXEC_SPACE) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CFLAGS) $(OMPFLAGS) -c $<

time_step_kokkos_rawpointers.o: time_step_kokkos_rawpointers.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) -DEXEC_SPACE=$(EXEC_SPACE) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CFLAGS) $(OMPFLAGS) -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ *.dat
	rm -f gnu_opt_report.txt *.optrpt
	rm -rf KokkosCore_* Makefile.kokkos.f90

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe fparser.log *.a
	rm -rf *_wpl_dir/
	${MAKE} -C ${COMMON_DIR} allclean
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean

summary:
	@echo -n "${CC}\t\tcpu\t"
	@./nemolite2d_kokkos.exe | awk 'BEGIN {ORS="\t"} \
		{if ($$2 == "checksum") {print $$4} } \
		{if ($$1 == "Time-stepping") {print $$5} } \
		END {ORS="\n"; print ""}'
