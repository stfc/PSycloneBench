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

# Additionally it needs the path to the Kokkos source code
KOKKOS_PATH ?= $(HOME)/tools/kokkos
KOKKOS_DEBUG ?= no  # Careful, 10x performance penalty in kernels.
CXXFLAGS = $(CFLAGS)  # Use same CFLAGS to compile Kokkos library.

# If no KOKKOS_DEVICES is specified, by default use the OpenMP
KOKKOS_DEVICES ?= OpenMP

# Configuration specific to the selected Kokkos device 
ifeq ($(KOKKOS_DEVICES),OpenMP)
$(info "Using OpenMP device")
KOKKOS_OPTIONS := "aggresive_vectorization"
else ifeq ($(KOKKOS_DEVICES),Serial)
$(info "Using Serial device")
KOKKOS_OPTIONS := "aggresive_vectorization"
else ifeq ($(KOKKOS_DEVICES),Pthreads)
$(error "Pthreads device is currently not supported")
else ifeq ($(KOKKOS_DEVICES),Cuda)
$(info "Using CUDA device")
KOKKOS_ARCH = Pascal61
# nvcc doesn't support all compiler flags, overwrite CFLAGS to just O3
CFLAGS := -O3
# Still use the selected compiler but using the Kokkos nvcc_wrapper
NVCC_WRAPPER_DEFAULT_COMPILER = $(CXX)
CXX := $(KOKKOS_PATH)/bin/nvcc_wrapper
# The enable lambda option is necessary for the nvcc compiler to recognise
# as CUDA kernels the lambda-inlined functions.
KOKKOS_CUDA_OPTIONS = "enable_lambda"
# If CUDA_LIB is not provided, infer path from the nvcc compiler location.
CUDA_LIB ?= $(shell echo $(shell which nvcc) | sed 's/bin\/nvcc/lib64/g')
else
$(error "Unrecognised KOKKOS_DEVICES value: $(KOKKOS_DEVICES)")
endif

.PHONY: nemolite2d_views_kokkos nemolite2d_rawpointers_kokkos clean allclean summary
.DEFAULT_GOAL := nemolite2d_views_kokkos
include ../common.mk

all: nemolite2d_views_kokkos nemolite2d_rawpointers_kokkos

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}
COMMON_MODULES = time_step_mod.o ${COMMON_LIB} ${INF_LIB}

# Normal targets
nemolite2d_views_kokkos: ${COMMON_LIB} timer_lib inf_lib time_step_views_kokkos.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_views_kokkos.o" \
		$@_$(KOKKOS_DEVICES).exe

nemolite2d_rawpointers_kokkos: ${COMMON_LIB} timer_lib inf_lib time_step_rawpointers_kokkos.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_rawpointers_kokkos.o" \
		$@_$(KOKKOS_DEVICES).exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order
time_step_mod.o: ${COMMON_LIB} inf_lib

# Import needed kokkos environment variables
include $(KOKKOS_PATH)/Makefile.kokkos

# Linking the executable needs the -lstdc++ to link also the std lib functions
# called from C++ code
%.exe: $(MODULE_LIST) $(KOKKOS_LINK_DEPENDS) timing.o
	$(CXX) -o $@ timing.o $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) \
		$(LDFLAGS) ${OMPFLAGS} $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) -lgfortran

# Kokkos files rules
%_kokkos.o: %_kokkos.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) -DEXEC_SPACE=$(KOKKOS_DEVICES) $(KOKKOS_CXXFLAGS) $(CFLAGS) -c $<

timing.o: timing.cpp
	$(CXX) -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ *.dat
	rm -f gnu_opt_report.txt *.optrpt
	rm -rf KokkosCore_* Makefile.kokkos.f90

allclean: clean
	rm -f *.exe fparser.log *.a
	rm -rf *_wpl_dir/
	${MAKE} -C ${COMMON_DIR} allclean
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean

summary:
	@echo -n "${CC}\t\tcpu\t"
	@./nemolite2d_views_kokkos_OpenMP.exe | awk 'BEGIN {ORS="\t"} \
		{if ($$2 == "checksum") {print $$4} } \
		{if ($$1 == "Time-stepping") {print $$5} } \
		END {ORS="\n"; print ""}'
