# Makefile for manual C++ version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler to use any flags from
# environment variables. e.g. to use gfortran:
# 
# export F90=gfortran
# export CXX=g++
# export CFLAGS="-O3"
# export F90FLAGS="-O3"
# export OMPFLAGS=" -fopenmp"
# export AR=ar
SYCL_COMP ?= dpcpp
SYCL_FLAGS ?= -fsycl -O3 -mtune=native

.PHONY: all nemolite2d_sycl
.DEFAULT_GOAL := nemolite2d_sycl

include ../common.mk

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}
COMMON_MODULES = time_step_mod.o ${COMMON_LIB} ${INF_LIB}

# Normal targets
all: nemolite2d_sycl

nemolite2d_sycl: ${COMMON_LIB} timer_lib inf_lib time_step_sycl.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_sycl.o" nemolite2d_sycl.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order
time_step_mod.o: ${COMMON_LIB} inf_lib

# Linking the executable needs the libgfortran runtime library installed in the system
# and available in the library path.
%.exe: $(MODULE_LIST) timing.o
	$(SYCL_COMP) -o $@ timing.o $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) \
		$(LDFLAGS) ${OMPFLAGS}  -lgfortran -lsycl

# C++ rules (be aware of filename clashes as there is already a %.o: %.f90 in
# the imported ../common.mk)
# Since kernels are included directly any modification in the kernel files
# must trigger a recompilation.
%.o: %.cpp ../../kernels/c_family/*.c
	$(SYCL_COMP) $(SYCL_FLAGS)  -fPIC -c $< -o $@

timing.o: timing.cpp
	$(CXX) -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ *.dat
	rm -f gnu_opt_report.txt *.optrpt

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/
	${MAKE} -C ${COMMON_DIR} allclean
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean

summary:
	@echo -n "${CC}\t\tcpu\t"
	@./nemolite2d.exe | awk 'BEGIN {ORS="\t"} \
		{if ($$2 == "checksum") {print $$4} } \
		{if ($$1 == "Time-stepping") {print $$5} } \
		END {ORS="\n"; print ""}'
