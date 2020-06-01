# Makefile for manual C++ kokkos version of PSyKAl'd NEMOLite2D.
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

.PHONY: nemolite2d_kokkos
.DEFAULT_GOAL := nemolite2d_kokkos

# Import Kokkos build infrastructure
KOKKOS_PATH=${HOME}/build/kokkos
include $(KOKKOS_PATH)/Makefile.kokkos


include ../common.mk

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}
COMMON_MODULES = time_step_mod.o ${COMMON_LIB} ${KERNELS} ${INF_LIB}

# Normal targets
nemolite2d_kokkos: ${COMMON_LIB} timer_lib inf_lib time_step_kokkos.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_kokkos.o" $@.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order
time_step_mod.o: ${COMMON_LIB} inf_lib momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o

# Linking the executable needs the -lstdc++ to link also the std lib functions
# called from C++ code
%.exe: $(MODULE_LIST) $(KOKKOS_LINK_DEPENDS) 
	$(F90) -o $@ $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) \
		$(LDFLAGS) ${OMPFLAGS} $(KOKKOS_LINK_DEPENDS) $(KOKKOS_LDFLAGS) $(KOKKOS_LIBS) -lstdc++

# C++ rules
time_step_kokkos.o: time_step_kokkos.cpp $(KOKKOS_CPP_DEPENDS)
	$(CXX) $(KOKKOS_CPPFLAGS) $(KOKKOS_CXXFLAGS) $(CFLAGS) -c $<

clean: 
	${MAKE} -C ${INF_DIR} clean
	rm -f *.o *.mod *.MOD *~ *.dat
	rm -f gnu_opt_report.txt *.optrpt
	rm -rf KokkosCore_*

allclean: clean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe fparser.log *.a
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
