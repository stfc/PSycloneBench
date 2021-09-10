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

.PHONY: all nemolite2d nemolite2d_omp nemolite2d_ompopt
.DEFAULT_GOAL := nemolite2d

include ../common.mk

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}
COMMON_MODULES = time_step_mod.o ${COMMON_LIB} ${KERNELS} ${INF_LIB}

# Normal targets
all: nemolite2d nemolite2d_omp

nemolite2d: ${COMMON_LIB} timer_lib inf_lib time_step.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step.o" nemolite2d.exe

nemolite2d_omp: ${COMMON_LIB} timer_lib inf_lib
	rm -f time_step_omp.o
	${MAKE} CFLAGS="${CFLAGS} ${OMPFLAGS}" LDFLAGS="${LDFLAGS} ${OMPFLAGS}" MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_omp.o" nemolite2d_omp.exe

nemolite2d_omp_gpu: ${COMMON_LIB} timer_lib inf_lib
	rm -f time_step_omp.o
	${MAKE} CFLAGS="${CFLAGS} ${OMPTARGETFLAGS}" LDFLAGS="${LDFLAGS} ${OMPTARGETFLAGS}" MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_omp.o" nemolite2d_omp_gpu.exe

nemolite2d_ompopt: ${COMMON_LIB} timer_lib inf_lib time_step_ompopt.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_ompopt.o" nemolite2d_ompopt.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order
time_step_mod.o: ${COMMON_LIB} inf_lib momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o

# Linking the executable needs the -lstdc++ to link also the std lib functions
# called from C++ code
%.exe: $(MODULE_LIST)
	$(F90) -o $@ $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) \
		$(LDFLAGS) -lstdc++

# C++ rules (be aware of filename clashes as there is already a %.o: %.f90 in
# the imported ../common.mk)
# Since kernels are included directly any modification in the kernel files
# must trigger a recompilation.
%.o: %.cpp ../../kernels/c_family/*.c
	$(CXX) $(CFLAGS) -c $< -o $@
	
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
