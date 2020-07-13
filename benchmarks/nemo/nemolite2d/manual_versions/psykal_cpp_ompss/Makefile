# Makefile for manual C++ OmpSs version of PSyKAl'd NEMOLite2D.
#
# This Makefile picks up the compiler to use any flags from
# environment variables. e.g. to use gfortran:
# 

.PHONY: all nemolite2d_ompss
.DEFAULT_GOAL := nemolite2d_ompss

include ../common.mk

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}
COMMON_MODULES = time_step_mod.o ${COMMON_LIB} ${KERNELS} ${INF_LIB}

# Normal targets
all: nemolite2d_ompss

nemolite2d_ompss: ${COMMON_LIB} timer_lib_parallel inf_lib time_step_ompss.o
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} time_step_ompss.o" nemolite2d_ompss.exe

nemolite2d.o: $(COMMON_MODULES)

# Interdependencies between modules, alphabetical order
time_step_mod.o: ${COMMON_LIB} inf_lib momentum_mod.o continuity_mod.o \
                 time_update_mod.o boundary_conditions_mod.o

# Linking the executable needs the -lstdc++ to link also the std lib functions
# called from C++ code

nemolite2d_ompss.exe: $(MODULE_LIST)
	mcxx --ompss -o $@ $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) \
              $(LDFLAGS) -lgfortran -lstdc++

%.exe: $(MODULE_LIST)
	$(CXX) -o $@ $(MODULE_LIST) ${COMMON_LIB} ${INF_LIB} $(TIMER_LIB) \
		$(LDFLAGS) ${OMPFLAGS} -lgfortran -lstdc++

# C++ rules (be aware of filename clashes as there is already a %.o: %.f90 in
# the imported ../common.mk)

time_step_ompss.o: time_step_ompss.cpp
	mcxx --ompss $(MCXXFLAGS) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CFLAGS) $(OMPFLAGS) -c $< -o $@

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
	@./nemolite2d_ompss.exe | awk 'BEGIN {ORS="\t"} \
		{if ($$2 == "checksum") {print $$4} } \
		{if ($$1 == "Time-stepping") {print $$5} } \
		END {ORS="\n"; print ""}'
