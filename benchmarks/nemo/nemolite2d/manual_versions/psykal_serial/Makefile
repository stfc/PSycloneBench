# This Makefile expects the following env. vars to be set:
#
# F90       - How to invoke the Fortran compiler
# F90FLAGS  - Flags to pass to the Fortran compiler
# OMPFLAGS  - Flags for compiling with OpenMP
# AR        - Command to use when creating an archive (.a)

.PHONY: all nemolite2d nemolite2d_align nemolite2d_cont_only nemolite2d_nodiv nemolite2d_preload nemolite2d_alignall nemolite2d_nopeel nemolite2d_8arrays
.DEFAULT_GOAL := nemolite2d

include ../common.mk

# INF_LIB is an archive that must come at the end of the list of objects
# passed to the linker
LIBS = ${COMMON_LIB} ${TIMER_LIB} ${INF_LIB}
COMMON_MODULES = ${COMMON_LIB} ${KERNELS} ${INF_LIB}

# Interdependencies between modules, alphabetical order
time_step_mod.o: ${COMMON_LIB} ${KERNELS} inf_lib

IMPLEMENTATIONS=nemolite2d nemolite2d_align nemolite2d_cont_only nemolite2d_nodiv \
nemolite2d_preload nemolite2d_alignall nemolite2d_nopeel nemolite2d_8arrays

# The targets that this Makefile supports
all: ${IMPLEMENTATIONS}

# Normal targets. Need to link with OpenMP in order to use OMP timer
nemolite2d: ${COMMON_MODULES} time_step_mod.o nemolite2d.o
	${MAKE} MODULE_LIST="nemolite2d.o time_step_mod.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d.exe

# Version of nemolite2d with the time-stepping routine hacked to
# encourage the Intel compiler to vectorise efficiently
nemolite2d_align: ${COMMON_MODULES} time_step_mod_align.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align
	rm -f time_step_mod.o time_step_mod.mod
	${MAKE} MODULE_LIST="time_step_mod_align.o nemolite2d.o ${COMMON_MODULES} ${LIBS}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_align.exe

# Time-stepping routine hacked to only call the Continuity kernel
nemolite2d_cont_only: ${COMMON_MODULES} time_step_mod_cont_only.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align
	rm -f time_step_mod.o time_step_mod.mod
	${MAKE} MODULE_LIST="time_step_mod_cont_only.o nemolite2d.o ${COMMON_MODULES} ${LIBS}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_cont_only.exe

# With division ops removed from Momentum kernels
nemolite2d_nodiv: ${COMMON_MODULES} time_step_mod_align_nodiv.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_nodiv
	rm -f time_step_mod.o time_step_mod.mod
	${MAKE} MODULE_LIST="time_step_mod_align_nodiv.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_nodiv.exe

nemolite2d_preload: ${COMMON_MODULES} time_step_mod_align_preload.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_preload
	rm -f time_step_mod.o time_step_mod.mod
	${MAKE} MODULE_LIST="time_step_mod_align_preload.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_preload.exe

nemolite2d_alignall: ${COMMON_MODULES} time_step_mod_alignall_preload.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_preload
	rm -f time_step_mod*.o time_step_mod*.mod
	${MAKE} MODULE_LIST="time_step_mod_alignall_preload.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_alignall.exe

nemolite2d_nopeel: ${COMMON_MODULES} time_step_mod_alignall_preload_nopeel.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_preload
	rm -f time_step_mod*.o time_step_mod*.mod
	${MAKE} MODULE_LIST="time_step_mod_alignall_preload_nopeel.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_nopeel.exe

nemolite2d_8arrays: ${COMMON_MODULES} time_step_mod_alignall_preload_nopeel_8arrays.o nemolite2d.o
        # Have to delete any existing time_step_mod because it contains a module
        # with the same name as that in time_step_mod_align_preload
	rm -f time_step_mod*.o time_step_mod*.mod
	${MAKE} MODULE_LIST="time_step_mod_alignall_preload_nopeel_8arrays.o nemolite2d.o ${COMMON_MODULES}" \
             LDFLAGS="${LDFLAGS} ${OMPFLAGS}" nemolite2d_8arrays.exe


clean: 
	rm -f *.o *.mod *.MOD *~ go2d_* tmask_*

libclean:
	${MAKE} -C ${INF_DIR} distclean
	${MAKE} -C ${TIMER_DIR} allclean

allclean: clean libclean
	rm -f ${KERNELS:.o=.f90}
	rm -f *.exe fparser.log
	rm -f *.lst *.cg *.opt *.optrpt gnu_opt_report.txt *.s
	rm -rf *_wpl_dir/

docs:
	doxygen gocean2d.doxy.config

summary:
	@echo -n "${CC}\t\tcpu\t"
	@./nemolite2d.exe | awk 'BEGIN {ORS="\t"} \
		{if ($$2 == "checksum") {print $$4} } \
		{if ($$1 == "Time-stepping") {print $$5} } \
		END {ORS="\n"; print ""}'

# Execute all enabled implementations with their best full node performance and
# report their ua and uv final checksums and the time-stepping time.
summary_table:
	@echo "Implementation\t\tCompiler\tArch\tua checksum\tuv checksum\ttime/step"
	@$(foreach var,$(IMPLEMENTATIONS), \
		echo -n "$(var)\t${CC}\t\tcpu\t" && \
		./$(var).exe | awk 'BEGIN {ORS="\t"} \
			{if ($$2 == "checksum") {print $$4} } \
			{if ($$1 == "Time-stepping") {print $$5} } \
			END {ORS="\n"; print ""}'; \
	)
