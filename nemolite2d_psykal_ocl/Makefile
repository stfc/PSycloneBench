# You must define:
#
#   CC - the C compiler
#   OPENCL_INCLUDE - root of path to search for OpenCL headers
#
# Optionally:
#
#   CFLAGS - flags to pass to the C compiler
#   LDFLAGS - flags to pass to the C linker

# The targets that this Makefile supports
EXECS = nemolite2d

all: $(EXECS)

COMMON_MODULES = opencl_utils.o timing.o

# Source files containing (OpenCL) kernels
KERNEL_MODULES = continuity_kern.o momentum_kern.o \
                 boundary_conditions_kern.o time_update_kern.o

KERNEL_BINARIES = $(KERNEL_MODULES:.o=.aocx)

AOC_FLAGS = -v -board=p385a_sch_ax115
#AOC_FLAGS += -march=emulator

# Normal targets
.PHONY: nemolite2d
nemolite2d: 
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} ${KERNEL_MODULES}" nemolite2d.exe

# Target for FPGA - we must compile the kernels using the off-line compiler
nemolite2d_fpga: nemolite2d_kernels.aocx nemolite2d

nemolite2d.o: $(COMMON_MODULES) ${KERNEL_MODULES}

# Compile all the kernels into a single file so that we can put them all on
# one device
nemolite2d_kernels.aocx: ${KERNEL_MODULES}
	aoc -D__OPENCL_VERSION__=100 ${AOC_FLAGS} ${KERNEL_MODULES:.o=.cl} -o $@

# Interdependencies between modules, alphabetical order

# Generic rules

%.exe: $(MODULE_LIST)
	$(CC) -o $@ $(MODULE_LIST) ${OPENCL_LIBS} $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) ${OPENCL_INCLUDE} -c $<

# Compile the kernels for emulation on the CPU
%.aocx: %.cl
	aoc -D__OPENCL_VERSION__=100 ${AOC_FLAGS} $< -o $@

clean: 
	rm -f *.o *~

allclean: clean
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/

# When we've gone through the pain of building kernels for real FPGA
# hardware then we're not often going to want to delete them - hence
# a separate target for that...
kernclean:
	rm -rf $(KERNEL_MODULES:.o=)
	rm -f ${KERNEL_BINARIES}
	rm -f $(KERNEL_MODULES:.o=.aoco)
	rm -f nemolite2d_kernels.aoc?
