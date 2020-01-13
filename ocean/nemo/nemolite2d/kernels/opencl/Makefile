# Makefile for the OpenCL kernels for NEMOLite2D.
#
# You must define:
#
#   CC - the C compiler
#   OPENCL_INCLUDE - root of path to search for OpenCL headers
#
# Optionally:
#
#   CFLAGS - flags to pass to the C compiler
#   LDFLAGS - flags to pass to the C linker

# Source files containing (OpenCL) kernels
ALL_KERNEL_MODULES = continuity_kern.o momentum_u_kern.o momentum_v_kern.o \
                 boundary_conditions_kern.o time_update_kern.o

# Flags for the Altera Offline Compiler. Need the -D__OPENCL_VERSION__
# because the kernel source files also include CPU code
AOC_FLAGS = -v -board=p385a_sch_ax115 -D__OPENCL_VERSION__=100
#AOC_FLAGS += -c
AOC_FLAGS += -march=emulator
#-emulator-channel-depth-model=strict
AOC_FLAGS += -I../c

# This causes us to use the chopped-down version of the u/v Momentum kernels
# which simplifies testing with the FPGA.
#AOC_FLAGS += "-DSIMPLE_MOMENTUM"
#CFLAGS += "-DSIMPLE_MOMENTUM"

all: fpga

# Compile all the kernels into a single file so that we can put them all on
# one device
fpga:
	${MAKE} KERNEL_MODULES="${ALL_KERNEL_MODULES}" nemolite2d_kernels.aocx

nemolite2d_kernels.aocx: ${KERNEL_MODULES:.o=.cl}
	aoc ${AOC_FLAGS} ${KERNEL_MODULES:.o=.cl} -o nemolite2d_kernels

# Generic rules

# The Altera Off-line Compiler expects sources to have the .cl suffix
%.cl: %.c
	ln -sf $< $@

%.o: %.c
	$(CC) $(CFLAGS) -I../c ${OPENCL_INCLUDE} -c $<

# Compiles individual kernels
%.aocx: %.cl
	aoc -D__OPENCL_VERSION__=100 ${AOC_FLAGS} $< -o $@

clean: 
	rm -f *.o *~ *.cl

# When we've gone through the pain of building kernels for real FPGA
# hardware then we're not often going to want to delete them - hence
# a separate target for that...
kernclean:
	rm -rf $(ALL_KERNEL_MODULES:.o=)
	rm -f *.aoc?
