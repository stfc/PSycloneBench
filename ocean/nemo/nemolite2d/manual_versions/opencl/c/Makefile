# You must define:
#
#   CC - the C compiler
#   OPENCL_INCLUDE - root of path to search for OpenCL headers
#   OPENCL_LIBS - necessary paths and libs to link to OpenCL
#
# Optionally:
#
#   CFLAGS - flags to pass to the C compiler
#   LDFLAGS - flags to pass to the C linker

# The targets that this Makefile supports
EXECS = nemolite2d

.PHONY: all nemolite2d fpga cont twokern twokernchannels

all: $(EXECS)

COMMON_MODULES = opencl_utils.o timing.o nemolite2d_utils.o

# Source files containing (OpenCL) kernels
KERNEL_DIR = ../../../kernels/opencl
ALL_KERNEL_MODULES = continuity_kern.o momentum_u_kern.o momentum_v_kern.o \
                 boundary_conditions_kern.o time_update_kern.o

# Flags for the Altera Offline Compiler. Need the -D__OPENCL_VERSION__
# because the kernel source files also include CPU code
AOC_FLAGS = -v -board=p385a_sch_ax115 -D__OPENCL_VERSION__=100
#AOC_FLAGS += -c
AOC_FLAGS += -march=emulator
#-emulator-channel-depth-model=strict

# Normal targets
nemolite2d: 
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES}" KERNEL_MODULES="${ALL_KERNEL_MODULES}" nemolite2d.exe

# Single-kernel version of the application
cont:
	${MAKE} CFLAGS="${CFLAGS} -DSINGLE_KERNEL" MODULE_LIST="cont_only.o ${COMMON_MODULES}" KERNEL_MODULES="continuity_kern.o" cont.exe

cont_kern: continuity_kern.cl
	aoc ${AOC_FLAGS} continuity_kern.cl

# Two-kernel version
twokern:
	${MAKE} KERNEL_MODULES="continuity_kern.o momentum_u_kern.o" fpga
	${MAKE} CFLAGS="${CFLAGS}" MODULE_LIST="twokern.o ${COMMON_MODULES}" KERNEL_MODULES="continuity_kern.o momentum_u_kern.o" twokern.exe

# Two-kernel version with streaming using FIFO buffers
twokernchannels:
	${MAKE} KERNEL_MODULES="continuity_kern_channel.o time_update_u_kern_channel.o" fpga
	${MAKE} CFLAGS="${CFLAGS}" MODULE_LIST="twokernchannels.o ${COMMON_MODULES}" KERNEL_MODULES="continuity_kern.o time_update_u_kern.o" twokernchannels.exe

simplechannel:
	${MAKE} KERNEL_MODULES="read_write_kern.o" fpga
	${MAKE} CFLAGS="${CFLAGS}" MODULE_LIST="simplechannel.o ${COMMON_MODULES}" KERNEL_MODULES="" simplechannel.exe

# Target for FPGA - we must compile the kernels using the off-line compiler
nemolite2d_fpga:
	${MAKE} KERNEL_MODULES="${ALL_KERNEL_MODULES}" fpga
	${MAKE} KERNEL_MODULES="${ALL_KERNEL_MODULES}" nemolite2d

nemolite2d.o: $(COMMON_MODULES) ${KERNEL_MODULES}

# Compile all the kernels into a single file so that we can put them all on
# one device
fpga:
	${MAKE} KERNEL_MODULES="${KERNEL_MODULES}" nemolite2d_kernels.aocx

nemolite2d_kernels.aocx: ${KERNEL_MODULES:.o=.cl}
	aoc ${AOC_FLAGS} ${KERNEL_MODULES:.o=.cl} -o nemolite2d_kernels

# Generic rules

%.exe: $(MODULE_LIST) ${KERNEL_MODULES}
	$(CC) -o $@ $(MODULE_LIST) ${KERNEL_MODULES} ${OPENCL_LIBS} $(LDFLAGS)

# The Altera Off-line Compiler expects sources to have the .cl suffix
%.cl: %.c
	ln -sf $< $@

# Kernel source is kept in the top-level kernels directory
%.c:
	ln -sf ${KERNEL_DIR}/$@ .

%.o: %.c
	$(CC) $(CFLAGS) ${OPENCL_INCLUDE} -c $<

%.o: %.cpp
	$(CC) $(CFLAGS) ${OPENCL_INCLUDE} -c $<

# Compile the kernels for emulation on the CPU
%.aocx: %.cl
	aoc -D__OPENCL_VERSION__=100 ${AOC_FLAGS} $< -o $@

clean: 
	rm -f *.o *~ *.cl

allclean: clean
	rm -f *.exe fparser.log gnu_opt_report.txt
	rm -rf *_wpl_dir/

# When we've gone through the pain of building kernels for real FPGA
# hardware then we're not often going to want to delete them - hence
# a separate target for that...
kernclean:
	rm -rf $(ALL_KERNEL_MODULES:.o=)
	rm -f *.aoc?
