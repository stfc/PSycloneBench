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

# Normal targets
.PHONY: nemolite2d
nemolite2d: 
	${MAKE} MODULE_LIST="nemolite2d.o ${COMMON_MODULES} ${KERNEL_MODULES}" nemolite2d.exe

nemolite2d.o: $(COMMON_MODULES) ${KERNEL_MODULES}

# Interdependencies between modules, alphabetical order

# Generic rules

%.exe: $(MODULE_LIST)
	$(CC) -o $@ $(MODULE_LIST) ${OPENCL_LIBS} $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) ${OPENCL_INCLUDE} -c $<

clean: 
	rm -f *.o *~

allclean: clean
	rm -f *.exe fparser.log
	rm -rf *_wpl_dir/
