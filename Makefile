# Top-level Makefile for PSycloneBench benchmarks.
# By default only builds those benchmarks that target the CPU (i.e. excluding
# OpenACC, OpenCL and Maxeler.) Separate targets for the OpenACC versions
# are provided.
#
# Picks-up the compiler and compiler flags from environment
# variables. See e.g. compiler_setup/gnu.sh

.PHONY: all shallow_cpu nemolite_cpu shallow_gen nemolite_gen

all: shallow_cpu nemolite_cpu

# All targets using PSyclone for code generation
all_gen: shallow_gen nemolite_gen nemolite_cpu

# All manual targets for CPU versions of Shallow
shallow_cpu:
	${MAKE} -C ./benchmarks/shallow/SEQ/original
	${MAKE} -C ./benchmarks/shallow/SEQ shallow_base
	${MAKE} -C ./benchmarks/shallow/OMP

# Requires PSyclone be installed
shallow_gen:
	${MAKE} -C ./benchmarks/shallow/SEQ shallow_gen

# All manual targets for CPU versions of NEMOLite2D
nemolite_cpu:
	${MAKE} -C ./benchmarks/nemo/nemolite2d/manual_versions

# Requires PSyclone be installed
nemolite_gen:
	${MAKE} -C ./benchmarks/nemo/nemolite2d/psykal

nemolite_acc:
	${MAKE} -C ./benchmarks/nemo/nemolite2d/manual_versions/single_file_acc
	${MAKE} -C ./benchmarks/nemo/nemolite2d/manual_versions/psykal_acc

clean allclean:
	${MAKE} -C ./benchmarks/nemo/nemolite2d/psykal $@
	${MAKE} -C ./benchmarks/nemo/nemolite2d/common $@
	${MAKE} -C ./benchmarks/nemo/nemolite2d/manual_versions $@
	${MAKE} -C ./benchmarks/nemo/nemolite2d/original $@
	${MAKE} -C ./benchmarks/shallow/SEQ $@
	${MAKE} -C ./benchmarks/shallow/SEQ/original $@
	${MAKE} -C ./benchmarks/shallow/OMP $@
	${MAKE} -C ./benchmarks/smallmatvec/manual_versions/openmp $@

