# Top-level Makefile for PSycloneBench benchmarks.
# Only supports those benchmarks that target the CPU (i.e. excluding
# OpenACC, OpenCL and Maxeler.)
#
# Picks-up the compiler and compiler flags from environment
# variables. See e.g. compiler_setup/gnu.sh

.PHONY: all shallow_cpu nemolite_cpu shallow_gen nemolite_gen

all: shallow_cpu nemolite_cpu

# All targets using PSyclone for code generation
all_gen: shallow_gen nemolite_gen

# All manual targets for CPU versions of Shallow
shallow_cpu:
	${MAKE} -C ./ocean/shallow/SEQ/original
	${MAKE} -C ./ocean/shallow/SEQ shallow_base
	${MAKE} -C ./ocean/shallow/OMP

# Requires PSyclone be installed
shallow_gen:
	${MAKE} -C ./ocean/shallow/SEQ shallow_gen

# All manual targets for CPU versions of NEMOLite2D
nemolite_cpu:
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/psykal_serial
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/psykal_omp
	${MAKE} -C ./ocean/nemo/nemolite2d/original

# Requires PSyclone be installed
nemolite_gen:
	${MAKE} -C ./ocean/nemo/nemolite2d/psykal

nemolite_acc:
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/single_file_acc
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/psykal_acc

clean:
	${MAKE} -C ./ocean/nemo/nemolite2d/psykal clean
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/psykal_serial clean
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/psykal_omp clean
	${MAKE} -C ./ocean/nemo/nemolite2d/original clean
	${MAKE} -C ./ocean/shallow/SEQ clean
	${MAKE} -C ./ocean/shallow/SEQ/original clean
	${MAKE} -C ./ocean/shallow/OMP clean

allclean:
	${MAKE} -C ./ocean/nemo/nemolite2d/psykal allclean
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/psykal_serial allclean
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/psykal_omp allclean
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/single_file_acc allclean
	${MAKE} -C ./ocean/nemo/nemolite2d/manual_versions/psykal_acc allclean
	${MAKE} -C ./ocean/nemo/nemolite2d/original allclean
	${MAKE} -C ./ocean/shallow/SEQ allclean
	${MAKE} -C ./ocean/shallow/SEQ/original allclean
	${MAKE} -C ./ocean/shallow/OMP allclean

