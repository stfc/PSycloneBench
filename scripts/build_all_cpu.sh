#!/bin/bash

# Simple script that builds all of the benchmark codes that run on CPU.
# Uses the Gnu compiler suite (needs gcc and gfortran).

. compiler_setup/gnu.sh

# The Shallow benchmarks
make -C ./ocean/shallow/SEQ/original
make -C ./ocean/shallow/SEQ shallow_base
make -C ./ocean/shallow/OMP
# The NEMOLite2D benchmarks
make -C  ./ocean/nemo/nemolite2d/manual_versions/psykal_serial
make -C  ./ocean/nemo/nemolite2d/manual_versions/psykal_omp
make -C  ./ocean/nemo/nemolite2d/original
# Versions which require that PSyclone be available
if hash psyclone 2>/dev/null; then
  make -C ./ocean/shallow/SEQ shallow_gen
  make -C ./ocean/nemo/nemolite2d/psykal
fi
