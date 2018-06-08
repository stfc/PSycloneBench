#!/bin/bash

# Simple script that builds all of the benchmark codes that run on CPU.
# Uses the Gnu compiler suite (needs gcc and gfortran).

. compiler_setup/gnu.sh

# Create list of locations containing benchmarks that compile for CPU
# (i.e. excluding those using OpenACC, OpenCL or Maxeler)
dirs="./ocean/shallow/SEQ/original"
dirs+=" ./ocean/shallow/SEQ"
dirs+=" ./ocean/shallow/OMP"
dirs+=" ./ocean/nemo/nemolite2d/manual_versions/psykal_serial"
dirs+=" ./ocean/nemo/nemolite2d/manual_versions/psykal_omp"
dirs+=" ./ocean/nemo/nemolite2d/original"
if hash psyclone 2>/dev/null; then
  # Versions which require that PSyclone be available
  dirs+=" ./ocean/nemo/nemolite2d/psykal"
else
  echo "PSyclone not found so skipping benchmarks requiring code generation"
fi

for dir in $dirs
do
  cd $dir
  make
  cd -
done
