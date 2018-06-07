# OpenMP, C version of Shallow #

This directory contains an OpenMP, C-version of the shallow-water
model 'shallow'.  It was provided by Graham Riley of the University
of Manchester.

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite.
If you are using Bash and the Gnu compiler then:

    > . ../../../compiler_setup/gnu.sh
    > make

will build the `shallow_omp.exe` executable.

Model parameters (including dimensions of the model grid) are set as
#define's in shallow_base_openmp_v3.c

## Running ##

Set the number of threads that the code will use as (in bash):
    > export OMP_NUM_THREADS=<num_threads>

and execute like so:

    > ./shallow_omp.exe

Note that for performance, you will probably need to ensure that threads
are bound to processor cores in a sensible fashion. See e.g. taskset, likwid
or Intels KMP_AFFINITY.

