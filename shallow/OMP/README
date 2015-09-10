This directory contains an OpenMP, C-version of the shallow-water model 'shallow'.
To build:
  ln -sf Makefile.include.<your compiler> Makefile.include
  make

You will probably have to edit/create Makefile.include.<your compiler> to match
your system.

Set the number of threads that the code will use as (in bash):
  export OMP_NUM_THREADS=<num_threads>

Model parameters (including dimensions of the model grid) are set as
#define's in shallow_base_openmp_v3.c

Andrew Porter
STFC Daresbury Laboratory
