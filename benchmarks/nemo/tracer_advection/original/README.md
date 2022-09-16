# README #

This directory contains the 'original' tra_adv benchmark as extracted from
the NEMO source code by Silvia Mocavero of CMCC.

## Building ##

The code uses the [dl_timer](https://bitbucket.org/apeg/dl_timer)
library for timing. Provided you followed the steps in the top-level
[README](../../../../README.md#obtaining-the-code) then the source code
for this library will have been obtained when you cloned the
repository. The library is built automatically as part of the
compilation process.

The Makefile supports the following targets:

* tra_adv_serial         - the original, sequential form of the benchmark.
* tra_adv_acc_kernels    - version transformed by PSyclone for OpenACC acceleration
                           on GPU using the kernels directive.
* tra_adv_acc_loops      - version transformed by PSyclone for OpenACC acceleration
                           on GPU using the explicit loop directive.
* tra_adv_omp_cpu_levels - version transformed by PSyclone for OpenMP threading
                           on CPU parallelising loops over the k-domain.
* tra_adv_omp_cpu        - version transformed by PSyclone for OpenMP threading
                           on CPU parallelising the outermost loop that is
                           parallelisable.
* tra_adv_omp_gpu        - version transformed by PSyclone for OpenMP offload to
                           GPU.

The compiler and compiler flags to use must be set through the following
environment variables:

* F90      - the command with which to invoke the Fortran compiler
* F90FLAGS - flags to pass to the compiler, e.g. -g
* OMPFLAGS - the flag(s) required to enable OpenMP with the chosen compiler
             (if desired)

e.g. to build with Gnu Fortran one might use:

    export F90=gfortran
    export F90FLAGS=-O3
    (export OMPFLAGS=-fopenmp)
    export MPIF90=mpif90

To use the NVIDIA compiler, OpenMP offload and managed memory:

    export F90=nvfortran
    export F90FLAGS="-O3"
    export OMPTARGETFLAGS="-mp=gpu"
    export UMEMFLAGS="-gpu=managed"

Scripts to do this for various compilers may be found in the `compiler_setup`
directory at the root of this repository.

Additionally, the optional `ENABLE_NVIDIA_PROFILE=yes` flag can be set to enable
GPU offloading versions to generate instrumented code, this will require the
`PSYCLONE_NVIDIA_LIB_DIR` to be set up to reference the path to PSyclone's nvidia
profiling lib (see
https://psyclone.readthedocs.io/en/stable/profiling.html#interface-to-third-party-profiling-tools).

## Running ##

The problem size and number of iterations/time-steps to perform are also
set using environment variables:

* JPI - x-dimension of the domain
* JPJ - y-dimension of the domain
* JPK - number of vertical levels
* IT - number of iterations to perform

e.g.:

    export JPI=130
    export JPJ=130
    export JPK=31
    export IT=100
