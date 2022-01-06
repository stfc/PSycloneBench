# Tracer-advection Benchmark(s) #

Below this directory are various versions of the tracer-advection
benchmark.

The `scripts` directory contains PSyclone transformation scripts that
are used to convert the benchmark into SIR-compliant form (see below)
and to add OpenACC Kernels directives. It also contains a simple
`problemsize.sh` script which automates the execution of the benchmark
for a range of problem sizes.

## original ##

The `original` directory contains the single-file version of the
benchmark based on that extracted from the NEMO code by Silvia
Mocavero of CMCC. The Makefile provides targets to generate serial,
OpenACC and OpenMP versions of this benchmark. The latter two use
PSyclone to transform the code.

## compute_in_subroutine ##

This directory contains a version of the benchmark where all of the
compute has been moved into a subroutine contained within a separate
Fortran source file. Two versions of this are then provided, which
differ in the way work arrays within the compute subroutine are
implemented. By default these arrays are 'automatic' (and thus
allocated by the compiler upon each call of the subroutine). In a
second version, this repeated allocation is avoided by promoting
the variables to have module scope. For each of these versions,
PSyclone is used to produce an OpenACC version, an SIR-compliant
version (with certain intrinsics and code structures replaced)
and an OpenACC+SIR-compliant version. (SIR is the Stencil
Intermediate Representation used by the Dawn/Dusk toolchain.)

## multi_kernel ##

This directory contains a version of the benchmark where the compute
has been further sub-divided into ten subroutines. This form can then
be processed by PSyclone to produce code that works with the Dawn/Dusk
toolchain (allowing CUDA to be produced). As with the `compute_in_subroutine`
version, PSyclone is used to produce OpenACC, SIR-compliant and
OpenACC+SIR-compliant versions.

