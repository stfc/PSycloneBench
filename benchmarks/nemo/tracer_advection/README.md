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
Fortran source file. The reason for this change is that the Stencil
Intermediate Representation (SIR) used by the Dawn/Dusk toolchain is
only able to represent computation and thus cannot be used for the
full application. With the compute isolated in a subroutine, the SIR
for that routine can be constructed.

Two versions of this form of the benchmark are then provided, which
differ in the way work arrays within the compute subroutine are
implemented. By default these arrays are 'automatic' (and thus
allocated by the compiler upon each call of the subroutine). In a
second version, this repeated allocation is avoided by promoting
the variables to have module scope. For each of these versions,
PSyclone is used to produce an OpenACC version, an SIR-compliant
version (with certain intrinsics and code structures replaced)
and an OpenACC+SIR-compliant version.

Note that at present the SIR produced for this version of the
benchmark does not result in correct CUDA being generated. The reasons
for this are under investigation in the ESIWACE2 project.  The
'multi-kernel' version of the benchmark described below is a workaround
that *does* permit correct CUDA code to be generated.

## multi_kernel ##

This directory contains a version of the benchmark where the compute
has been further sub-divided into ten subroutines. This form can then
be processed by PSyclone to produce SIR that Dawn is able to successfully
process allowing correct, optimised CUDA to be produced. As with the
`compute_in_subroutine` version, PSyclone is used to produce OpenACC,
SIR-compliant and OpenACC+SIR-compliant versions.

