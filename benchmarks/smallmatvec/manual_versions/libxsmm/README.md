
This directory allows the performance of the default matrix vector
kernel benchmark to be compared with a libxsmm version. In the libxsmm
version the matmul intrinsic is replaced and an explicit just-in-time
compilation call is placed outside of the k-loop.

This code was originally proposed, implemented and tested (on the full
LFRic code) by Wolfgang Hayek where he found a significant performance
improvement. This version differs slightly from his version as he also
included a special case when the number of dofs in both cases is 1. In
this benchmark we do not have this situation so can safely ignore this
additional optimisation. Optimising for particular numbers of dofs is
being done but is being treated as a separate optimisation.

Files are

matrix_vector_kernel_mod.f90: the original matrix vector kernel with
metadata removed to reduce the number of modules needed to compile the
example.

matrix_vector_kernel_xsmm_mod.f90: the libxsmm modified matrix vector
    kernel which also has the metadata removed.

kdriver.f90: a standalone driver for the benchmark which includes the
    coloured OpenMP parallelisation in the horizontal. This version
    generates its own data so can be used for different problem
    sizes. These sizes are set on the command line like this
    ... kdriver_version -g <nhoriz> <nvert>. This is a cut down
    version of code written by Sergi Siso.

constants_mod.f90: precision constants needed by the kernel and
    driver.

Makefile: a Makefile designed to use the intel compiler. The makefile
    relies on libxsmm being installed and pointed to appropriately in
    order to compile the codes (the original and the xsmm version).

runme.sh: a sample script to run the benchmark for different
    resolutions, different numbers of processors etc. This should be
    tailored to the particular cpu the code is being run on and the
    problem sizes being examined.

