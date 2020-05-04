# Manual versions of NEMOLite2D #

The manual implementation versions of NemoLite2D are divided in two subgroups:

1. The "PSyKAl" implementations, which have each application separated into
three layers (Algorithm, Parallel System and Kernel) of the PSyKAl model but
with the Parallel System layer implemented manually.
In these implementations usually the Algorithm can be found in the
`nemolite2d.f90` file, the manual PSy-layer is found in `time_step.f90` and
the kernels which are encoded as independent functions in `time_step.f90` are
found in the same directory or in the `../kernels` directory.
These are:

    - `psykal_serial`: Various serial versions of the PSyKAl form of
    NEMOLite2D with different optimisations applied to the PSy layer in order
    to test their effect on performance (and the ability of Habakkuk
    [https://github.com/arporter/habakkuk] to predict it).

    - `psykal_omp`: OpenMP version of the PSyKAl form of NEMOLite2D.

    - `psykal_dm`: Distributed Memory (MPI) version of the PSyKAl form of
    NEMOLite2D.

    - `psykal_acc`: OpenACC version of the PSyKAl form of NEMOLite2D.

    - `psykal_opencl`: OpenCL version of the PSyKAl form of NEMOLite2D.

    - `psykal_ompss`: OmpSs version of the PSyKAl form of NEMOLite2D.


2. The "single\_file" versions that due to either programming model limitations
or code simplification/optimisation do not follow the PSyKAl separation of
concerns. These are:

    - `single_file_acc`: OpenACC version of the single source-file form of
    NEMOLite2D with different optimisations applied.

    - `single_file_opencl_fortran`: OpenCL version of the single-file form of
    NEMOLite2D. This is a Fortran (using clFortran) data-parallel
    implementation of the non-PSyKAl form.

    - `single_file_opencl_c`: OpenCL version of the single-file form of
    NEMOLite2D. This is a C data-parallel implementation of the non-PSyKAl
    form.

## Makefile ##
A top-level Makefile is provided to recurse down to the multiple
implementations, compile them and execute them using a common set of
parameters (provided by the top-level namelist).
