# Manual versions of NEMOLite2D #

## single_file_acc ##

OpenACC versions of the single source-file form of NEMOLite2D with
different optimisations applied.

## psykal_serial ##

Various serial versions of the PSyKAl form of NEMOLite2D with different
optimisations applied to the PSy layer in order to test their effect
on performance (and the ability of Habakkuk
[https://github.com/arporter/habakkuk] to predict it).

## psykal_acc ##

OpenACC version of the PSyKAl form of NEMOLite2D.

## psykal_ocl ##

OpenCL versions of the PSyKAl form of NEMOLite2D. One with a C driver layer,
the second with a Fortran driver layer (using clFortran). The OpenCL kernels
for these benchmarks are in ../kernels/opencl.

## psykal_omp ##

OpenMP version of the PSyKAl form of NEMOLite2D.

