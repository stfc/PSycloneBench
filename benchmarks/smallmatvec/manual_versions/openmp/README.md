# Configurable OpenMP manual implementation


This implementation uses OpenMP and test the vectorization and threading
capabilities of multiple versions of the kernel and multiple parallel
traversing strategies.

Th most relevant files are:

- `kdriver.F90` multiple thread-level parallelism and iteration order
implementations.
- `matrix_vector_kernel_mod.F90` multiple implementation of the mv kernel.
- `reports/` Documents the mv kernel performance as well as full LFRic
scalability analysis made in the IPCC project.

# Limitations

Currently compiler and compiler flags are hardcoded into the Makefile, the
PSycloneBench compiler environment variables are ignored.

# Attributions 

It was initially created by Chris Maynard using the PSyKE and Dino tools to
extract the needed infrastructure and data from LFRic. Then extended by
Sergi Siso to make it configurable to multiple OpenMP implementations and
benchmark the KNL architecture for the IPCC at Hartree Centre project.

Some files from the LFRic infrastructure are needed to run the benchmark,
these are `argument_mod.F90`, `constants_mod.F90`, `dino_mod.F90` and
`kernel_mod.F90`. These files have been copied from
https://github.com/christophermaynard/LFRic-microbenchmarks/common
and are distributed with the attached licence file LICENNCE.MetOffice.
