# Configurable OpenMP manual implementation


This implementation uses OpenMP and test the vectorization and threading
capabilities of multiple versions of the kernel and multiple parallel
traversing strategies.

It was initially created by Chris Maynard using the PSyKE and Dino tools to
extract the needed infrastructure and data from LFRic. Then extended by
Sergi Siso to make it configurable to multiple OpenMP implementations and
benchmark the KNL architecture for the IPCC at Hartree Centre project.

Th most relevant files are:

- `kdriver.F90` multiple thread-level parallelism and iteration order
implementations.
- `matrix_vector_kernel_mod.F90` multiple implementation of the mv kernel.
- `reports/` Documents the mv kernel performance as well as full LFRic
scalability analysis made in the IPCC project.
