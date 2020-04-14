# LFRic matrix-vector multiplication kernel optimizations

This project contains multiple implementations of the LFRic Matrix-vector multiplication kernel.

Relevant files are:

- `kdriver.F90` multiple thread-level parallelism and iteration order implementations.
- `matrix_vector_kernel_mod.F90` multiple implementation of the mv kernel.
- `reports/` Documents the mv kernel performance as well as full LFRic scalability analysis made in the IPCC project.

Additionally, the Full LFRic MPI performance can be visualized [here](https://sergisiso.gitlab.io/multiperf/)

