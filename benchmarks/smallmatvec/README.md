# LFRic Small Matrix Vector Multiplications Benchmark

This project contains multiple implementations of the matrix-vector
multiplication operation done in the LFRic application. 

## Manual Versions

Manual implementations of the kernel are found in teh `manual_versions`
directory. At the moment the following versions are available:

    - `openmp`: Configurable implementation using OpenMP that test
    the vectorization and threading capabilities of multiple versions
    of the kernel and multiple parallel traversing strategies.
