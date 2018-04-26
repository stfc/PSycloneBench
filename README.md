# PSycloneBench

Various small benchmarks used to inform the development of the
[PSyclone](https://github.com/stfc/psyclone) Domain-Specific Compiler.

## Obtaining the Code ##

The benchmark codes contained in this project use the
[dl_timer](https://bitbucket.org/apeg/dl_timer) library for execution
timing and the [dl_esm_inf](https://github.com/stfc/dl_esm_inf)
infrastructure. These packages are included within the PSycloneBench
repository as git submodules.  As such, when cloning the repository
you need to use the `--recursive` flag in order to get the source code
of that submodule, e.g.:

    git clone --recursive https://github.com/stfc/PSycloneBench

If you forget to do this then your cloned repository will contain
empty `PSycloneBench/shared/dl_timer` and `dl_esm_inf` directories. To
populate them you can do:

    cd PSycloneBench
    git submodule init
    git submodule update

## PSyKAl Version of the Code ##

The `src/psykal` directory contains a version of the benchmark code
that has been re-structured following the Parallel System Kernel
Algorithm (PSyKAl) separation of concerns. It also makes use of the
dl_esm_inf infrastructure support for the creation of grid and field
objects.

