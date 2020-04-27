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

## Benchmarks ##

Are contained in directories beneath the `benchmarks` directory.

### NEMOLite2D ###

`benchmarks/nemo/nemolite2d` contains various versions of the NEMOLite2D
benchmark which is based upon the free-surface component of the NEMO
ocean model (www.nemo-ocean.eu). This was originally constructed
by Hedong Liu of the National Oceanography Centre, Liverpool, UK.

#### Manual implementations ####

The `manual_versions` directory contains various implementations of
NEMOLite2D where the application has been split-up into Algorithm PSy
and Kernel layers but the PSy layer has been written manually.
See the [README.md](ocean/nemo/nemolite2d/manual_versions/README.md)
file in the NemoLite2D manual`_versions directory for more information.

(Note that the OpenACC versions are based upon work by Jeremy
Appleyard of NVIDIA).

#### PSyclone code generation ####

The `psykal` directory contains the version of NEMOLite2D that uses
PSyclone to generate the PSy layer. Currently serial and OpenMP
versions may be generated.

### Shallow ###

`benchmarks/shallow` contains various versions of the Shallow benchmark,
originally developed by Paul Swarztrauber of NCAR. In summary these are:

* SEQ/original - original single-file serial version in Fortran.
* SEQ - PSyKAl version of Shallow including option to build with PSyclone.
* OMP - manual OpenMP implementation in C.

Please see the README.md files in the individual directories for more
information.

## Building ##

The Makefiles for all of the various benchmarks pick-up which compiler
to use (and associated flags) from environment variables. Bash scripts to
set these variables for a variety of compilers are in the compiler_setup
directory. When benchmarking you will need to edit the relevant file
and ensure that the compiler flags are set appropriately. The flags that
give the best performance for a particular compiler can vary from benchmark
to benchmark and even by problem size.

In order to use the settings for a particular compiler, you must source
the appropriate Bash script. e.g. for the Intel compiler (assuming you are
running a Bash shell):

    > . compiler_setup/intel.sh

Again, see the README.md file included with each benchmark.

A top-level Makefile is provided that has targets to build all of the
various NEMOLite2D and Shallow benchmarks for CPU:

 * all - build all targets that do not require PSyclone
 * all_gen - build all targets that use PSyclone code generation
 * shallow_cpu - build all CPU versions of Shallow that do not require
                 PSyclone
 * shallow_gen - build Shallow using PSyclone for code generation
 * nemolite_cpu - build all CPU versions of NEMOLite2D that do not require
                  PSyclone
 * nemolite_gen - build NEMOLite2D using PSyclone for code generation

Since the set-up for OpenACC and OpenCL is more involved, it does not
attempt to compile the benchmarks using those languages. See the Makefiles
in the individual benchmark folders for those cases.

### LFRic Small Matrix Vector Multiplication ###

`benchmarks/smallmatvec` contains multiple implementations of the
matrix-vector multiplication operation done in the LFRic application. 

Currently it only contains one manual implementation:

 * openmp - Configurable implementation using OpenMP that test the
            vectorization and threading capabilities of multiple versions
            of the kernel and multiple parallel traversing strategies.

Note that this benchmark does not use the `compiler_setup` settings yet.
