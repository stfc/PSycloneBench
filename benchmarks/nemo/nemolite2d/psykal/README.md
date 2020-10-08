# PSyKAl NEMOLite2D with Code Generation #

This directory contains the PSyKAl version of NEMOLite2D that uses
PSyclone to generate (and transform) the Parallel System (PSy) layer.

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite. In Bash you can source that file
to set things up for Gnu:

    > . ../../../../compiler_setup/gnu.sh

However, if the build target requires the MPI library, the `F90` environment
variable should be changed to the desired MPI compiler wrappers, for example:

    > export F90=mpif90

Note that the MPI compiler wrapper should call the same compiler family as that
loaded by the compiler setup script, otherwise the build system compiler flags
will be unrecognised. It is also advisable to perform a `make allclean` action
when requesting a new `F90` compiler in order to compile all dependencies with
the same parameters.

Since the Makefile uses PSyclone you must have PSyclone installed on your
system (see http://psyclone.readthedocs.io/en/stable/getting_going.html).

The Makefile supports many targets:

- nemolite2d_serial (default) - vanilla (serial) version
- nemolite2d_omp - OpenMP version with static scheduling
- nemolite2d_mpi - MPI version for distributed memory parallelism
- nemolite2d_hybrid - MPI and OpenMP combined version
- nemolite2d_acc - OpenACC (Not working - PSyclone will not generate it)
- nemolite2d_ocl - OpenCL version
- nemolite2_mpiocl - MPI and OpenCL combined version (Not working - Psyclone
will generate a OpenCL-only version)

Executing `make <target>` will generate a folder with the same name as the
target, containing at least a `alg.f90` and a `psy.f90` for the
PSyclone-generated Algorithm and PSy-layer components respectively.
It will also generate kernel source files when necessary. Finally, it will
build a `nemolite2d.exe` with the source files in that target folder.


## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

There are differences in how to execute the generated `nemolite2d.exe`,
notably:

- For the OpenMP version, you will need to set `OMP_NUM_THREADS` to the
desired number of threads to use. Note that if you are benchmarking
performance then you probably need to ensure that the threads are bound
sensibly to processor cores (see e.g. taskset, likwid or Intel's KMP_AFFINITY).

- For the MPI version, you will need to prefix the binary with the
`mpirun -n <#ranks>` command. See additional MPI parameters for the selected
implementation with `mpirun -h`.

- For the OpenCL version, you will need to specify which kernels file to use
with the `FORTCL_KERNELS_FILE` environment variable, it can be a source file
for JITed OpenCL or a ahead-of-time compiled binary object. Optionally you can
select which OpenCL platform to use with the `FORTCL_PLATFORM` environment
variable. For more information about the FortCL environment variables read
the [FortCL README](https://github.com/stfc/FortCL).
The target folder makefile also contains some execution helpers like:
`run_opencl_jit`, `run_xilinx_sw_emu` and `run_xilinx_sw_emu` to aid in the
execution of the OpenCL binary.

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.

