# Manual OpenCL version of PSyKAl NEMOLite2D #

This directory contains a version of the PSyKAl form of NEMOLite2D
that has had its PSy layer (in `time_step_mod.f90`) manually
parallelised using OpenCL.

## Compiling ##

The Makefile picks up the compiler and associated flags from environment
variables. See e.g. ../../../../../compiler_setup/gnu.sh for sample
settings for the Gnu compiler suite (with JIT OpenCL).
If you are using Bash and the Gnu compiler then:

    > source ../../../../../compiler_setup/gnu.sh
    > make

should build the `nemolite2d.exe` binary.

Note that in some platforms, like FPGAs, using a JIT solution for the
OpenCL kernels is not feasible. In this case the OpenCL kernel binary
has to be compiled ahead of time with the ``device_binary`` Makefile
target and the proper environment loaded. For example:

    > source ../../../../../compiler_setup/xilinx.sh
    > make device_binary

Since this process can take a long time it is convenient to launch
the building process detached from the terminal session with:

    > make device_binary_nohup

There are equivalent ``device_binary_tasks`` and ``device_binary_tasks_nohup``
to compile a task-based version of the OpenCL kernels. Currently this is
constrained to a single problem size set by a define pre-processor parameter
in the ``allkernel_tasks.cl`` (set by default to 250x250 problem sizes).

## Running ##

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.

Additionally the OpenCL execution requires that the environment variable
FORTCL_KERNELS_FILE is set up to reference the target device binary
or the target source code. If the source code is given, the OpenCL
runtime system will JIT compile the necessary device objects at runtime.
The following command shows an example of a JIT compiled execution:

    > FORTCL_KERNELS_FILE=allkernels.cl ./nemolite2d.exe 

Also, it is often necessary to specify the desired target platform and the
padding necessary (mandatory for OpenCL < 2) to launch the OpenCL application
as expected. These can be set respectively by the environment variables
``FORTCL_PLATFORM`` and `` DL_ESM_ALIGNMENT``. A complete example of an OpenCL
execution command could look like this:

    > DL_ESM_ALIGNMENT=64 FORTCL_PLATFORM=2 \
      FORTCL_KERNELS_FILE=allkernels.cl ./nemolite2d.exe

Note that if the ``task_optimization`` parameter has been turned true on the
``time_step_mod.f90``, the OpenCL kernels provided must also contain a task
implementation. For example:

    > FORTCL_KERNELS_FILE=allkernels_tasks.cl ./nemolite2d.exe 

## Output ##

If output is enabled (`irecord` < `nitend` in the `namelist` file) then
field values are written to the ASCII file `go2d_<time-step>.dat`. This
is formatted for use with gnuplot's `splot` command with data in columns:

x-ordinate, y-ordinate, depth, sea-surface height, u, v

where u and v are the x and y components of velocity, respectively.


