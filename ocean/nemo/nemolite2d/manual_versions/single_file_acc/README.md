# Single source-file OpenACC version of NEMOLite2D #

The source for this version is in nemolite2d.f90. Although it is primarily
parallelised using OpenACC, some experiments were also done with CUDA Fortran
and thus it must be compiled with CUDA Fortran enabled.

In addition to the 'final' version in `nemolite2d.f90`, this directory also
includes other versions that show various stages in the optimisation
process:

collapse_loops: loops have the COLLAPSE clause applied to increase the
                parallelised iteration space.
cuda_mom: uses the CUDA-Fortran version of the Momentum kernels.
fused_mom: has the two (u and v) Momentum kernels fused.
indep_clause: has the INDEPENDENT clause applied to loops around the Flather
              boundary conditions (to assure the compiler that they can
              be parallelised).
mom_short_loop: has the SHORTLOOP clause applied to the loop directives
                around the momentum kernels.
vanilla: the basic, unoptimised OpenACC version.

In order to build any of these versions, simply link or re-name the
corresponding file to `nemolite2d.f90`.

## Compiling ##

To build this version of NEMOLite you will need a Fortran compiler
with OpenACC and CUDA Fortran support, i.e. PGI.
A script is provided to set the necessary compiler flags etc. in
../../../../../compiler_setup/pgi_acc.sh. Depending on your specific system
(e.g. model of GPU) you many need to edit this file. Assuming that it
matches your requirements then doing (in Bash):

    > . ../../../../../compiler_setup/pgi_acc.sh
    > make

should compile the application and create a `nemolite2d.exe` binary.

## Running ##

Is simply a matter of doing:

    > ./nemolite2d.exe

If you want to verify that computation is actually being performed on
the GPU then you can set PGI_ACC_TIME=1 before executing the code.

Model parameters (size of domain [jpiglo,jpjglo], number of time-steps
[nitend], whether or not and how often to do output [irecord]) may be
configured by editing the `namelist` file.
