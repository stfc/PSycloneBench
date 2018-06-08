# PSycloneBench Utility Scripts #

This directory contains some simple utility scripts.

## build_all_cpu.sh ##

Uses the Gnu compiler suite to build all of the benchmarks that
only target a CPU (i.e. excluding those that need OpenACC, OpenCL or
Maxeler). If `psyclone` is found on the user's PATH then those
benchmarks requiring code generation are also built.

This script must be run from the `PSycloneBench` directory:

    > scripts/build_all_cpu.sh
