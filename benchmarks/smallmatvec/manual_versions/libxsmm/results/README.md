This directory contains results from running this benchmark on a
single Skylake node of Scafell Pike (SfP). These nodes contain 2 Skylake
CPUs of type Xeon Gold 6142, which have 16 cores each and run at
2.6GHz. The results were compiled using the Intel compiler version
"19.0.3.199 20190206". Results were run at the end of May 2020. The
source of libxsmm version 1.15 release was downloaded and compiled on
SfP using the same intel compilers as above.

Result (.txt) files
===================

Each line in the results files contain the following information 1)
Problem size, 2) number of levels, 3) the chosen OpenMP schedule, 4)
the chosen OpenMP affinity, 5) the number of cores requested, 6) the
number of SMT threads requested, 7) the time taken in seconds for the
original kernel, 8) a checksum of the results of the original kernel,
9) the time taken in seconds for the libxsmm kernel, 10) a checksum of
the results of the libxsmm kernel, 11) the time in seconds for the
kinner_reordered kernel and 12) a checksum of the results of the
kinner_reordered kernel.

In all cases KMP_AFFINITY was set to "scatter", the OMP_SCHEDULE was
set to "static" KMP_HW_SUBSET was set to choose the number of cores
required and to specify that no SMT threads were required (as initial
tests showed that they slowed the code down).

32x32 columns, 32 levels varying number of threads
==================================================

Results libxsmm_results_h32_v32_schstatic_affscatter_smt1_r1.txt and
libxsmm_results_h32_v32_schstatic_affscatter_smt1_r2.txt capture the
timing results when running the benchmark with 32x32 columns, 32
levels and varying the number of threads used from 1 to
32.

The gnuplot file "results_v32_threads.plot" will plot these results:
> gnuplot
gnuplot> load "results_v32_threads.plot"

The resultant plot is provided in "results_v32_threads.png"

32x32 columns, 100 levels varying number of threads
==================================================

This test is the same as the previous one except that the number of
levels is now 100 rather than 32. All associated files have the same
names as the previous test apart from "v32" being replaced with
"v100".

32x32 columns, 32 threads on 32 cores, varying number of levels
===============================================================

Result "libxsmm_results_h32_c32_schstatic_affscatter_smt1_r1.txt"
captures the timing results when running the benchmark with 32x32
columns, 32 threads on 32 cores (fully populated node) and varying the
number of levels to see how that affects the performance of the
different versions.

The gnuplot file "results_h32_vertical.plot" will plot these results:
> gnuplot
gnuplot> load "results_h32_vertical.plot"

The resultant plot is provided in "results_h32_vertical.png"

32 levels, 32 threads on 32 cores, varying number of cell columns
=================================================================

Result "libxsmm_results_v32_c32_schstatic_affscatter_smt1_r1.txt"
captures the timing results when running the benchmark with 32 levels,
32 threads on 32 cores (fully populated node) and varying the number
of columns to see how that affects the performance of the different
versions.

The gnuplot file "results_v32_horizontal.plot" will plot these results:
> gnuplot
gnuplot> load "results_v32_horizontal.plot"

The resultant plot is provided in "results_v32_horizontal.png"

100 levels, 32 threads on 32 cores, varying number of cell columns
=================================================================

This test is the same as the previous test except that the number of
levels is 100 rather than 32. All the resultant files have the same
names as the previous test apart from "v32" being replaced with
"v100".
