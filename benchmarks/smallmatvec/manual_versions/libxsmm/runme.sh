#!/usr/bin/env bash

# This script has Intel-specific KMP environment variables in it

# This script executes the selected original and xsmm matrix_vector
# multiplication implementation with a variety of OpenMP
# configurations to analyse its performance.

# The configurable parameters are:
#   - $v: matrix_vector implementation
#   - $ps: problem size (horizontal)
#   - $omp_sch: OpenMP scheduler [static, dynamic, guided]
#   - $omp_aff: OpenMP affinity [none, compact, scatter]
#   - $smt: Number of SMT
#   - $cores: Number of cores

# Update the for statements to configure each parameter search space.

echo "Don't forget to set the parameters to the required values"
exit

output=results.txt

echo "Original vs libxsmm plot" | tee $output

for v in orig xsmm ; do
    for ps in 32; do
	for vertical in 32 100; do
            for omp_sch in static ; do
            for omp_aff in scatter ;  do
            for smt in 1 2; do
            for cores in 1 2 4 6 8; do
                echo $v $smt $cores
                OMP_NUM_THREADS=$(( ${cores} * ${smt} )) \
                    KMP_HW_SUBSET=${cores}c,${smt}t \
                    KMP_AFFINITY=$omp_aff OMP_SCHEDULE=$omp_sch \
                    ./kdriver_$v -g $ps $vertical > output
                t=$(cat output | grep "Loop time" | awk '{print $5}')
                check=$(cat output | grep "Reduction value" | awk '{print $3}')
                echo $v $ps $vertical $omp_sch $omp_aff $cores $smt $t $check >> $output
            done
            done
            done
            done
        done
    done
done
