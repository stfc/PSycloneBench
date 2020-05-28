#!/usr/bin/env bash

# This script executes the selected matrix_vector multiplication implementation
# with a variety of OpenMP configurations to analyse the threading scalability.

# The configurable parameters are:
#   - $v: matrix_vector implementation
#   - $ps: problem size (horizontal)
#   - $omp_s: OpenMP implementation
#   - $omp_sch: OpenMP scheduler
#   - $omp_aff: OpenMP affinity
#   - $smt: Number of SMT
#   - $cores: Number of cores
# Update the for statements to configure each parameter search space.

output=scalability.txt

echo "Scalability plot" | tee $output

for v in orig xsmm ; do
    for ps in 32; do
        for omp_s in colouring ; do
            for omp_sch in static ; do
            for omp_aff in scatter ;  do
            for smt in 1 ; do
            for cores in 1 2 4 6; do
                echo $v $omp_s $smt $cores
                OMP_NUM_THREADS=$(( ${cores} * ${smt} )) \
                    KMP_HW_SUBSET=${cores}c,${smt}t \
                    KMP_AFFINITY=$omp_aff OMP_SCHEDULE=$omp_sch \
                    ./kdriver_$v -g $ps 16 > output
                t=$(cat output | grep "Loop time" | awk '{print $5}')
                check=$(cat output | grep "Reduction value" | awk '{print $3}')
                echo $v $ps $omp_s $omp_sch $omp_aff $cores $smt $t $check >> $output
            done
            done
            done
            done
        done
    done
done
