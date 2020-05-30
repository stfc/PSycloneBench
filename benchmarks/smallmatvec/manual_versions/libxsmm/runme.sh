#!/usr/bin/env bash

# Note, this script has Intel-specific KMP environment variables in it

# This script executes the selected original and xsmm matrix_vector
# multiplication implementation with a variety of OpenMP
# configurations to analyse its performance.

# Update the for statements to configure each parameter search space.

# The configurable parameters are:
#   - $ps: problem size (horizontal)
#   - $vertical: problem size (vertical)
#   - $omp_sch: OpenMP scheduler [static, dynamic, guided]
#   - $omp_aff: OpenMP affinity [none, compact, scatter]
#   - $smt: Number of SMT
#   - $cores: Number of cores
#   - $v: matrix_vector implementation [orig, xsmm]

# echo "Don't forget to set the parameters to the required values"
# exit

echo "Original vs libxsmm timing results"

for ps in 32; do
    for vertical in 32 100; do
        for omp_sch in static ; do
            for omp_aff in scatter ;  do
		for smt in 1; do
		    output=libxsmm_results_h${ps}_v${vertical}_sch${omp_sch}_aff${omp_aff}_smt${smt}.txt
                    for cores in 1 2 4 8 12 16 24 32; do
			echo $v $ps $vertical $omp_sch $omp_aff $cores $smt
			echo -n $ps $vertical $omp_sch $omp_aff $cores $smt >> $output
			for v in orig xsmm ; do
			    OMP_NUM_THREADS=$(( ${cores} * ${smt} )) \
			      KMP_HW_SUBSET=${cores}c,${smt}t \
			      KMP_AFFINITY=$omp_aff OMP_SCHEDULE=$omp_sch \
			      ./kdriver_$v -g $ps $vertical > output
			    t=$(cat output | grep "Loop time" | awk '{print $5}')
			    check=$(cat output | grep "Reduction value" | awk '{print $3}')
			    echo -n " " $t $check >> $output
			done
			echo >> $output
		    done
		done
            done
        done
    done
done
exit
