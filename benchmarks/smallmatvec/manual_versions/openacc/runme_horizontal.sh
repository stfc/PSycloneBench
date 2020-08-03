#!/usr/bin/env bash

# This script executes the selected xxxxx matrix
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
#   - $v: matrix_vector implementation [orig, xsmm, reorder]

# echo "Don't forget to set the parameters to the required values"
# exit

echo "v100 openacc timing results"

for vertical in 32 100; do
    #output=v100_reorder_horizontal_v${vertical}.txt
    output=v100_nvidia_horizontal_v${vertical}.txt
    for ps in 8 16 32 64 96 128 ; do
	echo $vertical $ps
	echo -n $vertical $ps >> $output
	#for version in reorder1 reorder2; do
	for version in nvidia4; do
            ./kdriver_openacc_${version} ${ps} ${vertical}> output
	    t=$(cat output | grep "Loop time" | awk '{print $5}')
	    check=$(cat output | grep "Reduction value" | awk '{print $3}')
	    echo -n " " $t $check >> $output
	done
	echo >> $output
    done
done
exit
