#!/usr/bin/env bash

output=scalability.txt

echo "Scalability plot" | tee $output

#make && make kdriver.novec

#for v in `make versions`; do
for v in nlayersf ; do
    for ps in 32; do
        for omp_s in omp-locking colouring colouring2 colouring-rows ; do # ompall
            for omp_sch in static dynamic guided ; do
		for omp_aff in none compact scatter ;  do
		    for smt in 1 2 4; do
			for cores in 1 4 8 `seq 16 16 64`; do
			    echo $v $omp_s $smt $cores
			    KMP_HW_SUBSET=${cores}c,${smt}t KMP_AFFINITY=$omp_aff OMP_SCHEDULE=$omp_sch ./kdriver.avx512.$v -g $ps 16 -t $omp_s > output
			    t=$(cat output | grep "Loop time" | awk '{print $5}')
			    check=$(cat output | grep "Reduction value" | awk '{print $3}')
			    #t=$(OMP_NUM_THREADS=${cores} OMP_SCHEDULE=$omp_sch ./kdriver.avx512.$v -g $ps 16 -t $omp_s | grep "Loop time" | awk '{print $5}')
			    echo $v $ps $omp_s $omp_sch $omp_aff $cores $smt $t $check >> $output
		        done
		    done
		done
	    done
        done
    done
done
