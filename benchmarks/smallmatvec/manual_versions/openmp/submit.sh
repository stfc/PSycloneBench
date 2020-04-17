#!/bin/bash

#PBS -l select=1:ncpus=68:mpiprocs=68:mem=86GB:mcdram=cache:numa=quadrant
#PBS -l walltime=10:00:00 
#PBS -A cin_external_0

numactl -H
cd /marconi/home/userexternal/ssiso000/gungho-mv
module load intel
./scalability.sh

