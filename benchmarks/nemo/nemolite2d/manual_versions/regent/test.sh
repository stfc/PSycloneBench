#!/bin/bash
#This script will run the regent implementation of NEMOLite2D with 
#different numbers of threads, tilesizes etc. to enable easy benchmarking.
#The individual results will be output to files, and results/res.txt will
#contain a markdown compatible table of results. 
echo "| cpu | util | runtime (s) |" >> results/res.txt
echo "|-----|------|-------------|" >> results/res.txt

for i in 24 28 #32
do
  for j in 4 8 
  do
  for t1 in 256
  do
  for t2 in 1024 #512 #256 512
  do
  for k in 0 1 2 
  do
    legion/language/regent.py algorithm.rg -ll:cpu $i -ll:util $j -ll:csize 8000 -t1 ${t1} -t2 ${t2} > results/results_${i}_${j}_${t1}_${t2}_${k}.txt
    b=`grep "Runtime" results/results_${i}_${j}_${t1}_${t2}_${k}.txt | awk '{print $3}'`
    echo "| $i | $j | $b |" >> results/res.txt
  done
  done
  done
  done
done
