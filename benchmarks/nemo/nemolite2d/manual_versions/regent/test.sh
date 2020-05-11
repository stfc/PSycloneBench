

for i in 24 28 32
do
  for j in 4 8 
  do
  for t1 in 84 #256
  do
  for t2 in 128 #512 #256 512
  do
    #legion/language/regent.py algorithm.rg -ll:cpu $i -ll:util $j -ll:ocpu 2 -ll:othr 6 -ll:ht_sharing 0  > results/results_${i}_${j}.txt
    #legion/language/regent.py algorithm.rg -ll:cpu $i -ll:util $j -ll:ocpu 2 -ll:othr 6 > results/results_${i}_${j}.txt
    #echo "legion/language/regent.py algorithm.rg -ll:cpu $i -ll:util $j -ll:csize 1024 -t1 ${t1} -t2 ${t2}"
    legion/language/regent.py algorithm.rg -ll:cpu $i -ll:util $j -ll:csize 1024 -t1 ${t1} -t2 ${t2} > results/results_${i}_${j}_${t1}_${t2}.txt
    b=`grep "Runtime" results/results_${i}_${j}_${t1}_${t2}.txt | awk '{print $3}'`
    echo $i $j ${t1} ${t2} $b >> results/res.txt
  done
  done
  done
done
