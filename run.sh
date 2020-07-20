#!/bin/bash
data=$1
result="results/results.txt"
scr=scripts
jstream=$scr/jstreamex.py
taco=$scr/tacoex.py
mkl=$scr/mklex.py
csb=$scr/csbex.py
collect=$scr/collect.sh
rcsv="results/results.csv"
summarize=$scr/summarize.py
summary="summary.csv"

#make -j8
#cd CSB
#make spmm_dall
#cd ..

bash $collect $data >>  $result


for script in  $jstream $taco $mkl $csb
do
	python3 $script $result >> $rcsv
done


python3 $summarize $rcsv $summary

