#!/bin/bash
data=$1
resultdir="results"
result="$resultdir/results.txt"
scr=scripts
jstream=$scr/jstreamex.py
taco=$scr/tacoex.py
mkl=$scr/mklex.py
csb=$scr/csbex.py
collect=$scr/collect.sh
rcsv="$resultdir/results.csv"
summarize=$scr/summarize.py
summary="summary.csv"

#make -j8
#cd CSB
#make spmm_dall
#cd ..

if [ ! -d $resultdir ]; then
  mkdir -p $resultdir;
fi


bash $collect $data >>  $result


for script in  $jstream $taco $mkl $csb
do
	python3 $script $result >> $rcsv
done

echo "Summarizing..."
python3 $summarize $rcsv $summary

