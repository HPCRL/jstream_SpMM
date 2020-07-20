

jspmm=bin/j_spmm_model.exe
jsddmm=bin/j_sddmm_model.exe
tspmm=bin/taco_spmm.exe
tsddmm=bin/taco_sddmm.exe
mkl=bin/mkl.exe
csb=CSB/spmm_d

folder=$1

cache=$(cat /sys/devices/system/cpu/cpu0/cache/index2/size)

#echo $cache

cache=$(echo $cache| cut -d'K' -f 1)

#echo $cache

cache=$(expr $cache \* 128)

echo Cache size for L2 is $cache doubles 

for i in $folder/*mtx
do
	for K in 128 1024
	do
		./$jspmm $i cache $cache Nk $K
		./$jsddmm $i cache $cache Nk $K
		./$tspmm $i $K
		./$tsddmm $i $K
		./$csb$K $i
		./$mkl $i $K
	done
echo Result collection for $i finished
done

