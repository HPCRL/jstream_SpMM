#!/bin/perl

#@matrices = ("largebasis","asic_320k","cont11","rucci","wiki2007","vanHeukelum-cage15","tsopf1","tsopf2","rajat30");

@matrices = ("tsopf2","rajat30");
@rblock = ("2","4","8");

foreach $matrices (@matrices)
{
	foreach $rblock (@rblock)
	{
		$str = "CILK_NPROC=1 numactl --cpunodebind 0 ./spmv_r$rblock ../BinaryMatrices/$matrices.bin nosym binary";
		print $str;
		system($str);

		$str = "CILK_NPROC=2 numactl --cpunodebind 0 ./spmv_r$rblock ../BinaryMatrices/$matrices.bin nosym binary";
		print $str;
		system($str);

		$str = "CILK_NPROC=4 numactl --cpunodebind 0 ./spmv_r$rblock ../BinaryMatrices/$matrices.bin nosym binary";
		print $str;
		system($str);

		$str = "CILK_NPROC=8 numactl --cpunodebind 0 ./spmv_r$rblock ../BinaryMatrices/$matrices.bin nosym binary";
		print $str;
		system($str);

		$str = "CILK_NPROC=16 numactl --cpunodebind 0 ./spmv_r$rblock ../BinaryMatrices/$matrices.bin nosym binary";
		print $str;
		system($str);
	}
}
