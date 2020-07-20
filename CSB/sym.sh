#!/bin/perl

#!/bin/perl

@matrices = ("wind","bone010","GHS_psdef-ldoor","kkt_power","grid3d-200","nlpkkt200");
@rblock = ("2","4","8");
@procs = ("1","2","4","8","16");

foreach $matrices (@matrices)
{
        foreach $rblock (@rblock)
        {
                foreach $procs (@procs)
                {
                        $str = "CILK_NPROC=$procs numactl --cpunodebind 0 ./spmv_r$rblock ../BinaryMatrices/$matrices.bin nosym binary";
                        print $str;
                        system($str);

                        $str = "CILK_NPROC=$procs numactl --cpunodebind 0 ./spmv_r$rblock ../BinaryMatrices/$matrices\_rcm.bin nosym binary";
                        print $str;
                        system($str);

                        $str = "CILK_NPROC=$procs numactl --cpunodebind 0 ./sym_r$rblock ../BinaryMatrices/$matrices\_rcm.bin nosym binary";
                        print $str;
                        system($str);
                }
        }
}


