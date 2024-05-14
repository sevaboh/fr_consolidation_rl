#!/bin/sh
ulimit -c unlimited
export OMP_THREAD_LIMIT=8
#g++ -O3 -fopenmp fr_consolidation_rl.cpp libfftpack/libfftpack.a -o fr_cons_rl
./fr_cons_rl debug_level 0 beta 2.0 ls_max_iter 100 testing 2 ls_max_tau 0.2 NB 25; mv out_H.txt out_H2_25.txt
./fr_cons_rl debug_level 0 beta 2.0 ls_max_iter 100 testing 2 ls_max_tau 0.1 NB 50; mv out_H.txt out_H2_50.txt
./fr_cons_rl debug_level 0 beta 2.0 ls_max_iter 100 testing 2 ls_max_tau 0.05 NB 100; mv out_H.txt out_H2_100.txt
./fr_cons_rl debug_level 0 beta 2.0 ls_max_iter 100 testing 2 ls_max_tau 0.025 NB 200; mv out_H.txt out_H2_200.txt
./fr_cons_rl debug_level 0 beta 2.0 ls_max_iter 100 testing 2 ls_max_tau 0.0125 NB 400; mv out_H.txt out_H2_400.txt
./fr_cons_rl debug_level 0 beta 2.0 ls_max_iter 100 testing 2 ls_max_tau 0.00625 NB 800; mv out_H.txt out_H2_800.txt
./fr_cons_rl debug_level 0 beta 2.0 ls_max_iter 100 testing 2 ls_max_tau 0.003125 NB 1600; mv out_H.txt out_H2_1600.txt

./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.2 NB 25; mv out_H.txt out_H1_25.txt
./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.1 NB 50; mv out_H.txt out_H1_50.txt
./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.05 NB 100; mv out_H.txt out_H1_100.txt
./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.025 NB 200; mv out_H.txt out_H1_200.txt
./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.0125 NB 400; mv out_H.txt out_H1_400.txt
./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.00625 NB 800; mv out_H.txt out_H1_800.txt
./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.003125 NB 1600; mv out_H.txt out_H1_1600.txt

# test constructed solution
time -o times.txt ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.1 NB 50; mv out_H.txt out_H_1_50.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.01 NB 50; mv out_H.txt out_H_01_50.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.001 NB 50; mv out_H.txt out_H_001_50.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.0001 NB 50; mv out_H.txt out_H_0001_50.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.1 NB 200; mv out_H.txt out_H_1_200.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.01 NB 200; mv out_H.txt out_H_01_200.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.001 NB 200; mv out_H.txt out_H_001_200.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.0001 NB 200; mv out_H.txt out_H_0001_200.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.1 NB 400; mv out_H.txt out_H_1_400.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.01 NB 400; mv out_H.txt out_H_01_400.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.001 NB 400; mv out_H.txt out_H_001_400.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.0001 NB 400; mv out_H.txt out_H_0001_400.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.1 NB 800; mv out_H.txt out_H_1_800.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.01 NB 800; mv out_H.txt out_H_01_800.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.001 NB 800; mv out_H.txt out_H_001_800.txt
time -o times.txt -a ./fr_cons_rl debug_level 0 testing 1 ls_max_tau 0.0001 NB 800; mv out_H.txt out_H_0001_800.txt
# test analytic solution
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.1 NB 50; mv out_H.txt out_H2_1_50.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.01 NB 50; mv out_H.txt out_H2_01_50.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.001 NB 50; mv out_H.txt out_H2_001_50.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.0001 NB 50; mv out_H.txt out_H2_0001_50.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.1 NB 200; mv out_H.txt out_H2_1_200.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.01 NB 200; mv out_H.txt out_H2_01_200.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.001 NB 200; mv out_H.txt out_H2_001_200.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.0001 NB 200; mv out_H.txt out_H2_0001_200.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.1 NB 400; mv out_H.txt out_H2_1_400.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.01 NB 400; mv out_H.txt out_H2_01_400.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.001 NB 400; mv out_H.txt out_H2_001_400.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.0001 NB 400; mv out_H.txt out_H2_0001_400.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.1 NB 800; mv out_H.txt out_H2_1_800.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.01 NB 800; mv out_H.txt out_H2_01_800.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.001 NB 800; mv out_H.txt out_H2_001_800.txt
time -o times.txt -a ./fr_cons_rl beta 2.0 debug_level 0 ls_max_iter 100 testing 2 ls_max_tau 0.0001 NB 800; mv out_H.txt out_H2_0001_800.txt
