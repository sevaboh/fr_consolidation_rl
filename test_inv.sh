#!/bin/sh
ulimit -c unlimited
#export OMP_THREAD_LIMIT=1
g++ -O3 -fopenmp fr_consolidation_rl.cpp libfftpack/libfftpack.a -o fr_cons_rl
#exit
echo > inv_log.txt
for noise in 0 0.05 0.1
do
for input in 1 2 3 4 5 6 7 8 9 10 11
do
for div in 2 5 10
do
for run in 0 1 2 3 4 5 6 7 8 9
do
echo $noise $input $div $run
time ./fr_cons_rl Tm 5.0 debug_level 0 inverse out_H_test${input}.txt ls_max_tau 0.025 NB 200 pso_rp 0.01 pso_n 30 pso_max_iter 25 noise $noise pso_o 0.6 pso_fi_p 0.3 pso_fi_g 0.3 div_t $div div_x $div | tee log_inv_${input}_${div}_${noise}_${run}.txt
mv out_H.txt out_H_inv_${input}_${div}_${noise}_${run}.txt
./inv_compare.py out_H_inv_${input}_${div}_${noise}_${run}.txt out_H_testL${input}.txt >> log_inv_${input}_${div}_${noise}_${run}.txt
tail -n 5 log_inv_${input}_${div}_${noise}_${run}.txt | tr '\n' ' ' >> inv_log.txt
echo "" >> inv_log.txt
exit
done
done
done
done