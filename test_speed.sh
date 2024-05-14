#!/bin/sh
ulimit -c unlimited
# test openmp
g++ -O3 -fopenmp fr_consolidation_rl.cpp libfftpack/libfftpack.a -o fr_cons_rl
g++ -O3 -fopenmp fr_consolidation_rl.cpp libfftpack_orig/libfftpack.a -o fr_cons_rl_orig
echo "" > times.txt
for s in 400 600 800 1000 1200 1400 1600 2000 2400 2800 3200 4000 4800 5600 6200 6600 7000
do
echo $s >> times.txt
export OMP_THREAD_LIMIT=1
time -o times.txt -a ./fr_cons_rl_orig Tm 0.1 debug_level 0 testing 1 ls_max_tau 0.01 NB $s
time -o times.txt -a ./fr_cons_rl Tm 0.1 debug_level 0 testing 1 ls_max_tau 0.01 NB $s
export OMP_THREAD_LIMIT=2
time -o times.txt -a ./fr_cons_rl Tm 0.1 debug_level 0 testing 1 ls_max_tau 0.01 NB $s
export OMP_THREAD_LIMIT=4
time -o times.txt -a ./fr_cons_rl Tm 0.1 debug_level 0 testing 1 ls_max_tau 0.01 NB $s
export OMP_THREAD_LIMIT=6
time -o times.txt -a ./fr_cons_rl Tm 0.1 debug_level 0 testing 1 ls_max_tau 0.01 NB $s
export OMP_THREAD_LIMIT=8
time -o times.txt -a ./fr_cons_rl Tm 0.1 debug_level 0 testing 1 ls_max_tau 0.01 NB $s
done