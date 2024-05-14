#!/bin/sh
ulimit -c unlimited
#g++ -O3 -fopenmp fr_consolidation_rl.cpp libfftpack/libfftpack.a -o fr_cons_rl

./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test1.txt

./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.85 alpha 0.7 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test2.txt
./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 2.0 alpha 0.7 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test3.txt
./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.55 alpha 0.7 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test4.txt

./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.7 alpha 0.85 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test5.txt
./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.7 alpha 1.0 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test6.txt
./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.7 alpha 0.55 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test7.txt

./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.25 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test8.txt
./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.75 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test9.txt

./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.5 Cv 0.005 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test10.txt
./fr_cons_rl debug_level 0 Tm 1.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.5 Cv 0.08 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_test11.txt


./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL1.txt

./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.85 alpha 0.7 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL2.txt
./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 2.0 alpha 0.7 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL3.txt
./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.55 alpha 0.7 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL4.txt

./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.7 alpha 0.85 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL5.txt
./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.7 alpha 1.0 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL6.txt
./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.7 alpha 0.55 zeta 0.5 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL7.txt

./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.25 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL8.txt
./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.75 Cv 0.02 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL9.txt

./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.5 Cv 0.005 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL10.txt
./fr_cons_rl debug_level 0 Tm 5.0 Om 0.05 beta 1.7 alpha 0.7 zeta 0.5 Cv 0.08 testing 0 ls_max_tau 0.025 NB 200; mv out_H.txt out_H_testL11.txt
