# fr_consolidation_rl
Solver for space-time-fractional filtration-consolidation equation with Caputo-Fabrizio and Riemann-Liouville derivatives.

Source code:
fr_consolidation_rl.cpp - main source for the solver;
gamma.h - code to calculate Gamma function values;
inverse.h - code for PSO inverse problem solver;
rl_findiff.h - code for the calculation of the Riemann-Liouville derivative w.r.t. the space variable representing it through Toeplitz matrices and using libfftpack to perform the multiplications of Toeplitz matrices on vectors;
toeplitz_mult.h - code for Toplitz matrix on vector multiplication;
libfftpack - OpenMP-parallelized version of libfftpack;
libfftpack_orig - original version of libfftpack.

Scripts:
test_accuracy.sh - runs tests for direct problem solution accuracy;
test_inv_gen_input.sh - generates input files for the testing of inverse problem solution;
test_inv.sh - runs tests of inverse problem solution;
test_speed.sh - runs tests for multithreaded performance of direct problem solver.

