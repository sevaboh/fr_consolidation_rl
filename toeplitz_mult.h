#include "libfftpack/fftpack.h"
#include "libfftpack/ls_fft.h"
#include "libfftpack/bluestein.h"
real_plan plan1;
int old_n1 = -1;
real_plan plan2;
int old_n2 = -1;
void fft(double *v, int n)
{
	if (old_n1 != n)
	{
		if (old_n1!=-1)
			kill_real_plan(plan1);
		plan1 = make_real_plan(n);
		old_n1 = n;
	}
    real_plan_forward_c(plan1, v);
    
}
void ifft(double *v, int n)
{
	if (old_n2 != n)
	{
		if (old_n2 != -1)
			kill_real_plan(plan2);
		plan2 = make_real_plan(n);
		old_n2 = n;
	}
	real_plan_backward_c(plan2, v);
}
// z=T~v 
//full=1, space_der=1 (T~=(-T[N-1],...,-T[1],T[0],T[1],...,T[N-1])
//full=0, space_der=1 (T~=(-T[N-1],...,-T[1],0.0,T[1],...,T[N-1])
//full=1, space_der=0 (T~=(0.0,...,0.0,T[0],T[1],...,T[N-1])
//full=0, space_der=0 (T~=(0.0,...,0.0,0.0,T[1],...,T[N-1])
double *tmps[4] = { NULL,NULL,NULL,NULL }; // 2n
double *old_T=NULL;
int tmp_n = 0;
void Toeplitz_mult(double *T, double *v, double *z,int N,int full,int space_der)
{
	if (tmp_n != N)
	{		
		for (int i = 0;i < 4;i++)
		{
			if (tmps[i]) delete[] tmps[i];
			tmps[i] = new double[5*N];
			memset(tmps[i], 0, 5 * N*sizeof(double));
		}
		tmp_n = N;
	}
	if (old_T!=T)
	{
	    for (int i = 0;i < N;i++)
	    {
			tmps[full][2 * i] = T[i];
			tmps[full][2 * i+1] = 0.0;
			if (space_der == 1)
				tmps[full][2 * (N + i)] = -T[N - i];
			else
				tmps[full][2 * (N + i)] = 0.0;
			tmps[full][2 * (N + i)+1] = 0.0;
	    }
	    if (full == 0)
			tmps[full][0] = 0.0;
	    tmps[full][2 * N] = 0.0;
	    tmps[full][2 * N+1] = 0.0;
	    // c=FFT(c)
	    fft(tmps[full], 2 * N );
	}
	// t1=(v,0)
	for (int i = 0;i < N;i++)
	{
		tmps[2][2*i] = v[i];
		tmps[2][2*i+1] = 0;
		tmps[2][2*(N + i)] = 0;
		tmps[2][2*(N + i)+1] = 0;
	}
	// v=FFT(v)
	fft(tmps[2], 2* N);
	// t2=c x v
#pragma omp parallel for
	for (int i = 0;i < 2 * N;i++)
	{
		tmps[3][2*i] = tmps[full][2*i] * tmps[2][2 * i] - tmps[full][2*i+1] * tmps[2][2*i+1]; // Re
		tmps[3][2*i+1] = tmps[full][2 * i] * tmps[2][2*i+1] + tmps[full][2*i+1] * tmps[2][2*i]; //Im
	}
	// t2=IFFT(t2)
	ifft(tmps[3],  2*N);
#pragma omp parallel for
	for (int i = 0;i < 2*2 * N;i++)
		tmps[3][i] /= 2 * N;
	// z[i]=tmps[3][i],i=0,...,N-1
	for (int i = 0;i < N;i++)
		z[i] = tmps[3][2*i];
}
