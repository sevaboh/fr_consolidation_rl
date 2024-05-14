#include "toeplitz_mult.h"
double delta2(double j,double s,double a)
{
    if (fabs(a-2.0)<1e-6) return (((j-1)==s)?1.0:0.0);
    return -(pow(j-s-1.0,2.0-a)-pow(j-s,2.0-a));
}
double *_T=NULL,*M=NULL,*pD1=NULL,*pD2=NULL;
int sn=0;
void rl_der2(double *F,double *D,int n,double a,double dx,int toepl=0) // in openmp mode must be run inside "parallel" block
{
    double mult=pow(dx,-a)/Gamma(3.0-a);
    if (sn!=n)
    {
        if (_T) delete [] _T;
        if (M) delete [] M;
        if (pD1) delete [] pD1;
        if (pD2) delete [] pD2;
        _T=new double[n];
        M=new double[n];
        pD1=new double[n];
        pD2=new double[n];
    }
    if (sn!=n)
#pragma omp parallel for
    for (int i=1;i<n-1;i++)
    {
        _T[n-i-1]=(delta2(n+1,i,a)-2.0*delta2(n,i,a)+delta2(n-1,i,a));
	pD1[i]=delta2(i+1,i,a);
	pD2[i]=(delta2(i+1,i-1,a)-2.0*delta2(i,i-1,a));
    }
    if (sn!=n)
    	    sn=n;
    if (toepl)
	Toeplitz_mult(_T,F,M,n,1,0);
    for (int i=1;i<n-1;i++)
    {
	double ret=0.0;
	if ((i-2)>=0)
	{
	    if (toepl==0)
	    for (int s=0;s<=i-2;s++)
		ret+=F[s+1]*_T[i-s-1];
	    else
		ret+=M[i];
	}
	ret+=F[i+1]*pD1[i];
	ret+=F[i]*pD2[i];
	D[i]=ret*mult;
    }
}
