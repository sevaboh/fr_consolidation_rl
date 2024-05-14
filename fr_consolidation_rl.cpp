#include "gamma.h"
#include "inverse.h"
#include <map>
#include <vector>
#include <omp.h>
#include "libfftpack/fftpack.h"
#include "libfftpack/ls_fft.h"
#include "libfftpack/bluestein.h"
//#define OMP_HERE
#define OMP_NO_NESTED
//////////////////////////////////////
///// global variables ///////////////
//////////////////////////////////////
int N=50;  // number of nodes
int testing_mode=0; // testing with constructed solution
int debug_level=0;
int zoutstep=1;
// linear solver (TFQMR)
int ls_max_iter=30; // maximal number of iterations
double ls_eps=1e-12; // accuracy threshold
double ls_min_tau=0.00001; // minimal time step length, s
double max_tau=0.01; // maximal time step length, s
double ls_mult=1.25; // time step multiplier/divisor
double ls_percent=0.66; // threshold in per cent of max.iterations used to make a decision to change time step length
int total_iterations=0;
///////////////////////////////////////
////////// basic solver class /////////
///////////////////////////////////////

class basic_solver;
class basic_solver {
public:
	// storage
	double *b_U; // solution
	double *sb_U; // to save solution
	double *MM,*RP; // multiplication results and right part for TFQMR
	// matrix parts
	double *pr_A=NULL,*pr_B=NULL,*pr_R=NULL;

	// steps and auxiliary
	double tau; // current time step length 
	double T; // current time
	double L,dL; // domain depth and space variable step length
	// output
	FILE *out;
	// linear solver
	// tridiagonal part: linear system coefficients - A*U(i-1)+R*U(i)+B*U(i+1)=Rp
	virtual double A_(int i)=0;
	virtual double B(int i)=0;
	virtual double R(int i)=0;
	virtual double Rp(int i,double *U,double *pickard_prevU)=0;
	// pre-step action
	virtual void prestep()=0;
	// matrix multiplication of vector UU
	virtual void mmult(double *UU)=0;
	void mmult_main(double *UU)
	{
#ifdef OMP_HERE
#pragma omp for nowait
#endif
		for (int i = 1;i < N ;i++)
		{
			double rr=pr_R[i];
			MM[i]=UU[i-1]*pr_A[i]+UU[i+1]*pr_B[i]+UU[i]*rr;
			// normalization
			if (rr!=0.0)
				MM[i]/=rr;
		}
	}
	// precalculations for current time step
	virtual void __precalc(double *pickard_prevU)=0;
	void __precalc_main(double *pickard_prevU)
	{
#ifdef OMP_HERE
#pragma omp for nowait
#endif
		for (int i = 0;i <= N ;i++)
		{
			pr_A[i]=A_(i);
			pr_B[i]=B(i);
			pr_R[i]=R(i);
			// right part
			RP[i]=Rp(i,b_U,pickard_prevU);
			// normalization
			if (pr_R[i]!=0.0)
				RP[i]/=pr_R[i];
		}
	}
	// linear solver (TFQMR) (returns 1 if tau was increased at the end)
	int calc_step(double *pickard_prevU=NULL,int fix_tau=0)
	{
		if (pickard_prevU==NULL) pickard_prevU=b_U;
		// first step
		if (T==tau)
		{
		    for (int i =0;i <= N+1 ;i++)
		    {
			MM[i]=0.0;
			RP[i]=0.0;
		    }
		}
		//////////
		double *w=new double[(N+2)];
		double *y[2];
		y[0]=new double[(N+2)];
		y[1]=new double[(N+2)];
		double *rr=new double[(N+2)];
		double *v=new double[(N+2)];
		double *d=new double[(N+2)];
		double *x=new double[(N+2)];
		double theta=0,nu=0,tt=0,ro=0,ro1=0,c=0,b=0;
		double sgm,aa,rv;
		int n,m;
start:
		theta=nu=tt=ro=ro1=c=b=0.0;
		memcpy(x,b_U,(N+2)*sizeof(double));
#ifdef OMP_HERE
#pragma omp parallel
#endif
	    {
		////////// precalculations
		__precalc(pickard_prevU);
		mmult(x);
#ifdef OMP_HERE
#pragma omp for reduction(+:tt)
#endif
		for (int i=0;i<(N+2);i++)
		{
			// w_1=y_1=r_0=b-Ax_0
			// rr0: ro=rr_0*r_0!=0 - rr0=r0
			w[i]=y[0][i]=RP[i]-MM[i];
			// d=0
			d[i]=0;
			// tau=||r0||
			tt+=w[i]*w[i];
		}
#ifdef OMP_HERE
#pragma omp single
#endif
		tt=sqrt(tt);
#ifdef OMP_HERE
#pragma omp barrier
#endif
		// random rr0, ro0=rr_0*r_0
#ifdef OMP_HERE
#pragma omp for reduction(+:ro)
#endif
		for (int i=0;i<(N+2);i++)
		{
			rr[i]=tt*((rand() % 10000) / (10000.0 - 1.0));
			ro+=rr[i]*w[i];
		}
		// v=Ay_1
		mmult(y[0]);
#ifdef OMP_HERE
#pragma omp for
#endif
		for (int i=0;i<(N+2);i++)
			v[i]=MM[i];
#ifdef OMP_HERE
#pragma omp single
#endif
		n=1;
#ifdef OMP_HERE
#pragma omp barrier
#endif
loop1:
		{
			int br=0;
			// sigma_n-1 - rr_0*v_n-1
#ifdef OMP_HERE
#pragma omp single
#endif
			sgm=0;
#ifdef OMP_HERE
#pragma omp barrier
#pragma omp for reduction(+:sgm)
#endif
			for (int i=0;i<(N+2);i++)
				sgm+=rr[i]*v[i];
			// a_n-1=po_n-1/sigma_n-1
#ifdef OMP_HERE
#pragma omp single
#endif
			aa=ro/sgm;
#ifdef OMP_HERE
#pragma omp barrier
#endif
			// y_2n=y_2n-1 - a_n-1 * v_n-1
#ifdef OMP_HERE
#pragma omp for
#endif
			for (int i=0;i<(N+2);i++)
				y[1][i]=y[0][i]-aa*v[i];
#ifdef OMP_HERE
#pragma omp single
#endif
			m=2*n-1;
#ifdef OMP_HERE
#pragma omp barrier
#endif
loop2:
			{
				double ot=theta,onu=nu;
				// w_m+1=w_m-a_n-1Ay_m
				mmult(y[m-(2*n-1)]);
#ifdef OMP_HERE
#pragma omp for
#endif
				for (int i=0;i<(N+2);i++)
					w[i]=w[i]-aa*MM[i];
				// theta_m=||w_m+1||/tau_m-1; c_m=1/sqrt(1+theta_m^2)
#ifdef OMP_HERE
#pragma omp single
#endif
				theta=0;
#ifdef OMP_HERE
#pragma omp barrier
#pragma omp for reduction(+:theta)
#endif
				for (int i=0;i<(N+2);i++)
					theta+=w[i]*w[i];
#ifdef OMP_HERE
#pragma omp single
#endif
			    {
				theta=sqrt(theta)/tt;
				c=1.0/sqrt(1.0+theta*theta);
				// tau_m=tau_m-1 * theta_m *c_m; nu_m=c_m^2 *a_n-1
				tt=tt*theta*c;
				nu=c*c*aa;
				rv=0.0;
			    }
#ifdef OMP_HERE
#pragma omp barrier
#endif
				// d_m = y_m+(theta_m-1^2 nu_m-1 / a_n-1)*d_m-1
				// x_m=x_m-1 + nu_m *d_m
#ifdef OMP_HERE
#pragma omp for
#endif
				for (int i=0;i<(N+2);i++)
				{
					d[i]=y[m-(2*n-1)][i]+d[i]*(ot*ot*onu/aa);
					x[i]=x[i]+nu*d[i];
				}
				mmult(x);
#ifdef OMP_HERE
#pragma omp for reduction(+:rv)
#endif
				for (int i=0;i<(N+2);i++)
					rv+=(RP[i]-MM[i])*(RP[i]-MM[i]);
#ifdef OMP_HERE
#pragma omp single
#endif
				rv=sqrt(rv)/((N+2));
#ifdef OMP_HERE
#pragma omp barrier
#endif
				if (rv<ls_eps)
				{
				    br=1;
				    goto eloop2;
				}
			}
#ifdef OMP_HERE
#pragma omp single
#endif
			m++;
#ifdef OMP_HERE
#pragma omp barrier
#endif
			if (m<=2*n)
			    goto loop2;
eloop2:
			if (br==1)
				goto eloop1;
			// ro_n=rr0*w_2n+1, beta_n=ro_n/ro_n-1
#ifdef OMP_HERE
#pragma omp single
#endif
			ro1=0;
#ifdef OMP_HERE
#pragma omp barrier
#pragma omp for reduction(+:ro1)
#endif
			for (int i=0;i<(N+2);i++)
				ro1+=rr[i]*w[i];
#ifdef OMP_HERE
#pragma omp single
#endif
		    {
			b=ro1/ro;
			ro=ro1;
		    }
#ifdef OMP_HERE
#pragma omp barrier
#endif
			// y_2n+1 = w_2n+1+beta_n*y_2n
#ifdef OMP_HERE
#pragma omp for
#endif
			for (int i=0;i<(N+2);i++)
				y[0][i]=w[i]+b*y[1][i];
			// v_n=Ay_2n+1+b*(Ay_2n+b*v_n-1)
			mmult(y[1]);
#ifdef OMP_HERE
#pragma omp for
#endif
			for (int i=0;i<(N+2);i++)
				v[i]=b*(MM[i]+b*v[i]);
			mmult(y[0]);
#ifdef OMP_HERE
#pragma omp for
#endif
			for (int i=0;i<(N+2);i++)
				v[i]=MM[i]+v[i];
		}
#ifdef OMP_HERE
#pragma omp single
#endif
		n++;
#ifdef OMP_HERE
#pragma omp barrier
#endif
		if (n<ls_max_iter)
		    goto loop1;
eloop1:;
	    }
		total_iterations+=n;
		// change tau and recalc if needed
		if (fix_tau==0)
		if (n==ls_max_iter)
		{
		     if (tau<ls_min_tau)
		     {
			 x[0]=NAN;
			 printf("minimal tau value reached\n");
			 exit(2);
		     }
		     else
		     {
		    	T-=tau;
		    	tau/=ls_mult;
		    	T+=tau;
			// debug output
			if (debug_level==2) printf("%p r n %d r %g T %g tau %g sol %g %g %g %g %g %g\n",this,n,rv,T,tau,x[0],x[1],x[2],b_U[0],b_U[1],b_U[2]);
		        goto start;
		     }
		}
		// save solution and free memory
		memcpy(b_U,x,(N+2)*sizeof(double));
		delete [] w;
		delete [] y[0];
		delete [] y[1];
		delete [] rr;
		delete [] v;
		delete [] d;
		delete [] x;
		// increase tau if needed and increase T
		int ret=0;
		if (fix_tau==0)
		if (n<ls_max_iter*ls_percent)
			if (tau<max_tau)
			{
				 tau*=ls_mult;
				 ret=1;
			}
		T+=tau;
		// debug output
		if (debug_level==2) printf("%p e n %d rv %g T %g tau %g\n",this,n,rv,T,tau);
		return ret;
	}
	// perform Pickard iteration 
	int pickard_calc_step(basic_solver **other=NULL,int nother=0)
	{
	    double *prev_U=new double[N+2];
	    double *init_U=new double[N+2];
	    double **prev_other_U=NULL,**init_other_U=NULL;
	    double diff=0;
	    int iter=0;
	    int mul=0;
	    int *other_muls=NULL;
	    prestep();
	    memcpy(init_U,b_U,sizeof(double)*(N+2));
restart:
	    memcpy(b_U,init_U,sizeof(double)*(N+2));
	    if (init_other_U)
	    for (int i=0;i<nother;i++)
		memcpy(other[i]->b_U,init_other_U[i],sizeof(double)*(N+2));
	    iter=0;
	    mul=calc_step();
	    if (mul==-1) return 0;
	    if (mul)
	    {	
		T-=tau;
		tau/=ls_mult;
		T+=tau;
	    }
	    if (nother) // run other solvers
	    {
		if (prev_other_U==NULL) { prev_other_U=new double *[nother]; memset(prev_other_U,0,nother*sizeof(void *));}
		if (init_other_U==NULL) { init_other_U=new double *[nother]; memset(init_other_U,0,nother*sizeof(void *));}
		if (other_muls==NULL) other_muls=new int[nother];
		for (int i=0;i<nother;i++)
		{
			if (prev_other_U[i]==NULL) prev_other_U[i]=new double [N+2];
			if (init_other_U[i]==NULL) init_other_U[i]=new double [N+2];
			memcpy(init_other_U[i],other[i]->b_U,sizeof(double)*(N+2));
			other_muls[i]=other[i]->calc_step();
			if (other_muls[i]==-1) return 0;
			if (other_muls[i])
			{	
				other[i]->T-=other[i]->tau;
				other[i]->tau/=ls_mult;
				other[i]->T+=other[i]->tau;
		        }
		}
		// find minimal tau
		double min_tau=tau;
		for (int i=0;i<nother;i++)
			if (other[i]->tau<min_tau)
				min_tau=other[i]->tau;
		// solve all once more with fixed minimal tau
		memcpy(b_U,init_U,sizeof(double)*(N+2));
		for (int i=0;i<nother;i++)
			memcpy(other[i]->b_U,init_other_U[i],sizeof(double)*(N+2));
		T-=2.0*tau;
		tau=min_tau;
		T+=tau;
		mul=calc_step(NULL,1);
		if (mul==-1) return 0;
		for (int i=0;i<nother;i++)
		{
			other[i]->T-=2.0*other[i]->tau;
			other[i]->tau=min_tau;
			other[i]->T+=other[i]->tau;
			other[i]->calc_step(NULL,1);
		}
	    }
	    // debug
	    if (debug_level==1)
	    {
		 printf("%p pickard init tau %g T %g other %d ",this,tau,T-tau,nother);
		for (int i=0;i<nother;i++) printf(" %d tau %g T %g ",i,other[i]->tau,other[i]->T-other[i]->tau);
		printf("\n");
	    }
	    do
	    {
		// save solution on previous iteration and restore initial values
		memcpy(prev_U,b_U,sizeof(double)*(N+2));
		memcpy(b_U,init_U,sizeof(double)*(N+2));
		for (int i=0;i<nother;i++)
		{
			memcpy(prev_other_U[i],other[i]->b_U,sizeof(double)*(N+2));
			memcpy(other[i]->b_U,init_other_U[i],sizeof(double)*(N+2));
		}
		// solve on next iteration
		T-=tau;
		int mm=calc_step(prev_U,1);
		if (mm==-1) return 0;
		for (int i=0;i<nother;i++)
		{
		    other[i]->T-=other[i]->tau;
		    other[i]->calc_step(prev_other_U[i],1);
		}
		// calculate difference
		diff=0;
		for (int i=0;i<N+1;i++)
		    diff+=(b_U[i]-prev_U[i])*(b_U[i]-prev_U[i]);
		for (int j=0;j<nother;j++)
		    for (int i=0;i<N+1;i++)
			diff+=(other[j]->b_U[i]-prev_other_U[j][i])*(other[j]->b_U[i]-prev_other_U[j][i]);
		diff=diff/(N*(nother?nother:1));
		// debug
		if (debug_level==1) printf("%p pickard iter %d diff %g T %g S %g %g %g - %g %g %g\n",this,iter,diff,T-tau,b_U[0],b_U[1],b_U[2],b_U[N-2],b_U[N-1],b_U[N]);
		if ((iter++)> ls_max_iter) 
			break;
	    }
	    while (diff>ls_eps);
	    if (!isfinite(diff))
	    {
		ls_max_iter/=2.0;
		goto restart;
	    }
	    if (iter>=ls_max_iter) // decrease step
	    {
		 if (tau<ls_min_tau)
		 {
			 b_U[0]=NAN;
			 printf("minimal tau value reached\n");
			 exit(2);
		 }
		 else
		 {
		    T-=tau;
		    tau/=ls_mult;
		    T+=tau;
		    if (debug_level==1) printf("%p pickard restart niter %d diff %g tau %g T %g\n",this,iter,diff,tau,T-tau);
		    goto restart;
		}
	    }
	    // increase tau if needed
	    if (mul)
	    {
		T-=tau;
		tau*=ls_mult;
		T+=tau;
	    }
	    if (nother)
		for (int i=0;i<nother;i++)
		    if (other_muls[i])
		    {
			other[i]->T-=other[i]->tau;
			other[i]->tau*=ls_mult;
			other[i]->T+=other[i]->tau;
		    }
	    // debug
	    if (debug_level==1) printf("%p pickard niter %d diff %g tau %g T %g\n",this,iter,diff,tau,T-tau);
	    // clean up
	    delete [] prev_U;
	    delete [] init_U;
	    if (nother)
	    {
		for (int i=0;i<nother;i++)
		{
		    delete [] prev_other_U[i];
		    delete [] init_other_U[i];
		}
		delete [] prev_other_U;
		delete [] init_other_U;
		delete [] other_muls;
	    }
	    return 1;
	}
	// constructor
	basic_solver(double _L)
	{
		b_U=new double[(N+2)];
		sb_U=new double[(N+2)];
		MM=new double[(N+2)];
		RP=new double[(N+2)];

		L = _L;
		dL = L/N;
		tau = max_tau;
		T=tau;
	}
	// desctructor
	~basic_solver()
	{
		delete [] b_U;
		delete [] sb_U;
		delete [] MM;
		delete [] RP;
		if (pr_A) delete [] pr_A;
		if (pr_B) delete [] pr_B;
		if (pr_R) delete [] pr_R;
	}
	// output
	virtual void output_solution()=0;
};
//////////////////////////////////////
// Solver for consolidation equation 
// dH/dt + zeta * D[CF][t,alpha]H=Cv * D[RL][x,beta]H
// 0<alpha<1, 1<beta<2
//////////////////////////////////////
class H_solver: public basic_solver {
public:
#include "rl_findiff.h"
	// parameters
	double Cv; 
	double zeta;
	double alpha,beta;
	// auxiliary
	double *b_prevH;
	double *b_RL;
	double *b_CF;
	double pr_etau;
	double pr_g;
	// testing - constructed solution - testing_mode==1
	double constructed_solution(double z,double t) // solution for testing mode
	{
		return t*t*z*z*z;
	}
	double testing_F(double z,double t)
	{
		double ret=-Cv*t*t*pr_g*pow(z,3.0-beta);
		ret+=2*t*z*z*z;
		ret+=zeta*z*z*z*(2/(alpha*alpha))*((1.0-alpha)*exp(-(alpha/(1.0-alpha))*t)+alpha*t+alpha-1);
		return ret;
	}
	// testing - analytical solution for beta=2 - testing_mode==2
	static double analytic_h(double x) // initial
	{
		return 4.0*x*(1.0-x);
	}
	static double lambda_n(int n)
	{
		return -n*n*M_PI*M_PI;
	}
	static double w_n(int n,double x)
	{
		if ((n%2)==0) return 0.0;
		return sqrt(2.0)*sin(n*M_PI*x);
	}
	std::map<int, double> save_h_n;
	// integrate by linear approximation and recursive subdivision
	double integrate(double (*func)(double a,double b,double z,double r),double a,double b,double z,double r0,double r1,double eps,int sd=1,double prev=0.0)
	{
		double v1=0.5*(r1-r0)*(func(a,b,z,r0)+func(a,b,z,r1));
		double v=0.0;
		for (int i=1;i<=sd;i++)
			v+=func(a,b,z,r0+(i/((double)(sd+1)))*(r1-r0));
		v*=(r1-r0);
		v+=v1;
		v/=(double)(sd+1);
		if (sd<50000) // recursive subdivition
		{
			if (sd==1)
				prev=v1;
			if (fabs(v-prev)<eps)
				return v;
			else
				return integrate(func,a,b,z,r0,r1,eps,sd*2.0,v);
		}
		return v;
	}
	static double f_n_func(double a,double b,double z,double r)
	{
		return analytic_h(r)*w_n((int)z,r);
	}
	double h_n(int n)
	{
		double ret=0.0;
		std::map<int,double>::iterator ii;
		if ((ii=save_h_n.find(n))!=save_h_n.end()) return ii->second; 
		ret=integrate(f_n_func,0.0,0.0,n,0,1,ls_eps,fabs(n)+1);
		save_h_n.insert(std::pair<int,double>(n,ret));
		return ret;
	}
	double delta_n(int n)
	{
		return pow(alpha+zeta-Cv*lambda_n(n)*(1.0-alpha),2.0)+4*alpha*Cv*lambda_n(n)*(1.0-alpha);
	}
	double s_n(int n,int i)
	{
		double ret=-(alpha+zeta-Cv*lambda_n(n)*(1.0-alpha));
		if (i==1) ret+=sqrt(delta_n(n)); else ret-=sqrt(delta_n(n));
		return ret/(2.0*(1.0-alpha));
	}
	double fi_n(int n,int i)
	{
		return s_n(n,i)+(alpha+zeta)/(1.0-alpha);
	}
	double u_n(int n,double t)
	{
		return (h_n(n)/(fi_n(n,1)-fi_n(n,2)))*(fi_n(n,1)*exp(s_n(n,1)*t)-fi_n(n,2)*exp(s_n(n,2)*t));
	}
	double analytic_solution(double z,double t)
	{
		int n=1;
		double ret=0;
		do
		{
		    double v=u_n(n,t)*w_n(n,z)+u_n(n+1,t)*w_n(n+1,z);
		    ret+=v;
		    if (fabs(v)<ls_eps)
			break;
		    n+=2;
		}
		while(1);
		return ret;
	}
	// testing - main solution function
	double testing_f(double z,double t) // solution for testing mode
	{
		if (testing_mode==1) return constructed_solution(z,t);
		if (testing_mode==2) return analytic_solution(z,t);
		return 0.0;
	}
	// upper boundary condition (first order)
	double Uc()
	{
	    if (testing_mode==1)
		return testing_f(0,T);
	    return 0.0;
	}
	// coefficients of three-diagonal linear equations system
	double A_(int i)
	{
		return 0.0;
	}
	double B(int i)
	{
		return 0.0;
	}
	double R(int i)
	{
		return (1.0/tau)+(zeta/(alpha*tau))*(1.0-pr_etau);
	}
	// right part 
	double Rp(int ii,double *b_Uold,double *pickard_prevU)
	{
		double ret=R(ii)*b_Uold[ii];
		ret-=zeta*pr_etau*b_CF[ii];
		if (testing_mode==1) ret+=testing_F(ii*dL,T);
		return ret;
	}
	// linear solver
	// matrix multiplication of vector UU
	void mmult(double *UU)
	{
		mmult_main(UU);
		// Riemann-Liouville derivative upon z
#ifdef OMP_HERE
#pragma omp single
#endif
		rl_der2(UU,b_RL,N+1,beta,dL,1);
#ifdef OMP_HERE
#pragma omp for
#endif
		for (int i=1;i<N;i++)
			MM[i]-=Cv*b_RL[i]/((pr_R[i]==0.0)?1.0:pr_R[i]);
		// first order boundary conditions
#ifdef OMP_HERE
#pragma omp single
#endif
		{
			MM[N]=UU[N];
			MM[0]=UU[0];
		}
	}
	void prestep()
	{
	    // update CF, save prevU
	    double etau=exp(-alpha*tau/(1.0-alpha));
	    for (int i=0;i<N;i++)
	    {
	    	if (T!=tau) b_CF[i]=etau*b_CF[i]+(1.0/(alpha*tau))*(1.0-etau)*(b_U[i]-b_prevH[i]); else b_CF[i]=0.0;
		b_prevH[i]=b_U[i];
	    }
	}
	// precalculations for current time step
	void __precalc(double *pickard_prevU)
	{
#ifdef OMP_HERE
#pragma omp single
#endif
	    {
	    memcpy(sb_U,b_U,sizeof(double)*(N+2));
	    memcpy(b_U,pickard_prevU,sizeof(double)*(N+2));
	    if (pr_A==NULL) pr_A=new double[N+2];
	    if (pr_B==NULL) pr_B=new double[N+2];
	    if (pr_R==NULL) pr_R=new double[N+2];
	    pr_etau=exp(-alpha*tau/(1.0-alpha));
	    memcpy(b_U,sb_U,sizeof(double)*(N+2));
	    }
	    // precalc linear system values
	    __precalc_main(pickard_prevU);
	    // boundary conditions
#ifdef OMP_HERE
#pragma omp single
#endif
	    {
	    RP[N]=0;
	    RP[0]=0;
	    if (testing_mode==1)
	        RP[N]=testing_f(N*dL,T);
	    }
	}
	// constructor
	H_solver(double *p) : Cv(p[3]),zeta(p[4]),alpha(p[5]),beta(p[6]),basic_solver(p[1])
	{
		
		b_prevH=new double[N+2];
		b_RL=new double[N+2];
		b_CF=new double[N+2];
		pr_g=Gamma(4)/Gamma(4-beta);
		// initial conditions
		for (int i = 0;i < N + 1;i++)
		{
			b_U[i] = p[2];
			if (testing_mode==1) b_U[i]=testing_f(i*dL,0);
			if (testing_mode==2) b_U[i]=analytic_h(i*dL);
		}
		if (testing_mode==0) b_U[0]=b_U[N]=0.0;
		//printf("init cv %g z %g a %g b %g L %g h0 %g\n",Cv,zeta,alpha,beta,L,p[1]);
		out=fopen("out_H.txt","wt");
	}
	// desctructor
	~H_solver()
	{
		delete [] b_prevH;
		delete [] b_RL;
		delete [] b_CF;
		fclose(out);
	}
	int process(double *&xs,double *&ts,double *&hs,int &nv,double eT,int print)
	{
		std::vector<double> xsv,tsv,hsv;
		nv=0;
		while (T<eT)
		{
			if (pickard_calc_step()==0) return 0;
			for (int i=0;i<=N;i++)
			{
				xsv.push_back(i*dL);
				tsv.push_back(T-tau);
				hsv.push_back(b_U[i]);
			}
			if (print) output_solution();
		}
		xs=new double [xsv.size()];
		ts=new double [xsv.size()];
		hs=new double [xsv.size()];
		nv=xsv.size();
		for (int i=0;i<xsv.size();i++)
		{
			xs[i]=xsv[i];
			ts[i]=tsv[i];
			hs[i]=hsv[i];
		}
		if (print) output_solution();
		return 1;
	}
	// output
	void output_solution()
	{
		fprintf(out,"t %g tau %g zeta %g cv %g a %g b %g - H: ",(T-tau),tau,zeta,Cv,alpha,beta);
		if (T==tau)
		{
		    for (int i=0;i<N+1;i+=zoutstep)
			fprintf(out,"%g ",i*dL);
		    fprintf(out,"\n");
		    fprintf(out,"t %g tau %g zeta %g cv %g a %g b %g - H: ",(T-tau),tau,zeta,Cv,alpha,beta);
		}
		for (int i=0;i<N+1;i+=zoutstep)
			fprintf(out,"%g ",b_U[i]);
		if (testing_mode)
		{
		    double diff=0,sq1=0,sq2=0;
		    for (int i=0;i<N+1;i++)
		    {
			diff+=(b_U[i]-testing_f(i*dL,T-tau))*(b_U[i]-testing_f(i*dL,T-tau));
			sq1+=b_U[i]*b_U[i];
			sq2+=testing_f(i*dL,T-tau)*testing_f(i*dL,T-tau);
//			fprintf(out,"%d %g %g\n",i,b_U[i],testing_f(i*dL,T-tau));
		    }
		    diff=sqrt(diff/N);
		    sq1=sqrt(sq1);
		    sq2=sqrt(sq2);
		    fprintf(out,"testing diff %g rel %g",diff,fabs(sq1-sq2)/sq2);
		    printf ("testing diff %g rel %g\n",diff,fabs(sq1-sq2)/sq2);
		}
		fprintf(out,"\n");
	}
};
// for inverse
int nparams=6;
int mask=0x1F<<2;
double bv[6][3]={{1,1,1},{1,1,1},{0.001,0.1,0.1},{0.1,5,2.5},{0.5,1,0.5},{1.5,2,1.5}};
////////////////////////////////////////
////////////////// main/////////////////
////////////////////////////////////////
int main(int argc,char **argv)
{
	double Tm = 5; // ending time, s
	double Om = 0.1; // interval for solution output
	double H0=1; // initial H, m water
	double L=3; // domain depth, m
	double zeta=0.5;
	double Cv=0.02;
	double alpha=0.7;
	double beta=1.7;
	int inverse=0;
	char *inv_file=NULL;
	int pso_n=10;
	double pso_o=0.3;
	double pso_fi_p=0.3;
	double pso_fi_g=0.3;
	double pso_eps=1e-10;
	int pso_max_iter=10;
	double rp=0.2;
	int div_t=10;
	int div_x=10;
	double noise=0.05;
#ifndef OMP_NO_NESTED
	omp_set_nested(1);
#endif
	// reading basic parameters
	for (int i=1;i<argc;i+=2)
	{
		// basic parameters
		if (strcmp(argv[i],"Tm")==0)
			Tm=atof(argv[i+1]);
		if (strcmp(argv[i],"Om")==0)
			Om=atof(argv[i+1]);
		if (strcmp(argv[i],"NB")==0)
			N=atoi(argv[i+1]);
		if (strcmp(argv[i],"H0")==0)
			H0=atof(argv[i+1]);
		if (strcmp(argv[i],"LB")==0)
			L=atof(argv[i+1]);
		if (strcmp(argv[i],"testing")==0)
			testing_mode=atoi(argv[i+1]);
		if (strcmp(argv[i],"debug_level")==0)
			debug_level=atoi(argv[i+1]);
		if (strcmp(argv[i],"Zoutstep")==0)
			zoutstep=atoi(argv[i+1]);
		// linear solver parameters
		if (strstr(argv[i],"ls_eps")!=NULL)
		    ls_eps=atof(argv[i+1]);
		if (strstr(argv[i],"ls_max_iter")!=NULL)
		    ls_max_iter=atoi(argv[i+1]);
		if (strstr(argv[i],"ls_percent")!=NULL)
		    ls_percent=atof(argv[i+1]);
		if (strstr(argv[i],"ls_min_tau")!=NULL) //  in seconds
		    ls_min_tau=atof(argv[i+1]);
		if (strstr(argv[i],"ls_mult")!=NULL)
		    ls_mult=atof(argv[i+1]);
		if (strstr(argv[i],"max_tau")!=NULL) // in seconds
		    max_tau=atof(argv[i+1]);
		// model parameters
		if (strcmp(argv[i],"zeta")==0)
			zeta=atof(argv[i+1]);
		if (strcmp(argv[i],"Cv")==0)
			Cv=atof(argv[i+1]);
		if (strcmp(argv[i],"alpha")==0)
			alpha=atof(argv[i+1]);
		if (strcmp(argv[i],"beta")==0)
			beta=atof(argv[i+1]);
		// inverse
		if (strcmp(argv[i],"inverse")==0)
		{
			inverse=1;
			inv_file=argv[i+1];
		}
		if (strcmp(argv[i],"pso_n")==0)
			pso_n=atoi(argv[i+1]);
		if (strcmp(argv[i],"pso_o")==0)
			pso_o=atof(argv[i+1]);
		if (strcmp(argv[i],"pso_fi_p")==0)
			pso_fi_p=atof(argv[i+1]);
		if (strcmp(argv[i],"pso_fi_g")==0)
			pso_fi_g=atof(argv[i+1]);
		if (strcmp(argv[i],"pso_max_iter")==0)
			pso_max_iter=atoi(argv[i+1]);
		if (strcmp(argv[i],"pso_eps")==0)
			pso_eps=atof(argv[i+1]);
		if (strcmp(argv[i],"pso_rp")==0)
			rp=atof(argv[i+1]);
		if (strcmp(argv[i],"pso_mask")==0)
			mask=atoi(argv[i+1]);
		if (strcmp(argv[i],"div_t")==0)
			div_t=atoi(argv[i+1]);
		if (strcmp(argv[i],"div_x")==0)
			div_x=atoi(argv[i+1]);
		if (strcmp(argv[i],"noise")==0)
			noise=atof(argv[i+1]);
	}
	printf("Grid size - %d max_tau %g, Tend %g Tsave %g div_t %d div_x %d noise %g\n",N,max_tau,Tm,Om,div_t,div_x,noise);
	// solve inverse problem
	if (inverse==1)
	{
		bv[0][0]=bv[0][1]=bv[0][2]=L;
		bv[1][0]=bv[1][1]=bv[1][2]=H0;
		H_solver *s=run_inverse<H_solver>(inv_file,pso_n,pso_o,pso_fi_p,pso_fi_g,pso_eps,pso_max_iter,rp,nparams,mask,bv,div_x,div_t,noise,Tm);
		if (s) delete s;
		return 0;
	}
	// solve direct problem
	// creating solvers
	if (testing_mode==2) L=1.0;
	double p[]={0,L,H0,Cv,zeta,alpha,beta};
	H_solver *solver=new H_solver(p);
	// run simulation
	double last_o=0;
	solver->output_solution();
	while (solver->T<Tm)
	{
		if (solver->pickard_calc_step()==0) return 0;
		if ((solver->T-last_o)>=Om)
		{
		    last_o=solver->T;
		    solver->output_solution();
		    printf("T %g total_iters %d\n",solver->T-solver->tau,total_iterations);
		}
	}
	solver->output_solution();
	return 0;
}
