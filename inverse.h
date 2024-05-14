#include <sys/time.h>
#include <ctime>
#include <unistd.h>
#include <sys/times.h>
double rnd()
{
	return ((rand() % 10000) / (10000.0 - 1.0));
}
unsigned int GetTickCount()
{
   struct tms t;
   long long time=times(&t);
   int clk_tck=sysconf(_SC_CLK_TCK);
   return (unsigned int)(((long long)(time*(1000.0/clk_tck)))%0xFFFFFFFF);    
}
//////inverse PSO solver//////////////////
template <class T> class PSO_solver
{
public:
	double restart_prox=0.2; // probability of particles "restart" - random initialization
	int n_fit_parameters; // total number of parameters
	int mask_fit; // bit mask of parameters to fit
	double **bounds_fit; // [min,max,default]

	PSO_solver(double rp,int nfp,int mf,double bf[][3])
	{
		restart_prox=rp;
		n_fit_parameters=nfp;
		mask_fit=mf;
		bounds_fit=new double*[nfp];
		for (int i=0;i<nfp;i++)
		{
			bounds_fit[i]=new double[3];
			for (int j=0;j<3;j++) bounds_fit[i][j]=bf[i][j];
		}
	}
	void init_particle(double *particle)
	{
		for (int i=0;i<n_fit_parameters;i++)
			if (mask_fit & (1<<i))
			{
			    if ((bounds_fit[i][0]>0)&&(bounds_fit[i][1]>0))
				particle[1+i]=exp(log(bounds_fit[i][0])+(log(bounds_fit[i][1])-log(bounds_fit[i][0]))*rnd());
			    else
				particle[1+i]=bounds_fit[i][0]+(bounds_fit[i][1]-bounds_fit[i][0])*rnd();
			}
			else
				particle[1+i]=bounds_fit[i][2];
	}
	T *init_solver(double *p)
	{
		for (int i=0;i<n_fit_parameters;i++)
		{
			if (p[i+1]<bounds_fit[i][0])
				p[i+1]=bounds_fit[i][0];
			if (p[i+1]>bounds_fit[i][1])
				p[i+1]=bounds_fit[i][1];
		}
		return new T(p);
	}
	// goal function - average absolute difference between modelled and experimental data
	// - p - array of parameters' values
	// - if solv!=NULL returns in it the created solver with the given parameters' values 
	double calc_err(double *p,double *values_x,double *values_t,double *values_h,int nvalues,int size,int print,double _maxT=-1.0,T**solv=NULL)
	{
		T *ss=NULL;
		double *xs=NULL,*ts=NULL,*hs=NULL;
		double err=0.0;
		int nfound=0;
		int nv=0,runs;
		unsigned int t1=GetTickCount();
		FILE *fi;
		if (print) fi=fopen("inverse_log.txt","wt");
		// solve
		ss=init_solver(p);
		// calc max T
		double maxT=0.0;
		for (int i=0;i<nvalues;i++) if (values_t[i]>maxT) maxT=values_t[i];
		if (_maxT!=-1.0) maxT=_maxT;
		maxT+=ss->tau;
		runs=ss->process(xs,ts,hs,nv,maxT,print);
		if (nv==0)
		{
			err=1e300;
			goto end;
		}
		// calc error
		// loop on inputs
		for (int j=0;j<nvalues;j++)
		{
			// loop on simulated data - find closest x,t
			int cl_k=-1;
			double md=1e300;
			for (int k=0;k<nv;k++)
			{
				double d=(values_t[j]-ts[k])*(values_t[j]-ts[k])+(values_x[j]-xs[k])*(values_x[j]-xs[k]);
				if (d<md)
				{
				    md=d;
				    cl_k=k;
				}
			}
			if (cl_k!=-1)
			{
				err+=(values_h[j]-hs[cl_k])*(values_h[j]-hs[cl_k]);
				if (print) 
					fprintf(fi,"%d %d x %g t %g h %g simh %g x %g t %g\n",j,cl_k,xs[cl_k],ts[cl_k],hs[cl_k],values_h[j],values_x[j],values_t[j]);
				//printf("%d %d x %g t %g h %g simh %g x %g t %g\n",j,cl_k,xs[cl_k],ts[cl_k],hs[cl_k],values_h[j],values_x[j],values_t[j]);
				nfound++;
			}
			if (!finite(err))
			{
				err=1e300;
				break;
			}
		}
end:
		delete [] xs;
		delete [] ts;
		delete [] hs;
		if (solv==NULL)
			delete ss;
		else
			solv[0]=ss;
		if ((err!=1e300)&&(nfound))
			err=sqrt(err/nfound); // average absolute error in terms of current I
		if (print)
			fclose(fi);
		// log
		char str[4096]="";
		char str2[1028]="";
		strcat(str,"param values: ");
		for (int i=0;i<n_fit_parameters;i++)
		if (mask_fit & (1<<i))
		{
		    if (size)
			sprintf(str2,"(%g %g) ",p[1+i],p[1+size+i]);
		    else
			sprintf(str2,"%g ",p[1+i]);
		    strcat(str,str2);
		}
		sprintf(str2,"err - %g time %d runs %d init_time %d thread %d\n",err,GetTickCount()-t1,runs,t1,omp_get_thread_num());
		strcat(str,str2);
		//printf("%s",str);
		return err;
	}
	void solve_all_and_test(double **p,int nparticles,double *values_x,double *values_t,double *values_h,int nvalues,int size)
	{
#pragma omp parallel for
		for (int i=0;i<nparticles;i++)
			p[i][0]=calc_err(p[i],values_x,values_t,values_h,nvalues,size,0);
	}
	void update_particle(int size,double **particles,int i,double *best_p,int pso_n,double pso_o,double pso_fi_p,double pso_fi_g,double pso_eps,int iter,double u)
	{
		// update velocity
		for (int j=0;j<size;j++)
		if (mask_fit & (1<<j))
		{
			double rp=(rand()%10000)/10000.0;
			double rg=(rand()%10000)/10000.0;
			particles[i][j+1+size]=pso_o*particles[i][j+1+size]+pso_fi_p*rp*(particles[i][j+1+2*size]-particles[i][j+1])+pso_fi_g*rg*(best_p[j+1]-particles[i][j+1]);    
		}
		// update position
		double psi=exp(particles[i][0]/u)/pow(1.0+exp(-particles[i][0]/u),iter); // DW-PSO
		for (int j=0;j<size;j++)
		if (mask_fit & (1<<j))
		{
			if (restart_prox>=0) 
				particles[i][1+j]+=particles[i][j+1+size];
			else
			{
				// DW-PSO
				double r=(rand()%10000)/10000.0;
				particles[i][1+j]=psi*particles[i][1+j]+(1-psi)*particles[i][j+1+size]+r*psi*best_p[j+1];
			}
		}
		// restart particle
		double rp=(rand()%10000)/10000.0;	
		if (rp<((pso_o<1)?(1-pso_o):0.0)*fabs(restart_prox))
		{
			init_particle(particles[i]);
			printf("r");
		}
	}
	double update_global(int size,double **particles,double *best_p,int pso_n,int iter,int print=1)
	{
		for (int i=0;i<pso_n;i++)
		{
			// update bests
			if (particles[i][0]<particles[i][1+3*size])
			{
				for (int j=0;j<size;j++)
					particles[i][j+1+2*size]=particles[i][j+1];
				particles[i][1+3*size]=particles[i][0];
			}
			if (particles[i][0]<best_p[0])
			{
				for (int j=0;j<size;j++)
					best_p[j+1]=particles[i][j+1];
				best_p[0]=particles[i][0];
			}
		}
		// check best-worst difference
		double max = 0.0;
		double avg=0.0;
		for (int i = 0;i < pso_n;i++)
		{
			if (max < particles[i][0])
				max = particles[i][0];
			avg+= particles[i][0];
		}
		avg/=pso_n;
		if (print)
		{
			printf("%d avg %g best: ", iter,avg);
			for (int j = 0;j <= size;j++)
				printf("%g ", best_p[j]);
			printf("\n");
		}
		return max;
	}
	void async_particle(double *values_x,double *values_t,double *values_h,int nvalues,
		int pso_n,double pso_o,double pso_fi_p,double pso_fi_g,double pso_eps,int pso_max_iter,
		double **particles,int i,double *best_p,int *total_niter,int size,double u)
	{
		int iter=0;
		double max=0;
		double o_init=pso_o;
		do
		{
			update_particle(size,particles,i,best_p,pso_n,pso_o,pso_fi_p,pso_fi_g,pso_eps,iter,u);
			particles[i][0]=calc_err(particles[i],values_x,values_t,values_h,nvalues,size,0);
#pragma omp critical
			max=update_global(size,particles,best_p,pso_n,total_niter[0],1);
			if ((((max - best_p[0]) < pso_eps)&&(iter>(pso_max_iter/2)))||(best_p[0]<pso_eps))
				break;
#pragma omp critical
			total_niter[0]++;
			if (total_niter[0]>pso_max_iter*pso_n)
				break;
			iter++;
			// C-PSO
			if (restart_prox<0) pso_o=o_init*sin(M_PI*pso_o);
		}
		while (1);
		printf("particle %d thread %d niters %d\n",i,omp_get_thread_num(),iter);
	}
	// main fitting routine
	// - values_x, values_t, values_h - arrays of size nvalues of experimental data
	T* fit_and_solve(double *values_x,double *values_t,double *values_h,int nvalues,
		int pso_n,double pso_o,double pso_fi_p,double pso_fi_g,double pso_eps,int pso_max_iter,double maxT)
	{
		int size=n_fit_parameters;
		double **particles;
		int best;
		double *best_p;
		int iter=0;
		// do initial randomize
		srand(time(NULL)); 
		// output
		printf("nvalues %d\n",nvalues);
		printf("n %d omega %g fi_p %g fi_g %g eps %g max_iter %d restart %g\n",pso_n, pso_o, pso_fi_p, pso_fi_g, pso_eps, pso_max_iter,restart_prox);
		printf("fit parameters:\n");
		for (int i=0;i<n_fit_parameters;i++)
			printf("%d min %g max %g default %g mask %d\n",i,bounds_fit[i][0],bounds_fit[i][1],bounds_fit[i][2],(mask_fit>>i)&1);
		// alloc
		best_p=new double[2*size+1];
		particles=new double *[pso_n];
		for (int i=0;i<pso_n;i++)
			particles[i]=new double[3*size+2]; // particles[0] contains f value, particles[1+size].. contain velocities, particles[1+2*size]... contain particle's best known position, particles[1+3*size] contains best known f value
		// initialize
		for (int i=0;i<pso_n;i++)
		{
			init_particle(particles[i]);
			for (int j=0;j<size;j++)
			{
				particles[i][1+size+j]=0;
				particles[i][1+2*size+j]=particles[i][1+j];
			}
		}
		solve_all_and_test(particles,pso_n,values_x,values_t,values_h,nvalues,size);
		for (int i=0;i<pso_n;i++)
			particles[i][1+3*size]=particles[i][0];
		best=0;
		double u=particles[0][0]; //average initial f
		for (int i=1;i<pso_n;i++)
		{
			if (particles[i][0]<particles[best][0])
				best=i;
			u+=particles[i][0];
		}
		u/=pso_n;
		// save best known position
		for (int j=0;j<=size;j++)
			best_p[j]=particles[best][j];
		printf("initial best: ");
		for (int j = 0;j <= size;j++)
			printf("%2.2g ", best_p[j]);
		printf("\n");
		fflush(stdout);
		fprintf(stderr,"0/%d\n",pso_max_iter);
		// process
		if (pso_max_iter>=1)
		{
			int total_niter=0;
#pragma omp parallel for num_threads(pso_n)
			for (int i=0;i<pso_n;i++)
				async_particle(values_x,values_t,values_h,nvalues,pso_n,pso_o,pso_fi_p,pso_fi_g,pso_eps,pso_max_iter,
						particles,i,best_p,&total_niter,size,u);
		}
		// solve with best parameters values
		T *solv;
		double err=calc_err(best_p,values_x,values_t,values_h,nvalues,size,1,maxT,&solv);
		printf("avg abs err %g param values ",err);
		FILE *f1=fopen("inverse_values.txt","wt");
		for (int i=0;i<n_fit_parameters;i++)
		{
			printf("%g ",best_p[i+1]);
			fprintf(f1,"%g ",best_p[i+1]);
		}
		fclose(f1);
		printf("\n");
		// free
		delete [] best_p;
		for (int i=0;i<pso_n;i++)
			delete [] particles[i];
		delete [] particles;
		return solv;
	}
};
// reads experimental data from fn into values_x,values_t, values_h,nvalues. Allocates needed memory
int read_dat_file (char *fn,double *&values_x,double *&values_t,double *&values_h,int &nvalues,int div_x,int div_t,double noise)
{
	printf("data file - %s\n",fn);
	if (FILE *fi=fopen(fn,"rt"))
	{
		char str[100000];
		char *dummy,*dummy2;
		std::vector<double> xvals;
		nvalues=0;
		while (fgets(str,100000,fi)) nvalues++;
		fseek(fi,0,SEEK_SET);
		nvalues--; // first line - x values
		nvalues/=div_t; // take only each div_t row
		printf("n rows %d\n",nvalues);
		// parse first line with x values
		dummy=fgets(str,100000,fi);
		while (dummy[0]!=':') dummy++;
		dummy+=2;
		dummy2=dummy;
		int id=1;
		printf("x values: ");
		while (dummy[0])
		{
			char s;
			while (dummy[0]&&(dummy[0]!=' ')) dummy++;
			s=dummy[0];
			dummy[0]=0;
			double v=atof(dummy2);
			dummy[0]=s;
			dummy2=dummy;
			if ((id%div_x)==0)
			{
				xvals.push_back(v);
				printf("%g ",v);
			}
			if (s==0)
			    break;
			dummy++;
			dummy2++;
			id++;
		}
		nvalues*=xvals.size();
		printf("\ntotal values %d\n",nvalues);
		if (nvalues<=0) return 0;
		// alloc memory
		values_x=new double[nvalues];
		values_t=new double[nvalues];
		values_h=new double[nvalues];
		// parse data
		nvalues=0;
		int row=1;
		while (dummy=fgets(str,100000,fi))
		if (((row++)%div_t)==0)
		{
			while (dummy[0]!=' ') dummy++;
			dummy++;
			dummy2=dummy;
			while (dummy[0]!=' ') dummy++;
			dummy[0]=0;
			double T=atof(dummy2);
			dummy[0]=' ';
			while (dummy[0]!=':') dummy++;
			dummy+=2;
			dummy2=dummy;
			int id=1;
			int xv=0;
			while (dummy[0])
			{
				char s;
				while (dummy[0]&&(dummy[0]!=' ')) dummy++;
				s=dummy[0];
				dummy[0]=0;
				double v=atof(dummy2);
				dummy[0]=s;
				dummy2=dummy;
				if ((id%div_x)==0)
				{
					values_x[nvalues]=xvals[xv++];
					values_t[nvalues]=T;
					values_h[nvalues]=v+noise*(rnd()-0.5);
					printf("%d %g %g %g\n",nvalues,values_t[nvalues],values_x[nvalues],values_h[nvalues]);
					nvalues++;
				}
				if (s==0)
				    break;
				dummy++;
				dummy2++;
				id++;
			}
		}
		fclose(fi);
		return 1;
	}
	else
	{
		printf("data file not found\n");
		exit(1);
	}
	exit(0);
	return 0;
}
// runs fitting for experimental data from fn
// - pso parameters in n,o,fi_p,fi_g
// - iteration procedure fininsh condition parameters - eps, max_iter
// - rp - probability of particle "restart"
// - nfit - number of parameters
// - mask - bit mask of the parameters to fit
// - bv[nfit][3] - array of [min,max,default] values of parameters
template<class T> T* run_inverse(char *fn,int pso_n,double pso_o,double pso_fi_p,double pso_fi_g,double pso_eps,int pso_max_iter,double rp,int nfit,int mask,double bv[][3],int div_x,int div_t,double noise,double maxT)
{
	double *vx,*vt,*vh;
	int nv;
	if (read_dat_file(fn,vx,vt,vh,nv,div_x,div_t,noise))
	{
		PSO_solver<T> *ps= new PSO_solver<T>(rp,nfit,mask,bv);
		T *s=ps->fit_and_solve(vx,vt,vh,nv,pso_n, pso_o, pso_fi_p, pso_fi_g, pso_eps, pso_max_iter,maxT);
		return s;
	}
	return NULL;
}
