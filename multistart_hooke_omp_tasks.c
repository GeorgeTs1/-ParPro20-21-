#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>

#define MAXVARS		(250)	/* max # of variables	     */
#define RHO_BEGIN	(0.5)	/* stepsize geometric shrink */
#define EPSMIN		(1E-6)	/* ending value of stepsize  */
#define IMAX		(5000)	/* max # of iterations	     */
#define NUM_THREADS     (3)

/* global variables */
unsigned long funevals = 0;


/* Rosenbrocks classic parabolic valley ("banana") function */
double f(double *x, int n,int id)
{
    double fv;
    int i;
	
	funevals++;
    fv = 0.0;
    //#pragma omp parallel for firstprivate(n,x,id) shared(fv)
    for (i=0; i<n-1; i++)   /* rosenbrock */
        fv = fv + 100.0* pow((x[i+1]-x[i]*x[i]),2) + pow((x[i]-1.0),2);
   
   // #pragma omp  nowait

    return fv;
}

/* given a point, look for a better one nearby, one coord at a time */
double best_nearby(double delta[MAXVARS], double point[MAXVARS], double prevbest, int nvars,int id)
{
	double z[MAXVARS];
	double minf, ftmp;
	int i;
	minf = prevbest;

	#pragma omp parallel 
	{
	
	#pragma omp task firstprivate(delta,prevbest,nvars,id) shared(minf,point) private(ftmp,z)
	{
	for (i = 0; i < nvars; i++)
		//#pragma omp task firstprivate(i,point) shared(z)
		z[i] = point[i];
	
	#pragma omp nowait
	
	//#pragma omp  for
	for (i = 0; i < nvars; i++) {
		//#pragma omp task firstprivate(i,point,nvars,id,delta) private(ftmp) shared(minf,z)
		//{ 
		z[i] = point[i] + delta[i];
		ftmp = f(z, nvars,id);
		if (ftmp < minf)
			minf = ftmp;
		else {
			delta[i] = 0.0 - delta[i];
			z[i] = point[i] + delta[i];
			ftmp = f(z, nvars,id);
			if (ftmp < minf)
				minf = ftmp;
			else
				z[i] = point[i];
		}
	}
	
		//}

	#pragma omp nowait

	//#pragma omp  for
	for (i = 0; i < nvars; i++)
		//#pragma omp task firstprivate(i,z) shared(point)
		point[i] = z[i];
    	#pragma omp nowait
}
	
}
	return (minf);
}


int hooke(int nvars, double startpt[MAXVARS], double endpt[MAXVARS], double rho, double epsilon, int itermax,int id)
{
	double delta[MAXVARS];
	double newf, fbefore, steplength, tmp;
	double xbefore[MAXVARS], newx[MAXVARS];
	int i, j, keep;
	int iters, iadj;

	#pragma omp parallel 
	{
	#pragma omp task firstprivate(nvars,startpt,rho,epsilon,itermax,id,endpt) shared(iters) private(delta,newf,fbefore,steplength,tmp,xbefore,newx,keep,iadj) 
	{
	for (i = 0; i < nvars; i++) {
		newx[i] = xbefore[i] = startpt[i];
		delta[i] = fabs(startpt[i] * rho);
		if (delta[i] == 0.0)
			delta[i] = rho;
	}
	
        #pragma omp nowait
	
	iadj = 0;
	steplength = rho;
	//#pragma omp critical
	iters = 0;
	fbefore = f(newx, nvars,id);
	newf = fbefore;


	while ((iters < itermax) && (steplength > epsilon)) {
		//#pragma omp atomic
		iters++;
		iadj++;
#if DEBUG
		printf("\nAfter %5d funevals, f(x) =  %.4le at\n", funevals, fbefore);
		//#pragma omp  for
		for (j = 0; j < nvars; j++)
			printf("   x[%2d] = %.4le\n", j, xbefore[j]);
#endif
		//#pragma omp nowait

		/* find best new point, one coord at a time */
		//#pragma omp  for
		for (i = 0; i < nvars; i++) {
			newx[i] = xbefore[i];
		}

		#pragma omp nowait
		newf = best_nearby(delta, newx, fbefore, nvars,id);

		//#pragma omp barrier
		/* if we made some improvements, pursue that direction */
		keep = 1;

		

		while ((newf < fbefore) && (keep == 1)) {
			iadj = 0;
			//#pragma omp for
			for (i = 0; i < nvars; i++) {
				/* firstly, arrange the sign of delta[] */
				if (newx[i] <= xbefore[i])
					delta[i] = 0.0 - fabs(delta[i]);
				else
					delta[i] = fabs(delta[i]);
				/* now, move further in this direction */
				tmp = xbefore[i];
				xbefore[i] = newx[i];
				newx[i] = newx[i] + newx[i] - tmp;
			}

			#pragma omp nowait

			fbefore = newf;
			newf = best_nearby(delta, newx, fbefore, nvars,id);
			/* if the further (optimistic) move was bad.... */
			if (newf >= fbefore)
				break;

			/* make sure that the differences between the new */
			/* and the old points are due to actual */
			/* displacements; beware of roundoff errors that */
			/* might cause newf < fbefore */
			keep = 0;
			for (i = 0; i < nvars; i++) {
				keep = 1;
				if (fabs(newx[i] - xbefore[i]) > (0.5 * fabs(delta[i])))
					break;
				else
					keep = 0;
			}
		}
		if ((steplength >= epsilon) && (newf >= fbefore)) {
			steplength = steplength * rho;
			//#pragma omp  for
			for (i = 0; i < nvars; i++) {
				delta[i] *= rho;
			}
			#pragma omp nowait
		}
	


	}
	//#pragma omp  for 
	for (i = 0; i < nvars; i++)
		endpt[i] = xbefore[i];

	#pragma omp nowait
}	

}

	return (iters);
}


double get_wtime(void)
{
    struct timeval t;

    gettimeofday(&t, NULL);

    return (double)t.tv_sec + (double)t.tv_usec*1.0e-6;
}

int main(int argc, char *argv[])
{
	double startpt[MAXVARS], endpt[MAXVARS];
	int itermax = IMAX;
	double rho = RHO_BEGIN;
	double epsilon = EPSMIN;
	int nvars;
	int trial, ntrials;
	double fx;
	int i, jj;
	double t0, t1;

	double best_fx = 1e10;
	double best_pt[MAXVARS];
	int best_trial = -1;
	int best_jj = -1;
	int id;
        int Nthrds;
        int istart; 
        int iend;  



	for (i = 0; i < MAXVARS; i++) best_pt[i] = 0.0;

	ntrials = 128*1024;	/* number of trials */
	nvars = 16;		/* number of variables (problem dimension) */
	srand48(1);

	t0 = get_wtime();

	omp_set_num_threads(NUM_THREADS);

	#pragma omp parallel firstprivate(ntrials,nvars)  private(id)
	{
		    id = omp_get_thread_num();
		    Nthrds = omp_get_num_threads();
		    istart = id * ntrials / Nthrds; 
		    iend = (id+1) * ntrials / Nthrds; 
		    if (id == omp_get_num_threads()-1) iend = ntrials;
		    //#pragma omp set_dynamic(0)
		   
			
	          // #pragma omp taskwait


		

		#pragma omp single nowait
		{
		#pragma omp task private(jj,fx,startpt,endpt) firstprivate(ntrials,nvars,rho,epsilon,itermax,id,istart,iend) shared(best_fx,best_trial,best_jj,best_pt) 
		{
		
		for (trial = istart; trial < iend; trial++) {
		/* starting guess for rosenbrock test function, search space in [-4, 4) */
		

		//#pragma omp single nowait
		//{
		for (i = 0; i < nvars; i++) {
			#pragma omp task firstprivate(i) shared(startpt)
			startpt[i] = 4.0*drand48()-4.0;
		//}
		}
	
		//#pragma omp taskwait

		jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax,id);
#if DEBUG
		printf("\n\n\nHOOKE %d USED %d ITERATIONS, AND RETURNED\n", trial, jj);

		for (i = 0; i < nvars; i++)
			printf("x[%3d] = %15.7le \n", i, endpt[i]);
		//#pragma omp nowait
#endif

		
		fx = f(endpt, nvars,id);
#if DEBUG
		printf("f(x) = %15.7le\n", fx);
#endif

			
		//#pragma omp taskwait

		#pragma omp critical
		{
		if (fx < best_fx) {
			best_trial = trial;
			best_jj = jj;
			best_fx = fx;
			for (i = 0; i < nvars; i++)
				best_pt[i] = endpt[i];  
		
		}

		}
   	
}


}

}

	#pragma omp taskwait

}

	t1 = get_wtime();

	printf("\n\nFINAL RESULTS:\n");
	printf("Elapsed time = %.3lf s\n", t1-t0);
	printf("Total number of trials = %d\n", ntrials);
	printf("Total number of function evaluations = %ld\n", funevals);
	printf("Best result at trial %d used %d iterations, and returned\n", best_trial, best_jj);
	for (i = 0; i < nvars; i++) {
		printf("x[%3d] = %15.7le \n", i, best_pt[i]);
	}
	printf("f(x) = %15.7le\n", best_fx);

	return 0;
}
