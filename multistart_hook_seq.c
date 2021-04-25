#include <stdio.h>
#include <math.h>
#include<time.h>
#include<omp.h>

#define      VARS		(250)	/* max # of variables	     */
#define      RHO_BEGIN		(0.5)	/* stepsize geometric shrink */
#define      EPSMIN		(1E-6)	/* ending value of stepsize  */
#define      IMAX		(5000)	/* max # of iterations	     */
#define      NTRAILS    (10) /*max number of points     */

/* global variables */
/* global variables */
int		funevals = 0;

#ifdef Woods
double f();
#else


/* Rosenbrocks classic parabolic valley ("banana") function */
double
f(x, n)
	   double	   x[VARS];
	   int		   n;
{
	   double	   a, b, c;
	   funevals++;
	   a = x[0];
	   b = x[1];
	   c = 100.0 * (b - (a * a)) * (b - (a * a));
	   return (c + ((1.0 - a) * (1.0 - a)));
}

#endif


/* given a point, look for a better one nearby, one coord at a time */
double
best_nearby(delta, point, prevbest, nvars)
	   double	   delta[VARS], point[VARS];
	   double	   prevbest;
	   int		   nvars;
{
	   double	   z[VARS];
	   double	   minf, ftmp;
	   int		   i;
	   minf = prevbest;
	   for (i = 0; i < nvars; i++)
		   z[i] = point[i];
	   for (i = 0; i < nvars; i++) {
		   z[i] = point[i] + delta[i];
		   
		   if(z[i]<-4 || z[i]>4)
		   {
		   		z[i]=point[i];
		   }
		
		   ftmp = f(z, nvars);
		   if (ftmp < minf)
			   minf = ftmp;
		   else {
			   delta[i] = 0.0 - delta[i];
			   z[i] = point[i] + delta[i];
			   ftmp = f(z, nvars);
			   if (ftmp < minf)
				   minf = ftmp;
			   else
				   z[i] = point[i];
		   }
	   }
	   for (i = 0; i < nvars; i++)
		   point[i] = z[i];
	   return (minf);
}


int
hooke(nvars, startpt, endpt, rho, epsilon, itermax)
	   double	   startpt[VARS], endpt[VARS];
	   int		   nvars, itermax;
	   double	   rho, epsilon;
{
	   double	   delta[VARS];
	   double	   newf, fbefore, steplength, tmp;
	   double	   xbefore[VARS], newx[VARS];
	   int		   i, j, keep;
	   int		   iters, iadj;
	   for (i = 0; i < nvars; i++) {
		   newx[i] = xbefore[i] = startpt[i];
		   delta[i] = fabs(startpt[i] * rho);
		   if (delta[i] == 0.0)
			   delta[i] = rho;
	   }
	   iadj = 0;
	   steplength = rho;
	   iters = 0;
	   fbefore = f(newx, nvars);
	   newf = fbefore;
	   while ((iters < itermax) && (steplength > epsilon)) {
		   iters++;
		   iadj++;
		   printf("\nAfter %5d funevals, f(x) =  %.4le at\n", funevals, fbefore);
		   for (j = 0; j < nvars; j++)
			   printf("   x[%2d] = %.4le\n", j, xbefore[j]);
		   /* find best new point, one coord at a time */
		   for (i = 0; i < nvars; i++) {
			   newx[i] = xbefore[i];
		   }
		   newf = best_nearby(delta, newx, fbefore, nvars);
		   /* if we made some improvements, pursue that direction */
		   keep = 1;
		   while ((newf < fbefore) && (keep == 1)) {
			   iadj = 0;
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
			   fbefore = newf;
			   newf = best_nearby(delta, newx, fbefore, nvars);
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
				   if (fabs(newx[i] - xbefore[i]) >
				       (0.5 * fabs(delta[i])))
					   break;
				   else
					   keep = 0;
			   }
		   }
		   if ((steplength >= epsilon) && (newf >= fbefore)) {
			   steplength = steplength * rho;
			   for (i = 0; i < nvars; i++) {
				   delta[i] *= rho;
			   }
		   }
	   }
	   for (i = 0; i < nvars; i++)
		   endpt[i] = xbefore[i];
	   return (iters);
}

#ifndef Woods
main()
{
	   int id;
	   double	   startpt[4][VARS],endpt[4][VARS],tmp_min_f,curr_f;
	   int		   nvars, itermax;
	   double	   rho, epsilon;
	   int		   i, jj;
	   int         point,j;
	   double      tmp_min_x[VARS];
	   int 		   final_point;
	   nvars=16;
	 
	   /* starting guess for rosenbrock test function */
	   
	   
	   itermax = IMAX;
	   rho = RHO_BEGIN;
	   epsilon = EPSMIN;
	   
	   srand(1);
	   
	   omp_set_num_threads(4);
	
	   	#pragma omp parallel private(i) shared(j) 
	   	{
			
		 #pragma omp for
	      for(i=0; i<NTRAILS; i++)
			{	
				id =  omp_get_thread_num();
				
				#pragma omp for 
			    for(j=0; j<nvars; j++)
			      startpt[id][j]=rand();
			
		        jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax);
				
		    	printf("\nThread with id = %d finished the method\n",id);
		    
				 	
				if(i==0) 
				{
					tmp_min_f = f(endpt[id],nvars);
					for(j=0; j<nvars; j++)
					{
						tmp_min_x[j]=endpt[id][j];
					}		
				}
			
				else
				{	
					curr_f = f(endpt[id],nvars);
					if(tmp_min_f>curr_f)
					{
						tmp_min_f = curr_f;
						for(j=0; j<nvars; j++)
						{
							tmp_min_x[j]=endpt[id][j];	
						}
			 		}
	  		    }
		    }
		}
	   	
	   
	   for (i = 0; i < nvars; i++)
		   printf("x[%3d] = %15.7le \n", i,tmp_min_x[i]);

	 	
	  
	   printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj);
}
	   
#else
/* The Hooke & Jeeves algorithm works reasonably well on
 * Rosenbrock's function, but can fare worse on some
 * standard test functions, depending on rho.  Here is an
 * example that works well when rho = 0.5, but fares poorly
 * with rho = 0.6, and better again with rho = 0.8.
 */

#ifndef RHO_WOODS
#define RHO_WOODS 0.6
#endif

/* Woods -- a la More, Garbow & Hillstrom (TOMS algorithm 566) */

/*double
f(x, n)
	double	x[VARS];
	int	n;
{
	double	s1, s2, s3, t1, t2, t3, t4, t5;
	funevals++;
	s1 = x[1] - x[0]*x[0];
	s2 = 1 - x[0];
	s3 = x[1] - 1;
	t1 = x[3] - x[2]*x[2];
	t2 = 1 - x[2];
	t3 = x[3] - 1;
	t4 = s3 + t3;
	t5 = s3 - t3;
	return 100*(s1*s1) + s2*s2 + 90*(t1*t1) + t2*t2
		+ 10*(t4*t4) + t5*t5/10.;
}

main()
{
	double	startpt[VARS], endpt[VARS];
	int	nvars, itermax;
	double	rho, epsilon;
	int	i, jj;

	/* starting guess test problem "Woods" */
	/*nvars = 4;
	startpt[0] = -3;
	startpt[1] = -1;
	startpt[2] = -3;
	startpt[3] = -1;

	itermax = IMAX;
	rho = RHO_WOODS;
	epsilon = EPSMIN;
	jj = hooke(nvars, startpt, endpt, rho, epsilon, itermax);
	printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj);
	for (i = 0; i < nvars; i++)
		printf("x[%3d] = %15.7le \n", i, endpt[i]);
	printf("True answer: f(1, 1, 1, 1) = 0.\n");
}*/
#endif  
