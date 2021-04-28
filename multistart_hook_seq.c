#include <stdio.h>
#include <math.h>
#include <time.h>

#define      VARS		(250)	/* max # of variables	     */
#define      RHO_BEGIN	(0.5)	/* stepsize geometric shrink */
#define      EPSMIN		(1E-6)	/* ending value of stepsize  */
#define      IMAX		(5000)	/* max # of iterations	     */

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
{	   int         m;
	   double	   a, b, c;
	   double      value=0;;
	   funevals++;
	   a = x[0];
	   b = x[1];
	   for(m=1; m<n-1; m++)
	   {
	   		c = 100.0 * (b - (a * a)) * (b - (a * a));
			value += 	(c + ((1.0 - a) * (1.0 - a)));
			a=x[m];
			b=x[m+1];   		
	   }
	   
	   return value;
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

		   /*if(z[i]<-4.0 || z[i]>4.0)
		   {
		   		z[i]= point[i] - delta[i];
		   		delta[i] = 0.0 - delta[i];
		   		z[i] = point[i] + delta[i];
		   		ftmp = f(z, nvars);
		   		if (ftmp < minf)
			        minf = ftmp;
		   	else {
			   z[i] =  + delta[i];
			   ftmp = f(z, nvars);
			   if (ftmp < minf)
				   minf = ftmp;
			   else
				   z[i] = point[i];
		   }*/
		   		
		   
		   
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
	   double	   startpt[VARS], endpt[VARS],tmp_min[VARS],prev_f,curr_f;
	   int		   nvars, itermax;
	   double	   rho, epsilon;
	   int		   i, j,jj;
	   int         k;
			
	   /* starting guess for rosenbrock test function */
	   nvars = 16;
		
		itermax = IMAX;
	   rho = RHO_BEGIN;
	   epsilon = EPSMIN;
		
		
		srand(1);
		

		for(i=0; i<10; i++)	
		{
			for(j=0; j<16; j++)
			{
				startpt[j]= rand()%8 -4;
				printf("\nStartpt[%d]=%lf\n",j,startpt[j]);
			}
				
			jj =hooke(nvars, startpt, endpt, rho, epsilon, itermax);
			
			if(i==0)
			{
				prev_f = f(endpt,nvars);
				
				printf("\nPREV_F=%lf\n",prev_f);
				
				for(k=0; k<16; k++)
				{
					tmp_min[k] = endpt[k];
					
				}
				
			}
			
			else
			{
				curr_f = f(endpt,nvars);
				
				for(k=0; k<16; k++)
					printf("   x[%2d] = %.4le\n", k, endpt[k]);
				
				if(curr_f < prev_f)
				{
					for(k=0; k<16; k++)
					{
					     tmp_min[k] = endpt[k];
					}
					
					prev_f = curr_f;
					
				}
				
				else if(curr_f == prev_f)
				{		
						prev_f=curr_f;
						continue;
				}	
			
			}
		}


	   
	   
	   printf("\n\n\nHOOKE USED %d ITERATIONS, AND RETURNED\n", jj);
	   for (i = 0; i < nvars; i++)
		   printf("x[%3d] = %15.7le \n f=%lf\n", i, tmp_min[i],prev_f);
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


#endif
