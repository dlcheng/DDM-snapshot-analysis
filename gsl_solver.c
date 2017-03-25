#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_errno.h>


double cNFW_from_vel(double v_ratio)
{
  /* c=2.15 is the minimum value, the method can get */
  int status;
  int iter = 0;     
  const gsl_root_fdfsolver_type *T;
  gsl_root_fdfsolver *s;
  double x0;
  double x = 100000;              /* initial guess, it is not possible to have such high c*/
  gsl_function_fdf FDF;
  
  FDF.f = &function_c;
  FDF.df = &d_function_c;
  FDF.fdf = &fd_function_c;
  FDF.params = &v_ratio;

  T = gsl_root_fdfsolver_newton;
  s = gsl_root_fdfsolver_alloc(T);
  gsl_root_fdfsolver_set(s, &FDF, x);

  do{
    iter++;
    status = gsl_root_fdfsolver_iterate(s);
    x0 = x;
    x = gsl_root_fdfsolver_root(s);
    status = gsl_root_test_delta(x, x0, 0, 1e-3);
    }
    while(status == GSL_CONTINUE && iter < 200);

  gsl_root_fdfsolver_free(s);

  return x;	
	
}   /* end cNFW_from_vel */		


double function_c(double c, void *param)
{
  double Rv = *(double *) param;	

  double result;
  
  result = c * (1 + c) - Rv * Rv / 0.216217 * ((1 + c) * log(1 + c) - c);
 
  return result;	
}	    /* end function_c */

double d_function_c(double c, void *param)
{
  double Rv = *(double *) param;	
  double result;
  
  result = 2 * c + 1 - Rv * Rv / 0.216217 * log(1 + c);
  
  return result;
}	/* end d_function_c */

void fd_function_c(double c, void *param, double *f, double *df)
{
  *f = function_c(c, param);
  *df = d_function_c(c, param);	
}    /* end fd_function_c */	
