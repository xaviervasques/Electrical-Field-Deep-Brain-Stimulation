/*  libc/recipe2double/recipe.frprmn.c                                        */
/** MATRIX INDICES FROM 0 to DIM-1                                           **/
/** MATRIX : DOUBLE version                                                  **/

/**  Fletcher Reeves Polak Ribiere minimization **/

#include <math.h>
#include "recipe.nrutil.h"

#include "utiAlloc.h"
#include "utistdIO.h"

#include "recipe.frprmn.h"
#include "recipe.linmin.h"

#define  EPS     1.0e-10
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      frprmn(double p[], int n, double ftol, int *iter, double *fretp,
           double (*func)(double []), void (*dfunc)(double [], double []), int itmax)
{ int       j, it;
  double   *xi, fc, fret;
  double    gg, gam, dgg;
  double   *g, *h;
  double    lintol = 2.0e-12 ;
  char      converg;
  static    char     prognamp[] = "recipe.frprmn::frprmn";
  static    char     form1p[] = "%s  Initial Function value = %g \n";
  static    char     form2p[] = "%s"
                         "  Initial First Gradient values = %g, %g, ..., %g ,%g \n";

  g  = dblAlloc((size_t)n, "recipe.frprmn::frprmn  g");
  h  = dblAlloc((size_t)n, "recipe.frprmn::frprmn  h");
  xi = dblAlloc((size_t)n, "recipe.frprmn::frprmn  xi");
  linminAlloc(n, func);

fPrintF(stderr,"%s  pointer p=%p  \n", prognamp, (void*)p);
  fc = (*func)(p);          /** Function **/
fPrintF(stderr,form1p, prognamp,fc);
  (*dfunc)(p,xi);           /** Gradient **/
fPrintF(stderr,form2p, prognamp, xi[0], xi[1], xi[n-2], xi[n-1]);

  fret = fc;
  for (j = 0;  j <n;  j++){ g[j] = -xi[j];  xi[j] = h[j] = g[j];}
  converg = 0;
  for (it = 1;  it <= itmax;  it++)
  { *iter = it;
    linmin(p, xi, n, &fret, lintol);
    if(it%100 == 0)
    { fPrintF(stderr,"frprmn::  it = %d;  after linmin best f= %f\n", it, fret);
    }
    if(2.0*fabs(fret-fc)  <=  ftol*(fabs(fret)+fabs(fc) + EPS)) 
    { converg = 1;  break;                       /** Normal RETURN if convergence **/
    }
    fc = (*func)(p);
    (*dfunc)(p, xi);
    dgg = gg = 0.0;
    for(j = 0;  j < n;  j++) { gg += g[j]*g[j];  dgg += (xi[j]+g[j])*xi[j];}
    if(gg == 0.0)
    { converg = 1;  break;
    }
    gam = dgg/gg;
    for (j = 0;  j < n;  j++){ g[j] = -xi[j];  xi[j] = h[j] = g[j]+gam*h[j];}
  }
                                      /* nrerror("Too many iterations in frprmn"); */
  fPrintF(stderr,"%s  pointer p=%p \n", prognamp, (void*)p);
  linminFree();  free(xi);  free(h);  free(g);
  *fretp = fret;
  fPrintF(stderr,"%s  pointer p=%p\n", prognamp, (void*)p);
  return;
}
/******************************************************************************/
/******************************************************************************/
#undef EPS
/* (C) Copr. 1986-92 Numerical Recipes Software 'i<. */
