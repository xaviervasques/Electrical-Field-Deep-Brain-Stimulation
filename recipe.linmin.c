/*  libc/recipe2double/recipe.linmin.c                                        */
/** MATRIX INDICES FROM 0 to DIM-1                                           **/
/** MATRIX : DOUBLE version                                                  **/

#include <math.h>
#include "recipe.nrutil.h"
#include "utiAlloc.h"
#include "utistdIO.h"

#include "recipe.linmin.h"
#include "recipe.linmin.statalloc.h"
#include "recipe.mnbrak.h"
#include "recipe.brent.h"

/** #define TOL 2.0e-12 **/

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      linmin(double p[], double xi[], int n, double *fretp, double tol)
{ int       j;
  double    xx,bx,ax, xmin;
  double    fx,fb,fa, fret;

  for(j = 0;  j < n;  j++){ pcom[j] = p[j];  xicom[j] = xi[j];}
  ax = 0.0;
  xx = 1.0;
  mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, linminf1dim);
/*
fPrintF(stderr,"linmin::linmin  after mnbrak  x's=%f,%f,%f  f's=%f,%f,%f \n",
                                                                  ax,xx,bx,fa,fx,fb);
*/
  fret = brent(ax, xx, bx, linminf1dim, tol, &xmin);
/*
fPrintF(stderr,"linmin::linmin  after brent   xmin=%f fret=%f \n", xmin,fret);
*/
  for(j = 0;  j < n;  j++){ xi[j] *= xmin;  p[j] += xi[j];}
  *fretp = fret;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double    linminf1dim(double x)
{ int       j;
  double    f;

  for(j = 0;  j < ncom;  j++){ xtcom[j] = pcom[j] + x * xicom[j];}
  f = (*nrfunc)(xtcom);
  return f;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      linminAlloc(int n, double (*func)(double []))
{ 
  ncom = n;
  nrfunc = func;
  pcom  = dblAlloc((size_t)n, "recipe.linmin::linminAlloc  pcom");
  xicom = dblAlloc((size_t)n, "recipe.linmin::linminAlloc  xicom");
  xtcom = dblAlloc((size_t)n, "recipe.linmin::linminAlloc  xtcom");
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      linminFree()
{ free(xtcom);  free(xicom);  free(pcom);
  return;
}
/******************************************************************************/
/******************************************************************************/
#undef TOL
/* (C) Copr. 1986-92 Numerical Recipes Software 'i<. */
