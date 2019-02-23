/*  libc/recipe2double/recipe.amotry.c                                        */
/*       Revised : G.M.               20030522                                */
/*  Last Revised : G.M.               20051028                                */

/** VECTOR and MATRIX INDICES FROM 0 to DIM-1                                **/
/** VECTOR and MATRIX : DOUBLE version                                       **/

#include  "utistdIO.h"
#include  "utiAlloc.h"
#include  "recipe.amotry.h"

/******************************************************************************/
double    amotry(double **p, double y[], double psum[], int ndim,
                                      double (*funk)(double []), int ihi, double fac)
{ int      j;
  double   fac1, fac2, ytry, *ptry;
/*  static char   prognamp[] = "recipe.amotry::amotry"; */

  ptry = dblAlloc((size_t)ndim, "recipe.amotry");
  fac1 = (1.0 - fac)/ndim;
  fac2 = fac1 - fac;
  for(j = 0;  j < ndim;  j++){ ptry[j] = psum[j]*fac1 - p[ihi][j]*fac2;}

/*
fPrintF(stderr, "%s  ihi=%d,  fac=%f,fac1=%f,fac2=%f \n", prognamp, ihi, fac, fac1, fac2);
for(j = 0;  j < ndim;  j++) fPrintF(stderr," %10.3e", psum[j]);
fPrintF(stderr, "\n");
*/

  ytry = (*funk)(ptry);
/* fPrintF(stderr, "%s  ihi=%d,  fac=%f, ytry=%f \n", prognamp, ihi, fac, ytry); */
  if (ytry < y[ihi])
  { y[ihi] = ytry;
    for(j = 0;  j < ndim;  j++)
    { psum[j]  += ptry[j] - p[ihi][j];
      p[ihi][j] = ptry[j];
    }
  }
  free(ptry);
  return ytry;
}
/******************************************************************************/
/******************************************************************************/
/* (C) Copr. 1986-92 Numerical Recipes Software 'i<. */
/******************************************************************************/
