/*  libc/recipe2double/recipe.amoeba.c                                        */
/*  Last Revised : G.M.               20030523                                */

/** VECTOR and MATRIX INDICES FROM 0 to DIM-1                                **/
/** VECTOR and MATRIX : DOUBLE version                                       **/

#include <math.h>
#include  "recipe.nrutil.h"
#include  "utistdIO.h"
#include  "utiAlloc.h"
#include  "recipe.amoeba.h"
#include  "recipe.amotry.h"

#define GET_PSUM   \
                for (j = 0;  j < ndim;  j++)  \
                { for (sum =0.0, i = 0;  i < mpts;  i++){ sum += p[i][j];}  \
                  psum[j] = sum;  \
                }

#define SWAP(a,b) { swap=(a); (a)=(b); (b)=swap;}

void      simplexReduce(double (*funk)(double []), double *psum, double y[], 
                                                     double **p, int ilow, int ndim);
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      amoeba(double **p, double y[], int ndim, double ftol,
                                   double (*funk)(double []), int *nfunk, int nfunkx)
{ int       i, ihi, ilo, inhi, j, mpts;
 /** ilo=index of lowest, ihi=index of highest, inhi=index of previous to highest **/
  double    rtol, sum, swap, ysave, ytry, *psum;
  static char    prognamp[] = "recipe.amoeba::amoeba";

  mpts = ndim + 1;
  psum = dblAlloc(ndim, prognamp);
  *nfunk = 0;
  GET_PSUM
  for( ; ; ) 
  { ilo = 0;
                                          /** define initial values for ihi, inhi **/
    ihi = y[0] > y[1] ? (inhi = 1, 0) : (inhi = 0,1);
    for (i = 0;  i < mpts;  i++)
    { if (y[i] <= y[ilo]) ilo = i;
      if (y[i] > y[ihi])                   { inhi = ihi;  ihi = i;}
      else if (y[i] > y[inhi]  &&  i != ihi) inhi = i;
    }
/*
fPrintF(stderr, "amoeba 1 ilo=%d,ihi=%d,inhi=%d, ylo=%f,yhi=%f,ynhi=%f\n",
                                                ilo,ihi,inhi, y[ilo],y[ihi],y[inhi]);
*/
/*
fPrintF(stderr, "      ");
for(j = 0;  j < ndim;  j++) fPrintF(stderr," %10.3e", p[ilo][j]);
fPrintF(stderr, "\n");
*/
    rtol = 2.0 * fabs(y[ihi] - y[ilo])/(fabs(y[ihi]) + fabs(y[ilo]));
    if (rtol < ftol) break;
    if (*nfunk >= nfunkx){ fPrintF(stderr,"recipe.amoeba : nfunkx exceeded\n"); break;}
    
    ytry = amotry(p, y, psum, ndim, funk, ihi, -1.0);  (*nfunk)++;
/*
fPrintF(stderr, "amoeba 2 ilo=%d,ihi=%d,inhi=%d, ylo=%f,yhi=%f,ynhi=%f, ytry=%f\n",
                                          ilo,ihi,inhi, y[ilo],y[ihi],y[inhi], ytry);
*/
    if (ytry <= y[ilo])
    { 
/*
fPrintF(stderr, "amoeba 3 TRY with -1 did improved ytry=%f \n", ytry);
*/
      ytry = amotry(p, y, psum, ndim, funk, ihi, 2.0);  (*nfunk)++;
/*
fPrintF(stderr, "amoeba 4 ilo=%d,ihi=%d,inhi=%d, ylo=%f,yhi=%f,ynhi=%f, ytry=%f\n",
                                          ilo,ihi,inhi, y[ilo],y[ihi],y[inhi], ytry);
*/
    }
    else if (ytry >= y[inhi])
    { ysave = y[ihi];
      ytry = amotry(p, y, psum, ndim, funk, ihi, 0.5);  (*nfunk)++;
/*
fPrintF(stderr, "amoeba 5 ilo=%d,ihi=%d,inhi=%d, ylo=%f,yhi=%f,ynhi=%f, ytry=%f\n",
                                          ilo,ihi,inhi, y[ilo],y[ihi],y[inhi], ytry);
*/
            /** cannot improve. Expect we are in the hole ==> reduce simplex size **/
      if (ytry >= ysave)
      { 
/* fPrintF(stderr, "amoeba 6 NO IMPROVMENT reduce simplex SIZE \n"); */
        simplexReduce( funk, psum, y, p, ilo, ndim);
        *nfunk += ndim;
        GET_PSUM
      }
    }
  }

  if(ilo != 0)
  { SWAP(y[0], y[ilo])
    for (j = 0;  j < ndim;  j++) SWAP(p[0][j], p[ilo][j])
  }
  free(psum);
  return;
}
/******************************************************************************/
void      simplexReduce(double (*funk)(double []), double *psum, double y[], 
                                                      double **p, int ilow, int ndim)
{ int       i, j, mpts;

  mpts = ndim + 1;
  for (i = 0;  i < mpts;  i++)
  { if (i == ilow) continue;
    
    for (j = 0;  j < ndim;  j++)
    { p[i][j] = psum[j] = 0.5 * (p[i][j] + p[ilow][j]);
    }
    y[i] = (*funk)(psum);
  }
  return;
}
/******************************************************************************/
/******************************************************************************/
#undef SWAP
#undef GET_PSUM
/******************************************************************************/
/* (C) Copr. 1986-92 Numerical Recipes Software 'i<. */
/******************************************************************************/
