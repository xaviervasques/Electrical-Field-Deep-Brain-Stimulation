/*  ../myGauss/syslin.GaussSeidel.c                                           */
/*                                                                            */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include  "dataPoints.glob.h"

#include  "syslin.GaussSeidel.h"
#include  "funcBasis.h"

/******************************************************************************/
/*                                                                            */
/* sl_gauss_seidel                                                            */
/*                                                                            */
/******************************************************************************/
void      sl_gauss_seidel(size_t n, double *ap, double *bp, double *xp,
		                            int iterx, double **tpp, double eps)
{ int       i, j, k;
  double    alfa = 1.0;
  double    d, s;
  double   *currentptp;
  double   *currentbasisp;
  double   *xnewp;

  xnewp = malloc(n * sizeof(double));

  k = 0;

  while(k < iterx  &&  alfa > eps)
  { alfa = 0.0;
    currentptp = ptsCentresVec.p;
    for(i=0;  i < n;  i++, currentptp += 3)
    { d = xp[i];
      s = bp[i];
      currentbasisp = ptsCentresVec.p; 
      for(j=0;  j < n;  j++, currentbasisp += 3)
      { ap[j] = fbasis(currentptp, currentbasisp);

if(ap[j] != 0.0)
{ fprintf(stderr, "GaussSeidel i=%d j=%d, Pt=(%f,%f,%f) Centre=(%f,%f,%f) f=%f\n",
         i,j, currentptp[0], currentptp[1], currentptp[2],
              currentbasisp[0], currentbasisp[1], currentbasisp[2], ap[j]); 
}

        if(i!=j) s -= ap[j] * xp[j];
      }
      xnewp[i] = s/ap[i];
      d -= xnewp[i];
      alfa += d*d;
    }

    for(i=0;  i < n;  i++){ xp[i] = xnewp[i];}
    
fprintf(stderr,"GaussSeidel k=%d\n", k);
for(i=0;  i < n;  i++){ fprintf(stderr, " %f", xp[i]);}
fprintf(stderr,"\n");

    k++;
  }
fprintf(stdout, "GaussSeidel k=%d,  alfa =%12.4e, eps =%12.4e \n", k, alfa, eps);

  return;
}
/******************************************************************************/
/******************************************************************************/
