/*  URMAE/reconstruction.GM/V1/syslin.LU.c                                    */
/*  Mennessier Gerard                 20060321                                */
/*  Last Revised : G.M.               20060325                                */
/*                                                                            */

#include  <stddef.h>
#include  <stdio.h>

#include  "utiTnsr.h"

#include  "recipe.lubksb.h"
#include  "recipe.ludcmp.h"

#include  "dataPoints.glob.h"
#include  "syslin.LU.h"
#include  "funcBasis.h"

/******************************************************************************/
/*                                                                            */
/* LU                                                                         */
/*                                                                            */
/******************************************************************************/
void      sl_LU(size_t n, double *bp, double *xp)
{ tnsr2dbl  sat;
  double  **sapp, *sap;
  double    det;
  int      *indexsavp;
  int       i, j;
  double   *currentptp;
  double   *currentbasisp;


  indexsavp = malloc(n*sizeof(int));
fprintf(stderr, "syslin.LU::sl_LU     indexsavp ALLOCATED \n");
  tnsr2dblPAlloc(&sat, n, n);
fprintf(stderr, "syslin.LU::sl_LU     sat ALLOCATED \n");

  sapp = sat.pp;
  currentptp = ptsCentresVec.p;
  for(i=0;  i < n;  i++, currentptp += 3, sapp++)
  { sap = *sapp;
    currentbasisp = ptsCentresVec.p; 
    for(j=0;  j < n;  j++, currentbasisp += 3, sap++)
    { *sap = fbasis(currentptp, currentbasisp);
    }
  }
  sat.ev1x = sat.ev2x = n;
fprintf(stderr, "syslin.LU::sl_LU    sat COMPUTED \n");
tnsr2dblPrint(stderr, &sat);
/* tnsr2dblPPrint(stderr, &sat); */
fprintf(stderr, "syslin.LU::sl_LU    sat PRINTED \n");


  sapp = sat.pp;
  ludcmp(sapp, n, indexsavp, &det);                           /** LU tiangulation **/
fprintf(stderr, "syslin.LU::sl_LU    LU triangulation DONE \n");


  lubksb(sapp, n, indexsavp, bp);                                     /** SOLVE **/
  for(i=0;  i < n;  i++){ xp[i] = bp[i];}
  free(indexsavp);
  free(sat.p); free(sat.pp);
  return;
}
/******************************************************************************/
/******************************************************************************/
