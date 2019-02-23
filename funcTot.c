#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include   "dataPoints.glob.h"
#include   "funcBasis.h"
#include   "funcTot.h"

/******************************************************************************/
/*                                                                            */
/*  currentxyzp ~                                                             */
/*  X=currentposition[0], Y=currentposition[1], Z=currentposition[2],         */
/*                                                                            */
/******************************************************************************/
double    ftot(double *currentxyzp)
{ int       i;
  double    f = 0.0;

  for (i=0;  i < NPTS ;  i++)
  { f  += weightVec.p[i] * fbasis(currentxyzp, ptsCentresVec.p + 3*i);
  }
  return f;
}
/******************************************************************************/
/*                                                                            */
/*  currentxyzp ~                                                             */
/*  X=currentposition[0], Y=currentposition[1], Z=currentposition[2],         */
/*                                                                            */
/******************************************************************************/
void      ftotGrad(double *fgradp, double *currentxyzp)
{ int       i;
  double    gradBp[3] = {0.0, 0.0, 0.0};

  fgradp[0] = fgradp[1] = fgradp[2] = 0.0;
  for (i=0;  i < NPTS ;  i++)
  { 
    fbasisGrad(gradBp, currentxyzp, ptsCentresVec.p + 3*i);
    fgradp[0] += weightVec.p[i] * gradBp[0];
    fgradp[1] += weightVec.p[i] * gradBp[1];
    fgradp[2] += weightVec.p[i] * gradBp[2];
/*
fprintf(stderr, "ftotGrad, iB=%2d, w=%f, gradB=(%f,%f,%f), partiel Tot=(%f,%f,%f)\n",
                i, weightVec.p[i], gradBp[0], gradBp[1], gradBp[2],
                fgradp[0], fgradp[1], fgradp[2] );
*/
  }

  return;
}
/******************************************************************************/
/*                                                                            */
/*  currentxyzp ~                                                             */
/*  X=currentposition[0], Y=currentposition[1], Z=currentposition[2],         */
/*                                                                            */
/******************************************************************************/
double    ftotValGrad(double *fgradp, double *currentxyzp)
{ int       i, id;
  double    f = 0.0;
  double   *gradtotp;
  double    gradp[3] = {0.0, 0.0, 0.0};

  for (i=0;  i < NPTS ;  i++)
  { f  += weightVec.p[i] * fbasisValGrad(gradp, currentxyzp, ptsCentresVec.p + 3*i);
    gradtotp = fgradp;
    for(id = 0;  id < 2;  id++, gradtotp++)
    { *gradtotp  += weightVec.p[i] * gradp[id];
    }
  }

  return f;
}
/******************************************************************************/
/******************************************************************************/
