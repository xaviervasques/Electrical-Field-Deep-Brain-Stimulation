#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include   "funcBasis.globalloc.h"
#include   "funcBasis.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      fbasisSetR(double rayon)
{ 
  R = rayon;
  R2 = R * R;
  return;
}
/******************************************************************************/
/*                                                                            */
/* centerxyzp ~ centerxyz[0],centerxyz[1],centerxyz[2]                        */
/*                                                                            */
/******************************************************************************/
double    fbasis(double *currentxyzp, double *centerxyzp)
{ double    dx, dy, dz;
  double    d2, d, dm1;
  double    f = 0.0;
/*
  double   *dp, *dcenterp;
  dp = currentxyzp;  dcenterp = centerxyzp;
*/

  dx = *currentxyzp++ - *centerxyzp++;
  dy = *currentxyzp++ - *centerxyzp++;
  dz = *currentxyzp   - *centerxyzp;

  d2 = dx*dx + dy*dy + dz*dz;
  if(d2 < R2)
  { 
/* */
    d = sqrt(d2);  dm1 = 1.0 - d/R;
/* */
/*
    dm1 = 1.0 - d2/R2;
*/
    f = dm1 * dm1;
  }
  else
  { f = 0.0;
  }
  return f;
}
/******************************************************************************/
/*                                                                            */
/*  given  fgradp[0-2]  which must be ALREADY allocated                       */
/*         current point *currentxyzp and reference center *centerxyzp        */
/*  compute gradient and store it in fgradp[0-2]                              */
/*                                                                            */
/******************************************************************************/
void      fbasisGrad(double *fgradp, double *currentxyzp, double *centerxyzp)
{ double    dx, dy, dz;
  double    d2, dm1;
  double    d;
  double    dfsd2;
  double   *gradp;

  dx = *currentxyzp++ - *centerxyzp++;
  dy = *currentxyzp++ - *centerxyzp++;
  dz = *currentxyzp   - *centerxyzp;
  gradp = fgradp;

  d2 = dx*dx + dy*dy + dz*dz;
  if(d2 < R2)
  { 
/* */
    d = sqrt(d2);  dm1 = 1.0 - d/R;
    dfsd2 = -2.0 * dm1/R;
    if(d == 0.0)
    { *gradp++ = 0.0;  *gradp++ = 0.0;  *gradp = 0.0;}
    else        
    { dfsd2 = -2.0 * dm1/(R * d);
      *gradp++ = dfsd2 * dx;  *gradp++ = dfsd2 * dy;  *gradp = dfsd2 * dz;
    }
/* */

/*
    dm1 = 1.0 - d2/R2;
    dfsd2 = - 2.0 * dm1/R2;
    *gradp++ = dfsd2 * dx;  *gradp++ = dfsd2 * dy;  *gradp = dfsd2 * dz;
 */

  }
  else
  { *gradp++ = 0.0;  *gradp++ = 0.0;  *gradp = 0.0;
  }
  return;
}
/******************************************************************************/
/*                                                                            */
/*  given  fgradp[0-2]  which must be ALREADY allocated                       */
/*         current point *currentxyzp and reference center *centerxyzp        */
/*  compute and return function value f                                       */
/*          gradient and store it in fgradp[0-2]                              */
/*                                                                            */
/******************************************************************************/
double    fbasisValGrad(double *fgradp, double *currentxyzp, double *centerxyzp)
{ double    dx, dy, dz;
  double    d2, d, dm1;
  double    dfsd2 = 0.0;
  double    f = 0.0;
  double   *gradp;

  dx = *currentxyzp++ - *centerxyzp++;
  dy = *currentxyzp++ - *centerxyzp++;
  dz = *currentxyzp   - *centerxyzp;
  gradp = fgradp;

  d2 = dx*dx + dy*dy + dz*dz;
  if(d2 < R2)
  { 
/* */
    d = sqrt(d2);  dm1 = 1.0 - d/R;
    f = dm1 * dm1;
    dfsd2 = -2.0 * dm1/R;
    if(d == 0.0)
    { *gradp++ = 0.0;  *gradp++ = 0.0;  *gradp = 0.0;}
    else        
    { dfsd2 = -2.0 * dm1/(R * d);
      *gradp++ = dfsd2 * dx;  *gradp++ = dfsd2 * dy;  *gradp = dfsd2 * dz;
    }
/* */

/*
    dm1 = 1.0 - d2/R2;
    f = dm1 * dm1;
    dfsd2 = - 2.0 * dm1/R2;
    *gradp++ = dfsd2 * dx;  *gradp++ = dfsd2 * dy;  *gradp = dfsd2 * dz;
*/
  }
  else
  { f = 0.0;
    *gradp++ = 0.0;  *gradp++ = 0.0;  *gradp = 0.0;
  }

  return f;
}
/******************************************************************************/
/******************************************************************************/
