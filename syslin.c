/*  URMAE/reconstruction.GM/V1/syslin.c                                       */
/*  Mennessier Gerard                 20060321                                */
/*  Last Revised : G.M.               20060325                                */
/*                                                                            */

#include  <stddef.h>
#include  <stdio.h>
#include  "utiVecDbl.h"
#include  "utiTnsr.h"

#include  "syslin.h"
#include  "syslin.GaussSeidel.h"
#include  "syslin.LU.h"
#include  "funcTot.h"

/******************************************************************************/
/*                                                                            */
/* syslin                                                                     */
/*                                                                            */
/******************************************************************************/
void      syslin(size_t npts, double *centerxyzp, double *weightp)
{ int       i;
  int       itermax = 50;
  tnsr2dbl  xt = {0,0,0,0,0, NULL,NULL};
  dblVec    bVec = {0,0, NULL}, aVec = {0,0, NULL}, xVec = {0,0, NULL};
  double   *dp;
  double    eps = 1.0e-7;
  double   *lpcentrep;
  double    f, gradp[3];

  tnsr2dblPRealloc(&xt, itermax, npts);

  dblPVecRealloc(&bVec, npts, 10);
  dp = bVec.p;
  for(i = 0;  i < npts;  i++) *dp++ = 1.0;

  dblPVecRealloc(&aVec, npts, 10);
  dp = aVec.p;
  for(i = 0;  i < npts;  i++) *dp++ = 0.0;

  dblPVecRealloc(&xVec, npts, 10);
  dp = xVec.p;
  for(i = 0;  i < npts;  i++) *dp++ = 0.0;


                                                  /** Resolution systeme lineaire **/
/*
fprintf(stderr, "syslin::syslin Calling sl_gauss_seidel \n");
  sl_gauss_seidel(npts, aVec.p, bVec.p, xVec.p, itermax, xt.pp, eps);
 */

fprintf(stderr, "syslin::syslin Calling sl_LU \n");
  sl_LU(npts, bVec.p, xVec.p);


fprintf(stderr, "\n WEIGHT \n");
  for(i = 0;  i < npts;  i++)
  { weightp[i] = xVec.p[i];
fprintf(stderr, " i=%4d, w= %f \n", i, weightp[i]);
  }
fprintf(stderr, "\n");

                                                               /** CHECK SOLUTION **/
fprintf(stderr, "  CHECK RESOLUTION SYSTEME \n");
  lpcentrep = centerxyzp;
  for(i = 0;  i < npts;  i++, lpcentrep +=3)
  { f = ftot(lpcentrep);
fprintf(stderr, "syslin  i=%4d, point=(%f, %f, %f) fTot=%12.4e \n",
               i, *lpcentrep, *(lpcentrep+1), *(lpcentrep+2), f );
    ftotGrad(gradp, lpcentrep);
fprintf(stderr, "syslin  gradient=(%f, %f, %f) \n",  gradp[0], gradp[1], gradp[2]);
  }
fprintf(stderr, "  CHECK RESOLUTION  FIN \n\n");

  return;
}
/******************************************************************************/
/******************************************************************************/
