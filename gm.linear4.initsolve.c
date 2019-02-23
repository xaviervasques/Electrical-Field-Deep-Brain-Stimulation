/*  URMAE/orientHaut/linear4.GL.V4/gm.linear4.initsolve.c                     */
/*  Mennessier Gerard                 20010613                                */
/*  Last Revised : G.M.               20040210                                */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiTnsr.h"

#include  "gm.linear4.initsolve.statalloc.h"
#include  "gm.linear4.initsolve.globalloc.h"
#include  "gm.linear4.initsolve.h"

#include  "constant.phys.URMAE.def.h"
#include  "utiMath.constant.def.h"

#include  "command.linear4.gl.h"
#include  "base.linear4.h"

#include  "gridZ.linear4.h"
#include  "gridR.linear4.h"
#include  "gridRZ.linear4.h"
#include  "varGrid.linear4.h"

#include  "cylPot.linear4.h"
#include  "eltrdPot.linear4.h"
#include  "pot.linear4.h"
#include  "varPot.linear4.h"

#include  "f2gradf.h"

static    char *filenamepp[10] = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
static    FILE *bufpp[10]      = {NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};

/******************************************************************************/
/*                                                                            */
/*  gmGetFilenamepp()                                                         */
/******************************************************************************/
char    **gmGetFilenamepp()
{ return    filenamepp;
}
/******************************************************************************/
/*                                                                            */
/*  gmGetStreambufpp()                                                        */
/******************************************************************************/
FILE    **gmGetStreambufpp()
{ return    bufpp;
}
/******************************************************************************/
/*                                                                            */
/*  gmGetCylPfp()                                                             */
/******************************************************************************/
gmCylPot      *gmGetCylPfp()
{ return    &cylPf;
}
/******************************************************************************/
/*                                                                            */
/*  gmGetEltrdPfp()                                                           */
/******************************************************************************/
gmEltrdPot    *gmGetEltrdPfp()
{ return    &eltrdPf;
}
/******************************************************************************/
/*                                                                            */
/*  gmInitBase()                                                              */
/******************************************************************************/
void      gmInitBase(int argx, char *argvpp[])
{ int       i, ok;
  FILE     *bufOutp;

                                                                 /** Great Matrix **/
#define     GAL4X_   10
#define     GAC4X_    5
  double  **gapp;
                                                                 /** Small Matrix **/
  double  **sApp, **sBpp, *sBp;

  int       bi;
  char      newpinmodsp[6];                /** such as mp0p           + null char **/
  char      newneckmodp[6];                /** MONO or BI             + null char **/
  static    char      prognamp[] = "gm.linear4.initsolve::gmInitBase";

  fPrintF(stderr,"%s  BEGIN, argx=%d, argvpp[0]=%s\n", prognamp, argx, argvpp[0]);
  fPrintF(stderr,"Going to command_line \n");
  commandLine_linear_gl(bufpp, filenamepp, newpinmodsp, newneckmodp, &bi,
                                                                       argx, argvpp);
  fPrintF(stderr,"\n");
  for(i = 0;  i < 10;  i++)
  { fPrintF(stderr, "i= %2d, filenamep = %p, bufp = %p ", i,filenamepp[i], (void*)bufpp[i]);
    if(filenamepp[i] != NULL) fPrintF(stderr, ", %s;",  filenamepp[i]);
    fPrintF(stderr, "\n");
  }
  fPrintF(stderr,"\n");
  fPrintF(stderr,"newpinmodes=%s, newneckmode=%s, bi=%d \n",
                                                       newpinmodsp, newneckmodp, bi);
  fPrintF(stderr,"\n");

  bufOutp = bufpp[4];

  ok = basisGet(varVbasp, Vbasp, eltrdPbasp, cylPbasp, &varGridi, &gridRi, &gridZi,
                                             &gridRi_mm, &gridZi_mm, bufOutp, bufpp);

                                        /** allocate and fill in the GREAT MATRIX **/
  tnsr2dblPAlloc(&gAt, GAL4X_, GAC4X_);
  gAt.ev1x = GAL4X_ ;  gAt.ev2x = GAC4X_ ;  gapp = gAt.pp;
  basisGetGA4(varVbasp, &gAt);
  tnsr2dblPrint(bufOutp, &gAt);
  tnsr2dblPPrint(bufOutp, &gAt);
                                                            /** final  allocation **/
  gmEltrdPotZero(&eltrdPf);
  gmCylPotZero(&cylPf);
  gmVarPotPAllocZero(&varVf, &varGridi, &eltrdPf, &cylPf);
  gmPotPAllocZero(&Vf, &gridRi, &gridZi);
  Vf.varPotp = &varVf;
  gmPotPAllocZero(&gradVf, &gridRi, &gridZi);             /** allocation for grad **/
                                                  /** allocation for SMALL MATRIX **/
  tnsr2dblPAlloc(&sAt, GAC4X_, GAC4X_);
  sAt.ev1x = GAC4X_ ;  sAt.ev2x = GAC4X_ ;  sApp = sAt.pp;
  tnsr2dblPAlloc(&sBt, 1, GAC4X_);
  sBt.ev1x = 1;  sBt.ev2x = GAC4X_ ;  sBpp = sBt.pp;  sBp = sBt.p;
  fflush(bufOutp);

  fPrintF(stderr,"%s  END\n", prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmSolveFromBasis()                                                        */
/******************************************************************************/
void      gmSolveFromBasis(char *voltstrp, char *pinmodstrp)
{ int       bi;
  int       pinmodsip[5];
  double    veff, veffpin;
  char     *cp = NULL;
  double  **gapp, **sApp, **sBpp, *sBp;
  FILE     *bufOutp;
  static    char      prognamp[] = "gm.linear4.initsolve::gmSolveFromBasis";

fPrintF(stderr, "%s  pinmod=%s, voltage=%s \n", prognamp, pinmodstrp, voltstrp);
  bi = gmBiFromStrMod(pinmodsip, pinmodstrp);
  cylPf.neckMod = bi;

fPrintF(stderr, "  bi=%d, pinmodi=%d %d %d %d \n",
                            bi, pinmodsip[0],pinmodsip[1],pinmodsip[2],pinmodsip[3]);
fPrintF(stderr, "%s   voltstrp=%p, cpp=%p \n", prognamp, voltstrp, (void*)&cp);

  veff = strtod(voltstrp, &cp);
fPrintF(stderr, "%s  veff=%f, voltstrp=%p, cpp=%p, *cpp=%p \n",
                                          prognamp, veff,  voltstrp, (void*)&cp, cp);
  bufOutp = stderr;

  veffpin = veff;
  if(bi){ veffpin = veff * .5;}
  gapp = gAt.pp;
  sApp = sAt.pp;
  sBpp = sBt.pp;  sBp = sBt.p;
  basisGetSA4fromPINS(bufOutp, gapp, sApp, sBp, bi, pinmodstrp, veffpin);
  basisSolve4(bufOutp, &sAt, &sBt, pinmodsip, &cylPf, &eltrdPf, &Vf,
                                                        cylPbasp, eltrdPbasp, Vbasp);
  getGeo2DGradNorm(gridRi.rhondp, gridRi.rhotndx, gridZi.zndp, gridZi.zndx, 
                    gridZi.irGeond_firstp, gridZi.irGeond_lastp, Vf.vpp, gradVf.vpp);

  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmGetImpedance()                                                          */
/******************************************************************************/
double    gmGetImpedance()
{ double    kappa;
  double    elcur, predictedZ = 10000.;
  double   *pinVp, neckV, voltage;
  gmCylPot      *cylPfp;
  gmEltrdPot    *eltrdPfp;
  static    char      prognamp[] = "gm.linear4.initsolve::gmGetImpedance";

  kappa = KAPPA_ ;
  cylPfp = gmGetCylPfp();  eltrdPfp = gmGetEltrdPfp();
  elcur = fabs(cylPfp->neckI) + fabs(eltrdPfp->effI);  elcur *= 0.5;
  elcur *= (2. * myPI * ELTRD_R_m * kappa);
  neckV = cylPfp->infiniV;  pinVp = eltrdPfp->pinVp;
  voltage = getDiffPot(pinVp, neckV);
  if(elcur) predictedZ = voltage/elcur;
fPrintF(stderr, "%s cur=%f, v=%f, neckV=%f, Z=%f\n",
                                        prognamp, elcur, voltage, neckV, predictedZ);
  return  fabs(predictedZ);
}
/******************************************************************************/
/******************************************************************************/
