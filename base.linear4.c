/*  URMAE/orientHaut/linear4.Common.V1/base.linear4.c                         */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20030515                                */

#include  <stddef.h>
#include  <math.h>
#include  <time.h>                                         /** for time and ctime **/

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiTnsr.h"

#include  "eltrd.def.h"
#include  "gridZ.linear4.delta.def.h"
#include  "gridR.linear4.delta.def.h"



#include  "base.linear4.h"

#include  "gridZ.linear4.h"
#include  "gridR.linear4.h"
#include  "gridRZ.linear4.h"
#include  "varGrid.linear4.h"

#include  "cylPot.linear4.h"
#include  "eltrdPot.linear4.h"
#include  "pot.linear4.h"
#include  "varPot.linear4.h"
#include  "copy.linear4.h"

#include  "recipe.lubksb.h"
#include  "recipe.ludcmp.h"
/******************************************************************************/
/*                                                                            */
/**  bufpp[10] =                                                             **/
/**  Basisp[0], ...,Basisp[3], any, any, GridZ, GridR,                       **/
/******************************************************************************/
int       basisGet(gmVarPot *varVbsp, gmPot *Vbsp, 
                   gmEltrdPot *eltrdPbsp, gmCylPot *cylPbsp,
                   gmVarGrid *varGridip, gmGridR *gridRip, gmGridZ *gridZip,
                   gmGridR *gridRi_mmp, gmGridZ *gridZi_mmp,
                                                        FILE *bufOutp, FILE *bufpp[])
{ int       ok = 1;
  time_t    timesec;

  double    infozpap    [MYDELTAZMN_ +1];
  int       infopinzpaxp[MYDELTAZMN_ +1];
  double    inforpap[MYDELTARN_ +1];
  FILE     *bufgriZIp, *bufgriRIp, *bufReadBasep;
  int       imb;

  int       nVarNodes;
  int       boolPrintIni;

  int       cylcond;                    /**  -1 for CYLI, 0 for NECK, +1 for CYLC **/
  int       ordersp[NPEX_] = {0,0,0,0,0,0,0,0};
  static    char    timeformp[] = "     %s \n";
  static    char    prognamp[] = "linear4/base.linear4::basisGet";


fPrintF(stderr," %s  BEGIN\n", prognamp);
                                                              /** Global Geometry **/
                             /** Z grid initialisation  imposed UNIT = millimeter **/
  bufgriZIp = bufpp[6];
fPrintF(stderr," %s  reading Z info\n", prognamp);
  gmGridZreadINI(bufgriZIp, infozpap, infopinzpaxp);
  fclose(bufgriZIp);
fPrintF(stderr," %s  reading Z info DONE\n", prognamp);   fflush(stderr);
/* fPrintF(stderr," gridZi_mm init  BEGIN\n");  fflush(stderr);  */
  gmGridZinit(bufOutp, gridZi_mmp, infozpap, infopinzpaxp);
/* */
  gmGridZprint(bufOutp, gridZi_mmp);     fPrintF(bufOutp,"\n");    fflush(bufOutp);
/* */
fPrintF(stderr," gridZi_mm init  DONE\n");  fflush(stderr);
fPrintF(stderr,"\n");

                             /** R grid initialisation  imposed UNIT = millimeter **/
  bufgriRIp = bufpp[7];
fPrintF(stderr," %s  reading R info\n", prognamp);
  gmGridRreadINI(bufgriRIp, inforpap, &cylcond);
  fclose(bufgriRIp);
fPrintF(stderr," %s  reading R info DONE\n", prognamp);   fflush(stderr);
fPrintF(stderr," cylcond = %d \n", cylcond);  
fPrintF(stderr," neckRadius=%f\n", *inforpap);
  gmGridRinit(bufOutp, gridRi_mmp, inforpap, cylcond);
  gmGridRprint(bufOutp, gridRi_mmp);     fPrintF(bufOutp,"\n");    fflush(bufOutp);
fPrintF(stderr," gridRi_mm init  DONE\n");
fPrintF(stderr,"\n");



                                                  /** R Z geometry initialisation **/
  gmGridRZGeoinit(bufOutp, gridRi_mmp, gridZi_mmp);
  gmGridZprint(bufOutp, gridZi_mmp);     fPrintF(bufOutp,"\n");    fflush(bufOutp);


                                /** Z grid initialisation  UNIT = electrod radius **/
  gmGridZPAlloc(gridZi_mmp->zllx, gridZip);
  gmGridZCpy("eltrdRadius", gridZi_mmp->eltrdR, gridZip, gridZi_mmp);
  gmGridZprint(bufOutp, gridZip);        fPrintF(bufOutp,"\n");    fflush(bufOutp);
fPrintF(stderr," gridZi  DONE\n");
fPrintF(stderr,"\n");
                                /** R grid initialisation  UNIT = electrod radius **/
  gmGridRPAlloc(gridRi_mmp->rhotllx, gridRip);
  gmGridRCpy("eltrdRadius", gridRi_mmp->eltrdR, gridRip, gridRi_mmp);
  gmGridRprint(bufOutp, gridRip);        fPrintF(bufOutp,"\n");    fflush(bufOutp);
fPrintF(stderr," gridRi  DONE\n");
fPrintF(stderr,"\n");
                                          /** Geometry for minimization variables **/
fPrintF(stderr," varGrid init CALLING\n");
  nVarNodes = gmVarGridInit(varGridip, gridRip, gridZip, ordersp);
  gmVarGridPrint(bufOutp, varGridip);    fPrintF(bufOutp,"\n");    fflush(bufOutp);
fPrintF(stderr," varGrid print DONE\n");


                                                  /** Potential INITIALISATION **/
fPrintF(stderr,"\n READING VARPOT BASIS BEGIN\n");
  boolPrintIni = 1;
  boolPrintIni = 0;
  /*int i;
  scanf("%i\n",i);*/

  for(imb = 0;  imb < 4;  imb++)
  { 
    gmEltrdPotZero(eltrdPbsp + imb);
    gmCylPotZero(cylPbsp + imb);

    gmVarPotPAllocZero(varVbsp + imb, varGridip, eltrdPbsp + imb, cylPbsp + imb);
    bufReadBasep = bufpp[imb];

    gmVarPotRead(bufReadBasep, varVbsp + imb);  fclose(bufReadBasep);
    if(boolPrintIni)
    { fPrintF(bufOutp,"MAIN: Printing  Initial varPot %d\n", imb);
      gmVarPotPrint(bufOutp, varVbsp + imb);
      fPrintF(bufOutp,"\n");  fflush(bufOutp);
    }
    gmPotPAllocZero(Vbsp + imb, gridRip, gridZip);
    gmVarP2P(Vbsp + imb, varVbsp + imb);
    fPrintF(stderr,"\n");
  }
fPrintF(stderr," READING VARPOT BASIS DONE\n");
  /*timesec = time(&timesec);
  fPrintF(bufOutp, timeformp, ctime(&timesec));  fflush(bufOutp);
fPrintF(stderr," %s  END\n", prognamp);
  */
  /*int i;
  scanf("%i\n",i);*/
  return ok;
}

/******************************************************************************/
/*  double  **basisGetGA4(gmVarPot *varVbp)                                   */
/*                                                                            */
/*                  p000  0p00  00p0  000p    V
          pin0V
          pin1V  
          pin2V
          pin3V
          neckV
          pin0I
          pin1I
          pin2I
          pin3I
          neckI                                                               */
/*                                                                            */
/******************************************************************************/
void      basisGetGA4(gmVarPot *varVbp, tnsr2dbl *gAtp)
{ size_t    ev1x, ev2x;
  int       ic, is;
  double   *gap, **gApp;
  gmVarPot *varVp;
  gmEltrdPot    *eltrdp;
  gmCylPot      *cylp;

  ev1x = gAtp->ev1x;  ev2x = gAtp->ev2x;  gApp = gAtp->pp;
  for(ic = 0, varVp = varVbp;  ic < 4; ic++, varVp++)
  { eltrdp = varVp->eltrdPotp;
    cylp = varVp->cylPotp;
    for(is = 0;  is < 4; is++){ gap = gApp[is];      gap[ic] = eltrdp->pinVp[is];}
                                gap = gApp[4];       gap[ic] = cylp->infiniV;
    for(is = 0;  is < 4; is++){ gap = gApp[5 + is];  gap[ic] = eltrdp->pinIp[is];}
                                gap = gApp[9];       gap[ic] = cylp->neckI;
  }
                                                                 /** set V column **/
  for(is = 0;  is <= 4;  is++){ gap = gApp[is];  gap[4] = 1.0;}
  for(is = 5;  is <= 9;  is++){ gap = gApp[is];  gap[4] = 0.0;}

  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      basisGetSA4fromPINS(FILE *bufp, double **gapp, double **sapp, double *sbp,
                                               int bi, char *pinmodstrp, double veff)
{ int       il = 0, ic, ip;
  double   *sap, *gap;
  char      c;

  for(ip = 0;  ip < 4; ip++, pinmodstrp++)
  { c = *pinmodstrp;
    if(c == '0'){ il = 5+ip;  sbp[ip] = 0.0;}
    else if(c == 'p')            { il = ip;  sbp[ip] =  veff;}
    else if(c == 'm' || c == 'n'){ il = ip;  sbp[ip] = -veff;}
    gap = gapp[il];  sap = sapp[ip];
    for(ic = 0;  ic < 5; ic++){ *sap++ = *gap++;}
  }
  sap = sapp[4];
  if(bi){ il = 9;}
  else  { il = 4;}
  gap = gapp[il];
  for(ic = 0;  ic < 5; ic++){ *sap++ = *gap++;}
  sbp[4] = 0.0;

  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      basisSolve4(FILE *bufp, tnsr2dbl *satp, tnsr2dbl *sbtp, int *newpinmodsip, 
                      gmCylPot *cylPfp, gmEltrdPot *eltrdPfp, gmPot *Vfp,
                      gmCylPot *cylPbsp, gmEltrdPot *eltrdPbsp, gmPot *Vbsp)
{ size_t    evx;
  int       indexsavp[5];
  double  **sApp, **sBpp, *sBp;
  double    det;

  evx = satp->ev1x;  sApp = satp->pp;
  sBpp = sbtp->pp;  sBp = sbtp->p;
                                                                        /** solve **/
  ludcmp(sApp, evx, indexsavp, &det);
  fPrintF(bufp,"\n");  fPrintF(bufp," sA matrix after triangle\n");
  tnsr2dblPPrint(bufp, satp);

  lubksb(sApp, evx, indexsavp, sBp);
  fPrintF(bufp,"\n");  fPrintF(bufp," sB after solve \n");
  tnsr2dblPPrint(bufp, sbtp);
                                          /** get Potential by linear combination **/
  copyLinearCyl(sBp, 4, sBp[4], cylPfp, cylPbsp);
  copyLinearEltrd(sBp, 4, sBp[4], eltrdPfp, eltrdPbsp);
  gmEltrdPotSetModes(eltrdPfp, newpinmodsip);
  copyLinearPot(sBp, 4, sBp[4], Vfp, Vbsp);
/*  gmPotPrint(bufp, Vfp); */
  return;
}
/******************************************************************************/
/******************************************************************************/
