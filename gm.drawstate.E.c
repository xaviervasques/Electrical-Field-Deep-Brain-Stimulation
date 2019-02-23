/*  URMAE/orientHaut/linear4.GL.V4/gm.drawstate.E.c                           */
/*  Mennessier Gerard                 20010613                                */
/*  Last Revised : G.M.               20040511                                */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utiMath.constant.def.h"

#include  "utistdIO.h"
#include  "utistdErr.h"

#include  "utiVecChr.h"
#include  "utiCurve.def.h"
#include  "utiCurve.set.h"
#include  "utiCurve.string.h"
#include  "utiCurve.2PS.h"

#include  "gm.drawstate.E.h"
#include  "gm.drawstate.E.statalloc.h"
#include  "gm.drawstate.glob.h"
#include  "gm.linear4.initsolve.glob.h"

#include  "pallidus.draw.h"
#include  "pallidus.geom.h"
#include  "gm.draw.E.h"
#include  "gm.drawstate.h"

/*
static    double    defaultELimRho[2] = {-4.65, +4.65};
static    double    defaultELimZ[2]   = {-6.10, +6.10};
static    double    defaultELimX[2] = {-4.05, +4.05};
static    double    defaultELimY[2] = {-4.05, +4.05};

static    double    defaultELimRho[2] = {-6.1, +6.1};
static    double    defaultELimZ[2]   = {-8.6, +8.6};
static    double    defaultELimX[2] = {-6.1, +6.1};
static    double    defaultELimY[2] = {-6.1, +6.1};

static    double    defaultELimRho[2] = { -8.1,  +8.1};
static    double    defaultELimZ[2]   = {-11.5, +11.5};
static    double    defaultELimX[2] = {-8.1, +8.1};
static    double    defaultELimY[2] = {-8.1, +8.1};

static    double    defaultELimRho[2] = {-12. , +12.};
static    double    defaultELimZ[2]   = {-17. , +12.};
static    double    defaultELimX[2] = {-12. , +12.};
static    double    defaultELimY[2] = {-12. , +12.};

static    double    defaultELimRho[2] = {-42., +42.};
static    double    defaultELimZ[2]   = {-59., +59.};
static    double    defaultELimX[2] = {-42., +42.};
static    double    defaultELimY[2] = {-42., +42.};
*/

static    double    defaultELimRho[2] = {-12. , +12.};
static    double    defaultELimZ[2]   = {-17. , +12.};
static    double    defaultELimX[2] = {-12. , +12.};
static    double    defaultELimY[2] = {-12. , +12.};
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateEInit()
{ static    char      prognamp[] = "gm.drawstate.E::gmDrawStateEInit";

  fPrintF(stderr,"%s  BEGIN\n", prognamp);

  chrPVecAlloc(&stateEpsiRotV, 12);
  gmDrawStateESetpsiRot(0.0);
  
  chrPVecAlloc(&stateEzhV, 12);
  gmDrawStateESetzh(-1.0);

  stateETitlesp[4] = stateEpsiRotV.p;
  stateETitlesp[5] = stateEzhV.p;

  gmDrawStateESetdefLimZ();
  gmDrawStateESetLimfromRhoZScaleRatio(1.0);

  fPrintF(stderr,"%s  END\n", prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateEGetpsiRot()
{ return   &stateEpsiRotV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateESetpsiRot(double psiDegree)
{ static    char    psistrp[8];
  double    statePsi = 0.0;              /** always 0. since one rotates pallidus **/
  lcSegVec *pallidusLcSgVp;
  static    char      prognamp[] = "gm.drawstate.E::gmDrawStateESetpsiRot";

  stateEpsiRotDeg = psiDegree;
  stateEpsiRotRad = stateEpsiRotDeg *(myPI/180.);
  psistrp[3] = psistrp[4] = psistrp[5] = psistrp[6] = psistrp[7] = 0;
  sprintf(psistrp, "%+.2f", stateEpsiRotDeg);

  stateEpsiRotV.x = 0;
  chrVecIncStr(&stateEpsiRotV, psistrp);
fPrintF(stderr, "     %s str psi=%s; num psi=%f degree\n",
                                        prognamp, stateEpsiRotV.p, stateEpsiRotDeg);
  stateETitlesp[4] = stateEpsiRotV.p;                    /** may have been change **/

  pallidusEBoundPsiRotate(stateEpsiRotRad);

  pallidusLcSgVp = pallidusERZGetlcSeg();
  pallidusEComputeRZSection(pallidusLcSgVp, statePsi, 0);

  pallidusLcSgVp = pallidusEXYGetlcSeg();
  pallidusEComputeXYSectionRotated(pallidusLcSgVp, stateEzh, 0);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateEGetzh()
{ return   &stateEzhV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateESetzh(double zh)
{ lcSegVec *pallidusLcSgVp;
  static    char    zhstrp[8];
  static    char    prognamp[] = "gm.drawstate.E::gmDrawStateESetzh";

  stateEzh = zh;
  zhstrp[3] = zhstrp[4] = zhstrp[5] = zhstrp[6] = zhstrp[7] = 0;
  sprintf(zhstrp, "%+.2f", stateEzh);

  stateEzhV.x = 0;
  chrVecIncStr(&stateEzhV, zhstrp);
fPrintF(stderr, "     %s str zh=%s; num zh=%f\n", prognamp, stateEzhV.p, stateEzh);
  stateETitlesp[5] = stateEzhV.p;                        /** may have been change **/

  pallidusLcSgVp = pallidusEXYGetlcSeg();
  pallidusEComputeXYSectionRotated(pallidusLcSgVp, stateEzh, 0);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateESetdefLimRho()
{ stateELimRho[0] = defaultELimRho[0];
  stateELimRho[1] = defaultELimRho[1];
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateESetdefLimZ()
{ stateELimZ[0] = defaultELimZ[0];
  stateELimZ[1] = defaultELimZ[1];
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateESetdefLimX()
{ stateELimX[0] = defaultELimX[0];
  stateELimX[1] = defaultELimX[1];
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateESetdefLimY()
{ stateELimY[0] = defaultELimY[0];
  stateELimY[1] = defaultELimY[1];
  return;
}
/******************************************************************************/
/*                                                                            */
/* given the E lim Z, compute symmetric Elim Rho,                             */
/* so that one get the given scale ratio h/w                                  */
/******************************************************************************/
void      gmDrawStateESetLimfromRhoZScaleRatio(double hwscaleratio)
{ double    deltah, deltaw;
  static    char    forma1p[] = "   %s\n"
         "   stateELimRho=%f %f, stateELimZ=%f %f,"
         "   stateELimX=%f %f, stateELimY=%f %f\n";
  static    char prognamp[] = "gm.drawstate.E::gmDrawStateESetLimfromRhoZScaleRatio";

  deltah = stateELimZ[1] - stateELimZ[0];
  deltaw = gmDrawStateHWScaleRatioA4Height2Width(hwscaleratio, deltah);
  stateELimRho[0] = - (deltaw/2);  stateELimRho[1] = (deltaw/2);

  stateELimX[0] = - stateELimRho[1];  stateELimX[1] =  stateELimRho[1];
  deltah = gmDrawStateHWScaleRatioA4Width2Height(hwscaleratio, deltaw);
  stateELimY[0] = - (deltah/2);  stateELimY[1] = (deltah/2);

fPrintF(stderr, forma1p, prognamp, stateELimRho[0], stateELimRho[1],
                         stateELimZ[0], stateELimZ[1],
                         stateELimX[0], stateELimX[1], stateELimY[0], stateELimY[1]);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double   *gmDrawStateEGetLimRho()
{ return    stateELimRho;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double   *gmDrawStateEGetLimZ()
{ return    stateELimZ;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double   *gmDrawStateEGetLimX()
{ return    stateELimX;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double   *gmDrawStateEGetLimY()
{ return    stateELimY;
}
/******************************************************************************/


/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cSetVec  *gmDrawStateEGetRZlevelCSetp()
{ 
  if(ERZlevelCSetv.p==NULL)
  {
   gmDrawStateESetRZallCSet();
  }
  return   &ERZlevelCSetv;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cSetVec  *gmDrawStateEGetRZallCSetp()
{ 
  return   &ERZallCSetv;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cSetVec  *gmDrawStateEGetXYallCSetp()
{ 
  return   &EXYallCSetv;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateESetRZallCSet()
{ cStrVec  *cstrTitlevp, *cstrRZPinIndicesvp;
  cSegVec  *csegsectionvp;

  cstrTitlevp = gmDrawEGetRZCStrTitle(stateTitlesp, stateVEI, stateETitlesp);

  ERZeltrdAxisCSetv.x = 0;
  gmDrawEGetRZeltrdAxis(&ERZeltrdAxisCSetv, &gridZi_mm, &gridRi_mm);
  cstrRZPinIndicesvp = gmDrawEGetCStrPinIndices();
  cSetVecInc1ptype(&ERZeltrdAxisCSetv, cstrRZPinIndicesvp, MY_CSTRING);

  ERZtitleCSetv.x =0;
  cSetVecInc1ptype(&ERZtitleCSetv, cstrTitlevp, MY_CSTRING);

  ERZlevelCSetv.x = 0;
  gmDrawEComputeRZLevels(&ERZlevelCSetv, &gridZi_mm, &gridRi_mm, stateFpp,
                                                                      &stateLCrange);
  ERZallCSetv.x = 0;
  csegsectionvp = gmDrawEGetRZtraceXYsection(stateEzh);
  cSetVecInc1ptype(&ERZeltrdAxisCSetv, csegsectionvp, MY_CSEGV);
  cSetVecAddSetVec(&ERZallCSetv, &ERZeltrdAxisCSetv);
  cSetVecAddSetVec(&ERZallCSetv, &ERZtitleCSetv);
  cSetVecAddSetVec(&ERZallCSetv, &ERZlevelCSetv);

  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateESetXYallCSet()
{ cStrVec  *cstrTitlevp;

  cstrTitlevp = gmDrawEGetXYCStrTitle(stateTitlesp, stateVEI, stateETitlesp);

  EXYeltrdAxisCSetv.x = 0;
  gmDrawEGetXYeltrdAxis(&EXYeltrdAxisCSetv, &gridZi_mm, &gridRi_mm, stateEzh);

  EXYtitleCSetv.x =0;
  cSetVecInc1ptype(&EXYtitleCSetv, cstrTitlevp, MY_CSTRING);

  EXYlevelCSetv.x = 0;
  gmDrawEComputeXYLevels(&EXYlevelCSetv, &gridZi_mm, &gridRi_mm, stateFpp,
                                                             &stateLCrange, stateEzh);
  EXYallCSetv.x = 0;
  cSetVecAddSetVec(&EXYallCSetv, &EXYeltrdAxisCSetv);
  cSetVecAddSetVec(&EXYallCSetv, &EXYtitleCSetv);
  cSetVecAddSetVec(&EXYallCSetv, &EXYlevelCSetv);

  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateESetPallidusCSet()
{ lcSegVec *lcSegRZvp, *lcSegXYvp;
  double    y = 0.0;
  int       iy = 0,  iz = 0;

  lcSegRZvp = pallidusERZGetlcSeg();
  pallidusEComputeRZSection(lcSegRZvp, y, iy);

  lcSegXYvp = pallidusEXYGetlcSeg();
  pallidusEComputeXYSectionRotated(lcSegXYvp, stateEzh, iz);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateEPSprint(FILE *psFilep)
{ lcSegVec *lcSegRZvp, *lcSegXYvp;

  PSgraphInit(psFilep);
               /** RZ page **/
  PSgraphPageInit(psFilep);
  PSgraphSetInitCTM(psFilep, stateELimRho, stateELimZ);
  PSgraphSet(psFilep, &ERZallCSetv);
  lcSegRZvp = pallidusERZGetlcSeg();
/*
  PSgraphLcSeg(psFilep, lcSegRZvp);
*/
  PSgraphPageEnd(psFilep);

               /** XY page **/
  PSgraphPageInit(psFilep);
  PSgraphSetInitCTM(psFilep, stateELimX, stateELimY);
  PSgraphSet(psFilep, &EXYallCSetv);
  lcSegXYvp = pallidusEXYGetlcSeg();
/*
  PSgraphLcSeg(psFilep, lcSegXYvp);
*/
  PSgraphPageEnd(psFilep);

  PSgraphEnd(psFilep);
  return;
}
/******************************************************************************/
/******************************************************************************/
