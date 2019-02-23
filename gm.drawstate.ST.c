/*  URMAE/orientHaut/linear4.GL.V4/gm.drawstate.ST.c                          */
/*  Mennessier Gerard                 20010613                                */
/*  Last Revised : G.M.               20040511                                */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utistdIO.h"
#include  "utistdErr.h"

#include  "utiMath.constant.def.h"

#include  "utiVecChr.h"
#include  "utiCurve.def.h"
#include  "utiCurve.set.h"
#include  "utiCurve.string.h"
#include  "utiCurve.2PS.h"

#include  "gm.drawstate.ST.h"
#include  "gm.drawstate.ST.statalloc.h"
#include  "gm.drawstate.glob.h"
#include  "gm.linear4.initsolve.glob.h"


#include  "eltrdPot.linear4.h"           /** for PinMod string to int and Neckmod **/
#include  "grid_pot.ST.linear4.h"
#include  "gm.linear4.initsolve.h"
#include  "pallidus.draw.h"
#include  "pallidus.geom.h"
#include  "gm.draw.ST.h"
#include  "gm.drawstate.h"

/*
static    double    defaultSTLimRho[2] = { -8.1,  +8.1};
static    double    defaultSTLimZ[2]   = {-11.5, +11.5};
static    double    defaultSTLimX[2] = {-8.1, +8.1};
static    double    defaultSTLimY[2] = {-8.1, +8.1};

static    double    defaultSTLimRho[2] = {-12. , +12.};
static    double    defaultSTLimZ[2]   = {-17. , +12.};
static    double    defaultSTLimX[2] = {-12. , +12.};
static    double    defaultSTLimY[2] = {-12. , +12.};

static    double    defaultSTLimRho[2] = {-42., +42.};
static    double    defaultSTLimZ[2]   = {-59., +59.};
static    double    defaultSTLimX[2] = {-42., +42.};
static    double    defaultSTLimY[2] = {-42., +42.};
*/

static    double    defaultSTLimRho[2] = {-12. , +12.};
static    double    defaultSTLimZ[2]   = {-17. , +12.};
static    double    defaultSTLimX[2] = {-12. , +12.};
static    double    defaultSTLimY[2] = {-12. , +12.};
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTInit()
{ static    char      prognamp[] = "gm.drawstate.ST::gmDrawStateSTInit";
  double    zst;

  fPrintF(stderr,"%s  BEGIN\n", prognamp);

  chrPVecAlloc(&stateSTpsiRotV, 12);
  gmDrawStateSTSetpsiRot(0.0);
  
  gmGridPotSTXYPAlloc(101, gridRi_mm.rhotndx, &stateSTGridPotalphsax);
  gmGridPotSTXYsetAlpha(&stateSTGridPotalphsax);
  gmGridPotSTXYsetSaxis(gridRi_mm.rhotndx, gridRi_mm.rhondp, &stateSTGridPotalphsax);

  chrPVecAlloc(&stateSTzhV, 12);

  zst = -1.0;
  gmDrawStateSTSetzh(zst);

  stateSTTitlesp[4] = stateSTpsiRotV.p;
  stateSTTitlesp[5] = stateSTzhV.p;

  gmDrawStateSTSetdefLimZ();
  gmDrawStateSTSetdefLimRho();
  gmDrawStateSTSetdefLimX();
  gmDrawStateSTSetdefLimY();
  gmDrawStateSTSetLimfromRhoZScaleRatio(1.0);

  fPrintF(stderr,"%s  END\n", prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateSTGetpsiRotV()
{ return   &stateSTpsiRotV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double    gmDrawStateSTGetpsiRot()
{ return    stateSTpsiRotDeg;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTSetpsiRot(double psiDegree)
{ static    char    psistrp[8];
  lcSegVec *pallidusLcSgVp;
  static    char      prognamp[] = "gm.drawstate.ST::gmDrawStateSTSetpsiRot";

  stateSTpsiRotDeg = psiDegree;
  stateSTpsiRotRad = stateSTpsiRotDeg *(myPI/180.);
  psistrp[3] = psistrp[4] = psistrp[5] = psistrp[6] = psistrp[7] = 0;
  sprintf(psistrp, "%+.2f", stateSTpsiRotDeg);

  stateSTpsiRotV.x = 0;
  chrVecIncStr(&stateSTpsiRotV, psistrp);
fPrintF(stderr, "     %s str psi=%s; num psi=%f degree\n",
                                      prognamp, stateSTpsiRotV.p, stateSTpsiRotDeg);
  stateSTTitlesp[4] = stateSTpsiRotV.p;                  /** may have been change **/

  pallidusSTBoundPsiRotate( - stateSTpsiRotRad);

  pallidusLcSgVp = pallidusSTRZGetlcSeg();
  pallidusSTComputeRZSection(pallidusLcSgVp, 0.0, 0);          /** always y = 0.0 **/

  pallidusLcSgVp = pallidusSTXYGetlcSeg();
  pallidusSTComputeXYSection(pallidusLcSgVp, stateSTzh, 0);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateSTGetzhV()
{ return   &stateSTzhV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double    gmDrawStateSTGetzh()
{ return    stateSTzh;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTSetzh(double zh)
{ lcSegVec *pallidusLcSgVp;
  static    char    zhstrp[8];
  static    char    prognamp[] = "gm.drawstate.ST::gmDrawStateSTSetzh";

  stateSTzh = zh;
  zhstrp[3] = zhstrp[4] = zhstrp[5] = zhstrp[6] = zhstrp[7] = 0;
  sprintf(zhstrp, "%+.2f", stateSTzh);

  stateSTzhV.x = 0;
  chrVecIncStr(&stateSTzhV, zhstrp);
fPrintF(stderr, "     %s str zh=%s; num zh=%f\n", prognamp, stateSTzhV.p, stateSTzh);
  stateSTTitlesp[5] = stateSTzhV.p;                      /** may have been change **/

  gmGridPotSTXYsetPot(stateFpp, stateSTzh, &stateSTGridPotalphsax);
/*
fPrintF(stderr,  "     %s stateSTGridPotST alpha saxis\n", prognamp);
gmGridPotSTXYPrint(stderr, &stateSTGridPotalphsax); 
*/  
  pallidusLcSgVp = pallidusSTXYGetlcSeg();
  pallidusSTComputeXYSection(pallidusLcSgVp, stateSTzh, 0);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTSetdefLimRho()
{ stateSTLimRho[0] = defaultSTLimRho[0];
  stateSTLimRho[1] = defaultSTLimRho[1];
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTSetdefLimZ()
{ stateSTLimZ[0] = defaultSTLimZ[0];
  stateSTLimZ[1] = defaultSTLimZ[1];
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTSetdefLimX()
{ stateSTLimX[0] = defaultSTLimX[0];
  stateSTLimX[1] = defaultSTLimX[1];
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTSetdefLimY()
{ stateSTLimY[0] = defaultSTLimY[0];
  stateSTLimY[1] = defaultSTLimY[1];
  return;
}
/******************************************************************************/
/*                                                                            */
/* given the ST lim Z, compute symmetric STlim Rho,                           */
/* so that one get the given scale ratio h/w                                  */
/******************************************************************************/
void      gmDrawStateSTSetLimfromRhoZScaleRatio(double hwscaleratio)
{ double    deltah, deltaw;
  static    char    forma1p[] = "   %s\n"
         "   stateSTLimRho=%f %f, stateSTLimZ=%f %f,"
         "   stateSTLimX=%f %f, stateSTLimY=%f %f\n";
  static    char prognamp[] = "gm.drawstate.ST::gmDrawStateSTSetLimfromRhoZScaleRatio";

  deltah = stateSTLimZ[1] - stateSTLimZ[0];
  deltaw = gmDrawStateHWScaleRatioA4Height2Width(hwscaleratio, deltah);
  stateSTLimRho[0] = - (deltaw/2);  stateSTLimRho[1] = (deltaw/2);

  stateSTLimX[0] = - stateSTLimRho[1];  stateSTLimX[1] =  stateSTLimRho[1];
  deltah = gmDrawStateHWScaleRatioA4Width2Height(hwscaleratio, deltaw);
  stateSTLimY[0] = - (deltah/2);  stateSTLimY[1] = (deltah/2);

fPrintF(stderr, forma1p, prognamp, stateSTLimRho[0], stateSTLimRho[1],
                     stateSTLimZ[0], stateSTLimZ[1],
                     stateSTLimX[0], stateSTLimX[1], stateSTLimY[0], stateSTLimY[1]);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double   *gmDrawStateSTGetLimRho()
{ return    stateSTLimRho;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double   *gmDrawStateSTGetLimZ()
{ return    stateSTLimZ;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double   *gmDrawStateSTGetLimX()
{ return    stateSTLimX;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double   *gmDrawStateSTGetLimY()
{ return    stateSTLimY;
}
/******************************************************************************/


/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cSetVec  *gmDrawStateSTGetRZallCSetp()
{ 
  return   &STRZallCSetv;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cSetVec  *gmDrawStateSTGetXYallCSetp()
{ 
  return   &STXYallCSetv;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTSetRZallCSet()
{ cSegVec  *csegsectionvp;
  static    char      prognamp[] = "gm.drawstate.ST::gmDrawStateSTSetRZallCSet";

fPrintF(stderr, "%s begin\n", prognamp);

  STRZallCSetv.x = 0;
  
  STRZeltrdAxisCSetv.x = 0;
  gmDrawSTGetRZeltrdAxis(&STRZeltrdAxisCSetv, &gridZi_mm, &gridRi_mm);
  
  fPrintF(stderr, "%s gm 1****************************\n", prognamp);
  
  csegsectionvp = gmDrawSTGetRZtraceXYsection(stateSTzh);
  
    fPrintF(stderr, "%s gm 2****************************\n", prognamp);
  
  cSetVecInc1ptype(&STRZeltrdAxisCSetv, csegsectionvp, MY_CSEGV);
  
    fPrintF(stderr, "%s gm 3****************************\n", prognamp);
    
  cSetVecAddSetVec(&STRZallCSetv, &STRZeltrdAxisCSetv);

  STRZtitleCSetv.x =0;
  STRZlevelCSetv.x = 0;
  
  fPrintF(stderr, "%s gm 4****************************\n", prognamp);
  gmDrawSTComputeRZLevels(&STRZlevelCSetv);
  
    fPrintF(stderr, "%s gm 5****************************\n", prognamp);
  
  cSetVecAddSetVec(&STRZallCSetv, &STRZlevelCSetv);
  
  fPrintF(stderr, "%s end\n", prognamp);

  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTSetXYallCSet()
{ cSegVec  *csegsectionvp;
  static    char    prognamp[] = "gm.drawstate.ST::gmDrawStateSTSetXYallCSet";

  STXYallCSetv.x = 0;

  STXYeltrdAxisCSetv.x = 0;
  gmDrawSTGetXYeltrdAxis(&STXYeltrdAxisCSetv, &gridZi_mm, &gridRi_mm, stateSTzh);
fPrintF(stdout, "%s  GetXYeltrdAxis DONE\n", prognamp);

  csegsectionvp = gmDrawSTGetXYtraceRZsection(stateSTpsiRotRad);
  cSetVecInc1ptype(&STXYeltrdAxisCSetv, csegsectionvp, MY_CSEGV);
fPrintF(stdout, "%s  GetXYtraceRZsection DONE\n", prognamp);
  cSetVecAddSetVec(&STXYallCSetv, &STXYeltrdAxisCSetv);

  STXYtitleCSetv.x = 0;

  STXYlevelCSetv.x = 0;
  gmDrawSTComputeXYLevels(&STXYlevelCSetv, &stateLCrange, &stateSTGridPotalphsax);
fPrintF(stdout, "%s  XY Level Curves Computed\n", prognamp);
/*
cSetVecPrint(stdout, &STXYlevelCSetv);
lcSegVecPrint(stdout, (STXYlevelCSetv.p)->p );
*/
  cSetVecAddSetVec(&STXYallCSetv, &STXYlevelCSetv);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTSetPallidusCSet()
{ lcSegVec *lcSegRZvp, *lcSegXYvp;
  double    y = 0.0;
  int       iy = 0,  iz = 0;

  lcSegRZvp = pallidusSTRZGetlcSeg();
  pallidusSTComputeRZSection(lcSegRZvp, y, iy);

  lcSegXYvp = pallidusSTXYGetlcSeg();
  pallidusSTComputeXYSection(lcSegXYvp, stateSTzh, iz);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSTPSprint(FILE *psFilep)
{ lcSegVec *lcSegRZvp, *lcSegXYvp;

  PSgraphInit(psFilep);
               /** RZ page **/
  PSgraphPageInit(psFilep);
  PSgraphSetInitCTM(psFilep, stateSTLimRho, stateSTLimZ);
  PSgraphSet(psFilep, &STRZallCSetv);
  lcSegRZvp = pallidusSTRZGetlcSeg();
/*
  PSgraphLcSeg(psFilep, lcSegRZvp);
*/
  PSgraphPageEnd(psFilep);

               /** XY page **/
  PSgraphPageInit(psFilep);
  PSgraphSetInitCTM(psFilep, stateSTLimX, stateSTLimY);
  PSgraphSet(psFilep, &STXYallCSetv);
  lcSegXYvp = pallidusSTXYGetlcSeg();
/*
  PSgraphLcSeg(psFilep, lcSegXYvp);
*/
  PSgraphPageEnd(psFilep);

  PSgraphEnd(psFilep);
  return;
}
/******************************************************************************/
/******************************************************************************/



