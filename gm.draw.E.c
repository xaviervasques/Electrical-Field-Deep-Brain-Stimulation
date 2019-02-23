/*  URMAE/orientHaut/linear4.GL.V4/gm.draw.E.c                                */
/*  Mennessier Gerard                 20010613                                */
/*  Last Revised : G.M.               20040210                                */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utiMath.type.def.h"

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiTnsr.h"
#include  "utiVecChr.h"
#include  "utiVecDbl.h"
#include  "utiVecInt.h"
#include  "utiClassRangeDbl.h"

#include  "utiCurve.def.h"
#include  "utiCurve.GC.h"
#include  "utiCurve.circle.h"
#include  "utiCurve.ellipse.h"
#include  "utiCurve.poly.h"
#include  "utiCurve.seg.h"
#include  "utiCurve.set.h"
#include  "utiCurve.string.h"

#include  "gm.draw.E.h"

#include  "gm.linear4.initsolve.glob.h"

#include  "gridZ.linear4.h"
#include  "gridR.linear4.h"
#include  "gridRZ.linear4.h"
#include  "varGrid.linear4.h"

#include  "cylPot.linear4.h"
#include  "eltrdPot.linear4.h"
#include  "pot.linear4.h"
#include  "varPot.linear4.h"

#include  "f2gradf.h"
#include  "gm.drawstate.h"
#include  "gm.drawstate.E.h"

/******************************************************************************/
/*                                                                            */
/*  gmDrawEGetRZeltrdAxis()                                                   */
/*  Electrode frame: drawing of elelectrode, pins and axis                    */
/******************************************************************************/
void      gmDrawEGetRZeltrdAxis(cSetVec *csetvp,
                                              gmGridZ *gridZp_mm, gmGridR *gridRp_mm)
{ static    cPoly     cpoly;
  static    cSegVec   csegvec1, csegvec2;   /** csegvec1: pins;  csegvec2: z axis **/
  static    cSegVec   csegvec3;      /** csegvec3: fixed radius part of electrode **/
  static    cEllipse  cellips;                    /** rounded electrode extremity **/
  static    curveGC   gc1, gc2;
  static    intVec    irGeofV = {0,0,NULL};
  size_t    zndx;
  double   *zndp, z, zext0, zext1;
  size_t    rhondx;
  double   *rhondp, rho;                                       /** origin on axis **/
  double   *znxp, *rhonxp;                                 /**  z and rho Elimits **/
  int       iz, iz0,iz1;
  int       ir, ir2, irEl;
  int      *irGeoip, *irGeofp;
  double    rEl,  ticksize = 0.08;

  cSetPVecRealloc(csetvp, 5, 1);  csetvp->x = 0;

  cPolyPRealloc(&cpoly, 10, 1);
  cSegPVecRealloc(&csegvec1, 10, 1);
  cSegPVecRealloc(&csegvec2, 10, 1);

  gc1.lineWidth = 1;
  gc2.lineWidth = 3;

  cpoly.gcp = &gc2;
  csegvec3.gcp = &gc2;  cellips.gcp = &gc2;

  csegvec1.gcp = &gc1;
  csegvec2.gcp = &gc1;

  zndx = gridZp_mm->zndx;
  zndp = gridZp_mm->zndp;
  rhondx =gridRp_mm->rhotndx;
  rhondp =gridRp_mm->rhondp;

  intPVecRealloc(&irGeofV, zndx, 10);
  irGeofp = irGeofV.p;
  for(iz=0;  iz < zndx;  iz++){ irGeofp[iz] = rhondx-1;}
  irGeofV.x = zndx;
  irGeoip = gridZp_mm->irGeond_firstp;
  irEl   = gridRp_mm->eltrdRhondI;
  rEl    = gridRp_mm->eltrdR;

  iz0 = gridZp_mm->eltrdExtIp[0];  
  iz1 = gridZp_mm->eltrdExtIp[1];
  zext0 = gridZp_mm->eltrdExtZp[0];
  zext1 = gridZp_mm->eltrdExtZp[1];
                                                            /** electrode contour **/
                                                            /** fixed radius part **/
  csegvec3.x = 0;
  ir = irGeoip[iz1];
  z = zndp[0];                                         /** lower cylindrical part **/
  cSegVecInc1CoorData(&csegvec3, rEl, zext1, rEl, z);
  cSegVecInc1CoorData(&csegvec3,-rEl, zext1,-rEl, z);
  cSetVecInc1ptype(csetvp, &csegvec3, MY_CSEGV);
                                                        /** ellipsoidal extremity **/
  cellips.xp[0] = 0.0;  cellips.xp[1] = zext1;
  cellips.orient = +90.0;
  cellips.axisp[0] = fabs(zext1 - zext0);
  cellips.axisp[1] = rEl;
  cellips.anglDp[0] = -90.;  cellips.anglDp[1] = 90.0;
  cellips.closed = 0;
  cSetVecInc1ptype(csetvp, &cellips, MY_CELLIPSE);

                                                                         /** pins **/
  for(iz=0;  iz < 8;  iz++)
  { z = gridZp_mm->pinExtZp[iz];
    cSegVecInc1CoorData(&csegvec1, rhondp[0], z, rhondp[irEl], z);
  }
  cSetVecInc1ptype(csetvp, &csegvec1, MY_CSEGV);
/*
cSegVecPrint(stderr, &csegvec1);
*/
                                                              /** revolution axis **/
  znxp = gmDrawStateEGetLimZ();
  csegvec2.x = 0;
  rho = rhondp[0];
  cSegVecInc1CoorData(&csegvec2, rho, znxp[0], rho, znxp[1]);
                                                                      /** z ticks **/
  iz0 = (int)znxp[0];  iz1 = (int)znxp[1];
  for(iz = iz0;  iz <= iz1;  iz++)
  { z = iz * 1.0;  
    cSegVecInc1CoorData(&csegvec2, rho, z, rho - ticksize, z);
  }
                                                         /** x axis at z = 0. mm **/
  rhonxp = gmDrawStateEGetLimRho();
  z = 0. ;
  cSegVecInc1CoorData(&csegvec2, rhonxp[0], z, rhonxp[1], z);
  ir2 = (int)rhonxp[1];
  for(ir = 0;  ir <= ir2; ir++)
  { rho = ir * 1.0;
    cSegVecInc1CoorData(&csegvec2, rho, z, rho, z - ticksize);      /** x ticks **/
  }
  cSetVecInc1ptype(csetvp, &csegvec2, MY_CSEGV);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawEGetRZtraceXYsection()                                              */
/*  Electrode frame: drawing of trace in RZplane of the XYsection altitude zh */
/******************************************************************************/
cSegVec  *gmDrawEGetRZtraceXYsection(double zh)
{ static    cSeg      cseg;
  static    cSegVec   csegvec;
  static    curveGC   gc;
  double   *rhonxp;

  rhonxp = gmDrawStateEGetLimRho();
  gc.lineWidth = 2;
  csegvec.gcp = &gc;
  csegvec.p = &cseg;
  csegvec.z = 1;
  csegvec.x = 0;
  cSegVecInc1CoorData(&csegvec, 0.0, zh, rhonxp[1], zh);
  return  &csegvec;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawEGetXYeltrdAxis()                                                   */
/*  Electrode frame: drawing of electrode and axis                            */
/******************************************************************************/
void      gmDrawEGetXYeltrdAxis(cSetVec *csetvp,
                                   gmGridZ *gridZp_mm, gmGridR *gridRp_mm, double zh)
{ static    cCircle   celtrd;
  static    cSegVec   csegvec2;
  static    cCircleVec  cflevelv = {0,0,NULL,NULL};
  static    cCircle   cflevel;
  static    curveGC   gc1, gc2;
  size_t    zndx;
  double   *zndp;
  double    x, y;
  double   *Xlimp, *Ylimp;
                                                         /** rho : origin on axis **/
  size_t    rhondx;
  double   *rhondp, rho;
  int       izh, cmp;
  int       irEl, i;
  int      *irGeoip;
  double    rEl,  ticksize = 0.08;
/*  static    char      prognamp[] = "gm.draw.E::gmDrawEGetXYEltrdAxis"; */

  cSetPVecRealloc(csetvp, 5, 1);  csetvp->x = 0;
  cSegPVecRealloc(&csegvec2, 10, 1);
  cCirclePVecRealloc(&cflevelv, 40, 2);

  gc1.lineWidth = 1;
  gc2.lineWidth = 3;

  celtrd.gcp = &gc2;
  csegvec2.gcp = &gc1;
  cflevel.gcp = &gc1;

  zndx = gridZp_mm->zndx;
  zndp = gridZp_mm->zndp;
  rhondx =gridRp_mm->rhotndx;
  rhondp =gridRp_mm->rhondp;

  irGeoip = gridZp_mm->irGeond_firstp;
  irEl   = gridRp_mm->eltrdRhondI;
  rEl    = gridRp_mm->eltrdR;

  cmp = classDbl(&izh, zh, zndp, zndx);
                                             /** electrode contour at altitude zh **/
  celtrd.xp[0] = 0.0;  celtrd.xp[1] = 0.0;
  rho = gmGridZGeoZ2R(zh, gridZp_mm);
  celtrd.r = rho;
/* fPrintF(stderr,"%s  zh=%f, izh=%d, radius rho=%f \n", prognamp, zh,izh,rho); */
/*  celtrd.anglp[0] = -90.0;  celtrd.anglp[1] = +90.0; */
  celtrd.anglp[0] = -180.0;  celtrd.anglp[1] = +180.0;
  cSetVecInc1ptype(csetvp, &celtrd, MY_CCIRCLE);

  Xlimp = gmDrawStateEGetLimX();
  Ylimp = gmDrawStateEGetLimY();
  csegvec2.x = 0;
                                                              /** y axis at x = 0 **/
  x = 0.0;
  cSegVecInc1CoorData(&csegvec2, x, Ylimp[0], x, Ylimp[1]);
                                                                      /** y ticks **/
  for(i = (int)Ylimp[0];  i <= (int)Ylimp[1];  i++)
  { y = i * 1.0;
    cSegVecInc1CoorData(&csegvec2, x, y, x - ticksize, y);
  }
                                                              /** x axis at y = 0 **/
  y = 0.0 ;
  cSegVecInc1CoorData(&csegvec2, Xlimp[0], y, Xlimp[1], y);
                                                                      /** x ticks **/
  for(i = (int)Xlimp[0];  i <= (int)Xlimp[1];  i++)
  { x = i * 1.0;
    cSegVecInc1CoorData(&csegvec2, x, y, x, y - ticksize);
  }
  cSetVecInc1ptype(csetvp, &csegvec2, MY_CSEGV);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawEComputeRZLevels()                                                  */
/*  Electrode frame: drawing of field level curve in RZ plane                 */
/******************************************************************************/
void      gmDrawEComputeRZLevels(cSetVec *csetvp, gmGridZ *gridZp_mm,
                                   gmGridR *gridRp_mm, double **fpp, lcZrange *lczrp)
{ static    lcSegVec  lcsegvec;
  static    curveGC   gc;
  static    intVec    irGeofV = {0,0,NULL};
  size_t    zndx;
  double   *zndp;
  size_t    rhondx;
  double   *rhondp;                                            /** origin on axis **/
  int      *irGeoip, *irGeofp;
  int       iz;
  myBOOL    bRect = 0;
  short     boolPrint = 1;
  static char    form1p[] = "     %s: NUMBER OF LEVEL SEG lcsegvec.x=%d \n";
  static char    prognamp[] = "gm.draw.E::gmDrawEComputeRZLevels";

  
  lcSegPVecRealloc(&lcsegvec, 100, 10);
  lcsegvec.gcp = &gc;
  gc.lineWidth = 1;

  zndx = gridZp_mm->zndx;
  zndp = gridZp_mm->zndp;
  rhondx =gridRp_mm->rhotndx;
  rhondp =gridRp_mm->rhondp;

  intPVecRealloc(&irGeofV, zndx, 10);
  irGeofp = irGeofV.p;
  for(iz = 0;  iz < zndx;  iz++){ irGeofp[iz] = rhondx -1;}
  irGeofV.x = zndx;
  irGeoip = gridZp_mm->irGeond_firstp;

  lcsegvec.x = 0;
  lcSegAddM(rhondp,rhondx, zndp,zndx, bRect,irGeoip,irGeofp, fpp, lczrp, &lcsegvec);
   if(boolPrint) fPrintF(stderr, form1p, prognamp, lcsegvec.x);
  
  fprintf(stderr, "gmDrawEComputeRZLevels*************1\n");
  cSetVecPrint(stderr, csetvp);
  csetvp->x = 0;
  
  cSetVecInc1ptype(csetvp, &lcsegvec, MY_LCSEGV);
  fprintf(stderr, "gmDrawEComputeRZLevels*************2\n");
  cSetVecPrint(stderr, csetvp);
  
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawEComputeXYLevels()                                                  */
/*  Electrode frame: drawing of field level curve in XY plane altitude zh     */
/******************************************************************************/
void      gmDrawEComputeXYLevels(cSetVec *csetvp, gmGridZ *gridZp_mm, 
                        gmGridR *gridRp_mm, double **fpp, lcZrange *lczrp, double zh)
{ static    cCircleVec  cflevelv = {0,0,NULL,NULL};
  static    curveGC   gc;
  static    dblVec    fzhV = {0,0,NULL};     /** interpollated field values at zh **/
  double    centerp[2] = {0.0, 0.0};
  double    anglp[2]   = {-180.0, +180.0};
  size_t    zndx;
  double   *zndp;
  size_t    rhondx;
  double   *rhondp;
  int      *irGeoip;
  double    rho;
  double    wei1, wei2, weir0, weir1;
  int       ir, irgeoi, izh, cmp;
  int       ifcn, ifcx, if0, if1, i, ifn, ifx;
  double   *fp1, *fp2, f0, f1, flev0 = 0.0, f, fdif, df;
  short     boolPrint = 1;
  static char    form1p[] = "     %s: NUMBER OF LEVEL CIRCLE cflevelv.x=%d \n";
  static char    form2p[] = "     %s: zh=%f, izh=%d \n";
  static char    form3p[] = "     %s  flev0=%f, df=%f, ifn=%d, ifx=%d\n";
  static char    form4p[] = "     %s  zh=%f, izh=%d, ir=%d, if0=%d, if1=%d, f0=%f, f1=%f\n";
  static char    prognamp[] = "gm.draw.E::gmDrawEComputeXYLevels";

  cCirclePVecRealloc(&cflevelv, 40, 4);  cflevelv.x = 0;
  gc.lineWidth = 1;
  cflevelv.gcp = &gc;

  zndx = gridZp_mm->zndx;
  zndp = gridZp_mm->zndp;
  irGeoip = gridZp_mm->irGeond_firstp;
  rhondx =gridRp_mm->rhotndx;
  rhondp =gridRp_mm->rhondp;

  cmp = classDbl(&izh, zh, zndp, zndx);        /** zh -> izh, among the zndp list **/
fPrintF(stderr, form2p, prognamp, zh, izh);
  wei1 = (zndp[izh +1]- zh)/(zndp[izh +1] - zndp[izh]);
  wei2 = (zh -   zndp[izh])/(zndp[izh +1] - zndp[izh]);
  irgeoi = irGeoip[izh];

  dblPVecRealloc(&fzhV, zndx, 2);
  fp1 = *(fpp + izh);
  fp2 = *(fpp + izh +1);
                                  /** compute z-interpollated values for each rho **/
  fzhV.x = 0;
  for(ir = 0;  ir < rhondx;  ir++)
  { f = *fp1++ * wei1 + *fp2++ * wei2;
    dblVecInc1(&fzhV, f);
/*
fPrintF(stderr,"%s zh=%f, izh=%d, ir=%d, f=%d \n", prognamp, zh, izh, ir, f);
*/
  }

  flev0 = lczrp->zo;
  df  = lczrp->zpa;
  ifn = lczrp->in;
  ifx = lczrp->ix;
if(boolPrint) fPrintF(stderr, form3p, prognamp, flev0, df, ifn,ifx);

  fp1 = fzhV.p + irgeoi;                      /** shift to avoid electrode inside **/
  f1 = *fp1++;
  fdif = (f1 - flev0)/df;
  if1 = (int)fdif;  if(fdif < 0.0) if1--;
  for(ir = 1 + irgeoi;  ir < rhondx;  ir++)
  { f0 = f1;  if0 = if1;
    f1 = *fp1++;
    fdif = (f1 - flev0)/df;
    if1 = (int)fdif;  if(fdif < 0.0) if1--;
if(boolPrint) fPrintF(stderr, form4p, prognamp, zh, izh, ir, if0,if1, f0,f1);
    if(if0 == if1) continue;
    if(if0 >= if1){ ifcx = if0;  ifcn = if1;}
    else          { ifcx = if1;  ifcn = if0;}
    if(ifcn > ifx) continue;
    if(ifcx < ifn) continue;
    ifcn++;
    if(ifcx > ifx) ifcx = ifx;
    if(ifcn < ifn) ifcn = ifn;
/*
fPrintF(stderr,"%s zh=%f, izh=%d, ir=%d, ifcn=%d, ifcx=%d, f0=%f, f1=%f\n",
                                            prognamp, zh, izh, ir, ifcn,ifcx, f0,f1);
*/
    for(i = ifcn;  i <= ifcx;  i++)
    { f = flev0 + i * df;
      weir0 = (f1 - f)/(f1 - f0);
      weir1 = (f - f0)/(f1 - f0);
      rho = rhondp[ir -1] * weir0 + rhondp[ir] * weir1;
      cCircleVecInc1Data(&cflevelv, centerp, rho, anglp);
    }
  }
if(boolPrint) fPrintF(stderr, form1p, prognamp, cflevelv.x);
  csetvp->x = 0;
  cSetVecInc1ptype(csetvp, &cflevelv, MY_CCIRCLEV);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawEGetRZCStrTitle(char **titlespp, int vei)                           */
/*                                                                            */
/* titlesp[0] = patient code                                                  */
/*        [1]   string pin mode                                               */
/*        [2]   string voltage                                                */
/*        [3]   string impedance                                              */
/*        [4]   psi rotation                                                  */
/*        [5]   zh                                                            */
/******************************************************************************/
cStrVec  *gmDrawEGetRZCStrTitle(char **titlespp, int vei, char **Etitlespp)
{ static    cStrVec   cstrv;
  static    cStr      cstrp[10];
  static    curveGC   gc1;
  static    char     *grveistrpp[3] = {"V" , "E", "I"};
  static    chrVec    titleZHV  = {0,0,NULL};
  static    chrVec    titlePSIV = {0,0,NULL};
  double    fontscale1 = 0.2;
  double    xposbaseL = -2.5, xposbaseR = 1.0;
  cStr     *curstrp;
  static    char      form1p[] = "     %s : %s\n";
  static    char      prognamp[] = "gm.draw.E::gmDrawEGetRZCStrTitle";

  chrPVecRealloc(&titleZHV, 10, 2);
  titleZHV.x = 0;
  chrVecIncStr(&titleZHV, "zh = ");  titleZHV.x-- ;
  chrVecIncStr(&titleZHV, Etitlespp[5]);
fPrintF(stderr, form1p, prognamp, titleZHV.p);

  chrPVecRealloc(&titlePSIV, 10, 2);
  titlePSIV.x = 0;
  chrVecIncStr(&titlePSIV, "psi = ");  titlePSIV.x-- ;
  chrVecIncStr(&titlePSIV, Etitlespp[4]);
fPrintF(stderr, form1p, prognamp, titlePSIV.p);

  cstrv.p = cstrp;  cstrv.z = 8;  cstrv.x = 0;
  cstrv.gcp = &gc1;
  gc1.lineWidth = 2;
                                                                    /** Font info **/
  cstrv.fscalp[0] = fontscale1;  cstrv.fscalp[1] = fontscale1;
  cstrv.findex = 0;
                                                                 /** Patient Code **/
  curstrp = cstrp;
  curstrp->strp = titlespp[0];
  curstrp->xyp[0] = xposbaseL;  curstrp->xyp[1] = 5.6;
  cstrv.x++;
                                                                     /** Pin mode **/
  curstrp++;
  curstrp->strp = titlespp[1];
  curstrp->xyp[0] = xposbaseL;  curstrp->xyp[1] = 5.2;
  cstrv.x++;
                                                                      /** Voltage **/
  curstrp++;
  curstrp->strp = titlespp[2];
  curstrp->xyp[0] = xposbaseL + 1.;  curstrp->xyp[1] = 5.2;
  cstrv.x++;

                                                                   /** VEI choice **/
  curstrp++;
  curstrp->strp = grveistrpp[vei];
  curstrp->xyp[0] = xposbaseR;  curstrp->xyp[1] = 5.2;
  cstrv.x++;
                                                                    /** Impedance **/
  if(vei == 2)
  { curstrp++;
    curstrp->strp = titlespp[3];
    curstrp->xyp[0] = xposbaseR + 0.5;  curstrp->xyp[1] = 5.2;
    cstrv.x++;
  }
                                                            /** ZH and PSI values **/
  curstrp++;
  curstrp->strp = titlePSIV.p;
  curstrp->xyp[0] = xposbaseL;  curstrp->xyp[1] = 4.8;
  cstrv.x++;

  curstrp++;
  curstrp->strp = titleZHV.p;
  curstrp->xyp[0] = xposbaseR;  curstrp->xyp[1] = 4.8;
  cstrv.x++;

  return  &cstrv;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawEGetXYCStrTitle(char **titlespp, int vei)                           */
/*                                                                            */
/* titlesp[0] = patient code                                                  */
/*        [1]   string pin mode                                               */
/*        [2]   string voltage                                                */
/*        [3]   string impedance                                              */
/*        [4]   psi rotation                                                  */
/*        [5]   zh                                                            */
/******************************************************************************/
cStrVec  *gmDrawEGetXYCStrTitle(char **titlespp, int vei, char **Etitlespp)
{ static    cStrVec   cstrv;
  static    cStr      cstrp[10];
  static    curveGC   gc1;
  double    fontscale1 = 0.2;
  double    xposbaseL = -2.5, xposbaseR = 0.5;
  cStr     *curstrp;
  double   *Xlimp, *Ylimp;
  static    char     *grveistrpp[3] = {"V" , "E", "I"};
  static    chrVec    titleZHV  = {0,0,NULL};
  static    chrVec    titlePSIV = {0,0,NULL};
  static    char      form1p[] = "     %s : %s\n";
  static    char      prognamp[] = "gm.draw::gmDrawGetEXYCStrTitle";

  chrPVecRealloc(&titleZHV, 10, 2);
  titleZHV.x = 0;
  chrVecIncStr(&titleZHV, "zh = ");  titleZHV.x-- ;
  chrVecIncStr(&titleZHV, Etitlespp[5]);
fPrintF(stderr, form1p, prognamp, titleZHV.p);

  chrPVecRealloc(&titlePSIV, 10, 2);
  titlePSIV.x = 0;
  chrVecIncStr(&titlePSIV, "psi = ");  titlePSIV.x-- ;
  chrVecIncStr(&titlePSIV, Etitlespp[4]);
fPrintF(stderr, form1p, prognamp, titlePSIV.p);

  cstrv.p = cstrp;  cstrv.z = 8;  cstrv.x = 0;
  cstrv.gcp = &gc1;
  gc1.lineWidth = 2;
                                                                    /** Font info **/
  cstrv.fscalp[0] = fontscale1;  cstrv.fscalp[1] = fontscale1;
  cstrv.findex = 0;

  Xlimp = gmDrawStateEGetLimX();
  Ylimp = gmDrawStateEGetLimY();
  xposbaseL = Xlimp[0] + 2. ;
  xposbaseR = 0.5;
                                                                 /** Patient Code **/
  curstrp = cstrp;
  curstrp->strp = titlespp[0];
  curstrp->xyp[0] = xposbaseL;  curstrp->xyp[1] = Ylimp[1] - 2. * fontscale1;
  cstrv.x++;

                                                                     /** Pin mode **/
  curstrp++;
  curstrp->strp = titlespp[1];
  curstrp->xyp[0] = xposbaseL;  curstrp->xyp[1] = Ylimp[1] - 4. * fontscale1;
  cstrv.x++;
                                                                      /** Voltage **/
  curstrp++;
  curstrp->strp = titlespp[2];
  curstrp->xyp[0] = xposbaseL + 1.0;  curstrp->xyp[1] = Ylimp[1] - 4. * fontscale1;
  cstrv.x++;
                                                                   /** VEI choice **/
  curstrp++;
  curstrp->strp = grveistrpp[vei];
  curstrp->xyp[0] = xposbaseR;  curstrp->xyp[1] = Ylimp[1] - 4. * fontscale1;
  cstrv.x++;
                                                                    /** Impedance **/
  if(vei == 2)
  { curstrp++;
    curstrp->strp = titlespp[3];
    curstrp->xyp[0] = xposbaseR + 0.5;  curstrp->xyp[1] = Ylimp[1] - 4. * fontscale1;
    cstrv.x++;
  }

                                                            /** ZH and PSI values **/
  curstrp++;
  curstrp->strp = titleZHV.p;
  curstrp->xyp[0] = xposbaseL;  curstrp->xyp[1] = Ylimp[1] - 6. * fontscale1;
  cstrv.x++;

  curstrp++;
  curstrp->strp = titlePSIV.p;
  curstrp->xyp[0] = xposbaseR;  curstrp->xyp[1] = Ylimp[1] - 6. * fontscale1;
  cstrv.x++;

  return  &cstrv;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawEGetCStrPinIndices()                                                */
/*                                                                            */
/* titlesp[0] = patient code                                                  */
/*        [1]   string pin mode                                               */
/*        [2]   string voltage                                                */
/*        [3]   string impedance                                              */
/******************************************************************************/
cStrVec  *gmDrawEGetCStrPinIndices()
{ static    cStrVec   cstrvpinindex;
  static    cStr      cstrpinindexp[4];
  static    curveGC   gc2;
  static    char     *pinindexstrpp[4] = {"0", "1", "2", "3"};
  cStr     *curstrp;
/*  static    char      prognamp[] = "gm.draw::gmGetCStrPinIndices"; */

                                                        /** Pin Geometric Indices **/
  cstrvpinindex.p = cstrpinindexp;
  cstrvpinindex.z = 4;  cstrvpinindex.x = 0;
  cstrvpinindex.gcp = &gc2;
  gc2.lineWidth = 3;
                                                                    /** Font info **/
  cstrvpinindex.fscalp[0] = 0.4;  cstrvpinindex.fscalp[1] = 0.4;
  cstrvpinindex.findex = 0;
  curstrp = cstrpinindexp;
  curstrp->strp = pinindexstrpp[3];
  curstrp->xyp[0] = 0.15;  curstrp->xyp[1] = -3.2;
  cstrvpinindex.x++;

  curstrp++;
  curstrp->strp = pinindexstrpp[2];
  curstrp->xyp[0] = 0.15;  curstrp->xyp[1] = -1.2;
  cstrvpinindex.x++;

  curstrp++;
  curstrp->strp = pinindexstrpp[1];
  curstrp->xyp[0] = 0.15;  curstrp->xyp[1] = +0.8;
  cstrvpinindex.x++;

  curstrp++;
  curstrp->strp = pinindexstrpp[0];
  curstrp->xyp[0] = 0.15;  curstrp->xyp[1] = +2.8;
  cstrvpinindex.x++;

  return &cstrvpinindex;
}
/******************************************************************************/
/******************************************************************************/
