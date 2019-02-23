/*  URMAE/orientHaut/linear4.GL.V4/gm.draw.ST.c                               */
/*  Mennessier Gerard                 20010613                                */
/*  Last Revised : G.M.               20040210                                */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utiMath.constant.def.h"

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

#include  "gm.draw.ST.h"

#include  "gm.linear4.initsolve.glob.h"
#include  "pallidus.geom.glob.h"

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
#include  "gm.drawstate.ST.h"
#include  "gm.drawstate.E.h"   /** to get Electrode RZ level cylindrical surfaces **/
#include  "grid_pot.ST.linear4.h"

static double  deltaDeg, deltaRad, cdelta, sdelta,
               sdeltastheta, s2deltas2theta, cdeltastheta,
               sorientm, corientm, rorientm, dorientm, 
               discr1, sqdiscr1;
/******************************************************************************/
/*                                                                            */
/*  gmDrawSTGetRZeltrdAxis()                                                  */
/*  ScannerTranslated frame: drawing of electrode, pins and axis              */
/******************************************************************************/
void      gmDrawSTGetRZeltrdAxis(cSetVec *csetvp,
                                              gmGridZ *gridZp_mm, gmGridR *gridRp_mm)
{ static    cPoly     cpoly;
                                   /** csegvec1: pins;   csegvec2: z and rho axis **/
  static    cSegVec   csegvec1, csegvec2;
  static    cSegVec   csegvec3;      /** csegvec3: fixed radius part of electrode **/
  static    cEllipse  cellips1;                            /** electrode cylinder **/
  static    cEllipse  cellips2;               /** ellipsoidal electrode extremity **/
  static    curveGC   gc1, gc2;
  static    intVec    irGeofV = {0,0,NULL};
  size_t    zndx;
  double   *zndp;
  double    z, zext0, zext1, zextL;
  double    rEl, rEl2;
  double    z12prEl2, z12prEl2mzL2;
  double    discr2, discr3, discr4;
  double    sqdiscr2 = 0.0, sqdiscr3 = 0.0, sqdiscr4 = 0.0;

                        /** Boolean. If points Q exist, boolQ=1=T; else boolQ=0=F **/
  short     boolQ = 0;
            /** If intersection with ellipsoid non void boolE=1=T; else boolE=0=F **/
  short     boolE = 0;
  double    szetaC, czetaC, rzetaC = 0.0, dzetaC = 0.0;
  double    szetaE, czetaE, rzetaE = 0.0, dzetaE = 0.0;
  double    vpQ, vmQ;
  double    centerp[2], axisp[2], anglp[2];
  size_t    rhondx;
  double   *rhondp, rho;                                       /** origin on axis **/
  double   *znxp, *rhonxp;                                 /**  z and rho Elimits **/
  int       iz, iz0,iz1;
  int       ir, ir2, irEl;
  int      *irGeoip, *irGeofp;
  double    ticksize = 0.08;
  short     boolPrint = 1;
  static char    prognamp[] = "gm.draw.ST::gmDrawSTGetRZeltrdAxis";

if(boolPrint) fPrintF(stderr, "%s  BEGIN\n", prognamp);
  cSetPVecRealloc(csetvp, 5, 1);  csetvp->x = 0;

  cPolyPRealloc(&cpoly, 10, 1);
  cSegPVecRealloc(&csegvec1, 10, 1);
  cSegPVecRealloc(&csegvec2, 10, 1);
  cEllipseZero(&cellips1);  cEllipseZero(&cellips2);

  gc1.lineWidth = 1;
  gc2.lineWidth = 3;

  cpoly.gcp = &gc2;
  cellips1.gcp = &gc2;  cellips2.gcp = &gc2;
  csegvec3.gcp = &gc2;
  csegvec3.x = 0;

  csegvec1.gcp = &gc1;
  csegvec1.x = 0;
  csegvec2.gcp = &gc1;
  csegvec2.x = 0;

  zndx = gridZp_mm->zndx;
  zndp = gridZp_mm->zndp;
  rhondx =gridRp_mm->rhotndx;
  rhondp =gridRp_mm->rhondp;

  intPVecRealloc(&irGeofV, zndx, 10);
  irGeofp = irGeofV.p;
  for(iz=0;  iz < zndx; iz++){ irGeofp[iz] = rhondx-1;}
  irGeofV.x = zndx;
  irGeoip = gridZp_mm->irGeond_firstp;
  irEl   = gridRp_mm->eltrdRhondI;
  rEl    = gridRp_mm->eltrdR;
  rEl2 = rEl*rEl;

  iz0 = gridZp_mm->eltrdExtIp[0];  
  iz1 = gridZp_mm->eltrdExtIp[1];
  zext0 = gridZp_mm->eltrdExtZp[0];
  zext1 = gridZp_mm->eltrdExtZp[1];
  zextL = fabs(zext1 - zext0);
  z12prEl2 = zext1*zext1 + rEl2;
  z12prEl2mzL2 = z12prEl2 - zextL*zextL;

  gmDrawSTSetGeomRZ();

  discr2 = rEl2 - z12prEl2mzL2 * s2deltas2theta;
  discr3 = rEl2 - z12prEl2 * s2deltas2theta;
  discr4 = rEl2 + (zextL*zextL - rEl2) * s2deltas2theta;
  sqdiscr4 = sqrt(discr4);
  boolE = (discr2 >= 0.0)? 1 : 0;
  if(boolE)
  { sqdiscr2 = sqrt(discr2);
  }
  boolQ = (discr3 >= 0.0)? 1 : 0;
  if(boolQ)
  { szetaC = (zext1 * sdeltastheta)/(rEl * sqdiscr1);
    szetaC = fabs(szetaC);
    sqdiscr3 = sqrt(discr3);
    czetaC = sqdiscr3/(rEl * sqdiscr1);
    rzetaC = atan2(szetaC, czetaC);
    dzetaC = rzetaC * myRAD2DEG;
    dzetaC = 90. - dzetaC;
    vpQ = rEl * czetaC;
    vmQ = zext1/sqdiscr1;
    sqdiscr2 = sqrt(discr2);
    szetaE = sqdiscr3 * sqdiscr4/(rEl * sqdiscr1 * sqdiscr2);
    czetaE = zext1 * zextL * s2deltas2theta/(rEl * sqdiscr1 * sqdiscr2);
    rzetaE = atan2(szetaE, czetaE);
    dzetaE = rzetaE * myRAD2DEG;
  }
  if(boolPrint)
  { fPrintF(stderr,"     %s  deltaDeg=%f, boolE=%d, boolQ=%d, sTheta=%f, cTheta=%f, orient=%f, zetaC=%f\n",
                 prognamp, deltaDeg, boolE, boolQ, stheta, ctheta, dorientm, dzetaC);
  }
                                                            /** electrode contour **/
                                              /** fixed radius = cylindrical part **/
  { centerp[0] = 0.0;  centerp[1] = 0.0;
    axisp[0] = (fabs(sdeltastheta) > 0.0001)? rEl/fabs(sdeltastheta) : rEl/0.0001;
    axisp[1] = rEl;
    cEllipseSetCAO(&cellips1, centerp, axisp, dorientm);
    if(boolQ)
    { anglp[0] = dzetaC;
      anglp[1] = 360. - dzetaC;
    }
    else
    { anglp[0] = -180.0;
      anglp[1] =  180.0;
    }
    cEllipseSetDir(&cellips1, anglp);
    cSetVecInc1ptype(csetvp, &cellips1, MY_CELLIPSE);
  }
                                                        /** ellipsoidal extremity **/
  if(boolE)
  { centerp[0] = (zext1 * rEl2 / discr4) * cdeltastheta;
    centerp[1] = (zext1 * rEl2 / discr4) * ctheta;
    axisp[0] = rEl * zextL * sqdiscr2/discr4;
    axisp[1] = rEl * sqdiscr2/sqdiscr4;
    cEllipseSetCAO(&cellips2, centerp, axisp, dorientm);
    if(boolQ)
    { anglp[0] =   dzetaE;
      anglp[1] = - dzetaE;
    }
    else
    { anglp[0] = 0.0;
      anglp[1] = 0.0;
    }
    cEllipseSetDir(&cellips2, anglp);
    if(boolE) cSetVecInc1ptype(csetvp, &cellips2, MY_CELLIPSE);
  }
                                                                         /** pins **/
/*
  for(iz=0;  iz < 8;  iz++)
  { z = gridZp_mm->pinExtZp[iz];
    cSegVecInc1CoorData(&csegvec1, rhondp[0], z, rhondp[irEl], z);
  }
  cSetVecInc1ptype(csetvp, &csegvec1, MY_CSEGV);
*/
                                                                       /** z axis **/
  znxp = gmDrawStateSTGetLimZ();
  rho = 0.0;
  cSegVecInc1CoorData(&csegvec2, rho, znxp[0], rho, znxp[1]);
                                                                      /** z ticks **/
  iz0 = (int)znxp[0];  iz1 = (int)znxp[1];
  for(iz = iz0;  iz <= iz1;  iz++)
  { z = iz * 1.0;
    cSegVecInc1CoorData(&csegvec2, rho, z, rho - ticksize, z);
  }
                                                          /** x axis at z = 0. mm **/
  rhonxp = gmDrawStateSTGetLimRho();
  z = 0.0 ;
  cSegVecInc1CoorData(&csegvec2, rhonxp[0], z, rhonxp[1], z);
  ir2 = (int)rhonxp[1];
  for(ir = -ir2;  ir <= ir2;  ir++)
  { rho = ir * 1.0;
    cSegVecInc1CoorData(&csegvec2, rho, z, rho, z - ticksize);        /** x ticks **/
  }
  cSetVecInc1ptype(csetvp, &csegvec2, MY_CSEGV);

if(boolPrint) fPrintF(stderr, "     %s  Calling cEllipsePrint\n", prognamp);
if(boolPrint){ cEllipsePrint(stderr, &cellips1); cEllipsePrint(stderr, &cellips2);}
if(boolPrint) fPrintF(stderr, "%s  END\n", prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawSTGetRZtraceXYsection()                                             */
/*  ScannerTranslated frame:                                                  */
/*    drawing in RZplane of trace of the XYsection altitude zh                */
/*                                              (i.e. a single straight line) */
/******************************************************************************/
cSegVec  *gmDrawSTGetRZtraceXYsection(double zh)
{ static    cSeg      cseg;
  static    cSegVec   csegvec;
  static    curveGC   gc;
  double   *rhonxp;

  rhonxp = gmDrawStateSTGetLimRho();
  gc.lineWidth = 2;
  csegvec.gcp = &gc;
  csegvec.p = &cseg;
  csegvec.z = 1;
  csegvec.x = 0;
  cSegVecInc1CoorData(&csegvec, rhonxp[0], zh, rhonxp[1], zh);
  return  &csegvec;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawSTGetXYeltrdAxis()                                                  */
/*  ScannerTranslated frame: drawing of electrode and axis                    */
/******************************************************************************/
void      gmDrawSTGetXYeltrdAxis(cSetVec *csetvp,
                                   gmGridZ *gridZp_mm, gmGridR *gridRp_mm, double zh)
{ static    cSegVec   csegvec2;
  static    cEllipse  cellips1;                            /** electrode cylinder **/
  static    cEllipse  cellips2;               /** ellipsoidal electrode extremity **/
  static    curveGC   gc1, gc2;
  size_t    zndx;
  double   *zndp;
  double    zext0, zext1, zextL,  z1ct, zLct, rst, z1cmzh;
  double    zcrit1, zcrit2, zcrit3;
  double    dban, disc1, disc2, disc3;
  double    sqdisc1 = 0.0, sqdisc2 = 0.0, sqdisc3 = 0.0;
                                 /** If points P exist, boolP=1=T; else boolP=0=F **/
  short     boolP = 0;
            /** If intersection with ellipsoid non void boolE=1=T; else boolE=0=F **/
  short     boolE = 0; 
  double    calphaP, salphaP, ralphaP, dalphaP = 0.0;
  double    x, y;
  double   *Xlimp, *Ylimp;
                                                         /** rho : origin on axis **/
  size_t    rhondx;
  double   *rhondp, rho;
  int       izh, cmp;
  int       irEl, i;
  int      *irGeoip;
  double    rEl, rEl2, ticksize = 0.08;
  short     boolPrint = 1;
  static    char      prognamp[] = "gm.draw.ST::gmDrawSTGetXYeltrdAxis";

if(boolPrint) fPrintF(stderr, "%s  BEGIN\n", prognamp);
  cSetPVecRealloc(csetvp, 5, 1);  csetvp->x = 0;
  cSegPVecRealloc(&csegvec2, 10, 1);

  gc1.lineWidth = 1;
  gc2.lineWidth = 3;

  cellips1.gcp = &gc2;  cellips2.gcp = &gc2;

  csegvec2.gcp = &gc1;
  csegvec2.x = 0;

  zndx = gridZp_mm->zndx;
  zndp = gridZp_mm->zndp;
  rhondx =gridRp_mm->rhotndx;
  rhondp =gridRp_mm->rhondp;

  irGeoip = gridZp_mm->irGeond_firstp;

  zext0 = gridZp_mm->eltrdExtZp[0];
  zext1 = gridZp_mm->eltrdExtZp[1];
  zextL = fabs(zext1 - zext0);
  z1ct  = zext1 * ctheta;
  zLct  = zextL * ctheta;
  z1cmzh = z1ct - zh;

  irEl = gridRp_mm->eltrdRhondI;
  rEl  = gridRp_mm->eltrdR;
  rEl2 = rEl*rEl;
  rst  = rEl * stheta;

  disc1 = rst*rst - z1cmzh*z1cmzh;
  boolP = (fabs(z1cmzh) < rst)? 1 : 0;
  if(boolP)
  { sqdisc1 = sqrt(disc1);
    calphaP = z1cmzh/rst;
    salphaP = sqdisc1/rst;
    ralphaP = atan2(salphaP, calphaP);
    dalphaP = ralphaP * myRAD2DEG;
  }

  disc2 = zLct*zLct + rst*rst;
  sqdisc2 = sqrt(disc2);

  disc3 = disc2 - z1cmzh*z1cmzh;
  boolE = (fabs(z1cmzh) < sqdisc2)? 1 : 0;
  if(boolE)
  { sqdisc3 = sqrt(disc3);
  }

  zcrit1 = z1ct - rst;
  zcrit2 = z1ct + rst;
  zcrit3 = z1ct + sqdisc2;

  cmp = classDbl(&izh, zh, zndp, zndx);
  rho = 0.0;
                                             /** electrode contour at altitude zh **/
                                                   /** Intersection with cylinder **/
  if(zh <= zcrit2)
  { rho = rEl;
    cellips1.xp[0] = eltrdSVecp[0] *(zh/ctheta);
    cellips1.xp[1] = eltrdSVecp[1] *(zh/ctheta);
    cellips1.orient = phiDeg;
    cellips1.axisp[0] = rEl/ctheta;
    cellips1.axisp[1] = rEl;
    if(zh <= zcrit1)                                            /** ONLY cylinder **/
    { cellips1.anglDp[0] = -180.0;  cellips1.anglDp[1] = 180.0;
    }
    else                                             /** limit cylinder/ellipsoid **/
    { cellips1.anglDp[0] = dalphaP;  cellips1.anglDp[1] = 360.0 - dalphaP;
    }
    cSetVecInc1ptype(csetvp, &cellips1, MY_CELLIPSE);
  }

                                                  /** Intersection with ellipsoid **/
  if(boolE && zh >= zcrit1)
  { dban = ((zextL*zextL - rEl2)* zh * ctheta + zext1 * rEl2) * stheta/disc2;
    cellips2.xp[0] = dban * cphi;
    cellips2.xp[1] = dban * sphi;
    cellips2.orient = phiDeg;
    dban = rEl * sqdisc3;
    cellips2.axisp[0] = dban * zextL/disc2;
    cellips2.axisp[1] = dban/sqdisc2;
    if(zh <= zcrit2)                                 /** limit cylinder/ellipsoid **/
    { cellips2.anglDp[0] = -dalphaP;  cellips2.anglDp[1] = dalphaP;
    }
    else                                                       /** ONLY ellipsoid **/
    { cellips2.anglDp[0] = -180.0;  cellips2.anglDp[1] = 180.0;
    }
    cSetVecInc1ptype(csetvp, &cellips2, MY_CELLIPSE);
  }

  Xlimp = gmDrawStateSTGetLimX();
  Ylimp = gmDrawStateSTGetLimY();
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

if(boolPrint) fPrintF(stderr, "     %s  Calling cEllipsePrint\n", prognamp);
if(boolPrint){ cEllipsePrint(stderr, &cellips1); cEllipsePrint(stderr, &cellips2);}
if(boolPrint) fPrintF(stderr, "%s  END\n", prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawSTGetXYtraceRZsection()                                             */
/*  ScannerTranslated frame:                                                  */
/*    drawing in XYplane of trace of the RZsection angle psi                  */
/*                                         (i.e. a single half straight line) */
/******************************************************************************/
cSegVec  *gmDrawSTGetXYtraceRZsection(double psi)
{ static    cSeg      cseg;
  static    cSegVec   csegvec;
  static    curveGC   gc;
  double   *Xlimp, *Ylimp, rhon, rhox, cpsi, spsi;

  cpsi = cos(psi);  spsi = sin(psi);
  Xlimp = gmDrawStateSTGetLimX();
  Ylimp = gmDrawStateSTGetLimY();
  rhon = sqrt(Xlimp[0]*Xlimp[0] + Ylimp[0]*Ylimp[0]);
  rhox = sqrt(Xlimp[1]*Xlimp[1] + Ylimp[1]*Ylimp[1]);
  gc.lineWidth = 2;
  csegvec.gcp = &gc;
  csegvec.p = &cseg;
  csegvec.z = 1;
  csegvec.x = 0;
  rhon = 0.0;
  cSegVecInc1CoorData(&csegvec,-rhon*cpsi, -rhon*spsi, rhox*cpsi, rhox*spsi);
  return  &csegvec;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawSTComputeRZLevels()                                                 */
/*  ScannerTranslated frame: drawing of field level curve in RZ plane         */
/******************************************************************************/
void      gmDrawSTComputeRZLevels(cSetVec *csetvp)
{ static    lcSegVec  lcsegV1, lcsegV2;
  static    curveGC   gc;
  cSetVec  *csetvERZp;
  cSet     *csetERZp;
  lcSegVec *lcsegvERZp;
  lcSeg    *lcsegip, *lcsegf1p, *lcsegf2p;
  size_t    nsegx;
  int       i, nf1 =0, nf2 = 0;
  double    rhoE, zE, discr3, sqdiscr3;
  double    v0m, v0p, v1m, v1p;
  double    sv0m = 0.0, cv0m = 0.0, sv0p = 0.0, cv0p = 0.0;
  double    sv1m = 0.0, cv1m = 0.0, sv1p = 0.0, cv1p = 0.0;
  double   *xip, *yip, *xf1p, *yf1p, *xf2p, *yf2p;
  short     exist0, exist1;
  static char progmp[]="gm.draw.ST::RZlevel";
  
  
  fprintf(stderr, "%s ***********************RZ 0********\n", progmp);
  
  csetvERZp  = gmDrawStateEGetRZlevelCSetp();
  
  cSetVecPrint(stderr, csetvERZp);
  fprintf(stderr, "%s  %d *********************** csetvERZp RZ 1********\n", progmp,csetvERZp);
  
  csetERZp   = csetvERZp->p;
  
  cSetPrint(stderr, csetERZp); 
  fprintf(stderr, "%s  %d *********************** csetERZp  RZ 2 x********\n", progmp,csetERZp);
  
  
  lcsegvERZp = (lcSegVec *)csetERZp->p;
  lcSegVecPrint(stderr, lcsegvERZp);
  
  fprintf(stderr, "%s  %d ***********************lcsegvERZp RZ 3********\n", progmp,lcsegvERZp);

  
  nsegx = lcsegvERZp->x;
  fprintf(stderr, "%s  %d ***********************nsegx RZ 4********\n", progmp,nsegx);

  
  lcsegip = lcsegvERZp->p;

  gc.lineWidth = 1;

  lcSegPVecRealloc(&lcsegV1, nsegx, 10);
  lcsegV1.x = 0;
  lcsegV1.gcp = &gc;
  lcsegf1p = lcsegV1.p;
  lcSegPVecRealloc(&lcsegV2, nsegx, 10);
  lcsegV2.x = 0;
  lcsegV2.gcp = &gc;
  lcsegf2p = lcsegV2.p;
   
  fprintf(stderr, "%s ***********************RZ 1********\n", progmp);
  
  gmDrawSTSetGeomRZ();
  
  fprintf(stderr, "%s ***********************RZ 2********\n", progmp);
  
  for(i = 0;  i < nsegx;  i++, lcsegip++)                 /** Level Segments loop **/
  { xip = lcsegip->x;  yip = lcsegip->y;
                                                          /** segment extremity 0 **/
    rhoE = *xip++;
    zE   = *yip++;
    discr3 = rhoE*rhoE*discr1 - zE*zE*s2deltas2theta;
    if(discr3 >= 0.0)
    { exist0 = 1;
      sqdiscr3 = sqrt(discr3);
      v0m = zE/sqdiscr1;
      v0p = sqdiscr3/sqdiscr1;
      cv0m = v0m * corientm;
      sv0m = v0m * sorientm;
      cv0p = v0p * corientm;
      sv0p = v0p * sorientm;
    }
    else
    { exist0 = 0;  v0m = v0p = 0.0;
    }
                                                          /** segment extremity 1 **/
    rhoE = *xip;
    zE   = *yip;
    discr3 = rhoE*rhoE*discr1 - zE*zE*s2deltas2theta;
    exist1 = (discr3 >= 0.0);
    if(discr3 >= 0.0)
    { exist1 = 1;
      sqdiscr3 = sqrt(discr3);
      v1m = zE/sqdiscr1;
      v1p = sqdiscr3/sqdiscr1;
      cv1m = v1m * corientm;
      sv1m = v1m * sorientm;
      cv1p = v1p * corientm;
      sv1p = v1p * sorientm;
    }
    else
    { exist1 = 0;  v1m = v1p = 0.0;
    }
                   /** if both circle intersect store segment and its symmetrical **/
    if(exist0 && exist1)
    { lcsegf1p->iz = lcsegip->iz;
      xf1p = lcsegf1p->x;  yf1p = lcsegf1p->y;
      lcsegf2p->iz = lcsegip->iz;                                 /** symmetrical **/
      xf2p = lcsegf2p->x;  yf2p = lcsegf2p->y;

                                                          /** segment extremity 0 **/
      *xf1p++ = cv0m + sv0p;
      *yf1p++ = sv0m - cv0p;
                                                          /** segment extremity 0 **/
      *xf2p++ = cv0m - sv0p;
      *yf2p++ = sv0m + cv0p;
                                                          /** segment extremity 1 **/
      *xf1p   = cv1m + sv1p;
      *yf1p   = sv1m - cv1p;

      *xf2p   = cv1m - sv1p;
      *yf2p   = sv1m + cv1p;
      nf1++;
      nf2++;
      lcsegf1p++;
      lcsegf2p++;
    }
                 /** if only 1 circle intersect, add a symmetric segment in set 1 **/
    else if(exist0)
    { lcsegf1p->iz = lcsegip->iz;
      xf1p = lcsegf1p->x;  yf1p = lcsegf1p->y;
      *xf1p++ = cv0m + sv0p;
      *yf1p++ = sv0m - cv0p;
      *xf1p   = cv0m - sv0p;
      *yf1p   = sv0m + cv0p;
      nf1++;
      lcsegf1p++;
    }
    else if(exist1)
    { lcsegf1p->iz = lcsegip->iz;
      xf1p = lcsegf1p->x;  yf1p = lcsegf1p->y;
      *xf1p++ = cv1m + sv1p;
      *yf1p++ = sv1m - cv1p;
      *xf1p   = cv1m - sv1p;
      *yf1p   = sv1m + cv1p;
      nf1++;
      lcsegf1p++;
    }
  }
  lcsegV1.x = nf1;  lcsegV2.x = nf2;

  csetvp->x = 0;
  cSetVecInc1ptype(csetvp, &lcsegV1, MY_LCSEGV);
  cSetVecInc1ptype(csetvp, &lcsegV2, MY_LCSEGV);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawSTComputeXYLevels()                                                 */
/*  ScannerTranslated frame:                                                  */
/*                   drawing of field level curve in XY plane altitude zh     */
/******************************************************************************/
void      gmDrawSTComputeXYLevels(cSetVec *csetvp, lcZrange *lczrp, gmGridPotSTXY *p)
{ static    lcSegVec  lcsegvec;
  lcSeg    *lcsegp;
  static    curveGC   gc;
  static    intVec    aiV = {0,0,NULL},  afV = {0,0,NULL};
  size_t    alphax;
  double   *alphap;
  size_t    saxx;
  double   *saxp;
  int      *saxip, *saxfp, i;
  double  **fpp;
  double    xyzp[3], *ap, *sap, zst;
  myBOOL    brect = 1;
  static    char      prognamp[] = "gm.draw.ST::gmDrawSTComputeXYLevels";

  lcSegPVecRealloc(&lcsegvec, 100, 10);
  lcsegvec.gcp = &gc;
  gc.lineWidth = 1;

  alphax = p->alphax;
  alphap = p->alphap;
  saxx = p->saxx;
  saxp = p->saxp;
  zst = p->zst;
  fpp = p->potpp;

  saxip = aiV.p;  saxfp = afV.p;     /** NOT really used because rectangular grid **/

  lcsegvec.x = 0;
  lcSegAddM(alphap,alphax, saxp,saxx, brect,saxip,saxfp, fpp, lczrp, &lcsegvec);

fPrintF(stdout, "%s  ST (saxis,alpha) level curves\n", prognamp);
/* lcSegVecPrint(stdout, &lcsegvec); */

                            /** from  (saxis,alpha) coordinate to (X,Y) ST frame  **/
  lcsegp = lcsegvec.p;
  for(i = 0;  i < lcsegvec.x;  i++, lcsegp++)
  { ap = lcsegp->x;  sap = lcsegp->y;  
    gmGridPotSTalphaSaxisZ2xy(xyzp, *ap, *sap, zst);
    *ap = xyzp[0];  *sap = xyzp[1];
    ap++;  sap++;
    gmGridPotSTalphaSaxisZ2xy(xyzp, *ap, *sap, zst);
    *ap = xyzp[0];  *sap = xyzp[1];
  }
fPrintF(stdout, "%s  ST (X,Y) level curves\n", prognamp);
/* lcSegVecPrint(stdout, &lcsegvec); */

  csetvp->x = 0;
  cSetVecInc1ptype(csetvp, &lcsegvec, MY_LCSEGV);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawSTGetRZCStrTitle(char **titlespp, int vei)                          */
/*                                                                            */
/* titlesp[0] = patient code                                                  */
/*        [1]   string pin mode                                               */
/*        [2]   string voltage                                                */
/*        [3]   string impedance                                              */
/*        [4]   psi rotation                                                  */
/*        [5]   zh                                                            */
/******************************************************************************/
cStrVec  *gmDrawSTGetRZCStrTitle(char **titlespp, int vei, char **STtitlespp)
{ static    cStrVec   cstrv;
  static    cStr      cstrp[10];
  static    curveGC   gc1;
  double    fontscale1 = 0.2;
  double    xposbaseL = -2.5, xposbaseR = 1.0;
  cStr     *curstrp;
  static    char     *grveistrpp[3] = {"V" , "E", "I"};
  static    chrVec    titleZHV  = {0,0,NULL};
  static    chrVec    titlePSIV = {0,0,NULL};
  static    char      prognamp[] = "gm.draw.ST::gmDrawSTGetRZCStrTitle";

  chrPVecRealloc(&titleZHV, 10, 2);
  titleZHV.x = 0;
  chrVecIncStr(&titleZHV, "zh = ");  titleZHV.x-- ;
  chrVecIncStr(&titleZHV, STtitlespp[5]);
fPrintF(stderr, "     %s : %s\n", prognamp, titleZHV.p);

  chrPVecRealloc(&titlePSIV, 10, 2);
  titlePSIV.x = 0;
  chrVecIncStr(&titlePSIV, "psi = ");  titlePSIV.x-- ;
  chrVecIncStr(&titlePSIV, STtitlespp[4]);
fPrintF(stderr, "     %s : %s\n", prognamp, titlePSIV.p);

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
/*  gmDrawSTGetXYCStrTitle(char **titlespp, int vei)                          */
/*                                                                            */
/* titlesp[0] = patient code                                                  */
/*        [1]   string pin mode                                               */
/*        [2]   string voltage                                                */
/*        [3]   string impedance                                              */
/*        [4]   psi rotation                                                  */
/*        [5]   zh                                                            */
/******************************************************************************/
cStrVec  *gmDrawSTGetXYCStrTitle(char **titlespp, int vei, char **STtitlespp)
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
  static    char      prognamp[] = "gm.draw.ST::gmDrawSTGetXYCStrTitle";

  chrPVecRealloc(&titleZHV, 10, 2);
  titleZHV.x = 0;
  chrVecIncStr(&titleZHV, "zh = ");  titleZHV.x-- ;
  chrVecIncStr(&titleZHV, STtitlespp[5]);
fPrintF(stderr, "     %s : %s\n", prognamp, titleZHV.p);

  chrPVecRealloc(&titlePSIV, 10, 2);
  titlePSIV.x = 0;
  chrVecIncStr(&titlePSIV, "psi = ");  titlePSIV.x-- ;
  chrVecIncStr(&titlePSIV, STtitlespp[4]);
fPrintF(stderr, "     %s : %s\n", prognamp, titlePSIV.p);

  cstrv.p = cstrp;  cstrv.z = 8;  cstrv.x = 0;
  cstrv.gcp = &gc1;
  gc1.lineWidth = 2;
                                                                    /** Font info **/
  cstrv.fscalp[0] = fontscale1;  cstrv.fscalp[1] = fontscale1;
  cstrv.findex = 0;

  Xlimp = gmDrawStateSTGetLimX();
  Ylimp = gmDrawStateSTGetLimY();
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
/*  gmDrawSTGetCStrPinIndices()                                               */
/*                                                                            */
/* titlesp[0] = patient code                                                  */
/*        [1]   string pin mode                                               */
/*        [2]   string voltage                                                */
/*        [3]   string impedance                                              */
/******************************************************************************/
cStrVec  *gmDrawSTGetCStrPinIndices()
{ static    cStrVec   cstrvpinindex;
  static    cStr      cstrpinindexp[4];
  static    curveGC   gc2;
  static    char     *pinindexstrpp[4] = {"0", "1", "2", "3"};
  cStr     *curstrp;
/*  static    char      prognamp[] = "gm.draw.ST::gmDrawSTGetCStrPinIndices"; */

                                                        /** Pin Geometric Indices **/
  cstrvpinindex.p = cstrpinindexp;
  cstrvpinindex.z = 4;  cstrvpinindex.x = 0;
  cstrvpinindex.gcp = &gc2;
  gc2.lineWidth = 3;
                                                                    /** Font info **/
  cstrvpinindex.fscalp[0] = 0.4;  cstrvpinindex.fscalp[1] = 0.4;
  cstrvpinindex.findex = 0;
  curstrp = cstrpinindexp;
  curstrp->strp = pinindexstrpp[0];
  curstrp->xyp[0] = 0.15;  curstrp->xyp[1] = -3.2;
  cstrvpinindex.x++;

  curstrp++;
  curstrp->strp = pinindexstrpp[1];
  curstrp->xyp[0] = 0.15;  curstrp->xyp[1] = -1.2;
  cstrvpinindex.x++;

  curstrp++;
  curstrp->strp = pinindexstrpp[2];
  curstrp->xyp[0] = 0.15;  curstrp->xyp[1] = +0.8;
  cstrvpinindex.x++;

  curstrp++;
  curstrp->strp = pinindexstrpp[3];
  curstrp->xyp[0] = 0.15;  curstrp->xyp[1] = +2.8;
  cstrvpinindex.x++;

  return &cstrvpinindex;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawSTfromEtoRZpsi(Xpsi, Zpsi, Vp, Vm, rhoE, zE)                        */
/*                                                                            */
/* Compute intersection of circle (rhoE,zE) with ST plane (R,Z)psi            */
/* return exist = 0 if no intersection                                        */
/*        exist = 1 if they intersect                                         */
/* store in coorST[] the 4 coordinates vm, vp, Zpsi, Xpsi                     */
/******************************************************************************/
short     gmDrawSTfromEtoRZpsi(double *coorST, double rhoE, double zE)
{ double    discr3, sqdiscr3;

  discr3 = rhoE*rhoE*discr1 - zE*zE*s2deltas2theta;
  if(discr3 < 0.0)
  { *coorST++ = 0.0;  *coorST++ = 0.0;  *coorST++ = 0.0;  *coorST = 0.0;
    return 0;
  }
  sqdiscr3 = sqrt(discr3);
  coorST[0] = zE/sqdiscr1;                                                 /** v- **/
  coorST[1] = sqdiscr3/sqdiscr1;                                           /** v+ **/
  coorST[2] = coorST[0] * sorientm - coorST[1] * corientm;         /** Zpsi = zST **/
  coorST[3] = coorST[0] * corientm + coorST[1] * sorientm;          /** Xpsi=Rpsi **/
  return 1;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawSTSetGeomRZ()                                                       */
/*                                                                            */
/******************************************************************************/
void      gmDrawSTSetGeomRZ()
{ 

  deltaDeg = gmDrawStateSTGetpsiRot() - phiDeg;
  deltaRad = deltaDeg * myDEG2RAD ;
  cdelta   = cos(deltaRad);
  sdelta   = sin(deltaRad);

  sdeltastheta   = stheta * sdelta;
  s2deltas2theta = sdeltastheta * sdeltastheta;
  cdeltastheta   = stheta * cdelta;
  discr1 = 1.0 - s2deltas2theta;
  sqdiscr1 = sqrt(discr1);
                              /** Orientation choice V- with respect to Xpsi axis **/
  corientm = cdeltastheta/sqdiscr1;
  sorientm = ctheta/sqdiscr1;

  rorientm = atan2(-cdeltastheta, ctheta);
  dorientm = rorientm * myRAD2DEG;
  rorientm = (myPI * 0.5) + rorientm;
  dorientm = 90.0 + dorientm;
  return;
}
/******************************************************************************/
/******************************************************************************/

