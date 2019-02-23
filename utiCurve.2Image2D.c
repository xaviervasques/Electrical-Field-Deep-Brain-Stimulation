/*  ../libmy/utiCurve.2Image2D.c                                              */
/*  Mennessier Gerard                   20030612                              */
/*       revised M.G.                   20041121                              */
/*  Last revised M.G.                   20050208                              */
/*                                                                            */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utistdIO.h"
#include  "utiMath.constant.def.h"
#include  "utiMath.type.def.h"
#include  "utiMath.macro.def.h"

#include  "utiImage2D.h"
#include  "utiImage2D.Phys.h"

#include  "utiCurve.def.h"
#include  "utiCurve.set.h"
#include  "utiCurve.poly.h"
#include  "utiCurve.seg.h"
#include  "utiCurve.circle.h"
#include  "utiCurve.ellipse.h"
#include  "utiCurve.level.h"
#include  "utiCurve.2Image2D.h"

#ifndef   myPIXELxySET
#define   myPIXELxySET(IX, IY)                             \
          cp = ecp + (IY)*widthbz + (IX)*ibx;              \
          inkp = inkip;                                    \
          for(ib = 0;  ib < ibx;  ib++) *cp++ = *inkp++;
#endif

static char    srcfilenamp[] = "utiCurve.2Image2D";
static short   boolPrintDebug = 0;
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      c2Image2DgraphSetInk(unsigned char *inkp)
{
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphSet(utiImage2D *ep, unsigned char *inkip,              */
/*                                        utiImage2DPhys *physp, cSeg *segp)  */
/*                                                                            */
/******************************************************************************/
void      c2Image2DgraphSet(utiImage2D *ep, unsigned char *inkip, 
                                             utiImage2DPhys *physp, cSetVec *setvecp)
{ cSet     *setp;
  void     *p;
  int       i, type;
  static  char   form3p[] = "%s::%s  string or unknown CurveSet type = %d\n";
  static  char   prognamp[] = "c2Image2DgraphSet";

  for(i = 0, setp = setvecp->p;  i < setvecp->x;  i++, setp++)
  { type = setp->type;  p = setp->p;
    switch(type)
    { case MY_CPOLY:    { c2Image2DgraphPoly    (ep, inkip, physp, (cPoly*)p);    break;}
      case MY_CPOLYV:   { c2Image2DgraphPolyV   (ep, inkip, physp, (cPolyVec*)p); break;}
      case MY_CSEGV:    { c2Image2DgraphSegV    (ep, inkip, physp, (cSegVec*)p);  break;}
      case MY_LCSEGV:   { c2Image2DgraphLcSeg   (ep, inkip, physp, (lcSegVec*)p); break;}
      case MY_CCIRCLE:  { c2Image2DgraphCircle  (ep, inkip, physp, (cCircle*)p);  break;}
      case MY_CCIRCLEV: { c2Image2DgraphCircleV (ep, inkip, physp, (cCircleVec*)p); break;}
      case MY_CELLIPSE: { c2Image2DgraphEllipse (ep, inkip, physp, (cEllipse*)p);   break;}
      case MY_CELLIPSEV:{ c2Image2DgraphEllipseV(ep, inkip, physp, (cEllipseVec*)p);break;}
      case MY_CPOINTS:  { c2Image2DgraphPoints  (ep, inkip, physp, (cPoints*)p);  break;}
      default:          { fPrintF(stderr, form3p, srcfilenamp,prognamp, type); break;}
    }
  }
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphSeg(utiImage2D *ep, unsigned char *inkip,              */
/*                                        utiImage2DPhys *physp, cSeg *segp)  */
/*                                                                            */
/* set all pixels crossed by the segment to value (color or gray)             */
/*     defined by the ibx (ep->bytesPerPixel) pointed by inkip                */
/******************************************************************************/
void      c2Image2DgraphSeg(utiImage2D *ep, unsigned char *inkip, 
                                                   utiImage2DPhys *physp, cSeg *segp)
{ int       j0, j1, j0i, j0f, j1i, j1f;
  double    dban, slope;
  double    deltax, deltay, dpx, dpy;
  myBOOL    exch01 = 0;
  double    effsegx0p[2], effsegx1p[2];
  double   *effwinx0p, *effwinx1p, effpixdx0, effpixdx1, effscale0, effscale1;
  int       j1ci, j1cf, j0ci, j0cf, iban;
  int       j1min, j1max;
  double    x0ci, x1ci, x0cf, x1cf;
  size_t    widthbz;
  int       x0xex, x1xex;
  int       ib, ibx;
  unsigned char *cp, *ecp;
  unsigned char *inkp;
  short     boolPrint = 0;
  static char    form1p[]   = "%s::%s  BEGIN\n";
  static char    form2p[]   = "%s::%s  END\n";
  static char    form3p[]   = "  %s  exch01=%d, j0i=%d, j0f=%d, j1i=%d, j1f=%d \n";
  static char    form4p[]   = "  %s %s  j0=%d, j1min=%d, j1max=%d \n";

  static char    form7p[]   = "  %s comment=%s j0i=%d, j0f=%d, j1i=%d, j1f=%d \n";
  static char    form8p[]   = "  %s GENERIC DEBUT \n";
  static char    form9p[]   = "  %s x0xex=%d, x1xex=%d, boolPrintDebug=%d,  boolPrint=%d\n";
  static char    prognamp[] = "c2Image2DgraphSeg";

if(boolPrintDebug) fPrintF(stderr, form1p, srcfilenamp, prognamp);
  ecp = (unsigned char *)(ep->dataV.p);
  ibx = ep->im0D.bytesPerPixel;
  widthbz = ep->wbz;

  deltax = segp->x[1] - segp->x[0];  deltay = segp->y[1] - segp->y[0];
  dpx = deltax * physp->scalex;      dpy = deltay * physp->scaley;

                                /** choose effx0 to be the lowest (in pixel unit) **/
  if(fabs(dpx) <= fabs(dpy))
  { exch01 = 0;
    effwinx0p = physp->xnxp;    effwinx1p = physp->ynxp;  
    effpixdx0 = physp->pixdx;   effpixdx1 = physp->pixdy;
    effscale0 = physp->scalex;  effscale1 = physp->scaley;
    effsegx0p[0] = segp->x[0];  effsegx0p[1] = segp->x[1];
    effsegx1p[0] = segp->y[0];  effsegx1p[1] = segp->y[1];
    x0xex = ep->wpx;  x1xex = ep->hpx;
  }
  else
  { exch01 = 1;
    effwinx0p = physp->ynxp;    effwinx1p = physp->xnxp;  
    effpixdx0 = physp->pixdy;   effpixdx1 = physp->pixdx;
    effscale0 = physp->scaley;  effscale1 = physp->scalex;
    effsegx0p[0] = segp->y[0];  effsegx0p[1] = segp->y[1];
    effsegx1p[0] = segp->x[0];  effsegx1p[1] = segp->x[1];
    x0xex = ep->hpx;  x1xex = ep->wpx;
  }
                                         /** exchange 0 and 1 segment extremities **/
                                   /** so that 0 index is increasing when drawing **/
  dban = (effsegx0p[0] - effwinx0p[0]) * effscale0;    j0i = myFLOORi(dban);
  dban = (effsegx0p[1] - effwinx0p[0]) * effscale0;    j0f = myFLOORi(dban);
  if(j0f < j0i)
  { mySWAP(effsegx0p[0], effsegx0p[1], dban);
    mySWAP(effsegx1p[0], effsegx1p[1], dban);
    mySWAP(j0i, j0f, iban);
  }

  dban = (effsegx1p[0] - effwinx1p[0]) * effscale1;    j1i = myFLOORi(dban);
  dban = (effsegx1p[1] - effwinx1p[0]) * effscale1;    j1f = myFLOORi(dban);
if(boolPrintDebug) fPrintF(stderr, form3p, prognamp, exch01, j0i, j0f, j1i, j1f);
if(j0i < 0      || j0f < 0      || j1i < 0      || j1f < 0)      boolPrint = 1;
if(j0i >= x0xex || j0f >= x0xex || j1i >= x1xex || j1f >= x1xex) boolPrint = 1;
if(boolPrintDebug) fPrintF(stderr, form9p, prognamp, x0xex, x1xex, boolPrintDebug, boolPrint);

  if(j0f < 0)
  { if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp,prognamp);  return;
  }
  if(j0i >=  x0xex)
  { if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp,prognamp);  return;
  }

  if(myMAX(j1i, j1f) < 0)
  { if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp,prognamp);  return;
  }
  if(myMIN(j1i, j1f) >= x1xex)
  { if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp,prognamp);  return;
  }

  slope = (effsegx1p[1] - effsegx1p[0]) / (effsegx0p[1] - effsegx0p[0]);

                                                                  /** FIRST pixel **/
  if(j0i == j0f)
  { x0cf = effsegx0p[1];  
    x1cf = effsegx1p[1];
  }
  else
  { x0cf = effwinx0p[0] + (j0i +1) * effpixdx0;
    x1cf = effsegx1p[0] + (x0cf - effsegx0p[0]) * slope;
  }
  j1ci = j1i;
  dban = (x1cf - effwinx1p[0]) * effscale1;
  j1cf = myFLOORi(dban);
                                                           /** check inside image **/
  if(j0i >= x0xex)
  { if(boolPrintDebug) fPrintF(stderr, form7p, prognamp, "j0i >= x0xex", j0i, j0f, j1i, j1f);
    if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp,prognamp);  return;
  }
  if(j0i >= 0)
  {
    j1min = myMIN(j1ci, j1cf);  j1max = myMAX(j1ci, j1cf);
if(boolPrintDebug) fPrintF(stderr, form4p, prognamp, "FIRST", j0i, j1min, j1max);
                                                        /** limiting inside image **/
    if(j1min < 0) j1min = 0;  if(j1max >= x1xex) j1max = x1xex -1;
    if(exch01)
    { for(j1 = j1min;  j1 <= j1max;  j1++){ myPIXELxySET(j1, j0i)}
    }
    else
    { for(j1 = j1min;  j1 <= j1max;  j1++){ myPIXELxySET(j0i, j1)}
    }
  }
  if(j0i == j0f)
  { if(boolPrintDebug) fPrintF(stderr, form7p, prognamp, "j0i == j0f", j0i, j0f, j1i, j1f);
    if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp,prognamp);  return;
  }
                                                        /** GENERIC along segment **/
if(boolPrint)      fPrintF(stderr, form8p, prognamp);
if(boolPrintDebug) fPrintF(stderr, form8p, prognamp);
  j0cf = j0i + 1;
  if(j0cf < 0) j0cf = 0;
  x0cf = effwinx0p[0] + j0cf * effpixdx0;
  x1cf = effsegx1p[0] + (x0cf - effsegx0p[0]) * slope;
  dban = (x1cf - effwinx1p[0]) * effscale1;
  j1cf = myFLOORi(dban);
  j0ci = j0cf;
                                                           /** check inside image **/
  for(j0 = j0ci;  j0 < myMIN(j0f, x0xex);  j0++)
  { j1ci = j1cf;
    x0cf = effwinx0p[0] + (j0 +1) * effpixdx0;
    x1cf = effsegx1p[0] + (x0cf - effsegx0p[0]) * slope;
    dban = (x1cf - effwinx1p[0]) * effscale1;
    j1cf = myFLOORi(dban);

    j1min = myMIN(j1ci, j1cf);  j1max = myMAX(j1ci, j1cf);
if(boolPrint) fPrintF(stderr, form4p, prognamp, "GENERIC", j0, j1min, j1max);
                                                        /** limiting inside image **/
    if(j1min < 0) j1min = 0;  if(j1max >= x1xex) j1max = x1xex -1;
    if(exch01)
    { for(j1 = j1min;  j1 <= j1max;  j1++){ myPIXELxySET(j1, j0)}
    }
    else
    { for(j1 = j1min;  j1 <= j1max;  j1++){ myPIXELxySET(j0, j1)}
    }
  }

                                                                   /** LAST pixel **/
  x0cf = effsegx0p[1];
  x1cf = effsegx1p[1];
  j0ci = j0f;
  j1cf = j1f;
  x0ci = effwinx0p[0] + j0ci * effpixdx0;
  x1ci = effsegx1p[0] + (x0ci - effsegx0p[0]) * slope;
  dban = (x1ci - effwinx1p[0]) * effscale1;
  j1ci = myFLOORi(dban);
                                                           /** check inside image **/
  if(j0f < x0xex)
  {
    j1min = myMIN(j1ci, j1cf);  j1max = myMAX(j1ci, j1cf);
if(boolPrintDebug)fPrintF(stderr, form4p, prognamp, "LAST", j0f, j1min, j1max);
                                                        /** limiting inside image **/
    if(j1min < 0) j1min = 0;  if(j1max >= x1xex) j1max = x1xex -1;

    if(exch01)
    { for(j1 = j1min;  j1 <= j1max;  j1++){ myPIXELxySET(j1, j0f)}
    }
    else
    { for(j1 = j1min;  j1 <= j1max;  j1++){ myPIXELxySET(j0f, j1)}
    }
  }
if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp,prognamp);
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphSegV(utiImage2D *ep, unsigned char *inkip,             */
/*                                  utiImage2DPhys *physp, cSegVec *segvecp)  */
/*                                                                            */
/* set all pixels crossed by the segment to value (color or gray)             */
/*     defined by the ibx (ep->bytesPerPixel) pointed by inkip                */
/******************************************************************************/
void      c2Image2DgraphSegV(utiImage2D *ep, unsigned char *inkip, 
                                             utiImage2DPhys *physp, cSegVec *segvecp)
{ int       nsg;
  cSeg     *segp;

  segp = segvecp->p;
  for(nsg = 0;  nsg < segvecp->x;  nsg++, segp++)
  { c2Image2DgraphSeg(ep, inkip, physp, segp);
  }
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphPoly(utiImage2D *ep, unsigned char *inkip,             */
/*                                          utiImage2DPhys *physp, cPoly *p)  */
/*                                                                            */
/* set all pixels crossed by the segment to value *inkip (color or gray)      */
/******************************************************************************/
void      c2Image2DgraphPoly(utiImage2D *ep, unsigned char *inkip,
                                                     utiImage2DPhys *physp, cPoly *p)
{ int       i, ix;
  cSeg      seg = {{0.0, 0.0}, {0.0, 0.0}};
  double   *xp, *yp;
  static char    form1p[]   = "%s::%s  BEGIN\n";
  static char    form2p[]   = "%s::%s  END\n";
  static char    prognamp[] = "c2Image2DgraphPoly";

if(boolPrintDebug) fPrintF(stderr, form1p, srcfilenamp, prognamp);
  ix = p->x;
  xp = p->xp;
  yp = p->yp;

  seg.x[1] = *xp++;  seg.y[1] = *yp++;
  for(i = 1;  i < ix;  i++)
  { seg.x[0 ] = seg.x[1];  seg.y[0 ] = seg.y[1];
    seg.x[1] = *xp++;  seg.y[1] = *yp++;
    c2Image2DgraphSeg(ep, inkip, physp, &seg);
  }
  if(p->closed)
  { seg.x[0 ] = seg.x[1];  seg.y[0 ] = seg.y[1];
    seg.x[1] = *(p->xp);  seg.y[1] = *(p->yp);
    c2Image2DgraphSeg(ep, inkip, physp, &seg);
  }
if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp, prognamp);
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphCircle(utiImage2D *ep, unsigned char *inkip,           */
/*                                     utiImage2DPhys *physp, cCircle *ccp)   */
/*                                                                            */
/* set all pixels crossed by the arc of circle to *inkip  (color or gray)     */
/******************************************************************************/
void      c2Image2DgraphCircle(utiImage2D *ep, unsigned char *inkip, 
                                                 utiImage2DPhys *physp, cCircle *ccp)
{ int       i, ix;
  double    xc, yc, r;
  double    dangli, danglf;                               /** variables in DEGREE **/
  double    rstep, rangli, ranglf, rangl, rstepif;        /** variables in RADIAN **/
  double    pixmin;
  myBOOL    sign;
  cSeg      seg = {{0.0, 0.0}, {0.0, 0.0}};
  static char    form1p[]   = "%s::%s  BEGIN\n";
  static char    form2p[]   = "%s::%s  END\n";
  static char    form3p[]   = "  center=%f, %f;  r=%f; anglesD=%f, %f; anglesR=%f, %f\n";
  static char    form4p[]   = "  ix=%d, rstep=%f, rstepif=%f\n";
  static char    form5p[]   = "  %s  i=%d, rangl=%f\n";
  static char    prognamp[] = "c2Image2DgraphCircle";

if(boolPrintDebug) fPrintF(stderr, form1p, srcfilenamp, prognamp);
  xc = ccp->xp[0];  yc = ccp->xp[1];
  r = fabs(ccp->r);
  dangli = ccp->anglp[0];  danglf = ccp->anglp[1];
  rangli = dangli * myDEG2RAD;  ranglf = danglf * myDEG2RAD;
if(boolPrintDebug) fPrintF(stderr, form3p, xc, yc, r, dangli, danglf, rangli, ranglf);
  rangl = ranglf - rangli;
  sign = (rangl>= 0.0)? 1: 0;
  rangl = fabs(rangl);
          /** angular step choosen as a function of pixel_size compared to radius **/
  pixmin = myMIN(fabs(physp->pixdx), fabs(physp->pixdy));
  if(r < pixmin){ rstep = myPI * 0.25;}
  else          { rstep = pixmin/r * 0.5;}
  ix = myMAX(4, (int)(rangl/rstep));
  rstep = rangl/(ix +1);
  rstepif = (rangl - ix * rstep + rstep) * 0.5;
  if(!sign){ rstepif = - rstepif;  rstep = - rstep;}
if(boolPrintDebug) fPrintF(stderr, form4p, ix, rstep, rstepif);

  rangl = rangli;
  seg.x[1] = xc + r * cos(rangl);  seg.y[1] = yc + r * sin(rangl);  

  for(i = 0;  i < ix;  i++)
  { seg.x[0] = seg.x[1];  seg.y[0] = seg.y[1];
    rangl = (rangli + rstepif + i * rstep);
if(boolPrintDebug) fPrintF(stderr, form5p, prognamp, i, rangl);
    seg.x[1] = xc + r * cos(rangl);  seg.y[1] = yc + r * sin(rangl); 
    c2Image2DgraphSeg(ep, inkip, physp, &seg);
  }
  
  seg.x[0] = seg.x[1];  seg.y[0] = seg.y[1];
  rangl = ranglf;
if(boolPrintDebug) fPrintF(stderr, form5p, prognamp, ix, rangl);
  seg.x[1] = xc + r * cos(rangl);  seg.y[1] = yc + r * sin(rangl); 
  c2Image2DgraphSeg(ep, inkip, physp, &seg);

if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp, prognamp);
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphEllipse(utiImage2D *ep, unsigned char *inkip,          */
/*                                       utiImage2DPhys *physp, cEllipse *p)  */
/*                                                                            */
/* set all pixels crossed by the arc of circle to *inkip  (color or gray)     */
/******************************************************************************/
void      c2Image2DgraphEllipse(utiImage2D *ep, unsigned char *inkip, 
                                                  utiImage2DPhys *physp, cEllipse *p)
{ int       i, ix;
  double    xc, yc, a, b, rx, ry, x, y;
  double    dangli, danglf;                               /** variables in DEGREE **/
  double    dorient, rorient, co, so;
  double    rstep, rangli, ranglf, rangl, rstepif;        /** variables in RADIAN **/
  double    pixx, pixy;
  myBOOL    sign;
  cSeg      seg = {{0.0, 0.0}, {0.0, 0.0}};
  static char    form1p[]   = "%s::%s  BEGIN\n";
  static char    form2p[]   = "%s::%s  END\n";
  static char    form3p[]   = "  center=%f, %f;  a=%f, b=%f; anglesDir=(%f, %f)Deg=(%f, %f)Rad; orient=%f Deg=%f Rad\n";
  static char    form4p[]   = "  ix=%d, rstep=%f, rstepif=%f\n";
  static char    form5p[]   = "  %s  i=%d, rangl=%f \n";
  static char    prognamp[] = "c2Image2DgraphEllipse";

if(boolPrintDebug) fPrintF(stderr, form1p, srcfilenamp, prognamp);
  xc = p->xp[0];    yc = p->xp[1];
  a = fabs(p->axisp[0]);  b = fabs(p->axisp[1]); 
  dorient = p->orient;
  rorient = dorient * myDEG2RAD;  co = cos(rorient);  so = sin(rorient);
  dangli = p->anglDp[0];  danglf = p->anglDp[1];              /** director angles **/
  rangli = dangli * myDEG2RAD;  ranglf = danglf * myDEG2RAD;
if(boolPrintDebug) fPrintF(stderr, form3p, xc, yc, a,b, dangli,danglf, rangli,ranglf, dorient,rorient);
  rangl = ranglf - rangli;
  sign = (rangl>= 0.0)? 1: 0;
  rangl = fabs(rangl);
  rx = a * fabs(co) + b * fabs(so);  ry = a * fabs(so) + b * fabs(co);
          /** angular step choosen as a function of pixel_size compared to radius **/
  pixx = fabs(physp->pixdx);  pixy = fabs(physp->pixdy);
  if(rx < pixx){ rstep = myPI * 0.5;}
  else         { rstep = pixx/rx;}
  if(ry < pixy){ rstep = myMIN(myPI * 0.5, rstep);}
  else         { rstep = myMIN(pixy/ry,    rstep);}

  ix = myMAX(4, (int)(rangl/rstep));
  rstep = rangl/(ix +1);
  rstepif = (rangl - ix * rstep + rstep) * 0.5;
  if(!sign){ rstepif = - rstepif;  rstep = - rstep;}
if(boolPrintDebug) fPrintF(stderr, form4p, ix, rstep, rstepif);

  rangl = rangli;
  x = a * cos(rangl);  y = b * sin(rangl);
  seg.x[1] = xc + (co * x - so * y);  seg.y[1] = yc + (so * x + co * y);
  for(i = 0;  i < ix;  i++)
  { seg.x[0] = seg.x[1];  seg.y[0] = seg.y[1];
    rangl = (rangli + rstepif + i * rstep);
if(boolPrintDebug) fPrintF(stderr, form5p, prognamp, i, rangl);
    x = a * cos(rangl);  y = b * sin(rangl);
    seg.x[1] = xc + (co * x - so * y);  seg.y[1] = yc + (so * x + co * y);
    c2Image2DgraphSeg(ep, inkip, physp, &seg);
  }

  seg.x[0] = seg.x[1];  seg.y[0] = seg.y[1];
  rangl = ranglf;
if(boolPrintDebug) fPrintF(stderr, form5p, prognamp, ix, rangl);
  x = a * cos(rangl);  y = b * sin(rangl);
  seg.x[1] = xc + (co * x - so * y);  seg.y[1] = yc + (so * x + co * y);
  c2Image2DgraphSeg(ep, inkip, physp, &seg);

if(boolPrintDebug) fPrintF(stderr, form2p, srcfilenamp, prognamp);
  return;
}

/******************************************************************************/
/*  void c2Image2DgraphLcSeg(utiImage2D *ep, unsigned char *inkip,            */
/*                                   utiImage2DPhys *physp, lcSegVec *lcsvp)  */
/*                                                                            */
/* set all pixels crossed by the segment to value (color or gray)             */
/*     defined by the ibx (ep->bytesPerPixel) pointed by inkip                */
/******************************************************************************/
void      c2Image2DgraphLcSeg(utiImage2D *ep, unsigned char *inkip, 
                                              utiImage2DPhys *physp, lcSegVec *lcsvp)
{ int       nsg;
  lcSeg    *lcsp;
  cSeg      seg;

  lcsp = lcsvp->p;
  for(nsg = 0;  nsg < lcsvp->x;  nsg++, lcsp++)
  { seg.x[0] = lcsp->x[0];  seg.x[1] = lcsp->x[1];
    seg.y[0] = lcsp->y[0];  seg.y[1] = lcsp->y[1];
    c2Image2DgraphSeg(ep, inkip, physp, &seg);
  }
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphPolyV(utiImage2D *ep, unsigned char *inkip,            */
/*                                    utiImage2DPhys *physp, cPolyVec *ccvp)  */
/*                                                                            */
/******************************************************************************/
void      c2Image2DgraphPolyV(utiImage2D *ep, unsigned char *inkip, 
                                             utiImage2DPhys *physp, cPolyVec *vecp)
{ int       n;
  cPoly    *ccp;

  ccp = vecp->p;
  for(n = 0;  n < vecp->x;  n++, ccp++) c2Image2DgraphPoly(ep, inkip, physp, ccp);
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphCircleV(utiImage2D *ep, unsigned char *inkip,          */
/*                                  utiImage2DPhys *physp, cCircleVec *ccvp)  */
/*                                                                            */
/******************************************************************************/
void      c2Image2DgraphCircleV(utiImage2D *ep, unsigned char *inkip, 
                                             utiImage2DPhys *physp, cCircleVec *vecp)
{ int       n;
  cCircle  *ccp;

  ccp = vecp->p;
  for(n = 0;  n < vecp->x;  n++, ccp++) c2Image2DgraphCircle(ep, inkip, physp, ccp);
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphEllipseV(utiImage2D *ep, unsigned char *inkip,         */
/*                                  utiImage2DPhys *physp, cEllipseVec *ccvp) */
/*                                                                            */
/******************************************************************************/
void      c2Image2DgraphEllipseV(utiImage2D *ep, unsigned char *inkip, 
                                            utiImage2DPhys *physp, cEllipseVec *vecp)
{ int       n;
  cEllipse *ccp;

  ccp = vecp->p;
  for(n = 0;  n < vecp->x;  n++, ccp++) c2Image2DgraphEllipse(ep, inkip, physp, ccp);
  return;
}
/******************************************************************************/
/*  void c2Image2DgraphPoints(utiImage2D *ep, unsigned char *inkip,           */
/*                                        utiImage2DPhys *physp, cPoints *cp) */
/*                                                                            */
/* VERIFIER LES ARRONDIS SELON SIGNE DE scaleX(Y) et (y)xnxp[1]-(y)xnxp[0]    */

/******************************************************************************/
void      c2Image2DgraphPoints(utiImage2D *ep, unsigned char *inkip, 
                                                   utiImage2DPhys *physp, cPoints *p)
{ unsigned char *ecp, *cp;
  int       ib, ibx;
  size_t    widthbz;
  double    xn, yn;                  /** xn = left extremity; yn = y of first row **/
  int       i;
  int       j0, j1;
  double   *dp, x, y, dban;
  unsigned char *inkp;

  ecp = (unsigned char *)(ep->dataV.p);
  ibx = ep->im0D.bytesPerPixel;
  widthbz = ep->wbz;

  xn = physp->xnxp[0];  yn = physp->ynxp[0];

  dp = p->xp;
  for(i = 0;  i < p->x;  i++)
  { x = *dp++;  y = *dp++;
    dban = (x - xn) * physp->scalex;
    j0 = myFLOORi(dban);
    dban = (y - yn) * physp->scaley;
    j1 = myFLOORi(dban);
    myPIXELxySET(j0, j1)
  }
  return;
}
/******************************************************************************/
#undef  myPIXELxySET

#include  "utiMath.macro.undef.h"
/******************************************************************************/
/******************************************************************************/
