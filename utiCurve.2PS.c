/*  essai/C/utiCurve.2PS.c                                                    */
/*  Mennessier Gerard                   20010502                              */
/*  Last revised M.G.                   20030624                              */
/*                                                                            */

#include  <stddef.h>
#include  "utiMath.constant.def.h"

#include  "utistdIO.h"
#include  "utiCurve.2PS.h"
#include  "utiCurve.def.h"
#include  "utiCurve.paperFormat.h"

static    double    xscale = 1.0, yscale = 1.0;
static    int       currentpage = 0;
/******************************************************************************/
/*  PSgraphInit(FILE *psfile, double xnx[2], double ynx[2])                   */
/*                                                                            */
/*  WARNING % as special meaning inside format                                */
/*          ==> print %% as string, not inside format                         */
/******************************************************************************/
void      PSgraphInit(FILE *psfile)
{ static    char pshead1[] = "%!PS-Adobe-3.0\n", 
                 pshead2[] = "%%Creator esLevelCurve, author Mennessier 2001 May\n",
                 pshead3[] = "%%Pages: (atend)\n%%PageOrder: Ascend\n",
                 pshead4[] = "%%EndComments \n\n";
/*
  static    char psinitmat[] = "\n  initmatrix\n";
*/

                                /** Default unit is 1 pt = 1/72 inch = 25.4/72 mm **/
  fPrintF(psfile, "%s%s%s%s", pshead1,pshead2,pshead3,pshead4);

/*  fPrintF(psfile, "%s", psinitmat); */  /* Set init matrix which is DEVICE dependant */
  currentpage = 0;
  return;
}
/******************************************************************************/
/*  PSgraphPageInit(FILE *psfile)                                             */
/*                                                                            */
/******************************************************************************/
void      PSgraphPageInit(FILE *psfile)
{ static    char      psdeb1[] = "%%Page:";
  static    char      psdeb2[] = "\nsave\n";
  static    char psinitmat[] = "\n  initmatrix\n";

  currentpage++;
  fPrintF(psfile, "\n%s %d %d\n", psdeb1, currentpage, currentpage);
  fPrintF(psfile, "%s", psinitmat);
  fPrintF(psfile, psdeb2);

  return;
}
/******************************************************************************/
/*  PSgraphPageEnd(FILE *psfile)                                              */
/*                                                                            */
/******************************************************************************/
void      PSgraphPageEnd(FILE *psfile)
{ static    char      pstail[] = "\n  restore  showpage \n";

  fPrintF(psfile, "%s", pstail);
  return;
}
/******************************************************************************/
/*  PSgraphEnd(FILE *psfile)                                                  */
/*                                                                            */
/******************************************************************************/
void      PSgraphEnd(FILE *psfile)
{ static    char      pstrail1[] = "%%Trailer";
  static    char      pstrail2[] = "%%Pages:";
  static    char      pstail[]   = "%%EOF";

  fPrintF(psfile, "\n%s\n%s %d\n", pstrail1, pstrail2, currentpage);
  fPrintF(psfile, "\n%s\n", pstail);
  return;
}

/******************************************************************************/
/*  PSgraphSetInitCTM(FILE *psfile, double xnx[2], double ynx[2])             */
/*                                                                            */
/*  set the initial Current Transformation Matrix                             */
/*  which map  the rectangle xnx[2],ynx[2]  onto A4 ~ 210*290                 */
/******************************************************************************/
void      PSgraphSetInitCTM(FILE *psfile, double xnx[2], double ynx[2])
{ static    char formvarfdef[] =  "/%s  %f def \n";
  double    marge = 0.03, reff;
  double    xt, yt, linewidth = 1.;
  double    mm2pt = 72.0/25.4;                  /** 1 pt = 1/72 inch = 25.4/72 mm **/

  reff = 1.0 - 2. * marge;
  xscale = reff * A4W_ * mm2pt / (xnx[1] - xnx[0]);
  yscale = reff * A4H_ * mm2pt / (ynx[1] - ynx[0]);
  fPrintF(psfile, formvarfdef, "myxscale", xscale);
  fPrintF(psfile, formvarfdef, "myyscale", yscale);

  linewidth = 0.9 /xscale;                             /* will be affected by scale */
  fPrintF(psfile, formvarfdef, "mylinewidthunit", linewidth);
  xt = -xnx[0] + (xnx[1] - xnx[0]) * marge;
  yt = -ynx[0] + (ynx[1] - ynx[0]) * marge;
  fPrintF(psfile, "%f  %f  scale\n", xscale, yscale);
  fPrintF(psfile, "%f  setlinewidth\n", linewidth);
  fPrintF(psfile, "%f  %f  translate\n", xt, yt);
  return;
}
/******************************************************************************/
/*  PSgraphSetInitPortCTM(FILE *psfile, double xnx[2], double ynx[2])         */
/*                                                                            */
/*  set the initial Current Transformation Matrix                             */
/*  which map  the rectangle xnx[2],ynx[2]  onto A4 ~ 210*290                 */
/*  in PORTRAIT mode                                                          */
/******************************************************************************/
void      PSgraphSetInitPortCTM(FILE *psfile, double xnx[2], double ynx[2])
{ PSgraphSetInitCTM(psfile, xnx, ynx);  return;
}
/******************************************************************************/
/*  PSgraphSetInitLandCTM(FILE *psfile, double xnx[2], double ynx[2])         */
/*                                                                            */
/*  set the initial Current Transformation Matrix                             */
/*  which map  the rectangle xnx[2],ynx[2]  onto A4 ~ 210*290                 */
/*  in LANDSCAPE mode                                                         */
/******************************************************************************/
void      PSgraphSetInitLandCTM(FILE *psfile, double xnx[2], double ynx[2])
{ static    char formvarfdef[] =  "/%s  %f def \n";
  double    marge = 0.03, reff;
  double    xt, yt, linewidth = 1.;
  double    angledeg = 90.0;
  double    mm2pt = 72.0/25.4;                  /** 1 pt = 1/72 inch = 25.4/72 mm **/

  xt = (1.0 - marge) * A4W_ * mm2pt;  yt = marge * A4H_ * mm2pt;
  fPrintF(psfile, "%f  %f  translate\n", xt, yt);
  fPrintF(psfile, "%f  rotate\n", angledeg);

  reff = 1.0 - 2. * marge;
  xscale = reff * A4H_ * mm2pt / (xnx[1] - xnx[0]);
  yscale = reff * A4W_ * mm2pt / (ynx[1] - ynx[0]);
  fPrintF(psfile, formvarfdef, "myxscale", xscale);
  fPrintF(psfile, formvarfdef, "myyscale", yscale);

  linewidth = 0.9 /xscale;                             /* will be affected by scale */
  fPrintF(psfile, formvarfdef, "mylinewidthunit", linewidth);
  xt = -xnx[0];
  yt = -ynx[0];
  fPrintF(psfile, "%f  %f  scale\n", xscale, yscale);
  fPrintF(psfile, "%f  setlinewidth\n", linewidth);
  fPrintF(psfile, "%f  %f  translate\n", xt, yt);
  return;
}

/******************************************************************************/
/*  PSgraphSet(FILE *psfile, cSetVec *setvecp)                                */
/*                                                                            */
/******************************************************************************/
void      PSgraphSet(FILE *psfile, cSetVec *setvecp)
{ cSet     *setp;
  void     *p;
  int       i, type;
  static    char      prognamp[] = "utiCurve.2PS::PSgraphSet";

  for(i = 0, setp = setvecp->p;  i < setvecp->x;  i++, setp++)
  { type = setp->type;  p = setp->p;
    switch(type)
    { case MY_CPOLY:    { PSgraphPoly    (psfile, (cPoly*)p);    break;}
      case MY_CPOLYV:   { PSgraphPolyV   (psfile, (cPolyVec*)p); break;}
      case MY_CSEGV:    { PSgraphSegV    (psfile, (cSegVec*)p);  break;}
      case MY_LCSEGV:   { PSgraphLcSeg   (psfile, (lcSegVec*)p); break;}
      case MY_CSTRING:  { PSgraphStr     (psfile, (cStrVec*)p);  break;}
      case MY_CCIRCLE:  { PSgraphCircle  (psfile, (cCircle*)p, NULL);   break;}
      case MY_CCIRCLEV: { PSgraphCircleV (psfile, (cCircleVec*)p);      break;}
      case MY_CELLIPSE: { PSgraphEllipse (psfile, (cEllipse*)p, NULL);  break;}
      case MY_CELLIPSEV:{ PSgraphEllipseV(psfile, (cEllipseVec*)p);     break;}

      default: fPrintF(stderr,"%s  unknown CurveSet type = %d\n", prognamp,type);
    }
    fPrintF(psfile, "\n");
  }
  return;
}

/******************************************************************************/
/*  PSgraphLcSeg(FILE *psfile, lcSegVec *segvecp)                             */
/*                                                                            */
/******************************************************************************/
void      PSgraphLcSeg(FILE *psfile, lcSegVec *segvecp)
{ lcSeg    *segp;
  int       nsg;
  curveGC  *gcp;

  fPrintF(psfile, "gsave  ");
  fPrintF(psfile, "newpath\n");

  gcp = segvecp->gcp;
  if(gcp == NULL) fPrintF(psfile, "mylinewidthunit  setlinewidth\n");
  else            
  { fPrintF(psfile, "mylinewidthunit  %d  mul setlinewidth\n", gcp->lineWidth);
  }

  for(nsg = 0, segp = segvecp->p;  nsg < segvecp->x;  nsg++, segp++)
  { fPrintF(psfile, "%f  %f  moveto\n", segp->x[0], segp->y[0]);
    fPrintF(psfile, "%f  %f  lineto\n", segp->x[1], segp->y[1]);
  }
  fPrintF(psfile, "stroke  ");
  fPrintF(psfile, "grestore\n");

  return;
}
/******************************************************************************/
/*  PSgraphSegV(FILE *psfile, cSegVec *segvecp)                               */
/*                                                                            */
/******************************************************************************/
void      PSgraphSegV(FILE *psfile, cSegVec *segvecp)
{ cSeg     *segp;
  int       nsg;
  curveGC  *gcp;

  fPrintF(psfile, "gsave  ");
  fPrintF(psfile, "newpath\n");

  gcp = segvecp->gcp;
  if(gcp == NULL) fPrintF(psfile, "mylinewidthunit  setlinewidth\n");
  else            
  { fPrintF(psfile, "mylinewidthunit  %d  mul setlinewidth\n", gcp->lineWidth);
  }

  for(nsg = 0, segp = segvecp->p;  nsg < segvecp->x;  nsg++, segp++)
  { fPrintF(psfile, "%f  %f  moveto\n", segp->x[0], segp->y[0]);
    fPrintF(psfile, "%f  %f  lineto\n", segp->x[1], segp->y[1]);
  }
  fPrintF(psfile, "stroke  ");
  fPrintF(psfile, "grestore\n");

  return;
}
/******************************************************************************/
/*  PSgraphPoly(FILE *psfile, cPoly *p)                                       */
/*                                                                            */
/******************************************************************************/
void      PSgraphPoly(FILE *psfile, cPoly *p)
{ double   *xp, *yp, *xip, *yip;
  int       nsg;
  curveGC  *gcp;

  fPrintF(psfile, "gsave  ");
  fPrintF(psfile, "newpath\n");

  gcp = p->gcp;
  if(gcp == NULL) fPrintF(psfile, "mylinewidthunit  setlinewidth\n");
  else            
  { fPrintF(psfile, "mylinewidthunit  %d  mul setlinewidth\n", gcp->lineWidth);
  }

  xip = xp = p->xp;  yip = yp = p->yp;
  fPrintF(psfile, "%f  %f  moveto\n", *xp++, *yp++);

  for(nsg = 1; nsg < p->x;  nsg++)
  { fPrintF(psfile, "%f  %f  lineto\n", *xp++, *yp++);
  }
  if(p->closed) fPrintF(psfile, "%f  %f  lineto\n", *xip, *yip);
  fPrintF(psfile, "stroke  ");
  fPrintF(psfile, "grestore\n");

  return;
}
/******************************************************************************/
/*  PSgraphPolyV(FILE *psfile, cPolyVec *p)                                   */
/*                                                                            */
/******************************************************************************/
void      PSgraphPolyV(FILE *psfile, cPolyVec *vecp)
{ int       n;
  curveGC  *gcp;
  cPoly    *p;

  gcp = vecp->gcp;
  p = vecp->p;
  for(n = 0;  n < vecp->x;  n++, p++)
  { if(p->gcp == NULL) p->gcp = gcp;
    PSgraphPoly(psfile, p);
  }
  return;
}
/******************************************************************************/
/*  PSgraphStr(FILE *psfile, cStrVec *strvecp)                                */
/*                                                                            */
/******************************************************************************/
void      PSgraphStr(FILE *psfile, cStrVec *strvecp)
{ cStr     *strp;
/*
  curveGC  *gcp;
  double    fontxyscale = 15.;
*/
  int       fontindex;
  int       n;
  static    char      timesp[] = "Times-Roman",  helvep[] = "Helvetica";
  char     *fontnamespp[2] = {timesp, helvep},  *fontnamesp;

  fontindex = strvecp->findex;
  fontnamesp = fontnamespp[fontindex];

  fPrintF(psfile, "\n");
/*
  fPrintF(psfile, "gsave\n");
  fPrintF(psfile, "  initmatrix\n");
  fPrintF(psfile, "  %f  %f  scale\n", strvecp->fscalp[0], strvecp->fscalp[1]);
*/

  fPrintF(psfile, "  /myfonta  ");
  fPrintF(psfile, "  /%s findfont \n", fontnamesp);
  fPrintF(psfile, "  [%f  0.  0.  %f  0.  0.] makefont \n",
                                             strvecp->fscalp[0], strvecp->fscalp[1]);
  fPrintF(psfile, "  def \n");

/* 
  fPrintF(psfile, "  /myfonta  "
  fPrintF(psfile, "  /%s findfont \n", fontnamesp);
  fPrintF(psfile, "  %f  scalefont \n",  fontxyscale);
  fPrintF(psfile, "  def \n");
  fPrintF(psfile, "grestore\n");
*/


  fPrintF(psfile, "  myfonta  setfont \n");
  for(n = 0, strp = strvecp->p;  n < strvecp->x;  n++, strp++)
  { fPrintF(psfile, "%f  %f  moveto\n", strp->xyp[0], strp->xyp[1]);
    fPrintF(psfile, "  (%s) show\n", strp->strp);
  }
  return;
}
/******************************************************************************/
/*  PSgraphCircle(FILE *psfile, cCircle *p, curveGC *vgcp)                    */
/*                                                                            */
/*  PS primitive arc                                                          */
/*  always draw in trigonometric direction (anticlockwise)                    */
/*  (adds k* 2 PI  to final angle to have value larger than initial one)      */
/******************************************************************************/
void      PSgraphCircle(FILE *psfile, cCircle *p, curveGC *vgcp)
{ double    xc, yc, r, angl, angli, anglf, c, s;
  double    deg2rad = myPI / 180.;
  curveGC  *gcp;

  xc = p->xp[0];  yc = p->xp[1];
  r = p->r;
  angli = p->anglp[0];  anglf = p->anglp[1];                 /** angles in DEGREE **/
  angl = angli * deg2rad;
  c = cos(angl);  s = sin(angl);

  gcp = p->gcp;                                 /** priority : circle GC if !NULL **/
  if(gcp == NULL) gcp = vgcp;                              /** else circle Vec GC **/

  fPrintF(psfile, "gsave  ");
  fPrintF(psfile, "newpath\n");

  if(gcp == NULL) fPrintF(psfile, "mylinewidthunit  setlinewidth\n");
  else            
  { fPrintF(psfile, "mylinewidthunit  %d  mul setlinewidth\n", gcp->lineWidth);
  }

  fPrintF(psfile, "%f  %f  moveto\n", xc + r * c, yc + r * s);
  fPrintF(psfile, "%f %f %f %f %f arc\n",  xc, yc, r, angli, anglf);

  fPrintF(psfile, "stroke  ");
  fPrintF(psfile, "grestore\n");
  return;
}
/******************************************************************************/
/*   PSgraphCircleV(FILE *psfile, cCircleVec *vecp)                           */
/*                                                                            */
/******************************************************************************/
void      PSgraphCircleV(FILE *psfile, cCircleVec *vecp)
{ curveGC  *gcp;
  cCircle  *ccp;
  int       n;

  gcp = vecp->gcp;
  ccp = vecp->p;
  for(n = 0;  n < vecp->x;  n++, ccp++){ PSgraphCircle(psfile, ccp, gcp);}
  return;
}
/******************************************************************************/
/*  PSgraphEllipse(FILE *psfile, cEllipse *p, curveGC *vgcp)                  */
/*                                                                            */
/*  PS primitive arc                                                          */
/*  always draw in trigonometric direction (anticlockwise)                    */
/*  (adds k* 2 PI  to final angle to have value larger than initial one)      */
/*                                                                            */
/*  here                                                                      */
/*  exchange of initial and final angles (if final < initial)                 */
/*  ==> draw from smallest angle to largest one                               */
/******************************************************************************/
void      PSgraphEllipse(FILE *psfile, cEllipse *p, curveGC *vgcp)
{ double    xc, yc, a,b;
  double    dangli, danglf, orient;                       /** variables in DEGREE **/
  double    rangli, ranglf;                               /** variables in RADIAN **/
  double    deg2rad = myPI / 180.;
  curveGC  *gcp;

  gcp = p->gcp;                                /** priority : ellipse GC if !NULL **/
  if(gcp == NULL) gcp = vgcp;                             /** else ellipse Vec GC **/

  xc = p->xp[0];  yc = p->xp[1];
  a = p->axisp[0];  b = p->axisp[1];
  orient = p->orient;   
  dangli = p->anglDp[0];      danglf = p->anglDp[1]; 

if(danglf < dangli){ danglf = p->anglDp[0];  dangli = p->anglDp[1];}

  rangli = dangli * deg2rad;  ranglf = danglf * deg2rad;
  fPrintF(psfile, "gsave  ");
  fPrintF(psfile, "newpath\n");
  fPrintF(psfile, "%f  %f  translate\n", xc, yc);
  fPrintF(psfile, "%f rotate\n", orient);
  fPrintF(psfile, "%f  %f  moveto\n", a * cos(rangli), b * sin(rangli));

  if(gcp == NULL) fPrintF(psfile, "mylinewidthunit  setlinewidth\n");
  else            
  { fPrintF(psfile, "mylinewidthunit  %d  mul setlinewidth\n", gcp->lineWidth);
  }
  fPrintF(psfile, "1.0   %f  scale\n", b/a);

  fPrintF(psfile, "%f %f %f %f %f arc\n",  0.0, 0.0, a, dangli, danglf);

  fPrintF(psfile, "stroke  ");
  fPrintF(psfile, "grestore\n");
/*
fPrintF(stderr,
  "\n\n utiCurve.2PS::PSgraphEllipse orient=%f,dangli=%f,danglf=%f,a=%f,b=%f,linewidth=%d \n\n",
                                                  orient, dangli, danglf, a, b, gcp->lineWidth);
*/
  return;
}
/******************************************************************************/
/*   PSgraphEllipseV(FILE *psfile, cEllipseVec *vecp)                         */
/*                                                                            */
/******************************************************************************/
void      PSgraphEllipseV(FILE *psfile, cEllipseVec *vecp)
{ curveGC  *gcp;
  cEllipse *ccp;
  int       n;

  gcp = vecp->gcp;
  ccp = vecp->p;
  for(n = 0;  n < vecp->x;  n++, ccp++){ PSgraphEllipse(psfile, ccp, gcp);}
  return;
}
/******************************************************************************/
/******************************************************************************/
