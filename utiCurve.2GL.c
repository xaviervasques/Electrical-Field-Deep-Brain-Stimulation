/*  essai/C/utiCurve.2GL.c                                                    */
/*  Mennessier Gerard                   20010629                              */
/*  Last revised M.G.                   20030626                              */
/*                                                                            */

#include  <stddef.h>
#include  <math.h>
#include  <GL/glut.h>
#include  "utiMath.constant.def.h"

#include  "utistdIO.h"
#include  "utiCurve.2GL.h"
#include  "utiCurve.def.h"

/******************************************************************************/
/*  GLgraphInit(int win, double xnx[2], double ynx[2])                        */
/*                                                                            */
/******************************************************************************/
void      GLgraphInit(int win, double xnx[2], double ynx[2])
{ 
  double    xscale = 1.0, yscale = 1.0, xt, yt, linewidth = 1.;
  double    mm2pt = 72.0/25.4;                  /** 1 pt = 1/72 inch = 25.4/72 mm **/

  xscale = 0.95 * 200. * mm2pt / (xnx[1] - xnx[0]);
  yscale = 0.95 * 250. * mm2pt / (ynx[1] - ynx[0]);
  linewidth = .1 /xscale;                             /* will be affected by scale */
  xt = -xnx[0] + (xnx[1] - xnx[0]) * 0.05;
  yt = -ynx[0] + (ynx[1] - ynx[0]) * 0.05;

  return;
}
/******************************************************************************/
/*  GLgraphEnd(int win)                                                       */
/*                                                                            */
/******************************************************************************/
void      GLgraphEnd(int win)
{ 
  return;
}

/******************************************************************************/
/*  GLgraphSet(int win, cSetVec *setvecp)                                     */
/*                                                                            */
/******************************************************************************/
void      GLgraphSet(int win, cSetVec *setvecp)
{ cSet     *setp;
  void     *p;
  int       i, type;
  static    char      prognamp[] = "utiCurve.2GL::GLgraphSet";

  for(i = 0, setp = setvecp->p;  i < setvecp->x;  i++, setp++)
  { type = setp->type;  p = setp->p;
    switch(type)
    { case MY_CPOLY:    { GLgraphPoly    (win, (cPoly*)p);    break;}
      case MY_CPOLYV:   { GLgraphPolyV   (win, (cPolyVec*)p); break;}
      case MY_CSEGV:    { GLgraphSegV    (win, (cSegVec*)p);  break;}
      case MY_LCSEGV:   { GLgraphLcSeg   (win, (lcSegVec*)p); break;}
      case MY_CSTRING:  { GLgraphStr     (win, (cStrVec*)p);  break;}
      case MY_CCIRCLE:  { GLgraphCircle  (win, (cCircle*)p, NULL);  break;}
      case MY_CCIRCLEV: { GLgraphCircleV (win, (cCircleVec*)p);     break;}
      case MY_CELLIPSE: { GLgraphEllipse (win, (cEllipse*)p, NULL);  break;}
      case MY_CELLIPSEV:{ GLgraphEllipseV(win, (cEllipseVec*)p);  break;}
      default: fPrintF(stderr,"%s  unknown CurveSet type = %d\n", prognamp,type);
    }
  }
  return;
}
/******************************************************************************/
/*  GLgraphSetGC(curveGC  *gcp)                                               */
/*                                                                            */
/******************************************************************************/
void      GLgraphSetGC(curveGC  *gcp)
{ 

  if(gcp == NULL)
  { glLineWidth(1.0);
    glColor3f(0.0, 0.0, 0.0);                                /** Default to Black **/
  }
  else
  { glLineWidth(gcp->lineWidth);
    glColor3f(gcp->color[0], gcp->color[1], gcp->color[2]);
  }

/*
  if(gcp != NULL) glLineWidth(gcp->lineWidth);
  else            glLineWidth(1.0);
*/
  return;
}
/******************************************************************************/
/*  GLgraphPoly(int win, cPoly *p)                                            */
/*                                                                            */
/*  GLfloat   lw;                                                             */
/*  lw = gcp->lineWidth;                                                      */
/*  glLineWidth(lw);                                                          */
/******************************************************************************/
void      GLgraphPoly(int win, cPoly *p)
{ double   *xp, *yp, *xip, *yip;
  int       nsg;

  GLgraphSetGC(p->gcp);
  if(p->closed) glBegin(GL_LINE_LOOP);
  else          glBegin(GL_LINE_STRIP);

  xip = xp = p->xp;  yip = yp = p->yp;
  for(nsg = 0;  nsg < p->x;  nsg++)
  { glVertex2f(*xp++, *yp++);
  }
  glEnd();
  return;
}
/******************************************************************************/
/*  GLgraphPolyV(int win, cPolyVec *p)                                        */
/*                                                                            */
/******************************************************************************/
void      GLgraphPolyV(int win, cPolyVec *vecp)
{ int       n;
  curveGC  *gcp;
  cPoly    *p;

  gcp = vecp->gcp;
  p = vecp->p;
  for(n = 0;  n < vecp->x;  n++, p++)
  { if(p->gcp == NULL) p->gcp = gcp;
    GLgraphPoly(win, p);
  }
  return;
}
/******************************************************************************/
/*  GLgraphSeg(int win, cSegVec *segvecp)                                     */
/*                                                                            */
/******************************************************************************/
void      GLgraphSegV(int win, cSegVec *segvecp)
{ cSeg     *segp;
  int       nsg;

  GLgraphSetGC(segvecp->gcp);
  glBegin(GL_LINES);
  for(nsg = 0, segp = segvecp->p;  nsg < segvecp->x;  nsg++, segp++)
  { glVertex2f(segp->x[0], segp->y[0]);
    glVertex2f(segp->x[1], segp->y[1]);
  }
  glEnd();
  return;
}
/******************************************************************************/
/*  GLgraphLcSeg(int win, lcSegVec *segvecp)                                  */
/*                                                                            */
/******************************************************************************/
void      GLgraphLcSeg(int win, lcSegVec *segvecp)
{ lcSeg    *segp;
  int       nsg;

  GLgraphSetGC(segvecp->gcp);
  glBegin(GL_LINES);
  for(nsg = 0, segp = segvecp->p; nsg < segvecp->x; nsg++, segp++)
  { glVertex2f(segp->x[0], segp->y[0]);
    glVertex2f(segp->x[1], segp->y[1]);
  }
  glEnd();
  return;
}
/******************************************************************************/
/*  GLgraphStr(int win, cStrVec *strvecp)                                     */
/*                                                                            */
/******************************************************************************/
void      GLgraphStr(int win, cStrVec *strvecp)
{ cStr     *strp;
  char     *cp, c;
  int       fontindex;
  int       n;
  void     *font;
  GLfloat   xof, yof;
  GLfloat   fxscal, fyscal;

  GLgraphSetGC(strvecp->gcp);

  fontindex = strvecp->findex;
  if(fontindex == 1) font = GLUT_STROKE_MONO_ROMAN;
  else               font = GLUT_STROKE_ROMAN;
  fxscal = strvecp->fscalp[0];  fyscal = strvecp->fscalp[1];
  
  for(n = 0, strp = strvecp->p;  n < strvecp->x;  n++, strp++)
  { cp = strp->strp;
    xof = strp->xyp[0];  yof = strp->xyp[1];
    glPushMatrix();
    glTranslatef(xof, yof, 0.0);
     /** def size of GLUT_STROKE..  char is ~ 100, thus * 0.01 to scale them to 1 **/
    glScalef(fxscal * 0.01, fyscal * 0.01, 1.0);
    while( (c = *cp++) ){ glutStrokeCharacter(font, c);}
    glPopMatrix();
  }
  return;
}
/******************************************************************************/
/*  GLgraphCircle(int win, cCircle *p)                                        */
/*                                                                            */
/******************************************************************************/
void      GLgraphCircle(int win, cCircle *p, curveGC *vgcp)
{ double    xc, yc, r, dx;
  double    dstep, dangli, danglf;                        /** variables in DEGREE **/
  double    rstep, rangli, ranglf, rangl, rstepif;        /** variables in RADIAN **/
  double    deg2rad = myPI / 180.;
  int       ix, i;
  short     sign;
  curveGC  *gcp;
/*  static    char      prognamp[] = "utiCurve.2GL::GLgraphCircle"; */

  gcp = p->gcp;                                 /** priority : circle GC if !NULL **/
  if(gcp == NULL) gcp = vgcp;                              /** else circle Vec GC **/
  GLgraphSetGC(gcp);
  
  xc = p->xp[0];  yc = p->xp[1];
  r = p->r;
  dangli = p->anglp[0];  danglf = p->anglp[1];
  rangli = dangli * deg2rad;  ranglf = danglf * deg2rad;
  rangl = fabs(ranglf - rangli);
  sign = (ranglf - rangli >= 0.0)? 1: 0;
  dstep = 2.0;                                /** default angle step is 2 degrees **/
  rstep = dstep * deg2rad;
  dx = rangl/rstep - 0.000001;
  ix = dx;
  rstepif = (rangl - ix * rstep) *0.5;
  if(!sign){ rstepif = - rstepif;  rstep = - rstep;}

/*
fPrintF(stderr,"%s ix=%d, dstep=%f, dstepif=%f, deg2rad=%f\n", 
                                              prognamp, ix, dstep, dstepif, deg2rad);
*/
  if(p->closed) glBegin(GL_LINE_LOOP);
  else          glBegin(GL_LINE_STRIP);
  rangl = rangli;
  glVertex2f(xc + r * cos(rangl), yc + r * sin(rangl));
  for(i = 0;  i < ix;  i++)
  { rangl = (rangli + rstepif + i * rstep);
    glVertex2f(xc + r * cos(rangl), yc + r * sin(rangl));
  }
  rangl = ranglf;
  glVertex2f(xc + r * cos(rangl), yc + r * sin(rangl));
  glEnd();

  return;
}
/******************************************************************************/
/*  GLgraphCircleV(int win, cCircleVec *vecp)                                 */
/*                                                                            */
/******************************************************************************/
void      GLgraphCircleV(int win, cCircleVec *vecp)
{ curveGC  *gcp;
  cCircle  *ccp;
  int       n;

  gcp = vecp->gcp;
  ccp = vecp->p;
  for(n = 0;  n < vecp->x;  n++, ccp++){ GLgraphCircle(win, ccp, gcp);}
  return;
}
/******************************************************************************/
/*  GLgraphEllipse(int win, cEllipse *p)                                      */
/*                                                                            */
/******************************************************************************/
void      GLgraphEllipse(int win, cEllipse *p, curveGC *vgcp)
{ double    xc, yc, a, b, dx;
  double    dstep, dangli, danglf, orient;                /** variables in DEGREE **/
  double    rstep, rangli, ranglf, rangl, rstepif;        /** variables in RADIAN **/
  double    deg2rad = myPI / 180. ;
  int       ix, i;
  short     sign;
  curveGC  *gcp;
/*  static    char      prognamp[] = "utiCurve.2GL::GLgraphEllipse"; */

  gcp = p->gcp;                                /** priority : ellipse GC if !NULL **/
  if(gcp == NULL) gcp = vgcp;                             /** else ellipse Vec GC **/
  GLgraphSetGC(gcp);
  
  xc = p->xp[0];    yc = p->xp[1];
  a = p->axisp[0];  b = p->axisp[1];
  orient = p->orient;                                               /** in DEGREE **/
  dstep = 2.0;
  rstep = dstep * deg2rad;
  dangli = p->anglDp[0];  danglf = p->anglDp[1];
                                                     /** rangli, ranglf in RADIAN **/
  rangli = dangli * deg2rad;  ranglf = danglf * deg2rad;
                              
  rangl = fabs(ranglf - rangli);
  sign = (ranglf - rangli >= 0.0)? 1: 0;
  dx = rangl/rstep  -0.000001;
  ix = dx;
  rstepif = (rangl - ix * rstep) * 0.5;
  if(!sign){ rstepif = - rstepif;  rstep = - rstep;}

/*
fPrintF(stderr,"%s ix=%d, rstep=%f, rstepif=%f\n", prognamp, ix, rstep, rstepif);
*/

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glTranslatef(xc, yc, 0.0);
  glRotatef(orient, 0.0, 0.0, 1.0);                   /** needs orient  in DEGREE **/

  if(p->closed) glBegin(GL_LINE_LOOP);
  else          glBegin(GL_LINE_STRIP);
  rangl = rangli;
  glVertex2f(a * cos(rangl), b * sin(rangl));
  for(i = 0;  i < ix;  i++)
  { rangl = (rangli + rstepif + i * rstep);
    glVertex2f(a * cos(rangl), b * sin(rangl));
  }
  rangl = ranglf;
  glVertex2f(a * cos(rangl), b * sin(rangl));
  glEnd();

  glPopMatrix();
  return;
}
/******************************************************************************/
/*  GLgraphEllipseV(int win, cCircleVec *vecp)                                */
/*                                                                            */
/******************************************************************************/
void      GLgraphEllipseV(int win, cEllipseVec *vecp)
{ curveGC  *gcp;
  cEllipse *ccp;
  int       n;

  gcp = vecp->gcp;
  ccp = vecp->p;
  for(n = 0;  n < vecp->x;  n++, ccp++) GLgraphEllipse(win, ccp, gcp);
  return;
}
/******************************************************************************/
/******************************************************************************/
