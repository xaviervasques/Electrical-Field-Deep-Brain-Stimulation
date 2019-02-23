/*  URMAE/orientHaut/linear4.GL.V1/pallidus.draw.c                            */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20020103                                */

#include  <stddef.h>
#include  <math.h>
#include  <GL/glut.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utistdStr.h"
#include  "utiMem.h"
#include  "utiStr.h"
#include  "utiAlloc.h"
#include  "utiVecInt.h"
#include  "utiVecDbl.h"
#include  "utiBookChr.h"
#include  "utiLec.h"
#include  "utiClassRangeDbl.h"
#include  "utiClassRangeInt.h"
#include  "utiCurve.def.h"
#include  "utiCurve.level.h"
#include  "utiCurve.poly.h"
#include  "utiCurve.set.h"

#include  "pallidus.draw.statalloc.h"
#include  "pallidus.geom.glob.h"
#include  "pallidus.draw.h"

/**
      S  scanner Frame
      ST scanner Translated Frame, same origin as electrode Frame
      E  electrode Frame
**/

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      pallidusInitDraw()
{ 

fPrintF(stderr, "            Computing pallidus E RZ section; psi=0\n");
  pallidusEComputeRZSection(&pallidusERZSegV, 0.0, 0);
  lcSegVecPrint(stderr, &pallidusERZSegV);
fPrintF(stderr, "            Computing pallidus E XY section; z=0\n");
  pallidusEComputeXYSectionRotated(&pallidusEXYSegV, 0.0, 0);
  lcSegVecPrint(stderr, &pallidusEXYSegV);
fPrintF(stderr, "\n");

fPrintF(stderr, "            Computing pallidus ST XZ section; psi=0\n");
  pallidusSTComputeRZSection(&pallidusSTRZSegV, 0.0, 0);
  lcSegVecPrint(stderr, &pallidusSTRZSegV);
fPrintF(stderr, "            Computing pallidus ST XY section; z=0\n");
  pallidusSTComputeXYSection(&pallidusSTXYSegV, 0.0, 0);
  lcSegVecPrint(stderr, &pallidusSTXYSegV);
fPrintF(stderr, "\n");

/*fPrintF(stderr, "            Showing pallidus all ST XY section for data z\n");
  pallidusSTallXYdataZsection(&pallidusSTallXYSetV);
  */
  return;
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
lcSegVec *pallidusERZGetlcSeg()
{ 
  return  &pallidusERZSegV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
lcSegVec *pallidusEXYGetlcSeg()
{ 
  return  &pallidusEXYSegV;
}
/******************************************************************************/
/*                                                                            */
/*  Electrode frame                                                           */
/*  draw  rho z Section at fixed y = 0.0, AFTER psi rotation                  */
/*  Use psi rotatation rather than y translation to keep electrode axis       */
/******************************************************************************/
void      pallidusEComputeRZSection(lcSegVec *lcsegvecp, double y, int iy)
{ int       n0, n1;
  int      *ip, i, ix;
  int       lp[3] = {0,2,1};
  short     closed = 1;
  double   *d0p, *d1p;
  static    curveGC   gc;
  static    char      prognamp[] = "pallidus.draw::pallidusEComputeRZSection";

  if(eltrdSPosVd.x == 0) return;

  gc.lineWidth = 3;
  lcsegvecp->gcp = &gc;
  lcsegvecp->x = 0;

  ip = pallidusNptPlaneVi.p;
  ix = pallidusNptPlaneVi.x;                                 /** number of planes **/
  n0 = *ip++;                                  /** number of pints in first plane **/
  d0p = pallidusEBoundRotatedVd.p;
  for(i = 1;  i < ix;  i++)
  { n1 = *ip++;
    if(n1 != n0) fPrintF(stderr,"%s, WARNING n0=%d, n1=%d\n", prognamp, n0, n1);
    d1p = d0p + 3 * n0;
    lcSegAddSFromQuadrangles(d0p, d1p, (size_t)n0, closed, lp, y, iy, lcsegvecp);
    d0p = d1p;  n0 = n1;
  }
lcSegVecPrint(stderr, lcsegvecp);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  Electrode frame                                                           */
/*  x y  Section at fixed altitude z, AFTER psi rotation                      */
/*  Use psi rotatation rather than y translation to keep electrode axis       */
/******************************************************************************/
void      pallidusEComputeXYSectionRotated(lcSegVec *lcsegvecp, double z, int iz)
{ int       n0, n1;
  int      *ip, i, ix;
  int       lp[3] = {0,1,2};
  short     closed = 1;
  double   *d0p, *d1p;
  static    curveGC   gc;
  static    char    prognamp[] = "pallidus.draw::pallidusEComputeXYSectionRotated";

  if(eltrdSPosVd.x == 0) return;

  gc.lineWidth = 3;
  lcsegvecp->gcp = &gc;
  lcsegvecp->x = 0;

  ip = pallidusNptPlaneVi.p;
  ix = pallidusNptPlaneVi.x;                                 /** number of planes **/
  n0 = *ip++;                                  /** number of pints in first plane **/
  d0p = pallidusEBoundRotatedVd.p;
  for(i = 1;  i < ix;  i++)
  { n1 = *ip++;
    if(n1 != n0) fPrintF(stderr,"%s, WARNING n0=%d, n1=%d\n", prognamp, n0, n1);
    d1p = d0p + 3 * n0;
    lcSegAddSFromQuadrangles(d0p, d1p, (size_t)n0, closed, lp, z, iz, lcsegvecp);
    d0p = d1p;  n0 = n1;
  }
lcSegVecPrint(stderr, lcsegvecp);
  return;
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
lcSegVec *pallidusSTRZGetlcSeg()
{
  return  &pallidusSTRZSegV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
lcSegVec *pallidusSTXYGetlcSeg()
{
  return  &pallidusSTXYSegV;
}
/******************************************************************************/
/*                                                                            */
/*  Scanner Translated frame                                                  */
/*  x z  Sections at given value y (usually 0)                                */
/******************************************************************************/
void      pallidusSTComputeRZSection(lcSegVec *lcsegvecp, double y, int iy)
{ int       n0, n1;
  int      *ip, i, ix;
  int       lp[3] = {0,2,1};
  short     closed = 1;
  double   *d0p, *d1p;
  static    curveGC   gc;
  static    char    prognamp[] = "pallidus.draw::pallidusSTComputeXZSection";

  if(eltrdSPosVd.x == 0) return;

  gc.lineWidth = 3;
  lcsegvecp->gcp = &gc;
  lcsegvecp->x = 0;

  ip = pallidusNptPlaneVi.p;
  ix = pallidusNptPlaneVi.x;                                 /** number of planes **/
  n0 = *ip++;                                  /** number of pints in first plane **/
  d0p = pallidusSTBoundRotatedVd.p;
  for(i = 1;  i < ix;  i++)
  { n1 = *ip++;
    if(n1 != n0) fPrintF(stderr,"%s, WARNING n0=%d, n1=%d\n", prognamp, n0, n1);
    d1p = d0p + 3 * n0;
    lcSegAddSFromQuadrangles(d0p, d1p, (size_t)n0, closed, lp, y, iy, lcsegvecp);
    d0p = d1p;  n0 = n1;
  }
lcSegVecPrint(stderr, lcsegvecp);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  Scanner Translated frame                                                  */
/*  x y  Sections at altitude z                                               */
/******************************************************************************/
void      pallidusSTComputeXYSection(lcSegVec *lcsegvecp, double z, int iz)
{ int       n0, n1;
  int      *ip, i, ix;
  int       lp[3] = {0,1,2};
  short     closed = 1;
  double   *d0p, *d1p;
  static    curveGC   gc;
  static    char    prognamp[] = "pallidus.draw::pallidusSTComputeXYSection";

  if(eltrdSPosVd.x == 0) return;

  gc.lineWidth = 3;
  lcsegvecp->gcp = &gc;
  lcsegvecp->x = 0;

  ip = pallidusNptPlaneVi.p;
  ix = pallidusNptPlaneVi.x;                                 /** number of planes **/
  n0 = *ip++;                                  /** number of pints in first plane **/
  d0p = pallidusSTBoundVd.p;
  for(i = 1;  i < ix;  i++)
  { n1 = *ip++;
    if(n1 != n0) fPrintF(stderr,"%s, WARNING n0=%d, n1=%d\n", prognamp, n0, n1);
    d1p = d0p + 3 * n0;
    lcSegAddSFromQuadrangles(d0p, d1p, (size_t)n0, closed, lp, z, iz, lcsegvecp);
    d0p = d1p;  n0 = n1;
  }
/* lcSegVecPrint(stderr, lcsegvecp); */
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cSetVec  *pallidusSTXYallVGetSet()
{ 
  return  &pallidusSTallXYSetV;
}
/******************************************************************************/
/*                                                                            */
/*  Scanner Translated frame                                                  */
/*  Sections in all data z planes                                             */
/******************************************************************************/
void      pallidusSTallXYdataZsection(cSetVec *csetvecp)
{ int       n;
  int      *ip, i, ix,  j;
  short     closed = 1;
  double   *dp;
  static    curveGC   gcColorp[20];           /** ASSUME <= 20  z planes for DATA **/
  static    curveGC  *gcp;
  int       lineWidth = 3;
  float     color[3] = {1.0 , 0.0, 0.0};
  static    size_t    polypz = 0;
  static    cPolyp    polyp = NULL, polycp;
  static    char    prognamp[] = "pallidus.draw::pallidusSTallXYdataZsection";

  if(eltrdSPosVd.x == 0) return;

  ix = pallidusNptPlaneVi.x;                                 /** number of planes **/
  if(polyp == NULL){ polypz = ix;  polyp = cPolyAlloc(polypz, prognamp);}
  else           { polyp = cPolyChkRealloc(polyp, &polypz, (size_t)ix, 1, prognamp);}
  ip = pallidusNptPlaneVi.p;
  dp = pallidusSTBoundVd.p;  polycp = polyp;  gcp = gcColorp;
  for(i = 0;  i < ix;  i++, polycp++, gcp++)
  { n = *ip++;                              /** number of points in current plane **/
    gcp->lineWidth = lineWidth;
    gcp->color[0] = 1.  - i * 0.8 /ix;
    gcp->color[1] = color[1];  gcp->color[2] = color[2];
    polycp->gcp = gcp;
    polycp->closed = closed;
    for(j = 1;  j <= n;  j++){ cPolyPInc1p(polycp, dp);  dp += 3;}
    cSetVecInc1ptype(csetvecp, polycp, MY_CPOLY);
  }
  return;
}

/******************************************************************************/
/*                                                                            */
/*                    pure  GL  routines                                      */
/*                                                                            */
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/*  Sections in all data z planes                                             */
/******************************************************************************/
void      pallidusSTallXYZaxisGL()
{ int       i, in, ix;
  double    x, y, z, ticksize = 0.05;
  double    dmin = -7.1, dmax = +7.1;
  int       lineWidth = 1;
  float     color[3] = {0.0 , 0.0, 0.0};
/*  static    char    prognamp[] = "pallidus.draw::pallidusSTallXYZaxisGL";  */

  glColor3f(color[0], color[1], color[2]);
  glLineWidth(lineWidth);
  in = dmin;  ix = dmax;

  glBegin(GL_LINES);
                                                       /** x axis at y = 0, z = 0 **/
  y = 0.0;  z = 0.0;
  glVertex3f(dmin, y, z);  glVertex3f(dmax, y, z);
                                                                        /** ticks **/
  for(i = in;  i <= ix;  i++)
  { x = i * 1.0;
    glVertex3f(x, y, z);   glVertex3f(x, y - ticksize, z);
  }
                                                       /** y axis at x = 0, z = 0 **/
  x = 0.0;  z = 0.0;
  glVertex3f(x, dmin, z);  glVertex3f(x, dmax, z);
                                                                        /** ticks **/
  for(i = in;  i <= ix;  i++)
  { y = i * 1.0;
    glVertex3f(x, y, z);   glVertex3f(x - ticksize, y, z);
  }
                                                       /** z axis at x = 0, y = 0 **/
  x = 0.0;  y = 0.0;
  glVertex3f(x, y, dmin);  glVertex3f(x, y, dmax);
                                                                      /** x ticks **/
  for(i = in;  i <= ix;  i++)
  { z = i * 1.0;
    glVertex3f(x, y, z);   glVertex3f(x - ticksize, y, z);
  }
                                                                      /** y ticks **/
  for(i = in;  i <= ix;  i++)
  { z = i * 1.0;
    glVertex3f(x, y, z);   glVertex3f(x, y - ticksize, z);  
  }
  glEnd();
  return;
}
/******************************************************************************/
/*                                                                            */
/*  Sections in all data z planes                                             */
/******************************************************************************/
void      pallidusSTallXYZaxisGLlist()
{ int       Pallidus3DaxisListIndex = 1;
/*  static    char    prognamp[] = "pallidus.draw::pallidusSTallXYZaxisGLlist"; */

  glNewList(Pallidus3DaxisListIndex, GL_COMPILE);
  pallidusSTallXYZaxisGL();
  glEndList();

  return;
}
/******************************************************************************/
/*                                                                            */
/*  Sections in all data z planes                                             */
/******************************************************************************/
void      pallidusSTallXYdataZsectionGL()
{ int       n, ni;
  int      *ip, i, ix,  j;
  double   *dp;
  int       lineWidth = 3;
  float     color[3] = {1.0 , 0.0, 0.0};
  static    char    form1p[] = "  %s\n";
  static    char    form2p[] = "  plane index=%d, number of points=%d\n";
  static    char    form3p[] = "  point= %f  %f  %f\n";
  static    char    prognamp[] = "pallidus.c::pallidusSTallXYdataZsectionGL";

  if(eltrdSPosVd.x == 0) return;
fPrintF(stderr,form1p, prognamp);

  glLineWidth(lineWidth);
  ix = pallidusNptPlaneVi.x;                                 /** number of planes **/
  ip = pallidusNptPlaneVi.p;
  ni = *ip;                                        /** number of points per plane **/
  dp = pallidusSTBoundVd.p;
  for(i = 0;  i < ix;  i++, ip++)
  { n = *ip;
    if(n != ni) fPrintF(stderr,"%s, WARNING ni=%d NotEqual n=%d\n", prognamp, ni, n);
fPrintF(stderr,form2p, i, n);
    color[0] = 0.2 + i * 0.8 /ix;
    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_LINE_LOOP);
    for(j = 1;  j <= n;  j++, dp +=3)
    { glVertex3f(*dp, *(dp+1), *(dp+2) );
fPrintF(stderr,form3p, *dp, *(dp+1), *(dp+2) );
    }
    glEnd();
  }

  glPointSize(7.0);
  glColor3f(0.0, 0.0, 1.0);
  ip = pallidusNptPlaneVi.p;
  dp = pallidusSTBoundVd.p;
  glBegin(GL_POINTS);
  for(i = 0;  i < ix;  i++, ip++)
  { n = *ip;
    for(j = 1;  j <= n;  j++, dp +=3){ glVertex3f(*dp, *(dp+1), *(dp+2) );}
  }
  glEnd();
  return;
}
/******************************************************************************/
/*                                                                            */
/*  Sections in all data z planes                                             */
/******************************************************************************/
void      pallidusSTallXYdataZsectionGLlist()
{ int       Pallidus3DsectionListIndex = 2;
/*  static    char    prognamp[] = "pallidus.draw::pallidusSTallXYdataZsectionGLlist"; */

  glNewList(Pallidus3DsectionListIndex, GL_COMPILE);
  pallidusSTallXYdataZsectionGL();
  glEndList();

  return;
}
/******************************************************************************/
/*                                                                            */
/*  Quadrangles strips in all data z planes                                   */
/******************************************************************************/
void      pallidusSTallXYdataZquadGL()
{ int       n, ni;
  int      *ip, i, ix,  j;
  double   *dep, *dop, *decp, *docp;
  int       lineWidth = 3;
  float     colorf[3] = {0.0 , 1.0, 0.0};
  static    char    prognamp[] = "pallidus.draw::pallidusSTallXYdataZquadGL";

  if(eltrdSPosVd.x == 0) return;

  glLineWidth(lineWidth);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  ix = pallidusNptPlaneVi.x;                                 /** number of planes **/
  ip = pallidusNptPlaneVi.p;
  ni = *ip;                                        /** number of points per plane **/
  dep = pallidusSTBoundVd.p;  dop = dep + 3 * ni;
  for(i = 1;  i < ix;  i++, ip++)
  { n = *ip;
    if(n != ni) fPrintF(stderr,"%s, WARNING ni=%d NotEqual n=%d\n", prognamp, ni, n);
    colorf[1] =  0.2 + i * 0.8/ix;
    glColor3f(colorf[0], colorf[1], colorf[2]);

    decp = dep;  docp = dop;
    glBegin(GL_QUAD_STRIP);
    for(j = 1;  j <= n;  j++, decp +=3, docp +=3)
    { glVertex3f(*decp, *(decp+1), *(decp+2) ); 
      glVertex3f(*docp, *(docp+1), *(docp+2) );
    }
                                                              /** close the strip **/
    glVertex3f(*dep, *(dep+1), *(dep+2) ); 
    glVertex3f(*dop, *(dop+1), *(dop+2) );
    dep = dop;  dop += 3 * ni;
    glEnd();
  }
  return;
}
/******************************************************************************/
/*                                                                            */
/*  Quadrangles strips in all data z planes                                   */
/******************************************************************************/
void      pallidusSTallXYdataZquadGLlist()
{ int       Pallidus3DquadListIndex = 3;
/*  static    char    prognamp[] = "pallidus.draw::pallidusSTallXYdataZquadGLlist"; */

  glNewList(Pallidus3DquadListIndex, GL_COMPILE);
  pallidusSTallXYdataZquadGL();
  glEndList();

  return;
}
/******************************************************************************/
/*                                                                            */
/*  Filled Quadrangles strips in all data z planes                            */
/******************************************************************************/
void      pallidusSTallXYdataZquadfillGL()
{ int       n, ni;
  int      *ip, i, ix,  j;
  double   *dep, *dop, *decp, *docp;
  int       lineWidth = 3;
  float     colorb[3] = {1.0 , 0.0, 0.0};
  float     colorf[3] = {0.0 , 1.0, 0.0};
  static    char    prognamp[] = "pallidus.draw::pallidusSTallXYdataZquadfillGL";

  if(eltrdSPosVd.x == 0) return;

  glLineWidth(lineWidth);
/*
  glPolygonMode(GL_BACK,  GL_LINE);
  glPolygonMode(GL_BACK,  GL_FILL);
  glPolygonMode(GL_FRONT, GL_FILL);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
*/

/*  glEnable(GL_CULL_FACE); */
/*  glCullFace(GL_BACK); */

  ix = pallidusNptPlaneVi.x;                                 /** number of planes **/
  ip = pallidusNptPlaneVi.p;
  ni = *ip;                                        /** number of points per plane **/
  dep = pallidusSTBoundVd.p;  dop = dep + 3 * ni;
  for(i = 1;  i < ix;  i++, ip++)
  { n = *ip;
    if(n != ni) fPrintF(stderr,"%s, WARNING ni=%d NotEqual n=%d\n", prognamp, ni, n);
/*
    glPolygonMode(GL_FRONT, GL_FILL);
    colorb[0] = (0.5 + i * 0.5/ix);
    glColor3f(colorb[0], colorb[1], colorb[2]);
*/
    glPolygonMode(GL_BACK,  GL_FILL);
    colorf[1] = (0.5 + i * 0.5/ix);
    glColor3f(colorf[0], colorf[1], colorf[2]);

    glPolygonMode(GL_FRONT, GL_FILL);
    colorb[0] = (0.5 + i * 0.5/ix);
    glColor3f(colorb[0], colorb[1], colorb[2]);

    decp = dep;  docp = dop;
    glBegin(GL_QUAD_STRIP);
    for(j = 1;  j <= n;  j++, decp +=3, docp +=3)
    { glVertex3f(*decp, *(decp+1), *(decp+2) ); 
      glVertex3f(*docp, *(docp+1), *(docp+2) );
    }
                                                              /** close the strip **/
    glVertex3f(*dep, *(dep+1), *(dep+2) ); 
    glVertex3f(*dop, *(dop+1), *(dop+2) );
    dep = dop;  dop += 3 * ni;
    glEnd();
  }
  return;
}
/******************************************************************************/
/*                                                                            */
/*  Filled Quadrangles strips in all data z planes                            */
/******************************************************************************/
void      pallidusSTallXYdataZquadfillGLlist()
{ int       Pallidus3DquadListIndex = 4;
/*  static    char    prognamp[] = "pallidus.draw::pallidusSTallXYdataZquadfillGLlist"; */

  glNewList(Pallidus3DquadListIndex, GL_COMPILE);
  pallidusSTallXYdataZquadfillGL();
  glEndList();

  return;
}
/******************************************************************************/
/******************************************************************************/
