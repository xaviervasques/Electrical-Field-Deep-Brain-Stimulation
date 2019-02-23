/*  essai/C/utiCurve.ellipse.c                                                */
/*  Mennessier Gerard                   20011212                              */
/*  Last revised M.G.                   20030624                              */
/*                                                                            */

#include  <stddef.h>
#include  <math.h>                                       /** for function  fabs() **/

#include  "utiMath.constant.def.h"

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

#include  "utiCurve.ellipse.h"
#include  "utiCurve.GC.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cEllipse  *cEllipseAlloc(size_t  nz,char *progcallp)
{ cEllipse *p;
  static    char prognamp[] = "utiCurve.ellipse::cEllipseAlloc";
  static    char form1p[] = "called from %s. malloc failed; cEllipse size nz=%d\n";

  p = (cEllipse *)malloc( nz*sizeof(cEllipse) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  return p;
}
/******************************************************************************/
/*  cEllipse *cEllipseChkRealloc(cEllipse,nzp,neednz,incrnz,progcallp)        */
/*                                                                            */
/* check if there is enough memory allocated for the neednz needed cEllipse   */
/* old and possibly new size  stored as  *nzp                                 */
/* if reallocation is needed, newsize =  neednz + incrnz                      */
/******************************************************************************/
cEllipse  *cEllipseChkRealloc(cEllipse *pi,size_t *nzp,size_t neednz,size_t incrnz,
                                                                     char *progcallp)
{ cEllipse *p;
  size_t    nfz;
  static    char prognamp[] = "utiCurve.ellipse::cEllipseChkRealloc";
  static    char form1p[] = "called from %s. realloc failed;"
                                " initial pointer=%p, desired cEllipse size nz=%d\n";
  static    char form2p[] = "called from %s. malloc failed;"
                                                    " desired cEllipse size nz=%d\n";
  
  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (cEllipse *)malloc( nfz*sizeof(cEllipse) );
fPrintF(stderr,"%s called from %s; realloc p=%p\n", prognamp, progcallp, (void*)pi);
    if(p == NULL) myErr1(-1,stderr,prognamp, form2p, progcallp, neednz);
  }
  else
  { if(neednz  <=  *nzp) return pi;
    p = (cEllipse *)realloc(pi, nfz*sizeof(cEllipse) );
    if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,pi,neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*  void  cEllipsePrint(FILE *bufp,cEllipse *p)                               */
/******************************************************************************/
void      cEllipsePrint(FILE  *bufp, cEllipse *p)
{ static    char form1p[] = 
  "%s  pointer=%p, center={%f, %f}, axis={%f, %f}, orient=%f, angleTs={%f, %f}, angleDs={%f, %f}, closed=%d\n";
  static    char prognamp[] = "utiCurve.ellipse::cEllipsePrint";

  fPrintF(bufp,form1p, prognamp, (void*)p, 
                               p->xp[0],p->xp[1], p->axisp[0],p->axisp[1], p->orient,
                      p->anglTp[0],p->anglTp[1], p->anglDp[0],p->anglDp[1], p->closed);
  return;
}
/******************************************************************************/
/*  void  cEllipseArrayPrint(FILE *bufp,cSeg *p,size_t n)                     */
/******************************************************************************/
void      cEllipseArrayPrint(FILE *bufp,cEllipse *p,size_t n)
{ n++;
  while(--n){ cEllipsePrint(bufp, p);  p++;}
  return;
}
/******************************************************************************/
/*  void  cEllipseZero(cEllipse *p)                                           */
/******************************************************************************/
void      cEllipseZero(cEllipse *p)
{ double   *dp;

  dp = p->xp;       dp[0] = dp[1] = 0.0;
  dp = p->axisp;    dp[0] = dp[1] = 0.0;
  p->orient = 0.0;
  dp = p->anglTp;   dp[0] = dp[1] = 0.0;
  dp = p->anglDp;   dp[0] = dp[1] = 0.0;
  p->gcp = NULL;
  p->closed = 0;
  return;
}
/******************************************************************************/
/*        cEllipseSetCAO(p,xp,axisp,orient)                                   */
/*                                                                            */
/* Set Center, Axises algebraic length, Orientation                           */
/******************************************************************************/
void      cEllipseSetCAO(cEllipse *p, double *xp, double *axisp, double orient)
{ double   *dp;

  dp = p->xp;     *dp++ = xp[0];     *dp = xp[1];
  dp = p->axisp;  *dp++ = axisp[0];  *dp = axisp[1];
  p->orient = orient;
  return ;
}
/******************************************************************************/
/*        cEllipseSetTru(p,anglTrup[2])                                       */
/*                                                                            */
/* Set anglTru in DEGREE                                                      */
/******************************************************************************/
void      cEllipseSetTru(cEllipse *p, double *anglTrup)
{ double   *dp, *axisp;

  dp = p->anglTp;
  *dp++ = anglTrup[0];  *dp = anglTrup[1];
  axisp = p->axisp;
  p->anglDp[0] = cEllipseAnglTru2Dir(axisp, anglTrup[0]);
  p->anglDp[1] = cEllipseAnglTru2Dir(axisp, anglTrup[1]);
  p->closed = cEllipseAngleIsClosed(anglTrup);
  return ;
}
/******************************************************************************/
/*        cEllipseSetDir(p,anglDirp[2])                                       */
/*                                                                            */
/* Set anglDir in DEGREE                                                      */
/******************************************************************************/
void      cEllipseSetDir(cEllipse *p, double *anglDirp)
{ double   *dp, *axisp;

  dp = p->anglDp;
  *dp++ = anglDirp[0];  *dp = anglDirp[1];
  axisp = p->axisp;
  p->anglTp[0] = cEllipseAnglDir2Tru(axisp, anglDirp[0]);
  p->anglTp[1] = cEllipseAnglDir2Tru(axisp, anglDirp[1]);
  p->closed = cEllipseAngleIsClosed(anglDirp);
  return ;
}
/******************************************************************************/
/*        cEllipseSetCAOTru(p,xp,axisp,orient,anglTrup)                       */
/*                                                                            */
/* Set Center, Axises algebraic length, Orientation and Initial/Final angles  */
/*  anglTru[] in DEGREE                                                       */
/******************************************************************************/
void      cEllipseSetCAOTru(cEllipse *p, double *xp, double *axisp,
                                                     double orient, double *anglTrup)
{ double   *dp;

  dp = p->xp;     *dp++ = xp[0];     *dp = xp[1];
  dp = p->axisp;  *dp++ = axisp[0];  *dp = axisp[1];
  p->orient = orient;
  dp = p->anglTp;
  *dp++ = anglTrup[0];  *dp = anglTrup[1];
  p->anglDp[0] = cEllipseAnglTru2Dir(axisp, anglTrup[0]);
  p->anglDp[1] = cEllipseAnglTru2Dir(axisp, anglTrup[1]);
  p->closed = cEllipseAngleIsClosed(anglTrup);
  return ;
}

/******************************************************************************/
/*        cEllipseSetCAODir(p,xp,axisp,orient,anglDirp)                       */
/*                                                                            */
/* Set Center, Axises algebraic length, Orientation and Initial/Final angles  */
/*  anglDir[] in DEGREE                                                       */
/******************************************************************************/
void      cEllipseSetCAODir(cEllipse *p, double *xp, double *axisp,
                                                     double orient, double *anglDirp)
{ double   *dp;

  dp = p->xp;   
  *dp++ = xp[0];     *dp = xp[1];
  dp = p->axisp;
  *dp++ = axisp[0];  *dp = axisp[1];
  p->orient = orient;
  dp = p->anglDp;
  *dp++ = anglDirp[0];  *dp = anglDirp[1];
  p->anglTp[0] = cEllipseAnglDir2Tru(axisp, anglDirp[0]);
  p->anglTp[1] = cEllipseAnglDir2Tru(axisp, anglDirp[1]);
  p->closed = cEllipseAngleIsClosed(anglDirp);
  return ;
}
/******************************************************************************/
/*        cEllipseCopypp(pf, pi)                                              */
/******************************************************************************/
void      cEllipseCopypp(cEllipse *pf, cEllipse *pi)
{ double   *dfp, *dip;
  
  dip = pi->xp;  dfp = pf->xp;        *dfp++ = *dip++;  *dfp = *dip;
  dip = pi->axisp;  dfp = pf->axisp;  *dfp++ = *dip++;  *dfp = *dip;
  pf->orient = pi->orient;
  dip = pi->anglTp;  dfp = pf->anglTp;  *dfp++ = *dip++;  *dfp = *dip;
  dip = pi->anglDp;  dfp = pf->anglDp;  *dfp++ = *dip++;  *dfp = *dip;
  pf->closed = pi->closed;
  pf->gcp = pi->gcp;
  return;
}
/******************************************************************************/
/*        cEllipseCopyp(pf, cc)                                               */
/*                                                                            */
/******************************************************************************/
void      cEllipseCopyp(cEllipse *pf, cEllipse cc)
{ double   *dfp, *dip;
  
  dip = cc.xp;  dfp = pf->xp;        *dfp++ = *dip++;  *dfp = *dip;
  dip = cc.axisp;  dfp = pf->axisp;  *dfp++ = *dip++;  *dfp = *dip;
  pf->orient = cc.orient;
  dip = cc.anglTp;  dfp = pf->anglTp;  *dfp++ = *dip++;  *dfp = *dip;
  dip = cc.anglDp;  dfp = pf->anglDp;  *dfp++ = *dip++;  *dfp = *dip;
  pf->closed = cc.closed;
  pf->gcp = cc.gcp;
  return;
}
/******************************************************************************/
/*        cEllipseAnglDir2Tru(axisp, danglDir)                                */
/*                                                                            */
/*  Compute and return "true" angle (in DEGREE)                               */
/*  from Director angle (in DEGREE)   and axises length                       */
/******************************************************************************/
double    cEllipseAnglDir2Tru(double *axisp, double danglDir)
{ double    x0, x1, raD, raT, daT;
  double    deg2rad = myPI/180.0, rad2deg = 180.0/myPI ;

  raD = danglDir * deg2rad;
  x0 = axisp[0] * cos(raD);
  if(x0 == 0.0) return danglDir;
  x1 = axisp[1] * sin(raD);
  if(x1 == 0.) return danglDir;
  raT = atan2(x1, x0);
  daT = raT * rad2deg;
  daT += 360.0 * cEllipseAngleDetermination360(danglDir);
  return daT;
}
/******************************************************************************/
/*        cEllipseAnglTru2Dir(axisp, danglTru)                                */
/*                                                                            */
/*  Compute and return  Director angle (in DEGREE)                            */
/*  from "true" angle (in DEGREE)   and axises length                         */
/******************************************************************************/
double    cEllipseAnglTru2Dir(double *axisp, double danglTru)
{ double    caD, saD, raT, raD, daD;
  double    deg2rad = myPI/180.0, rad2deg = 180.0/myPI ;

  raT = danglTru * deg2rad;
  caD = cos(raT);
  if(caD == 0.0) return danglTru;
  caD = caD/axisp[0];
  saD = sin(raT);
  if(saD == 0.0) return danglTru;
  saD = saD/axisp[1];

  raD = atan2(saD, caD);
  daD = raD * rad2deg;
  daD += 360.0 * cEllipseAngleDetermination360(danglTru);
  return daD;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       cEllipseAngleDetermination360(double alphaDeg)
{ int       n = 0;

/*
  n = ( (int)(alphaDeg) )/180;
  if(alphaDeg == n * 180.) return n;
*/
  n = (alphaDeg >= 0.0)? ((int)(alphaDeg + 180))/360 : ((int)(alphaDeg - 180))/360;
  return n;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       cEllipseAngleIsClosed(double *alphaDegp)
{ int       n = 0, boolIsClosed = 0;
  double    diff, pre = 0.0000001;

  diff = fabs(alphaDegp[1] - alphaDegp[0]);
  n = (int)(diff + pre) /360;
  if( fabs(diff - n*360.) <= pre ) boolIsClosed = 1;
  return boolIsClosed;
}
/******************************************************************************/
/******************************************************************************/



/******************************************************************************/
/*        cEllipseVecAlloc(nz,progcallp)                                      */
/******************************************************************************/

/******************************************************************************/
/*        cEllipsePVecAlloc(vecp,nz)                                          */
/*                                                                            */
/******************************************************************************/
void      cEllipsePVecAlloc(cEllipseVec *vecp,size_t  nz)
{ 
  vecp->p = cEllipseAlloc(nz,"cEllipsePVecAlloc");
  vecp->z = nz;  vecp->x = 0;
  cEllipseZero(vecp->p);
  return ;
}
/******************************************************************************/
/*        cEllipsePVecRealloc(vecp,neednz,incrnz)                             */
/******************************************************************************/
void      cEllipsePVecRealloc(cEllipseVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = cEllipseChkRealloc(vecp->p, &(vecp->z),neednz,incrnz,
                                                              "cEllipsePVecRealloc");
  return ;
}
/******************************************************************************/
/*        cEllipsePVecFree(vecp)                                              */
/******************************************************************************/
void      cEllipsePVecFree(cEllipseVec *vecp)
{ 
  free(vecp->p);  vecp->p = NULL;  vecp->z = 0;  vecp->x = 0;  return;
}
/******************************************************************************/
/*        cEllipseVecFree(vecp)                                               */
/******************************************************************************/
void      cEllipseVecFree(cEllipseVec *vecp)
{ 
  if(vecp->p != NULL) cEllipsePVecFree(vecp);
  free(vecp);  return;
}
/******************************************************************************/
/*        cEllipseVecPrint(vecp)                                              */
/******************************************************************************/
void      cEllipseVecPrint(FILE  *bufp, cEllipseVec *vecp)
{ cEllipse  *p;
  static    char form1[] = "%s  cEllipseVec.p=%p, vec.z=%d, vec.x=%d \n";
  static    char prognamp[] = "utiCurve.ellipse::cEllipseVecPrint";

  p = vecp->p;
  fPrintF(bufp,form1, prognamp, p,vecp->z,vecp->x);
  if(p == NULL) return;
  cEllipseArrayPrint(bufp,p,vecp->x);
  fPrintF(bufp, "\n");  return ;
}
/******************************************************************************/
/*        cEllipseVecInc1(vecp,y)                                             */
/* add 1 cEllipse into a cEllipseVec structure                                */
/******************************************************************************/
void      cEllipseVecInc1(cEllipseVec *vecp, cEllipse cc)
{ cEllipse  *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cEllipsePVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  cEllipseCopyp(p, cc);
  return ;
}
/******************************************************************************/
/*        cEllipseVecInc1p(vecp,yp)                                           */
/* add 1 *cEllipse into a cEllipseVec structure                               */
/******************************************************************************/
void      cEllipseVecInc1p(cEllipseVec *vecp, cEllipse *ccp)
{ cEllipse  *p;
  size_t    newz,oldx;

  oldx = vecp->x ; newz = oldx+1;
  cEllipsePVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  cEllipseCopypp(p, ccp);
  return;
}
/******************************************************************************/
/*        cEllipseVecInc1DataTru(vecp, xp, axisp, orient, anglTp)             */
/* add 1 cEllipse into a cEllipseVec structure                                */
/*               from center, axis, radius, angles                            */
/******************************************************************************/
void      cEllipseVecInc1DataTru(cEllipseVec *vecp, double *xp, double *axisp,
                                                        double orient, double *anglTp)
{ cEllipse  *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cEllipsePVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  cEllipseSetCAOTru(p, xp, axisp, orient, anglTp);
  p->gcp = NULL;
  return;
}
/******************************************************************************/
/******************************************************************************/
