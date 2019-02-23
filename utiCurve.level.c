/*  essai/C/utiCurve.level.c                                                  */
/*  Mennessier Gerard                   20010502                              */
/*  Last revised M.G.                   20040331                              */
/*                                                                            */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"
#include  "utiMath.type.def.h"
#include  "utiCurve.level.h"

static myBOOL  bzcmp[4], bcross[4];
static double  xycross[4], xcrossp[4], ycrossp[4];

static char    srcfilenamp[] = "utiCurve.level";
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
lcSeg     *lcSegAlloc(size_t  nz,char *progcallp)
{ lcSeg     *p;
  static char    form1p[] = "called from %s. malloc failed; lcSeg size nz=%d\n";
  static char    prognamp[] = "lcSegAlloc";

  p = (lcSeg *)malloc( nz*sizeof(lcSeg) );
  if(p == NULL) myErr2(-1, stderr, srcfilenamp, prognamp,  form1p, progcallp, nz);
  return p;
}
/******************************************************************************/
/*  lcSeg *lcSegChkRealloc(lcSegp,nzp,neednz,incrnz,progcallp)                */
/*                                                                            */
/* check if there is enough memory allocated for the neednz needed lcSeg      */
/* old and possibly new size  stored as  *nzp                                 */
/* if reallocation is needed, newsize =  neednz + incrnz                      */
/******************************************************************************/
lcSeg     *lcSegChkRealloc(lcSeg *pi,size_t  *nzp,size_t neednz,size_t incrnz,
                                                                     char *progcallp)
{ lcSeg    *p;
  size_t    nfz;
  static char    form1p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired lcSeg size nz=%d\n";
  static char    form2p[]= "called from %s. malloc failed;"
                                                       " desired lcSeg size nz=%d\n";
  static char    prognamp[] = "utiCurve.level::lcSegChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (lcSeg *)malloc( nfz*sizeof(lcSeg) );
fPrintF(stderr,"%s called from %s; realloc p=%p\n", prognamp, progcallp, (void*)pi);
    if(p == NULL) myErr2(-1, stderr, srcfilenamp, prognamp, form2p, progcallp, neednz);
  }
  else
  { if(neednz  <=  *nzp) return pi;
    p = (lcSeg *)realloc(pi, nfz*sizeof(lcSeg) );
    if(p == NULL) myErr2(-1,stderr, srcfilenamp,prognamp, form1p,progcallp,pi,neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*  void  lcSegPrint(FILE *bufp,lcSeg *p)                                     */
/******************************************************************************/
void      lcSegPrint(FILE *bufp,lcSeg *p)
{ 
  fPrintF(bufp,"P1=(%f, %f), P2=(%f, %f), level number=%d",
                                          p->x[0], p->y[0], p->x[1], p->y[1], p->iz);
  return;
}
/******************************************************************************/
/*  void  lcSegArrayPrint(FILE *bufp,lcSeg *p,size_t n)                       */
/******************************************************************************/
void      lcSegArrayPrint(FILE *bufp,lcSeg *p,size_t n)
{ n++;
  while(--n)
  { fPrintF(bufp,"P1=(%f, %f), P2=(%f, %f), level number=%d\n",
                                          p->x[0], p->y[0], p->x[1], p->y[1], p->iz);
    p++;
  }
  return;
}
/******************************************************************************/
/*  void  lcSegZero(lcSeg *p)                                                 */
/******************************************************************************/
void      lcSegZero(lcSeg *p)
{ 
  p->x[0] = 0.; p->x[1] = 0.;  p->y[0] = 0.; p->y[1] = 0.;  p->iz = 0;  return;
}
/******************************************************************************/
/*  void  lcSegCpy1p(lcSeg *sfp,lcSeg *sip)                                   */
/*                                                                            */
/* copy 1 lcSeg  *sip into  *sfp                                              */
/******************************************************************************/
void      lcSegCpy1p(lcSeg *sfp,lcSeg *sip)
{ double    *dip, *dfp;

  dip = sip->x;  dfp = sfp->x;  *dfp++ = *dip++;  *dfp = *dip;
  dip = sip->y;  dfp = sfp->y;  *dfp++ = *dip++;  *dfp = *dip;
  sfp->iz = sip->iz;  return;
}
/******************************************************************************/



/******************************************************************************/
/*        lcZrangeZero(lcZrange *p)                                           */
/******************************************************************************/
void      lcZrangeZero(lcZrange *p)
{
  p->zo  = 0.0;
  p->zpa = 1.0;
  p->zn  = 1.0001;  p->zx = -1.0001;
  p->in = 1;  p->ix = -1;
  return;
}
/******************************************************************************/
/*        lcZrangePrint(FILE *bufp, lcZrange *p)                              */
/******************************************************************************/
void      lcZrangePrint(FILE *bufp, lcZrange *p)
{ static char    form1p[] = "%s::%s  LevelCurve Z range\n";
  static char    form2p[] = "  Z origin=%f, step=%f, min z=%f, max z=%f" \
                                                   ", min index=%d, max index=%d \n";
  static char    prognamp[] = "lcZrangePrint";

  fPrintF(bufp, form1p, srcfilenamp, prognamp);
  fPrintF(bufp, form2p, p->zo, p->zpa, p->zn, p->zx, p->in, p->ix);
  return;
}
/******************************************************************************/



/******************************************************************************/
/*        lcSegVecAlloc(nz,progcallp)                                         */
/******************************************************************************/

/******************************************************************************/
/*        lcSegPVecAlloc(vecp,nz)                                             */
/******************************************************************************/
void      lcSegPVecAlloc(lcSegVec *vecp,size_t  nz)
{ 
  vecp->p = lcSegAlloc(nz,"lcSegPVecAlloc");
  vecp->z = nz;  vecp->x = 0;
  lcSegZero(vecp->p);
  return ;
}
/******************************************************************************/
/*        lcSegPVecRealloc(vecp,neednz,incrnz)                                */
/******************************************************************************/
void      lcSegPVecRealloc(lcSegVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = lcSegChkRealloc(vecp->p, &(vecp->z), neednz, incrnz, "lcSegPVecRealloc");
  return ;
}
/******************************************************************************/
/*        lcSegPVecFree(vecp)                                                 */
/******************************************************************************/
void      lcSegPVecFree(lcSegVec *vecp)
{ 
  free(vecp->p);  vecp->p = NULL;  vecp->z = 0;  vecp->x = 0;  return ;
}
/******************************************************************************/
/*        lcSegVecFree(vecp)                                                  */
/******************************************************************************/
void      lcSegVecFree(lcSegVec *vecp)
{ 
  if(vecp->p != NULL) lcSegPVecFree(vecp);
  free(vecp);  return ;
}
/******************************************************************************/
/*        lcSegVecPrint(vecp)                                                 */
/******************************************************************************/
void      lcSegVecPrint(FILE  *bufp, lcSegVec *vecp)
{ lcSeg    *p;
  static    char    form1p[] = "%s::%s lcSegVec.p=%p, vec.z=%d, vec.x=%d \n" ;
  static    char    prognamp[] = "lcSegVecPrint";

  p = vecp->p;
  fPrintF(bufp, form1p, srcfilenamp, prognamp, p, vecp->z, vecp->x);
  if(p == NULL) return;
  lcSegArrayPrint(bufp, p, vecp->x);
  fPrintF(bufp, "\n");  return ;
}
/******************************************************************************/
/*        lcSegVecInc1(vecp,y)                                                */
/* add 1 lcSeg into a lcSegVec structure                                      */
/******************************************************************************/
void      lcSegVecInc1(lcSegVec *vecp,lcSeg s)
{ lcSeg    *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  lcSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  p->x[0] = s.x[0];  p->x[1] = s.x[1];  p->y[0] = s.y[0];  p->y[1] = s.y[1];
  p->iz = s.iz;
  return ;
}
/******************************************************************************/
/*        lcSegVecInc1p(vecp,yp)                                              */
/* add 1 *lcSeg into a lcSegVec structure                                     */
/******************************************************************************/
void      lcSegVecInc1p(lcSegVec *vecp,lcSeg *sp)
{ lcSeg    *p;
  size_t    newz,oldx;

  oldx = vecp->x ; newz = oldx+1;
  lcSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;  lcSegCpy1p(p,sp);  return ;
}
/******************************************************************************/
/*        lcSegVecInc1PointData(vecp,*p1,*p2,iz)                              */
/* add 1 segment, point p1 to p2, level number iz into a lcSegVec structure   */
/* *p1 = (x1,y1) , *p2 = (x2,y2)                                              */
/******************************************************************************/
void      lcSegVecInc1PointData(lcSegVec *vecp,double *p1p, double *p2p, int iz)
{ lcSeg    *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  lcSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  p->x[0] = *p1p++;  p->x[1] = *p2p++;  p->y[0] = *p1p;  p->y[1] = *p2p;
  p->iz = iz;
  return ;
}
/******************************************************************************/
/*        lcSegVecInc1CoorData(vecp,x0,y0,x1,y1,iz)                           */
/* add 1 segment, point (x0,y0) to (x1,y1), level number iz                   */
/*                                                into a lcSegVec structure   */
/******************************************************************************/
void      lcSegVecInc1CoorData(lcSegVec *vecp,double x0, double y0,
                                                        double x1, double y1, int iz)
{ lcSeg    *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  lcSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  p->x[0] = x0;  p->x[1] = x1;  p->y[0] = y0;  p->y[1] = y1;
  p->iz = iz;
  return ;
}
/******************************************************************************/
/*        lcSegVecInc1CoorPData(vecp,xp,yp,iz)                                */
/* add 1 segment, point (xp[0],yp[0]) to (xp[1],yp[1]), level number iz       */
/*                                                into a lcSegVec structure   */
/******************************************************************************/
void      lcSegVecInc1CoorPData(lcSegVec *vecp,double *xp, double *yp, int iz)
{ lcSeg    *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  lcSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  p->x[0] = *xp++;  p->x[1] = *xp;  p->y[0] = *yp++;  p->y[1] = *yp;
  p->iz = iz;
  return ;
}
/******************************************************************************/



/******************************************************************************/
/*        lcSegSFromLineInit                                                  */
/*                                                                            */
/*   y1         node 1 (level z1)                                             */
/*                  |                                                         */
/*                 |                                                          */
/*           link 0                                                           */
/*               |                                                            */
/*              |                                                             */
/*   y0   node 0       (level z0)                                             */
/*                                                                            */
/*             x0   x1                                                        */
/******************************************************************************/
void      lcSegSFromLineInit(double x0, double x1, double *yp, 
                                                      double z0, double z1, double z)
{ 
  bzcmp[0] = (z >= z0);  bzcmp[1] = (z >= z1);
  bcross[0] = (bzcmp[0] != bzcmp[1]) ;

/** fPrintF(stderr,"INIT: P0=(%f, %f, %f), P1=(%f, %f, %f), z=%f\n", 
                                                    x0,yp[0],z0, x1,yp[1],z1, z); **/
  return ;
}
/******************************************************************************/
/*        lcSegAddSFromLineRectangle                                          */
/*                                                                            */
/*    y1       1---1---2      levels *zp1                                     */
/*             |       |                                                      */
/*             |       |                                                      */
/*             0       2                                                      */
/*             |       |                                                      */
/*             |       |                                                      */
/*    y0       0---3---3      levels *zp0                                     */
/*                                                                            */
/*             x0     x1                                                      */
/*                                                                            */
/*  nx   rectangles                                                           */
/*  nx+1 points xp[0], xp[nx]                                                 */
/******************************************************************************/
int       lcSegAddSFromLineRectangle(double *xp, size_t nx, double *yp,
                        double *z0p, double *z1p, double z, int iz, lcSegVec *lcsegp)
{ int       ok = 1;
  int       ix, ncross;
  double    x0, x1, deltax, y0, y1, deltay;
  double    deltaz, zweight1 = 0.0, zmean, zcol, zsum0, zsum1;
  double    z00, z01, z10, z11;
  static    char prognamp[] = "utiCurve.level::lcSegAddSFromLineRectangle";

  if(nx <= 0) return ok;

  y0 = *yp++;  y1 = *yp;  deltay = y1 - y0;
/* fPrintF(stderr,"y0=%f, y1=%f\n", y0,y1); */
  lcSegPVecRealloc(lcsegp, lcsegp->x + nx, nx/20);
                                                          /** Loop Initialisation **/
  x1 = *xp++;
  z01 = *z0p++;  z11 = *z1p++;  deltaz = z11 - z01;  zsum1 = z11 + z01;
  bzcmp[2] = (z >= z11);  bzcmp[3] = (z >= z01);
  bcross[2] = (bzcmp[2] != bzcmp[3]);
  if(bcross[2])
  { deltaz = z11 - z01;   zweight1 = (z - z01)/deltaz;
    xycross[2] = y0 + deltay * zweight1;
  }
                                                                       /** x Loop **/
  for(ix = 0;  ix < nx;  ix++)
  { x0 = x1;  x1 = *xp++;  deltax = x1 - x0;
    z00 = z01;  z10 = z11;
    z01 = *z0p++;  z11 = *z1p++;
    zsum0 = zsum1;  zsum1 = z11 + z01;
    zmean = (zsum0 + zsum1)/4 ;

    bzcmp[0] = bzcmp[3];
    bzcmp[1] = bzcmp[2];
    bcross[0] = bcross[2];
    xycross[0] = xycross[2];
                                   /** bzcmp[0] and bzcmp[1] already computed **/
    bzcmp[2] = (z >= z11);  bzcmp[3] = (z >= z01);
    if(bcross[0])
    { ncross = 0;
      xcrossp[ncross] = x0;  ycrossp[ncross] = xycross[0];
    }
    else
    { ncross = -1;
    }
    bcross[1] = bzcmp[1] != bzcmp[2];
    if(bcross[1])
    { ncross++;
      deltaz = z11 - z10;     zweight1 = (z - z10)/deltaz;
      xycross[1] = x0 + deltax * zweight1;
      xcrossp[ncross] = xycross[1];  ycrossp[ncross] = y1;
    }
    bcross[2] = bzcmp[2] != bzcmp[3];
    if(bcross[2])
    { ncross++;
      deltaz = z11 - z01;   zweight1 = (z - z01)/deltaz;
      xycross[2] = y0 + deltay * zweight1;
      xcrossp[ncross] = x1;  ycrossp[ncross] = xycross[2];
    }
    bcross[3] = bzcmp[3] != bzcmp[0];
    if(bcross[3])
    { ncross++;
      deltaz = z01 - z00;     zweight1 = (z - z00)/deltaz;
      xycross[3] = x0 + deltax * zweight1;
      xcrossp[ncross] = xycross[3];  ycrossp[ncross] = y0;
    }
    if(ncross < 0) continue;

if((ncross % 2) == 0) fPrintF(stderr,"%s : BIZARRE, ncross = %d\n", prognamp,ncross);
                                                      /** case : 2 intersections. **/
    if(ncross < 2)
    { lcSegVecInc1CoorPData(lcsegp, xcrossp, ycrossp, iz);
      continue;
    }
                        /** case : 4 intersections. Compare to saddle point value **/
    zcol = (z00 - zmean)*(z11 - zmean) - (z01 - zmean)*(z10 - zmean);
    zcol = zmean + zcol/(z00 + z11 - z01 - z10);
    if(z >= zcol && bzcmp[1])
    { lcSegVecInc1CoorPData(lcsegp,xcrossp,ycrossp, iz);
      lcSegVecInc1CoorPData(lcsegp,xcrossp+2,ycrossp+2, iz);
    }
    else
    { lcSegVecInc1CoorPData(lcsegp,xcrossp+1,ycrossp+1, iz);
      lcSegVecInc1CoorData(lcsegp,xcrossp[3],ycrossp[3],xcrossp[0],ycrossp[0], iz);
    }
  }
  return ok;
}
/******************************************************************************/
/*        lcSegAddSFromLineTriangle                                           */
/*                                                                            */
/*                            xi                                              */
/*    y1                      i~2   level zi                                  */
/*                          | |                                               */
/*                        |  |                                                */
/*                     0    1                                                 */
/*                  |      |                                                  */
/*                |       |                                                   */
/*    y0       0----2----1         levels *zp0                                */
/*                                                                            */
/*         x0=xp[ix]  x1=xp[ix+1]                                             */
/*                                                                            */
/*  nx+1 points xp[0],..., xp[nx]  donc nx triangles                          */
/******************************************************************************/
int       lcSegAddSFromLineTriangle(double *xp, size_t nx, double xi,
                                                              double y0, double y1,
                          double *z0p, double zi, double z, int iz, lcSegVec *lcsegp)
{ int       ok = 1;
  int       ix, ncross;
  double    x0, x1, deltaxi0, deltaxi1, deltax10;
  double    deltay;
  double    z0, z1, deltaz, zweight = 0.0;
  static char    form1p[]   = "%s::%s  BIZARRE, ncross = %d\n";
  static char    prognamp[] = "lcSegAddSFromLineTriangle";

  if(nx <= 0) return ok;

/*fPrintF(stdout,"Triangle: y0=%f, y1=%f, xi=%f, zi=%f, nx=%d \n", y0,y1,xi,zi,nx); */
  lcSegPVecRealloc(lcsegp, lcsegp->x + nx, nx/20);
  deltay = y1 - y0;
  bzcmp[2] = (z >= zi);
                                                          /** Loop Initialisation **/
  x1 = *xp++;   deltaxi1 = xi - x1;
  z1 = *z0p++;
  bzcmp[1] = (z >= z1);
  bcross[1] = (bzcmp[1] != bzcmp[2]);
  ncross = -1;
  if(bcross[1])
  { deltaz = zi - z1;     zweight = (z - z1)/deltaz;
    ncross = 0;
    xcrossp[ncross] = x1 + deltaxi1 * zweight;
    ycrossp[ncross] = y0 + deltay * zweight;
  }
                                                                       /** x Loop **/
  for(ix = 1;  ix <= nx;  ix++)
  { x0 = x1;  deltaxi0 = deltaxi1;
    bzcmp[0] = bzcmp[1];
    bcross[0] = bcross[1];
    z0 = z1;
                                       /** bzcmp[0] and bzcmp[2] already computed **/
    x1 = *xp++;  deltaxi1 = xi - x1;  deltax10 = x1 - x0;
    z1 = *z0p++;
    bzcmp[1] = (z >= z1);
    if(bcross[0])
    {                         /** previous one was stored in highest ncross value **/
      xcrossp[0] = xcrossp[ncross];
      ycrossp[0] = ycrossp[ncross];
      ncross = 0;
    }
    else
    { ncross = -1;
    }
    bcross[2] = bzcmp[1] != bzcmp[0];
    if(bcross[2])
    { ncross++;
      deltaz = z1 - z0;     zweight = (z - z0)/deltaz;
      xcrossp[ncross] = x0 + deltax10 * zweight;
      ycrossp[ncross] = y0;
    }
    bcross[1] = (bzcmp[1] != bzcmp[2]);   
    if(bcross[1])                        
    { ncross++;
      deltaz = zi - z1;     zweight = (z - z1)/deltaz;
      xcrossp[ncross] = x1 + deltaxi1 * zweight;
      ycrossp[ncross] = y0 + deltay * zweight;
    }
    if(ncross < 0) continue;

if((ncross % 2) == 0) fPrintF(stderr, form1p, srcfilenamp, prognamp, ncross);
                                                       /** case : 2 intersections **/
    lcSegVecInc1CoorPData(lcsegp,xcrossp,ycrossp, iz);
    continue;
  }
  return ok;
}
/******************************************************************************/
/*                                                                            */
/* Domain :                                                                   */
/*          ny points  yp[0],...,yp[ny-1]                                     */
/*            i.e. iy in [0       , ny-1]                                     */
/*          nx points  xp[0],...,xp[nx-1] restricted to                       */
/*            ix in [ ixiGeop[iy], ixfGeop[iy]]                               */
/*                                                                            */
/*          z = f(x,y) in  zfpp[iy][ix]                                       */
/******************************************************************************/
int       lcSegAddS(double *xp, size_t nx, double *yp, size_t ny, 
                                       myBOOL brectangle, int *ixiGeop, int *ixfGeop,
                               double **zfpp,  double z, int iz, lcSegVec *lcsegvecp)
{ int       ok = 1;
  int       iy;
  int       ixi0, ixi1, ixf0,ixf1;
  int       ixix, ixin, ntix, ixfx, ixfn, ntfx;
  double   *zf0p, *zf1p;
  myBOOL    bi, bf;

  zf1p = *zfpp++ ;
  if(brectangle)
  { for(iy = 0;  iy < ny - 1;  iy++, yp++)
    { zf0p = zf1p;
      zf1p = *zfpp++ ;
      ok = lcSegAddSFromLineRectangle(xp,nx -1, yp, zf0p,zf1p, z,iz, lcsegvecp);
    }
    return ok;
  }
  else
  { for(iy = 0;  iy < ny - 1;  iy++, yp++)
    { zf0p = zf1p;
      zf1p = *zfpp++ ;
      ixi0 = ixiGeop[iy];  ixi1 = ixiGeop[iy+1];
      if((bi = ixi1 >= ixi0)){ ixix = ixi1; ixin = ixi0;}
      else                   { ixix = ixi0; ixin = ixi1;}
      ntix = ixix-ixin;
      ixf0 = ixfGeop[iy];  ixf1 = ixfGeop[iy+1];
      if((bf = ixf1 <= ixf0)){ ixfn = ixf1;  ixfx = ixf0;}
      else                   { ixfn = ixf0;  ixfx = ixf1;}
      ntfx = ixfx-ixfn;

      if(bi)
      { ok = lcSegAddSFromLineTriangle(xp+ixin, ntix, xp[ixix], yp[0], yp[1],
                                             zf0p+ixin, zf1p[ixix], z,iz, lcsegvecp);
      }
      else
      { ok = lcSegAddSFromLineTriangle(xp+ixin, ntix, xp[ixix], yp[1], yp[0],
                                             zf1p+ixin, zf0p[ixix], z,iz, lcsegvecp);
      }
      ok = lcSegAddSFromLineRectangle(xp+ixix, ixfn-ixix, yp,
                                              zf0p+ixix, zf1p+ixix, z,iz, lcsegvecp);
      if(bf)
      { ok = lcSegAddSFromLineTriangle(xp+ixfn, ntfx, xp[ixfn], yp[0], yp[1],
                                             zf0p+ixfn, zf1p[ixfn], z,iz, lcsegvecp);
      }
      else
      { ok = lcSegAddSFromLineTriangle(xp+ixfn, ntfx, xp[ixfn], yp[1], yp[0],
                                             zf1p+ixfn, zf0p[ixfn], z,iz, lcsegvecp);
      }
    }
  }
  return ok;
}
/******************************************************************************/
/*        lcSegAddMFromLineRectangle                                          */
/*                                                                            */
/*    y1       1---1---2      levels *zp1                                     */
/*             |       |                                                      */
/*             |       |                                                      */
/*             0       2                                                      */
/*             |       |                                                      */
/*             |       |                                                      */
/*    y0       0---3---3      levels *zp0                                     */
/*                                                                            */
/*             x0     x1                                                      */
/*                                                                            */
/*  nx   rectangles                                                           */
/*  nx+1 points xp[0], xp[nx]                                                 */
/******************************************************************************/
int       lcSegAddMFromLineRectangle(double *xp, size_t nx, double *yp,
                          double *z0p, double *z1p, lcZrange *lcZp, lcSegVec *lcsegp)
{ int       ok = 1;
  int       ix, ncross;
  double    x0, x1, deltax, y0, y1, deltay;
  int       iz, izn, izx;
  double    zweight = 0.0, zmean, zcol, zsum0, zsum1;
  double    z00, z01, z10, z11, z0n, z0x, z1n, z1x, zcn, zcx;
  double    delta0z, delta1z, delta2z, delta3z;
  double    zo, zn, zx, zpa, z;
  static char    form1p[]   = "%s::%s  BIZARRE, ncross = %d\n";
  static char    prognamp[] = "utiCurve.level::lcSegAddMFromLineRectangle";

  y0 = *yp++;  y1 = *yp;  deltay = y1 - y0;
/* fPrintF(stderr,"Rectangle: y0=%f, y1=%f, deltay = %f\n", y0,y1,deltay); */

  lcSegPVecRealloc(lcsegp, lcsegp->x + nx, nx/20);
  zo = lcZp->zo;  zn = lcZp->zn;  zx = lcZp->zx;  zpa = lcZp->zpa;

                                                          /** Loop Initialisation **/
  x1 = xp[0];
  z01 = z0p[0];  z11 = z1p[0];  delta2z = z11 - z01;  zsum1 = z11 + z01;
  if(z01 > z11){ z1n = z11;  z1x = z01;}
  else         { z1n = z01;  z1x = z11;}
                                                                       /** x Loop **/
  for(ix = 1;  ix <= nx;  ix++)
  { x0 = x1;  x1 = xp[ix];  deltax = x1 - x0;
    z00 = z01;  z10 = z11;  z0n = z1n;  z0x = z1x;
    z01 = z0p[ix];  z11 = z1p[ix];
    delta0z = delta2z;    delta1z = z11 - z10;     /* delta0z = z10 - z00 */
    delta2z = z11 - z01;  delta3z = z01 - z00;
    zsum0 = zsum1;  zsum1 = z11 + z01;
    zmean = (zsum0 + zsum1)/4 ;
    zcol = (z00 - zmean)*(z11 - zmean) - (z01 - zmean)*(z10 - zmean);
    zcol = zmean + zcol/(z00 + z11 - z01 - z10);

    if(z01 > z11){ z1n = z11;  z1x = z01;}
    else         { z1n = z01;  z1x = z11;}
    if(z0x > z1x){ zcx = z0x;}
    else         { zcx = z1x;}
    if(z0n < z1n){ zcn = z0n;}
    else         { zcn = z1n;}
    if(zcx > zx) zcx = zx;
    if(zcn < zn) zcn = zn;
    izn = (zcn - zo)/zpa;  izx = (zcx - zo)/zpa;
                                                                   /** level Loop **/
    for(iz = izn;  iz <= izx;  iz++)
    { z = zo + iz * zpa;
      bzcmp[0] = (z >= z00);  bzcmp[1] = (z >= z10);
      bzcmp[2] = (z >= z11);  bzcmp[3] = (z >= z01);
      bcross[0] = bzcmp[0] != bzcmp[1];
      bcross[1] = bzcmp[1] != bzcmp[2];
      bcross[2] = bzcmp[2] != bzcmp[3];
      bcross[3] = bzcmp[3] != bzcmp[0];
      ncross = -1;
      if(bcross[0])
      { ncross++;  xcrossp[ncross] = x0;
        zweight = (z - z00)/delta0z;  ycrossp[ncross] = y0 + deltay * zweight;
      }
      if(bcross[2])
      { ncross++;  xcrossp[ncross] = x1;
        zweight = (z - z01)/delta2z;  ycrossp[ncross] = y0 + deltay * zweight;
      }
      if(bcross[1])
      { ncross++;  ycrossp[ncross] = y1;
        zweight = (z - z10)/delta1z;  xcrossp[ncross] = x0 + deltax * zweight;
      }
      if(bcross[3])
      { ncross++;  ycrossp[ncross] = y0;
        zweight = (z - z00)/delta3z;  xcrossp[ncross] = x0 + deltax * zweight;
      }
      if(ncross < 0) continue;

if((ncross % 2) == 0) fPrintF(stderr,form1p, srcfilenamp, prognamp,ncross);
                                                       /** case : 2 intersections **/
      if(ncross < 2)
      { lcSegVecInc1CoorPData(lcsegp, xcrossp, ycrossp, iz);
        continue;
      }
                       /** case : 4 intersections. Compare to saddle point value **/
      if(z >= zcol && bzcmp[1])
      { lcSegVecInc1CoorPData(lcsegp,xcrossp,ycrossp, iz);
        lcSegVecInc1CoorPData(lcsegp,xcrossp+2,ycrossp+2, iz);
      }
      else
      { lcSegVecInc1CoorPData(lcsegp,xcrossp+1,ycrossp+1, iz);
        lcSegVecInc1CoorData(lcsegp,xcrossp[3],ycrossp[3],xcrossp[0],ycrossp[0], iz);
      }
    }
  }
  return ok;
}
/******************************************************************************/
/*        lcSegAddMFromLineTriangle                                           */
/*                                                                            */
/*                            xi                                              */
/*    y1                      i~2   level zi                                  */
/*                          | |                                               */
/*                        |  |                                                */
/*                     0    1                                                 */
/*                  |      |                                                  */
/*                |       |                                                   */
/*    y0       0----2----1         levels *zp0                                */
/*                                                                            */
/*         x0=xp[ix]  x1=xp[ix+1]                                             */
/*                                                                            */
/*  nx   triangles                                                            */
/*  nx+1 points xp[0],..., xp[nx]                                             */
/******************************************************************************/
int       lcSegAddMFromLineTriangle(double *xp, size_t nx, double xi,
                                                              double y0, double y1,
                            double *z0p, double zi, lcZrange *lcZp, lcSegVec *lcsegp)
{ int       ok = 1;
  int       ix, ncross, izn, izx, iz;
  double    x1, x0, deltax10, deltaxi0, deltaxi1;
  double    deltay;
  double    z0, z1, deltaz, zweight = 0.0;
  double    z0n, z0x, z1n, z1x, zcn, zcx;
  double    zo, zn, zx, zpa, z;
  static char    form1p[]   = "%s::%s  BIZARRE, ncross = %d\n";
  static char    prognamp[] = "lcSegAddMFromLineTriangle";

  if(nx <= 0) return ok;

  deltay = y1 - y0;
/*fPrintF(stdout,"Triangle: y0=%f, y1=%f, xi=%f, zi=%f, nx=%d \n", y0,y1,xi,zi,nx); */

  lcSegPVecRealloc(lcsegp, lcsegp->x + nx, nx/20);
  zo = lcZp->zo;  zn = lcZp->zn;  zx = lcZp->zx;  zpa = lcZp->zpa;

                                                          /** Loop Initialisation **/
  x1 = *xp++;  deltaxi1 = xi - x1;
  z1 = *z0p++;
  if(zi > z1){ z1x = zi;  z1n = z1;}
  else       { z1x = z1;  z1n = zi;}
                                                                       /** x Loop **/
  for(ix = 1;  ix <= nx;  ix++)
  { x0 = x1;      deltaxi0 = deltaxi1;
    x1 = *xp++;   deltaxi1 = xi - x1;  deltax10 = x1 - x0;
    z0 = z1;     z0x = z1x;  z0n = z1n;
    z1 = *z0p++;
    if(zi > z1){ z1x = zi;  z1n = z1;}
    else       { z1x = z1;  z1n = zi;}
    if(z0x > z1x){ zcx = z0x;}
    else         { zcx = z1x;}
    if(z0n < z1n){ zcn = z0n;}
    else         { zcn = z1n;}
    if(zcx > zx) zcx = zx;
    if(zcn < zn) zcn = zn;
    izn = (zcn - zo)/zpa;  izx = (zcx - zo)/zpa;
                                                                   /** level Loop **/
    for(iz = izn;  iz <= izx;  iz++)
    { z = zo + iz * zpa;
      bzcmp[0] = z >= z0;
      bzcmp[1] = z >= z1;
      bzcmp[2] = z >= zi;
      bcross[0] = (bzcmp[0] != bzcmp[2]);
      bcross[1] = (bzcmp[1] != bzcmp[2]);
      bcross[2] = (bzcmp[0] != bzcmp[1]);
      ncross = -1;
      if(bcross[0])
      { ncross++;
        deltaz = zi - z0;     zweight = (z - z0)/deltaz;
        xcrossp[ncross] = x0 + deltaxi0 * zweight;
        ycrossp[ncross] = y0 + deltay * zweight;
      }
      if(bcross[1])
      { ncross++;
        deltaz = zi - z1;     zweight = (z - z1)/deltaz;
        xcrossp[ncross] = x1 + deltaxi1 * zweight;
        ycrossp[ncross] = y0 + deltay * zweight;
      }
      if(bcross[2])
      { ncross++;
        deltaz = z1 - z0;     zweight = (z - z0)/deltaz;
        ycrossp[ncross] = y0;
        xcrossp[ncross] = x0 + deltax10 * zweight;
      }
      if(ncross < 0) continue;

if((ncross % 2) == 0) fPrintF(stderr, form1p, srcfilenamp, prognamp, ncross);
                                                       /** case : 2 intersections **/
      lcSegVecInc1CoorPData(lcsegp,xcrossp,ycrossp, iz);
    }
  }
  return ok;
}
/******************************************************************************/
/*                                                                            */
/* Domain :                                                                   */
/*          ny points  yp[0],...,yp[ny-1]                                     */
/*            i.e. iy in [0       , ny-1]                                     */
/*          nx points  xp[0],...,xp[nx-1] restricted to                       */
/*            ix in [ ixiGeop[iy], ixfGeop[iy]]                               */
/*                                                                            */
/*          z = f(x,y) in  zfpp[iy][ix]                                       */
/******************************************************************************/
int       lcSegAddM(double *xp, size_t nx, double *yp, size_t ny, 
                                       myBOOL brectangle, int *ixiGeop, int *ixfGeop,
                                  double **zfpp, lcZrange *lcZp, lcSegVec *lcsegvecp)
{ int       ok = 1;
  int       iy;
  int       ixi0, ixi1, ixf0,ixf1;
  int       ixix, ixin, ntix, ixfx, ixfn, ntfx;
  double   *zf0p, *zf1p;
  myBOOL    bi, bf;

  zf1p = *zfpp++ ;
  if(brectangle)
  { for(iy = 0;  iy < ny - 1;  iy++, yp++)
    { zf0p = zf1p;
      zf1p = *zfpp++ ;
      ok = lcSegAddMFromLineRectangle(xp,nx -1, yp, zf0p,zf1p, lcZp, lcsegvecp);
    }
    return ok;
  }
  else
  { for(iy = 0;  iy < ny - 1;  iy++, yp++)
    { zf0p = zf1p;
      zf1p = *zfpp++ ;
      ixi0 = ixiGeop[iy];  ixi1 = ixiGeop[iy+1];
      if((bi = ixi1 >= ixi0)){ ixix = ixi1; ixin = ixi0;}
      else                   { ixix = ixi0; ixin = ixi1;}
      ntix = ixix-ixin;
      ixf0 = ixfGeop[iy];  ixf1 = ixfGeop[iy+1];
      if((bf = ixf1 <= ixf0)){ ixfn = ixf1;  ixfx = ixf0;}
      else                   { ixfn = ixf0;  ixfx = ixf1;}
      ntfx = ixfx-ixfn;

      if(bi)
      { ok = lcSegAddMFromLineTriangle(xp+ixin, ntix, xp[ixix], yp[0], yp[1],
                                             zf0p+ixin, zf1p[ixix], lcZp, lcsegvecp);
      }
      else
      { ok = lcSegAddMFromLineTriangle(xp+ixin, ntix, xp[ixix], yp[1], yp[0],
                                             zf1p+ixin, zf0p[ixix], lcZp, lcsegvecp);
      }
      ok = lcSegAddMFromLineRectangle(xp+ixix, ixfn-ixix, yp,
                                              zf0p+ixix, zf1p+ixix, lcZp, lcsegvecp);
      if(bf)
      { ok = lcSegAddMFromLineTriangle(xp+ixfn, ntfx, xp[ixfn], yp[0], yp[1],
                                             zf0p+ixfn, zf1p[ixfn], lcZp, lcsegvecp);
      }
      else
      { ok = lcSegAddMFromLineTriangle(xp+ixfn, ntfx, xp[ixfn], yp[1], yp[0],
                                             zf1p+ixfn, zf0p[ixfn], lcZp, lcsegvecp);
      }
    }
  }
  return ok;
}

/******************************************************************************/
/*                                                                            */
/* Given 2 sets of n points {P00, P01,...,P0n}, {P10, P11,...,P1n}            */
/*    where each point is a (double valued) triplet (x0,x1,x2)~(x,y,z)        */
/*                                                                            */
/* Consider the quadrangles (P00,P01,P11,P10), (P01,P02,P12,P11),             */
/*                          (P02,P03,P13,P12), ...                            */
/* with possibly P.n+1 = P.0 if closed                                        */
/*                                                                            */
/* Compute the 2 coordinates index lp[0], lp[1],                              */
/*    for level h, index ih, of coordinate lp[2]                              */
/******************************************************************************/
int       lcSegAddSFromQuadrangles(double *x0p, double *x1p, size_t nx, short close,
                                        int *lp, double h, int ih, lcSegVec *lcsegvp)
{ int       ok = 1, lx, ly, lz;
  int       i, ix, ncross;
  double   *x0cp, *x1cp;
  double    x00, x01, x10, x11, y00, y01, y10, y11, z00, z01, z10, z11;
  double    deltaz, zweight1 = 0.0, zmean, zcol, zsum0, zsum1;
  double    x1inter = 0.0,  y1inter = 0.0;

  if(nx <= 0) return ok;
  lcSegPVecRealloc(lcsegvp, lcsegvp->x + nx/100 +10, nx/200);

  lx = lp[0];  ly = lp[1];  lz = lp[2];

  x0cp = x0p;  x1cp = x1p;
  x01 = *(x0cp +lx);  x11 = *(x1cp +lx);
  y01 = *(x0cp +ly);  y11 = *(x1cp +ly);
  z01 = *(x0cp +lz);  z11 = *(x1cp +lz);
  zsum1 = z11 + z01;
  bzcmp[2] = (h >= z11);  bzcmp[3] = (h >= z01);
  bcross[2] = (bzcmp[2] != bzcmp[3]);
  if(bcross[2])
  { deltaz = z11 - z01;   zweight1 = (h - z01)/deltaz;
    x1inter = x01 + (x11 - x01) * zweight1;
    y1inter = y01 + (y11 - y01) * zweight1;
  }
                                                                       /** i Loop **/
  ix = (close)? nx -1: nx -2;
  for(i = 0;  i <= ix;  i++)
  { x00 = x01;   x10 = x11;
    y00 = y01;   y10 = y11;
    z00 = z01;   z10 = z11;  zsum0 = zsum1;
    bzcmp[0] = bzcmp[3];
    bzcmp[1] = bzcmp[2];
    bcross[0] = bcross[2];
    xcrossp[0] = xcrossp[2];
    ycrossp[0] = ycrossp[2];
    x0cp += 3;   x1cp += 3;
    if(i == nx -1){ x0cp = x0p;  x1cp = x1p;}

    x01 = *(x0cp +lx);  x11 = *(x1cp +lx);
    y01 = *(x0cp +ly);  y11 = *(x1cp +ly);
    z01 = *(x0cp +lz);  z11 = *(x1cp +lz);
    zsum1 = z11 + z01;
    zmean = (zsum0 + zsum1)/4 ;
                                       /** bzcmp[0] and bzcmp[1] already computed **/
    bzcmp[2] = (h >= z11);  bzcmp[3] = (h >= z01);
    if(bcross[0])
    { ncross = 0;
      xcrossp[ncross] = x1inter;  ycrossp[ncross] = y1inter;
    }
    else
    { ncross = -1;
    }
    bcross[1] = bzcmp[1] != bzcmp[2];
    if(bcross[1])
    { ncross++;
      deltaz = z11 - z10;     zweight1 = (h - z10)/deltaz;
      xcrossp[ncross] = x10 + (x11 - x10) * zweight1;
      ycrossp[ncross] = y10 + (y11 - y10) * zweight1;
    }
    bcross[2] = bzcmp[2] != bzcmp[3];
    if(bcross[2])
    { ncross++;
      deltaz = z11 - z01;   zweight1 = (h - z01)/deltaz;
      x1inter = x01 + (x11 - x01) * zweight1;
      y1inter = y01 + (y11 - y01) * zweight1;
      xcrossp[ncross] = x1inter;
      ycrossp[ncross] = y1inter;
    }
    bcross[3] = bzcmp[3] != bzcmp[0];
    if(bcross[3])
    { ncross++;
      deltaz = z01 - z00;     zweight1 = (h - z00)/deltaz;
      xcrossp[ncross] = x00 + (x01 - x00) * zweight1;
      ycrossp[ncross] = y00 + (y01 - y00) * zweight1;
    }
    if(ncross < 0) continue;

                                                      /** case : 2 intersections. **/
    if(ncross < 2)
    { lcSegVecInc1CoorPData(lcsegvp, xcrossp, ycrossp, ih);
      continue;
    }
                        /** case : 4 intersections. Compare to saddle point value **/
    zcol = (z00 - zmean)*(z11 - zmean) - (z01 - zmean)*(z10 - zmean);
    zcol = zmean + zcol/(z00 + z11 - z01 - z10);
    if(h >= zcol && bzcmp[1])
    { lcSegVecInc1CoorPData(lcsegvp, xcrossp, ycrossp, ih);
      lcSegVecInc1CoorPData(lcsegvp, xcrossp+2, ycrossp+2, ih);
    }
    else
    { lcSegVecInc1CoorPData(lcsegvp, xcrossp+1, ycrossp+1, ih);
      lcSegVecInc1CoorData(lcsegvp, xcrossp[3],ycrossp[3],xcrossp[0],ycrossp[0], ih);
    }
  }
  return ok;
}
/******************************************************************************/
/******************************************************************************/
