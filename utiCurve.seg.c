/*  essai/C/utiCurve.seg.c                                                    */
/*  Mennessier Gerard                   20010507                              */
/*  Last revised M.G.                   20020904                              */
/*                                                                            */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

#include  "utiCurve.seg.h"
#include  "utiCurve.GC.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cSeg     *cSegAlloc(size_t  nz,char *progcallp)
{ cSeg     *p;
  static    char      prognamp[] = "utiCurve.seg::cSegAlloc";
  static    char      form1p[] = "called from %s. malloc failed; cSeg size nz=%d\n";

  p = (cSeg *)malloc( nz*sizeof(cSeg) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  return p;
}
/******************************************************************************/
/*  cSeg *cSegChkRealloc(cSegp,nzp,neednz,incrnz,progcallp)                   */
/*                                                                            */
/* check if there is enough memory allocated for the neednz needed cSeg       */
/* old and possibly new size  stored as  *nzp                                 */
/* if reallocation is needed, newsize =  neednz + incrnz                      */
/******************************************************************************/
cSeg     *cSegChkRealloc(cSeg *pi,size_t  *nzp,size_t neednz,size_t incrnz,
                                                                      char *progcallp)
{ cSeg     *p;
  size_t    nfz;
  static    char      prognamp[] = "utiCurve.seg::cSegChkRealloc";
  static    char      form1p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired cSeg size nz=%d\n";
  static    char      form2p[]= "called from %s. malloc failed;"
                                                       " desired cSeg size nz=%d\n";
  
  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (cSeg *)malloc( nfz*sizeof(cSeg) );
fPrintF(stderr,"%s called from %s; realloc p=%p\n", prognamp, progcallp, (void*)pi);
    if(p == NULL) myErr1(-1,stderr,prognamp, form2p,progcallp,neednz);
  }
  else
  { if(neednz  <=  *nzp) return pi;
    p = (cSeg *)realloc(pi, nfz*sizeof(cSeg) );
    if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,pi,neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*  void  cSegPrint(FILE *bufp,cSeg *p)                                       */
/******************************************************************************/
void      cSegPrint(FILE *bufp,cSeg *p)
{ 
  fPrintF(bufp,"P1=(%f, %f), P2=(%f, %f)\n",  p->x[0], p->y[0], p->x[1], p->y[1]);
  return;
}
/******************************************************************************/
/*  void  cSegArrayPrint(FILE *bufp,cSeg *p,size_t n)                         */
/******************************************************************************/
void      cSegArrayPrint(FILE *bufp,cSeg *p,size_t n)
{ n++;
  while(--n)
  { fPrintF(bufp,"P1=(%f, %f), P2=(%f, %f)\n",  p->x[0], p->y[0], p->x[1], p->y[1]);
    p++;
  }
  return;
}
/******************************************************************************/
/*  void  cSegZero(cSeg *p)                                                   */
/******************************************************************************/
void      cSegZero(cSeg *p)
{ 
  p->x[0] = 0.; p->x[1] = 0.;  p->y[0] = 0.; p->y[1] = 0.;  return;
}
/******************************************************************************/
/*  void  cSegCpy1p(cSeg *sfp,cSeg *sip)                                      */
/*                                                                            */
/* copy 1 cSeg  *sip into  *sfp                                               */
/******************************************************************************/
void      cSegCpy1p(cSeg *sfp,cSeg *sip)
{ double    *dip, *dfp;

  dip = sip->x;  dfp = sfp->x;  *dfp++ = *dip++;  *dfp = *dip;
  dip = sip->y;  dfp = sfp->y;  *dfp++ = *dip++;  *dfp = *dip;
  return;
}
/******************************************************************************/


/******************************************************************************/
/*        cSegVecAlloc(nz,progcallp)                                          */
/******************************************************************************/

/******************************************************************************/
/*        cSegPVecAlloc(vecp,nz)                                              */
/******************************************************************************/
void      cSegPVecAlloc(cSegVec *vecp,size_t  nz)
{ 
  vecp->p = cSegAlloc(nz,"cSegPVecAlloc");
  vecp->z = nz;  vecp->x = 0;
  cSegZero(vecp->p);
  return ;
}
/******************************************************************************/
/*        cSegPVecRealloc(vecp,neednz,incrnz)                                 */
/******************************************************************************/
void      cSegPVecRealloc(cSegVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = cSegChkRealloc(vecp->p, &(vecp->z),neednz,incrnz,"cSegPVecRealloc");
  return ;
}
/******************************************************************************/
/*        cSegPVecFree(vecp)                                                  */
/******************************************************************************/
void      cSegPVecFree(cSegVec *vecp)
{ 
  free(vecp->p);  vecp->p = NULL;  vecp->z = 0;  vecp->x = 0;  return ;
}
/******************************************************************************/
/*        cSegVecFree(vecp)                                                   */
/******************************************************************************/
void      cSegVecFree(cSegVec *vecp)
{ 
  if(vecp->p != NULL) cSegPVecFree(vecp);
  free(vecp);  return ;
}
/******************************************************************************/
/*        cSegVecPrint(vecp)                                                  */
/******************************************************************************/
void      cSegVecPrint(FILE  *bufp, cSegVec *vecp)
{ cSeg     *p;
  static    char    form1[] = "%s  cSegVec.p=%p, vec.z=%d, vec.x=%d \n" ;
  static    char    prognamp[] = "uitCurve.seg::cSegVecPrint";

  p = vecp->p;
  fPrintF(bufp,form1, prognamp, p,vecp->z,vecp->x);
  if(p == NULL) return;
  cSegArrayPrint(bufp,p,vecp->x);
  fPrintF(bufp, "\n");  return ;
}
/******************************************************************************/
/*        cSegVecInc1(vecp,y)                                                 */
/* add 1 cSeg into a cSegVec structure                                        */
/******************************************************************************/
void      cSegVecInc1(cSegVec *vecp,cSeg s)
{ cSeg     *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  p->x[0] = s.x[0];  p->x[1] = s.x[1];  p->y[0] = s.y[0];  p->y[1] = s.y[1];
  return ;
}
/******************************************************************************/
/*        cSegVecInc1p(vecp,yp)                                               */
/* add 1 *cSeg into a cSegVec structure                                       */
/******************************************************************************/
void      cSegVecInc1p(cSegVec *vecp,cSeg *sp)
{ cSeg     *p;
  size_t    newz,oldx;

  oldx = vecp->x ; newz = oldx+1;
  cSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;  cSegCpy1p(p,sp);  return ;
}
/******************************************************************************/
/*        cSegVecInc1PointData(vecp,*p1,*p2)                                  */
/* add 1 segment, point p1 to p2                  into a cSegVec structure    */
/* *p1 = (x1,y1) , *p2 = (x2,y2)                                              */
/******************************************************************************/
void      cSegVecInc1PointData(cSegVec *vecp,double *p1p, double *p2p)
{ cSeg     *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  p->x[0] = *p1p++;  p->x[1] = *p2p++;  p->y[0] = *p1p;  p->y[1] = *p2p;
  return ;
}
/******************************************************************************/
/*        cSegVecInc1CoorData(vecp,x0,y0,x1,y1,iz)                            */
/* add 1 segment, point (x0,y0) to (x1,y1),                                   */
/*                                                into a cSegVec structure    */
/******************************************************************************/
void      cSegVecInc1CoorData(cSegVec *vecp,double x0, double y0,
                                                                double x1, double y1)
{ cSeg    *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  p->x[0] = x0;  p->x[1] = x1;  p->y[0] = y0;  p->y[1] = y1;
  return ;
}
/******************************************************************************/
/*        cSegVecInc1CoorPData(vecp,xp,yp,iz)                                 */
/* add 1 segment, point (xp[0],yp[0]) to (xp[1],yp[1]),                       */
/*                                                into a cSegVec structure    */
/******************************************************************************/
void      cSegVecInc1CoorPData(cSegVec *vecp,double *xp, double *yp)
{ cSeg    *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cSegPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  p->x[0] = *xp++;  p->x[1] = *xp;  p->y[0] = *yp++;  p->y[1] = *yp;
  return ;
}
/******************************************************************************/
/******************************************************************************/
