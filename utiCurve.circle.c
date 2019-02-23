/*  essai/C/utiCurve.circle.c                                                 */
/*  Mennessier Gerard                   20010507                              */
/*  Last revised M.G.                   20021027                              */
/*                                                                            */

#include  <stddef.h>
#include  <math.h>                                       /** for function  fabs() **/
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

#include  "utiCurve.circle.h"
#include  "utiCurve.GC.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cCircle  *cCircleAlloc(size_t  nz,char *progcallp)
{ cCircle     *p;
  static    char      prognamp[] = "utiCurve.circle::cCircleAlloc";
  static    char      form1p[] = "called from %s. malloc failed; cCircle size nz=%d\n";

  p = (cCircle *)malloc( nz*sizeof(cCircle) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  return p;
}
/******************************************************************************/
/*  cCircle *cCircleChkRealloc(cCircle,nzp,neednz,incrnz,progcallp)           */
/*                                                                            */
/* check if there is enough memory allocated for the neednz needed cCircle    */
/* old and possibly new size  stored as  *nzp                                 */
/* if reallocation is needed, newsize =  neednz + incrnz                      */
/******************************************************************************/
cCircle  *cCircleChkRealloc(cCircle *pi, size_t *nzp, size_t neednz, size_t incrnz,
                                                                     char *progcallp)
{ cCircle    *p;
  size_t    nfz;
  static    char      prognamp[] = "utiCurve.circle::cCircleChkRealloc";
  static    char      form1p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired cCircle size nz=%d\n";
  static    char      form2p[]= "called from %s. malloc failed;"
                                                       " desired cCircle size nz=%d\n";
  
  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (cCircle *)malloc( nfz*sizeof(cCircle) );
fPrintF(stderr,"%s called from %s; realloc p=%p\n", prognamp, progcallp, (void*)pi);
    if(p == NULL) myErr1(-1,stderr,prognamp, form2p, progcallp, neednz);
  }
  else
  { if(neednz  <=  *nzp) return pi;
    p = (cCircle *)realloc(pi, nfz*sizeof(cCircle) );
    if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,pi,neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*  void  cCirclePrint(FILE *bufp,cCircle *p)                                 */
/******************************************************************************/
void      cCirclePrint(FILE  *bufp, cCircle *p)
{ static    char form1[] = "%s  center={%f, %f}, radius=%f, angles={%f, %f},"
                                                                      " closed=%d\n";
  static    char prognamp[] = "utiCurve.circle::cCirclePrint";

  fPrintF(bufp,form1, prognamp, p->xp[0],p->xp[1], p->r, p->anglp[0],p->anglp[1],
                                                                          p->closed);
  return;
}
/******************************************************************************/
/*  void  cCircleArrayPrint(FILE *bufp,cSeg *p,size_t n)                      */
/******************************************************************************/
void      cCircleArrayPrint(FILE *bufp,cCircle *p,size_t n)
{ n++;
  while(--n){ cCirclePrint(bufp, p);  p++;}
  return;
}

/******************************************************************************/
/*  void  cCircleZero(cCircle *p)                                             */
/******************************************************************************/
void      cCircleZero(cCircle *p)
{ double   *dp;

  dp = p->xp;       dp[0] = dp[1] = 0.0;
  p->r = 0.0;
  dp = p->anglp;    dp[0] = dp[1] = 0.0;
  p->gcp = NULL;
  p->closed = 0;
  return;
}
/******************************************************************************/
/*        cCircleSet(p,xp,r,anglp)                                            */
/******************************************************************************/
void      cCircleSet(cCircle *p, double *xp, double r, double *anglp)
{ double   *dp, diffangl;

  dp = p->xp;  *dp++ = xp[0];  *dp = xp[1];
  p->r = r;
  diffangl = fabs(anglp[1] - anglp[0]) + 0.0000001;
  dp = p->anglp;  *dp++ = anglp[0];  *dp = anglp[1];
  if( (int)(diffangl) % 360 == 0.0) p->closed = 1;
  else                              p->closed = 0;
  return ;
}
/******************************************************************************/
/*        cCircleCopypp(p,xp,r,anglp)                                         */
/******************************************************************************/
void      cCircleCopypp(cCircle *pf, cCircle *pi)
{ double   *dfp, *dip;
  
  dip = pi->xp;  dfp = pf->xp;  *dfp++ = *dip++;  *dfp = *dip;
  pf->r = pi->r;
  dip = pi->anglp;  dfp = pf->anglp;  *dfp++ = *dip++;  *dfp = *dip;
  pf->closed = pi->closed;
  pf->gcp = pi->gcp;
  return;
}
/******************************************************************************/
/*        cCircleCopyp(p,xp,r,anglp)                                          */
/******************************************************************************/
void      cCircleCopyp(cCircle *pf, cCircle cc)
{ double   *dfp, *dip;
  
  dip = cc.xp;  dfp = pf->xp;  *dfp++ = *dip++;  *dfp = *dip;
  pf->r = cc.r;
  dip = cc.anglp;  dfp = pf->anglp;  *dfp++ = *dip++;  *dfp = *dip;
  pf->closed = cc.closed;
  pf->gcp = cc.gcp;
  return;
}
/******************************************************************************/


/******************************************************************************/
/*        cCircleVecAlloc(nz,progcallp)                                       */
/******************************************************************************/

/******************************************************************************/
/*        cCirclePVecAlloc(vecp,nz)                                           */
/******************************************************************************/
void      cCirclePVecAlloc(cCircleVec *vecp,size_t  nz)
{ 
  vecp->p = cCircleAlloc(nz,"cCirclePVecAlloc");
  vecp->z = nz;  vecp->x = 0;
  cCircleZero(vecp->p);
  return ;
}
/******************************************************************************/
/*        cCirclePVecRealloc(vecp,neednz,incrnz)                              */
/******************************************************************************/
void      cCirclePVecRealloc(cCircleVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = cCircleChkRealloc(vecp->p, &(vecp->z),neednz,incrnz,"cCirclePVecRealloc");
  return ;
}
/******************************************************************************/
/*        cCirclePVecFree(vecp)                                               */
/******************************************************************************/
void      cCirclePVecFree(cCircleVec *vecp)
{ 
  free(vecp->p);  vecp->p = NULL;  vecp->z = 0;  vecp->x = 0;  return ;
}
/******************************************************************************/
/*        cCircleVecFree(vecp)                                                */
/******************************************************************************/
void      cCircleVecFree(cCircleVec *vecp)
{ 
  if(vecp->p != NULL) cCirclePVecFree(vecp);
  free(vecp);  return ;
}
/******************************************************************************/
/*        cCircleVecPrint(vecp)                                               */
/******************************************************************************/
void      cCircleVecPrint(FILE  *bufp, cCircleVec *vecp)
{ cCircle  *p;
  static    char    form1[] = "%s  cCircleVec.p=%p, vec.z=%d, vec.x=%d \n" ;
  static    char    prognamp[] = "utiCurve.circle::cCircleVecPrint";

  p = vecp->p;
  fPrintF(bufp,form1, prognamp, p,vecp->z,vecp->x);
  if(p == NULL) return;
  cCircleArrayPrint(bufp,p,vecp->x);
  fPrintF(bufp, "\n");  return ;
}
/******************************************************************************/
/*        cCircleVecInc1(vecp,y)                                              */
/* add 1 cCircle into a cCircleVec structure                                  */
/******************************************************************************/
void      cCircleVecInc1(cCircleVec *vecp, cCircle cc)
{ cCircle  *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cCirclePVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  cCircleCopyp(p, cc);
  return ;
}
/******************************************************************************/
/*        cCircleVecInc1p(vecp,yp)                                            */
/* add 1 *cCircle into a cCircleVec structure                                 */
/******************************************************************************/
void      cCircleVecInc1p(cCircleVec *vecp, cCircle *ccp)
{ cCircle  *p;
  size_t    newz,oldx;

  oldx = vecp->x ; newz = oldx+1;
  cCirclePVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  cCircleCopypp(p, ccp);
  return;
}
/******************************************************************************/
/*        cCircleVecInc1Data(vecp,y)                                          */
/* add 1 cCircle into a cCircleVec structure                                  */
/*               from center, radius, angles                                  */
/******************************************************************************/
void      cCircleVecInc1Data(cCircleVec *vecp, double *xp, double r, double *anglp)
{ cCircle  *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cCirclePVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  cCircleSet(p, xp, r, anglp);
  p->gcp = NULL;
  return;
}
/******************************************************************************/
/******************************************************************************/
