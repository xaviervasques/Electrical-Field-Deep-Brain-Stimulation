/*  essai/C/utiCurve.poly.c                                                   */
/*  Mennessier Gerard                   20010507                              */
/*  Last revised M.G.                   20020904                              */
/*                                                                            */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

#include  "utiCurve.poly.h"
#include  "utiCurve.GC.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cPoly     *cPolyAlloc(size_t  nz,char *progcallp)
{ cPoly     *p;
  static    char      prognamp[] = "utiCurve.poly::cPolyAlloc";
  static    char      form1p[] = "called from %s. malloc failed; cPoly size nz=%d\n";

  p = (cPoly *)malloc( nz*sizeof(cPoly) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  return p;
}
/******************************************************************************/
/*  cPoly *cPolyChkRealloc(cPoly,nzp,neednz,incrnz,progcallp)                 */
/*                                                                            */
/* check if there is enough memory allocated for the neednz needed cPoly      */
/* old and possibly new size  stored as  *nzp                                 */
/* if reallocation is needed, newsize =  neednz + incrnz                      */
/******************************************************************************/
cPoly     *cPolyChkRealloc(cPoly *pi, size_t *nzp, size_t neednz, size_t incrnz,
                                                                     char *progcallp)
{ cPoly    *p;
  size_t    nfz;
  static    char      prognamp[] = "utiCurve.poly::cPolyChkRealloc";
  static    char      form1p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired cPoly size nz=%d\n";
  static    char      form2p[]= "called from %s. malloc failed;"
                                                       " desired cPoly size nz=%d\n";
  
  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (cPoly *)malloc( nfz*sizeof(cPoly) );
fPrintF(stderr,"%s called from %s; realloc p=%p\n", prognamp, progcallp, (void*)pi);
    if(p == NULL) myErr1(-1,stderr,prognamp, form2p,progcallp,neednz);
  }
  else
  { if(neednz  <=  *nzp) return pi;
    p = (cPoly *)realloc(pi, nfz*sizeof(cPoly) );
    if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,pi,neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*  void  cPolyPrint(FILE *bufp,cPoly *p)                                     */
/******************************************************************************/
void      cPolyPrint(FILE *bufp,cPoly *p)
{ static    char      prognamp[] = "utiCurve.poly::cPolyPrint";

  fPrintF(bufp,"%s  size(alloc, used)=(%d,%d), pointer(x,y)=(%p,%p), closed=%d \n",
         prognamp, (int)p->z, (int)p->x, (void*)p->xp, (void*)p->yp, (int)p->closed);
  return;
}
/******************************************************************************/
/*  void  cPolyZero(cPoly *p)                                                 */
/******************************************************************************/
void      cPolyZero(cPoly *p)
{ 
  p->z = 0; p->x = 0;  p->xp = NULL; p->yp = NULL;  p->closed = 0;
  return;
}
/******************************************************************************/
/*        cPolyCopyp(pf, cc)                                                  */
/*                                                                            */
/******************************************************************************/
void      cPolyCopyp(cPoly *pf, cPoly cc)
{ 
  pf->z = cc.z;  pf->x = cc.x;
  pf->xp = cc.xp;  pf->yp = cc.yp;
  pf->gcp = cc.gcp;
  pf->closed = cc.closed;
  return;
}
/******************************************************************************/
/*        cPolyCopypp(pf, pi)                                                 */
/*                                                                            */
/******************************************************************************/
void      cPolyCopypp(cPoly *pf, cPoly *pi)
{  
  pf->z = pi->z;  pf->x = pi->x;
  pf->xp = pi->xp;  pf->yp = pi->yp;
  pf->gcp = pi->gcp;
  pf->closed = pi->closed;
  return;
}
/******************************************************************************/
/******************************************************************************/

/******************************************************************************/
/*        cPolyPAlloc(p,nz)                                                   */
/******************************************************************************/
void      cPolyPAlloc(cPoly *p,size_t  nz)
{ 
  p->xp = dblAlloc(nz,"cPolyPAlloc");
  p->yp = dblAlloc(nz,"cPolyPAlloc");  
  p->z = nz;  p->x = 0;
  p->closed = 0;
  return ;
}
/******************************************************************************/
/*        cPolyPRealloc(p, neednz, incrnz)                                    */
/******************************************************************************/
void      cPolyPRealloc(cPoly *p,size_t neednz,size_t incrnz)
{ size_t    z;

  z = p->z;
  p->xp = dblChkRealloc(p->xp,      &z,neednz,incrnz,"cPolyPRealloc");
  p->yp = dblChkRealloc(p->yp, &(p->z),neednz,incrnz,"cPolyPRealloc");
  return ;
}
/******************************************************************************/
/*        cPolyPFree(p)                                                       */
/******************************************************************************/
void      cPolyPFree(cPoly *p)
{ 
  if(p->xp != NULL)free(p->xp);  if(p->yp != NULL)free(p->yp);
  p->xp = NULL;  p->yp = NULL;
  p->z = 0;  p->x = 0;
  return ;
}
/******************************************************************************/
/*        cPolyPPrint(p)                                                      */
/******************************************************************************/
void      cPolyPPrint(FILE  *bufp, cPoly *p)
{ int       i;
  double   *xp, *yp;
  static    char    form1[] = "%s  xp=%p, yp=%p, z=%d, x=%d, closed=%d\n";
  static    char    prognamp[] = "utiCurve.poly::cPolyPPrint";

  fPrintF(bufp,form1, prognamp, p->xp,p->yp, p->z,p->x, p->closed);
  xp = p->xp;  yp = p->yp;
  if(xp == NULL) return;
  for(i = 0;  i < p->x;  i++, xp++,yp++){ fPrintF(bufp,"(%f, %f) ", *xp, *yp);}
  fPrintF(bufp, "\n");  return ;
}
/******************************************************************************/
/*        cPolyPInc1(p,x,y)                                                   */
/* add 1 point (x,y) into a cPoly structure                                   */
/******************************************************************************/
void      cPolyPInc1(cPoly *p, double x, double y)
{ size_t    newz, oldx;
  double   *dp;

  oldx = p->x;  newz = oldx +1;
  cPolyPRealloc(p,newz, newz/8);
  dp = p->xp; *(dp+oldx) = x;
  dp = p->yp; *(dp+oldx) = y;
  (p->x)++;
  return ;
}
/******************************************************************************/
/*        cPolyPInc1p(p,xy)                                                   */
/* add 1 point (xy[0],xy[1]) into a cPoly structure                           */
/******************************************************************************/
void      cPolyPInc1p(cPoly *p, double xy[2])
{ size_t    newz, oldx;
  double   *dp;

  oldx = p->x;  newz = oldx +1;
  cPolyPRealloc(p, newz, newz/8);
  dp = p->xp + oldx; *dp = xy[0];  dp = p->yp + oldx; *dp = xy[1];
  (p->x)++ ;
  return ;
}
/******************************************************************************/
/*        cPolyPIncN(p,xp,yp,n)                                               */
/* add N points into a cPoly structure                                        */
/******************************************************************************/
void      cPolyPIncN(cPoly *p, double *xp, double *yp, size_t n)
{ size_t    newz, oldx;
  double   *xfp, *yfp;
  int       i;

  oldx = p->x;  newz = oldx + n;
  cPolyPRealloc(p,newz, newz/8);
  xfp = p->xp;  xfp += oldx;  yfp = p->yp;  yfp += oldx;
  for(i = 0;  i < n;  i++){ *xfp++ = *xp++;  *yfp++ = *yp++;}
  p->x += n;
  return ;
}
/******************************************************************************/


/******************************************************************************/
/*        cPolyVecAlloc(nz,progcallp)                                         */
/******************************************************************************/

/******************************************************************************/
/*        cPolyPVecAlloc(vecp,nz)                                             */
/*                                                                            */
/******************************************************************************/
void      cPolyPVecAlloc(cPolyVec *vecp,size_t  nz)
{ 
  vecp->p = cPolyAlloc(nz,"cPolyPVecAlloc");
  vecp->z = nz;  vecp->x = 0;
  cPolyZero(vecp->p);
  return ;
}
/******************************************************************************/
/*        cPolyPVecRealloc(vecp,neednz,incrnz)                                */
/******************************************************************************/
void      cPolyPVecRealloc(cPolyVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = cPolyChkRealloc(vecp->p, &(vecp->z), neednz, incrnz, "cPolyPVecRealloc");
  return ;
}
/******************************************************************************/
/*        cPolyPVecFree(vecp)                                                 */
/******************************************************************************/
void      cPolyPVecFree(cPolyVec *vecp)
{ 
  free(vecp->p);  vecp->p = NULL;  vecp->z = 0;  vecp->x = 0;  vecp->gcp = NULL;
  return;
}
/******************************************************************************/
/*        cPolyVecFree(vecp)                                                  */
/******************************************************************************/
void      cPolyVecFree(cPolyVec *vecp)
{ 
  if(vecp->p != NULL) cPolyPVecFree(vecp);
  free(vecp);  return;
}
/******************************************************************************/
/*        cPolyVecInc1(vecp,y)                                                */
/* add 1 cPoly into a cPolyVec structure                                      */
/******************************************************************************/
void      cPolyVecInc1(cPolyVec *vecp, cPoly cc)
{ cPoly  *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cPolyPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  cPolyCopyp(p, cc);
  return ;
}
/******************************************************************************/
/*        cPolyVecInc1p(vecp,yp)                                              */
/* add 1 *cPoly into a cPolyVec structure                                     */
/******************************************************************************/
void      cPolyVecInc1p(cPolyVec *vecp, cPoly *ccp)
{ cPoly  *p;
  size_t    newz,oldx;

  oldx = vecp->x ; newz = oldx+1;
  cPolyPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx;  vecp->x = newz;
  cPolyCopypp(p, ccp);
  return;
}
/******************************************************************************/
/*        cPolyVecPrint(p)                                                    */
/******************************************************************************/
void      cPolyVecPrint(FILE *bufp, cPolyVec *p)
{ static    char    form1[] = "%s  cPolyVecp=%p, size(alloc, used)=(%d,%d), cPolyp=%p, gcp=%d\n";
  static    char    prognamp[] = "utiCurve.poly::cPolyVecPrint";

  fPrintF(bufp,form1, prognamp, (void*)p, p->z, p->x, (void*)p->p, (void*)p->gcp);
  return ;
}
/******************************************************************************/
/*        cPolyPVecPrint(p)                                                   */
/******************************************************************************/
void      cPolyPVecPrint(FILE *bufp, cPolyVec *vecp)
{ int       i;
  cPoly    *p;
  static    char    form1[] = "%s  cPolyVecp=%p, size(alloc, used)=(%d,%d), cPolyp=%p\n";
  static    char    prognamp[] = "utiCurve.poly::cPolyPVecPrint";

  fPrintF(bufp,form1, prognamp, (void*)vecp, vecp->z, vecp->x, (void*)vecp->p);
  p = vecp->p;
  for(i = 0;  i < vecp->x;  i++, p++){ cPolyPrint(bufp, p);}
  return ;
}
/******************************************************************************/
/******************************************************************************/
