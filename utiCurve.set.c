/*  essai/C/utiCurve.set.c                                                    */
/*  Mennessier Gerard                   20010507                              */
/*  Last revised M.G.                   20020911                              */
/*                                                                            */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

#include  "utiCurve.set.h"
#include  "utiCurve.GC.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cSet     *cSetAlloc(size_t nz, char *progcallp)
{ cSet     *p;
  static char    prognamp[] = "utiCurve.set::cSetAlloc";
  static char    form1p[] = "called from %s. malloc failed; cSet size nz=%d\n";

  p = (cSet *)malloc( nz*sizeof(cSet) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  return p;
}
/******************************************************************************/
/*  cSet *cSetChkRealloc(cSetp,nzp,neednz,incrnz,progcallp)                   */
/*                                                                            */
/* check if there is enough memory allocated for the neednz needed cSet       */
/* old and possibly new size  stored as  *nzp                                 */
/* if reallocation is needed, newsize =  neednz + incrnz                      */
/******************************************************************************/
cSet     *cSetChkRealloc(cSet *pi, size_t *nzp, size_t neednz, size_t incrnz,
                                                                      char *progcallp)
{ cSet     *p;
  size_t    nfz;
  static char    prognamp[] = "utiCurve.set::cSetChkRealloc";
  static char    form1p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired cSet size nz=%d\n";
  static char    form2p[]= "called from %s. malloc failed;"
                                                       " desired cSet size nz=%d\n";
  
  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (cSet *)malloc( nfz*sizeof(cSet) );
fPrintF(stderr,"%s called from %s; realloc p=%p\n", prognamp, progcallp, (void*)pi);
    if(p == NULL) myErr1(-1,stderr,prognamp, form2p, progcallp, neednz);
  }
  else
  { if(neednz  <=  *nzp) return pi;
    p = (cSet *)realloc(pi, nfz*sizeof(cSet) );
    if(p == NULL) myErr1(-1,stderr,prognamp, form1p, progcallp, pi, neednz);
  }
  *nzp = nfz;
  return p;
}
/******************************************************************************/
/*  void  cSetPrint(FILE *bufp,cSet *p)                                       */
/******************************************************************************/
void      cSetPrint(FILE *bufp, cSet *p)
{ static char    prognamp[] = "utiCurve.set::cSetPrint";

  fPrintF(bufp,"%s cSet address=%p, pointer=%p, type=%d\n",
                                                  prognamp, (void*)p, p->p, p->type);
  return;
}
/******************************************************************************/
/*  void  cSetArrayPrint(FILE *bufp,cSet *p,size_t n)                         */
/******************************************************************************/
void      cSetArrayPrint(FILE *bufp, cSet *p, size_t n)
{ n++;
  while(--n)
  { fPrintF(bufp,"pointer=%p, type=%d\n", p->p, p->type);
    p++;
  }
  return;
}
/******************************************************************************/
/*  void  cSetZero(cSet *p)                                                   */
/******************************************************************************/
void      cSetZero(cSet *p)
{ 
  p->p = NULL; p->type = 0;  return;
}
/******************************************************************************/
/*  void  cSetCpy1p(cSet *sfp, cSet *sip)                                     */
/*                                                                            */
/* copy 1 cSet  *sip into  *sfp                                               */
/******************************************************************************/
void      cSetCpy1p(cSet *sfp, cSet *sip)
{ 
  sfp->p = sip->p;  sfp->type = sip->type;   return;
}
/******************************************************************************/


/******************************************************************************/
/*        cSetVecAlloc(nz,progcallp)                                          */
/******************************************************************************/

/******************************************************************************/
/*        cSetPVecAlloc(vecp,nz)                                              */
/******************************************************************************/
void      cSetPVecAlloc(cSetVec *vecp, size_t  nz)
{ 
  vecp->p = cSetAlloc(nz,"cSetPVecAlloc");
  vecp->z = nz;  vecp->x = 0;
  cSetZero(vecp->p);
  return ;
}
/******************************************************************************/
/*        cSetPVecRealloc(vecp,neednz,incrnz)                                 */
/******************************************************************************/
void      cSetPVecRealloc(cSetVec *vecp, size_t neednz, size_t incrnz)
{ 
  vecp->p = cSetChkRealloc(vecp->p, &(vecp->z), neednz, incrnz, "cSetPVecRealloc");
  return ;
}
/******************************************************************************/
/*        cSetPVecFree(vecp)                                                  */
/******************************************************************************/
void      cSetPVecFree(cSetVec *vecp)
{ 
  free(vecp->p);  vecp->p = NULL;  vecp->z = 0;  vecp->x = 0;  return ;
}
/******************************************************************************/
/*        cSetVecFree(vecp)                                                   */
/******************************************************************************/
void      cSetVecFree(cSetVec *vecp)
{ 
  if(vecp->p != NULL) cSetPVecFree(vecp);
  free(vecp);  return ;
}
/******************************************************************************/
/*        cSetVecPrint(vecp)                                                  */
/******************************************************************************/
void      cSetVecPrint(FILE *bufp, cSetVec *vecp)
{ cSet     *p;
  static char    form1[] = "%s  cSetVecp=%p, cSetVec.p=%p, vec.z=%d, vec.x=%d \n" ;
  static char    prognamp[] = "utiCurve.set::cSetVecPrint";

  p = vecp->p;
  fPrintF(bufp, form1, prognamp, vecp, p, vecp->z, vecp->x);
  if(p == NULL) return;
  cSetArrayPrint(bufp, p, vecp->x);
  fPrintF(bufp, "\n");  return ;
}
/******************************************************************************/
/*        cSetVecInc1(vecp,y)                                                 */
/* add 1 cSet into a cSetVec structure                                        */
/******************************************************************************/
void      cSetVecInc1(cSetVec *vecp, cSet s)
{ cSet     *csetp;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  cSetPVecRealloc(vecp, newz, newz/8);
  csetp = vecp->p + oldx;  vecp->x = newz;
  csetp->p = s.p;  csetp->type = s.type;
  return ;
}
/******************************************************************************/
/*        cSetVecInc1p(vecp,yp)                                               */
/* add 1 *cSet into a cSetVec structure                                       */
/******************************************************************************/
void      cSetVecInc1p(cSetVec *vecp, cSet *sp)
{ cSet     *csetp;
  size_t    newz,oldx;

  oldx = vecp->x ; newz = oldx+1;
  cSetPVecRealloc(vecp,newz, newz/8);
  csetp = vecp->p + oldx;  vecp->x = newz;
  cSetCpy1p(csetp, sp);  return;
}
/******************************************************************************/
/*        cSetVecInc1ptype(cSetVec *vecp, void *p, int type)                  */
/* add the cSet = {p,type} into a cSetVec structure                           */
/******************************************************************************/
void      cSetVecInc1ptype(cSetVec *vecp, void *p, int type)
{ cSet     *csetp;
  size_t    newz, oldx;

  oldx = vecp->x ; newz = oldx + 1;
  cSetPVecRealloc(vecp,newz, newz/8);
  csetp = vecp->p + oldx;  vecp->x = newz;  
  csetp->p = p;  csetp->type = type;
  return;
}
/******************************************************************************/
/*        cSetVecAddSetVec(cSetVec *vecp, cSetVec *vecnewp)                   */
/* add sets pointed by vecnewp into vecp                                      */
/******************************************************************************/
void      cSetVecAddSetVec(cSetVec *vecp, cSetVec *vecnewp)
{ cSet     *csfp, *csip;
  size_t    newz, oldx;
  int       i, ix;

  oldx = vecp->x ; newz = oldx + vecnewp->x;
  cSetPVecRealloc(vecp,newz, newz/8);
  csfp = vecp->p + oldx;  csip = vecnewp->p;
  ix = vecnewp->x;
  for(i = 0;  i < ix;  i++,  csfp++, csip++)
  { csfp->p = csip->p;  csfp->type = csip->type;
  }
  vecp->x = newz;
  return;
}
/******************************************************************************/
/******************************************************************************/
