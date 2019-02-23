/*  essai/C/utiCurve.string.c                                                 */
/*  Mennessier Gerard                   20010507                              */
/*  Last revised M.G.                   20020904                              */
/*                                                                            */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

#include  "utiCurve.string.h"
#include  "utiCurve.GC.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
cStr     *cStrAlloc(size_t  nz,char *progcallp)
{ cStr     *p;
  static    char      prognamp[] = "utiCurve.string::cStrAlloc";
  static    char      form1p[] = "called from %s. malloc failed; cStr size nz=%d\n";

  p = (cStr *)malloc( nz*sizeof(cStr) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p, progcallp, nz);
  return p;
}
/******************************************************************************/
/*  cStr *cStrChkRealloc(cStr,nzp,neednz,incrnz,progcallp)                   */
/*                                                                            */
/* check if there is enough memory allocated for the neednz needed cPoly      */
/* old and possibly new size  stored as  *nzp                                 */
/* if reallocation is needed, newsize =  neednz + incrnz                      */
/******************************************************************************/
cStr     *cStrChkRealloc(cStr *pi, size_t *nzp, size_t neednz, size_t incrnz,
                                                                     char *progcallp)
{ cStr     *p;
  size_t    nfz;
  static    char      prognamp[] = "utiCurve.string::cStrChkRealloc";
  static    char      form1p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired cStr size nz=%d\n";
  static    char      form2p[]= "called from %s. malloc failed;"
                                                       " desired cStr size nz=%d\n";
  
  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (cStr *)malloc( nfz*sizeof(cStr) );
fPrintF(stderr,"%s called from %s; realloc p=%p\n", prognamp, progcallp, (void*)pi);
    if(p == NULL) myErr1(-1,stderr,prognamp, form2p, progcallp, neednz);
  }
  else
  { if(neednz  <=  *nzp) return pi;
    p = (cStr *)realloc(pi, nfz*sizeof(cStr) );
    if(p == NULL) myErr1(-1,stderr,prognamp, form1p, progcallp, pi, neednz);
  } 
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*  void  cStrZero(cStr *p)                                                   */
/******************************************************************************/
void      cStrZero(cStr *p)
{ 
  p->xyp[0] = 0.; p->xyp[1] = 0.;
  return;
}
/******************************************************************************/


/******************************************************************************/
/*        cStrVecAlloc(nz,progcallp)                                          */
/******************************************************************************/
cStrVec  *cStrVecAlloc(size_t nz, char *progcallp)
{ cStrVec  *p, *cp;
  int       i;
  static char    form1p[] = "called from %s. malloc failed; cStrVec size nz=%d\n";
  static char    prognamp[] = "utiCurve.string::cStrVecAlloc";

  cp = p = (cStrVec *)malloc( nz*sizeof(cStrVec) );
  if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, nz);
  for(i = nz;  i > 0;  i--, cp++)
  { cp->z = 0;  cp->x = 0;  cp->p = NULL;
    cp->gcp = NULL;
  }
  return p;
}
/******************************************************************************/
/*        cStrVecChkRealloc(p,nzp,neednz,incrnz,prognamp)                     */
/******************************************************************************/
cStrVec  *cStrVecChkRealloc(cStrVec *pi, size_t *nzp, size_t neednz, size_t incrnz, 
                                                                     char *progcallp)
{ cStrVec  *cp, *p;
  size_t    nfz;
  int       i;
  static char    form2p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired cStrVec size nz=%d\n";
  static char    prognamp[] = "utiCurve.string::cStrVecChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL) p = cStrVecAlloc(nfz, prognamp);
  else
  { if(neednz <= *nzp) return pi;
    p = (cStrVec *)realloc(pi, nfz*sizeof(cStrVec) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
    cp = p + *nzp;
    for(i = (nfz - *nzp);  i > 0;  i--, cp++)
    { cp->z = 0;  cp->x = 0;  cp->p = NULL;
      cp->gcp = NULL;
    }
  }
  *nzp = nfz;  return p;
}

/******************************************************************************/
/*        cStrPVecAlloc(vecp,nz)                                              */
/******************************************************************************/
void      cStrPVecAlloc(cStrVec *vecp,size_t  nz)
{ 
  vecp->p = cStrAlloc(nz,"cStrPVecAlloc");
  vecp->z = nz;  vecp->x = 0;
  cStrZero(vecp->p);
  return ;
}
/******************************************************************************/
/*        cStrPVecRealloc(vecp,neednz,incrnz)                                 */
/******************************************************************************/
void      cStrPVecRealloc(cStrVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = cStrChkRealloc(vecp->p, &(vecp->z),neednz,incrnz,"cStrPVecRealloc");
  return ;
}
/******************************************************************************/
/*        cStrPVecFree(vecp)                                                  */
/******************************************************************************/
void      cStrPVecFree(cStrVec *vecp)
{ 
  free(vecp->p);  vecp->p = NULL;  vecp->z = 0;  vecp->x = 0;  return ;
}
/******************************************************************************/
/*        cStrVecFree(vecp)                                                   */
/******************************************************************************/
void      cStrVecFree(cStrVec *vecp)
{ 
  if(vecp->p != NULL) cStrPVecFree(vecp);
  free(vecp);  return ;
}
/******************************************************************************/
/******************************************************************************/
