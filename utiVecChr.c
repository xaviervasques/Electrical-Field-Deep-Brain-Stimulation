/*  ../libmy/utiVecChr.c                                                      */
/*  Mennessier Gerard                   19940822                              */
/*  Last revised : M.G.                 20020213                              */

#include  "utiVecChr.h"

#include  "utistdErr.h"
#include  "utistdMem.h"                                            /** for memCpy **/
#include  "utistdStr.h"                                            /** for strLen **/
#include  "utiAlloc.h"

/******************************************************************************/
/*        chrVecAlloc(nz,prognamp)                                            */
/******************************************************************************/
/**
chrVec   *chrVecAlloc(size_t  nz,char *progcallp)
{ chrVec   *vecp , *cvecp; 
  int       i;
  static    char      prognamp[] = "utiVecChr::chrVecAlloc";
  static    char      form1p[]= "called from %s. malloc failed; chrVec size nz=%d\n";

  cvecp = vecp = (chrVec *)malloc( nz*sizeof(chrVec) );
  if(vecp == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  for(i=0; i<nz; i++){ cvecp->z=0; cvecp->x=0; cvecp->p=NULL; cvecp++; }
  return vecp;
}
**/
/******************************************************************************/
/*        chrPVecAlloc(vecp,nz)                                               */
/******************************************************************************/
void      chrPVecAlloc(chrVec *vecp,size_t nz)
{ char     *p;

  if(nz <= 0) nz = 1;
  p = chrAlloc(nz, "chrPVecAlloc");  vecp->z = nz;  vecp->x = 0;
  vecp->p = p;
  if(p != NULL) *p = 0;
  return ;
}
/******************************************************************************/
/*        chrPVecRealloc(vecp,neednz,incrnz)                                  */
/******************************************************************************/
void      chrPVecRealloc(chrVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = chrChkRealloc(vecp->p, &(vecp->z),neednz,incrnz,"chrPVecRealloc");
  return ;
}
/******************************************************************************/
/*        chrPVecFree(vecp)                                                   */
/******************************************************************************/
void      chrPVecFree(chrVec *vecp)
{  
  free(vecp->p);  vecp->p=NULL;  vecp->z=0;  vecp->x=0; 
  return ;
}

/******************************************************************************/
/*        chrVecFree(vecp)                                                    */
/******************************************************************************/
void      chrVecFree(chrVecp vecp)
{
  if(vecp->p != NULL){ chrPVecFree(vecp); }
  free(vecp);  vecp = NULL;
  return ;
}
/******************************************************************************/
/*        chrVecPrint(vecp)                                                   */
/******************************************************************************/
void      chrVecPrint(FILE  *bufp, chrVec *vecp)
{ char     *p, *px, c;
  char      form1[]="chrVec at %p, chrVec.p=%p, vec.z=%d, vec.x=%d \n" ;

  fPrintF(bufp,form1, vecp, vecp->p, vecp->z, vecp->x);
  if(vecp->p == NULL){ return; }
  p = vecp->p;  px = p + vecp->x;
  while(p < px)
  { c = *p++;
    if(c) fPrintF(bufp, "%c", c);
    else  break;
  }
  fPrintF(bufp, "\n");
  return ;
}
/******************************************************************************/
/*        chrVecXPrint(vecp)                                                  */
/******************************************************************************/
void      chrVecXPrint(FILE  *bufp, chrVec *vecp)
{ char     *p, *px;
  char      form1[]="chrVec at %p, chrVec.p=%p, vec.z=%d, vec.x=%d \n" ;

  fPrintF(bufp,form1, vecp, vecp->p,vecp->z,vecp->x);
  if(vecp->p == NULL){ return;}
  p = vecp->p;  px = p + vecp->x;
  for(;  p < px;  p++){ fPrintF(bufp, "%2x ",*p);}
  fPrintF(bufp, "\n");
  return ;
}
/******************************************************************************/
/*        chrVecInc1(vecp,c)                                                  */
/* add 1 char into a chrVec structure                                         */
/******************************************************************************/
void      chrVecInc1(chrVec *vecp,char c)
{ char     *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  chrPVecRealloc(vecp,newz, newz/8);
  p=vecp->p;  *(p+oldx)=c;  vecp->x=newz ;
  return ;
}
/******************************************************************************/
/*        chrVecIncN(vecp,memp,n)                                             */
/* add n char into a chrVec structure                                         */
/******************************************************************************/
void      chrVecIncN(chrVec *vecp,char *memp,size_t n)
{ char     *cp;
  size_t    newz,oldx;
  
  oldx = vecp->x;  newz = oldx +n;
  chrPVecRealloc(vecp,newz, newz/8);
  cp = vecp->p + oldx;  memCpy(cp, memp, n);  vecp->x = newz;
  return ;
}
/******************************************************************************/
/*        chrVecEnd(vecp,c)                                                   */
/* write int c (converted to char) in the last used memory                    */
/******************************************************************************/
void      chrVecEnd(chrVec *vecp,int c)
{ *(vecp->p + vecp->x -1) = (char)c;
  return ;
}
/******************************************************************************/
/*        chrVecIncStr(vecp,cp)                                               */
/* add a string into a chrVec structure                                       */
/******************************************************************************/
void      chrVecIncStr(chrVec *vecp,char *cp)
{ size_t    n;

  n = strLen(cp) +1;
  chrVecIncN(vecp,cp,n);
  return ;
}
/******************************************************************************/
/******************************************************************************/
