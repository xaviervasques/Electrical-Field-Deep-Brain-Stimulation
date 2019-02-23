/*  ../libmy/utiVecUInt.c                                                     */
/*  Mennessier Gerard                   20030926                              */
/*  Last revised : M.G.                 20030926                              */

#include  "utiVecUInt.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

/******************************************************************************/
/*        uintVecAlloc(nz,prognamp)                                            */
/******************************************************************************/
/**
uintVec   *uintVecAlloc(size_t nz, char *progcallp)
{ uintVec   *vecp, *cvecp; 
  int       i;
  static    char      prognamp[] = "utiVecUInt::uintVecAlloc";
  static    char      form1p[]= "called from %s. malloc failed; uintVec size nz=%d\n";

  cvecp = vecp = (uintVec *)malloc(nz*sizeof(uintVec));
  if(vecp == NULL) myErr1(-1,stderr,prognamp, form1p,progcall,nz);
  for(i = 0;  i < nz;  i++){ cvecp->z = 0;  cvecp->x = 0; cvecp->p = NULL;  cvecp++;}
  return vecp;
}
**/
/******************************************************************************/
/*        uintPVecAlloc(vecp,nz)                                              */
/******************************************************************************/
void      uintPVecAlloc(uintVec *vecp, size_t nz)
{ static    char      prognamp[] = "utiVecUInt::uintPVecAlloc";

  vecp->p = (unsigned int *)intAlloc(nz, prognamp);
  vecp->z = nz;  vecp->x = 0;  *(vecp->p) = 0;
  return ;
}
/******************************************************************************/
/*        uintPVecRealloc(vecp,neednz,incrnz)                                 */
/******************************************************************************/
void      uintPVecRealloc(uintVec *vecp, size_t neednz, size_t incrnz)
{ static    char      prognamp[] = "utiVecUInt::uintPVecRealloc";

  vecp->p = (unsigned int *)intChkRealloc((int *)vecp->p, &(vecp->z),
                                                           neednz, incrnz, prognamp);
  return ;
}
/******************************************************************************/
/*        uintPVecFree(vecp)                                                  */
/******************************************************************************/
void      uintPVecFree(uintVec *vecp)
{  
  free(vecp->p);  vecp->p = NULL;  vecp->z = 0;  vecp->x = 0; 
  return ;
}
/******************************************************************************/
/*        uintVecFree(vecp)                                                   */
/******************************************************************************/
void      uintVecFree(uintVec *vecp)
{ 
  if(vecp == NULL) return;
  if(vecp->p != NULL) uintPVecFree(vecp);
  free(vecp);  vecp = NULL;
  return ;
}
/******************************************************************************/
/*        uintVecPrint(vecp)                                                  */
/******************************************************************************/
void      uintVecPrint(FILE  *bufp, uintVec *vecp)
{ unsigned int  *p, *px;
  static    char      form1p[] ="uintVec.p=%p, vec.z=%d, vec.x=%d \n" ;

  fPrintF(bufp,form1p, vecp->p,vecp->z,vecp->x);
  if(vecp->p == NULL) return;
  p = vecp->p;  px = p + vecp->x;
  for(;  p < px;  p++){ fPrintF(bufp, "%u ",*p);}
  fPrintF(bufp, "\n");
  return ;
}
/******************************************************************************/
/*        uintVecInc1(vecp,y)                                                 */
/* add 1 int into a uintVec structure                                         */
/******************************************************************************/
void      uintVecInc1(uintVec *vecp, unsigned int y)
{ unsigned int  *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx + 1;
  uintPVecRealloc(vecp,newz, newz/8);
  p = vecp->p;  *(p+oldx) = y;  vecp->x = newz;
  return;
}
/******************************************************************************/
/*        uintVecIncN(vecp,yp,n)                                              */
/* add n int into a uintVec structure                                         */
/******************************************************************************/
void      uintVecIncN(uintVec *vecp, unsigned int *yp, size_t n)
{ unsigned int  *p;
  size_t    newz,i,oldx;
  
  oldx = vecp->x;  newz = oldx + n;
  uintPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx; 
  for(i = 0;  i < n;  i++){ *p = *yp;  p++;  yp++;}
  vecp->x = newz;
  return;
}
/******************************************************************************/
/******************************************************************************/
