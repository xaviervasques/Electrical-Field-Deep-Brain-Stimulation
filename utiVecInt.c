/*  ../libmy/utiVecInt.c                                                      */
/*  Mennessier Gerard                   940822                                */
/*  Last revised : M.G.                 991026                                */

#include  "utiVecInt.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

/******************************************************************************/
/*        intVecAlloc(nz,prognamp)                                            */
/******************************************************************************/
/**
intVec   *intVecAlloc(size_t  nz,char *progcallp)
{ intVec   *vecp, *cvecp; 
  int       i;
  static    char      prognamp[] = "utiVecInt::intVecAlloc";
  static    char      form1p[]= "called from %s. malloc failed; intVec size nz=%d\n";

  cvecp = vecp = (intVec *)malloc(nz*sizeof(intVec));
  if(vecp == NULL) myErr1(-1,stderr,prognamp, form1p,progcall,nz);
  for(i = 0;  i < nz;  i++){ cvecp->z = 0;  cvecp->x = 0; cvecp->p = NULL; cvecp++;}
  return vecp;
}
**/
/******************************************************************************/
/*        intPVecAlloc(vecp,nz)                                               */
/******************************************************************************/
void      intPVecAlloc(intVec *vecp,size_t  nz)
{ static    char      prognamp[] = "utiVecInt::intPVecAlloc";

  vecp->p = intAlloc(nz,prognamp);
  vecp->z = nz;  vecp->x = 0;  *(vecp->p) = 0;
  return ;
}
/******************************************************************************/
/*        intPVecRealloc(vecp,neednz,incrnz)                                  */
/******************************************************************************/
void      intPVecRealloc(intVec *vecp,size_t neednz,size_t incrnz)
{ static    char      prognamp[] = "utiVecInt::intPVecRealloc";

  vecp->p = intChkRealloc(vecp->p, &(vecp->z),neednz,incrnz,prognamp);
  return ;
}
/******************************************************************************/
/*        intPVecFree(vecp)                                                   */
/******************************************************************************/
void      intPVecFree(intVec *vecp)
{  
  free(vecp->p);  vecp->p = NULL;  vecp->z = 0;  vecp->x = 0; 
  return ;
}
/******************************************************************************/
/*        intVecFree(vecp)                                                    */
/******************************************************************************/
void      intVecFree(intVec *vecp)
{ 
  if(vecp == NULL) return;
  if(vecp->p != NULL) intPVecFree(vecp);
  free(vecp);  vecp = NULL;
  return ;
}
/******************************************************************************/
/*        intVecPrint(vecp)                                                   */
/******************************************************************************/
void      intVecPrint(FILE  *bufp, intVec *vecp)
{ int      *p, *px;
  static    char      form1p[] ="intVec.p=%p, vec.z=%d, vec.x=%d \n" ;

  fPrintF(bufp,form1p, vecp->p,vecp->z,vecp->x);
  if(vecp->p == NULL) return;
  p = vecp->p;  px = p + vecp->x;
  for(;  p < px;  p++){ fPrintF(bufp, "%d ",*p);}
  fPrintF(bufp, "\n");
  return ;
}
/******************************************************************************/
/*        intVecInc1(vecp,y)                                                  */
/* add 1 int into a intVec structure                                          */
/******************************************************************************/
void      intVecInc1(intVec *vecp,int y)
{ int       *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx + 1;
  intPVecRealloc(vecp,newz, newz/8);
  p = vecp->p;  *(p+oldx) = y;  vecp->x = newz;
  return;
}
/******************************************************************************/
/*        intVecIncN(vecp,yp,n)                                               */
/* add n int into a intVec structure                                          */
/******************************************************************************/
void      intVecIncN(intVec *vecp,int *yp,size_t n)
{ int     *p;
  size_t    newz,i,oldx;
  
  oldx = vecp->x;  newz = oldx + n;
  intPVecRealloc(vecp,newz, newz/8);
  p = vecp->p + oldx; 
  for(i = 0;  i < n;  i++){ *p = *yp;  p++;  yp++;}
  vecp->x = newz;
  return;
}
/******************************************************************************/
/******************************************************************************/
