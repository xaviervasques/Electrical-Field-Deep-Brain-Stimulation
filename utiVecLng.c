/*  ../libmy/utiVecLng.c                                                      */
/*  Mennessier Gerard                   940822                                */
/*  Last revised : M.G.                 991026                                */

#include  "utiVecLng.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

/******************************************************************************/
/*        lngVecAlloc(nz,prognamp)  ->  (lngVec *)myVecAlloc()                */
/******************************************************************************/
/**
lngVec   *lngVecAlloc(size_t  nz,char *prognamp)
{ lngVec   *vecp, *cvecp; 
  size_t    i;
  char      form1[]= "lngVecAlloc ERROR; lngVec size nz=%d \n" ;

  cvecp = vecp = (lngVec *)malloc(nz*sizeof(lngVec));
  if(vec == NULL){ myErr1(-1,stderr,prognamp, form1, nz);}
  for(i=0; i<nz; i++){ cvecp->z=0; cvecp->x=0; cvecp->p=NULL; cvecp++;}
  return vecp;
}
**/
/******************************************************************************/
/*        lngPVecAlloc(vecp,nz)                                               */
/******************************************************************************/
void      lngPVecAlloc(lngVec *vecp,size_t  nz)
{  
  vecp->p = lngAlloc( nz,"lngPVecAlloc");  vecp->z=nz; vecp->x=0; *(vecp->p)=0;
  return ;
}
/******************************************************************************/
/*        lngPVecRealloc(vecp,neednz,incrnz)                                  */
/******************************************************************************/
void      lngPVecRealloc(lngVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = lngChkRealloc(vecp->p, &(vecp->z),neednz,incrnz,"lngPVecRealloc");
  return ;
}
/******************************************************************************/
/*        lngPVecFree(vecp)                                                   */
/******************************************************************************/
void      lngPVecFree(lngVec *vecp)
{  
  free(vecp->p);  vecp->p=NULL; vecp->z=0; vecp->x=0; 
  return ;
}
/******************************************************************************/
/*        lngVecFree(vecp)                                                    */
/******************************************************************************/
void      lngVecFree(lngVec *vecp)
{ if(vecp == NULL){ return;}
  if(vecp->p != NULL){ lngPVecFree(vecp); }
  free(vecp);  vecp = NULL;
  return ;
}
/******************************************************************************/
/*        lngVecPrint(vecp)                                                   */
/******************************************************************************/
void      lngVecPrint(FILE  *bufp, lngVec *vecp)
{ long     *p, *px;
  char    form1[]="lngVec.p=%p, vec.z=%d, vec.x=%d \n" ;

  fPrintF(bufp,form1, vecp->p,vecp->z,vecp->x);
  if(vecp->p == NULL){ return;}
  p = vecp->p;  px = p + vecp->x;
  for(; p<px; p++){ fPrintF(bufp, "%ld ",*p);}
  fPrintF(bufp, "\n");
  return ;
}
/******************************************************************************/
/*        lngVecInc1(vecp,y)                                                  */
/* add 1 long into a lngVec structure                                         */
/******************************************************************************/
void      lngVecInc1(lngVec *vecp,long y)
{ long       *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx+1;
  lngPVecRealloc(vecp,newz, newz/8);
  p=vecp->p;  *(p+oldx)=y; vecp->x=newz;
  return ;
}
/******************************************************************************/
/*        lngVecIncN(vecp,yp,n)                                               */
/* add n long into a lngVec structure                                         */
/******************************************************************************/
void      lngVecIncN(lngVec *vecp,long *yp,size_t n)
{ long     *p;
  size_t    newz,i,oldx;
  
  oldx = vecp->x ; newz = oldx +n;
  lngPVecRealloc(vecp,newz, newz/8);
  p=vecp->p + oldx ; 
  for(i=0; i<n; i++){ *p=*yp; p++; yp++; }
  vecp->x = newz ;
  return ;
}
/******************************************************************************/
/******************************************************************************/
