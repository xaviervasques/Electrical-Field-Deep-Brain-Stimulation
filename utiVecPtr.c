/*  ../libmy/utiVecPtr.c                                                      */
/*  Mennessier Gerard                   940822                                */
/*  Last revised : M.G.                 991026                                */

#include  "utiVecPtr.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

/******************************************************************************/
/*        ptrVecAlloc(nz,prognamp)                                            */
/******************************************************************************/
/**
ptrVec   *ptrVecAlloc(size_t  nz,char *prognamp)
{ ptrVec   *vecp, *cvecp; 
  int       i;
  char      form1[]= "ptrVecAlloc ERROR; ptrVec size nz=%d \n" ;

  cvecp = vecp = (ptrVec *)malloc(nz*sizeof(ptrVec));
  if(vecp == NULL){ myErr1(-1,stderr,prognamp, form1, nz);}
  for(i=0; i<nz; i++){ cvecp->z=0; cvecp->x=0; cvecp->p=NULL; cvecp++;}
  return vecp;
}
**/
/******************************************************************************/
/*        ptrVecChkRealloc(p,nzp,neednz,incrnz,prognamp)                      */
/******************************************************************************/
/**
ptrVec   *ptrVecChkRealloc(void *p,size_t *nzp,size_t neednz,size_t incrnz, 
                                                                      char *prognamp)
{ char      form1[]= "ptrVecChkRealloc ERROR; desired ptrVec size nz=%d \n";

  if( neednz <= *nzp ){ return p;}
  p = (ptrVec *)realloc(p, (neednz+incrnz)*sizeof(ptrVec) );
  if(p == NULL){ myErr1(-1,stderr,prognamp, form1, neednz);}
  *nzp = (neednz+incrnz);  return p;
}
**/
/******************************************************************************/
/*        ptrPVecAlloc(vecp,nz)                                               */
/******************************************************************************/
void      ptrPVecAlloc(ptrVec *vecp,size_t  nz)
{  
  vecp->p = ptrAlloc( nz,"ptrPVecAlloc");  vecp->z=nz;  vecp->x=0;
  return ;
}
/******************************************************************************/
/*        ptrPVecRealloc(vecp,nz)                                             */
/******************************************************************************/
void      ptrPVecRealloc(ptrVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = ptrChkRealloc(vecp->p, &(vecp->z),neednz,incrnz, "ptrPVecRealloc");
  return ;
}
/******************************************************************************/
/*        ptrPVecFree(vecp)                                                   */
/******************************************************************************/
void      ptrPVecFree(ptrVec *vecp)
{  
  free(vecp->p);  vecp->p=NULL; vecp->z=0; vecp->x=0; 
  return ;
}
/******************************************************************************/
/*        ptrVecFree(vecp)                                                    */
/******************************************************************************/
void      ptrVecFree(ptrVec *vecp)
{ if(vecp == NULL){ return;}
  if(vecp->p != NULL){ ptrPVecFree(vecp); }
  free(vecp);  vecp = NULL;
  return ;
}
/******************************************************************************/
/*        ptrVecPrint(bufp,vecp)                                              */
/******************************************************************************/
void      ptrVecPrint(FILE  *bufp, ptrVec *vecp)
{ void    **pp, **ppx;
  char      form1[]="ptrVec.p=%p, vec.z=%d, vec.x=%d \n" ;

  fPrintF(bufp,form1, vecp->p,vecp->z,vecp->x);
  if(vecp->p == NULL){ return; }
  pp = vecp->p;  ppx = pp + vecp->x;
  for(; pp<ppx; pp++){ fPrintF(bufp, "%p ",*pp);}
  fPrintF(bufp, "\n");
  return ;
}
/******************************************************************************/
/*        ptrVecStrPrint(bufp,vecp)                                           */
/******************************************************************************/
void      ptrVecStrPrint(FILE  *bufp, ptrVec *vecp)
{ void    **pp, **ppx;

  if(vecp->p == NULL){ return; }
  pp = vecp->p;  ppx = pp + vecp->x;
  for(; pp<ppx; pp++){ fPrintF(bufp, "%s\n",(char*)*pp);}
  return ;
}
/******************************************************************************/
/*        ptrVecInc1(vecp,yp)                                                 */
/* add 1 ptr into a ptrVec structure                                          */
/******************************************************************************/
void      ptrVecInc1(ptrVec *vecp,void *yp)
{ void       **pp;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  ptrPVecRealloc(vecp,newz, newz/8); 
  pp = vecp->p;  *(pp+oldx)=yp; vecp->x=newz ;
  return ;
}
/******************************************************************************/
/*        ptrVecIncN(vecp,ypp,n)                                              */
/* add n ptr into a ptrVec structure                                          */
/******************************************************************************/
void      ptrVecIncN(ptrVec *vecp,void **ypp,size_t n)
{ void     **pp;
  size_t    newz,oldx;
  int       i;
  
  oldx = vecp->x;  newz = oldx +n;
  ptrPVecRealloc(vecp,newz, newz/8);
  pp = vecp->p + oldx; 
  for(i=n; i>0; i--,pp++,ypp++){ *pp = *ypp;}
  vecp->x = newz ;
  return ;
}
/******************************************************************************/
/******************************************************************************/
