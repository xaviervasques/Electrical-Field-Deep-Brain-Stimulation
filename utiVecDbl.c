/*  ../libmy/utiVecDbl.c                                                      */
/*  Mennessier Gerard                   940822                                */
/*  Last revised              M.G.      970912                                */

#include  "utiVecDbl.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

/******************************************************************************/
/*        dblVecAlloc(nz,prognamp)                                            */
/******************************************************************************/
/**
dblVec   *dblVecAlloc(size_t  nz,char *prognamp)
{ dblVec   *vecp, *cvecp; 
  size_t    i;
  char      form1[]= "dblVecAlloc ERROR; dblVec size nz=%d \n" ;

  cvecp = vecp = (dblVec *)malloc( nz*sizeof(dblVec) );
  if(vecp==NULL){ myErr1(-1,stderr,prognamp, form1, nz);}
  for(i=0;i<nz;i++){ cvecp->z=0; cvecp->x=0; cvecp->p=NULL; cvecp++;}
  return vecp;
}
**/
/******************************************************************************/
/*        dblPVecAlloc(vecp,nz)                                               */
/******************************************************************************/
void      dblPVecAlloc(dblVec *vecp,size_t  nz)
{  
  vecp->p = dblAlloc( nz, "dblPVecAlloc" ); vecp->z=nz; vecp->x=0; *(vecp->p)=0;
  return ;
}
/******************************************************************************/
/*        dblPVecRealloc(vecp,neednz,incrnz)                                  */
/******************************************************************************/
void      dblPVecRealloc(dblVec *vecp,size_t neednz,size_t incrnz)
{ 
  vecp->p = dblChkRealloc(vecp->p, &(vecp->z),neednz,incrnz, "dblPVecRealloc");
  return ;
}
/******************************************************************************/
/*        dblPVecFree(vecp)                                                   */
/******************************************************************************/
void      dblPVecFree(dblVec *vecp)
{  
  free(vecp->p);  vecp->p=NULL; vecp->z=0; vecp->x=0; 
  return ;
}
/******************************************************************************/
/*        dblVecFree(vecp)                                                    */
/******************************************************************************/
void      dblVecFree(dblVec *vecp)
{
  if(vecp->p != NULL){ dblPVecFree(vecp); }
  free(vecp);  vecp = NULL;
  return ;
}
/******************************************************************************/
/*        dblVecPrint(vecp)                                                   */
/******************************************************************************/
void      dblVecPrint(FILE  *bufp, dblVec *vecp)
{ double   *p, *px;
  char    form1[]="dblVec.p=%p, vec.z=%d, vec.x=%d \n" ;

  fPrintF(bufp,form1, vecp->p,vecp->z,vecp->x);
  if(vecp->p == NULL){ return;}
  p=vecp->p;  px = p + vecp->x;
  for(; p<px; p++){ fPrintF(bufp, "%e ",*p);}
  fPrintF(bufp, "\n");
  return ;
}
/******************************************************************************/
/*        dblVecInc1(vecp,y)                                                  */
/* add 1 double into a dblVec structure                                       */
/******************************************************************************/
void      dblVecInc1(dblVec *vecp,double y)
{ double   *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx+1;
  dblPVecRealloc(vecp,newz, newz/8);
  p = vecp->p;  *(p+oldx) = y;  vecp->x = newz;
  return ;
}
/******************************************************************************/
/*        dblVecIncN(vecp,yp,n)                                               */
/* add n double into a dblVec structure                                       */
/******************************************************************************/
void      dblVecIncN(dblVec *vecp,double *yp,size_t n)
{ double   *p;
  size_t    newz,i,oldx;
  
  oldx = vecp->x ; newz = oldx +n;
  dblPVecRealloc(vecp,newz, newz/8);
  p=vecp->p + oldx ; 
  for(i=0; i<n; i++){ *p=*yp; p++; yp++; }
  vecp->x = newz ;
  return ;
}
/******************************************************************************/
/******************************************************************************/
