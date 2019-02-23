/*  ../libmy/utiVecMyVec.c                                                    */
/*  Mennessier Gerard                   970513                                */
/*  Last revised : M.G.                 991026                                */

#include  "utiVecMyVec.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

                    /************************************/
                    /**       myVec structure          **/
                    /************************************/

/******************************************************************************/
/*  void  myVecPrint(FILE *bufp,myVec *p)                                     */
/******************************************************************************/
void      myVecPrint(FILE *bufp,myVec *p)
{ char      form1p[]= "(%p,%d,%d)";

  fPrintF(bufp,form1p, p->p, p->z, p->x);  return;
}
/******************************************************************************/
/*  void  myVecArrayPrint(FILE *bufp,myVec *p,size_t n)                       */
/******************************************************************************/
void      myVecArrayPrint(FILE *bufp,myVec *p,size_t n)
{ char      form1p[]= "(%p,%d,%d)";

  n++;
  while(--n){ fPrintF(bufp,form1p, p->p, p->z, p->x);  p++;}
  return;
}
/******************************************************************************/
/*  void  myVecZero(myVec *p)                                                 */
/******************************************************************************/
void      myVecZero(myVec *p)
{ p->p=NULL;  p->z=0;  p->x=0;  return;
}
/******************************************************************************/
/*  void  myVecCpy1F(myVecp xfp,void *p,size_t z,size_t x)                    */
/*                                                                            */
/******************************************************************************/
void      myVecCpy1F(myVec *xfp,void *p,size_t z,size_t x)
{ xfp->p=p;  xfp->z=z;  xfp->x=x;  return;
}
/******************************************************************************/
/*  void  myVecCpy1(myVecp xfp,myVec xi)                                      */
/*                                                                            */
/* copy 1 myVec  xi into  *xfp                                                */
/******************************************************************************/
void      myVecCpy1(myVec *xfp,myVec xi)
{ xfp->p = xi.p;  xfp->z = xi.z;  xfp->x = xi.x;  return;
}
/******************************************************************************/
/*  void  myVecCpy1p(myVec *xfp,myVec *xip)                                   */
/*                                                                            */
/* copy 1 myVec  *xip into  *xfp                                              */
/******************************************************************************/
void      myVecCpy1p(myVec *xfp,myVec *xip)
{ xfp->p = xip->p;  xfp->z = xip->z;  xfp->x = xip->x;  return;
}
/******************************************************************************/


                    /************************************/
                    /**       myVec pointing on myVec  **/
                    /************************************/

/******************************************************************************/
/*        myVecPVecAlloc(vecp,nz)                                             */
/******************************************************************************/
void      myVecPVecAlloc(myVec *vecp,size_t  nz)
{ vecp->p = (void*)myVecAlloc(nz, "myVecPVecAlloc");
  vecp->z = nz;  vecp->x = 0;  myVecZero(vecp->p);
  return ;
}
/******************************************************************************/
/*        myVecPVecRealloc(vecp,neednz,incrnz)                                */
/******************************************************************************/
void      myVecPVecRealloc(myVec *vecp,size_t neednz,size_t incrnz)
{ vecp->p = (void*)myVecChkRealloc(vecp->p,&(vecp->z),neednz,incrnz,
                                                                 "myVecPVecRealloc");
  return ;
}
/******************************************************************************/
/*        myVecPVecFree(vecp)                                                 */
/******************************************************************************/
void      myVecPVecFree(myVec *vecp)
{ void     *p;

  p = vecp->p;  if(p != NULL){ free(p);  vecp->p = NULL;}
  vecp->z = 0;  vecp->x = 0; 
  return ;
}
/******************************************************************************/
/*        myVecVecInc1(vecp,y)                                                */
/* add 1 myVec into a myVecVec structure                                      */
/******************************************************************************/
void      myVecVecInc1(myVec *vecp,myVec y)
{ myVec    *p;
  size_t    newz,oldx;

  oldx = vecp->x;  newz = oldx +1;
  myVecPVecRealloc(vecp,newz, newz/8);
  p = (myVec*)vecp->p;
  p += oldx;
  p->p = y.p;  p->z = y.z;  p->x = y.x;
  vecp->x = newz;
  return ;
}
/******************************************************************************/
/*        myVecVecIncN(vecp,yp,n)                                             */
/* add n myVec into a myVecVec structure                                      */
/******************************************************************************/
void      myVecVecIncN(myVec *vecp,myVec *yp,size_t n)
{ myVec    *p;
  size_t    newz,i,oldx;
  
  oldx = vecp->x ; newz = oldx +n;
  myVecPVecRealloc(vecp,newz, newz/8);
  p=(myVec*)(vecp->p) + oldx ; 
  for(i=0; i<n; i++){ p->p = yp->p;  p->z = yp->z;  p->x = yp->x;  p++; yp++; }
  vecp->x = newz ;
  return ;
}
/******************************************************************************/
/******************************************************************************/
