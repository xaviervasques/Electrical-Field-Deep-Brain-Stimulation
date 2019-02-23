/*   ../libmy/utiAlloc.c                                                      */
/*   Mennessier Gerard                  19940822                              */
/*   Last revised : M.G.                20030128                              */

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

/* static char    form10p[] = "%s called from %s; realloc p=%p\n"; */
/******************************************************************************/
/*        chrAlloc(nz,prognamp)                                               */
/******************************************************************************/
char     *chrAlloc(size_t nz, char *progcallp)
{ char     *p;
  static char    form1p[] = "called from %s. malloc failed;"
                                                        " desired char size nz=%d\n";
  static char    prognamp[] = "utiAlloc::chrAlloc";

  p = (char *)malloc( nz*sizeof(char) );
  if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, nz);
  return p;
}
/******************************************************************************/
/*        chrChkRealloc(p,nzp,neednz,incrnz,prognamp)                         */
/******************************************************************************/
char     *chrChkRealloc(char *pi, size_t *nzp, size_t neednz, size_t incrnz,
                                                                     char *progcallp)
{ char     *p;
  size_t    nfz;
  static char    form2p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired char size nz=%d\n";
  static char    form1p[]= "called from %s. malloc failed;"
                                                       " desired char size nz=%d\n";
  static char    prognamp[] = "utiAlloc::chrChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (char *)malloc( nfz*sizeof(char) );
/* fPrintF(stderr, form10p, prognamp, progcallp, (void*)pi); */
    if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, neednz);
  }
  else
  { if(neednz <= *nzp ) return pi;
    p = (char *)realloc(pi, nfz*sizeof(char) );
    if(p == NULL) myErr1(-1,stderr,prognamp, form2p,progcallp,pi,neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*        dblAlloc(nz,prognamp)                                               */
/******************************************************************************/
double   *dblAlloc(size_t nz, char *progcallp)
{ double   *p;
  static char    form1p[] = "called from %s. malloc failed; double size nz=%d\n";
  static char    prognamp[] = "utiAlloc::dblAlloc";

  p = (double *)malloc( nz*sizeof(double) );
  if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, nz);
  return p;
}
/******************************************************************************/
/*        dblChkRealloc(p,nzp,neednz,incrnz,prognamp)                         */
/******************************************************************************/
double   *dblChkRealloc(double *pi, size_t *nzp, size_t neednz, size_t incrnz,
                                                                     char *progcallp)
{ double   *p;
  size_t    nfz;
  static char    form2p[] = "called from %s. realloc failed;"
                                  " initial pointer=%p, desired double size nz=%d\n";
  static char    form1p[] = "called from %s. malloc failed;"
                                                      " desired double size nz=%d\n";
  static char    prognamp[] = "utiAlloc::dblChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (double *)malloc( nfz*sizeof(double) );
/* fPrintF(stderr, form10p, prognamp, progcallp, (void*)pi); */
    if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, neednz);
  }
  else
  { if(neednz <= *nzp) return pi;
    p = (double *)realloc(pi, nfz*sizeof(double) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*        fltAlloc(nz,prognamp)                                               */
/******************************************************************************/
float    *fltAlloc(size_t nz, char *progcallp)
{ float    *p;
  static char    form1p[] = "called from %s. malloc failed; float size nz=%d\n";
  static char    prognamp[] = "utiAlloc::fltAlloc";

  p = (float *)malloc( nz*sizeof(float) );
  if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, nz);
  return p;
}
/******************************************************************************/
/*        fltChkRealloc(p,nzp,neednz,incrnz,prognamp)                         */
/******************************************************************************/
float    *fltChkRealloc(float *pi, size_t *nzp, size_t neednz, size_t incrnz,
                                                                     char *progcallp)
{ float   *p;
  size_t    nfz;
  static char    form2p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired float size nz=%d\n";
  static char    form1p[] = "called from %s. malloc failed;"
                                                       " desired float size nz=%d\n";
  static char    prognamp[] = "utiAlloc::fltChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (float *)malloc( nfz*sizeof(float) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, neednz);
  }
  else
  { if(neednz <= *nzp) return pi;
    p = (float *)realloc(pi, nfz*sizeof(float) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*        intAlloc(nz,prognamp)                                               */
/******************************************************************************/
int      *intAlloc(size_t nz, char *progcallp)
{ int      *p;
  static char    prognamp[] = "utiAlloc::intAlloc";
  static char    form1p[] = "called from %s. malloc failed; int size nz=%d\n";

  p = (int *)malloc( nz*sizeof(int) );
  if(p == NULL) myErr1(-1,stderr, prognamp, form1p, progcallp, nz);
  return p;
}
/******************************************************************************/
/*        intChkRealloc(p,nzp,neednz,incrnz,prognamp)                         */
/******************************************************************************/
int      *intChkRealloc(int *pi, size_t *nzp, size_t neednz, size_t incrnz,
                                                                     char *progcallp)
{ int      *p;
  size_t    nfz;
  static char    form2p[] = "called from %s. realloc failed;"
                                     " initial pointer=%p, desired int size nz=%d\n";
  static char    form1p[]= "called from %s. malloc failed;"
                                                         " desired int size nz=%d\n";
  static char    prognamp[] = "utiAlloc::intChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (int *)malloc( nfz*sizeof(int) );
/* fPrintF(stderr, form10p, prognamp, progcallp, (void*)pi); */
    if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, neednz);
  }
  else
  { if(neednz <= *nzp ) return pi;
    p = (int *)realloc(pi, nfz*sizeof(int) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*        uintAlloc(nz,prognamp)                                              */
/******************************************************************************/
unsigned int  *uintAlloc(size_t nz, char *progcallp)
{ unsigned int  *p;
  static char    form1p[] = "called from %s. malloc failed; int size nz=%d\n";
  static char    prognamp[] = "utiAlloc::uintAlloc";

  p = (unsigned int *)malloc( nz*sizeof(unsigned int) );
  if(p == NULL) myErr1(-1,stderr, prognamp, form1p, progcallp, nz);
  return p;
}
/******************************************************************************/
/*        uintChkRealloc(p,nzp,neednz,incrnz,prognamp)                        */
/******************************************************************************/
unsigned int  *uintChkRealloc(unsigned int *pi, size_t *nzp, 
                                       size_t neednz, size_t incrnz, char *progcallp)
{ unsigned int  *p;
  size_t    nfz;
  static char    form2p[] = "called from %s. realloc failed;"
                                     " initial pointer=%p, desired int size nz=%d\n";
  static char    form1p[]= "called from %s. malloc failed;"
                                                         " desired int size nz=%d\n";
  static char    prognamp[] = "utiAlloc::uintChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (unsigned int *)malloc( nfz*sizeof(unsigned int) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, neednz);
  }
  else
  { if(neednz <= *nzp ) return pi;
    p = (unsigned int *)realloc(pi, nfz*sizeof(unsigned int) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*        lngAlloc(nz,prognamp)                                               */
/******************************************************************************/
long     *lngAlloc(size_t nz, char *progcallp)
{ long     *p; 
  static char    prognamp[] = "utiAlloc::lngAlloc";
  static char    form1p[] = "called from %s. malloc failed; long size nz=%d\n";

  p = (long *)malloc( nz*sizeof(long) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p, progcallp, nz);
  return p;
}
/******************************************************************************/
/*        lngChkRealloc(p,nzp,neednz,incrnz,prognamp)                         */
/******************************************************************************/
long     *lngChkRealloc(long *pi, size_t *nzp, size_t neednz, size_t incrnz,
                                                                     char *progcallp)
{ long     *p;
  size_t    nfz;
  static char    form2p[] = "called from %s. realloc failed;"
                                    " initial pointer=%p, desired long size nz=%d\n";
  static char    form1p[]= "called from %s. malloc failed;"
                                                        " desired long size nz=%d\n";
  static char    prognamp[] = "utiAlloc::lngChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (long *)malloc( nfz*sizeof(long) );
/* fPrintF(stderr, form10p, prognamp, progcallp, (void*)pi); */
    if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, neednz);
  }
  else
  { if(neednz <= *nzp ) return pi;
    p = (long *)realloc(pi, nfz*sizeof(long) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*        ptrAlloc(nz,prognamp)                                               */
/******************************************************************************/
void     *ptrAlloc(size_t nz, char *progcallp)
{ void     *p;
  static    char      prognamp[] = "utiAlloc::ptrAlloc";
  static    char      form1p[] = "called from %s. malloc failed; ptr size nz=%d\n";

  p = (void *)malloc( nz*sizeof(p) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  return p;
}
/******************************************************************************/
/*        ptrChkRealloc(p,nzp,neednz,incrnz,prognamp)                         */
/******************************************************************************/
void     *ptrChkRealloc(void *pi,size_t *nzp,size_t neednz,size_t incrnz, 
                                                                     char *progcallp)
{ void     *p;
  size_t    nfz;
  static char    form2p[] = "called from %s. realloc failed;"
                                     " initial pointer=%p, desired ptr size nz=%d\n";
  static char    form1p[]= "called from %s. malloc failed;"
                                                         " desired ptr size nz=%d\n";
  static char    prognamp[] = "utiAlloc::ptrChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (void *)malloc( nfz*sizeof(p) );
/* fPrintF(stderr, form10p, prognamp, progcallp, (void*)pi); */
    if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, neednz);
  }
  else
  { if(neednz <= *nzp ) return pi;
    p = (void *)realloc(pi, nfz*sizeof(p) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/*        ptrdfAlloc(nz,prognamp)                                             */
/******************************************************************************/
ptrdiff_t  *ptrdfAlloc(size_t nz, char *progcallp)
{ ptrdiff_t  *p;
  static char    form1p[] = "called from %s. malloc failed;"
                                                           " ptrdiff_t size nz=%d\n";
  static char    prognamp[] = "utiAlloc::ptrdfAlloc";

  p = (ptrdiff_t *)malloc( nz*sizeof(ptrdiff_t) );
  if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, nz);
  return p;
}
/******************************************************************************/
/*        ptrdfChkRealloc(p,nzp,neednz,incrnz,prognamp)                       */
/******************************************************************************/
ptrdiff_t  *ptrdfChkRealloc(ptrdiff_t *pi, size_t  *nzp,size_t neednz, size_t incrnz,
                                                                     char *progcallp)
{ ptrdiff_t  *p;
  size_t    nfz;
  static char    form2p[] = "called from %s. realloc failed;"
                                     " initial pointer=%p, desired ptr size nz=%d\n";
  static char    form1p[]= "called from %s. malloc failed;"
                                                         " desired ptr size nz=%d\n";
  static char    prognamp[] = "utiAlloc::ptrdfChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL)
  { p = (ptrdiff_t *)malloc( nfz*sizeof(ptrdiff_t) );
/* fPrintF(stderr, form10p, prognamp, progcallp, (void*)pi); */
    if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, neednz);
  }
  else
  { if(neednz <= *nzp ) return pi;
    p = (ptrdiff_t *)realloc(pi, nfz*sizeof(ptrdiff_t) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/

                    /**************************************/
                    /**         myVec structure          **/
                    /**************************************/
/******************************************************************************/
/*        myVecAlloc(nz,prognamp)                                             */
/******************************************************************************/
myVec    *myVecAlloc(size_t nz, char *progcallp)
{ myVec    *p, *cp;
  int       i;
  static char    form1p[] = "called from %s. malloc failed; myVec size nz=%d\n";
  static char    prognamp[] = "utiAlloc::myVecAlloc";

  cp = p = (myVec *)malloc( nz*sizeof(myVec) );
  if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, nz);
  for(i = nz;  i > 0;  i--, cp++){ cp->z = 0;  cp->x = 0;  cp->p = NULL;}
  return p;
}
/******************************************************************************/
/*        myVecChkRealloc(p,nzp,neednz,incrnz,prognamp)                       */
/******************************************************************************/
myVec    *myVecChkRealloc(myVec *pi, size_t *nzp, size_t neednz, size_t incrnz, 
                                                                     char *progcallp)
{ myVec    *cp, *p;
  size_t    nfz;
  int       i;
  static char    form2p[] = "called from %s. realloc failed;"
                                   " initial pointer=%p, desired myVec size nz=%d\n";
  static char    prognamp[] = "utiAlloc::myVecChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL) p = myVecAlloc(nfz, prognamp);
  else
  { if(neednz <= *nzp) return pi;
    p = (myVec *)realloc(pi, nfz*sizeof(myVec) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
    cp = p + *nzp;
    for(i = (nfz - *nzp);  i > 0;  i--, cp++){ cp->z = 0;  cp->x = 0;  cp->p = NULL;}
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/

                    /**************************************/
                    /**         myBook structure         **/
                    /**************************************/
/******************************************************************************/
/*        myBookAlloc(nz,prognamp)                                            */
/******************************************************************************/
myBook   *myBookAlloc(size_t nz, char *progcallp)
{ myBook   *p, *cp;
  int       i;
  static char    form1p[] = "called from %s. malloc failed; myBook size nz=%d\n";
  static char    prognamp[] = "utiAlloc::myBookAlloc";

  cp = p = (myBook *)malloc( nz*sizeof(myBook) );
  if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, nz);
  for(i = nz;  i > 0;  i--, cp++)
  { cp->bz = 0;  cp->bx = 0;  cp->bp = NULL;
    cp->iz = 0;  cp->ix = 0;  cp->ip = NULL;
  }
  return p;
}
/******************************************************************************/
/*        myBookChkRealloc(p,nzp,neednz,incrnz,prognamp)                       */
/******************************************************************************/
myBook   *myBookChkRealloc(myBook *pi, size_t *nzp, size_t neednz, size_t incrnz, 
                                                                     char *progcallp)
{ myBook   *p, *cp;
  size_t    nfz;
  int       i;
  static char    form2p[] = "called from %s. realloc failed;"
                                  " initial pointer=%p, desired myBook size nz=%d\n";
  static char    prognamp[] = "utiAlloc::myBookChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL) p = myBookAlloc(nfz, prognamp);
  else
  { if(neednz <= *nzp) return pi;
    nfz = neednz+incrnz;
    p = (myBook *)realloc(pi, nfz*sizeof(myBook) );
    if(p == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, pi, neednz);
    cp = p + *nzp;
    for(i =(nfz - *nzp);  i > 0;  i--, cp++)
    { cp->bz = 0;  cp->bx = 0;  cp->bp = NULL;
      cp->iz = 0;  cp->ix = 0;  cp->ip = NULL;
    }
  }
  *nzp = nfz;  return p;
}
/******************************************************************************/
/******************************************************************************/
