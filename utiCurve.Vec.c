/*   ../libmy/utiCurve.Vec.c                                                  */
/*   Mennessier Gerard                  20020902                              */
/*   Last revised : M.G.                20020904                              */

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"
#include  "utiCurve.Vec.h"

/******************************************************************************/
/*        myCurveVecAlloc(nz,prognamp)                                        */
/******************************************************************************/
myCurveVec    *myCurveVecAlloc(size_t nz, char *progcallp)
{ myCurveVec    *p, *cp;
  int       i;
  static char    form1p[] = "called from %s. malloc failed; myCurveVec size nz=%d\n";
  static char    prognamp[] = "utiCurve.Vec::myCurveVecAlloc";

  cp = p = (myCurveVec *)malloc( nz*sizeof(myCurveVec) );
  if(p == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, nz);
  for(i = nz;  i > 0;  i--, cp++)
  { cp->z = 0;  cp->x = 0;  cp->p = NULL;
    cp->gcp = NULL;
  }
  return p;
}
/******************************************************************************/
/*        myCurveVecChkRealloc(p,nzp,neednz,incrnz,prognamp)                  */
/******************************************************************************/
myCurveVec    *myCurveVecChkRealloc(myCurveVec *pi, size_t *nzp, size_t neednz,
                                                      size_t incrnz, char *progcallp)
{ myCurveVec    *cp, *p;
  size_t    nfz;
  int       i;
  static char    form2p[] = "called from %s. realloc failed;"
                              " initial pointer=%p, desired myCurveVec size nz=%d\n";
  static char    prognamp[] = "utiCurve.Vec::myCurveVecChkRealloc";

  nfz = neednz + incrnz;
  if(pi == NULL) p = myCurveVecAlloc(nfz, prognamp);
  else
  { if(neednz <= *nzp) return pi;
    p = (myCurveVec *)realloc(pi, nfz*sizeof(myCurveVec) );
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
/******************************************************************************/
