/*  ../libmy/utiTnsr.c                                                        */
/*  Mennessier Gerard                   20010611                              */
/*       revised : M.G.                 20030205                              */
/*       revised : M.G.                 20050706                              */
/*  Last revised : M.G.                 20051124                              */

#include  <stddef.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"
#include  "utiTnsr.h"

/******************************************************************************/
/*        tnsr2dblAlloc(nz,prognamp)                                          */
/*                                                                            */
/******************************************************************************/
tnsr2dbl *tnsr2dblAlloc(size_t nz, char *progcallp)
{ tnsr2dbl *p;
  static    char      prognamp[] = "utiTnsr::tnsr2dblAlloc";
  static    char      form1p[] = "called from %s. malloc failed; ptr size nz=%d\n";

  p = (tnsr2dbl *)malloc( nz*sizeof(tnsr2dbl) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  return p;
}
/******************************************************************************/
/*        tnsr2dblPAlloc(tnsr2dbl *p, size_t ev1z, size_t ev2z)               */
/******************************************************************************/
void      tnsr2dblPAlloc(tnsr2dbl *tp, size_t ev1z, size_t ev2z)
{ int       i;
  size_t    z;
  double  **pp, *p;
  static    char      prognamp[] = "utiTnsr::tnsr2dblPAlloc";

  z = ev1z * ev2z;
  p = dblAlloc(z, prognamp);
  tp->p = p;
  tp->pz = z;
  pp = (double **)ptrAlloc(ev1z, prognamp);
  tp->pp = pp;  tp->ev1z = ev1z;  tp->ev2z = ev2z;
  for(i = 0;  i < ev1z;  i++, pp++){ *pp = p;  p+= ev2z;}
  tp->ev1x = 0;  tp->ev2x = 0;          /** effective dimensions of Vector Spaces **/
  return;
}
/******************************************************************************/
/*        tnsr2dblPRealloc(tnsr2dbl *p, size_t ev1z, size_t ev2z)             */
/******************************************************************************/
void      tnsr2dblPRealloc(tnsr2dbl *tp, size_t needev1z, size_t needev2z)
{ int       i;
  size_t    z;
  double  **pp, *p;
  static    char      prognamp[] = "utiTnsr::tnsr2dblPAlloc";

  z = needev1z * needev2z;
  p = dblChkRealloc(tp->p, &(tp->pz), z, 0, prognamp);
  tp->p = p;
  
  pp = (double **)ptrChkRealloc(tp->pp, &(tp->ev1z), needev1z, 0, prognamp);
  tp->pp = pp;

  tp->ev2z = needev2z;
  for(i = 0;  i < needev1z;  i++, pp++){ *pp = p;  p+= needev2z;}
  for(i = needev1z;  i < tp->ev1z;  i++, pp++){ *pp = NULL;}

  tp->ev1x = 0;  tp->ev2x = 0;          /** effective dimensions of Vector Spaces **/
  return;
}
/******************************************************************************/
/*        tnsr2dblPrint(tnsr2dbl *p, size_t ev1z, size_t ev2z)                */
/******************************************************************************/
void      tnsr2dblPrint(FILE *bufp, tnsr2dbl *p)
{ static    char      form1p[] = "%s  printing %p\n";
  static    char      prognamp[] = "utiTnsr::tnsr2dblPrint";

  fPrintF(bufp,form1p, prognamp, p);
  fPrintF(bufp, "p allocated at %p, total size=%d\n", (void*)p->p, p->pz);
  fPrintF(bufp, "pp allocated at %p, ev1z=%d, ev1x=%d, ev2z=%d, ev2x=%d\n",
                                   (void*)p->pp, p->ev1z, p->ev1x, p->ev2z, p->ev2x);
  return;
}
/******************************************************************************/
/*        tnsr2dblPPrint(FILE *bufp, tnsr2dbl *p)                             */
/******************************************************************************/
void      tnsr2dblPPrint(FILE *bufp, tnsr2dbl *p)
{ int       i1, i2, i2x;
  double  **dpp, *dp;
  static    char      form1p[] = "%s  printing %p.p\n";
  static    char      prognamp[] = "utiTnsr::tnsr2dblPPrint";

  fPrintF(bufp,form1p, prognamp, p);
  dpp = p->pp;
  i2x = p->ev2x;
  for(i1 = 0;  i1 < p->ev1x;  i1++)
  { dp = *dpp++;
    fPrintF(bufp,"%5d: ", i1);
    for(i2 = 1;  i2 <= i2x;  i2++)
    { fPrintF(bufp,"%+11.4e ", *dp++);
      if(i2%10 == 0) fPrintF(bufp,"\n");
    }
    fPrintF(bufp,"\n");
  }
  return;
}
/******************************************************************************/
/*        tnsr2dblCopy(tnsr2dbl *pf, tnsr2dbl *pi)                            */
/******************************************************************************/
void      tnsr2dblCopy(tnsr2dbl *pf, tnsr2dbl *pi)
{ int       i1, i2, i2x;
  double  **dipp, *dip, **dfpp, *dfp;

  if(pf->p == NULL || pf->pp == NULL) tnsr2dblPAlloc(pf, pi->ev1x, pi->ev2x);
  dipp = pi->pp;  dfpp = pf->pp;
  i2x = pi->ev2x;
  for(i1 = 0;  i1 < pi->ev1x;  i1++)
  { dip = *dipp++;  dfp = *dfpp++;
    for(i2 = 1;  i2 <= i2x;  i2++){ *dfp++ = *dip++;}
  }
  pf->ev1x = pi->ev1x;  pf->ev2x = pi->ev2x;
  return;
}
/******************************************************************************/
/*        tnsr2dblAdd(tnsr2dbl *pi, tnsr2dbl *p)                              */
/*                                                                            */
/* add p into pi, result in pi                                                */
/* p and pi must have same values fot ev1x, ev2x                              */
/******************************************************************************/
void      tnsr2dblAdd(tnsr2dbl *pi, tnsr2dbl *p)
{ int       i1, i2, i2x;
  double  **dipp, *dip, **dpp, *dp;
  static    char      form1p[] = "%s  ERROR: NOT same dimensions, pi (%d,%d), p (%d,%d)\n";
  static    char      prognamp[] = "utiTnsr::tnsr2dblAdd";

  if(pi->ev1x != p->ev1x  ||  pi->ev2x != p->ev2x)
  { fPrintF(stderr, form1p, prognamp, pi->ev1x, pi->ev2x, p->ev1x, p->ev2x);
  }
  dipp = pi->pp;  dpp = p->pp;
  i2x = pi->ev2x;
  for(i1 = 0;  i1 < pi->ev1x;  i1++)
  { dip = *dipp++;  dp = *dpp++;
    for(i2 = 1;  i2 <= i2x;  i2++, dip++, dp++){ *dip += *dp;}
  }
  return;
}
/******************************************************************************/
/*        tnsr2dblShift(tnsr2dbl *pi, double d)                               */
/*                                                                            */
/* add constant  d  to all values in pi, result in pi                         */
/******************************************************************************/
void      tnsr2dblShift(tnsr2dbl *pi, double d)
{ int       i1, i2, i2x;
  double  **dipp, *dip;

  dipp = pi->pp;
  i2x = pi->ev2x;
  for(i1 = 0;  i1 < pi->ev1x;  i1++)
  { dip = *dipp++;
    for(i2 = 1;  i2 <= i2x;  i2++, dip++){ *dip += d;}
  }
  return;
}
/******************************************************************************/
/*        tnsr2dblShiftE1a(tnsr2dbl *pi, double *shiftp)                      */
/*                                                                            */
/* shift all "lines" by "line" shifp                                          */
/*                                                                            */
/* i.e.     pi[i1][i2] -> pi[i1][i2] +  shifp[i2]                             */
/*                                                                            */
/* shiftp MUST point to ALREADY allocated and defined doubles,                */
/*                                  at least pi->ev2x  double values          */
/* Method a                                                                   */
/******************************************************************************/
void      tnsr2dblShiftE1a(tnsr2dbl *pi, double *shiftp)
{ int       i1, i2, i2x;
  double  **dipp, *dip;
  double    d;

  i2x = pi->ev2x;

  for(i2 = 0;  i2 < i2x;  i2++)
  { d = *shiftp++;
    dipp = pi->pp;
    for(i1 = 0;  i1 < pi->ev1x;  i1++, dipp++)
    { dip = i2 + *dipp;
      *dip += d;
    }
  }
  return;
}
/******************************************************************************/
/*        tnsr2dblShiftE1b(tnsr2dbl *pi, double *shiftp)                      */
/*                                                                            */
/* shift all "lines" by "line" shifp                                          */
/*                                                                            */
/* i.e.     pi[i1][i2] -> pi[i1][i2] +  shifp[i2]                             */
/*                                                                            */
/* shiftp MUST point to ALREADY allocated and defined doubles,                */
/*                                  at least pi->ev2x  double values          */
/* Method b                                                                   */
/******************************************************************************/
void      tnsr2dblShiftE1b(tnsr2dbl *pi, double *shiftp)
{ int       i1, i2, i2x;
  double  **dipp, *dip;
  double   *dp;

  i2x = pi->ev2x;

  dipp = pi->pp;
  for(i1 = 0;  i1 < pi->ev1x;  i1++, dipp++)
  { dp = shiftp;
    dip = *dipp;
    for(i2 = 0;  i2 < i2x;  i2++, dip++){ *dip += *dp++;}
  }
  return;
}
/******************************************************************************/
/*        tnsr2dblSumE1a(tnsr2dbl *pi, double *sump)                          */
/*                                                                            */
/* add all "lines" into line sump                                             */
/*                                                                            */
/* i.e.     sump[i2] = SUM(over i1) pi[i1][i2]                                */
/*                                                                            */
/* sump MUST point to already allocated, at least pi->ev2x  doubles           */
/* Method a                                                                   */
/******************************************************************************/
void      tnsr2dblSumE1a(tnsr2dbl *pi, double *sump)
{ int       i1, i2, i2x;
  double  **dipp, *dip;
  double   *dp;

  i2x = pi->ev2x;
  dp = sump;
  for(i2 = 0;  i2 < i2x;  i2++){ *dp++ = 0.0;}

  dp = sump;
  for(i2 = 0;  i2 < i2x;  i2++, dp++)
  { dipp = pi->pp;
    for(i1 = 0;  i1 < pi->ev1x;  i1++){ dip = *dipp++;  *dp += *(dip + i2);}
  }
  return;
}
/******************************************************************************/
/*        tnsr2dblSumE1b(tnsr2dbl *pi, double *sump)                          */
/*                                                                            */
/* add all "lines" into line sump                                             */
/*                                                                            */
/* i.e.     sump[i2] = SUM(over i1) pi[i1][i2]                                */
/*                                                                            */
/* sump MUST point to ALREADY allocated, at least pi->ev2x  doubles           */
/* Method b                                                                   */
/******************************************************************************/
void      tnsr2dblSumE1b(tnsr2dbl *pi, double *sump)
{ int       i1, i2, i2x;
  double  **dipp, *dip;
  double   *dp;

  i2x = pi->ev2x;
  dp = sump;
  for(i2 = 0;  i2 < i2x;  i2++){ *dp++ = 0.0;}

  dipp = pi->pp;
  for(i1 = 0;  i1 < pi->ev1x;  i1++, dipp++)
  { dp = sump;
    dip = *dipp;
    for(i2 = 0;  i2 < i2x;  i2++, dp++){ *dp += *dip++;}
  }
  return;
}
/******************************************************************************/
/*        tnsr2dblMinMaxE1b(tnsr2dbl *pi, double *minp, double maxp)          */
/*                                                                            */
/* get min and max over all "lines"  into lines minp, maxp                    */
/*                                                                            */
/* i.e.      minp[i2] = MIN(over i1) pi[i1][i2]                               */
/*           maxp[i2] = MAX(over i1) pi[i1][i2]                               */
/* minp, maxp MUST point to ALREADY allocated, at least pi->ev2x  doubles     */
/* Method b                                                                   */
/******************************************************************************/
void      tnsr2dblMinMaxE1b(tnsr2dbl *pi, double *minp, double *maxp)
{int       i1, i2, i2x;
  double  **dipp, *dip, d;
  double   *dminp, *dmaxp;

  i2x = pi->ev2x;
  dipp = pi->pp;
  dip = *dipp;
  dminp = minp;  dmaxp = maxp;
  for(i2 = 0;  i2 < i2x;  i2++){ *dminp++ = *dmaxp++ = *dip++;}

  dipp++;
  for(i1 = 1;  i1 < pi->ev1x;  i1++, dipp++)
  { dip = *dipp;
    dminp = minp;  dmaxp = maxp;
    for(i2 = 0;  i2 < i2x;  i2++, dminp++, dmaxp++)
    { d = *dip++;
           if(d > *dmaxp) *dmaxp = d;
      else if(d < *dminp) *dminp = d;
    }
  }
  return;
}
/******************************************************************************/


/******************************************************************************/
/*        tnsr2fltAlloc(nz,prognamp)                                          */
/*                                                                            */
/******************************************************************************/
tnsr2flt *tnsr2fltAlloc(size_t nz, char *progcallp)
{ tnsr2flt *p;
  static    char      prognamp[] = "utiTnsr::tnsr2fltAlloc";
  static    char      form1p[] = "called from %s. malloc failed; ptr size nz=%d\n";

  p = (tnsr2flt *)malloc( nz*sizeof(tnsr2flt) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  return p;
}
/******************************************************************************/
/*        tnsr2fltPAlloc(tnsr2flt *p, size_t ev1z, size_t ev2z)               */
/******************************************************************************/
void      tnsr2fltPAlloc(tnsr2flt *tp, size_t ev1z, size_t ev2z)
{ int       i;
  size_t    z;
  float   **pp, *p;
  static    char      prognamp[] = "utiTnsr::tnsr2fltPAlloc";

  z = ev1z * ev2z;
  p = fltAlloc(z, prognamp);
  tp->p = p;
  tp->pz = z;
  pp = (float  **)ptrAlloc(ev1z, prognamp);
  tp->pp = pp;  tp->ev1z = ev1z;  tp->ev2z = ev2z;
  for(i = 0;  i < ev1z;  i++, pp++){ *pp = p;  p+= ev2z;}
  tp->ev1x = 0;  tp->ev2x = 0;          /** effective dimensions of Vector Spaces **/
  return;
}
/******************************************************************************/
/*        tnsr2fltPrint(tnsr2flt *p, size_t ev1z, size_t ev2z)                */
/******************************************************************************/
void      tnsr2fltPrint(FILE *bufp, tnsr2flt *p)
{ static    char      prognamp[] = "utiTnsr::tnsr2fltPrint";
  static    char      form1p[] = "%s  printing %p\n";

  fPrintF(bufp,form1p, prognamp, p);
  fPrintF(bufp, "p allocated at %p, total size=%d\n", (void*)p->p, p->pz);
  fPrintF(bufp, "pp allocated at %p, ev1z=%d, ev1x=%d, ev2z=%d, ev2x=%d\n",
                                   (void*)p->pp, p->ev1z, p->ev1x, p->ev2z, p->ev2x);
  return;
}
/******************************************************************************/
/*        tnsr2fltPPrint(FILE *bufp, tnsr2flt *p)                             */
/******************************************************************************/
void      tnsr2fltPPrint(FILE *bufp, tnsr2flt *p)
{ int       i1, i2, i2x;
  float   **dpp, *dp;
  static    char      prognamp[] = "utiTnsr::tnsr2fltPPrint";
  static    char      form1p[] = "%s  printing %p.p\n";

  fPrintF(bufp,form1p, prognamp, p);
  dpp = p->pp;
  i2x = p->ev2x;
  for(i1 = 0;  i1 < p->ev1x;  i1++)
  { dp = *dpp++;
    fPrintF(bufp,"%5d: ", i1);
    for(i2 = 1;  i2 <= i2x;  i2++)
    { fPrintF(bufp,"%+11.4f ", *dp++);
      if(i2%10 == 0) fPrintF(bufp,"\n");
    }
    fPrintF(bufp,"\n");
  }
  return;
}
/******************************************************************************/



/******************************************************************************/
/*        tnsr2intAlloc(nz,prognamp)                                          */
/******************************************************************************/
tnsr2int *tnsr2intAlloc(size_t nz, char *progcallp)
{ tnsr2int *p;
  static    char      prognamp[] = "utiTnsr::tnsr2intAlloc";
  static    char      form1p[] = "called from %s. malloc failed; ptr size nz=%d\n";

  p = (tnsr2int *)malloc( nz*sizeof(tnsr2int) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  return p;
}
/******************************************************************************/
/*        tnsr2intPAlloc(tnsr2int *p, size_t ev1z, size_t ev2z)               */
/******************************************************************************/
void      tnsr2intPAlloc(tnsr2int *tp, size_t ev1z, size_t ev2z)
{ int       i;
  size_t    z;
  int     **pp, *p;
  static    char      prognamp[] = "utiTnsr::tnsr2intPAlloc";

  z = ev1z * ev2z;
  p = intAlloc(z, prognamp);
  tp->p = p;
  tp->pz = z;
  pp = (int **)ptrAlloc(ev1z, prognamp);
  tp->pp = pp;  tp->ev1z = ev1z;  tp->ev2z = ev2z;
  for(i = 0;  i < ev1z;  i++, pp++){ *pp = p;  p+= ev2z;}
  tp->ev1x = 0;  tp->ev2x = 0;          /** effective dimensions of Vector Spaces **/
  return;
}
/******************************************************************************/
/*        tnsr2intPrint(tnsr2int *p, size_t ev1z, size_t ev2z)                */
/******************************************************************************/
void      tnsr2intPrint(FILE  *bufp, tnsr2int *p)
{ static    char      prognamp[] = "utiTnsr::tnsr2intPrint";
  static    char      form1p[] = "%s  printing %p\n";

  fPrintF(bufp,form1p, prognamp, p);
  fPrintF(bufp, "p allocated at %p, total size=%d\n", (void*)p->p, p->pz);
  fPrintF(bufp, "pp allocated at %p, ev1z=%d, ev1x=%d, ev2z=%d, ev2x=%d\n",
                                   (void*)p->pp, p->ev1z, p->ev1x, p->ev2z, p->ev2x);
  return;
}
/******************************************************************************/
/*        tnsr2intPPrint(FILE *bufp, tnsr2int *p)                             */
/******************************************************************************/
void      tnsr2intPPrint(FILE *bufp, tnsr2int *p)
{ int       i1, i2, i2x;
  int     **ipp, *ip;
  static    char      prognamp[] = "utiTnsr::tnsr2intPPrint";
  static    char      form1p[] = "%s  printing %p.p\n";

  fPrintF(bufp,form1p, prognamp, p);
  ipp = p->pp;
  i2x = p->ev2x;
  for(i1 = 0;  i1 < p->ev1x;  i1++)
  { ip = *ipp++;
    fPrintF(bufp,"%5d: ", i1);
    for(i2 = 1;  i2 <= i2x;  i2++)
    { fPrintF(bufp,"%+7d ", *ip++);
      if(i2%15 == 0) fPrintF(bufp,"\n");
    }
    fPrintF(bufp,"\n");
  }
  return;
}
/******************************************************************************/
/******************************************************************************/
