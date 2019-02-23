/*  libc/recipe2double/recipe.mnbrak.c                                        */
/** MATRIX : DOUBLE version                                                  **/

#include <math.h>
#include "recipe.nrutil.h"
#include "utiAlloc.h"
#include "recipe.mnbrak.h"

#define GOLD    1.618034
#define GLIMIT  100.0
#define TINY    1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

static    double dmaxarg1, dmaxarg2;
/******************************************************************************/
/*                                                                            */
/* Input : initial points *axp, *bxp                                          */
/*         function func                                                      */
/* Output : *axp, *bxp, *cxp, *fap, *fbp, *fcp                                */
/******************************************************************************/
void      mnbrak(double *axp, double *bxp, double *cxp,
                       double *fap, double *fbp, double *fcp, double (*func)(double))
{ double    ulim,u,r,q,fu,dum;
  double    ax,bx,cx, fa,fb,fc;
  double    bxmu;

  ax = *axp;         bx = *bxp;  
  fa = (*func)(ax);  fb = (*func)(bx);
  if(fb > fa){ SHFT(dum,ax,bx,dum)  SHFT(dum,fb,fa,dum)}
  cx = (bx) + GOLD*(bx-ax);
  fc = (*func)(cx);

  while (fb > fc)
  { r = (bx-ax)*(fb-fc);
    q = (bx-cx)*(fb-fa);
    bxmu = ((bx-cx)*q - (bx-ax)*r) / (2.0*SIGN(DMAX(fabs(q-r),TINY),q-r));
    u = bx - bxmu;
    ulim = bx + GLIMIT*(cx-bx);
    if(bxmu*(u-cx) > 0.0) 
    { fu = (*func)(u);
      if(fu < fc)
      {	*axp = bx;  *bxp = u;   *cxp = cx;
        *fap = fb;  *fbp = fu;  *fcp = fc;
        return;
      }
      else if(fu > fb)
      { *axp = ax;  *bxp = bx;  *cxp = u;
        *fap = fa;  *fbp = fb;  *fcp = fu;
        return;
      }
      u = cx + GOLD*(cx-bx);
      fu = (*func)(u);
    }
    else if((cx-u)*(u-ulim) > 0.0)
    { fu = (*func)(u);
      if(fu < fc)
      { SHFT(bx,cx,u, cx + GOLD*(cx-bx))
        SHFT(fb,fc,fu, (*func)(u))
      }
    }
    else if((u-ulim)*(ulim-cx) >= 0.0)
    { u = ulim;
      fu = (*func)(u);
    }
    else
    { u = cx + GOLD*(cx-bx);
      fu = (*func)(u);
    }
    SHFT(ax,bx,cx,u)
    SHFT(fa,fb,fc,fu)
  }
  *axp = ax;  *bxp = bx;  *cxp = cx;
  *fap = fa;  *fbp = fb;  *fcp = fc;
  return;
}
/******************************************************************************/
/******************************************************************************/
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
/* (C) Copr. 1986-92 Numerical Recipes Software 'i<. */
