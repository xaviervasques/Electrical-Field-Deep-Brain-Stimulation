/*  libc/recipe2double/recipe.brent.c                                         */
/** MATRIX : DOUBLE version                                                  **/

#include <math.h>
#include "recipe.nrutil.h"
#include "utiAlloc.h"
#include "recipe.brent.h"

#define ITMAX  100
#define CGOLD  0.3819660
#define ZEPS   1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

/******************************************************************************/
/*                                                                            */
/* Input : initial points *axp, *bxp                                          */
/*         function func                                                      */
/* Output : *axp, *bxp, *cxp, *fap, *fbp, *fcp                                */
/******************************************************************************/
double    brent(double ax, double bx, double cx, double (*f)(double),
                                                            double tol, double *xmin)
{ int       iter;
  double    a,b, p,q,r, tol1,tol2, xm;
  double    u,v,w,x, fu,fv,fw,fx;
  double    d = 0.0, e = 0.0, etemp;

  a = (ax < cx ? ax : cx);
  b = (ax > cx ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = (*f)(x);
  for(iter = 1;  iter <= ITMAX;  iter++)
  { xm = 0.5*(a+b);
    tol1 = tol*fabs(x) + ZEPS;
    tol2 = 2.0*tol1;
    if(fabs(x-xm) <= (tol2 - 0.5*(b-a)))
    { *xmin = x;
      return fx;
    }
    if(fabs(e) > tol1)
    { r = (x-w)*(fx-fv);
      q = (x-v)*(fx-fw);
      p = (x-v)*q-(x-w)*r;
      q = 2.0*(q-r);
      if(q > 0.0)  p = -p;
      q = fabs(q);
      etemp = e;
      e = d;
      if(fabs(p) >= fabs(0.5*q*etemp)  ||  p <= q*(a-x)  ||  p >= q*(b-x))
      { e= (x >= xm ? a-x : b-x);
        d = CGOLD*e;
      }
      else
      { d = p/q;
        u = x+d;
        if(u-a < tol2  ||  b-u < tol2) d = SIGN(tol1,xm-x);
      }
    }
    else
    { e=(x >= xm ? a-x : b-x);
      d = CGOLD*e;
    }
    u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
    fu = (*f)(u);
    if(fu <= fx)
    { if(u >= x) a=x;
      else       b=x;
      SHFT(v,w,x,u) SHFT(fv,fw,fx,fu)
    } 
    else
    { if(u < x) a=u;
      else      b=u;
      if(fu <= fw  ||  w == x)
      { v = w;
        w = u;
        fv = fw;
        fw = fu;
      }
      else if(fu <= fv  ||  v == x  ||  v == w)
      { v = u;
        fv = fu;
      }
    }
  }
  *xmin = x;
  return fx;
}
/******************************************************************************/
/******************************************************************************/
#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
/* (C) Copr. 1986-92 Numerical Recipes Software 'i<. */
