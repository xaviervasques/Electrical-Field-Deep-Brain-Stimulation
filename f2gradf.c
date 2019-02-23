/*  bio_phys/URMAE/numerical/linear4/f2gradf.c                                */
/*  Mennessier Gerard                 20010528                                */
/*  Last Revised : G.M.               20010528                                */

#include  <stddef.h>
#include  <math.h>

#include  "utistdIO.h"

/******************************************************************************/
/*                                                                            */
/* Given rectangular grid size xx * yx (stored in xp, yp)                     */
/*       values of function stored in fpp                                     */
/* compute the norm of gradient at each NODE                                  */
/*                                                                            */
/* ASSUME  fp and gNp are allocated to  xx * yx at least,                     */
/*         fpp and gNpp points to begining of each y=Cte , x variable         */
/******************************************************************************/
void      getGeo2DGradNorm(double *xp, size_t xx, double *yp, size_t yx,
                       int *iyGeo2ixfp, int *iyGeo2ixlp, double **fpp, double **gNpp)
{ int       ix, ixf, ixl, iy;
  double   *fp1p, *f0p, *fm1p, *gNp, f1, f2;
  double    gdef = 0.0, gr, gr2;
  double    dxp1, dxm1, dyp1, dym1;

  gNp = gNpp[0];
  for(ix = 0;  ix < xx;  ix++){ *gNp++ = gdef;}
  gNp = gNpp[yx-1];
  for(ix = 0;  ix < xx;  ix++){ *gNp++ = gdef;}

  for(iy = 1;  iy < yx-1;  iy++)
  { 
    ixf = iyGeo2ixfp[iy];  ixl = iyGeo2ixlp[iy];
    fp1p = fpp[iy+1];  f0p = fpp[iy];  fm1p = fpp[iy-1];
    dyp1 = yp[iy+1] - yp[iy];  dym1 = yp[iy-1] - yp[iy];
    gNp = gNpp[iy];
    for(ix = 0;      ix < ixf; ix++){ *gNp++ = gdef;}
    gNp = gNpp[iy] + ixl+1;
    for(ix = ixl+1;  ix < xx;  ix++){ *gNp++ = gdef;}
    gNp = gNpp[iy] + ixf;
    for(ix = ixf;  ix <= ixl; ix++)
    { dxp1 = xp[ix+1] - xp[ix];  dxm1 = xp[ix-1] - xp[ix];
      if(ix == ixf){ gr = (f0p[ix +1] - f0p[ix]) /dxp1;}
      else if(ix == ixl){ gr = (f0p[ix -1] - f0p[ix]) /dxm1;}
      else
      { f1 = f0p[ix +1] - f0p[ix -1];  f2 = f0p[ix +1] + f0p[ix -1] - 2.0 * f0p[ix];
        gr = f1 / (dxp1 - dxm1) * ( -(dxp1*dxp1+dxm1*dxm1)/(dxp1*dxm1) )
                                                       + f2 *(dxp1+dxm1)/(dxp1*dxm1);
        gr /= 2.0;
      }
      gr2 = gr*gr;
      f1 = fp1p[ix] - fm1p[ix];  f2 = fp1p[ix] + fm1p[ix] - 2.0 * f0p[ix];
      gr = f1 / (dyp1 - dym1) * ( -(dyp1*dyp1+dym1*dym1)/(dyp1*dym1) )
                                                       + f2 *(dyp1+dym1)/(dyp1*dym1);
      gr2 += gr*gr;
      *gNp++ = sqrt(gr2);
    }
  }
  return;
}
/******************************************************************************/
/******************************************************************************/
