/*  URMAE/orientHaut/linear4.GL.V4/grid_pot.ST.linear4.c                      */
/*  Mennessier Gerard                 20011228                                */
/*  Last Revised : G.M.               20040210                                */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utiMath.constant.def.h"

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"
#include  "utiClassRangeDbl.h"

#include  "grid_pot.ST.linear4.h"

#include  "gm.linear4.initsolve.glob.h"  /** for Vf, gradVf, gridZi_mm, gridRi_mm **/
#include  "pallidus.geom.glob.h"

#include  "gridZ.linear4.h"
#include  "gridR.linear4.h"
#include  "pot.linear4.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
gmGridPotSTXY *gmGridPotSTXYAlloc(size_t nz, char *progcallp)
{ gmGridPotSTXY *p;
  static    char      form1p[] = "called from %s. malloc failed;"
                                                       " gmGridPotSTXY size nz=%d\n";
  static    char      prognamp[] = "grid_pot.ST.linear4::gmGridPotSTXYAlloc";

  p = (gmGridPotSTXY *)malloc( nz*sizeof(gmGridPotSTXY) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  p->unitnamep = NULL;
  p->alphap = NULL;
  p->alphaz = 0;
  p->saxp = NULL;
  p->saxz = 0;
  return p;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridPotSTXYPAlloc(size_t az, size_t saxz, gmGridPotSTXY *p)
{
  p->alphap = dblAlloc(az, "grid_pot.ST.linear4::gmGridPotSTXYarPAlloc alphap");
  p->alphaz = az;  p->alphax = 0;
  p->saxp = dblAlloc(saxz, "grid_pot.ST.linear4::gmGridPotSTXYarPAlloc saxp");
  p->saxz = saxz;  p->saxx = 0;
  p->potp = dblAlloc(saxz*az, "grid_pot.ST.linear4::gmGridPotSTXYarPAlloc potp");
  p->potz = saxz*az;  p->potx = 0;
  p->potpp = ptrAlloc(saxz, "grid_pot.ST.linear4::gmGridPotSTXYarPAlloc potpp");
  p->potpz = saxz;  p->potpx = 0;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridPotSTXYsetAlpha(gmGridPotSTXY *p)
{ double    defstep = 2.0;                    /** default angle step is 2 degrees **/
  double    dstep, *dp;
  size_t    ax;
  int       i, ahx;
  static    char      prognamp[] = "grid_pot.ST.linear4::gmGridPotSTXYsetAlpha";
   
  ahx = 180.0/defstep;  ax = 2 * ahx +1;
  dstep = myPI /(double)ahx;
fPrintF(stderr, "     %s  dstep=%f, ahx=%d, ax=%d\n", prognamp, dstep, ahx, ax);
  p->alphap = dblChkRealloc(p->alphap, &(p->alphaz), ax, 1, prognamp);
  dp = p->alphap;
  for(i = -ahx;  i <= ahx;  i++){ *dp++ = i * dstep;}
  p->alphax = ax;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridPotSTXYsetSaxis(size_t rx, double *rhop, gmGridPotSTXY *p)
{ int       i;
  double   *dp;
  static    char      prognamp[] = "grid_pot.ST.linear4::gmGridPotSTXYsetSaxis";

  p->saxp = dblChkRealloc(p->saxp, &(p->saxz), rx, 0, prognamp);
  dp = p->saxp;
  for(i = 0;  i < rx;  i++){ *dp++ = *rhop++;}
  p->saxx = rx;
  return;
}
/******************************************************************************/
/*                                                                            */
/*  compute the values potpp[isax][ialpha]  at the nodes of the rectangle     */
/*          from values stored at fpp[iz][ir]                                 */
/******************************************************************************/
void      gmGridPotSTXYsetPot(double  **fpp, double zst, gmGridPotSTXY *p)
{ int       is, ia;
  double  **dpp, *dp;
  double    saxis, zel, alphar, *saxp, *alphap;
  size_t    ax, sax, potz;
  static    char      prognamp[] = "grid_pot.ST.linear4::gmGridPotSTXYsetPot";

  p->zst = zst;
  ax = p->alphax;  sax = p->saxx;
  potz = ax * sax;
  p->potp = dblChkRealloc(p->potp, &(p->potz), potz, 0, prognamp);
  p->potx = 0;
  p->potpp = ptrChkRealloc(p->potpp, &(p->potpz), sax, 0, prognamp);
  dpp = p->potpp;  dp = p->potp;
  for(is = 0;  is < sax;  is++){ *dpp++ = dp;  dp += ax;}
  p->potpx = sax;

  if(fpp == NULL)
  { for(dp = p->potp, is = 0;  is < potz;  is++) *dp++ = 0.0;
      p->potx = potz;  return;
  }

  for(dpp = p->potpp, saxp = p->saxp, is = 0;  is < sax;  is++)
  { saxis = *saxp++;
    dp = *dpp++;
    for(alphap = p->alphap, ia = 0;  ia < ax;  ia++)
    { alphar = *alphap++; 
      zel = (zst + saxis * cos(alphar) * stheta)/ctheta;
      *dp++ = gmPotInterpol(saxis, zel, fpp, &gridRi_mm, &gridZi_mm);
    }
  }
  p->potx = potz;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridPotSTXYPrint(FILE *bufp, gmGridPotSTXY *p)
{ int       i, j;
  double   *dp, **dpp;
  static    char      prognamp[] = "grid_pot.ST.linear4::gmGridPotSTXYPrint";

  fPrintF(bufp,"%s  p=%p \n", prognamp, (void*)p);
  if(p->unitnamep) fPrintF(bufp,"  unit = %s\n", p->unitnamep);
  fPrintF(bufp,"  alphap = %p, alphaz = %d, alphax = %d\n", 
                                             (void*)p->alphap, p->alphaz, p->alphax);
  for(dp = p->alphap, i=0;  i < p->alphax;  i++) fPrintF(bufp,"%f ", *dp++);
  fPrintF(bufp,"\n");

  fPrintF(bufp,"  saxp = %p, saxz = %d, saxx = %d\n",
                                                   (void*)p->saxp, p->saxz, p->saxx);
  for(dp = p->saxp, i=0;  i < p->saxx;  i++) fPrintF(bufp,"%f ", *dp++);
  fPrintF(bufp,"\n");

  fPrintF(bufp,"  potp = %p, potz = %d, potx = %d\n",
                                                   (void*)p->potp, p->potz, p->potx);

  fPrintF(bufp,"  potpp = %p, potpz = %d, potpx = %d\n",
                                                (void*)p->potpp, p->potpz, p->potpx);
  for(dpp = p->potpp,  i=0;  i < p->potpx;  i++, dpp++)
  { fPrintF(bufp,"  i = %d, potpp[i]= %p\n", i, (void*)dpp);
    for(dp = *dpp,  j=0;  j < p->alphax;  j++) fPrintF(bufp,"%f ", *dp++);
    fPrintF(bufp,"\n");
  }

  return;
}
/******************************************************************************/
/*                                                                            */
/*  Given director angle alpha, small axis saxis and z ST value               */
/*                 return the ST frame X and Y coordinates                    */
/******************************************************************************/
void      gmGridPotSTalphaSaxisZ2xy(double *xyp,
                                             double alphar, double saxis, double zst)
{ double    x, y;

                /** first compute the x,y coordinates in the (theta, phi=0) frame **/
                /**               from director angle alphar (radian),            **/
                /**               small axis saxis, and large axis saxis/ctheta   **/
  x = (zst * stheta + saxis * cos(alphar))/ctheta;
  y = saxis * sin(alphar) ;
                                                     /** then rotate by phi angle **/
  xyp[0] = x * cphi - y * sphi;
  xyp[1] = x * sphi + y * cphi;
  return;
}
/******************************************************************************/
/*                                                                            */
/*  Given director angle alpha, small axis saxis and z ST value               */
/*                 return the E frame X, Y and Z coordinates                  */
/******************************************************************************/
void      gmGridPotSTalphaSaxisZ2Exyz(double *xyzp,
                                             double alphar, double saxis, double zst)
{ double    rho, ca, sa;

  rho = saxis;
  ca = cos(alphar);  sa = sin(alphar);
  xyzp[0] = rho * ca;
  xyzp[1] = rho * sa;
  xyzp[2] = (zst + saxis * ca * stheta)/ctheta;
  return;
}
/******************************************************************************/
/******************************************************************************/

