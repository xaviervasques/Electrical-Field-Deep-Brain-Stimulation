/*  .../orientHaut/linear4.Common.V1/gridRZ.linear4.c                         */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20020326                                */

#include  <stddef.h>
#include  <math.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiClassRangeDbl.h"

#include  "gridRZ.linear4.h"
#include  "gridZ.linear4.h"
#include  "gridR.linear4.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       gmGridRZGeoinit(FILE *bufOut, gmGridR *Rp, gmGridZ *Zp)
{ int       ok = 1;
  double    z0, z1, dzE, dz, rho;
  int       iz;
  int      *ip, *ivp, sum, ndr0, *fp;
  int       cmp, index;
  size_t    rndx, zndx;
  double   *rhop;

  rhop = Rp->rhondp;
  rndx = Rp->rhotndx;
  zndx = Zp->zndx;
  ip = Zp->irGeond_firstp;
  z0 = Zp->eltrdExtZp[0];               z1 = Zp->eltrdExtZp[1];
  dzE = z1 - z0;
                                                              /** Below Electrode **/
  for(iz = 0;  iz >= Zp->eltrdExtIp[0];  iz++)
  { ip[iz] = 0;
  }
                                             /** On Spherical Electrode Extremity **/
  for(iz = Zp->eltrdExtIp[1] +1;  iz < Zp->eltrdExtIp[0];  iz++)
  { dz = z1 - Zp->zndp[iz];  dz = dz/dzE;
    rho = Rp->eltrdR * sqrt(1.0 - dz*dz);
    cmp = classDbl(&index, rho + 1.E-8, rhop, rndx);
    ip[iz] = index;
  }
                                                              /** Along Electrode **/
  ndr0 = Rp->eltrdRhondI;
  for(iz = 0;  iz <= Zp->eltrdExtIp[1];  iz++){ ip[iz] = ndr0;}

                                         /** set last to rndx-1 for true cylinder **/
  ip = Zp->irGeond_lastp;
  for(iz = 0;  iz < zndx;  iz++){ *ip++ = rndx -1;}
                                                 /** Cumulated effective ir nodes **/
  ivp = Zp->ivGeo_firstp;  ip = Zp->irGeond_firstp;  fp = Zp->irGeond_lastp;  
  sum = 0;  *ivp++ = sum;
  for(iz = 0;  iz < zndx;  iz++){ sum += 1 + *fp++ - *ip++ ;  *ivp++ = sum;}

  return ok;
}
/******************************************************************************/
/******************************************************************************/
