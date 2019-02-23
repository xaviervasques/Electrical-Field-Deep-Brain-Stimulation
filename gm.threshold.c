/*  URMAE/orientHaut/linear4.GL.V4/gm.threshold.c                             */
/*  Mennessier Gerard                 20030509                                */
/*  Last Revised : G.M.               20040210                                */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utiMath.constant.def.h"

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiVecChr.h"

#include  "gm.threshold.h"

#include  "gridZ.linear4.h"
#include  "gridR.linear4.h"
#include  "pot.linear4.h"
#include  "gm.linear4.initsolve.glob.h"                  /** for current E=grad V **/

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double    gmThresholdGetVolGlob(double vthr)
{ double     vol, voltot, voldelta;
  size_t    zllx;
  double   *zllp, *zndp, zll, z0, z1;
  int       iz;
  double   *rhotllp, *rhondp;
  int      *irGeoip, *irGeofp;
  int       ir;
  double  **vpp, *vp0, *vp1, vmean;
  short     boolPrint = 0;
  static char    form1p[] = "%s z0=%f, z1=%f, r0=%f, r1=%f, vmean=%f \n";
  static char    prognamp[] = "gm.threshold::gmThresholdGetVolGlob";


  zllx = gridZi_mm.zllx;
  zllp = gridZi_mm.zllp;
  zndp = gridZi_mm.zndp;
  rhotllp =gridRi_mm.rhotllp;
  rhondp =gridRi_mm.rhondp;
  irGeoip = gridZi_mm.irGeond_firstp;
  irGeofp = gridZi_mm.irGeond_lastp;
  vpp = gradVf.vpp;
  vp1 = *vpp++;

  voltot = 0.0;
  for(iz = 0;  iz < zllx;  iz++)
  { vol = 0.0;
    zll = zllp[iz];
    z0 = zndp[iz];  z1 = zndp[iz +1];
    vp0 = vp1;
    vp1 = *vpp++;
    for(ir = irGeoip[iz];  ir < irGeofp[iz] -1;  ir++)
    { 
      vmean = 0.25*( *(vp0 + ir) + *(vp0 + ir +1) + *(vp1 + ir) + *(vp1 + ir +1) );
      if(vmean > vthr)
      { voldelta = zll * rhotllp[ir] * 0.5*(rhondp[ir] + rhondp[ir+1]);
        vol += 2. * myPI * voldelta;
        if(boolPrint) fPrintF(stderr, form1p, prognamp, z0, z1, rhondp[ir], rhondp[ir+1], vmean);
      }
    }
    voltot += vol;
  }
  return voltot;
}
/******************************************************************************/
/******************************************************************************/
