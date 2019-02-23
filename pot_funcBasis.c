/*  bio_phys/URMAE/orientHaut/linear4.Common/pot_funcBasis.c                  */
/*  Mennessier Gerard                 20000415                                */
/*  Last Revised : G.M.               20040210                                */

#include  <stddef.h>
#include  <math.h>

#include  "utiMath.constant.def.h"

#include  "pot_funcBasis.h"
#include  "pot_funcBasis.def.h"

/******************************************************************************/
/*                                                                            */
/* Given  zeff = +- (z - zangle), algebraic z-ordinate from singular point    */
/*   and  rht                     rhot-ordinate from pin                      */
/* compute the Potential (normalized to 0 at singular point) at this point    */
/*   from coefs series                                                        */
/*                                                                            */
/* a -> alpha = 1/2;                                                          */
/* rta = rt^alpha;                                                            */
/* sa0   = sin(alpha psi),      ca0 = cos(alpha psi)                          */
/* sap1  = sin( (alpha+1) psi), cap1 = cos( (alpha+1) psi)                    */
/* sap2  = sin( (alpha+2) psi), cap2 = cos( (alpha+2) psi)                    */
/* sam1  = sin( (alpha-1) psi), cam1 = cos( (alpha-1) psi)                    */
/* sam2  = sin( (alpha-2) psi), cam2 = cos( (alpha-2) psi)                    */
/*                                                                            */
/* WARNING : zeff = +- (z - zangle);                                          */
/*                  + if Conductor above Insulator                            */
/*                  - if Conductor below Insulator                            */
/*  i.e.     zeff = z - znICa;                                                */
/*        or zeff = -(z - znH)                                                */
/*  so that Conductor for psi = 0,  Insulator for psi = Pi                    */
/*                                                                            */
/* par[0] = S00;  par[1] = S11;  par[2] = S22;  par[3] = S33;                 */
/******************************************************************************/
double    angle(double par[NCOEF_],int order,double rht,double zeff)
{ double    rt, rt2, rta, logrt;
  double    psi, s1,c1,c1p1,c1m1, c11sq;
  double    sa0, ca0;
  double    sap1= 0.0, cap1= 0.0, sam1= 0.0, cam1= 0.0,
            sap2= 0.0, cap2= 0.0, sam2= 0.0, cam2= 0.0,
            sap3= 0.0, cap3= 0.0, sam3= 0.0, cam3= 0.0;
  double    g;
  double    g00= 0., g10= 0., g11= 0., g20= 0., g21= 0., g22= 0.,
                                                  g30= 0., g31= 0., g32= 0., g33= 0.;
  double    s00= 0.0, s11= 0.0, s10= 0.0,
            s22= 0.0, s21= 0.0, s20= 0.0, s33= 0.0, s32= 0.0, s31= 0.0, s30= 0.0;
  double    s200p0= 0.0, c200p0= 0.0;

  if(order < 0) return 0.0;

  rt2 = rht*rht + zeff*zeff;
  rt = sqrt(rt2);   rta = sqrt(rt);
  logrt = log(rt);
                           /** rht = 0. **/
  if(rt2 == 0.0)
  { c1 = 1.0;   s1 = 0.0;
    ca0 = 1.0;  sa0 = 0.0;  psi = 0.0;
    logrt = 0.;
    return 0.0;
  }

                                                                 /** GENERIC rht  **/
  c1 = zeff/rt;  s1 = rht/rt;  psi = atan2(s1,c1);
                                           /** sa0^2 = 2(1-c1)/4 = s1^2/(2(1+c1))
                                               ca0^2 = 2(1+c1)/4 = s1^2/(2(1-c1)) **/
  if(zeff >= 0.0)
  { c1p1 = 1.0 + c1;  c11sq = sqrt(2.0*c1p1);
    ca0 = 0.5*c11sq;  sa0 = s1/c11sq;
  }
  else
  { c1m1 = 1.0 - c1;  c11sq = sqrt(2.0*c1m1);
    ca0 = s1/c11sq;  sa0 = 0.5*c11sq;
  }

  if(order == 0) goto ETIQ1 ;
  sap1 = sa0*c1 + ca0*s1;    sam1 = sa0*c1 - ca0*s1;
  cap1 = ca0*c1 - sa0*s1;    cam1 = ca0*c1 + sa0*s1;
  if(order == 1) goto ETIQ1 ;

  sap2 = sap1*c1 + cap1*s1;  sam2 = sap1*c1 - cap1*s1;
  cap2 = cap1*c1 - sap1*s1;  cam2 = cap1*c1 + sap1*s1;
  if(order == 2) goto ETIQ1 ;

  sap3 = sap2*c1 + cap2*s1;  sam3 = sap2*c1 - cap2*s1;
  cap3 = cap2*c1 - sap2*s1;  cam3 = cap2*c1 + sap2*s1;
ETIQ1:

  s00 = par[0];
    g00 = s00 * sa0;
  if(order == 0) goto ETIQ2 ;

  s11 = -s00 * 0.5/(1.5 * myPI);
    g11 = s11 * sap1;
  s10 = par[1];
    g10 = s11 * psi*cap1 + s10 * sap1 + s00 * 0.25*(cap1 - cam1);
  if(order == 1) goto ETIQ2 ;

  s22 = - 0.25/(2.5 * myPI) *s11;
    g22 = s22 * sap2;
  s21 = - 2.0/(2.5) *s22 - 0.5/(2.5 * myPI) *s10;
    g21 = s22 * 2.0*psi*cap2 + s21 * sap2 + s11 * 0.25*(cap2 - ca0);
  s22 = par[2];
  s200p0 = s00 /16. *3.5/1.5 ;  c200p0 = - s10 * 0.25;
    g20 = - s22 * psi*psi * sap2 - s11 * 0.25*psi*sap2 + s21 * psi*cap2
          + s20 * sap2 + s11 * 0.25*psi*sa0
          + s200p0 * sa0 + c200p0 * (ca0 - cap2) - s00 *3./32. *sam2;
  if(order == 2) goto ETIQ2 ;
/*
  s33 = - 1.0/6.0/(3.5 * myPI) *s20;
  s31 = -3.0/(3.5) * s33 -1.0/(3.5 * myPI) *c32 + 0.25/(3.5 * myPI) *s22
                                                            - 0.25/(3.5 * myPI) *s21;
*/
  s33 = s32 = s31 =0.0;
    g33 = 0.0;
    g32 = 0.0;
    g31 = 0.0;
  s30 = par[3];
    g30 = 0.0;
ETIQ2:

  g = 0.0;
  if(order >= 3) g += ((g33*logrt + g32)*logrt + g31)*logrt + g30;
  if(order >= 2) g += g*rt + (g22*logrt + g21)*logrt + g20;
  if(order >= 1) g += g*rt + g11*logrt + g10;
  g = g*rt + g00;

  return  rta*g;
}
/******************************************************************************/
/*                                                                            */
/* valuesp[0] = value at (drhot!=0, dz= 0.0 )                                 */
/* valuesp[1] = value at (drhot =0, dz=-|dz|)                                 */
/*                                                                            */
/* compute the first 2 coef from 2 values                                     */
/*    at (z=zs, rhot = drhot) and (z=zs-dz, rhot =0, i.e. along insulator)    */
/*    ( value is 0 at singular point)                                         */
/*                                                                            */
/******************************************************************************/
void      anglePerpTOcoefs(double parp[2], int order, double valuesp[2],
                                                          double dz, double drhot)
{ double    basisp[2][2];
  double    par[NCOEF_] = {0.0, 0.0, 0.0};
  double    det;

  dz = fabs(dz);  drhot = fabs(drhot);
  switch (order)
  { case (0):
    { par[0] = 1.0;
      basisp[0][0] = angle(par, 0,drhot, 0.0);
      basisp[0][1] = angle(par, 0,0.0  , -dz);
      parp[0] = 0.5 * (valuesp[0] / basisp[0][0] + valuesp[1] / basisp[0][1]);
      parp[1] = 0.0;
      break;
    }
    case (1):
    { par[0] = 1.0;  par[1] = 0.0;
      basisp[0][0] = angle(par, 1,drhot, 0.0);
      basisp[0][1] = angle(par, 1,0.0 , -dz);
      par[0] = 0.0;  par[1] = 1.0;
      basisp[1][0] = angle(par, 1,drhot, 0.0);
      basisp[1][1] = angle(par, 1,0.0 , -dz);
      det = basisp[0][0] * basisp[1][1] - basisp[0][1] * basisp[1][0];
      parp[0] = (valuesp[0] * basisp[1][1] - valuesp[1] * basisp[1][0]) / det;
      parp[1] = (valuesp[1] * basisp[0][0] - valuesp[0] * basisp[0][1]) / det;
      break;
    }
    default : parp[0] = parp[1] = 0.0;
  }
  return;
}
/******************************************************************************/
/*                                                                            */
/* Contribution of interval [zs, zs + deltaz] of the electrode (rhot = 0.0)   */
/******************************************************************************/
double    angleTOintensity0(double par[NCOEF_],int order, double deltaz)
{ double    Is = 0.0;
  double    rt, rta;                                  /** rta = rt^alpha = rt^1/2 **/
  double    s00 = 0.0, s10 = 0.0, s11 = 0.0;
  
  if(order < 0)     return 0.0;
  if(deltaz <= 0.0) return 0.0;

  s00 = par[0];
  if(order == 0) goto ETIQ1 ;

  s11 = -s00 * 0.5/(1.5 * myPI);
  s10 = par[1];
  if(order == 1) goto ETIQ1 ;
ETIQ1:

  rt = deltaz;                          rta = sqrt(rt);
  Is = s00;
  if(order == 0) goto ETIQ2 ;

  Is += ( s11 * log(rt) + s10) * rt;
  if(order == 0) goto ETIQ2 ;
ETIQ2:

  return  (rta * Is);
}
/******************************************************************************/
/*                                                                            */
/* Contribution of rectangle      [zs - deltaz, zs + deltaz] * [0, rhot]      */
/* to the action (= chi2)                                                     */
/******************************************************************************/
double    angleTOaction(double par[NCOEF_],int order, double deltaz, double rhot)
{ double    chi2 = 0.0;
  double    rt, rt2;
  double    spsi0,cpsi0;
  double    lps, lpc;
  double    s00 = 0.0, s10 = 0.0, s11 = 0.0;

  if(order < 0) return 0.0;

  if(deltaz < 0.0) deltaz = -deltaz;
  rt2 = deltaz*deltaz + rhot*rhot;      rt = sqrt(rt2);
  spsi0 = rhot / rt;                    cpsi0 = deltaz / rt;
  lps = log( (1. + spsi0)/cpsi0 );      lpc = log( (1. + cpsi0)/spsi0 );

  s00 = par[0];
  chi2 = s00 * 0.5 *(deltaz * lps + rhot * lpc);
  if(order == 0) goto ETIQ1 ;

  s11 = - s00/( 3. * myPI);
  s10 = par[1];
  chi2 += - 3. /8. * myPI *( deltaz*deltaz *(1./cpsi0 - 1.) + rhot*rhot * lpc);
  if(order == 1) goto ETIQ1 ;
ETIQ1:

  chi2 = chi2 * s00;
  return chi2;
}
/******************************************************************************/
/*                                                                            */
/* V = par[0]/r * (1 + par[1] *(h/r)^2 (3 cos^2 -1)/2 + O(h/r ^4) )           */
/* where cos = z/r                                                            */
/* WARNING :                                                                  */
/* r = rho = Reltrd + rht                                                     */
/******************************************************************************/
double    mono_quadru(double par[NINFINITYX_],double rho,double z)
{ double    v, r2, rinv, r2inv, z2;
  
  z2 = z*z;
  r2 = rho*rho + z2;
  r2inv = 1.0/r2;
  rinv = sqrt(r2inv);
  v = 1. + r2inv * par[1] *(1.5*z2*r2inv - 0.5);
  return  par[0]*rinv*v;
}
/******************************************************************************/
/******************************************************************************/
