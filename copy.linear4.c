/*  bio_phys/URMAE/numerical/linear4/copy.linear4.c                           */
/*  Mennessier Gerard                 20010528                                */
/*  Last Revised : G.M.               20010601                                */

#include  <stddef.h>
#include  <math.h>

#include  "utistdIO.h"

#include  "eltrd.def.h"
#include  "cylPot.linear4.h"
#include  "eltrdPot.linear4.h"
#include  "pot.linear4.h"
#include  "varPot.linear4.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void     copyLinearCyl(double *weightp, int weightx, double dV,
                                                 gmCylPot *cylPp, gmCylPot *cylPbasp)
{ int      i;
  double   infiniV, neckI = 0.0;
  double  *wp;

  wp = weightp;
  infiniV = dV;
  for(i = 0;  i < weightx;  i++)
  { infiniV += *wp   * (cylPbasp+i)->infiniV;
    neckI   += *wp++ * (cylPbasp+i)->neckI;
  }
  cylPp->infiniV = infiniV;  cylPp->neckI = neckI;
gmCylPotPrint(stderr, cylPp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void     copyLinearEltrd(double *weightp, int weightx, double dV,
                                         gmEltrdPot *eltrdPp, gmEltrdPot *eltrdPbasp)
{ int      i, is, j;
  double   power = 0.0, effI = 0.0;
  double  *wp, w;
  gmEltrdPot    *eltrdPip;

                                                                    /** Pure Copy **/
  eltrdPip = eltrdPbasp;
  eltrdPp->scoefx = eltrdPip->scoefx;
  for(is = 0;  is < NPEX_ ;  is++){ eltrdPp->sorderp[is] = eltrdPip->sorderp[is];}

                                                                  /** Set to zero **/
  for(is = 0;  is < NPX_ ;  is++)
  { eltrdPp->pinVp[is] = dV;  eltrdPp->pinIp[is] = 0.0;
  }
  for(j = 0;  j < NPEX_ *NCOEF_ ;  j++){ eltrdPp->scoefmatp[j] = 0.0;}

                                                                       /** Sum up **/
  wp = weightp;
  for(i = 0;  i < weightx;  i++, eltrdPip++)
  { w = *wp++;
    power += w * eltrdPip->power;
    for(is = 0;  is < NPX_ ;  is++)
    { eltrdPp->pinVp[is] += w * eltrdPip->pinVp[is];
      eltrdPp->pinIp[is] += w * eltrdPip->pinIp[is];
    }
    for(j = 0;  j < NPEX_ *NCOEF_ ;  j++)
    { eltrdPp->scoefmatp[j] += w * eltrdPip->scoefmatp[j];
    }
  }
  eltrdPp->power = power;
  for(is = 0;  is < NPX_ ;  is++){ effI += fabs(eltrdPp->pinIp[is]);}
  eltrdPp->effI = effI;
gmEltrdPotPrint(stderr, eltrdPp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void     copyLinearPot(double *weightp, int weightx, double dV,
                                                             gmPot *Vp, gmPot *Vbasp)
{ int      i, j;
  double  *vp, *vip;
  double  *wp, w;
  gmPot   *Vip;

  vp = Vp->vp;
  for(j = 0;  j < Vp->vpx;  j++){ *vp++ = dV;}

  Vip = Vbasp;
  wp = weightp;
  for(i = 0;  i < weightx;  i++, Vip++)
  { w = *wp++;
    vip = Vip->vp;    vp = Vp->vp;
    for(j = 0;  j < Vp->vpx;  j++){ *vp++ += w * (*vip++);}
  }
  return;
}
/******************************************************************************/
/******************************************************************************/
