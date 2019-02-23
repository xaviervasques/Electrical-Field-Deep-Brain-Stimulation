/*  .../orientHaut/linear4.Common.V1/varPot.linear4.c                         */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20020410                                */

#include  <stddef.h>
#include  <math.h>

#include  "utistdIO.h"
#include  "utiAlloc.h"

#include  "varPot.linear4.h"
#include  "gridZ.linear4.h"
#include  "gridR.linear4.h"
#include  "varGrid.linear4.h"
#include  "eltrdPot.linear4.h"
#include  "cylPot.linear4.h"
#include  "pot_funcBasis.h"
#include  "pot.linear4.h"

static    char    form7p[] = "  The %d nodes values :\n";
static    char    form8p[] = "  For izf = %d to izl = %d \n";

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmVarPotPAllocZero(gmVarPot *p, gmVarGrid *varGridp,
                                               gmEltrdPot *eltrdp, gmCylPot *cylp)
{ size_t    nz;
  static char    prognamp[] = "varPot.linear4::gmVarPotPAllocZero";
  
  p->cylPotp   = cylp;
  p->eltrdPotp = eltrdp;
  p->varGridp  = varGridp;
  nz = varGridp->ndx;
fPrintF(stderr,"%s  nz = %d, varp=%p \n", prognamp, nz, (void*)p->varp);
  p->varp = dblAlloc(nz, prognamp);
  p->ndz = nz;
fPrintF(stderr,"%s  Allocation Done\n", prognamp);
  p->ndx = 0;
  p->ndx = nz;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmVarPotSetSCoef(gmVarPot *p, double **scoefpp, int *sorderp)
{ gmEltrdPot    *eltrdp;

  eltrdp = p->eltrdPotp;
  gmEltrdPotSetSCoef(eltrdp, scoefpp, sorderp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmVarPotSetModes(gmVarPot *p, int neckMod, int *modp)
{ gmEltrdPot    *eltrdp;
  gmCylPot      *cylp;
  gmVarGrid     *vargridp;

  eltrdp = p->eltrdPotp;
  gmEltrdPotSetModes(eltrdp, modp);
  vargridp = p->varGridp;
  cylp = p->cylPotp;
  gmCylPotSetMode(cylp, vargridp->cylcond, neckMod);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmVarPotSetV(gmVarPot *p, double neckv, double *vp)
{ gmEltrdPot    *eltrdp;
  gmCylPot      *cylp;

  eltrdp = p->eltrdPotp;
  gmEltrdPotSetV(eltrdp, vp);
  cylp = p->cylPotp;
  gmCylPotSetV(cylp, neckv);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmVarPotPrint(FILE *bufOut, gmVarPot *p)
{ int       iz, ir, irf, irl, n;
  int      *irfp,  *irlp;
  double   *varp;
  gmGridZ  *gridZp;
  gmVarGrid    *varGridp;
  static    char    form1p[] = "  gmVarPotp=%p, eltrdPotp=%p, cylPotp=%p, varGridp=%p, gridRp=%p, gridZp=%p\n";
  static    char    form2p[] = "          varp=%p, ndz=%d, ndx = %d\n";
  static    char    form9p[] = "  iz=%3d,  irf=%3d, irl=%3d, nf=%4d \n";
  static    char    form10p[] = "%+9.6e ";
  static    char    form11p[] = "varPot.linear4::gmVarPotPrint  BEGIN\n";
  static    char    form12p[] = "varPot.linear4::gmVarPotPrint  END\n";

  fPrintF(bufOut,form11p);
  varGridp = p->varGridp;
  gridZp = varGridp->gridZp;
  fPrintF(bufOut,form1p, p, p->eltrdPotp, p->cylPotp,
                                    p->varGridp, varGridp->gridRp, varGridp->gridZp);
  fPrintF(bufOut,form2p, p->varp, (int)(p->ndz), (int)(p->ndx));
  fPrintF(bufOut,form7p, p->ndx);
  varp = p->varp;
  irfp = varGridp->irnd_firstp;  irlp = varGridp->irnd_lastp;
  n = 0;
  for(iz = varGridp->iznd_first;  iz <= varGridp->iznd_last;  iz++)
  { irf = irfp[iz];  irl = irlp[iz];
    if(irf <= irl)
    { fPrintF(bufOut,form9p, iz, irf, irl, n);
      for(ir = irf;  ir <= irl;  ir++){ fPrintF(bufOut,form10p, *varp++);}
      fPrintF(bufOut,"\n");
      n += irl -irf + 1;
    }
  }
  fPrintF(bufOut,"\n");

  gmEltrdPotPrint(bufOut, p->eltrdPotp);
  fPrintF(bufOut,"\n");
  gmCylPotPrint(bufOut, p->cylPotp);

  fPrintF(bufOut,form12p);
  fflush(bufOut);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmVarPotWrite(FILE *bufWritep, gmVarPot *p)
{ int       iz, ir, irf, irl;
  int      *irfp,  *irlp;
  double   *varp;
  gmGridZ  *gridZp;
  gmVarGrid    *varGridp;
  static    char    form1p[] = "%s\n";
  static    char    form9p[] = "  iz=%3d,  irf=%2d, irl=%2d\n";
  static    char    form10p[] = "%+9.6e ";
  static    char    form11p[] = "%s  BEGIN\n";
  static    char    form12p[] = "%s  END\n";
  static    char    prognamp[]= "varPot.linear4::gmVarPotWrite";

fPrintF(stderr,form11p, prognamp);

  fPrintF(bufWritep,form1p, prognamp);  fflush(bufWritep);
fPrintF(stderr, "varPot.linear4::gmVarPotWrite  TITLE PRINTED\n");
  varGridp = p->varGridp;
  gridZp = varGridp->gridZp;

  fPrintF(bufWritep,form7p, p->ndx);
  fPrintF(bufWritep,form8p, varGridp->iznd_first, varGridp->iznd_last);
  irfp = varGridp->irnd_firstp;  irlp = varGridp->irnd_lastp;

  varp = p->varp;

  for(iz = varGridp->iznd_first;  iz <= varGridp->iznd_last;  iz++)
  { irf = irfp[iz];  irl = irlp[iz];
    if(irf <= irl)
    { fPrintF(bufWritep,form9p, iz, irf, irl);
      for(ir = irf;  ir <= irl;  ir++){ fPrintF(bufWritep,form10p, *varp++);}
      fPrintF(bufWritep,"\n");
    }
  }
  fPrintF(bufWritep,"\n");

  gmEltrdPotWrite(bufWritep, p->eltrdPotp);
  fPrintF(bufWritep,"\n");
  gmCylPotWrite(bufWritep, p->cylPotp);   fflush(bufWritep);
  fflush(bufWritep);
fPrintF(stderr,form12p, prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmVarPotRead(FILE *bufReadp, gmVarPot *p)
{ int       iz, izz, izf, izl;
  int       ir, irf, irl;
  int      *irfp,  *irlp;
  int       ndx;
  double   *varp;
  gmGridZ  *gridZp;
  gmVarGrid    *varGridp;
  gmEltrdPot *eltrdp;
  char      cp[2048];
  static    char    form0p[] = "%le ";
  static    char    form9p[] = "  iz=%3d,  irf=%2d, irl=%2d\n";
  static    char    form11p[] = "%s  BEGIN\n";
  static    char    form12p[] = "%s  END\n";
  static    char    prognamp[]= "varPot.linear4::gmVarPotRead";
 
fPrintF(stderr,form11p, prognamp);

  fscanf(bufReadp,"%s", cp);
fPrintF(stderr," cp=%s \n", cp);
  varGridp = p->varGridp;
  gridZp = varGridp->gridZp;

  fscanf(bufReadp,form7p, &ndx);
fPrintF(stderr,form7p, ndx);
  fscanf(bufReadp,form8p, &izf, &izl);
fPrintF(stderr,form8p, izf,izl);
  varGridp->iznd_first = izf;   varGridp->iznd_last = izl;
  p->ndx = (size_t)ndx;
  varp = p->varp;
  irfp = varGridp->irnd_firstp;  irlp = varGridp->irnd_lastp;
  for(iz = izf;  iz <= izl;  iz++)
  { fscanf(bufReadp,form9p, &izz, &irf, &irl);
/* fPrintF(stderr,form9p, izz,irf,irl); */
    irfp[izz] = irf;   irlp[izz] = irl;
    for(ir = irf;  ir <= irl;  ir++){ fscanf(bufReadp,form0p, varp++);}
/* fPrintF(stderr,"izz=%d, last ir=%d, value=%+9.6e \n",  izz, irl, *(varp-1)); */
    fscanf(bufReadp,"\n");
    if(izz >= izl) break;
  }
  fscanf(bufReadp,"\n");
fPrintF(stderr,"izz=%d, last ir=%d, value=%+9.6e \n",  izz, irl, *(varp-1));

fPrintF(stderr,"\n");
  eltrdp = p->eltrdPotp;
  gmEltrdPotRead(bufReadp, p->eltrdPotp);
  fscanf(bufReadp,"\n");
fPrintF(stderr,"\n");
  gmCylPotRead(bufReadp, p->cylPotp);

fPrintF(stderr,form12p, prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/* assume    field pf->vpp   is already allocated large enough                */
/******************************************************************************/
void      gmVarPotCpy(gmVarPot *pf, gmVarPot *pi)
{ int       i;
  double   *vip, *vfp;
  gmGridZ  *gridZp;
  gmVarGrid     *varGridp;

  vip = pi->varp;     vfp = pf->varp;
  pf->ndx = pi->ndx;
  for(i = 0;  i < pi->ndx;  i++){ *vfp++ = *vip++;}

  varGridp = pi->varGridp;
  if(pf->varGridp == NULL) pf->varGridp = pi->varGridp;
  gridZp = varGridp->gridZp;

  gmEltrdPotCpy(pf->eltrdPotp, pi->eltrdPotp);
  gmCylPotCpy(pf->cylPotp, pi->cylPotp);

  return;
}
/******************************************************************************/
/*                                                                            */
/* assume    field Potp->vpp   is already allocated large enough              */
/******************************************************************************/
void      gmVarP2P(gmPot *Potp, gmVarPot *varPotp)
{ int       iz, izf, izl, izx, iz0, iz1;
  int       ir, irf, irl, irx;
  int      *irfp, *irlp;
  int      *irGeofp, neckIndx;
  int       is, isx, order, *orderp;
  int      *ip;
  double   *scoefp;
  double  **Vpp, *Vp, **vpp, *vp, v, Vinfini, Vltrd;
  double   *varp;
  double   *zp, zs, dz, *zpinEp;
  double   *rhotp;
  int      *izLimFp, *izLimLp;
  gmGridZ  *gridZp;
  gmGridR  *gridRp;
  gmVarGrid    *varGridp;
  gmEltrdPot   *eltrdp;
  gmCylPot     *cylp;
  int       cylcond, neckmode;
  static    char      form1p[] = "%s  z-Grid pointers differ  STOP\n";
  static    char      form2p[] = "%s  r-Grid pointers differ  STOP\n";
  static    char      prognamp[] = "varPot.linear4::gmVarP2P";

  eltrdp = varPotp->eltrdPotp;
  cylp   = varPotp->cylPotp;
  varGridp = varPotp->varGridp;
  gridZp   = varGridp->gridZp;
  gridRp   = varGridp->gridRp;
  if(gridZp != Potp->gridZp){ fPrintF(stderr,form1p,prognamp); exit(1);}
  if(gridRp != Potp->gridRp){ fPrintF(stderr,form2p,prognamp); exit(1);}
  Potp->varPotp = varPotp;

  varp = varPotp->varp;
  Vpp = Potp->vpp;
  izf  = varGridp->iznd_first;          izl  = varGridp->iznd_last;
  irfp = varGridp->irnd_firstp;         irlp = varGridp->irnd_lastp;
  for(iz = izf;  iz <= izl;  iz++)
  { irf = irfp[iz];  irl = irlp[iz];
    Vp = Vpp[iz] + irf;
    for(ir = irf;  ir <= irl;  ir++){ *Vp++ = *varp++;}
  }

                                                                   /** Boundaries **/
  cylcond = cylp->cylcond;
  neckmode = cylp->neckMod;
  Vinfini = cylp->infiniV;
  izx = gridZp->zndx -1;
  irGeofp = gridZp->irGeond_firstp;

  irx = gridRp->rhotndx -1;
  neckIndx = gridRp->neckIndx;
                                                 /** Neck and Upper cylinder face **/
  iz = izx;  Vp = Vpp[iz];
       if(cylcond == 1) irx = gridRp->rhotndx -1;
  else if(cylcond == 0) irx = neckIndx;
  for(ir = 0;  ir <= irx; ir++){ *Vp++ = Vinfini;}
                                                                      /** Large r **/
  if(cylcond == 1)
  { irx = gridRp->rhotndx -1;
    izx = gridZp->zndx -1;
    for(iz = 0;  iz < izx;  iz++){ Vp = Vpp[iz];  Vp[irx] = Vinfini;}
  }
                                                          /** Lower cylinder face **/
  if(cylcond == 1)
  { iz = 0;  Vp = Vpp[iz];
    irx = gridRp->rhotndx -1;
    for(ir = 0;  ir <= irx;  ir++){ *Vp++ = Vinfini;}
  }
                                                                         /** Pins **/
  isx = gridZp->pinX;
  ip = gridZp->pinExtIp;
  for(is = 0;  is < isx;  is++)
  { iz1 = *ip++;    iz0 = *ip++;
    Vltrd = eltrdp->pinVp[is];
    vpp = Vpp + iz0;
    for(iz = iz0 ;  iz <= iz1 ;  iz++){ Vp = *vpp++; *(Vp + irGeofp[iz]) = Vltrd;}
  }
                                                              /** Pin Extremities **/
  izLimFp = varGridp->izLimFp;        
  izLimLp = varGridp->izLimLp;
  zp = gridZp->zndp;
  rhotp = gridRp->rhotndp;
  isx = gridZp->pinExtX;
  zpinEp = gridZp->pinExtZp;
  orderp = eltrdp->sorderp;
  for(is = 0;  is < isx;  is++)
  { zs = *zpinEp++;
    Vltrd = eltrdp->pinVp[is/2];
    order = *orderp++;
    scoefp = (eltrdp->scoefpp)[is];
    iz1 = izLimFp[is];
    iz0 = izLimLp[is];
    vpp = Vpp + iz0;
    for(iz = iz0;  iz <= iz1;  iz++)
    { Vp = *vpp++;
      dz = zp[iz] - zs;
      if(is % 2  == 0) dz = -dz;
      irf = irfp[iz];
      vp = Vp + irGeofp[iz];
      for(ir = irGeofp[iz];  ir < irf;  ir++)
      { *vp++ = Vltrd + angle(scoefp, order, rhotp[ir], dz);
      }
    }
  }
                              /* Electrode interior : set boundary value at same z */
  vpp = Vpp;
  for(iz = 0;  iz <= izx;  iz++, vpp++)
  { irf = irGeofp[iz];
    vp = *vpp;  v = *(vp + irf);
    for(ir = 0;  ir < irf;  ir++){ *vp++ = v;}
  }
  return;
}
/******************************************************************************/
/*                                                                            */
/* assume    field varPotp->varp   is already allocated large enough          */
/******************************************************************************/
void      gmP2VarP(gmVarPot *varPotp, gmPot *Potp)
{ int       iz, izf, izl;
  int       ir, irf, irl;
  int      *irfp, *irlp, *irGeofp;
  int       is, isx, order, izs, *zpinIp;
  double   *scoefp, valuep[2];
  double  **Vpp, *Vp, Vltrd;
  double   *varp;
  double   *zp, zs, dz, *zpinEp;
  double   *rhotp;
  gmGridZ  *gridZp;
  gmGridR  *gridRp;
  gmVarGrid    *varGridp;
  gmEltrdPot   *eltrdp;
  static    char  form1p[] = "%s  z-Grid pointers differ  STOP";
  static    char  form2p[] = "%s  r-Grid pointers differ  STOP";
  static    char  prognamp[] = "varPot.linear4::gmP2VarP";

  eltrdp = varPotp->eltrdPotp;
  varGridp = varPotp->varGridp;
  gridZp   = varGridp->gridZp;
  gridRp   = varGridp->gridRp;
  if(gridZp != Potp->gridZp){ fPrintF(stderr, prognamp,form1p); exit(1);}
  if(gridRp != Potp->gridRp){ fPrintF(stderr, prognamp,form2p); exit(1);}

  varp = varPotp->varp;
  Vpp = Potp->vpp;
  izf  = varGridp->iznd_first;    izl  = varGridp->iznd_last;
  irfp = varGridp->irnd_firstp;   irlp = varGridp->irnd_lastp;
  for(iz = izf;  iz <= izl;  iz++)
  { 
    irf = irfp[iz];  irl = irlp[iz];
    Vp = Vpp[iz] + irf;
    for(ir = irf;  ir <= irl;  ir++){ *varp++ = *Vp++ ;}
  }
                                                               /** Pin Potentials **/
  isx = gridZp->pinX;
  for(is = 0;  is < isx;  is++)
  { iz = gridZp->pinExtIp[2*is];
    eltrdp->pinVp[is] = *(Vpp[iz]);
  }
                                                              /** Pin Extremities **/
  zp = gridZp->zndp;
  rhotp = gridRp->rhotndp;
  irGeofp = gridZp->irGeond_firstp;
  zpinEp = gridZp->pinExtZp;
  zpinIp = gridZp->pinExtIp;
  isx = gridZp->pinExtX;
  for(is = 0;  is < isx;  is++)
  { zs = zpinEp[is];  izs = zpinIp[is];
    ir = irGeofp[izs];
    Vltrd = *Vpp[izs];
    scoefp = eltrdp->scoefpp[is];
    order = eltrdp->sorderp[is];
    valuep[0] = *(Vpp[izs] + ir +1) - Vltrd;
    izf = (is % 2 == 0)?  izs + 1: izs - 1;
    dz = zp[izf] - zs; 
    valuep[1] = *(Vpp[izf] +ir)  - Vltrd;
    anglePerpTOcoefs(scoefp, order, valuep, dz, rhotp[ir +1]);
  }

  return;
}
/******************************************************************************/
/******************************************************************************/
