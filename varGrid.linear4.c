/*  .../orientHaut/linear4.Common.V1/varGrid.linear4.c                        */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20020327                                */

#include  <stddef.h>
#include  <math.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

#include  "gridZ.linear4.h"
#include  "gridR.linear4.h"
#include  "varGrid.linear4.h"
#include  "eltrd.def.h"
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmVarGridPAlloc(size_t ndz, gmVarGrid *p)
{
  p->irnd_firstp = intAlloc(ndz,"varGrid.linear4::gmVarGridPAlloc  irnd_firstp");
  p->irnd_lastp  = intAlloc(ndz,"varGrid.linear4::gmVarGridPAlloc  irnd_lastp");
  p->iv_firstp   = intAlloc(ndz +1,"varGrid.linear4::gmVarGridPAlloc  iv_firstp");
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       gmVarGridInit(gmVarGrid *p, gmGridR *gridRp, gmGridZ *gridZp, int *orderp)
{ size_t    zndx, rndx;
  int       is, isx;
  int       iz, izx, izn;
  int       iznd_first, iznd_last;
  int       varVndx;
  int      *irfp, *irlp, irf, ir, irx, nr, sum;
  int      *irGeop, *ivf;
  int      *pinEip;
  int      *izLimFp, *izLimLp, izLimF, izLimL;
  int       cylcond;
  double   *rLimsp,  *r2Limsp;
  double   *zp, *rp, *pinEzp;
  double    z, z2, zs, dz, dz2, rhot2;
  
  p->gridZp = gridZp;
  izLimFp = p->izLimFp;                 izLimLp = p->izLimLp;
  p->gridRp = gridRp;
  rLimsp = p->rLimsp;                   r2Limsp = p->r2Limsp;
  cylcond = gridRp->cylcond;

  isx = gridZp->pinExtX;
  pinEip = gridZp->pinExtIp;            pinEzp = gridZp->pinExtZp;
  zndx = gridZp->zndx;
  zp = gridZp->zndp;
  irGeop = gridZp->irGeond_firstp;

  rndx = gridRp->rhotndz;
  rp = gridRp->rhotndp;

  for(is = 0;  is < isx;  is++)
  { int       ns, order;
    double    rhot, rcoef;
    order = orderp[is];
    if(order < 0){ rLimsp[is] = 0.0;  continue;}
    if(order == 0) rcoef = 1.01;
    else           rcoef = 1.51;
    ns = gridZp->pinExtIp[is];
    if(is % 2 == 0){ rhot = (gridZp->zllp)[ns];}
    else           { rhot = (gridZp->zllp)[ns -1];}
/**    rhot = gridZp->insL / 10.0 ;                  **/
    rLimsp[is] = rhot * rcoef;
  }

  for(is = 0;  is < isx;  is++){ r2Limsp[is] = rLimsp[is]*rLimsp[is];}
  p->cylcond = cylcond;

                                           /** Limiting iz around pin extremities **/
  for(is = 0;  is < isx;  is++){ izLimFp[is] = izLimLp[is] = 0;}
  izx = zndx-1;
  for(is = 0;  is < isx;  is++)
  { zs = pinEzp[is];
    izLimF = pinEip[is];
    for(iz = izx;  iz >= izLimF;  iz--)
    { dz = zp[iz] - zs;
      if(dz*dz <= r2Limsp[is] + 1.E-10){ izLimF = iz;  break;}
    }
    izLimFp[is] = izLimF;

    izLimL = pinEip[is];
    for(iz = izLimL;  iz >= 0;  iz--)
    { dz = zp[iz] - zs;
      if(dz*dz >  r2Limsp[is] + 1.E-10){ izLimL = iz+1;  break;}
    }
    izLimLp[is] = izLimL;
    izx = izLimL;
  }

                       /** determination of first and last rho index for pot.var. **/
  gmVarGridPAlloc(zndx, p);
  irfp = p->irnd_firstp;                irlp = p->irnd_lastp;
  varVndx = 0;

  if(cylcond > 0){ izn = 1;  izx = zndx -2;  irx = rndx -2;}
  else           { izn = 0;  izx = zndx -1;  irx = rndx -1;}
  
  iznd_first = zndx;  iznd_last = 0;
  irfp[0] = rndx;                       irlp[0] = 0;
  irfp[zndx -1] = rndx;                 irlp[zndx -1] = 0;

                                            /**  low z part : electrode insulator **/
  for(iz = izn;  iz < izLimLp[isx -1];  iz++)
  { irfp[iz] = irGeop[iz];  irlp[iz] = irx;
    nr = irlp[iz] - irfp[iz] + 1;
    if(nr > 0)
    { if(iznd_first > iz) iznd_first = iz;
      if(iznd_last  < iz) iznd_last = iz;
      varVndx += nr;
    }
  }
  izn = izLimLp[isx -1];
  for(is = isx-1;  is > 0;  is -= 2)
  {                                                       /** Below Pin Extremity **/
    for(iz = izn;  iz < izLimLp[is];  iz++)
    { irfp[iz] = irGeop[iz];  irlp[iz] = irx;
      nr = irlp[iz] - irfp[iz] + 1;
      if(nr > 0)      
      { if(iznd_first > iz) iznd_first = iz;
        if(iznd_last  < iz) iznd_last = iz;
        varVndx += nr;
      }
    } 
                                                   /** Around Lower Pin Extremity **/
    for(iz = izLimLp[is];  iz <= izLimFp[is];  iz++)
    { z = zp[iz];  z2 = z * z;  dz2 = (z - pinEzp[is])*(z - pinEzp[is]);
      irf = irGeop[iz];
      for(ir = irGeop[iz];  ir <= irx;  ir++)
      { rhot2 = rp[ir] * rp[ir] + dz2;
        if(rhot2 <= r2Limsp[is]) irf = ir;
        else break;
      }
      irfp[iz] = irf + 1;
                                 /** determination of last rho index for pot.var. **/
      irlp[iz] = irx;
      nr = irlp[iz] - irfp[iz] + 1;
      if(nr > 0)
      { if(iznd_first > iz) iznd_first = iz;
        if(iznd_last  < iz) iznd_last = iz;
        varVndx += nr;
      }
    }
                                                                    /** Along Pin **/
    for(iz = izLimFp[is] +1;  iz < izLimLp[is-1];  iz++)
    { irfp[iz] = irGeop[iz] + 1;
                                 /** determination of last rho index for pot.var. **/
      irlp[iz] = irx;
      nr = irlp[iz] - irfp[iz] + 1;
      if(nr > 0)
      { if(iznd_first > iz) iznd_first = iz;
        if(iznd_last  < iz) iznd_last = iz;
        varVndx += nr;
      }
    }
                                                   /** Around Upper Pin Extremity **/
    for(iz = izLimLp[is-1];  iz <= izLimFp[is-1];  iz++)
    { z = zp[iz];  z2 = z * z;  dz2 = (z - pinEzp[is-1])*(z - pinEzp[is-1]);
      irf = irGeop[iz];
      for(ir = irGeop[iz];  ir <= irx;  ir++)
      { rhot2 = rp[ir] * rp[ir] + dz2;
        if(rhot2 <= r2Limsp[is]) irf = ir;
        else break;
      }
      irfp[iz] = irf + 1;
                                 /** determination of last rho index for pot.var. **/
      irlp[iz] = irx;
      nr = irlp[iz] - irfp[iz] + 1;
      if(nr > 0)
      { if(iznd_first > iz) iznd_first = iz;
        if(iznd_last  < iz) iznd_last = iz;
        varVndx += nr;
      }
    }
    izn = izLimFp[is-1] +1;
  }
                                 /**  high z part : electrode extremity and above **/
  for(iz = izn;  iz < izx;  iz++)
  { irfp[iz] = irGeop[iz];  irlp[iz] = irx;
    nr = irlp[iz] - irfp[iz] + 1;
    if(nr > 0)
    { if(iznd_first > iz) iznd_first = iz;
      if(iznd_last  < iz) iznd_last = iz;
      varVndx += nr;
    }
  }
                      /** case iz = izx (z = zmax) : neck and upper cylinder face **/
  iz = izx;
  irfp[iz] = rndx;  irlp[iz] = 0;
       if(cylcond == 0){ irfp[iz] = gridRp->neckIndx +1;  irlp[iz] = irx;}
  else if(cylcond <  0){ irfp[iz] = 0;                    irlp[iz] = irx;}
  nr = irlp[iz] - irfp[iz] +1;
  if(nr > 0)
  { if(iznd_first > iz) iznd_first = iz;
    if(iznd_last  < iz) iznd_last = iz;
    varVndx += nr;
  }

  p->iznd_first = iznd_first;
  p->iznd_last  = iznd_last;
  p->ndx = (size_t)varVndx;

  ivf = p->iv_firstp;
  sum = 0;  *ivf++ = sum;
  for(iz = 0;  iz <= izx;  iz++)
  { nr = irlp[iz] - irfp[iz] + 1;
    if(nr > 0) sum += nr;
    *ivf++ = sum;
  }
  return  varVndx;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmVarGridPrint(FILE *bufOut, gmVarGrid *p)
{ int       is, isx;
  int       iz, izx, izf, izl;
  int       n, nr;
  int      *ip;
  gmGridR  *gridRp;
  gmGridZ  *gridZp;
  static    char    prognamp[] = "varGrid.linear4::gmVarGridPrint";
  static    char    form1p[] = "%s  rho Unit = %s, z Unit = %s\n";
  static    char    form2p[] = "  Limiting radius for the %d singular points\n";
  static    char    form3p[] = "  First effective iz nodes = %d,"
                               "  Last effective iz nodes = %d, out of %d z nodes\n";
  static    char    form4p[] = "  List of the First ir variable nodes =\n";
  static    char    form5p[] = "  List of the Last  ir variable nodes =\n";
  static    char    form6p[] = "  Total number of variable nodes = %d, expected %d\n";
  static    char    form7p[] = "  List of the %d index of first varPot per line =\n";
  static    char    form8p[] = "  CYL CONDUCTANCE = %+d  \n";

  gridRp = p->gridRp;  gridZp = p->gridZp;
  isx = gridZp->pinExtX;
  fPrintF(bufOut,form1p, prognamp, gridRp->unitnamep, gridZp->unitnamep);
  fPrintF(bufOut,form8p, gridRp->cylcond);
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,form2p, isx);
  for(is = 0;  is < isx;  is++)
  { fPrintF(bufOut,"  pin index %d,  rLim=%f,  r2Lim=%f \n",
                                              is, (p->rLimsp)[is], (p->r2Limsp)[is]);
    fPrintF(bufOut,"    izLimF=%d, izLimL=%d \n", p->izLimFp[is], p->izLimLp[is]);
  }
  fPrintF(bufOut,"\n");  

  izf = p->iznd_first;  izl = p->iznd_last;
  fPrintF(bufOut,form3p, izf,izl, gridZp->zndx);
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,form4p);
  izl = 0;
  izx = gridZp->zndx -1;
  ip = p->irnd_firstp;
  for(iz = 0, n = 0;  iz <= izx;  iz++, n++)
  { if(n != 0  &&  n%20 == 0) fPrintF(bufOut,"\n");
    fPrintF(bufOut," %2d,", *ip++);
  }
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,form5p);
  ip = p->irnd_lastp;
  for(iz = 0, n = 0;  iz <= izx;  iz++, n++)
  { if(n != 0  &&  n%20 == 0) fPrintF(bufOut,"\n");
    fPrintF(bufOut," %2d,", *ip++);
  }
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,form7p, gridZp->zndx +1);
  ip = p->iv_firstp;
  for(iz = 0, n = 0;  iz <= gridZp->zndx;  iz++, n++)
  { if(n != 0  &&  n%20 == 0) fPrintF(bufOut,"\n");
    fPrintF(bufOut," %4d,", *ip++);
  }
  fPrintF(bufOut,"\n");
                                                                        /** check **/
  for(iz = p->iznd_first, n = 0;  iz <= p->iznd_last;  iz++)
  { nr = p->irnd_lastp[iz] - p->irnd_firstp[iz] +1;
    if(nr > 0) n += nr;
  }
  fPrintF(bufOut,form6p, n, p->ndx);

  fflush(bufOut);
  return;
}
/******************************************************************************/
/******************************************************************************/
