/*  URMAE/orientHaut/linear4.Common.V1/gridZ.linear4.c                        */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20020404                                */

#include  <stddef.h>
#include  <math.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"
#include  "utiClassRangeDbl.h"

#include  "gridZ.linear4.delta.def.h"
#include  "gridZ.linear4.h"
#include  "cyl.def.h"
#include  "eltrd.def.h"


/******************************************************************************/
/*                                                                            */
/******************************************************************************/
gmGridZ  *gmGridZAlloc(size_t nz, char *progcallp)
{ gmGridZ  *p;
  static    char      prognamp[] = "gridZ.linear4::gmGridZAlloc";
  static    char      form1p[] = "called from %s. malloc failed; gmGridZ size nz=%d\n";

  p = (gmGridZ *)malloc( nz*sizeof(gmGridZ) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  p->unitnamep = NULL;
  p->zndp = p->zllp = NULL;
  p->zndz = p->zndx = 0;
  p->zllz = p->zllx = 0;
  p->eltrdR = 0.0;
  p->zndMIN = p->zndMAX = 0.0;
  p->insL = p->cndL = 0.0;
  p->eltrdExtCyl = p->eltrdExtEll = 0.0;
  p->pinX = p->pinExtX = 0;
  p->symPlZ = 0.0;
  p->symPlI = 0;
  p->eltrdExtX = 0;
  p->irGeond_firstp = p->irGeond_lastp = p->ivGeo_firstp = 0;
  return p;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridZPAlloc(size_t llz, gmGridZ *p)
{ size_t    ndz;

  ndz = llz + 1;
  p->zndp = dblAlloc(ndz, "gridZ.linear4::gmGridZPAlloc  zndp");
  p->zndz = ndz;
  p->zllp = dblAlloc(llz, "gridZ.linear4::gmGridZPAlloc  zllp");
  p->zllz = llz;
  p->irGeond_firstp = intAlloc(ndz, "gridZ.linear4::gmGridZPAlloc irGeond_firstp");
  p->irGeond_lastp  = intAlloc(ndz, "gridZ.linear4::gmGridZPAlloc irGeond_lastp");
  p->ivGeo_firstp   = intAlloc(ndz +1, "gridZ.linear4::gmGridZPAlloc ivGeo_firstp");

  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       gmGridZinit(FILE *bufOut, gmGridZ *p, double *infozpap, int *infozpaxp)
{ int       ok = 1;
  double    eltrdR_mm;
  double    insL_mm, cndL_mm, eltrdExtCyl_mm, eltrdExtEll_mm, zndMIN_mm, zndMAX_mm;
  double    symPlZ_mm;
  int       pinX = NPX_, pinExtX = NPEX_;
  double   *pinExtZp, *eltrdExtZp;
  int       is;
  double    zc;
  double    deltazc, deltazcic, deltazci, deltazie;
  double    deltazp[MYDELTAZN_ ],  zpap[MYDELTAZN_ ], delta;
  int       ndeltazp[MYDELTAZN_ ], idz;
  int       iz, izt;
  int       zllhx;
  size_t    zllx, zndx;
  double   *zndp, *zllp;
  int       cmp, index;
  int      *ip;
  static    char    prognamp[]  = "gridZ.linear4::gmGridZinit";
  static    char    form1p[] = "%s; "
                " InsulatorPinLength =%f, ConductorPinLength =%f, \n"
                " EllipsoidalExtremityLength =%f, CylindricalExtremityLength =%f, \n"
                                                    " zndMIN =%f, zndMAX =%f (%s)\n";
  static    char    form3p[] = "    deltaz index %2d, %2d intv length %f, tot=%f\n";
  static    char    form4p[] = "    TOTAL z intv = %d, TOTAL z nodes = %d\n";
  static    char    form5p[] = "    izt = %d,  expected = %d\n";
  static    char    form6p[] = "    Pin Extremities \n";
  static    char    form7p[] = "    Pin number %d, z0 = %+f, z1 = %+f,"
                                                           " iz0 = %3d, iz1 = %3d\n";
  static    char    form8p[] = "    Ellipsoidal Extremity z0 = %f, z1 = %f,"
                                                             " iz0 = %d, iz1 = %d\n";
/** **************************************** ****************************************
                                   | approximate symetry axis
                                   | 
     D-4    D-3    D-2    D-1     D0    D+1    D+2    D+3    D+4
 ...      |      |      |      |      |      |      |      |        |          )
     insL   cndL   insL   cndL   insL   cndL   insL   cndL  eltrdCyl  eltrdEll
         pe0    pe1    pe2    pe3    pe4    pe5    pe6    pe7      
**************************************** **************************************** **/

  fPrintF(bufOut,"%s  BEGIN\n", prognamp);

  p->unitnamep = "millimetres";
  p->eltrdR = eltrdR_mm = ELTRD_R_mm;

  insL_mm = ELTRD_I_mm;     cndL_mm = ELTRD_C_mm;
  zndMIN_mm = CYL_Zmin_mm;  zndMAX_mm = CYL_Zmax_mm;
  symPlZ_mm = 0.0;
  eltrdExtCyl_mm = ELTRD_IX1_mm;
  eltrdExtEll_mm = ELTRD_IX0_mm;
  fPrintF(bufOut,form1p, prognamp, insL_mm, cndL_mm, 
                 eltrdExtEll_mm, eltrdExtCyl_mm, zndMIN_mm, zndMAX_mm, p->unitnamep);

  p->insL = insL_mm;  p->cndL = cndL_mm;
  p->eltrdExtEll = ELTRD_IX0_mm;
  p->eltrdExtCyl = ELTRD_IX1_mm;
  p->zndMAX = zndMAX_mm;  p->zndMIN = zndMIN_mm;
  p->symPlZ = symPlZ_mm;
                                                /** Pin Extremities z-coordinates **/
  p->pinX = pinX;
  p->pinExtX = pinExtX = 2 * pinX;
  pinExtZp = p->pinExtZp;
  pinExtZp[4] = symPlZ_mm - insL_mm/2;  pinExtZp[5] = pinExtZp[4] - cndL_mm;
  pinExtZp[6] = pinExtZp[5] - insL_mm;  pinExtZp[7] = pinExtZp[6] - cndL_mm;

  pinExtZp[3] = symPlZ_mm + insL_mm/2;  pinExtZp[2] = pinExtZp[3] + cndL_mm;
  pinExtZp[1] = pinExtZp[2] + insL_mm;  pinExtZp[0] = pinExtZp[1] + cndL_mm;

                                              /** Electrode Ellipsoidal Extremity **/
  p->eltrdExtX = 2;
  eltrdExtZp = p->eltrdExtZp;
                                      /** z value of the begining of rounded part **/
  eltrdExtZp[1] = pinExtZp[0] + eltrdExtCyl_mm;
                                                     /** z value of the extremity **/
  eltrdExtZp[0] = eltrdExtZp[1] + eltrdExtEll_mm;

                                                             /** Partition Choice **/
  deltazc = cndL_mm - insL_mm;
  deltazcic = 2.0 * insL_mm;
  deltazci  = insL_mm;
  deltazie = 2.0 * eltrdExtCyl_mm - insL_mm;
  deltazp[MYDELTAZMN_ -2] = deltazp[MYDELTAZMN_ +2] = deltazp[MYDELTAZMN_ ] = deltazcic;
  deltazp[MYDELTAZMN_ -3] = deltazp[MYDELTAZMN_ +3] = 
                           deltazp[MYDELTAZMN_ -1] = deltazp[MYDELTAZMN_ +1] = deltazc;
  deltazp[MYDELTAZMN_ -4] = deltazp[MYDELTAZMN_ +4] = deltazci;
  deltazp[MYDELTAZMN_ -5] = deltazp[MYDELTAZMN_ +5] = deltazie;
  deltazp[MYDELTAZMN_ -6] = deltazp[MYDELTAZMN_ +6] =
                                   2.0 * (eltrdExtEll_mm - eltrdExtCyl_mm) + insL_mm;

  deltazp[MYDELTAZMN_ -7] = deltazp[MYDELTAZMN_ +7] = 2.0 ;
  deltazp[MYDELTAZMN_ -8] = deltazp[MYDELTAZMN_ +8] = 4.0;
  deltazp[MYDELTAZMN_ -9] = deltazp[MYDELTAZMN_ +9] = 8.0;
  delta = 0.0;
  for(idz = 1;  idz < MYDELTAZMN_ ;  idz++){ delta += deltazp[MYDELTAZMN_ - idz];}
  deltazp[0]           = fabs(zndMIN_mm) - deltazp[MYDELTAZMN_ ]/2 - delta;
  delta = 0.0;
  for(idz = 1;  idz < MYDELTAZMN_ ;  idz++){ delta += deltazp[MYDELTAZMN_ + idz];}
  deltazp[2*MYDELTAZMN_ ] = fabs(zndMAX_mm) - deltazp[MYDELTAZMN_ ]/2 - delta;

                                                   /** Default step length choice **/
  zpap[MYDELTAZMN_ -2] = zpap[MYDELTAZMN_ +2] = zpap[MYDELTAZMN_ ] = deltazcic/20;
  zpap[MYDELTAZMN_ -3] = zpap[MYDELTAZMN_ +3] = zpap[MYDELTAZMN_ -1] = zpap[MYDELTAZMN_ +1]
                                                                         = deltazc/5;
  zpap[MYDELTAZMN_ -4] = zpap[MYDELTAZMN_ +4] = deltazci/10;
  zpap[MYDELTAZMN_ -5] = zpap[MYDELTAZMN_ +5] = deltazie/10;
  zpap[MYDELTAZMN_ -5] = zpap[MYDELTAZMN_ +6] = deltazp[MYDELTAZMN_ +6] /8;

  zpap[MYDELTAZMN_ -6] = zpap[MYDELTAZMN_ +7] = deltazp[MYDELTAZMN_ +7] /8;
  zpap[MYDELTAZMN_ -7] = zpap[MYDELTAZMN_ +8] = deltazp[MYDELTAZMN_ +8] /4;
  zpap[MYDELTAZMN_ -8] = zpap[MYDELTAZMN_ +9] = deltazp[MYDELTAZMN_ +9] /4;
  zpap[MYDELTAZMN_ -10] = zpap[MYDELTAZMN_ +10] = 3.0;

                                                 /** Effective Step length choice **/
                                       /** include==>read  "gridZ.linear4.INI..." **/
  zpap[MYDELTAZMN_] = deltazp[MYDELTAZMN_]/infozpaxp[0];
  for(idz = 1;  idz <= MYDELTAZMN_ ; idz++)
  { zpap[MYDELTAZMN_ -idz] = zpap[MYDELTAZMN_ +idz] =
                                             deltazp[MYDELTAZMN_ +idz]/infozpaxp[idz];
  }

                                                  /** Partition mesh point number **/
  for(idz = 0;  idz <= 2*MYDELTAZMN_ ; idz++)
  { ndeltazp[idz] = (deltazp[idz] + 1.E-8) / zpap[idz];
    zpap[idz] = deltazp[idz] / ndeltazp[idz];
  }
  zllx = 0;
  for(idz = 0;  idz <= 2*MYDELTAZMN_ ; idz++){ zllx += ndeltazp[idz];}
  zndx = zllx +1;

  fPrintF(bufOut,"\n");
  fPrintF(bufOut,"  Effective Values \n");
  for(idz = 0;  idz <= 2*MYDELTAZMN_ ; idz++)
  { fPrintF(bufOut,form3p, idz, ndeltazp[idz], zpap[idz], deltazp[idz]);
  }
  fPrintF(bufOut,form4p, (int)zllx, (int)zndx);
  
  p->zllx = zllx;
  p->zndx = zndx;
                                                           /** Memory Allocations **/
  gmGridZPAlloc(zllx, p);
  zndp = p->zndp;
  zllp = p->zllp;

  for(idz = 0, zllhx = 0;  idz < MYDELTAZMN_ ; idz++){ zllhx += ndeltazp[idz];}
  zllhx += ndeltazp[MYDELTAZMN_] / 2;

                                           /** Effective z nodes and link-lengths **/
  izt = 0;
  zc = zndMIN_mm;  zndp[0] = zndMIN_mm;
  for(idz = 0;  idz <= 2*MYDELTAZMN_ ; idz++)
  { for(iz = 1;  iz <= ndeltazp[idz];  iz++)
    { zllp[izt] = zpap[idz];  izt++;  zndp[izt] = zpap[idz] * iz + zc;
    }
    zc = zndp[izt];
  }
  fPrintF(bufOut,form5p, izt,zllx);
                                                      /** Pin Extremities z-index **/
  fPrintF(bufOut,form6p);
  for(is = 0;  is < p->pinX;  is++)
  { zc = p->pinExtZp[2*is];
    cmp = classDbl(&index, zc, zndp, zndx);
    if( fabs(zc - zndp[index+1])  <  fabs(zc - zndp[index])) index++;
    p->pinExtIp[2*is] = index;
    zc = p->pinExtZp[2*is +1];
    cmp = classDbl(&index, zc, zndp, zndx);
    if( fabs(zc - zndp[index+1])  <  fabs(zc - zndp[index])) index++;
    p->pinExtIp[2*is +1] = index;
    fPrintF(bufOut,form7p, is, p->pinExtZp[2*is], p->pinExtZp[2*is +1],
                                            p->pinExtIp[2*is], p->pinExtIp[2*is +1]);
  }
                                      /** Electrode Ellipsoidal Extremity z-index **/
  for(is = 0;  is < p->eltrdExtX;  is++)
  { cmp = classDbl(&index, p->eltrdExtZp[is] + 1.E-8, zndp, zndx);
    p->eltrdExtIp[is] = index;
  }
  fPrintF(bufOut,form8p, p->eltrdExtZp[0], p->eltrdExtZp[1],
                                                 p->eltrdExtIp[0], p->eltrdExtIp[1]);

                                                          /** Geometry allocation **/
  for(iz = 0, ip = p->irGeond_firstp;  iz < zndx;  iz++) *ip++ = 0;
  for(iz = 0, ip = p->irGeond_lastp;   iz < zndx;  iz++) *ip++ = 0;
  for(iz = 0, ip = p->ivGeo_firstp;   iz <= zndx;  iz++) *ip++ = 0;
  fPrintF(bufOut,"%s  END\n", prognamp);

  fPrintF(bufOut,"\n");
  return ok;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridZprint(FILE *bufOut, gmGridZ *p)
{ int       iz, izx, is, iz0, iz1;
  int      *ip;
  double   *dp;
  static    char    prognamp[] = "gridZ.linear4::gmGridZprint";
  static    char    form1p[] = "%s;  p=%p;  unit=%s \n"
                               "     radius=%f, node z_MIN=%f, node z_MAX=%f\n"
                               "     insulatorLength=%f, conductorLength=%f, %s\n";
  static    char    form3p[] = "%s;  z nodes : zndz=%d zndx=%d;  values(%s) =";
  static    char    form4p[] = "%s;  z links : zllz=%d zllx=%d; lengths(%s) =";
  static    char    form5p[] = "%s;  %d Pins (%s) =\n";
  static    char    form6p[] = "  pin number %d,"
                                             "  z0=%+f  iz0=%3d, z1=%+f  iz1=%3d \n";
  static    char    form7p[] = "%s;"
               "  %d Electrode extremities, \n"
               "  cylinderLength=%f, ellipsoidLength=%f, \n"
               "  z0=%+f  iz0=%3d  zndp[iz0]=%+f,  z1=%+f  iz1=%3d  zndp[iz1]=%+f\n";

  static    char    form8p[] = "%s;  irGeond_first, nx=%d, index values: ";
  static    char    form9p[] = "%s;  irGeond_last,  nz=%d, index values: ";
  static    char    form10p[]= "%s;  ivGeo_first, %d cumulated index values: "; 

  static    char    form20p[]= "%+9.3e ";
  static    char    form21p[] = "%3d ";
  static    char    form22p[] = "%4d ";

  fPrintF(bufOut,form1p, prognamp, p, p->unitnamep, p->eltrdR, p->zndMIN, p->zndMAX,
                                                     p->insL, p->cndL, p->unitnamep);
  fPrintF(bufOut,form3p, prognamp, p->zndz, p->zndx, p->unitnamep);
  izx = p->zndx;
  dp  = p->zndp;
  for(iz = 0;  iz < izx;  iz++, dp++)
  { if(iz % 10 == 0) fPrintF(bufOut,"\n");
    fPrintF(bufOut,form20p, *dp);
  }
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,form4p, prognamp, p->zllz, p->zllx, p->unitnamep);
  izx = p->zllx;
  dp  = p->zllp;
  for(iz = 0; iz < izx; iz++, dp++)
  { if(iz % 10 == 0) fPrintF(bufOut,"\n");
    fPrintF(bufOut,form20p, *dp);
  }
  fPrintF(bufOut,"\n");

                                                      /** Pin Extremities **/
  fPrintF(bufOut,"\n");
  fPrintF(bufOut,form5p, prognamp, p->pinX, p->unitnamep);
  dp = p->pinExtZp;  ip = p->pinExtIp;
  for(is = 0;  is < p->pinX;  is++)
  { fPrintF(bufOut,form6p, is, *dp, *ip,  *(dp+1), *(ip+1));
    dp++;  dp++;  ip++; ip++;
  }
  
                                                      /** Electrode Extremity **/
  iz0 = p->eltrdExtIp[0];  iz1 = p->eltrdExtIp[1];
  fPrintF(bufOut,form7p, prognamp, p->eltrdExtX, p->eltrdExtCyl, p->eltrdExtEll,
                                   p->eltrdExtZp[0], p->eltrdExtIp[0], p->zndp[iz0],
                                   p->eltrdExtZp[1], p->eltrdExtIp[1], p->zndp[iz1]);

                                                /** Print global geometry indices **/
  fPrintF(bufOut,"\n");
  izx = p->zndx;
  fPrintF(bufOut,form8p, prognamp, izx);
  ip = p->irGeond_firstp;
  if(izx > 0)
  { for(iz = 0;  iz < izx;  iz++)
    { if(iz % 20 == 0) fPrintF(bufOut,"\n");
      fPrintF(bufOut,form21p, *ip++);
    }
  }
  fPrintF(bufOut,"\n");
  fPrintF(bufOut,form9p, prognamp, izx);
  ip = p->irGeond_lastp;
  if(izx > 0)
  { for(iz = 0;  iz < izx;  iz++)
    { if(iz % 20 == 0) fPrintF(bufOut,"\n");
      fPrintF(bufOut,form21p, *ip++);
    }
  }
  fPrintF(bufOut,"\n");
  fPrintF(bufOut,form10p, prognamp, izx +1);
  ip = p->ivGeo_firstp;
  if(izx > 0)
  { for(iz = 0;  iz <= izx;  iz++)
    { if(iz % 20 == 0) fPrintF(bufOut,"\n");
      fPrintF(bufOut,form22p, *ip++);
    }
  }

  fPrintF(bufOut,"\n");
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridZCpy(char *unitp, double scale, gmGridZ *pf, gmGridZ *pi)
{ int       i, is;
  int      *iip, *ifp;
  double   *dip, *dfp;

  pf->unitnamep = unitp;
  pf->zndx = pi->zndx;                  pf->zllx = pi->zllx;

  dip = pi->zndp;  dfp = pf->zndp;
  for(i = 0;  i < pi->zndx;  i++){ *dfp++ = *dip++ /scale;}
  dip = pi->zllp;  dfp = pf->zllp;
  for(i = 0;  i < pi->zllx;  i++){ *dfp++ = *dip++ /scale;}

  pf->eltrdR = pi->eltrdR /scale;
  pf->zndMIN = pi->zndMIN /scale;       pf->zndMAX = pi->zndMAX /scale;
  pf->insL = pi->insL /scale;           pf->cndL = pi->cndL /scale;

  pf->pinX = pi->pinX;                  pf->pinExtX = pi->pinExtX;

  dip = pi->pinExtZp;  dfp = pf->pinExtZp;
  iip = pi->pinExtIp;  ifp = pf->pinExtIp;
  for(is = 0;  is < pi->pinExtX;  is++){ *dfp++ = *dip++ /scale;  *ifp++ = *iip++;}

  pf->symPlZ = pi->symPlZ /scale;       pf->symPlI = pi->symPlI;

  pf->eltrdExtX = pi->eltrdExtX;
  pf->eltrdExtCyl = pi->eltrdExtCyl /scale;
  pf->eltrdExtEll = pi->eltrdExtEll /scale;
  dip = pi->eltrdExtZp;  dfp = pf->eltrdExtZp;
  iip = pi->eltrdExtIp;  ifp = pf->eltrdExtIp;
  for(is = 0;  is < pi->eltrdExtX;  is++){ *dfp++ = *dip++ /scale;  *ifp++ = *iip++;}

  pf->irGeond_firstp = pi->irGeond_firstp;
  pf->irGeond_lastp  = pi->irGeond_lastp;
  pf->ivGeo_firstp   = pi->ivGeo_firstp;
  return;
}
/******************************************************************************/
/*                                                                            */
/* Read : 10 double                                                           */
/*         2 int                                                              */
/* Total  10 double into infozpap, 2 int into infopinzpaxp                    */
/******************************************************************************/
int       gmGridZreadINI(FILE *bufGridZinitp, double *infozpap, int *infozpaxp)
{ int       ok = 1;
  double   *dp;
  int      *ip, i;
  char      cp[512];
  static    char form1p[] =
             "   infozpap[0]=%f, infozpap[1]=%f, ..."
                         " infozpap[MYDELTAZMN_ -1]=%f, infozpap[MYDELTAZMN_]=%f \n";
  static    char form2p[] =
             "   infozpaxp[0]=%d, infozpaxp[1]=%d, ..."
                       " infozpaxp[MYDELTAZMN_ -1]=%d, infozpaxp[MYDELTAZMN_]=%d \n";
  static    char prognamp[] = "gridZ.linear4::gmGridZreadINI";

fPrintF(stderr,"%s  BEGIN\n", prognamp);
  dp = infozpap;  ip = infozpaxp;
                                                               /** default values **/
  *dp++ = 0.05;  *dp++ = 5;  *dp++ = 0.05;  *dp++ = 5;  *dp++ = 0.05;
  *dp++ = 0.1;  *dp++ = 0.5;  *dp++ = 1.0;  *dp++ = 2.0;  *dp = 3.0;

                                                          /** reading file values **/
  fscanf(bufGridZinitp, "%s", cp);  fPrintF(stderr, "  Reading %s\n", cp);
  fscanf(bufGridZinitp, "%s", cp);  fPrintF(stderr, "  Fisrt Line Header= %s \n", cp);
  dp = infozpap;
  for(i = 0;  i <= MYDELTAZMN_ ;  i++){ fscanf(bufGridZinitp,"%le ", dp++);  }
  fscanf(bufGridZinitp,"\n");
  fscanf(bufGridZinitp, "%s", cp);  fPrintF(stderr, "  Second Line Header= %s \n", cp);
  ip = infozpaxp;
  for(i = 0; i <= MYDELTAZMN_; i++){ fscanf(bufGridZinitp,"%d ", ip++);}
fPrintF(stderr, form1p, infozpap[0],  infozpap[1],
                                   infozpap[MYDELTAZMN_ -1],  infozpap[MYDELTAZMN_]);
fPrintF(stderr, form2p, infozpaxp[0], infozpaxp[1],
                                  infozpaxp[MYDELTAZMN_ -1], infozpaxp[MYDELTAZMN_]);
fPrintF(stderr,"%s  END\n", prognamp);
  return  ok;
}
/******************************************************************************/
/*                                                                            */
/*  given z , return the radius rho  of the electrode at this altitude        */
/*  unit is that of  gmGridZ *Zp                                              */
/******************************************************************************/
double    gmGridZGeoZ2R(double z, gmGridZ *Zp)
{ double    radius, rho;
  double    z0, z1, dzE, dz;

  z0 = Zp->eltrdExtZp[0];  z1 = Zp->eltrdExtZp[1];
  radius = Zp->eltrdR;
  if(z >= z0) return  0.0;
  if(z <= z1) return  radius;
  dzE = z1 - z0;  dz = (z1 - z)/dzE;
  rho = radius * sqrt(1.0 - dz*dz);
  return rho;
}
/******************************************************************************/
/******************************************************************************/
