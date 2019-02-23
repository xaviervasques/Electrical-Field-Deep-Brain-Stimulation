/*  bio_phys/URMAE/numerical/linear4/pot.linear4.c                            */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20011228                                */

#include  <stddef.h>
#include  <math.h>

#include  "utistdIO.h"
#include  "utiAlloc.h"
#include  "utiClassRangeDbl.h"
#include  "eltrdPot.linear4.h"
#include  "cylPot.linear4.h"
#include  "pot.linear4.h"
#include  "varPot.linear4.h"

double  **allocVp(double *vnp, size_t rhotndx, size_t zndx);
void      printV(FILE *bufOut, double **vnpp, int rhotndx, int zndx);
void      printVMath(FILE *bufOut, double **vnpp, int rhotndx, int zndx);
int       initVfromDataFile(FILE *bufOut,double **vnpp,double *varp);

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmPotPAllocZero(gmPot *p, gmGridR *gridRp, gmGridZ *gridZp)
{ size_t    nz;
  double   *vp;
  int       i;

  p->gridRp = gridRp;  p->gridZp = gridZp;
  nz = gridRp->rhotndx * gridZp->zndx;
  vp = dblAlloc(nz, "pot.linear4::gmPotInit : vp");
  p->vp = vp;                 p->vpx = nz;
  for(i = 0;  i < nz;  i++){ *vp++ = 99.0;}
  nz = gridZp->zndx;
  p->vpp = allocVp(p->vp, gridRp->rhotndx, nz);
  p->vppx = nz;
  p->varPotp = NULL;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmPotPrint(FILE *bufOut, gmPot *p)
{ gmGridR  *gridRp;
  gmGridZ  *gridZp;
  gmVarPot *varpotp;
  static    char    form1p[] = "%s;  BEGIN \n";
  static    char    form2p[] = "  gmPotp=%p,  gridRp=%p, gridZp=%p, "
                                                         "created from varPot=%p \n";
  static    char    form3p[] = "  vp=%p, vpx=%d,    vpp=%p, vppx=%d\n";
  static    char    form9p[] = "%s;  END \n";
  static    char    prognamp[] = "pot.linear4::gmPotPrint";

  fPrintF(bufOut,form1p, prognamp);
  gridRp = p->gridRp;  gridZp = p->gridZp;
  fPrintF(bufOut,form2p, p, gridRp, gridZp, p->varPotp);
  fPrintF(bufOut,form3p, p->vp, (int)(p->vpx), p->vpp, (int)(p->vppx));

  printV(bufOut, p->vpp, gridRp->rhotndx, gridZp->zndx);

  varpotp = p->varPotp;
  if(varpotp != NULL)
  { fPrintF(bufOut,"\n");
    gmEltrdPotPrint(bufOut, varpotp->eltrdPotp);
    fPrintF(bufOut,"\n");
    gmCylPotPrint(bufOut, varpotp->cylPotp);
  }

  fPrintF(bufOut,form9p, prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmPotPrintMath(FILE *bufOut, gmPot *p)
{ gmGridR  *gridRp;
  gmGridZ  *gridZp;
  static    char    form1p[] = "%s;  BEGIN \n";
  static    char    form2p[] = "  gmPotp=%p,  gridRp=%p, gridZp=%p, "
                                                         "created from varPot=%p \n";
  static    char    form3p[] = "  vp=%p, vpx=%d,    vpp=%p, vppx=%d\n";
  static    char    form9p[] = "%s;  END \n";
  static    char    prognamp[] = "pot.linear4::gmPotPrintMath";

  fPrintF(bufOut,form1p, prognamp);
  gridRp = p->gridRp;  gridZp = p->gridZp;
  fPrintF(bufOut,form2p, p, gridRp, gridZp, p->varPotp);
  fPrintF(bufOut,form3p, p->vp, (int)(p->vpx), p->vpp, (int)(p->vppx));

  printVMath(bufOut, p->vpp, gridRp->rhotndx, gridZp->zndx);

  fPrintF(bufOut,form9p, prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double  **allocVp(double *vnp, size_t rhotndx, size_t zndx)
{ int       iz;
  double  **vnpp;

  vnpp = ptrAlloc(zndx, "potentiel::allocVp : vnpp");
  for(iz = 0; iz < zndx; iz++){ vnpp[iz] = vnp + iz * rhotndx;}
  return  vnpp;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      printV(FILE *bufOut, double **vnpp, int rhotndx, int zndx)
{ int       ir, iz;
  double   *vnp;
/*  
  static    char    form1p[] = "pot.linear4::printV;  BEGIN, Vx= rndx*zndx =%d\n";
  static    char    form9p[] = "pot.linear4::printV   END\n";
*/
  static    char    form2p[] = "  %d rhot mesh * %d z mesh = %d nodes\n";
  static    char    form3p[] = "  iz=%d, ir=0, ... ,%d\n";
  static    char    form10p[] = "%+9.6e ";

/*  fPrintF(bufOut,form1p, rhotndx * zndx); */
  fPrintF(bufOut,form2p, rhotndx,zndx, rhotndx * zndx); 
  for(iz = 0; iz < zndx; iz++)
  { vnp = vnpp[iz];
    fPrintF(bufOut,form3p, iz,rhotndx-1);
    for(ir = 0; ir < rhotndx; ir++){ fPrintF(bufOut,form10p,vnp[ir]);}
    fPrintF(bufOut,"\n");
  }
/*  fPrintF(bufOut,form9p); */
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      printVMath(FILE *bufOut, double **vnpp, int rhotndx, int zndx)
{ int       ir, iz;
  double   *vnp;
  static    char    form2p[] = "  %d rhot mesh * %d z mesh = %d nodes\n";
  static    char    form10p[] = "%+9.6e ";

  fPrintF(bufOut,form2p, rhotndx,zndx, rhotndx * zndx); 
  fPrintF(bufOut,"pot = {\n");
  for(iz = 0; iz < zndx; iz++)
  { vnp = vnpp[iz];
    for(ir = 0; ir < rhotndx; ir++){ fPrintF(bufOut,form10p,vnp[ir]);}
    fPrintF(bufOut,"\n");
  }
  fPrintF(bufOut,"}\n");
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       initVfromDataFile(FILE *bufOut,double **vnpp,double *varp)
{ int       ok = 1;
  static    char    prognamp[] = "pot.linear4::initVfromDataFile";
  static    char    forma1p[] = "%s;  TO BE DONE\n";


  fPrintF(bufOut,forma1p, prognamp);
  return ok;
}
/******************************************************************************/
/*                                                                            */
/* rho = rhot + radius                                                        */
/******************************************************************************/
double    gmPotFunc(double rho, double z, gmPot *p)
{ int       cmpz, cmpr, iz, irh;
  gmGridR  *gridrp;
  gmGridZ  *gridzp;
  double  **vpp, *vp0, *vp1, v0, v1 , v00, v01, v10, v11, v;

  vpp = p->vpp;
  gridzp = p->gridZp;
  cmpz = classDbl(&iz, z, gridzp->zndp, gridzp->zndx);
  if(iz < 0) iz++;
  else if(iz == gridzp->zndx -1) iz--;
  vpp += iz;
  vp0 = *vpp++;  vp1 = *vpp;
/* 
static    char    prognamp[] = "pot.linear4::gmPotFunc";
fPrintF(stderr, "%s  z=%f, iz=%d, vp0=%p, vp1=%p, ", prognamp, z,iz,vp0,vp1); 
*/
  gridrp = p->gridRp;
  cmpr = classDbl(&irh, rho, gridrp->rhondp, gridrp->rhotndx);
  if(irh < 0) irh++;
  else if(irh == gridrp->rhotndx -1) irh--;
/* fPrintF(stderr, " rho=%f, irh=%d \n", rho,irh); */

  vp0 += irh;  vp1 += irh;
  if(cmpr == 0){ v0 = *vp0;  v1 = *vp1;}
  else
  { v00 = *vp0++;  v01 = *vp0; 
    v0 = ( (rho - gridrp->rhondp[irh]) * v01 - (rho - gridrp->rhondp[irh+1]) * v00 )
         /gridrp->rhotllp[irh];
    v10 = *vp1++;  v11 = *vp1; 
    v1 = ( (rho - gridrp->rhondp[irh]) * v11 - (rho - gridrp->rhondp[irh+1]) * v10 )
         /gridrp->rhotllp[irh];
  } 
  if(cmpz == 0){ v = v0;}
  else
  { v = ( (z - gridzp->zndp[iz]) * v1 - (z - gridzp->zndp[iz+1]) * v0) 
        /gridzp->zllp[iz];
  }
  return v;
}
/******************************************************************************/
/*                                                                            */
/* rho = rhot + radius                                                        */
/******************************************************************************/
double    gmPotFunc2(double rho, double z, gmPot *p)
{ 
  return gmPotInterpol(rho,z, p->vpp, p->gridRp, p->gridZp);
}
/******************************************************************************/
/*                                                                            */
/*  given the grids gridzp and gridrp,                                        */
/*        the values of a function at the corresponding vertices,             */
/*                   pointed by vpp  (vpp[iz][ir])                            */
/*  interpolated linearly in both variables to get value in (z,rho)           */
/******************************************************************************/
double    gmPotInterpol(double rho, double z, double **vpp,
                                                    gmGridR *gridrp, gmGridZ *gridzp)
{ int       cmpz, cmpr, iz, irh;
  double   *vp0, *vp1, v0, v1 , v00, v01, v10, v11, v = 0.0;

  cmpz = classDbl(&iz, z, gridzp->zndp, gridzp->zndx);
  if(iz < 0) iz++;
  else if(iz == gridzp->zndx -1) iz--;
  vpp += iz;
  vp0 = *vpp++;  vp1 = *vpp;

  cmpr = classDbl(&irh, rho, gridrp->rhondp, gridrp->rhotndx);
  if(irh < 0) irh++;
  else if(irh == gridrp->rhotndx -1) irh--;

  vp0 += irh;  vp1 += irh;
  if(cmpr == 0){ v0 = *vp0;  v1 = *vp1;}
  else
  { v00 = *vp0++;  v01 = *vp0; 
    v0 = ( (rho - gridrp->rhondp[irh]) * v01 - (rho - gridrp->rhondp[irh+1]) * v00)
                                                               /gridrp->rhotllp[irh];
    v10 = *vp1++;  v11 = *vp1; 
    v1 = ( (rho - gridrp->rhondp[irh]) * v11 - (rho - gridrp->rhondp[irh+1]) * v10)
                                                               /gridrp->rhotllp[irh];
  } 
  if(cmpz == 0){ v = v0;}
  else
  { v = ( (z - gridzp->zndp[iz]) * v1 - (z - gridzp->zndp[iz+1]) * v0) 
                                                                   /gridzp->zllp[iz];
  }
  return v;
}
/******************************************************************************/
/******************************************************************************/
