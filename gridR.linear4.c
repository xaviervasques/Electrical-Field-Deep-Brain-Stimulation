/*  URMAE/orientHaut/linear4.Common.V1/gridR.linear4.c                        */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20020404                                */

#include  <stddef.h>
#include  <math.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"
#include  "utiClassRangeDbl.h"

#include  "gridR.linear4.delta.def.h"
#include  "gridR.linear4.h"
#include  "cyl.def.h"
#include  "eltrd.def.h"


/******************************************************************************/
/*                                                                            */
/******************************************************************************/
gmGridR  *gmGridRAlloc(size_t nz, char *progcallp)
{ gmGridR  *p;
  static    char    prognamp[] = "gridR.linear4::gmGridRAlloc";
  static    char    form1p[] = "called from %s. malloc failed; gmGridR size nz=%d\n";

  p = (gmGridR *)malloc( nz*sizeof(gmGridR) );
  if(p == NULL) myErr1(-1,stderr,prognamp, form1p,progcallp,nz);
  p->rhotndp = p->rhondp = p->rhotllp = NULL;
  p->rhotndz = p->rhotndx = 0;
  p->rhotllz = p->rhotllx = 0;
  p->eltrdR = 0.0;
  p->rhotndMAX = p->rhondMAX = 0.0;
  p->neckRho = 0.0;
  p->neckIndx = 0;
  p->cylcond = 0;  p->cylgeom = 1;
  return p;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridRPAlloc(size_t rllz, gmGridR *p)
{ size_t    rndz;

  rndz = rllz + 1;
  p->rhotndp = dblAlloc(rndz,  "gridR.linear4::gmGridRPAlloc  rhotndp");
  p->rhotndz = rndz;
  p->rhondp  = dblAlloc(rndz,  "gridR.linear4::gmGridRPAlloc  rhondp");
  p->rhotllp = dblAlloc(rllz,  "gridR.linear4::gmGridRPAlloc  rhotllp");
  p->rhotllz = rllz;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       gmGridRinit(FILE *bufOut, gmGridR *p, double *inforpap, int cylcond)
{ int       ok = 1;
  double    eltrdR_mm;
  double    RhotndMAX_mm, RhondMAX_mm, neckRho_mm;
  double    rc;

  double    deltarp[MYDELTARN_ ],  rpap[MYDELTARN_ ], delta;
  int       ndeltarp[MYDELTARN_ ], idr;
  int       ir, irt;
  int       ndr0;
  size_t    rllx, rndx;
  double   *rndp, *rllp;
  int       cmp, index;
  static    char    prognamp[]  = "gridR.linear4::gmGridRinit";
  static    char    form1p[] = "%s; "
                          " Relectrode =%f, rhotndMAX =%f, rhondMAX =%f (%s)\n";
  static    char    form3p[] = "    deltar index %2d, %2d intv length %f, tot=%f\n";
  static    char    form4p[] = "    TOTAL r intv = %d, TOTAL rho nodes = %d\n";
  static    char    form5p[] = "    irt = %d,  expected = %d\n";


/** **************************************** ****************************************
**************************************** **************************************** **/
  fPrintF(bufOut,"%s BEGIN\n", prognamp);

  p->unitnamep = "millimetres";
  p->eltrdR = eltrdR_mm = ELTRD_R_mm;

  RhondMAX_mm = CYL_R_mm;
  RhotndMAX_mm = RhondMAX_mm - eltrdR_mm;
  neckRho_mm = NECK_R_mm;
  fPrintF(bufOut,form1p,prognamp, eltrdR_mm, RhotndMAX_mm, RhondMAX_mm, p->unitnamep);
  p->rhotndMAX = RhotndMAX_mm;          p->rhondMAX  = RhondMAX_mm;

                                                             /** Partition Choice **/
  deltarp[0] = eltrdR_mm;
  deltarp[1] = 1.0 * eltrdR_mm;
  deltarp[2] = 1.0 * eltrdR_mm;
  deltarp[3] = 2.0 * eltrdR_mm;
  deltarp[4] = 6.0 * eltrdR_mm;
  deltarp[5] = 15.0 * eltrdR_mm;
  delta = 0.0;
  for(idr = 0;  idr <= 5;  idr++)
  { delta += deltarp[idr];  ndeltarp[idr] = 1;
  }
  deltarp[6] = neckRho_mm - delta;
  deltarp[7] = deltarp[6];
  delta = 0.0;
  for(idr = 0;  idr < MYDELTARN_ -1;  idr++)
  { delta += deltarp[idr];  ndeltarp[idr] = 1;
  }
  deltarp[MYDELTARN_ -1] = RhondMAX_mm - delta;

                                                   /** Default step length choice **/

                                                 /** Effective Step length choice **/
                                       /** include==>read  "gridR.linear4.INI..." **/
  neckRho_mm = *inforpap++;
  for(idr = 0;  idr < MYDELTARN_ ;  idr++){ rpap[idr] = *inforpap++;}

                                                  /** Partition mesh point number **/
  for(idr = 0;  idr < MYDELTARN_ ;  idr++)
  { ndeltarp[idr] = (deltarp[idr] + 1.E-8) / rpap[idr];
    if(ndeltarp[idr] == 0) ndeltarp[idr] = 1;
    rpap[idr] = deltarp[idr] / ndeltarp[idr];
  }
  rllx = 0;
  for(idr = 0;  idr < MYDELTARN_ ;  idr++){ rllx += ndeltarp[idr];}
  rndx = rllx +1;

  fPrintF(bufOut,"\n");
  fPrintF(bufOut,"  Effective Values \n");
  for(idr = 0;  idr < MYDELTARN_ ;  idr++)
  { fPrintF(bufOut,form3p, idr, ndeltarp[idr], rpap[idr], deltarp[idr]);
  }
  fPrintF(bufOut,form4p, (int)rllx, (int)rndx);
  
  p->rhotllx = rllx;
  p->rhotndx = rndx;
  ndr0 = p->eltrdRhollI = p->eltrdRhondI = ndeltarp[0];
                                                           /** Memory Allocations **/
  gmGridRPAlloc(rllx, p);
                                         /** Effective rho nodes and link-lengths **/
  rndp = p->rhondp;
  rllp = p->rhotllp;
  irt = 0;
  rc = 0.0;  rndp[0] = 0.0;
  for(idr = 0;  idr < MYDELTARN_ ; idr++)
  { for(ir = 1;  ir <= ndeltarp[idr];  ir++)
    { rllp[irt] = rpap[idr];  irt++;  rndp[irt] = rc + rpap[idr] * ir;
    }
    rc = rndp[irt];
  }
  fPrintF(bufOut,form5p, irt,rllx);
                                                         /** Effective rhot nodes **/
  ndr0 = ndeltarp[0];
fPrintF(bufOut,"ndr0 = %d \n", ndr0);
  rndp = p->rhotndp;
  irt = ndr0;
  rc = 0.0;  rndp[irt] = 0.0;
  for(idr = 1;  idr < MYDELTARN_ ; idr++)
  { for(ir = 1;  ir <= ndeltarp[idr];  ir++)
    { irt++;  rndp[irt] = rc + rpap[idr] * ir;
    }
    rc = rndp[irt];
  }
  idr = 0;   irt = ndr0 -1;
  rc = 0.0;
  for(ir = 1;  ir <= ndr0;){ rndp[irt--] = rc - rpap[0] * ir++;}
                                                                         /** Neck **/
  p->cylcond = cylcond;
  if(cylcond == 0) p->cylgeom = 1;
  else             p->cylgeom = 0;
  p->neckRho = neckRho_mm;
  cmp = classDbl(&index, neckRho_mm + 1.E-8, p->rhondp, rndx);
  p->neckIndx = index;
  fPrintF(bufOut,"%s END\n", prognamp);

  fPrintF(bufOut,"\n");
  return ok;    
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridRprint(FILE *bufOut, gmGridR *p)
{ int       ir, irx, irneck;
  double   *dp;
  static    char    prognamp[] = "gridR.linear4::gmGridRprint";
  static    char    form1p[] = "%s;  p=%p;  unit=%s \n" 
                           "  radius=%f,\n  node rhot_MAX=%f, node rho_MAX=%f, %s\n";
  static    char    form3p[] = "%s;  rhot nodes : rndz=%d rndx=%d;  values(%s) =";
  static    char    form4p[] = "%s;  rho  nodes : rndz=%d rndx=%d;  values(%s) =";
  static    char    form5p[] = "%s;  rho  links : rllz=%d rllx=%d; lengths(%s) =";
  static    char    form6p[] = "%s;  cylcond=%d  cylgeom=%d  neck_radius=%+f  index=%2d  rho[%d]=%+f\n";
  static    char    form10p[] = "%+9.3e ";
  
  fPrintF(bufOut,form1p, prognamp, p, p->unitnamep,
                                 p->eltrdR, p->rhotndMAX, p->rhondMAX, p->unitnamep);
  fPrintF(bufOut,form3p, prognamp, p->rhotndz, p->rhotndx, p->unitnamep);
  irx = p->rhotndx;
  dp  = p->rhotndp;
  for(ir = 0; ir < irx; ir++, dp++)
  { if(ir % 10 == 0) fPrintF(bufOut,"\n");
    fPrintF(bufOut,form10p, *dp);
  }
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,form4p, prognamp, p->rhotndz, p->rhotndx, p->unitnamep);
  irx = p->rhotndx;
  dp  = p->rhondp;
  for(ir = 0;  ir < irx;  ir++, dp++)
  { if(ir % 10 == 0) fPrintF(bufOut,"\n");
    fPrintF(bufOut,form10p, *dp);
  }
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,form5p, prognamp, p->rhotllz, p->rhotllx, p->unitnamep);
  irx = p->rhotllx;
  dp  = p->rhotllp;
  for(ir = 0;  ir < irx;  ir++, dp++)
  { if(ir % 10 == 0) fPrintF(bufOut,"\n");
    fPrintF(bufOut,form10p, *dp);
  }
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,"\n");
  irneck = p->neckIndx;
  fPrintF(bufOut,form6p, prognamp, p->cylcond, p->cylgeom, p->neckRho, irneck, irneck, p->rhondp[irneck]);
  fPrintF(bufOut,"\n");
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmGridRCpy(char *unitp, double scale, gmGridR *pf, gmGridR *pi)
{ int       i;
  double   *dip, *dfp;

  pf->unitnamep = unitp;
  pf->rhotndz = pi->rhotndz;  pf->rhotndx = pi->rhotndx;
  pf->rhotllz = pi->rhotllz;  pf->rhotllx = pi->rhotllx;

  dip = pi->rhotndp;  dfp = pf->rhotndp;
  for(i = 0;  i < pi->rhotndx;  i++){ *dfp++ = *dip++ /scale;}
  dip = pi->rhondp;  dfp = pf->rhondp;
  for(i = 0;  i < pi->rhotndx;  i++){ *dfp++ = *dip++ /scale;}
  dip = pi->rhotllp;  dfp = pf->rhotllp;
  for(i = 0;  i < pi->rhotllx;  i++){ *dfp++ = *dip++ /scale;}

  pf->eltrdR = pi->eltrdR /scale;
  pf->eltrdRhondI = pi->eltrdRhondI;    pf->eltrdRhollI = pi->eltrdRhollI;

  pf->rhotndMAX = pi->rhotndMAX /scale; pf->rhondMAX = pi->rhondMAX /scale;

  pf->cylcond = pi->cylcond;            pf->cylgeom = pi->cylgeom;
  pf->neckIndx = pi->neckIndx;
  pf->neckRho = pi->neckRho /scale;
  return;
}
/******************************************************************************/
/*                                                                            */
/* Read : neckRho_mm                                                          */
/*        7 double                                                            */
/* Total  8 double into inforpap                                              */
/******************************************************************************/
int       gmGridRreadINI(FILE *bufGridRinitp, double *inforpap, int *cylcondp)
{ int       ok = 1;
  double   *dp, neckRho_mm;
  int       i, cylcond = 0;
  char      cp[512];
  static char    prognamp[] = "gridR.linear4::gmGridRreadINI";
  
fPrintF(stderr,"%s  BEGIN\n", prognamp);
  dp = inforpap;
                                                               /** default values **/
  *dp =  0.0;               /** BI | MONO CYL ==> CYLI or CYLC neckRho_mm =  0.0; **/
  *dp = 20.0;                                         /** NECK neckRho_mm = 20.0; **/
  dp++;
  *dp++ = 0.1;
  *dp++ = 0.05;
  *dp++ = 0.10;  *dp++ = 0.20;
  *dp++ = 1.0;  *dp++ = 3.0;  *dp++ = 5.0;

  *cylcondp = 0;                                        /** default is NECK model **/
                                                          /** reading file values **/
  fscanf(bufGridRinitp, "%s", cp);  fPrintF(stderr, "  Reading %s\n", cp);
  fscanf(bufGridRinitp, "%s", cp);  fPrintF(stderr, "  Fisrt Line Header= %s \n", cp);
  fscanf(bufGridRinitp,"%d ", &cylcond);
  *cylcondp = cylcond;

  fscanf(bufGridRinitp, "%s", cp);  fPrintF(stderr, "  Second Line Header= %s \n", cp);
  fscanf(bufGridRinitp,"%le ", &neckRho_mm);
  dp = inforpap;
  *dp = neckRho_mm;  dp++;
fPrintF(stderr, "..  cylConductivityModel=%d, neckRho_mm= %f\n", cylcond, neckRho_mm);
  fscanf(bufGridRinitp,"\n");

  fscanf(bufGridRinitp, "%s", cp);  fPrintF(stderr, "  Third Line Header= %s \n", cp);
  for(i = 1;  i <= MYDELTARN_ ;  i++){ fscanf(bufGridRinitp,"%le ", dp++);}
  fscanf(bufGridRinitp,"\n");
fPrintF(stderr, "..  cylConductivityModel=%d \n", cylcond);
fPrintF(stderr, "..  inforpap[0]=%e, inforpap[1]=%e, inforpap[2]=%e, "
                "..., inforpap[MYDELTARN_ -1]=%e, inforpap[MYDELTARN_]=%e\n",
                                      inforpap[0], inforpap[1], inforpap[2],
                                      inforpap[MYDELTARN_ -1], inforpap[MYDELTARN_]);
fPrintF(stderr,"%s  END\n", prognamp);
  return  ok;
}
/******************************************************************************/
/******************************************************************************/
