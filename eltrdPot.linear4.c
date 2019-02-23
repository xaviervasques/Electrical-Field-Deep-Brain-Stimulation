/*  bio_phys/URMAE/numerical/linear4/eltrdPot.linear4.c                       */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20030515                                */

#include  <stddef.h>
#include  <math.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"
#include  "eltrdPot.linear4.h"
#include  "pot_funcBasis.def.h"

static    char    form1p[] = "%s  BEGIN\n";
static    char    form2p[] = "%s  END\n";
static    char    form3p[] = "%s; %s = %p\n";
static    char    form4p[] = "%f ";
static    char    form5p[] = "%+9.6e ";

static    char    form21p[] = "  pin_extremity_number  %d ,total_coef_number %d\n";
static    char    form22p[] = "  pin_extremity_index   %d ,devel_order %d\n";
static    char    form23p[] = "  pin_number  %d ,variable_pin_number  %d\n";
static    char    form24p[] = "  pin_index   %d ,mode %d ,V %+9.6e ,I %+9.6e \n";
static    char    form25p[] = "  power %+9.6e  effI %+9.6e \n";

static    char    form30p[] = "%le ";
static    char    form31p[] = "  pin_extremity_number  %d ,total_coef_number %d\n";
static    char    form32p[] = "  pin_extremity_index   %d ,devel_order %d\n";
static    char    form33p[] = "  pin_number  %d ,variable_pin_number  %d\n";
static    char    form34p[] = "  pin_index   %d ,mode %d ,V %le ,I %le \n";
static    char    form35p[] = "  power %le  effI %le \n";

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmEltrdPotZero(gmEltrdPot *p)
{ int       is, i;
  int      *sorderp;
  double  **scoefpp, *scoefp;

  p->scoefpp[0] = p->scoefmatp;
  for(is = 1;  is < NPEX_;  is++){ p->scoefpp[is] = p->scoefpp[is -1] + NCOEF_;}

  p->scoefx = 0;
  sorderp = p->sorderp;
  scoefpp = p->scoefpp;
  for(is = 0;  is < NPEX_;  is++)
  { *sorderp++ = 0;
    scoefp = *scoefpp++;
    for(i = 0;  i < NCOEF_; i++){ *scoefp++ = 0.0 ;}
  }

  p->pinVarX = 0;
  for(is = 0;  is < NPX_;  is++)
  { p->pinModp[is] = -1;  p->pinVp[is] = 0.0;  p->pinIp[is] = 0.0;
  }
  p->power = 0.0;  p->effI = 0.0;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmEltrdPotPrint(FILE *bufOut, gmEltrdPot *p)
{ int       is, isx;
  int       io, iox;
  static    char    prognamp[] = "eltrdPot.linear4::eltrdPotPrint";

  fPrintF(bufOut,form3p, prognamp, "gmEltrdPot", p);
 
  isx = NPEX_;
  fPrintF(bufOut,form21p, isx, p->scoefx);
  for(is = 0;  is < isx;  is++)
  { iox = p->sorderp[is];
    fPrintF(bufOut,form22p, is, iox);
    for(io = 0;  io <= iox;  io++){ fPrintF(bufOut,form4p, p->scoefpp[is][io]);}
    fPrintF(bufOut,"\n");
  }
  fPrintF(bufOut,"\n");

  isx = NPX_;
  fPrintF(bufOut,form23p, isx, p->pinVarX);
  for(is = 0;  is < isx;  is++)
  { fPrintF(bufOut,form24p, is, p->pinModp[is], p->pinVp[is], p->pinIp[is]);
  }
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,form25p, p->power, p->effI);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmEltrdPotWrite(FILE *bufOut, gmEltrdPot *p)
{ int       is, isx;
  int       io, iox;
/*  static    char    prognamp[] = "eltrdPot.linear4::eltrdPotWrite"; */

  fPrintF(bufOut,"cyl_eltrdPot.linear4.INI.____.________.____.___\n");
  fPrintF(bufOut,"eltrdPot\n");
  isx = NPEX_;
  fPrintF(bufOut,form21p, isx, p->scoefx);
  for(is = 0;  is < isx;  is++)
  { iox = p->sorderp[is];
    fPrintF(bufOut,form22p, is, iox);
    for(io = 0;  io <= iox;  io++){ fPrintF(bufOut,form5p, p->scoefpp[is][io]);}
    fPrintF(bufOut,"\n");
  }

  isx = NPX_;
  fPrintF(bufOut,form23p, isx, p->pinVarX);
  for(is = 0;  is < isx;  is++)
  { fPrintF(bufOut,form24p, is, p->pinModp[is], p->pinVp[is], p->pinIp[is] );
  }
  fPrintF(bufOut,"\n");

  fPrintF(bufOut,form25p, p->power, p->effI);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmEltrdPotRead(FILE *bufReadp, gmEltrdPot *p)
{ int       is, iss, isx, io, iox, scoefx;
  int       mode, pinVarX;
  double    scoef, v, cur, power, effI;
  char      cp[512];
  static    char    prognamp[] = "eltrdPot.linear4::eltrdPotRead";

fPrintF(stderr,form1p, prognamp);

  fscanf(bufReadp,"%s", cp);
fPrintF(stderr," cp=%s \n", cp);
  fscanf(bufReadp,"%s", cp);
fPrintF(stderr," cp=%s \n", cp);

  fscanf(bufReadp,form31p, &isx, &scoefx);
fPrintF(stderr,form21p, isx, scoefx);
  if(isx != NPEX_)
  { fPrintF(stderr,"ERROR isx=%d, expected=%d\n", isx, NPEX_);
  }
  p->scoefx = scoefx;
  for(is = 0;  is <  NPEX_;  is++)
  { fscanf(bufReadp,form32p, &iss, &iox);
fPrintF(stderr,form22p, iss, iox);
    p->sorderp[is] = iox;
    for(io = 0;  io <= iox;  io++)
    { fscanf(bufReadp,form30p, &scoef);  p->scoefpp[is][io] = scoef;
    }
    fscanf(bufReadp,"\n");
for(io = 0;  io <= iox;  io++){ fPrintF(stderr,form4p, p->scoefpp[is][io]);}
fPrintF(stderr,"\n");
  }
  fscanf(bufReadp,"\n");

  fscanf(bufReadp,form33p, &isx, &pinVarX);
fscanf(bufReadp,"\n");
fPrintF(stderr,form23p, isx, pinVarX);
  p->pinVarX = pinVarX;
  for(is = 0;  is < NPX_;  is++)
  { fscanf(bufReadp,form34p, &iss, &mode, &v, &cur); fscanf(bufReadp,"\n");
    p->pinModp[is] = mode;  p->pinVp[is] = v;  p->pinIp[is] = cur;
fPrintF(stderr,form24p, iss, mode, v, cur);
  }
  fscanf(bufReadp,"\n");

  fscanf(bufReadp,form35p, &power, &effI);
  p->power = power;  p->effI = effI;

fPrintF(stderr,form25p, power, effI);
fPrintF(stderr,form2p, prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmEltrdPotSetSCoef(gmEltrdPot *p, double **scoefpp, int *sorderp)
{ int       is, io;
  int       scoefx = 0, order;
  double   *scoefp, *elscoefp;

  for(is = 0;  is < NPEX_;  is++)
  { order = sorderp[is];
    p->sorderp[is] = order;
    if(order < 0) continue;
    scoefp = scoefpp[is];  elscoefp = p->scoefpp[is];
    for(io = 0;  io <= order;  io++){ *elscoefp++ = *scoefp++;}
    scoefx += order + 1;
  }
  p->scoefx = scoefx;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmEltrdPotSetModes(gmEltrdPot *p, int *modp)
{ int       is, pinvx = 0;

  for(is = 0;  is < NPX_ ;  is++, modp++)
  { p->pinModp[is] = *modp;
    if(*modp > 0) pinvx++;
  }
  p->pinVarX = pinvx;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmEltrdPotSetV(gmEltrdPot *p, double *vp)
{ int       is;

  for(is = 0;  is < NPX_;  is++){ p->pinVp[is] = *vp++;}
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmEltrdPotSetI(gmEltrdPot *p, double *ip)
{ int       is;

  for(is = 0;  is < NPX_;  is++){ p->pinIp[is] = *ip++;}
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmEltrdPotCpy(gmEltrdPot *pf, gmEltrdPot *pi)
{ int       is, isx, io, order;
  int       scoefx;
  static    char   form11p[] = "find scoefx=%d, initial scoefx=%d DIFFERS STOP\n";
  static    char   prognamp[] = "eltrdPot.linear4::gmEltrdPotCpy";
 
  isx = NPEX_;
  scoefx = 0;
  for(is = 0;  is < isx;  is++)
  { order = pi->sorderp[is];
    pf->sorderp[is] = order;
    scoefx += order +1;
    for(io = 0;  io <= order;  io++){ pf->scoefpp[is][io] = pi->scoefpp[is][io];}
  }
  pf->scoefx = pi->scoefx;
  if(scoefx != pi->scoefx){ fPrintF(stderr,prognamp,form11p, scoefx, pi->scoefx);}

  isx = NPX_;
  for(is = 0;  is < isx;  is++)
  { pf->pinModp[is] = pi->pinModp[is];
    pf->pinVp[is] = pi->pinVp[is];
    pf->pinIp[is] = pi->pinIp[is];
  }
  pf->pinVarX = pi->pinVarX;
  pf->power = pi->power;  pf->effI = pi->effI;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       gmBiFromPinV(double *pinVp)
{ int       bi = 0, existp = 0, existm = 0;
  int       i;
  double    v;

  for(i = 0;  i < NPX_; i++)
  { v = *pinVp++;
    if(v > 0) existp = 1;
    else if(v < 0) existm = 1;
  }
  if(existp && existm) bi = 1;
  return bi;
}
/******************************************************************************/
/*                                                                            */
/*  set the int pinmodi[] values,                                             */
/*  and return the neckmode "bi"                                              */
/******************************************************************************/
int       gmBiFromStrMod(int *pinmodip, char *pinmodstrp)
{ int       bi = 0, existp = 0, existm = 0, modi, i;
  char      c;

  for(i = 0;  i < NPX_; i++)
  { c = *pinmodstrp++;
    modi = 1;
         if(c == 'p')            { existp = 1;  modi = 0;}
    else if(c == 'm' || c == 'n'){ existm = 1;  modi = 0;}
    *pinmodip++ = modi;
  }
  if(existp && existm) bi = 1;
  return bi;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double    getDiffPot(double  *pinVp, double neckV)
{ double    vn, vx, v;
  int       is;
  static    char   prognamp[] = "eltrdPot.linear4::getDiffPot";

  fPrintF(stderr,"%s  BEGIN\n", prognamp);
  vn = vx = neckV;
  for(is = 0;  is < NPX_;  is++)
  { v = *pinVp++;
         if(v > vx) vx = v;
    else if(v < vn) vn = v;
  }
  fPrintF(stderr,"%s  END\n", prognamp);
  return (vx - vn);
}
/******************************************************************************/
/******************************************************************************/
