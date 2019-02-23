/*  URMAE/orientHaut/linear4.GL.V2/gm.drawstate.c                             */
/*  Mennessier Gerard                 20010613                                */
/*  Last Revised : G.M.               20030512                                */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utistdIO.h"
#include  "utistdErr.h"

#include  "utiVecChr.h"
#include  "utiCurve.paperFormat.h"
#include  "utiCurve.level.h"

#include  "gm.drawstate.h"
#include  "gm.drawstate.globalloc.h"
#include  "gm.linear4.initsolve.glob.h"                        /** for Vf, gradVf **/


#include  "eltrdPot.linear4.h"           /** for PinMod string to int and Neckmod **/
#include  "gm.linear4.initsolve.h"
#include  "pallidus.draw.h"
#include  "pallidus.geom.h"

#include  "constant.phys.URMAE.def.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateInit()
{ char     *cp;
  static    char      prognamp[] = "gm.drawstate::gmDrawStateInit";

  fPrintF(stderr,"%s  BEGIN\n", prognamp);
  chrPVecAlloc(&stateTitlePatientV, 16);
  chrVecIncStr(&stateTitlePatientV, "patient_code");
/**  stateTitlePatientV.x--;                         suppress final string char 0 **/

  chrPVecAlloc(&stateTitlePinModV, 5);
  chrVecIncStr(&stateTitlePinModV, "p000");
/**  stateTitlePinModV.x--; **/
  stateNeckMod = gmBiFromStrMod(statePinModp, stateTitlePinModV.p);

  chrPVecAlloc(&stateTitleVoltageV, 12);
  chrVecIncStr(&stateTitleVoltageV, "1.0");
  stateVoltage = strtod(stateTitleVoltageV.p, &cp);
/**  stateTitleVoltageV.x--; **/

  chrPVecAlloc(&stateTitleImpedanceV, 12);
  chrVecIncStr(&stateTitleImpedanceV, "1000.0");
  stateImpedance = strtod(stateTitleImpedanceV.p, &cp);
/**  stateTitleImpedanceV.x--; **/

  stateKappaEff = KAPPA_ ;

  gmDrawStateSetVEI(0);

  stateTitlesp[0] = stateTitlePatientV.p;
  stateTitlesp[1] = stateTitlePinModV.p;
  stateTitlesp[2] = stateTitleVoltageV.p;
  stateTitlesp[3] = stateTitleImpedanceV.p;

  fPrintF(stderr,"%s  END\n", prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateGetPatient()
{
  return &stateTitlePatientV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSetPatient(chrVec *patientNameVp)
{
  stateTitlePatientV.x = 0;
  chrVecIncN(&stateTitlePatientV, patientNameVp->p, patientNameVp->x);
  stateTitlesp[0] = stateTitlePatientV.p;           /** may be changed by realloc **/
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateGetPinMode()
{
  return &stateTitlePinModV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSetPinMode(chrVec *pinModVp)
{
  stateTitlePinModV.x = 0;
  chrVecIncN(&stateTitlePinModV, pinModVp->p, pinModVp->x);
  stateTitlesp[1] = stateTitlePinModV.p;           /** may be changed by realloc **/
  stateNeckMod = gmBiFromStrMod(statePinModp, stateTitlePinModV.p);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateGetVoltage()
{
  return &stateTitleVoltageV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSetVoltage(chrVec *voltageVp)
{ char     *cp;

  stateTitleVoltageV.x = 0;
  chrVecIncN(&stateTitleVoltageV, voltageVp->p, voltageVp->x);
  stateTitlesp[2] = stateTitleVoltageV.p;
  stateVoltage = strtod(stateTitleVoltageV.p, &cp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateGetImpedance()
{
  return &stateTitleImpedanceV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateSetImpedance(chrVec *impedanceVp)
{ char     *cp;

  stateTitleImpedanceV.x = 0;
  chrVecIncN(&stateTitleImpedanceV, impedanceVp->p, impedanceVp->x);
  stateTitlesp[3] = stateTitleImpedanceV.p;
  stateImpedance = strtod(stateTitleImpedanceV.p, &cp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       gmDrawStateGetVEI()
{
  return  stateVEI;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawStateSetVEI(int vei)                                                */
/*  set stateVEI, stateFpp, stateLCrange                                      */
/******************************************************************************/
void      gmDrawStateSetVEI(int vei)
{ 
                                  
  stateVEI = vei;
  if(stateVEI == 0) stateFpp = Vf.vpp;
  else              stateFpp = gradVf.vpp;

  gmDrawStateSetLevelMult();
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawStateSetLevelMult()                                                 */
/*  set stateLCrange for multiple levels                                      */
/******************************************************************************/
void      gmDrawStateSetLevelMult()
{ double    dlevelorigin = 0.0;
  int       ilevelmin = -5, ilevelmax = 5;
  double    dlevelpas = 1.0;

  if(stateVEI == 0)
  { if(stateNeckMod){ dlevelpas = 0.050;  ilevelmin = -60;  ilevelmax = +60;}
    else            { dlevelpas = 0.050;  ilevelmin = -80;  ilevelmax = +80;}
  }
  else 
  {                                                                      /** V/mm **/
    if(stateVEI == 1){ dlevelpas = 0.025;  ilevelmin = 1;  ilevelmax = +40;}
                                                                      /** mA/mm^2 **/
    else             { dlevelpas = 0.005/stateKappaEff;  ilevelmin = 1;  ilevelmax = +80;}
  }
  stateLCrange.zo  = dlevelorigin;
  stateLCrange.zpa = dlevelpas;
  stateLCrange.zn  = dlevelorigin + ilevelmin * dlevelpas;
  stateLCrange.zx  = dlevelorigin + ilevelmax * dlevelpas;
  stateLCrange.in  = ilevelmin;
  stateLCrange.ix  = ilevelmax;
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawStateSetLevelSingle()                                               */
/*  set stateLCrange for multiple levels                                      */
/******************************************************************************/
void      gmDrawStateSetLevelSingle(int vei, double levelValue)
{ int       ilevelmin = 0,  ilevelmax = 0;
  double    dlevelpas = 0.1, dlevelorigin;

  gmDrawStateSetVEI(vei);
  dlevelorigin = levelValue;

  stateLCrange.zo  = dlevelorigin;
  stateLCrange.zpa = dlevelpas;
  stateLCrange.zn  = dlevelorigin + ilevelmin * dlevelpas;
  stateLCrange.zx  = dlevelorigin + ilevelmax * dlevelpas;
  stateLCrange.in  = ilevelmin;
  stateLCrange.ix  = ilevelmax;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
lcZrange *gmDrawStateGetLCrange()
{
  return  &stateLCrange;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStatePrintLCrange()
{ lcZrange *p;

  p = gmDrawStateGetLCrange();
  fPrintF(stderr, "LEVELS INFO o=%f, pas=%f, dmin=%f, dmax=%f, imin=%d, imax=%d\n",
                                          p->zo, p->zpa, p->zn, p->zx, p->in, p->ix);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
double    gmDrawStateSolveFromBasis()
{ double    kappa = KAPPA_ ;

  gmSolveFromBasis(stateTitleVoltageV.p, stateTitlePinModV.p);
  statePredictedImpedance =  gmGetImpedance();
  stateKappaEff = kappa * (statePredictedImpedance/stateImpedance);

  return  statePredictedImpedance;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
char    **gmDrawStateGetTitles()
{
  return  stateTitlesp;
}
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/*  gmDrawStateHWScaleRatioA4Width2Height(double ratiohw, double deltaw)      */
/*  given a scale ratio  height/width, and the actual width, get the height   */
/*    assuming A4 paper format                                                */
/******************************************************************************/
double    gmDrawStateHWScaleRatioA4Width2Height(double scratiohw, double deltaw)
{ double    deltah;

  if(scratiohw < .01) scratiohw = 0.01;
  deltah = A4H_ *(deltaw/A4W_) /scratiohw;
  return  deltah;
}
/******************************************************************************/
/*                                                                            */
/*  gmDrawStateHWScaleRatioA4Height2Width(double ratiohw, double deltaw)      */
/*  given a scale ratio  height/width, and the actual width, get the height   */
/*    assuming A4 paper format                                                */
/******************************************************************************/
double    gmDrawStateHWScaleRatioA4Height2Width(double scratiohw, double deltah)
{ double    deltaw;

  if(scratiohw > 100.) scratiohw = 100.;
  deltaw = A4W_ *(deltah/A4H_) *scratiohw;
  return  deltaw;
}
/******************************************************************************/
/******************************************************************************/
