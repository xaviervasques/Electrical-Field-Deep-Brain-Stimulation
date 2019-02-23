/*  URMAE/orientHaut/linear4.GL.V2/gm.drawstate.threshold.c                   */
/*  Mennessier Gerard                 20030509                                */
/*  Last Revised : G.M.               20030513                                */

#include  <stddef.h>
#include  <math.h>                                                   /** for fabs **/

#include  "utistdIO.h"
#include  "utistdErr.h"

#include  "utiVecChr.h"

#include  "gm.drawstate.threshold.h"
#include  "gm.drawstate.threshold.statalloc.h"

#include  "gm.threshold.h"                               /** for current E=grad V **/
#include  "gm.drawstate.h"
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateThresholdInit()
{ static    char      stateThresholdstrp[] = "1000.0";
  chrVec    cV = {0,0,NULL};
  static    char      prognamp[] = "gm.drawstate::gmDrawStateThresholdInit";

  fPrintF(stderr,"%s  BEGIN\n", prognamp);
  chrPVecAlloc(&stateThresholdV, 16);
  chrPVecAlloc(&stateVolumeGlobalV, 16);
  chrPVecAlloc(&stateVolumePallidusV, 16);

  chrVecIncStr(&cV, stateThresholdstrp);
  gmDrawStateThresholdSet(&cV);

  fPrintF(stderr,"%s  END\n", prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmDrawStateThresholdSet(chrVec *thresholdVp)
{ char     *cp;
  char      bufp[32];
  int       vei = 1;
  static    char      prognamp[] = "gm.drawstate.threshold::gmDrawStateSetThreshold";

  stateThresholdV.x = 0;
  chrVecIncN(&stateThresholdV, thresholdVp->p, thresholdVp->x);
  stateThreshold = strtod(stateThresholdV.p, &cp);

  stateVolumeGlobal = 0.0;
  stateVolumePallidus = 0.0;

/** 
New THRESHOLD : COMPUTE VOLUMES
**/
  stateVolumeGlobal = gmThresholdGetVolGlob(stateThreshold);
fprintf(stderr, "%s, NUM threshold=%12.2e, VolumeGlobal=%12.2e, VolumePallidus=%12.2e \n",
          prognamp, stateThreshold, stateVolumeGlobal, stateVolumePallidus);

  sprintf(bufp, "%7.1f", stateVolumeGlobal);
  stateVolumeGlobalV.x = 0;
  chrVecIncStr(&stateVolumeGlobalV, bufp);
  
  sprintf(bufp, "%7.1f", stateVolumePallidus);
  stateVolumePallidusV.x = 0;
  chrVecIncStr(&stateVolumePallidusV, bufp);

fprintf(stderr, "%s, STR threshold=%s, VolumeGlobal=%s, VolumePallidus=%s \n",
          prognamp, stateThresholdV.p, stateVolumeGlobalV.p, stateVolumePallidusV.p);

/** Set new level curves **/
  gmDrawStateSetLevelSingle(vei, stateThreshold);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateGetThreshold()
{
  return &stateThresholdV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateGetVolumeGlobal()
{
  return &stateVolumeGlobalV;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
chrVec   *gmDrawStateGetVolumePallidus()
{
  return &stateVolumePallidusV;
}
/******************************************************************************/
/******************************************************************************/
