/*  URMAE/orientHaut/linear4.GL.V6/gm.utiCurveLevel.dump.c                    */
/*  Mennessier Gerard                 20060213                                */
/*  Last Revised : G.M.               20060213                                */

#include  <stddef.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiCurve.def.h"
#include  "utiCurve.level.h"
#include  "utiCurve.set.h"

#include  "gm.utiCurveLevel.dump.h"

#include  "gm.drawstate.E.h"
#include  "gm.drawstate.glob.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void     gmLevelDump(char *filenamp)
{ int       i, ix;
  FILE     *levelFilep;
  cSetVec  *ERZlevelSetvp;
  cSet     *csetp;
  lcSegVec *lcsegvp;
  static char    form1p[] = "  Cannot open %s = %s file\n";
  static char    prognamp[] = "gm.utiCurveLevel.dump.c::gmLevelDump";

  levelFilep = fopen(filenamp, "w+");
  if(levelFilep == NULL) 
  { myErr1(-1,stderr, prognamp, form1p, "Level Curves Dump File", levelFilep);}
  else
  { fPrintF(stderr, "%s  %s opened \n", prognamp, filenamp);}

  ERZlevelSetvp = gmDrawStateEGetRZlevelCSetp();
  ix = ERZlevelSetvp->x;
  csetp = ERZlevelSetvp->p;
fPrintF(stderr, "%s set number in ERZlevelSetv: ix=%d, first cset pointer =%p \n",
                                                         prognamp, ix, (void*)csetp);
  for(i = 0;  i < ix;  i++, csetp++)
  { if(csetp->type != MY_LCSEGV) continue;
    lcsegvp = csetp->p;
    utiCurveLevelDump(levelFilep,&stateLCrange,lcsegvp);
  }
  fflush(levelFilep);
  fclose(levelFilep);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      utiCurveLevelDump(FILE *streamp, lcZrange *lczrp, lcSegVec *lcsegvecp)
{ int       i, ix;
  lcSeg    *lcp;
  double   zo, zpa, z;
  static char    form1p[] = "%6.3f,%8.3f,%8.3f,\n"; 
  static char    prognamp[] = "gm.utiCurveLevel.dump.c::utiCurveLevelDump";

  zo = lczrp->zo;
  zpa = lczrp->zpa;
fPrintF(stderr, "%s zo = %8.3f, zpa = %8.3f \n", prognamp, zo, zpa);
  ix = lcsegvecp->x;
  lcp = lcsegvecp->p;
fPrintF(stderr, "%s segment number in lcSegVec: ix = %d, first lc pointer =%p \n",
                                                           prognamp, ix, (void*)lcp);
  for(i = 0;  i < ix;  i++, lcp++)
  { z = zo + zpa * (lcp->iz);
    fPrintF(streamp, form1p, z, lcp->y[0], lcp->x[0]);
    fPrintF(streamp, form1p, z, lcp->y[1], lcp->x[1]);
  }
  return;
}
/******************************************************************************/
/******************************************************************************/

