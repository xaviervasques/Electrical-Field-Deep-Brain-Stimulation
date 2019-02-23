/*  URMAE/numerical/linear4/cylPot.linear4.c                                  */
/*  Mennessier Gerard                 20010522                                */
/*  Last Revised : G.M.               20010629                                */

#include  <stddef.h>
#include  "utistdIO.h"

#include  "cylPot.linear4.h"

static    char    form1p[] = "%s  BEGIN\n";
static    char    form2p[] = "%s  END\n";
static    char    form3p[] = "%s; %s = %p\n";

static    char    form21p[] = "  cylcond %d  neckMod %d  neckV %+9.6e  neckI %+9.6e \n";
static    char    form31p[] = "  cylcond %d  neckMod %d  neckV %le  neckI %le \n";

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmCylPotZero(gmCylPot *p)
{ p->infiniV = 0.0;  p->neckI = 0.0;
  p->neckMod = 0;  p->cylcond = 0;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmCylPotPrint(FILE *bufOut, gmCylPot *p)
{ static    char    prognamp[] = "cylPot.linear4::gmCylPotPrint";

  fPrintF(bufOut,form3p, prognamp, "gmCylPot pointer", p);
  fPrintF(bufOut,form21p, p->cylcond, p->neckMod, p->infiniV, p->neckI);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmCylPotWrite(FILE *bufOut, gmCylPot *p)
{
/*  static    char    prognamp[] = "cylPot.linear4::gmCylPotWrite"; */

  fPrintF(bufOut,"cylPot\n");
  fPrintF(bufOut,form21p, p->cylcond, p->neckMod, p->infiniV, p->neckI);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmCylPotRead(FILE *bufReadp, gmCylPot *p)
{ static    char    prognamp[]= "cylPot.linear4::gmCylPotRead";
  char      cp[512];

fPrintF(stderr,form1p, prognamp);

  fscanf(bufReadp,"%s", cp);
fPrintF(stderr," cp=%s \n", cp);
  fscanf(bufReadp,form31p, &(p->cylcond),&(p->neckMod), &(p->infiniV), &(p->neckI) );

fPrintF(stderr, form21p, p->cylcond, p->neckMod, p->infiniV, p->neckI);
fPrintF(stderr,form2p, prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmCylPotSetMode(gmCylPot *p, int cylcond, int neckmode)
{  p->cylcond = cylcond;   p->neckMod = neckmode;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmCylPotSetV(gmCylPot *p, double neckv)
{  p->infiniV = neckv;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmCylPotSetI(gmCylPot *p, double necki)
{ p->neckI = necki;
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      gmCylPotCpy(gmCylPot *pf, gmCylPot *pi)
{ 
  pf->infiniV = pi->infiniV;   pf->neckI = pi->neckI;
  pf->neckMod = pi->neckMod;   pf->cylcond = pi->cylcond;
  return;
}
/******************************************************************************/
/******************************************************************************/
