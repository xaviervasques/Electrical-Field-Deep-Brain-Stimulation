/*  ../libmy/utiImage0D.c                                                     */
/*  Mennessier Gerard                   20031007                              */
/*  Last revised M.G.                   20040510                              */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utiImage0D.h"

static char    srcfilenamp[] = "utiImage0D";
/******************************************************************************/
/*  void  utiImage0DSet(utiImage0D *ep, int byteOrder,                        */
/*                                              int depth, int bytesPerPixel) */
/*                                                                            */
/******************************************************************************/
void      utiImage0DSet(utiImage0D *ep, int byteOrder, int depth, int bytesPerPixel)
{
  ep->byteOrder     = byteOrder;
  ep->depth         = depth;
  ep->bytesPerPixel = bytesPerPixel;
  return;
}
/******************************************************************************/
/*  void  utiImage0DPrint(FILE *streamp, utiImage0D *ep)                      */
/*                                                                            */
/******************************************************************************/
void      utiImage0DPrint(FILE *streamp, utiImage0D *ep)
{ int       byteOrder;
  static char    lsbp[] = "LSBfirst",  msbp[] = "MSBFirst";
  char          *lmsbpp[2] = { lsbp, msbp };
  static char    form1p[] = "%s::%s  BEGIN pointer=%p \n";
  static char    form2p[] = "%s::%s  END\n";
  static char    form3p[] = "  byteOrder=%d =%s  depth=%d  bytesPerPixel=%d\n";
  static char    prognamp[] = "utiImage0DPrint";

  byteOrder = ep->byteOrder;
  fPrintF(streamp, form1p, srcfilenamp, prognamp, (void*)ep);
  fPrintF(streamp, form3p, byteOrder, lmsbpp[byteOrder],
                                                       ep->depth, ep->bytesPerPixel);
  fPrintF(streamp, form2p, srcfilenamp, prognamp);
  return;
}
/******************************************************************************/
/*  void  utiImage0DCpy(utiImage0D *efp, utiImage0D ei)                       */
/*                                                                            */
/******************************************************************************/
void      utiImage0DCpy(utiImage0D *efp, utiImage0D ei)
{
  efp->byteOrder     = ei.byteOrder;
  efp->depth         = ei.depth;
  efp->bytesPerPixel = ei.bytesPerPixel;
  return;
}
/******************************************************************************/
/*  void  utiImage0DCpyp(utiImage0D *efp, utiImage0D *eip)                    */
/*                                                                            */
/******************************************************************************/
void      utiImage0DCpyp(utiImage0D *efp, utiImage0D *eip)
{
  efp->byteOrder     = eip->byteOrder;
  efp->depth         = eip->depth;
  efp->bytesPerPixel = eip->bytesPerPixel;
  return;
}
/******************************************************************************/
/*  short  utiImage0DCompatible(utiImage0D *e0p, utiImage0D *e1p)             */
/*                                                                            */
/*  check byteOrder and depth, BUT NOT bytesPerPixel                          */
/*  if equal at e0p and e1p return OK = 1 = true                              */
/*  else                    return 0                                          */
/******************************************************************************/
short     utiImage0DCompatible(utiImage0D *e0p, utiImage0D *e1p)
{ 
  if(e0p->byteOrder != e1p->byteOrder) return 0;
  if(e0p->depth     != e1p->depth) return 0;
  return 1;
}
/******************************************************************************/
/******************************************************************************/
