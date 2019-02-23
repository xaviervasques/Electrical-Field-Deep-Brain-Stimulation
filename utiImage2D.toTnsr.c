/*  ../libmy/utiImage2D.toTnsr.c                                              */
/*  Mennessier Gerard                   20041011                              */
/*  Last revised M.G.                   20041012                              */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiVecChr.h"

#include  "utiByteOrderENDIAN.def.h"
#include  "utiByteOrderENDIAN.h"
#include  "utiImage0D.h"
#include  "utiImage2D.h"
#include  "utiTnsr.h"

#include  "utiImage2D.toTnsr.h"

static char    srcfilenamp[] = "utiImage2D.toTnsr.c";
/******************************************************************************/
/* It is implicitly assumed that pixels are coding GRAY levels, not colors    */
/* so that it make sense to convert pixel to integer or float or double       */
/*                                                                            */
/* It is also ASSUMED, but may cause ALIGNMENT problems, that                 */
/*               char pointer ep->dataV.p                                     */
/* can be cast onto a                                                         */
/*              short pointer  or  int pointer                                */
/*                                                                            */
/* in instructions                                                            */
/*         (unsigned short *)ucp and (unsigned int *)ucp                      */
/******************************************************************************/

/******************************************************************************/
/*  void  utiImage2DtoTnsrInt(tnsr2int *tp, utiImage2D *ep)                   */
/*                                                                            */
/******************************************************************************/
void      utiImage2DtoTnsrInt(tnsr2int *tp, utiImage2D *ep)
{ short    *firstp;
  int       machineByteOrder;

  int       eByteOrder, needSwap = 1;
  int       bytesPerPixel, oknbp = 0;  
  size_t    w, h, wbz, wbx;  
  unsigned char *ucip, *ucp;
  unsigned short ush, *ushp;
  unsigned int   ui, *uip;
  int     **ipp, *ip;
  int       i, j;
  static char    form1p[] = "UNEXPECTED case STOP bytesPerPixel=%d; expected 1,2,4 \n";
  static char    prognamp[] = "utiImage2DtoTnsrInt";

  firstp = utiByteOrderGet();
  if(*firstp == 0) machineByteOrder = myMSBFirst;
  else             machineByteOrder = myLSBFirst;

  eByteOrder = (ep->im0D).byteOrder;
  if(machineByteOrder == eByteOrder) needSwap = 0;

  bytesPerPixel = (ep->im0D).bytesPerPixel;
  if(bytesPerPixel == 1) needSwap = 0;

  if(bytesPerPixel == 1 || bytesPerPixel == 2 ||bytesPerPixel == 4) oknbp = 1;
  if(!oknbp)
  { myErr2(myILLEGAL_CASE, stderr, srcfilenamp, prognamp, form1p, bytesPerPixel);
  }

  w = ep->wpx;  h = ep->hpx;  wbz = ep->wbz;  wbx = w * bytesPerPixel;
  ucip = (unsigned char *)ep->dataV.p;

  if(tp->p == NULL || tp->ev1z < h  || tp->ev2z < w)
  { tnsr2intPAlloc(tp, h, w);
  }
  tp->ev1x = 0;  tp->ev2x = 0;

  ipp = tp->pp;
  for(i = 0;  i < h;  i++, ucip += wbz, ipp++)
  { ucp = ucip;
    ip = *ipp;
         if(bytesPerPixel == 1)
    { for(j = 0;  j < w;  j++){ ush = (unsigned short)*ucp++;  *ip++ = (int)ush;}
    }
    else if(bytesPerPixel == 2)
    { if(needSwap == 0)
      { for(j = 0;  j < w;  j++, ucp += 2)
        { ushp = (unsigned short *)ucp;  ush = *ushp;  *ip++ = (int)ush;}
      }
      else
      { for(j = 0;  j < w;  j++, ucp += 2)
        { ush = utiSwapB2ushort(ucp);  *ip++ = (int)ush;}
      }
    }
    else if(bytesPerPixel == 4)
    { if(needSwap == 0)
      { for(j = 0;  j < w;  j++, ucp += 4)
        { uip = (unsigned int *)ucp; ui = *uip;  *ip++ = (int)ui;}
      }
      else
      { for(j = 0;  j < w;  j++, ucp += 4)
        { ui = utiSwapB4uint(ucp);  *ip++ = (int)ui;}
      }
    }
  }

  tp->ev1x = h;  tp->ev2x = w;
  return;
}
/******************************************************************************/
/*  void  utiImage2DtoTnsrDbl(tnsr2dbl *tp, utiImage2D *ep)                   */
/*                                                                            */
/******************************************************************************/
void      utiImage2DtoTnsrDbl(tnsr2dbl *tp, utiImage2D *ep)
{ short    *firstp;
  int       machineByteOrder;

  int       eByteOrder, needSwap = 1;
  int       bytesPerPixel, oknbp = 0;  
  size_t    w, h, wbz, wbx;  
  unsigned char *ucip, *ucp;
  unsigned short ush, *ushp;
  unsigned int   ui, *uip;
  double  **dpp, *dp;
  int       i, j;
  static char    form1p[] = "UNEXPECTED case STOP bytesPerPixel=%d; expected 1,2,4 \n";
  static char    prognamp[] = "utiImage2DtoTnsrDbl";

  firstp = utiByteOrderGet();
  if(*firstp == 0) machineByteOrder = myMSBFirst;
  else             machineByteOrder = myLSBFirst;

  eByteOrder = (ep->im0D).byteOrder;
  if(machineByteOrder == eByteOrder) needSwap = 0;

  bytesPerPixel = (ep->im0D).bytesPerPixel;
  if(bytesPerPixel == 1) needSwap = 0;

  if(bytesPerPixel == 1 || bytesPerPixel == 2 ||bytesPerPixel == 4) oknbp = 1;
  if(!oknbp)
  { myErr2(myILLEGAL_CASE, stderr, srcfilenamp, prognamp, form1p, bytesPerPixel);
  }

  w = ep->wpx;  h = ep->hpx;  wbz = ep->wbz;  wbx = w * bytesPerPixel;
  ucip = (unsigned char *)ep->dataV.p;

  if(tp->p == NULL || tp->ev1z < h  || tp->ev2z < w)
  { tnsr2dblPAlloc(tp, h, w);
  }
  tp->ev1x = 0;  tp->ev2x = 0;

  dpp = tp->pp;
  for(i = 0;  i < h;  i++, ucip += wbz, dpp++)
  { ucp = ucip;
    dp = *dpp;
         if(bytesPerPixel == 1)
    { for(j = 0;  j < w;  j++){ ush = (unsigned short)*ucp++;  *dp++ = (double)ush;}
    }
    else if(bytesPerPixel == 2)
    { if(needSwap == 0)
      { for(j = 0;  j < w;  j++, ucp += 2)
        { ushp = (unsigned short *)ucp;  ush = *ushp;  *dp++ = (double)ush;}
      }
      else
      { for(j = 0;  j < w;  j++, ucp += 2)
        { ush = utiSwapB2ushort(ucp);  *dp++ = (double)ush;}
      }
    }
    else if(bytesPerPixel == 4)
    { if(needSwap == 0)
      { for(j = 0;  j < w;  j++, ucp += 4)
        { uip = (unsigned int *)ucp;  ui = *uip;  *dp++ = (double)ui;}
      }
      else
      { for(j = 0;  j < w;  j++, ucp += 4)
        { ui = utiSwapB4uint(ucp);  *dp++ = (double)ui;}
      }
    }
  }

  tp->ev1x = h;  tp->ev2x = w;
  return;
}
/******************************************************************************/
/******************************************************************************/
