/*  ../libmy/utiImage2D.c                                                     */
/*  Mennessier Gerard                   20031007                              */
/*  Last revised M.G.                   20040503                              */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utiVecChr.h"
#include  "utiByteOrderENDIAN.def.h"

#include  "utiImage0D.h"
#include  "utiImage2D.h"

static char    srcfilenamp[] = "utiImage2D";
/******************************************************************************/
/*  void  utiImage2DdataPVecAlloc(utiImage2D *ep, size_t nz)                  */
/*                                                                            */
/******************************************************************************/
void      utiImage2DdataPVecAlloc(utiImage2D *ep, size_t nz)
{ utiImage0D     im0D;
  chrVec        *dataVp;

  im0D = ep->im0D;
  dataVp = &(ep->dataV);
  chrPVecAlloc(dataVp, nz);
  return;
}
/******************************************************************************/
/*  void  utiImage2DdataPVecRealloc(utiImage2D *ep,                           */
/*                                             size_t neednz, size_t incrnz)  */
/*                                                                            */
/******************************************************************************/
void      utiImage2DdataPVecRealloc(utiImage2D *ep, size_t neednz, size_t incrnz)
{ chrVec        *dataVp;

  dataVp = &(ep->dataV);
  chrPVecRealloc(dataVp, neednz, incrnz);
  return;
}
/******************************************************************************/
/*  void  utiImage2DdataPVecFree(utiImage2D *ep)                              */
/*                                                                            */
/******************************************************************************/
void      utiImage2DdataPVecFree(utiImage2D *ep)
{
  chrPVecFree(&(ep->dataV) );
  return;
}
/******************************************************************************/
/*  void  utiImage2DZero(utiImage2D *ep)                                      */
/*                                                                            */
/******************************************************************************/
void      utiImage2DZero(utiImage2D *ep)
{

  (ep->im0D).bytesPerPixel = 0;
  (ep->dataV).x = 0;
  ep->wpx = 0;  ep->wbz = 0;
  ep->hpx = 0;
  ep->graymin = 1;  ep->graymax = 0;
  return;
}
/******************************************************************************/
/*  void  utiImage2DdataBytesPrint(FILE *streamp, utiImage2D *ep)             */
/*                                                                            */
/******************************************************************************/
void      utiImage2DdataBytesPrint(FILE *streamp, utiImage2D *ep)
{ int       m1, m1x, m1z, m2, m2x, ib, ibx;
  utiImage0D     im0D;
  unsigned char *cp, *cip;
  static char    form1p[] = "%s::%s  BEGIN pointer=%p\n";
  static char    form2p[] = "%s::%s  END \n";
  static char    form5p[] = "  utiImage2D line index=%d \n";
  static char    form6p[] = " %2x";
  static char    form7p[] = ",";
  static char    prognamp[] = "utiImage2DdataBytesPrint";

  fPrintF(streamp, form1p, srcfilenamp,prognamp, (void*)ep);
  im0D = ep->im0D;
  ibx = im0D.bytesPerPixel; 
  m1x = ep->wpx;  m1z = ep->wbz;
  m2x = ep->hpx;
  cip = (unsigned char *)( (ep->dataV).p);

  for(m2 = 0;  m2 < m2x;  m2++, cip += m1z)
  { fPrintF(streamp, form5p, m2);
    cp = cip;
    for(m1 = 0;  m1 < m1x;  m1++)
    { for(ib = 0;  ib < ibx;  ib++) fPrintF(streamp, form6p, *cp++);
      fPrintF(streamp, form7p);
    }
    fPrintF(streamp, "\n");
  }
  fPrintF(streamp, form2p, srcfilenamp,prognamp);
  return;
}
/******************************************************************************/
/*  void  utiImage2DPrint(FILE *streamp, utiImage2D *ep)                      */
/*                                                                            */
/******************************************************************************/
void      utiImage2DPrint(FILE *streamp, utiImage2D *ep)
{ 
  static char    form1p[] = "%s::%s  BEGIN pointer=%p\n";
  static char    form2p[] = "%s::%s  END\n";
  static char    form3p[] = "  dataV.p=%p  dataV.z=%d  dataV.x=%d\n";
  static char    form4p[] = "  pixel width=%d, allocatedByteWidth=%d\n";
  static char    form5p[] = "  pixel height=%d\n";
  static char    form7p[] = "  gray min value=%d, gray max value=%d \n";
  static char    prognamp[] = "utiImage2DPrint";

  fPrintF(streamp, form1p, srcfilenamp,prognamp, (void*)ep);
  fPrintF(streamp, "  ");  utiImage0DPrint(streamp, &(ep->im0D));
  fPrintF(streamp, form3p, (ep->dataV).p, (ep->dataV).z, (ep->dataV).x);
  fPrintF(streamp, form4p, ep->wpx, ep->wbz);
  fPrintF(streamp, form5p, ep->hpx);
  fPrintF(streamp, form7p, ep->graymin, ep->graymax);
  fPrintF(streamp, form2p, srcfilenamp,prognamp);
  return;
}
/******************************************************************************/
/*  void  utiImage2DSetAll(utiImage2D *ep, unsigned char *inkip)              */
/*                                                                            */
/*  set all pixels of the image to value pointed by inkip                     */
/*  it is NOT checked that inkip points to correct bytesPerPixel number       */
/******************************************************************************/
void      utiImage2DSetAll(utiImage2D *ep, unsigned char *inkip)
{ int       ib, ibx;
  int       m1, m1x, m1z, m2, m2x;
  unsigned char *cp, *cip;
  unsigned char *inkp;

  ibx = (ep->im0D).bytesPerPixel;
  m1x = ep->wpx;  m1z = ep->wbz;
  m2x = ep->hpx;
  cip = (unsigned char *)((ep->dataV).p);
  for(m2 = 0;  m2 < m2x;  m2++, cip += m1z)
  { cp = cip;
    for(m1 = 0;  m1 < m1x;  m1++)
    { inkp = inkip;
      for(ib = 0;  ib < ibx;  ib++) *cp++ = *inkp++;
    }
  }
  (ep->dataV).x = m1z * m2x;
  return;
}
/******************************************************************************/
/*  void  utiImage2DSet1(utiImage2D *ep, unsigned char *inkip, int ix,int iy) */
/*                                                                            */
/*  set 1 pixel of the image, coord (ix,iy) to value pointed by inkip         */
/*  it is NOT checked that inkip points to correct bytesPerPixel number       */
/******************************************************************************/
void      utiImage2DSet1(utiImage2D *ep, unsigned char *inkip, int ix,int iy)
{ int       ib, ibx;
  int       m1z;
  unsigned char *cp;
                                             /**  case OUT of image : do nothing  **/
  if(ix < 0)        return;
  if(ix >= ep->wpx) return;
  if(iy < 0)        return;
  if(iy >= ep->hpx) return;
                                                                 /**  normal case **/
  ibx = (ep->im0D).bytesPerPixel;
  m1z = ep->wbz;
  cp = (unsigned char *)((ep->dataV).p);
  cp += iy*m1z + ix*ibx;
  for(ib = 0;  ib < ibx;  ib++) *cp++ = *inkip++;
  return;
}
/******************************************************************************/
/*  void  utiImage2DMinMax(utiImage2D *ep)                                    */
/*                                                                            */
/******************************************************************************/
void      utiImage2DMinMax(utiImage2D *ep)
{ utiImage0D     im0D;
  int       m1, m1x, m1z, m2, m2x;
  int       bytesPerPixel;
  int       byteOrder;
  unsigned char *cp, *cip;
  unsigned short ugrayvalue = 0;
  unsigned int   grayvalue;
  unsigned int   graymin = 1, graymax = 0;
  short     boolPrintDebug = 1;
  static char    form3p[] = "%s::%s  ERROR NOT a GRAY 16, bytesPerPixel=%d\n";
  static char    form4p[] = "%s::%s  graymin=X%X =d %d, graymax=X%X =d %d\n";
  static char    prognamp[] = "utiImage2DMinMax";

  im0D = ep->im0D;
  bytesPerPixel = im0D.bytesPerPixel;
  if(bytesPerPixel != 2)
  { fPrintF(stderr, form3p, srcfilenamp,prognamp, bytesPerPixel);
    ep->graymin = 1;  ep->graymax = 0;
    return;
  }

  byteOrder = im0D.byteOrder;
  cip = (unsigned char *)( (ep->dataV).p);
  my2UBYTEStoUSHORT(byteOrder, cip, ugrayvalue);
  graymin = graymax = ugrayvalue;

  m1x = ep->wpx;  m1z = ep->wbz;
  m2x = ep->hpx;
  for(m2 = 0;  m2 < m2x;  m2++, cip += m1z)
  { cp = cip;
    for(m1 = 0;  m1 < m1x;  m1++, cp++, cp++)
    { my2UBYTEStoUSHORT(byteOrder, cp, ugrayvalue);
      grayvalue = ugrayvalue;
           if(grayvalue > graymax) graymax = grayvalue;
      else if(grayvalue < graymin) graymin = grayvalue;
    }
  }
if(boolPrintDebug)
{ fPrintF(stderr, form4p, srcfilenamp,prognamp, graymin, graymin, graymax, graymax);
}
  ep->graymin = graymin;  ep->graymax = graymax;
  return;
}
/******************************************************************************/
/******************************************************************************/
