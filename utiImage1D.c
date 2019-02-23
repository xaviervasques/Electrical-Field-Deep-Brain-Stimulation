/*  ../libmy/utiImage1D.c                                                     */
/*  Mennessier Gerard                   20031008                              */
/*  Last revised M.G.                   20040503                              */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utiVecChr.h"
#include  "utiByteOrderENDIAN.def.h"

#include  "utiImage0D.h"
#include  "utiImage1D.h"

static char    srcfilenamp[] = "utiImage1D";
/******************************************************************************/
/*  void  utiImage1DdataPVecAlloc(utiImage1D *ep, size_t nz)                  */
/*                                                                            */
/******************************************************************************/
void      utiImage1DdataPVecAlloc(utiImage1D *ep, size_t nz)
{ utiImage0D     im0D;
  chrVec        *dataVp;

  im0D = ep->im0D;
  dataVp = &(ep->dataV);
  chrPVecAlloc(dataVp, nz);
  return;
}
/******************************************************************************/
/*  void  utiImage1DdataPVecRealloc(utiImage1D *ep,                           */
/*                                             size_t neednz, size_t incrnz)  */
/*                                                                            */
/******************************************************************************/
void      utiImage1DdataPVecRealloc(utiImage1D *ep, size_t neednz, size_t incrnz)
{ chrVec        *dataVp;

  dataVp = &(ep->dataV);
  chrPVecRealloc(dataVp, neednz, incrnz);
  return;
}
/******************************************************************************/
/*  void  utiImage1DdataPVecFree(utiImage1D *ep)                              */
/*                                                                            */
/******************************************************************************/
void      utiImage1DdataPVecFree(utiImage1D *ep)
{
  chrPVecFree(&(ep->dataV) );
  return;
}
/******************************************************************************/
/*  void  utiImage1DdataBytesPrint(FILE *streamp, utiImage1D *ep)             */
/*                                                                            */
/******************************************************************************/
void      utiImage1DdataBytesPrint(FILE *streamp, utiImage1D *ep)
{ int       m1, m1x, ib, ibx;
  utiImage0D     im0D;
  unsigned char *cp, *cip;
  static char    form1p[] = "%s::%s  BEGIN pointer=%p\n";
  static char    form2p[] = "%s::%s  END \n";
  static char    form4p[] = " %2x";
  static char    form5p[] = ",";
  static char    prognamp[] = "utiImage1DdataBytesPrint";

  fPrintF(streamp, form1p, srcfilenamp,prognamp, (void*)ep);
  im0D = ep->im0D;
  ibx = im0D.bytesPerPixel; 
  m1x = ep->wpx;
  cip = (unsigned char *)( (ep->dataV).p);
  cp = cip;
  for(m1 = 0;  m1 < m1x;  m1++)
  { for(ib = 0;  ib < ibx;  ib++) fPrintF(streamp, form4p, *cp++);
    fPrintF(streamp, form5p);
  }
  fPrintF(streamp, "\n");
  fPrintF(streamp, form2p, srcfilenamp,prognamp);
  return;
}
/******************************************************************************/
/*  void  utiImage1DZero(utiImage1D *ep)                                      */
/*                                                                            */
/******************************************************************************/
void      utiImage1DZero(utiImage1D *ep)
{

  (ep->im0D).bytesPerPixel = 0;
  (ep->dataV).x = 0;
  ep->wpx = 0;
  ep->graymin = 1;  ep->graymax = 0;
  return;
}
/******************************************************************************/
/*  void  utiImage1DPrint(FILE *streamp, utiImage1D *ep)                      */
/*                                                                            */
/******************************************************************************/
void      utiImage1DPrint(FILE *streamp, utiImage1D *ep)
{ 
  static char    form1p[] = "%s::%s  BEGIN  pointer=%p\n";
  static char    form2p[] = "%s::%s  END\n";
  static char    form3p[] = "  dataVp=%p  dataVz=%d  dataVx=%d\n";
  static char    form4p[] = "  pixel width=%d ";
  static char    form5p[] = "  gray min value=%d, gray max value=%d \n";
  static char    prognamp[] = "utiImage1DPrint";

  fPrintF(streamp, form1p, srcfilenamp, prognamp, (void*)ep);
  fPrintF(streamp, "  ");  utiImage0DPrint(streamp, &(ep->im0D));
  fPrintF(streamp, form3p, (ep->dataV).p, (ep->dataV).z, (ep->dataV).x);
  fPrintF(streamp, form4p, ep->wpx);
  fPrintF(streamp, form5p, ep->graymin, ep->graymax);
  fPrintF(streamp, form2p, srcfilenamp, prognamp);
  return;
}
/******************************************************************************/
/*  void  utiImage1DMinMax(utiImage1D *ep)                                    */
/*                                                                            */
/******************************************************************************/
void      utiImage1DMinMax(utiImage1D *ep)
{ utiImage0D     im0D;
  int       m1, m1x;
  int       bytesPerPixel;
  int       byteOrder;
  unsigned char *cp, *cip;
  unsigned short ugrayvalue = 0;
  unsigned int   grayvalue;
  unsigned int   graymin = 1, graymax = 0;
  short     boolPrintDebug = 1;
  static char    form3p[] = "%s::%s  ERROR NOT a GRAY 16, bytesPerPixel=%d\n";
  static char    form4p[] = "%s::%s  graymin=X%X =d %d, graymax=X%X =d %d\n";
  static char    prognamp[] = "utiImage1DMinMax";

  im0D = ep->im0D;
  bytesPerPixel = im0D.bytesPerPixel;
  if(bytesPerPixel != 2)
  { fPrintF(stderr, form3p, srcfilenamp, prognamp, bytesPerPixel);
    ep->graymin = 1;  ep->graymax = 0;
    return;
  }

  byteOrder = im0D.byteOrder;
  cip = (unsigned char *)( (ep->dataV).p);
  my2UBYTEStoUSHORT(byteOrder, cip, ugrayvalue);
  graymin = graymax = ugrayvalue;

  m1x = ep->wpx;
  cp = cip;
  for(m1 = 0;  m1 < m1x;  m1++, cp++, cp++)
  { my2UBYTEStoUSHORT(byteOrder, cp, ugrayvalue);
    grayvalue = ugrayvalue;
         if(grayvalue > graymax) graymax = grayvalue;
    else if(grayvalue < graymin) graymin = grayvalue;
  }
if(boolPrintDebug)
{ fPrintF(stderr, form4p, srcfilenamp, prognamp, graymin, graymin, graymax, graymax);
}
  ep->graymin = graymin;  ep->graymax = graymax;
  return;
}
/******************************************************************************/
/******************************************************************************/
