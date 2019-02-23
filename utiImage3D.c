/*  ../libmy/utiImage3D.c                                                     */
/*  Mennessier Gerard                   20031010                              */
/*  Last revised M.G.                   20040510                              */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utiVecChr.h"
#include  "utiVecUInt.h"
#include  "utiByteOrderENDIAN.def.h"

#include  "utiImage0D.h"
#include  "utiImage3D.h"

static char    srcfilenamp[] = "utiImage3D";
/******************************************************************************/
/*  void  utiImage3DdataPVecAlloc(utiImage3D *ep, size_t dataz, size_t zz)    */
/*                                                                            */
/*  allocate dataz bytes for image data                                       */
/*        and   zz (unsigned)int for both min and max values per plane        */
/*                                                                            */
/*  dataz needed size for dataV;                                              */
/*  zz    needed size for graymin2DV, graymax2DV                              */
/******************************************************************************/
void      utiImage3DdataPVecAlloc(utiImage3D *ep, size_t dataz, size_t zz)
{ chrVec        *dataVp;
  uintVec       *uiVp;

  dataVp = &(ep->dataV);
  chrPVecAlloc(dataVp, dataz);

  uiVp = &(ep->graymin2DV);
  uintPVecAlloc(uiVp, zz);
  uiVp = &(ep->graymax2DV);
  uintPVecAlloc(uiVp, zz);

  return;
}
/******************************************************************************/
/*  void  utiImage3DdataPVecRealloc(utiImage3D *ep,                           */
/*         size_t needataz, size_t dataincrnz, size_t needzz, size_t incrnzz) */
/*                                                                            */
/*  needataz needed size for dataV;                                           */
/*  needzz   needed size for graymin2DV, graymax2DV                           */
/******************************************************************************/
void      utiImage3DdataPVecRealloc(utiImage3D *ep,
                   size_t needataz, size_t dataincrnz, size_t needzz, size_t incrnzz)
{ chrVec        *dataVp;
  uintVec       *uiVp;

  dataVp = &(ep->dataV);
  chrPVecRealloc(dataVp, needataz, dataincrnz);

  uiVp = &(ep->graymin2DV);
  uintPVecRealloc(uiVp, needzz, incrnzz);
  uiVp = &(ep->graymax2DV);
  uintPVecRealloc(uiVp, needzz, incrnzz);

  return;
}
/******************************************************************************/
/*  void  utiImage3DdataPVecFree(utiImage3D *ep)                              */
/*                                                                            */
/******************************************************************************/
void      utiImage3DdataPVecFree(utiImage3D *ep)
{
  chrPVecFree(&(ep->dataV) );
  uintPVecFree(&(ep->graymin2DV) );
  uintPVecFree(&(ep->graymax2DV) );

  return;
}
/******************************************************************************/
/*  void  utiImage3DZero(utiImage3D *ep)                                      */
/*                                                                            */
/******************************************************************************/
void      utiImage3DZero(utiImage3D *ep)
{

  (ep->im0D).bytesPerPixel = 0;
  (ep->dataV).x = 0;
  ep->wpx = 0;  ep->wbz = 0;
  ep->hpx = 0;  ep->recbz = 0;  ep->recpx = 0;
  ep->zpx = 0;

  (ep->graymin2DV).x = 0;  (ep->graymax2DV).x = 0;
  ep->graymin3D = 1;  ep->graymax3D = 0;
  return;
}
/******************************************************************************/
/*  void  utiImage3DdataBytesPrint(FILE *streamp, utiImage3D *ep)             */
/*                                                                            */
/******************************************************************************/
void      utiImage3DdataBytesPrint(FILE *streamp, utiImage3D *ep)
{ int       m1, m1x, m1z, m2, m2x, m2z, m3, m3x, ib, ibx;
  utiImage0D     im0D;
  unsigned char *cp, *ci1p, *ci2p;
  static char    form1p[] = "%s::%s  BEGIN pointer=%p \n";
  static char    form2p[] = "%s::%s  END \n";
  static char    form4p[] = "    utiImage3D slice index=%d \n";
  static char    form5p[] = "    utiImage3D line index=%d \n";
  static char    form6p[] = " %2x";
  static char    form7p[] = ",";
  static char    prognamp[] = "utiImage3DdataBytesPrint";

  fPrintF(streamp, form1p, srcfilenamp,prognamp, (void*)ep);
  im0D = ep->im0D;
  ibx = im0D.bytesPerPixel; 
  m1x = ep->wpx;  m1z = ep->wbz;
  m2x = ep->hpx;  m2z = ep->recbz;
  m3x = ep->zpx;
  ci2p = (unsigned char *)( (ep->dataV).p);

  for(m3 = 0;  m3 < m3x;  m3++, ci2p += m2z)
  { fPrintF(streamp, form4p, m3);
    ci1p = ci2p;
    for(m2 = 0;  m2 < m2x;  m2++, ci1p += m1z)
    { fPrintF(streamp, form5p, m2);
      cp = ci1p;
      for(m1 = 0;  m1 < m1x;  m1++)
      { for(ib = 0;  ib < ibx;  ib++) fPrintF(streamp, form6p, *cp++);
        fPrintF(streamp, form7p);
      }
      fPrintF(streamp, "\n");
    }
  }
  fPrintF(streamp, form2p, srcfilenamp,prognamp);
  return;
}
/******************************************************************************/
/*  void  utiImage3DPrint(FILE *streamp, utiImage3D *ep)                      */
/*                                                                            */
/******************************************************************************/
void      utiImage3DPrint(FILE *streamp, utiImage3D *ep)
{ 
  static char    form1p[] = "%s::%s  BEGIN pointer=%p \n";
  static char    form2p[] = "%s::%s  END\n";
  static char    form3p[] = "  dataVp=%p  dataVz=%d  dataVx=%d\n";
  static char    form4p[] = "  pixel width=%d, allocatedByteWidth=%d\n";
  static char    form5p[] = "  pixel height=%d, rectanglePixelNumber=%d, allocatedByteRectangle=%d\n";
  static char    form6p[] = "  2D rectangle number==slice number = %d \n";
  static char    form7p[] = "  graymin2DVp=%p  graymin2DVz=%d  graymin2DVx=%d\n";
  static char    form8p[] = "  graymax2DVp=%p  graymax2DVz=%d  graymax2DVx=%d\n";
  static char    form9p[] = "  gray min3D value=%d, gray max3D value=%d \n";
  static char    prognamp[] = "utiImage3DPrint";

  fPrintF(streamp, form1p, srcfilenamp,prognamp, (void*)ep);
  utiImage0DPrint(streamp, &(ep->im0D));
  fPrintF(streamp, form3p, (ep->dataV).p, (ep->dataV).z, (ep->dataV).x);
  fPrintF(streamp, form4p, ep->wpx, ep->wbz);
  fPrintF(streamp, form5p, ep->hpx, ep->recpx, ep->recbz);
  fPrintF(streamp, form6p, ep->zpx);
  fPrintF(streamp, form7p, (ep->graymin2DV).p, (ep->graymin2DV).z, (ep->graymin2DV).x);
  fPrintF(streamp, form8p, (ep->graymax2DV).p, (ep->graymax2DV).z, (ep->graymax2DV).x);
  fPrintF(streamp, form9p, ep->graymin3D, ep->graymax3D);
  fPrintF(streamp, form2p, srcfilenamp,prognamp);
  return;
}
/******************************************************************************/
/*  void  utiImage3DGet2D(utiImage3D *e3p, utiImage2D *e2p, int iz)           */
/*                                                                            */
/*  set *e2p from slice/rectangle index iz of e3p                             */
/*  pixel data are NOT COPIED, e2.dataV just points into e3 data              */
/******************************************************************************/
void      utiImage3DGet2D(utiImage3D *e3p, utiImage2D *e2p, int iz)
{

  utiImage2DZero(e2p);
  utiImage0DCpy( &(e2p->im0D), e3p->im0D);
  e2p->wbz = e3p->wbz;  e2p->wpx = e3p->wpx;  e2p->hpx = e3p->hpx;
  (e2p->dataV).p = (e3p->dataV).p + iz * e3p->recbz;
  (e2p->dataV).x = e3p->recbz;
  (e2p->dataV).z = e3p->recbz;     /** NOT significant **/
  if( (e3p->graymin2DV).p != NULL) e2p->graymin = *( (e3p->graymin2DV).p + iz);
  else                             e2p->graymin = 1;
  if( (e3p->graymax2DV).p != NULL) e2p->graymax = *( (e3p->graymax2DV).p + iz);
  else                             e2p->graymax = 0;
  return;
}
/******************************************************************************/
/*  short  utiImage3DAdd2D(utiImage3D *e3p, utiImage2D *e2p)                  */
/*                                                                            */
/*  add a utiImage2D to a utiImage3D at highest index value                   */
/*  return OK = 1 if possible                                                 */
/*         OK = 0 if incompatible                                             */
/******************************************************************************/
short     utiImage3DAdd2D(utiImage3D *e3p, utiImage2D *e2p)
{ short     OK = 0, OKC;
  int       bytesPerPixelMin = 0;
  int       bytesPerPixel2D, bytesPerPixel3D;
  size_t    wpx, hpx;
  size_t    zpxi, zpxf, zpxfC;
  size_t    recbz, dataix, datafx;
  int       m1, m1x, m2, m2x, ib;
  char     *cf1p, *ci1p, *cfp, *cip;
  static char    form1p[] = "%s::%s \n";
  static char    form3p[] = "  width  of new image2D differs from 3D. Width  2D=%d, 3D=%d\n";
  static char    form4p[] = "  height of new image2D differs from 3D. Height 2D=%d, 3D=%d\n";
  static char    prognamp[] = "utiImage3DAdd2D";

  if(e2p == NULL) return OK;
  OKC = utiImage0DCompatible( &(e3p->im0D), &(e2p->im0D) );
  if(! OKC) return OK;
  if( (e2p->dataV).p == NULL) return OK;

  zpxi = e3p->zpx;
                     /**  Case FIRST image added. Choose values from that image2D **/
  if(zpxi == 0)
  { if( (e3p->im0D).bytesPerPixel == 0) (e3p->im0D).bytesPerPixel = (e2p->im0D).bytesPerPixel;
    e3p->wpx = e2p->wpx;  e3p->hpx = e2p->hpx;
    e3p->wbz = e2p->wbz;
    e3p->recpx = e3p->wpx * e3p->hpx;
    e3p->recbz = e3p->recpx * (e3p->im0D).bytesPerPixel;
  }

  wpx = e3p->wpx;
  if(e2p->wpx != wpx)
  { fPrintF(stderr, form1p, srcfilenamp,prognamp);  fPrintF(stderr, form3p, e2p->wpx, wpx);
    return OK;
  }
  hpx = e3p->hpx;
  if(e2p->hpx != hpx)
  { fPrintF(stderr, form1p, srcfilenamp,prognamp);  fPrintF(stderr, form4p, e2p->hpx, hpx);
    return OK;
  }

  OK = 1;
  bytesPerPixelMin = 1 + ( (e2p->im0D).depth -1)/8;
  bytesPerPixel2D = (e2p->im0D).bytesPerPixel;
  bytesPerPixel3D = (e3p->im0D).bytesPerPixel;
  zpxf = zpxi + 1;  zpxfC = zpxi;
  recbz = e3p->recbz;
  dataix = zpxi * recbz;
  datafx = dataix + recbz;
  chrPVecRealloc( &(e3p->dataV),  datafx, recbz);
  cf1p = (e3p->dataV).p + zpxfC * recbz;
  ci1p = (e2p->dataV).p;  cfp = cf1p;
  m1x = wpx;  m2x = hpx;
  for(m2 = 0;  m2 < m2x;  m2++, ci1p += e2p->wbz, cf1p += e3p->wbz)
  { cip = ci1p;  cfp = cf1p;
    for(m1 = 0;  m1 < m1x;  m1++)
    { for(ib = 0;  ib < bytesPerPixelMin;  ib++) *cfp++ = *cip++;
      for(ib = 0;  ib < (bytesPerPixel3D - bytesPerPixelMin);  ib++) *cfp++ = 0;
      cip += (bytesPerPixel2D - bytesPerPixelMin);
    }
    for(ib = m1x * bytesPerPixel3D;  ib < e3p->wbz;  ib++) *cfp++ = 0;
  }
  for(ib = m2x * e3p->wbz;  ib < recbz;  ib++) *cfp++ = 0;

  e3p->zpx = zpxf;
  return OK;
}
/******************************************************************************/
/*  void  utiImage3DMinMax2D(utiImage3D *ep, int iz, unsigned int minmaxp[2]) */
/*                                                                            */
/*  compute min-max of slice indexC=iz, and store them into minmaxp           */
/******************************************************************************/
void      utiImage3DMinMax2D(utiImage3D *ep, int iz, unsigned int minmaxp[2])
{ utiImage2D     e2;

  utiImage2DZero(&e2);
  utiImage0DCpy( &(e2.im0D), ep->im0D );
  e2.wbz = ep->wbz;  e2.wpx = ep->wpx;  e2.hpx = ep->hpx;
  (e2.dataV).p = (ep->dataV).p + iz * ep->recbz;
  utiImage2DMinMax(&e2);
  minmaxp[0] = e2.graymin;  minmaxp[1] = e2.graymax;
  return;
}
/******************************************************************************/
/*  void  utiImage3DAdd2DMinMax(utiImage3D *e3p, utiImage2D *e2p)             */
/*                                                                            */
/*  add a utiImage2D min-max to a utiImage3D min-max at highest index value   */
/******************************************************************************/
void      utiImage3DAdd2DMinMax(utiImage3D *e3p, utiImage2D *e2p)
{ size_t    zpxfC, zpxf;

  zpxf = e3p->zpx;  zpxfC = zpxf - 1;
  uintPVecRealloc( &(e3p->graymin2DV), zpxf, 1);
  uintPVecRealloc( &(e3p->graymax2DV), zpxf, 1);
  if(e2p->graymin > e2p->graymax) utiImage2DMinMax(e2p);

  *( (e3p->graymin2DV).p + zpxfC) = e2p->graymin;
  if(e2p->graymin < e3p->graymin3D) e3p->graymin3D = e2p->graymin;
  (e3p->graymin2DV).x = zpxf;

  *( (e3p->graymax2DV).p + zpxfC) = e2p->graymax;
  if(e2p->graymax > e3p->graymax3D) e3p->graymax3D = e2p->graymax;
  (e3p->graymax2DV).x = zpxf;

  return;
}
/******************************************************************************/
/*  void  utiImage3DMinMax(utiImage3D *e3p)                                   */
/*                                                                            */
/*  set all slices and global min-max                                         */
/******************************************************************************/
void      utiImage3DMinMax(utiImage3D *e3p)
{ size_t    zpxf;
  int       m3;
  unsigned int   minmaxp[2], graymin, graymax;
  unsigned int  *uiminp, *uimaxp;

  zpxf = e3p->zpx;
  uintPVecRealloc( &(e3p->graymin2DV), zpxf, 1);
  uintPVecRealloc( &(e3p->graymax2DV), zpxf, 1);
  uiminp = (e3p->graymin2DV).p;  uimaxp = (e3p->graymax2DV).p;

  m3 = 0;
  utiImage3DMinMax2D(e3p, m3, minmaxp);
  graymin = minmaxp[0];  graymax = minmaxp[1];
  *uiminp++ = minmaxp[0];  *uimaxp++ = minmaxp[1];

  for(m3 = 1;  m3 < zpxf;  m3++)
  { utiImage3DMinMax2D(e3p, m3, minmaxp);
    *uiminp++ = minmaxp[0];  *uimaxp++ = minmaxp[1];
    if(minmaxp[0] < graymin) graymin = minmaxp[0];
    if(minmaxp[1] > graymax) graymax = minmaxp[1];
  }
  e3p->graymin3D = graymin;  e3p->graymax3D = graymax;
  (e3p->graymin2DV).x = (e3p->graymax2DV).x = zpxf;
  return; 
}
/******************************************************************************/
/******************************************************************************/
