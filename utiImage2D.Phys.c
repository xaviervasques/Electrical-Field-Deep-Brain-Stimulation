/*  ../libmy/utiImage2D.Phys.c                                                */
/*  Mennessier Gerard                   20030618                              */
/*  Last Revised : G.M.                 20040504                              */

#include  "utiImage2D.Phys.h"

/* static char    srcfilenamp[] = "utiImage2D.Phys"; */
/******************************************************************************/
/*  void utiImage2DPhysSetBoundary(utiImage2DPhys *physWinp,                  */
/*                                              double xnx[2], double ynx[2]) */
/*                                                                            */
/*  define rectangle window physical boundaries                               */
/*           x = {xnx[0] , xnx[1]}, y = {ynx[0] , ynx[1]}                     */
/******************************************************************************/
void      utiImage2DPhysSetBoundary(utiImage2DPhys *physWinp,
                                                        double xnx[2], double ynx[2])
{ 
  physWinp->xnxp[0] = xnx[0];  physWinp->xnxp[1] = xnx[1];
  physWinp->ynxp[0] = ynx[0];  physWinp->ynxp[1] = ynx[1];

  physWinp->deltax = xnx[1] - xnx[0];
  physWinp->deltay = ynx[1] - ynx[0];

  return;
}
/******************************************************************************/
/*  void utiImage2DPhysSetWH(utiImage2DPhys *physWinp, int w, int h)        */
/*                                                                            */
/*  define rectangle                                                          */
/*           x={xnx[0]-xnx[1]}, y={ynx[0]-ynx[1]}                             */
/*         to be mapped onto window,                                          */
/*           width =w (x direction),                                          */
/*           height=h (y direction)                                           */
/******************************************************************************/
void      utiImage2DPhysSetWH(utiImage2DPhys *physWinp, int w, int h)
{
  physWinp->pixdx  = physWinp->deltax / w;
  physWinp->scalex = w / physWinp->deltax;
  physWinp->pixdy  = physWinp->deltay / h;
  physWinp->scaley = h / physWinp->deltay;

  return;
}
/******************************************************************************/
/******************************************************************************/
