/*  URMAE/orientHaut/linear4.GL.V4/pallidus.geom.c                            */
/*  Mennessier Gerard                 20000414                                */
/*  Last Revised : G.M.               20040210                                */

#include  <stddef.h>
#include  <math.h>

#include  "utiMath.constant.def.h"

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utistdStr.h"
#include  "utiMem.h"
#include  "utiStr.h"
#include  "utiAlloc.h"
#include  "utiVecInt.h"
#include  "utiVecDbl.h"
#include  "utiBookChr.h"
#include  "utiLec.h"
#include  "utiClassRangeDbl.h"
#include  "utiClassRangeInt.h"

#include  "pallidus.geom.globalloc.h"
#include  "pallidus.geom.h"

/**
      S  scanner Frame
      ST scanner Translated Frame, same origin as electrode Frame
      E  electrode Frame
**/
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      pallidusInitGeom(char *pallidusFileNamep)
{ double    psi;
  size_t    coordnx, nptot;
  char     *positionFileNamep;
  static    char    positionFileNameDefp[] = "position.data";
  static    char    prognamp[] = "pallidus.geom::pallidusInitGeom";

  positionFileNamep = (pallidusFileNamep)? pallidusFileNamep : positionFileNameDefp;
fPrintF(stderr, "%s: calling pallidusGetData with file %s\n",
                                                    prognamp, positionFileNamep);
  coordnx = pallidusGetData(positionFileNamep);
  nptot = coordnx/3;
fPrintF(stderr, "            pallidusGetData DONE, %d pallidus points\n", nptot);

  pallidusSsortPosData();
fPrintF(stderr, "            pallidusS sort PosData DONE\n");

  pallidusGetAffineRotation();
fPrintF(stderr, "            pallidus Get AffineRotation Matrices DONE\n");
  pallidusS2ERotateCheck();  pallidusE2SRotateCheck();
fPrintF(stderr, "            pallidus Rotate Check DONE\n");

  pallidusSsortedData2STBound();
fPrintF(stderr, "            pallidus sortedData to Translated Boundary DONE\n");

  psi = 0.0;
  pallidusSTBoundPsiRotate(psi);
fPrintF(stderr, "            pallidusSTBoundPsiRotate DONE\n");

  pallidusST2EBound();
fPrintF(stderr, "            pallidus Boundary rotated to ElectrodeFrame DONE\n");

  pallidusS2EDataAffineRotate();
fPrintF(stderr, "            pallidusS2EDataAffineRotate DONE\n");

  psi = 0.0;
  pallidusEBoundPsiRotate(psi);
fPrintF(stderr, "            pallidusEBoundPsiRotate DONE\n");

  return;
}
/******************************************************************************/
/*                                                                            */
/*  return the number of pallidus data (3 * number of points)                 */
/**        electrode data excluded                                            */
/******************************************************************************/
size_t    pallidusGetData(char *pallidusDataFileNamp)
{ FILE     *pallidusDataStreamp = NULL;
  chrBook  *chrBp;
  char     *cp,  *cip, *cfp, *ciip;
  char      space = ' ' ;
  ptrdiff_t     *ip;
  int       i, plotind, xip[3], isplot = 0, isinsul = 0;
  int       izprev, iz, npt, nplane;
  int      *iip;
  double    d;
  static    intVec    eltrdPinPosVi = {0,0,NULL};
  static    intVec    eltrdInsulatorPosVi = {0,0,NULL};
  static    char      form1p[] = "  Cannot open %s = %s file\n";
  static    char      form2p[] = "%s ERROR while reading data. \n"
                                                       "  line at%p; cip=%p, cp=%p\n"
                                                              "  Current line =%s\n";
  static    char      prognamp[] = "pallidus.geom::pallidusGetData";

  intPVecRealloc(&eltrdPinPosVi, 4, 2);          eltrdPinPosVi.x = 0;
  intPVecRealloc(&eltrdInsulatorPosVi, 4, 2);    eltrdInsulatorPosVi.x = 0;
  intPVecRealloc(&pallidusPosDataVi, 30, 4);      pallidusPosDataVi.x = 0;
  
  intPVecRealloc(&eltrdPinIndVi, 4, 2);           eltrdPinIndVi.x = 0;
  intPVecRealloc(&pallidusNptPlaneDataVi, 8, 2);  pallidusNptPlaneDataVi.x = 0;

  pallidusDataStreamp = fopen(pallidusDataFileNamp, "r");
  if(pallidusDataStreamp == NULL)
  { myErr1(-1,stderr,prognamp, form1p, "pallidus data file", pallidusDataFileNamp);
  }
  chrBp = lec3all(pallidusDataStreamp);
  fclose(pallidusDataStreamp);

  ciip = chrBp->bp;
  for(i = 0,  ip = chrBp->ip;  i < chrBp->ix;  i++, ip++)
  { cip = ciip + *ip;
    cp = strNormalizeSpace(cip, cip);
    if(cp < cip + 7) continue; 
    cfp = strChr(cip, space);
    if( memEq("plot", 4, cip, cfp - cip) )
    { plotind = (int)strtol(cfp, &cp, 0);
      isplot = 1;
      cip = cp;
      intVecInc1(&eltrdPinIndVi, plotind);
    }
    else if( memEq("i(", 2, cip, 2) )
    { isinsul = 1;
      cip = cfp;
    }
    else
    { isplot = 0;  isinsul = 0;
      cip = cfp;
    }
    xip[0] = (int)(0.5 + 10. * strtod(cip, &cp) );
    if(cp == cip) fPrintF(stderr, form2p, prognamp, ciip + *ip, cip, cp, ciip + *ip);
    cip = cp;  xip[1] = (int)(0.5 + 10. * strtod(cip, &cp) );
    if(cp == cip) fPrintF(stderr, form2p, prognamp, ciip + *ip, cip, cp, ciip + *ip);
    cip = cp;  xip[2] = (int)(0.5 + 10. * strtod(cip, &cp) );
    if(cp == cip) fPrintF(stderr, form2p, prognamp, ciip + *ip, cip, cp, ciip + *ip);
         if(isplot)  intVecIncN(&eltrdPinPosVi, xip, 3);
    else if(isinsul) intVecIncN(&eltrdInsulatorPosVi, xip, 3);
    else             intVecIncN(&pallidusPosDataVi, xip, 3);
  }
  chrBookFree(chrBp);

fPrintF(stderr,"%s  eltrdPinIndVi\n", prognamp);
intVecPrint(stderr, &eltrdPinIndVi);
fPrintF(stderr,"%s  eltrdPinPosVi\n", prognamp);
intVecPrint(stderr, &eltrdPinPosVi);
fPrintF(stderr,"%s  eltrdInsulatorPosVi\n", prognamp);
intVecPrint(stderr, &eltrdInsulatorPosVi);
fPrintF(stderr,"%s  pallidusPosDataVi\n", prognamp);
intVecPrint(stderr, &pallidusPosDataVi);

                                                     /** compute number of planes **/
  pallidusNptPlaneDataVi.x = 0;
  iip = pallidusPosDataVi.p + 2;
  izprev = *iip;  npt = 0;  nplane = 0;
  for(i = 0;  i < pallidusPosDataVi.x;  i+=3)
  { iz = *iip;   iip += 3;
    if(iz == izprev){ npt++;}
    else            { nplane++;  intVecInc1(&pallidusNptPlaneDataVi, npt);  npt = 1;}
    izprev = iz;
  }
  nplane++;  intVecInc1(&pallidusNptPlaneDataVi, npt);

fPrintF(stderr,"%s  pallidusNptPlaneDataVi  nplane=%d\n", prognamp, nplane);
intVecPrint(stderr, &pallidusNptPlaneDataVi);

  dblPVecRealloc(&eltrdSPosVd, eltrdPinPosVi.x, 2);
  eltrdSPosVd.x = 0;
  for(i = 0, iip = eltrdPinPosVi.p;  i < eltrdPinPosVi.x;  i++)
  { d = (double)*iip++;
    dblVecInc1(&eltrdSPosVd, 0.1 * d);
  }
  return  pallidusPosDataVi.x;
}
/******************************************************************************/
/*                                                                            */
/* Get Pallidus Boundary from Position Data                                   */
/******************************************************************************/
void      pallidusSsortPosData()
{ static    intVec    pallidusNptPlaneOffsetVi = {0,0,NULL};
  static    intVec    sortindexV = {0,0,NULL};
  static    intVec    zplanelevelV = {0,0,NULL};
  static    dblVec    anglV = {0,0,NULL};
  size_t    nplane;
  int       i, j, ifirst;
  int      *inp, *ioffp, ioff, ix, ns, n, npt, nptx;
  int      *zplanelevelp, *sortindexp, *izp, *ip;
  int       idminp[2], idmaxp[2];
  double    xb, yb, xm, ym;
  double    d2min, d2max, d2;
  double   *dp, *dip, *djp, *dnp, *anglp, d, dx, dy, x, y;

  static    char      prognamp[] = "pallidus.geom::pallidusSsortPosData";

  nplane = pallidusNptPlaneDataVi.x;

/*  zplanelevelp = intAlloc(nplane, prognamp); */
  intPVecRealloc(&zplanelevelV, nplane, 2);
  zplanelevelp = zplanelevelV.p;  zplanelevelV.x = 0;

  intPVecRealloc(&sortindexV, nplane, 2);
                                        sortindexp = sortindexV.p;  sortindexV.x = 0;
  intPVecRealloc(&pallidusNptPlaneOffsetVi, nplane, 2);
                                                      pallidusNptPlaneOffsetVi.x = 0;
                                        /** compute offsets and nptx = max(npt[]) **/
  nptx = 0;
  inp = pallidusNptPlaneDataVi.p;   
  ioffp = pallidusNptPlaneOffsetVi.p;
  ioff = 0;  intVecInc1(&pallidusNptPlaneOffsetVi, ioff);
  for(n = 1;  n < nplane;  n++)
  { npt = *inp++;  if(npt > nptx){ nptx = npt;}
    ioff += npt;  intVecInc1(&pallidusNptPlaneOffsetVi, ioff);
  }
fPrintF(stderr,"%s  pallidusNptPlaneOffsetVi\n", prognamp);
intVecPrint(stderr, &pallidusNptPlaneOffsetVi);

                                                       /** SORTING SCANNER PLANES **/

                                               /** get z coordinate of each plane **/
  inp = zplanelevelp;
  ioffp = pallidusNptPlaneOffsetVi.p;
  izp = pallidusPosDataVi.p + 2;                       /** +2 to get z coordinate **/
  for(n = 0;  n < nplane;  n++){ ioff = *ioffp++;  *inp++ = *(izp + 3 * ioff);}
                                                                  /** sort planes **/
  ranindInt(zplanelevelp, sortindexp, nplane);  sortindexV.x = nplane;

                                      /** store sorted number of points per plane **/
  intPVecRealloc(&pallidusNptPlaneVi, nplane, 0); pallidusNptPlaneVi.x = 0;
  for(n = 0;  n < nplane;  n++)
  { ns = sortindexp[n];  npt = *(pallidusNptPlaneDataVi.p + ns);
    intVecInc1(&pallidusNptPlaneVi, npt);
  }
fPrintF(stderr,"%s  pallidusNptPlaneVi  nplane=%d\n", prognamp, nplane);
intVecPrint(stderr, &pallidusNptPlaneVi);

                             /** store sorted pallidus data coordinates as double **/
  dblPVecRealloc(&pallidusSPosVd, pallidusPosDataVi.x, 2);
                                                  pallidusSPosVd.x = 0;
  for(n = 0;  n < nplane;  n++)
  { ns = sortindexp[n];
    ioff = *(pallidusNptPlaneOffsetVi.p + ns);
    npt = *(pallidusNptPlaneDataVi.p + ns);
    ip = pallidusPosDataVi.p + 3 * ioff;        ix = 3 * npt;
    for(i = 0;  i < ix;  i++)
    { d = (double)(*ip++);  dblVecInc1(&pallidusSPosVd, 0.1 * d);
    }
  }
fPrintF(stderr,"     plane sorted position data\n");
dblVecPrint(stderr, &pallidusSPosVd);
/*  free(zplanelevelp);  */

                                                 /** SORTING POINTS IN EACH PLANE **/
  intPVecRealloc(&sortindexV, nptx, 2);  sortindexp = sortindexV.p;
  dblPVecRealloc(&anglV, nptx, 2);
  dblPVecRealloc(&pallidusSsortedPosVd, pallidusSPosVd.x, 2);
                                                  pallidusSsortedPosVd.x = 0;
  inp = pallidusNptPlaneVi.p;
  npt = *inp;
  dnp = pallidusSPosVd.p;
  for(n = 0;  n < nplane;  n++)
  { npt = *inp++;
                                        /** compute barycenter, skip z coordinate **/
    xb = yb = 0.0;
    for(i = 0, dip = dnp;  i < npt;  i++){ xb += *dip++;  yb += *dip;  dip += 2;}
    xb = xb/npt;  yb = yb/npt;

                           /** compute greatest and smallest distances in z plane **/
    dip = dnp;  djp = dnp + 3;
    dx = *djp++ - *dip++;  dy = *djp - *dip;
    d2min = d2max = dx*dx + dy*dy;
    idminp[0] = 0;  idminp[1] = 1;  idmaxp[0] = 0;  idmaxp[1] = 1;
    for(i = 0;  i < npt -1;  i++)
    { dip = dnp + 3*i;  x = *dip++;  y = *dip;
      for(j = i+1;  j < npt;  j++)
      { djp = dnp + 3*j;
        dx = *djp++ - x;  dy = *djp - y;
        d2 = dx*dx + dy*dy;
        if     (d2 < d2min){ d2min = d2;  idminp[0] = i;  idminp[1] = j;}
        else if(d2 > d2max){ d2max = d2;  idmaxp[0] = i;  idmaxp[1] = j;}
      }
    }
                             /** ifirst = index of distant points with smallest y **/
    dip = dnp + 3*idmaxp[0] +1;  djp = dnp + 3*idmaxp[1] +1;
    ifirst = (*dip <= *djp)?  idmaxp[0] : idmaxp[1];
fPrintF(stderr,"  plane number %d, ifirst = %d\n", n, ifirst);
                                             /** compute middle of nearest points **/
    dip = dnp + 3 * idminp[0];  djp = dnp + 3 * idminp[1];
    xm = 0.5 * (*dip++ + *djp++);  ym = 0.5 * (*dip++ + *djp++);
                                                             /** WEIGHTED AVERAGE **/
    xb = 0.25 * (xb + 3. * xm);  yb = 0.25 * (yb + 3. * ym);

                                                               /** compute angles **/
    for(i = 0, dp = dnp, anglp = anglV.p;  i <  npt;  i++)
    { dx = *dp++  - xb;  dy = *dp++  - yb;  dp++; 
      *anglp++ = atan2(dy, dx);
    }

    ranindDbl(anglV.p, sortindexp, npt);  sortindexV.x = npt;
fPrintF(stderr,"     angle sort indices\n");
intVecPrint(stderr, &sortindexV);
    pallidusCircPerm(sortindexp, npt, ifirst);
fPrintF(stderr,"     angle sort indices after permutation, ifirst = %d\n", ifirst);
intVecPrint(stderr, &sortindexV);

    for(i = 0;  i < npt;  i++)
    { ns = sortindexp[i];  dp = dnp + 3 * ns;
      dblVecIncN(&pallidusSsortedPosVd, dp, 3);
    }
    dnp += 3 * npt;
  }
fPrintF(stderr,"     angle sorted position data\n");
dblVecPrint(stderr, &pallidusSsortedPosVd);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  from electrode position data : eltrdPinIndVi, eltrdSPosVd                 */
/*  get translation and rotation parameters (from Scanner to electrode)       */
/******************************************************************************/
void      pallidusGetAffineRotation()
{ double    x1d, x2d, x3d, xn;
  double   *dip, *dfp;
  double    ct, st,  cph, sph;
  double    x1OE, x2OE, x3OE;
  int      *ip, i1, i2, i, wi, n;
  static    char    prognamp[] = "pallidus.geom::pallidusGetAffineRotation";

  n = eltrdPinIndVi.x;
  i1 = 0;  i2 = n - 1;
  dip = eltrdSPosVd.p;   dfp = dip + 3 * i2;
                                   /** orientation up z, from pin 0  toward pin 3 **/
  x1d = *(dfp)     - *(dip);
  x2d = *(dfp + 1) - *(dip + 1);
  x3d = *(dfp + 2) - *(dip + 2);
  x1d = -x1d;  x2d = -x2d;  x3d = -x3d;                  /** opposite orientation **/
  xn = sqrt( x1d*x1d + x2d*x2d + x3d*x3d);
  x1d = x1d/xn;  x2d = x2d/xn;  x3d = x3d/xn;
  ct = x3d;  st = sqrt(x1d*x1d + x2d*x2d);
  ctheta = ct;  stheta = st;
  thetaRad = atan2(stheta, ctheta);  thetaDeg = thetaRad * myRAD2DEG;
  cph = x1d/st;  sph = x2d/st;
  cphi =cph;  sphi = sph;
  phiRad = atan2(sphi, cphi);  phiDeg = phiRad * myRAD2DEG;
  x1OE = x2OE = x3OE = 0.0;  wi = 0;
  for(i = 0,  dip = eltrdSPosVd.p,  ip = eltrdPinIndVi.p;  i < n;  i++)
  { x1OE += *dip++;  x2OE += *dip++;  x3OE += *dip++;
    wi += 2 * (*ip++) -3;
  }
  x1OE = x1OE/n;  x2OE = x2OE/n;  x3OE = x3OE/n;
  x1OE += (wi * x1d)/n;  x2OE += (wi * x2d)/n;  x3OE += (wi * x3d)/n;
fPrintF(stderr,"%s  electrode vector(x1,x2,x3) = (%f, %f, %f)\n",
                                                              prognamp, x1d,x2d,x3d);
fPrintF(stderr,"%s  (ct,st)=(%f, %f), (cphi,sphi)=(%f, %f) \n",
                                                            prognamp, ct,st,cph,sph);
fPrintF(stderr,"     Coord of Electrode Center in Scanner frame = %f, %f, %f \n",
                                                                   x1OE, x2OE, x3OE);
  eltrdSVecp[0] = x1d;   eltrdSVecp[1] = x2d;   eltrdSVecp[2] = x3d;
  eltrdSOrip[0] = x1OE;  eltrdSOrip[1] = x2OE;  eltrdSOrip[2] = x3OE;
  S2E0p[0] = ct*cph;  S2E0p[1] = ct*sph;  S2E0p[2] = -st;
  S2E1p[0] =   -sph;  S2E1p[1] =    cph;  S2E1p[2] = 0.;
  S2E2p[0] = st*cph;  S2E2p[1] = st*sph;  S2E2p[2] =  ct;
fPrintF(stderr,"%s    Rotation Matrix ST -> E\n", prognamp);
fPrintF(stderr,"     %8f  %8f  %8f\n", S2E0p[0],S2E0p[1],S2E0p[2]);
fPrintF(stderr,"     %8f  %8f  %8f\n", S2E1p[0],S2E1p[1],S2E1p[2]);
fPrintF(stderr,"     %8f  %8f  %8f\n", S2E2p[0],S2E2p[1],S2E2p[2]);

  E2ST0p[0] = S2E0p[0];  E2ST0p[1] = S2E1p[0];  E2ST0p[2] = S2E2p[0];
  E2ST1p[0] = S2E0p[1];  E2ST1p[1] = S2E1p[1];  E2ST1p[2] = S2E2p[1];
  E2ST2p[0] = S2E0p[2];  E2ST2p[1] = S2E1p[2];  E2ST2p[2] = S2E2p[2];
fPrintF(stderr,"%s     Rotation Matrix E -> ST\n", prognamp);
fPrintF(stderr,"     %8f  %8f  %8f\n", E2ST0p[0],E2ST0p[1],E2ST0p[2]);
fPrintF(stderr,"     %8f  %8f  %8f\n", E2ST1p[0],E2ST1p[1],E2ST1p[2]);
fPrintF(stderr,"     %8f  %8f  %8f\n", E2ST2p[0],E2ST2p[1],E2ST2p[2]);

  return;
}
/******************************************************************************/
/*                                                                            */
/* given the 3 point coordinates sxip[3], in Scanner frame,                   */
/*   get the 3 point coordinates exfp[3], in Electrode frame.                 */
/* WARNING exfp should be already allocated                                   */
/******************************************************************************/
void      pallidusS2EAffineRotate1pt(double *exfp, double *sxip)
{ double    xp[3];

  xp[0] = *sxip++  - eltrdSOrip[0];
  xp[1] = *sxip++  - eltrdSOrip[1];
  xp[2] = *sxip    - eltrdSOrip[2];

  *exfp = S2E0p[0] * xp[0] + S2E0p[1] * xp[1] + S2E0p[2] * xp[2];
  exfp++;
  *exfp = S2E1p[0] * xp[0] + S2E1p[1] * xp[1] + S2E1p[2] * xp[2];
  exfp++;
  *exfp = S2E2p[0] * xp[0] + S2E2p[1] * xp[1] + S2E2p[2] * xp[2];
  return;
}
/******************************************************************************/
/*                                                                            */
/* given the 3 vector coordinates sxip[3], in Scanner frame,                  */
/*   get the 3 vector coordinates exfp[3], in Electrode frame.                */
/* WARNING exfp should be already allocated                                   */
/******************************************************************************/
void      pallidusST2ERotate1pt(double *exfp, double *sxip)
{ double    xp[3];

  xp[0] = *sxip++;
  xp[1] = *sxip++;
  xp[2] = *sxip;

  *exfp = S2E0p[0] * xp[0] + S2E0p[1] * xp[1] + S2E0p[2] * xp[2];
  exfp++;
  *exfp = S2E1p[0] * xp[0] + S2E1p[1] * xp[1] + S2E1p[2] * xp[2];
  exfp++;
  *exfp = S2E2p[0] * xp[0] + S2E2p[1] * xp[1] + S2E2p[2] * xp[2];
  return;
}
/******************************************************************************/
/*                                                                            */
/* given the 3 point coordinates sxip[3], in Electrode frame,                 */
/*   get the 3 point coordinates exfp[3], in Scanner Translated frame.        */
/*                                         (i.e. same origin for both frames) */
/* WARNING exfp should be already allocated                                   */
/******************************************************************************/
void      pallidusE2SAffineRotate1pt(double *exfp, double *sxip)
{ double    xp[3];

  xp[0] = *sxip++;
  xp[1] = *sxip++;
  xp[2] = *sxip;

  *exfp = E2ST0p[0] * xp[0] + E2ST0p[1] * xp[1] + E2ST0p[2] * xp[2];
  *exfp += eltrdSOrip[0];
  exfp++;
  *exfp = E2ST1p[0] * xp[0] + E2ST1p[1] * xp[1] + E2ST1p[2] * xp[2];
  *exfp += eltrdSOrip[1];
  exfp++;
  *exfp = E2ST2p[0] * xp[0] + E2ST2p[1] * xp[1] + E2ST2p[2] * xp[2];
  *exfp += eltrdSOrip[2];
  return;
}

/******************************************************************************/
/*                                                                            */
/* given the 3 vector coordinates sxip[3], in Electrode frame,                */
/*   get the 3 vector coordinates exfp[3], in Scanner Translated frame.       */
/*                                         (i.e. same origin for both frames) */
/* WARNING exfp should be already allocated                                   */
/******************************************************************************/
void      pallidusE2STRotate1pt(double *exfp, double *sxip)
{ double    xp[3];

  xp[0] = *sxip++;
  xp[1] = *sxip++;
  xp[2] = *sxip;

  *exfp = E2ST0p[0] * xp[0] + E2ST0p[1] * xp[1] + E2ST0p[2] * xp[2];
  exfp++;
  *exfp = E2ST1p[0] * xp[0] + E2ST1p[1] * xp[1] + E2ST1p[2] * xp[2];
  exfp++;
  *exfp = E2ST2p[0] * xp[0] + E2ST2p[1] * xp[1] + E2ST2p[2] * xp[2];
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      pallidusS2ERotateCheck()
{ double    eltrdVecEp[3], xfp[3], *xip;
  int       i;
  static    char    prognamp[] = "pallidus.geom::pallidusS2ERotateCheck";

  pallidusST2ERotate1pt(eltrdVecEp, eltrdSVecp);
fPrintF(stderr,"%s  CHECK eltrdVecEp = (%f, %f, %f)\n", 
                              prognamp, eltrdVecEp[0], eltrdVecEp[1], eltrdVecEp[2]);

  pallidusS2EAffineRotate1pt(xfp, eltrdSOrip);
fPrintF(stderr,"%s  CHECK eltrdOriginEp = (%f, %f, %f)\n", 
                                                   prognamp, xfp[0], xfp[1], xfp[2]);
  for(i = 0,  xip = eltrdSPosVd.p;  i < eltrdPinIndVi.x;  i++, xip += 3)
  { pallidusS2EAffineRotate1pt(xfp, xip);
fPrintF(stderr,"%s  CHECK eltrdPin %d = (%f, %f, %f)\n", 
                           prognamp, *(eltrdPinIndVi.p + i), xfp[0], xfp[1], xfp[2]);
  }
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      pallidusE2SRotateCheck()
{ double    eltrdVecEp[3] = {0.0, 0.0, 1.0},  eltrdOriEp[3] = {0.0, 0.0, 0.0};
  double    xfp[3], *xip;
  int       i;
  static    char    prognamp[] = "pallidus.geom::pallidusE2SRotateCheck";

  pallidusE2STRotate1pt(xfp, eltrdVecEp);
fPrintF(stderr,"%s  CHECK eltrdVecSTp = (%f, %f, %f)\n", 
                                                   prognamp, xfp[0], xfp[1], xfp[2]);

  pallidusE2SAffineRotate1pt(xfp, eltrdOriEp);
fPrintF(stderr,"%s  CHECK eltrdOriginSp = (%f, %f, %f)\n", 
                                                   prognamp, xfp[0], xfp[1], xfp[2]);
  xip = eltrdOriEp;
  for(i = 0;  i <= +3;  i++)
  { xip[2] = (double)(-3 + 2*i);
    pallidusE2STRotate1pt(xfp, xip);
fPrintF(stderr,"%s  CHECK eltrdPin %d = (%f, %f, %f) in ST\n", 
                                                prognamp, i, xfp[0], xfp[1], xfp[2]);
    pallidusE2SAffineRotate1pt(xfp, xip);
fPrintF(stderr,"%s  CHECK eltrdPin %d = (%f, %f, %f) in S\n", 
                                                prognamp, i, xfp[0], xfp[1], xfp[2]);
  }
  return;
}
/******************************************************************************/
/*                                                                            */
/* Get Pallidus Boundary in the Scanner Translated ST frame                   */
/*   from Scanner S frame sorted data                                         */
/******************************************************************************/
void      pallidusSsortedData2STBound()
{ int       i, ix, i3x;
  double   *dip, dp[3];
  static    char    prognamp[] = "pallidus.geom::pallidusSsortedData2STBound";

  dblPVecRealloc(&pallidusSTBoundVd, pallidusSsortedPosVd.x, 2);
                                                             pallidusSTBoundVd.x = 0;
  i3x = pallidusSsortedPosVd.x;  ix = i3x/3;
  dip = pallidusSsortedPosVd.p;
  for(i = 1;  i <= ix;  i++)
  { dp[0] = *dip++ - eltrdSOrip[0];  dp[1] = *dip++ - eltrdSOrip[1];
    dp[2] = *dip++ - eltrdSOrip[2];
    dblVecIncN(&pallidusSTBoundVd, dp, 3);
  }
fPrintF(stderr,"\n");
fPrintF(stderr,"%s  pallidus sorted data ScannerFrame\n", prognamp);
dblVecPrint(stderr, &pallidusSsortedPosVd);
fPrintF(stderr,"\n");
fPrintF(stderr,"%s  pallidus Boundary ScannerTranslatedFrame\n", prognamp);
dblVecPrint(stderr, &pallidusSTBoundVd);
  return;
}
/******************************************************************************/
/*                                                                            */
/* Rotate ST frame Pallidus Boundary by angle psi                             */
/******************************************************************************/
void      pallidusSTBoundPsiRotate(double psi)
{ int       i, ix, i3x;
  double    c, s, xi1, xi2;
  double   *dip, dp[3];
  static    char    prognamp[] = "pallidus.geom::pallidusSTBoundPsiRotate";

  dblPVecRealloc(&pallidusSTBoundRotatedVd, pallidusSTBoundVd.x, 2);
                                                      pallidusSTBoundRotatedVd.x = 0;
  c = cos(psi);  s = sin(psi);
  i3x = pallidusSTBoundVd.x;  ix = i3x/3;
  dip = pallidusSTBoundVd.p;
  for(i = 1;  i <= ix;  i++)
  { xi1 = *dip++;  xi2 = *dip++;  dp[2] = *dip++;
    dp[0] = xi1 * c - xi2 * s;
    dp[1] = xi1 * s + xi2 * c;
    dblVecIncN(&pallidusSTBoundRotatedVd, dp, 3);
  }
fPrintF(stderr,"%s  ST pallidus data rotated, psi=%f rad\n", prognamp, psi);
dblVecPrint(stderr, &pallidusSTBoundRotatedVd);
  return;
}
/******************************************************************************/
/*                                                                            */
/* Get Pallidus Boundary in the Electrode E frame                             */
/*   from ST frame data                                                       */
/******************************************************************************/
void      pallidusST2EBound()
{ int       i, ix, i3x;
  double   *dip, dp[3], dfp[3];
  static    char    prognamp[] = "pallidus.geom::pallidusST2EBound";

  dblPVecRealloc(&pallidusEBoundVd, pallidusSTBoundVd.x, 2);
                                                              pallidusEBoundVd.x = 0;
  i3x = pallidusSTBoundVd.x;  ix = i3x/3;
  dip = pallidusSTBoundVd.p;
  for(i = 1;  i <= ix;  i++)
  { dp[0] = *dip++;  dp[1] = *dip++;  dp[2] = *dip++;
    pallidusST2ERotate1pt(dfp, dp);
    dblVecIncN(&pallidusEBoundVd, dfp, 3);
  }
fPrintF(stderr,"\n");
fPrintF(stderr,"%s  pallidus Boundary ElectrodeFrame\n", prognamp);
dblVecPrint(stderr, &pallidusEBoundVd);
  return;
}
/******************************************************************************/
/*                                                                            */
/* Get Pallidus Boundary in the Electrode E frame                             */
/*   from S frame data  (translation then rotation)                           */
/******************************************************************************/
void      pallidusS2EDataAffineRotate()
{ int       i, ix, i3x;
  double   *dip, dp[3];
  static    char    prognamp[] = "pallidus.geom::pallidusS2EDataAffineRotate";

  dblPVecRealloc(&pallidusEBoundVd, pallidusSsortedPosVd.x, 2);
                                                              pallidusEBoundVd.x = 0;
  i3x = pallidusSsortedPosVd.x;  ix = i3x/3;
  dip = pallidusSsortedPosVd.p;
  for(i = 1;  i <= ix;  i++, dip +=3)
  { pallidusS2EAffineRotate1pt(dp, dip);
    dblVecIncN(&pallidusEBoundVd, dp, 3);
  }
fPrintF(stderr,"\n");
fPrintF(stderr,"%s  pallidus sorted data ScannerFrame\n", prognamp);
dblVecPrint(stderr, &pallidusSsortedPosVd);
fPrintF(stderr,"%s  pallidus data ElectrodeFrame\n", prognamp);
dblVecPrint(stderr, &pallidusEBoundVd);
  return;
}
/******************************************************************************/
/*                                                                            */
/* Rotate E frame Pallidus Boundary by angle psi                              */
/******************************************************************************/
void      pallidusEBoundPsiRotate(double psi)
{ int       i, ix, i3x;
  double    c, s, xi1, xi2;
  double   *dip, dp[3];
  static    char    prognamp[] = "pallidus.geom::pallidusEBoundPsiRotate";

  dblPVecRealloc(&pallidusEBoundRotatedVd, pallidusEBoundVd.x, 2);
                                                       pallidusEBoundRotatedVd.x = 0;
  c = cos(psi);  s = sin(psi);
  i3x = pallidusEBoundVd.x;  ix = i3x/3;
  dip = pallidusEBoundVd.p;
  for(i = 1;  i <= ix;  i++)
  { xi1 = *dip++;  xi2 = *dip++;  dp[2] = *dip++;
    dp[0] = xi1 * c - xi2 * s;
    dp[1] = xi1 * s + xi2 * c;
    dblVecIncN(&pallidusEBoundRotatedVd, dp, 3);
  }
fPrintF(stderr,"%s  E pallidus data rotated, psi=%f rad\n", prognamp, psi);
dblVecPrint(stderr, &pallidusEBoundRotatedVd);
  return;
}
/******************************************************************************/
/*                                                                            */
/*  used to put in third position the fixed section variable (z role)         */
/*  while the 2 others play te role of coordinates in the section (x,y role)  */
/*  NO MORE NECESSARY after change in "levelcurve from Quadr"                 */
/******************************************************************************/
/*
void      pallidusEQuadr()
{ int       i, ix, i3x;
  double   *dip, dp[3];
  static    char    prognamp[] = "pallidus.geom::pallidusEQuadr";

  dblPVecRealloc(&pallidusEQuadrVd, pallidusEBoundVd.x, 2);
  pallidusEQuadrVd.x = 0;
  i3x = pallidusEBoundVd.x;  ix = i3x/3;
  dip = pallidusEBoundRotatedVd.p;
  for(i = 1;  i <= ix;  i++)
  { dp[0] = *dip++;  dp[1] = *dip++;  dp[2] = *dip++;
    dblVecIncN(&pallidusEQuadrVd, dp, 3);
  }
fPrintF(stderr,"%s  pallidus Quadr coordinates\n", prognamp);
dblVecPrint(stderr, &pallidusEQuadrVd);
  return;
}
*/
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      pallidusCircPerm(int *ip, int ix, int i0)
{ int      *jp, j, is;

  is = *ip;
  while(is != i0)
  { jp = ip;
    for(j = 0;  j < ix -1;  j++){ *jp = *(jp+1);  jp++;}
    *jp = is;
    is = *ip;
  }
  return;
}
/******************************************************************************/
/******************************************************************************/
