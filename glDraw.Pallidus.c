#include  <stdio.h>
#include  <string.h>
#include  <stdlib.h>
#include  <math.h>
#include  <GL/glut.h>

#include  "grille.glob.h"
#include  "glDraw.glob.h"
#include  "marchCube.glob.h"

#include  "glDraw.Pallidus.h"
#include  "grille.h"
#include  "marchCube.h"
#include  "funcTot.h"

static char    srcfilenamp[] = "glDrawPallidus";
static int     DEBUG = 0;
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      glDrawTriangles() 
{ int       ix, iy, iz;
  int       TailleYm1, TailleZm1;
  int       nfaces, nbTotFaces = 0;
  static char    form1p[] = "%s::%s ix=%4d \n";
  static char    form2p[] = "%s::%s ix=%4d, nbTotFaces=%d"  
                                                 "  DONE for TailleY=%d, TailleZ=%d\n";
  static char    form3p[] = "%s::%s DONE  nbTotFaces=%d \n\n";
  static char    prognamp[] = "glDrawTriangles";


  glNewList( NumList , GL_COMPILE_AND_EXECUTE );
  glBegin(GL_TRIANGLES);

                                   /** Taille. sommets ==> Taille. -1 intervalles **/
  TailleYm1 = TailleY - 1;
  TailleZm1 = TailleZ - 1;
  for(ix=0;  ix < TailleX - 1;  ix++)
  { 
fprintf(stderr, form1p, srcfilenamp, prognamp, ix);
    for (iy=0;  iy < TailleYm1;  iy++)
    { for (iz=0;  iz < TailleZm1;  iz++)
      { nfaces = glDrawCell(ix, iy, iz);
        nbTotFaces += nfaces;
      }
    }
fprintf(stderr, form2p, srcfilenamp, prognamp, ix, nbTotFaces, TailleY, TailleZ);
  }

  glEnd();
  glEndList( );
fprintf(stderr, form3p, srcfilenamp, prognamp, nbTotFaces);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       glDrawCell(int ix, int iy, int iz)
{ int       nbfaces;
  int       fce, is;
  double   *posdp;
  double    normaldp[3];
  double    norm, invNorm;
  double    elemCubeNormalp[3];
  float     posfp[3];
  float     normalfp[3];
  static char    form1p[] = "%s::%s  ix iy iz=(%2d,%2d,%2d)\n";
  static char    form2p[] = "%s::%s     face index %2d \n";
  static char    form3p[] = "%s::%s     vertex %1d = (%f,%f,%f), normal = (%f,%f,%f)\n";
  static char    prognamp[] = "glDrawCell";

                          /** set Function Values and indexValues of Cube Summits **/
                                  /** Return Triangles number inside current Cube **/
  nbfaces = glDrawSetElemCube(ix, iy, iz);
  if(nbfaces <= 0) return 0;

if(DEBUG) fprintf(stderr, form1p, srcfilenamp, prognamp, ix,iy,iz);

  for(fce=0;  fce < nbfaces;  fce++)
  {
if(DEBUG) fprintf(stderr, form2p, srcfilenamp, prognamp, fce);

    for(is=0;  is < 3;  is++) 
    { posdp = Triangles[fce].Vertpp[is];
      ftotGrad(normaldp, posdp);
      norm = normaldp[0] * normaldp[0] + normaldp[1] * normaldp[1]  
                                       + normaldp[2] * normaldp[2];
      if(norm == 0.0){ invNorm = 1.0;}
      else           { invNorm = 1.0/sqrt(norm);}

      elemCubeNormalp[0] = normaldp[0] * invNorm;
      elemCubeNormalp[1] = normaldp[1] * invNorm;
      elemCubeNormalp[2] = normaldp[2] * invNorm;

                                         /* donne la normale et le vertex à OpenGL */
      posfp[0] = posdp[0];  posfp[1] = posdp[1];  posfp[2] = posdp[2];
      normalfp[0] = -elemCubeNormalp[0];  normalfp[1] = -elemCubeNormalp[1];
                                                   normalfp[2] = -elemCubeNormalp[2];
      glNormal3fv(normalfp); 
      glVertex3fv(posfp);
if(DEBUG) fprintf(stderr, form3p, srcfilenamp, prognamp, is,
                posfp[0], posfp[1], posfp[2], normalfp[0], normalfp[1], normalfp[2]);
    }
  }
  return nbfaces;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       glDrawSetElemCube(int ix, int iy, int iz)
{ double  **dpp, *dp;
  int       elemCubeIndexp[3];
  int       nbtriangles = 0;
  double    isolevel;
  int       boolZ = 1;

                                /** set Function Values at current 8 Cube Summits **/
  dpp = grilleValuespp;  dpp += TailleY * ix + iy;
  dp = *dpp + iz;
  elemCubeValues[0] = *dp;
  elemCubeValues[3] = *(dp +1);
  if(elemCubeValues[0] != 0.0) boolZ = 0;
  if(elemCubeValues[3] != 0.0) boolZ = 0;
  
  dpp = grilleValuespp;  dpp += TailleY * (ix + 1)  + iy;
  dp = *dpp + iz;
  elemCubeValues[1] = *dp;
  elemCubeValues[2] = *(dp +1);
  if(elemCubeValues[1] != 0.0) boolZ = 0;
  if(elemCubeValues[2] != 0.0) boolZ = 0;

  dpp = grilleValuespp;  dpp += TailleY * ix + (iy + 1);
  dp = *dpp + iz;
  elemCubeValues[4] = *dp;
  elemCubeValues[7] = *(dp +1);
  if(elemCubeValues[4] != 0.0) boolZ = 0;
  if(elemCubeValues[7] != 0.0) boolZ = 0;

  dpp = grilleValuespp;  dpp += TailleY * (ix + 1)  + (iy + 1);
  dp = *dpp + iz;
  elemCubeValues[5] = *dp;
  elemCubeValues[6] = *(dp +1);
  if(elemCubeValues[5] != 0.0) boolZ = 0;
  if(elemCubeValues[6] != 0.0) boolZ = 0;


                                              /* calcule les positions des sommets */
  elemCubeIndexp[0] = ix;  elemCubeIndexp[1] = iy;  elemCubeIndexp[2] = iz;  
  grilleIndex2Phys(elemCubePositions[0], elemCubeIndexp);

  elemCubeIndexp[0] = ix + 1;  elemCubeIndexp[1] = iy;  elemCubeIndexp[2] = iz;  
  grilleIndex2Phys(elemCubePositions[1], elemCubeIndexp);
  
  elemCubeIndexp[0] = ix + 1;  elemCubeIndexp[1] = iy;  elemCubeIndexp[2] = iz + 1;  
  grilleIndex2Phys(elemCubePositions[2], elemCubeIndexp);

  elemCubeIndexp[0] = ix;  elemCubeIndexp[1] = iy;  elemCubeIndexp[2] = iz + 1;  
  grilleIndex2Phys(elemCubePositions[3], elemCubeIndexp);

  elemCubeIndexp[0] = ix;  elemCubeIndexp[1] = iy + 1;  elemCubeIndexp[2] = iz;
  grilleIndex2Phys(elemCubePositions[4], elemCubeIndexp);

  elemCubeIndexp[0] = ix + 1;  elemCubeIndexp[1] = iy + 1;  elemCubeIndexp[2] = iz;
  grilleIndex2Phys(elemCubePositions[5], elemCubeIndexp);

  elemCubeIndexp[0] = ix + 1;  elemCubeIndexp[1] = iy + 1;  elemCubeIndexp[2] = iz + 1;
  grilleIndex2Phys(elemCubePositions[6], elemCubeIndexp);

  elemCubeIndexp[0] = ix;  elemCubeIndexp[1] = iy + 1;  elemCubeIndexp[2] = iz + 1;
  grilleIndex2Phys(elemCubePositions[7], elemCubeIndexp);


  isolevel = IsoLeveld;
  nbtriangles = marchCubePolygon(isolevel);
/*
fprintf(stderr, "glDrawSetElemCube ix=%d,iy=%d,iz=%d  nbtriangles=%d\n",
                                                            ix, iy, iz, nbtriangles);
*/
  return nbtriangles;
}
/******************************************************************************/
/******************************************************************************/
