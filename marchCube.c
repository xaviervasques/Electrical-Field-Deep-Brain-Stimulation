#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>

#include  "marchCube.h"
#include  "marchCube.globalloc.h"
#include  "marchCube.Definitions.General.h"

#include  "grille.h"

                   /* liste des vertex potentiels faisant partie du polygone final */
double    VertexList[12][3];
/******************************************************************************/
/*                                                                            */
/* calcule le polygone pour une case donnée,                                  */
/* et renvoie le nombre de triangles qui le composent                         */
/*                                                                            */
/******************************************************************************/
int       marchCubePolygon(double isolevel)
{ int       Index;
  int       nbtriangles = 0;
  int       i, is;

  Index = 0;
/* détermine l'index dans le EdgeArray qui nous dit quels côtés la surface intersecte */

  if (elemCubeValues[0] < isolevel) 
    Index |= 1;
  if (elemCubeValues[1] < isolevel)
    Index |= 2; 
  if (elemCubeValues[2] < isolevel) 
    Index |= 4; 
  if (elemCubeValues[3] < isolevel)
    Index |= 8; 
  if (elemCubeValues[4] < isolevel) 
    Index |= 16;
  if (elemCubeValues[5] < isolevel) 
    Index |= 32;
  if (elemCubeValues[6] < isolevel)
    Index |= 64;
  if (elemCubeValues[7] < isolevel)
    Index |= 128; 


                /* calcule la position des vertex où la surface intersecte la case */
                                     /* la case est entièrement hors de la surface */
  if ( (Index == 0) || (Index == 255) ) return nbtriangles;


  if (EdgeArray[Index] & 1)
    Interpolate(elemCubePositions[0], elemCubePositions[1],
                      elemCubeValues[0], elemCubeValues[1], isolevel, VertexList[0]);

  if (EdgeArray[Index] & 2)
    Interpolate(elemCubePositions[1], elemCubePositions[2],
                      elemCubeValues[1], elemCubeValues[2], isolevel, VertexList[1]);

  if (EdgeArray[Index] & 4) 
    Interpolate(elemCubePositions[2], elemCubePositions[3],
                      elemCubeValues[2], elemCubeValues[3], isolevel, VertexList[2]);

  if (EdgeArray[Index] & 8)
    Interpolate(elemCubePositions[3], elemCubePositions[0],
                      elemCubeValues[3], elemCubeValues[0], isolevel, VertexList[3]); 

  if (EdgeArray[Index] & 16)
    Interpolate(elemCubePositions[4], elemCubePositions[5],
                      elemCubeValues[4], elemCubeValues[5], isolevel, VertexList[4]) ;

  if (EdgeArray[Index] & 32) 
    Interpolate(elemCubePositions[5], elemCubePositions[6],
                      elemCubeValues[5], elemCubeValues[6], isolevel, VertexList[5]);

  if (EdgeArray[Index] & 64)
    Interpolate(elemCubePositions[6], elemCubePositions[7],
                      elemCubeValues[6], elemCubeValues[7], isolevel, VertexList[6]);

  if (EdgeArray[Index] & 128)
    Interpolate(elemCubePositions[7], elemCubePositions[4],
                      elemCubeValues[7], elemCubeValues[4], isolevel, VertexList[7]); 

  if (EdgeArray[Index] & 256)
    Interpolate(elemCubePositions[0], elemCubePositions[4],
                      elemCubeValues[0], elemCubeValues[4], isolevel, VertexList[8]);

  if (EdgeArray[Index] & 512) 
    Interpolate(elemCubePositions[1], elemCubePositions[5],
                      elemCubeValues[1], elemCubeValues[5], isolevel, VertexList[9]);

  if (EdgeArray[Index] & 1024) 
    Interpolate(elemCubePositions[2], elemCubePositions[6],
                      elemCubeValues[2], elemCubeValues[6], isolevel, VertexList[10]);

  if (EdgeArray[Index] & 2048)
    Interpolate(elemCubePositions[3], elemCubePositions[7],
                      elemCubeValues[3], elemCubeValues[7], isolevel, VertexList[11]);


/* calcule les triangles */
  nbtriangles = 0;
  for(i=0;  FaceArray[Index][i] !=-1;  i+=3) 
  { is = FaceArray[Index][i];
    Triangles[nbtriangles].Vertpp[0] = VertexList[is];
    Triangles[nbtriangles].Nump[0]   = is;
    is = FaceArray[Index][i+1];
    Triangles[nbtriangles].Vertpp[1] = VertexList[is]; 
    Triangles[nbtriangles].Nump[1]   = is;
    is = FaceArray[Index][i+2];
    Triangles[nbtriangles].Vertpp[2] = VertexList[is]; 
    Triangles[nbtriangles].Nump[2]   = is; 
    nbtriangles++; 
  }
  return nbtriangles;
}
/******************************************************************************/
/*                                                                            */
/*  vect1 & vect2 sont les extremités                                         */
/*  (avec val1 & val2 les valeurs correspondantes)                            */
/*  vect est la position recherchée                                           */
/******************************************************************************/
int       Interpolate(double vect1[3], double vect2[3], double val1, double val2,
                                                     double isolevel, double vect[3])
{ double     coef, valdif;
  double    dtest;

  dtest = isolevel - val1;
  if ( fabs(dtest) < 0.00001 )
  { vect[0] = vect1[0];
    vect[1] = vect1[1];
    vect[2] = vect1[2];
    return 0; 
  }
  dtest = isolevel - val2;
  if ( fabs(dtest) < 0.00001 )
  { vect[0] = vect2[0];
    vect[1] = vect2[1];
    vect[2] = vect2[2]; 
    return 0; 
  }

  valdif = val2 - val1;
  if( fabs(valdif) < 0.000001 ) coef = 0.5;
  else
  { coef = (isolevel - val1) / (val2 - val1);
  }
  vect[0] = vect1[0] + coef * (vect2[0] - vect1[0]);
  vect[1] = vect1[1] + coef * (vect2[1] - vect1[1]);
  vect[2] = vect1[2] + coef * (vect2[2] - vect1[2]); 
  return 0;
}
/******************************************************************************/
/******************************************************************************/
