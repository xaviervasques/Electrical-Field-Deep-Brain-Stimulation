/*******************************************************************************
URMA
Gerard Mennessier
Xavier Vasques
Last Revised : 29/03/2005

Problème corrigé au niveau de la taille de la grille dans grilleSetSizes
(Inversement de deux lignes sinon pas de prise en compte de la taille)
********************************************************************************/

#include  <math.h>
#include  <stdio.h>
#include  <stdlib.h>

#include  "grille.globalloc.h"
#include  "dataPoints.glob.h"

#include  "grille.h"
#include  "funcTot.h"

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      grilleSetDim(double dimen)
{
  grilleDimensiondp[0] = dimen;
  grilleDimensiondp[1] = dimen;
  grilleDimensiondp[2] = dimen;

  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      grilleSetSizes(int sizeX, int sizeY, int sizeZ)
{ size_t    totz;
  int       ix, iy;
  double  **dpp, *dp;


  
  /*grilleSizep[0] = sizeX;  grilleSizep[1] = sizeY;  grilleSizep[2] = sizeZ;*/
  
  sizeX=grilleSizep[0];  sizeY=grilleSizep[1] ;  sizeZ=grilleSizep[2] ;
  
  TailleX = grilleSizep[0];  TailleY = grilleSizep[1];  TailleZ = grilleSizep[2];

  grilleSizeM1dp[0] = sizeX - 1.0;
  grilleSizeM1dp[1] = sizeY - 1.0;
  grilleSizeM1dp[2] = sizeZ - 1.0;

  totz = sizeX * sizeY * sizeZ;
  dblPVecRealloc(&grilleValuesVec, totz, 10);
  grilleValuesp = grilleValuesVec.p;
  dp = grilleValuesp;

  grilleValuespp = (double **)ptrAlloc(TailleX * TailleY, "grille::grilleSetSizes");
  dpp = grilleValuespp;
  for(ix = 0;  ix < TailleX;  ix++)
  { for(iy = 0;  iy < TailleY;  iy++,  dp += TailleZ)
    { *dpp++ = dp;
    }
  }

  return; 
}
/******************************************************************************/
/*                                                                            */
/* from  grille point indices              posIndexp[3],                      */
/* compute      its   physical coordinates posPhysp[3]                        */
/*                                                                            */
/******************************************************************************/
void      grilleIndex2Phys(double *posPhysp, int *posIndexp)
{

  *posPhysp++ = ( (*posIndexp++)/grilleSizeM1dp[0] - 0.5) * grilleDimensiondp[0];
  *posPhysp++ = ( (*posIndexp++)/grilleSizeM1dp[1] - 0.5) * grilleDimensiondp[1];
  *posPhysp   = ( (*posIndexp  )/grilleSizeM1dp[2] - 0.5) * grilleDimensiondp[2];

  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      grilleValues()
{ int       ix, iy, iz;                                               /* compteurs */
  int       indexp[3];
  double    posp[3];               /* la position physique d'un point de la grille */
  double    fval = 0.0;
  double  **dpp, *dp;
    
fprintf(stderr, "grille::grilleValues  BEGIN \n");
 
  dpp = grilleValuespp;
  for(ix=0;  ix < TailleX;  ix++)
  { indexp[0] = ix;
    for(iy=0;  iy < TailleY;  iy++, dpp++)
    { indexp[1] = iy;
      dp = *dpp;
      for(iz=0;  iz < TailleZ;  iz++, dp++)
      { indexp[2] = iz;

        grilleIndex2Phys(posp, indexp);
        fval = ftot(posp);
        *dp = fval;
      }

    }
fprintf(stderr, "grille::grilleValues  ix=%4d DONE for TailleY, TailleZ \n", ix);
}
return;
}

                             


/******************************************************************************/
/******************************************************************************/
