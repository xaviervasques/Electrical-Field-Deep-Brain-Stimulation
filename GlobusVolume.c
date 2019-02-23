/** Xavier Vasques      01/06/2006 **/
/** Last Revised X.V. : 03/06/2006 **/
/** Last Revised X.V. : 01/08/2006 **/
/** Last Revised X.V. : 28/09/2006 **/
/** Last Revised X.V. : 18/10/2006 **/
/** Last Revised X.V. : 14/12/2006 **/

/** Calcul du volume et de la surface de la structure ainsi que le volume
   posterieur et anterieur de la structure **/

#include  <stddef.h>
#include  <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <GL/glut.h>

#include  "dataPoints.glob.h"
#include  "dataFile.h"
#include  "funcBasis.h"
#include  "syslin.h"
#include  "grille.h"
#include  "grille.glob.h"
#include  "funcTot.h"
#include  "main.Reconstruction.h"

/** Equation d'une droite passant par l'origine **/
double    droite(double x, double y, int Gpiside)
{
 double d=0.0;
  
  if(Gpiside == 0)
  {
  d = y - x;
  }
  else
  {
  d = y + x;
  }
return d;
}

/** Calcule des volumes et des surfaces **/
void GlobusVolume()
{
     /* Calcul de la dimension physique de la grille */
     double    typicalSizec;
     /* dimension physique de la grille */
     double    dimPhys; 
     
     typicalSizec = fabs(ptsCentresMin);
     if(ptsCentresMax > typicalSizec) typicalSizec = ptsCentresMax;
     typicalSizec *= 1.5;
     typicalSize = typicalSizec;
     dimPhys = typicalSize * 3.0; 
     
     /* Fichier sortant contenant surfaces et volumes calculés */
     static char Mesure[] = "MesureVolumeNoyau.h";
     FILE *f;

     int    ax,ay,az;
     float  delta=(dimPhys/grilleSize);
     
     /* Parcours de la grille pour le volume*/
     double xyzc[3]={0.0,0.0,0.0};
     double xyzmin[3]={0.0,0.0,0.0};
     double xyzmax[3]={0.0,0.0,0.0};
     
     /* Volume de la structure */
     float  Volume=0.0; 
     /* Volume de la moitié posterieur de la structure */
     float  Volumeposterieur = 0.0;
     /* Volume de la moitié antérieur de la structure */
     float  Volumeanterieur = 0.0;
     /* Surface de la structure */
     float  Surface = 0.0; 
     /* Surface de la moitié posterieur de la structure */
     float  SurfacePost = 0.0;
     /* Surface de la moitié antérieur de la structure */
     float  SurfaceAnt = 0.0;
     
     /* Calcul des volumes */
     printf(" CALCUL DES VOLUMES GPI \n");
     az = 1;
     while(az <= grilleSize)
     {
      ay = 1;
      while(ay <= grilleSize)
      {
       ax = 1;
       while(ax <= grilleSize)
       {
         xyzc[0]= -grilleDimensiondp[0]*0.5 + (delta/2) + ax*(delta);
         xyzc[1]= -grilleDimensiondp[1]*0.5 + (delta/2) + ay*(delta);
         xyzc[2]= -grilleDimensiondp[2]*0.5 + (delta/2) + az*(delta);
         
         xyzmin[0] = -grilleDimensiondp[0]*0.5 + (delta/2) + (ax - 1)*(delta);
         xyzmin[1] = xyzc[1];
         xyzmin[2] = xyzc[2];
                   
         xyzmax[0] = -grilleDimensiondp[0]*0.5 + (delta/2) + (ax + 1)*(delta);
         xyzmax[1] = xyzc[1];
         xyzmax[2] = xyzc[2];
         
         
         if(ftot(xyzc)>=1)  
         {
          Volume = Volume + delta*delta*delta;
          /*printf(" Le volume de la structure vaut : %f \n", Volume);*/
          if(ftot(xyzmin)<1)
          {
           Surface = Surface + delta*delta;
           if(droite(xyzc[0],xyzc[1],Gpiside)<=0)
           {
            SurfacePost = SurfacePost + delta*delta;
           }
           else
           {
            SurfaceAnt = SurfaceAnt + delta*delta;
           }
          }           
          if(ftot(xyzmax)<1)
          {
           Surface = Surface + delta*delta;
           if(droite(xyzc[0],xyzc[1],Gpiside)<=0)
           {
           SurfacePost = SurfacePost + delta*delta;
           }
           else
           {
           SurfaceAnt = SurfaceAnt + delta*delta;
           }
          }
         
         if(droite(xyzc[0],xyzc[1],Gpiside)<=0){Volumeposterieur=Volumeposterieur + delta*delta*delta;}
         if(droite(xyzc[0],xyzc[1],Gpiside)>0){Volumeanterieur=Volumeanterieur + delta*delta*delta;}
        }
        ax = ax + 1; 
       }
       ay = ay + 1;
      }
      az = az + 1;
     }
     
     az = 1;
     while(az <= grilleSize)
     {
      ax = 1;
      while(ax <= grilleSize)
      {
       ay = 1;
       while(ay <= grilleSize)
       {
         xyzc[0]= -grilleDimensiondp[0]*0.5 + (delta/2) + ax*(delta);
         xyzc[1]= -grilleDimensiondp[1]*0.5 + (delta/2) + ay*(delta);
         xyzc[2]= -grilleDimensiondp[2]*0.5 + (delta/2) + az*(delta);
                   
         xyzmin[0] = xyzc[0];
         xyzmin[1] = -grilleDimensiondp[1]*0.5 + (delta/2) + (ay - 1)*(delta);
         xyzmin[2] = xyzc[2];
                   
         xyzmax[0] = xyzc[0];
         xyzmax[1] = -grilleDimensiondp[1]*0.5 + (delta/2) + (ay + 1)*(delta);
         xyzmax[2] = xyzc[2];
                             
         if(ftot(xyzc)>=1)  
         {
          if(ftot(xyzmin)<1)
          {
            Surface = Surface + delta*delta;
            /*printf("La surface de la structure vaut : %f\n",Surface);*/
            if(droite(xyzc[0],xyzc[1],Gpiside)<=0)
            {
            SurfacePost = SurfacePost + delta*delta;
            }
            else
            {
            SurfaceAnt = SurfaceAnt + delta*delta;
            }
          }           
          if(ftot(xyzmax)<1)
          {
            Surface = Surface + delta*delta;
            if(droite(xyzc[0],xyzc[1],Gpiside)<=0)
            {
            SurfacePost = SurfacePost + delta*delta;
            }
            else
            {
            SurfaceAnt = SurfaceAnt + delta*delta;
            }
          }
       }
       ay = ay + 1; 
      }
      ax = ax + 1;
     }
     az = az + 1;
    }

/* Impression du fichier des résultats du calcul*/
     f = fopen(Mesure,"w");
     fPrintF(f," Le volume de la structure est de  : %f mm3\n", Volume);
     fPrintF(f,"\n");
     fPrintF(f," Le volume de la moitié posterieur de la structure est de : %f mm3\n",Volumeposterieur);
     fPrintF(f,"\n");
     fPrintF(f," Le volume de la moitié antérieur de la structure est de : %f mm3\n",Volumeanterieur);
     fPrintF(f,"\n");
     fPrintF(f," La surface de la structure est de : %f mm2\n", Surface);
     fPrintF(f,"\n");
     fPrintF(f," La surface de la moitié postérieur est de : %f mm2\n",SurfacePost);
     fPrintF(f,"\n");
     fPrintF(f," La surface de la moitié antérieur est de : %f mm2\n",SurfaceAnt);
     fPrintF(f,"\n"); 
     fflush(f);
     fclose(f);
     printf(" FIN CALCUL DES VOLUMES GPI \n");
     
return;
}
/******************************************************************************/
/******************************************************************************/
