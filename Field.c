/* Programme de visualisations des lignes ISO-CHAMPS et l'intersection avec   */
/* la surface donné en 3D : Calcul de volumes  et surfaces                    */
/* Xavier Vasques                                                             */
/* Last Revised : 10/04/2006                                                  */
/* L.R. : 01/08/2006                                                          */
/* L.R. : 22/08/2006                                                          */
/* L.R. : 28/09/2006                                                          */
/* L.R. : 18/10/2006                                                          */
/* L.R. : 14/12/2006                                                          */
/* L.R. : 06/05/2007                                                          */

/* Biblio standard*/
#include  <stdio.h>
#include  <string.h>
#include  <stdlib.h>
#include  <math.h>
#include  <GL/glut.h>


#include  "Field.h"
#include  "dataPoints.glob.h"
#include  "dataFile.h"
#include  "funcBasis.h"
#include  "syslin.h"
#include  "grille.h"
#include  "grille.glob.h"
#include  "funcTot.h"
#include  "main.Reconstruction.h"

/* Visualisation des lignes de champs autour de l'électrode*/
void IsoField(double theta, double phi)
{
     /* Declaration de variables */
     int NbLignes = 0;
     int i;
     float a,b,c;
     
     /* Ouverture du fichier */
     FILE* FieldPoint;
     FieldPoint = fopen("FieldPoint.h","r");

     /* Tester le fichier */     
     if(FieldPoint != NULL){}
     else{printf("le fichier n'a pu être ouvert\n");}
 
     /* Compteur du nombre de lignes du fichier */
     while(!feof(FieldPoint))
     {/*wchar_t *j=(wchar_t *) malloc(sizeof(wchar_t));*/
     a,b,c=0.0;
     fscanf(FieldPoint,"%f,%f,%f,\n",&a,&b,&c);
     NbLignes=NbLignes+1;
     }
     
     float champs[NbLignes][3];  
     
     /* On revient au début du fichiers */
     rewind(FieldPoint);
     
     /* On rempli le tableau champs[][] */
     i=0;
     while(!feof(FieldPoint) && i<NbLignes)
     {
     fscanf(FieldPoint,"%f,%f,%f,\n", &(champs[i][0]),&(champs[i][1]),&(champs[i][2]));
     i=i+1;
     }
     
     /* Plus besoin du fichier, on le ferme */
     fclose(FieldPoint);   
      
     /* Pour chaque z=champs[][1], on trace cercle de rayon=champs[][2] */
     static GLUquadricObj *isoline;   
     isoline = gluNewQuadric();
     gluQuadricDrawStyle(isoline, GLU_LINE);
     for (i=0;i<NbLignes;i++)
     {
     glPushMatrix();
     glRotatef ( (GLfloat)phi,   0.0, 1.0, 0.0);
     glRotatef ( (GLfloat)theta, 1.0, 0.0, 0.0);
     glTranslatef(0.0, 0.0, champs[i][1]);
     glColor3f(1.0, 1.0, 1.0);
     gluDisk(isoline, champs[i][2],champs[i][2] , 100, 100);
     glPopMatrix();
     }
}

/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/

/* Equation d'un cercle*/
double    fcercle(double currentx,double currenty,double rayon)
{
double f = 0.0;
f = (currentx)*(currentx) + (currenty)*(currenty)- (rayon)*(rayon);
return f;
}

/* Calcul de la plus petite valeur de la premiere colonne d'un tableau*/
double minvalue(float tab[][3])
{
double min = tab[0][0];
size_t Np = (sizeof tab) / (sizeof tab[2]);
int i;
for(i=1; i<Np; i++)
{
         if(tab[i][0] < min)
         { 
         min = tab[i][0];
         }
}
return min;
}

/* Calcul du nombre de ligne min dans le tableau*/
int numbermin(float tab[][3], double min, int Np)
{
int i, a;
i = 0;
a = 0;
while(i<Np)
{
if(tab[i][0] == min){a = a + 1;}
i=i+1;
}
return a;
}

/* Equation d'une droite passant par l'origine*/
double    droite1(double x, double y, int Gpiside)
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

/******************************************************************************/
/******************************************************************************/
/*  Calcul pour chaque lignes de champs le volume intersecté avec le pallidum */
/* Ainsi que le calcul de l'intersection avec la motié postérieur et antérieur*/
/******************************************************************************/
/******************************************************************************/
void InterVolumeFile(double theta, double phi)
{
/* Déclaration de la variable InterVolumeFile pour l'impression du fichier final*/
     static char InterVolumeFile[] = "MesureIntersectionChamp.h";
     FILE *o;
     
      /* Declaration de variables */
     int Np = 0;
     int i;
     float a,b,c;
     
     /* Ouverture du fichier */
     FILE* FieldPoint;
     FieldPoint = fopen("FieldPoint.h","r");

     /* Tester le fichier */     
     if(FieldPoint != NULL){}
     else{printf("le fichier n'a pu être ouvert\n");}
 
     /* Compteur du nombre de lignes du fichier */
     while(!feof(FieldPoint))
     {/*wchar_t *j=(wchar_t *) malloc(sizeof(wchar_t));*/
     a,b,c=0.0;
     fscanf(FieldPoint,"%f,%f,%f,\n",&a,&b,&c);
     Np=Np+1;
     }
     
     float Champs[Np][3];  
     
     /* On revient au début du fichiers */
     rewind(FieldPoint);
     
     /* On rempli le tableau champs[][] */
     i=0;
     while(!feof(FieldPoint) && i<Np)
     {
     fscanf(FieldPoint,"%f,%f,%f,\n", &(Champs[i][0]),&(Champs[i][1]),&(Champs[i][2]));
     i=i+1;
     }
     
     /* Plus besoin du fichier, on le ferme */
     fclose(FieldPoint);   

     printf("Le fichier Champ[][] est rempli\n");
     for(i=0; i<Np;i++)
     {
     printf("%f,%f,%f\n",Champs[i][0],Champs[i][1],Champs[i][2]);
     }
     
/*Declaration de variable pour la dimension de la grille*/
     double    typicalSizec;
     double    dimPhys;
     typicalSizec = fabs(ptsCentresMin);
     if(ptsCentresMax > typicalSizec) typicalSizec = ptsCentresMax;
     typicalSizec *= 1.5;
     typicalSize = typicalSizec;
     dimPhys = typicalSize * 3.0;
     float delta=((dimPhys) / (grilleSize));
     
/*Déclaration de constantes*/
     int j, k ,l , finalsize;
     double  minfield;
     int  M;
     double x,y,z= 0.0;
     double xyz[4];
     
/* Conversion des degrés en radians*/
     double thetaRad = theta * 0.017453293;
     double phiRad = phi * 0.017453293; 
          
/* On va trier par ordre croissant les lignes de champs */
     int mini;
     double temp[3];
     for(i=0; i<Np-1; i++)
     {
     mini = i;
              for(j=i+1; j<Np; j++)
              {
              if(Champs[j][0] < Champs[mini][0]){mini = j;}
              } 
                  temp[0] = Champs[i][0];
                  temp[1] = Champs[i][1];
                  temp[2] = Champs[i][2];
                  
                  Champs[i][0] = Champs[mini][0];
                  Champs[i][1] = Champs[mini][1];
                  Champs[i][2] = Champs[mini][2];
                  
                  Champs[mini][0] = temp[0];
                  Champs[mini][1] = temp[1];
                  Champs[mini][2] = temp[2];
     }

/* On commence l'impression du fichier*/
   o = fopen(InterVolumeFile,"w");
   fPrintF(o, "Ligne de Champs; Volume Intersection; Volume Inter 1/2 postérieur; Volume Inter 1/2 antérieur; \n");
   fPrintF(o, " \n");

   int p= 0;
   while(p < Np)
   {
     
     M = numbermin(Champs, Champs[p][0],Np);
     /*printf("le nombre de ligne de champs %f  est %i\n",Champs[p][0],M);*/

/* Newfield est un tableau ne contenant que la ligne de champs choisi ainsi que son z et son rayon*/
     double newfield[M][3];
     int e=0;
     for(i=0; i<Np; i++)
     {
             if(Champs[i][0] == Champs[p][0])
             {
              newfield[e][0] = Champs[i][0];
              newfield[e][1] = Champs[i][1];
              newfield[e][2] = Champs[i][2];
              /*printf("Position %i du tableau : (E,Z,R) = (%f,%f,%f)\n",a,newfield[a][0],newfield[a][1],newfield[a][2]);*/
              e = e + 1;
             }
     }
     
/* On trie le tableau newfield de sorte que les z  soit ordonné par ordre croissant*/
     size_t Np1 = (sizeof newfield) / (sizeof newfield[2]);
     double temp[3];
     int min;
     /*printf("NP1 vaut %i\n",Np1);*/
     
     for(i=0; i<Np1-1; i++)
     {
     min = i;
              for(j=i+1; j<Np1; j++)
              {
              if(newfield[j][1] < newfield[min][1]){min = j;}
              } 
                  temp[0] = newfield[i][0];
                  temp[1] = newfield[i][1];
                  temp[2] = newfield[i][2];
                  
                  newfield[i][0] = newfield[min][0];
                  newfield[i][1] = newfield[min][1];
                  newfield[i][2] = newfield[min][2];
                  
                  newfield[min][0] = temp[0];
                  newfield[min][1] = temp[1];
                  newfield[min][2] = temp[2];
     }
/* Nous devons ensuite supprimer les doublons qui ne servent à rien du tout*/
   double finalfield[Np1][3];
   finalsize=0;
   l=0;
   for(i=0;i<Np1;i++)
   {
   l = 0;
   for(j=i+1; j<Np1;j++)
   {
      if(newfield[i][1] == newfield[j][1] && newfield[i][2] == newfield[j][2])
      {
      l = l +1;
      }
   }
   if(l==0)
   {  
   finalfield[finalsize][0] = newfield[i][0];
   finalfield[finalsize][1] = newfield[i][1];
   finalfield[finalsize][2] = newfield[i][2];
   /*printf(" %f, %f, %f,\n",finalfield[finalsize][0],finalfield[finalsize][1],finalfield[finalsize][2]);*/
   finalsize=finalsize+1;
   }
   }

   /*printf("finalsize vaut %i\n",finalsize);*/

/* Calcul du Volume d'intersection en ne prenant pas en compte l'electrode et les creux*/
   float zi;
   float ri;
   float jymax;
   float jy;
   float jxmax;
   float jx;
   float InterVolume = 0.0;
   float InterVolPost = 0.0; /* Intersection avec la motié postérieur*/
   float InterVolAnt = 0.0; /* Intersection avec la moitié antérieur*/
   float Cube = 0.0;

   for(i = 0; i < finalsize;i++)
   {
      zi = finalfield[i][1];
      ri = finalfield[i][2];
      
      if(finalfield[i][1] == finalfield[i+1][1] && finalfield[i][2] != finalfield[i+1][2])
      {
      }
      else
      {
       j=0;
       while(j*(delta/2) <= 2*ri)
       {
        y = (ri - j*(delta/2));
        k=0;
        while(k*(delta/2) <= 2*sqrt(ri*ri - y*y))
        {
         x = sqrt(ri*ri - y*y) - k*(delta/2);
      
         xyz[0] = x * cos(phiRad) + sin(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad))+ (targetp[0] - cg[0]) ;
         xyz[1] = y*cos(thetaRad) - zi*sin(thetaRad)+ (targetp[1] - cg[1]);
         xyz[2] = (-x*sin(phiRad) + cos(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad)))+ (targetp[2] - cg[2]);
      
         if(finalfield[i][1] == finalfield[i-1][1] && finalfield[i][2] != finalfield[i-1][2])
         {
          if(zi <= 4.25)
          {
           if( ( sqrt(y*y + x*x) < finalfield[i-1][2]) || (sqrt(y*y + x*x) < 0.635))
           {
           }
           else
           {
            if(ftot(xyz)>=1)
            {
             Cube =(delta/2)*(delta/2)* (finalfield[i+1][1] - finalfield[i][1]);
             InterVolume = InterVolume + Cube;
             /*printf("Le volume ajouté de la ligne %f est de %f\n" ,Champs[p][0],InterVolume);*/
             if(droite1(xyz[0],xyz[1],Gpiside)<=0){InterVolPost = InterVolPost + Cube;}
             if(droite1(xyz[0],xyz[1],Gpiside)>0){InterVolAnt = InterVolAnt + Cube;}
            }
           }
          }
          else
          {
           if( ( sqrt(y*y + x*x) < finalfield[i-1][2]))
           {
           }
           else
           {
            if(ftot(xyz) >=1)
            {
            Cube =(delta/2)*(delta/2)* (finalfield[i+1][1] - finalfield[i][1]);
            InterVolume = InterVolume + Cube;
            /*printf("Le volume ajouté de la ligne %f est de %f\n" ,Champs[p][0],InterVolume);*/
            if(droite1(xyz[0],xyz[1],Gpiside)<=0){InterVolPost = InterVolPost + Cube;}
            if(droite1(xyz[0],xyz[1],Gpiside)>0){InterVolAnt = InterVolAnt + Cube;}
            }
           }
          }
         }
         else
         {
          if(zi<= 4.25)
          {
           if( (sqrt(y*y + x*x) < 0.635))
           {
           }
           else
           {
            if(ftot(xyz)>=1)
            {
            Cube =(delta/2)*(delta/2)* (finalfield[i+1][1] - finalfield[i][1]);
            InterVolume = InterVolume + Cube;
            /*printf("Le volume ajouté de la ligne %f est de %f\n" ,Champs[p][0],InterVolume);*/
            if(droite1(xyz[0],xyz[1],Gpiside)<=0){InterVolPost = InterVolPost + Cube;}
             if(droite1(xyz[0],xyz[1],Gpiside)>0){InterVolAnt = InterVolAnt + Cube;}
            }
           }
          }
          else
          {
           if(ftot(xyz)>=1)
           {
           Cube =(delta/2)*(delta/2)* (finalfield[i+1][1] - finalfield[i][1]);
           InterVolume = InterVolume + Cube;
           /*printf("Le volume ajouté de la ligne %f est de %f\n" ,Champs[p][0],InterVolume);*/
           if(droite1(xyz[0],xyz[1],Gpiside)<=0){InterVolPost = InterVolPost + Cube;}
           if(droite1(xyz[0],xyz[1],Gpiside)>0){InterVolAnt = InterVolAnt + Cube;}
           }
          }
         }
    k = k+ 1;
    }
    j= j + 1;
    }
   }
  }
    fPrintF(o, " %f, %f, %f, %f\n",Champs[p][0], InterVolume,InterVolPost,InterVolAnt);
    p = p + M;
    }
    fflush(o);
    fclose(o);
    return;
    }     
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/* Calcul du volume de chaque lignes de champs. On ne prend pas en compte     */
/* l'électrode et le creux aux extrmités du champ                             */
/******************************************************************************/
/******************************************************************************/

void VolumeChamps(double theta, double phi)
{
/* Déclaration de la variable VolumeFieldValue pour l'impression du fichier final*/     
     static char VolumeFieldValue[] = "MesureVolumeChamp.h";
     FILE *o;
     
 /* Declaration de variables */
     int Np = 0;
     int i;
     float a,b,c;
     
     /* Ouverture du fichier */
     FILE* FieldPoint;
     FieldPoint = fopen("FieldPoint.h","r");

     /* Tester le fichier */     
     if(FieldPoint != NULL){}
     else{printf("le fichier n'a pu être ouvert\n");}
 
     /* Compteur du nombre de lignes du fichier */
     while(!feof(FieldPoint))
     {/*wchar_t *j=(wchar_t *) malloc(sizeof(wchar_t));*/
     a,b,c=0.0;
     fscanf(FieldPoint,"%f,%f,%f,\n",&a,&b,&c);
     Np=Np+1;
     }
     
     float Champs[Np][3];  
     
     /* On revient au début du fichiers */
     rewind(FieldPoint);
     
     /* On rempli le tableau champs[][] */
     i=0;
     while(!feof(FieldPoint) && i<Np)
     {
     fscanf(FieldPoint,"%f,%f,%f,\n", &(Champs[i][0]),&(Champs[i][1]),&(Champs[i][2]));
     i=i+1;
     }
     
     /* Plus besoin du fichier, on le ferme */
     fclose(FieldPoint);   

     printf("Le fichier Champ[][] est rempli\n");
     for(i=0; i<Np;i++)
     {
     printf("%f,%f,%f\n",Champs[i][0],Champs[i][1],Champs[i][2]);
     }
     
/*Declaration de variable pour la dimension de la grille*/
     double    typicalSizec;
     double    dimPhys;
     typicalSizec = fabs(ptsCentresMin);
     if(ptsCentresMax > typicalSizec) typicalSizec = ptsCentresMax;
     typicalSizec *= 1.5;
     typicalSize = typicalSizec;
     dimPhys = typicalSize * 3.0;
     float delta=((dimPhys) / (grilleSize));
     
/*Déclaration de constantes*/
     int j, k ,l , finalsize;
     double  minfield;
     int  M;
     double x,y,z= 0.0;
     double xyz[4];
     double xyzco[4];

/* Conversion des degrés en radians */
     double thetaRad = theta * 0.017453293;
     double phiRad = phi * 0.017453293; 
     
/* On va trier par ordre croissant les lignes de champs */
     int mini;
     double temp[3];
     for(i=0; i<Np-1; i++)
     {
     mini = i;
              for(j=i+1; j<Np; j++)
              {
              if(Champs[j][0] < Champs[mini][0]){mini = j;}
              } 
                  temp[0] = Champs[i][0];
                  temp[1] = Champs[i][1];
                  temp[2] = Champs[i][2];
                  
                  Champs[i][0] = Champs[mini][0];
                  Champs[i][1] = Champs[mini][1];
                  Champs[i][2] = Champs[mini][2];
                  
                  Champs[mini][0] = temp[0];
                  Champs[mini][1] = temp[1];
                  Champs[mini][2] = temp[2];
     }
     
o = fopen(VolumeFieldValue,"w");
fPrintF(o, "Ligne de Champs; Volume complet;\n");
fPrintF(o, " \n");

int p= 0;
while(p < Np)
{
     
     M = numbermin(Champs, Champs[p][0],Np);
     /*printf("le nombre de ligne de champs %f  est %i\n",Champs[p][0],M);*/
     
/* newfield est un tableau ne contenant que la ligne de champs choisi ainsi   */
/* que son z et son rayon*/
     double newfield[M][3];
     int e=0;
     for(i=0; i<Np; i++)
     {
             if(Champs[i][0] == Champs[p][0])
             {
              
              newfield[e][0] = Champs[i][0];
              newfield[e][1] = Champs[i][1];
              newfield[e][2] = Champs[i][2];
              /*printf("Position %i du tableau : (E,Z,R) = (%f,%f,%f)\n",a,newfield[a][0],newfield[a][1],newfield[a][2]);*/
              e = e + 1;
             
             }
     }
     
/* On trie le tableau newfield de sorte que les z  soit ordonné par ordre croissant*/
     size_t Np1 = (sizeof newfield) / (sizeof newfield[2]);
     double temp[3];
     int min;
     /*printf("NP1 vaut %i\n",Np1);*/
     
     for(i=0; i<Np1-1; i++)
     {
     min = i;
              for(j=i+1; j<Np1; j++)
              {
              if(newfield[j][1] < newfield[min][1]){min = j;}
              } 
                  temp[0] = newfield[i][0];
                  temp[1] = newfield[i][1];
                  temp[2] = newfield[i][2];
                  
                  newfield[i][0] = newfield[min][0];
                  newfield[i][1] = newfield[min][1];
                  newfield[i][2] = newfield[min][2];
                  
                  newfield[min][0] = temp[0];
                  newfield[min][1] = temp[1];
                  newfield[min][2] = temp[2];
     }
/* Nous devons ensuite supprimer les doublons qui ne servent à rien du tout*/
   double finalfield[Np1][3];
   finalsize=0;
   l=0;

   for(i=0;i<Np1;i++)
   {
    l = 0;
    for(j=i+1; j<Np1;j++)
    {
      if(newfield[i][1] == newfield[j][1] && newfield[i][2] == newfield[j][2])
      {
      l = l +1;
      }
    }
    if(l==0)
    {  
    finalfield[finalsize][0] = newfield[i][0];
    finalfield[finalsize][1] = newfield[i][1];
    finalfield[finalsize][2] = newfield[i][2];
    /*printf(" %f, %f, %f,\n",finalfield[finalsize][0],finalfield[finalsize][1],finalfield[finalsize][2]);*/
    finalsize=finalsize+1;
    }
   }
/*printf("finalsize vaut %i\n",finalsize);*/

/* Calcul du volume du champs actuel en parcourant toute la grille */
   float zi;
   float ri;
   float jymax;
   float jy;
   float jxmax;
   float jx;
   float InterVolume = 0.0;
   float Cube = 0.0;

   for(i = 0; i < finalsize;i++)
   {
      zi = finalfield[i][1];
      ri = finalfield[i][2];
      
      if(finalfield[i][1] == finalfield[i+1][1] && finalfield[i][2] != finalfield[i+1][2])
      {
      }
      else
      {
      j=0;
      while(j*(delta/2) <= 2*ri)
      {
        y = (ri - j*(delta/2));
        k=0;
        while(k*(delta/2) <= 2*sqrt(ri*ri - y*y))
        {
         x = sqrt(ri*ri - y*y) - k*(delta/2);
      
         xyz[0] = x * cos(phiRad) + sin(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad))+ (targetp[0] - cg[0]) ;
         xyz[1] = y*cos(thetaRad) - zi*sin(thetaRad)+ (targetp[1] - cg[1]);
         xyz[2] = -x*sin(phiRad) + cos(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad))+ (targetp[2] - cg[2]);
      
      if(finalfield[i][1] == finalfield[i-1][1] && finalfield[i][2] != finalfield[i-1][2])
      {
       if(zi <=4.25)
       {
          if(( sqrt(y*y + x*x) < finalfield[i-1][2]) || (sqrt(y*y + x*x) < 0.635))
          {
          }
          else
          {
           Cube =(delta/2)*(delta/2)* (finalfield[i+1][1] - finalfield[i][1]);
           InterVolume = InterVolume + Cube;
           /*printf("le volume ajouté de la ligne de champ %f est %f \n",Champs[p][0],InterVolume);*/
          }
       }
       else
       {
        if( ( sqrt(y*y + x*x) < finalfield[i-1][2]))
        {
        }
        else
        {
        Cube =(delta/2)*(delta/2)* (finalfield[i+1][1] - finalfield[i][1]);
        InterVolume = InterVolume + Cube;
        /*printf("le volume ajouté de la ligne de champ %f est %f \n",Champs[p][0],InterVolume);*/
        }
       }
     }
     else
     {
       if(zi<=4.25)
       {
         if( ( sqrt(y*y + x*x) < 0.635))
         {
         }
         else
         {
           Cube =(delta/2)*(delta/2)* (finalfield[i+1][1] - finalfield[i][1]);
           InterVolume = InterVolume + Cube;
           /*printf("le volume ajouté de la ligne de champ %f est %f \n",Champs[p][0],InterVolume);*/
         }
       }
       else
       {
         Cube =(delta/2)*(delta/2)* (finalfield[i+1][1] - finalfield[i][1]);
         InterVolume = InterVolume + Cube;
         /*printf("le volume ajouté de la ligne de champ %f est %f \n",Champs[p][0],InterVolume);*/
       }
    }
    k = k+ 1;
    }
    j= j + 1;
    }
}
}
fPrintF(o, " %f, %f\n",Champs[p][0], InterVolume);
p = p + M;
}
fflush(o);
fclose(o);
return;
}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/*      Calcul des surfaces intersectées pour chaque lignes de champs         */
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
   
void InterSurfaceValue(double theta, double phi)
{
  
/* Declaration de la variable IntersurfaceValue pour l'imression du fichier   */
/* contenant les valeurs de surface intersectées pour chaque lignes de champs */                        
     static char InterSurfaceValue[] ="MesureSurfaceInOut.h";
     FILE *g;     
     
 /* Declaration de variables */
     int Np = 0;
     int i;
     float a,b,c;
     
     /* Ouverture du fichier */
     FILE* FieldPoint;
     FieldPoint = fopen("FieldPoint.h","r");

     /* Tester le fichier */     
     if(FieldPoint != NULL){}
     else{printf("le fichier n'a pu être ouvert\n");}
 
     /* Compteur du nombre de lignes du fichier */
     while(!feof(FieldPoint))
     {/*wchar_t *j=(wchar_t *) malloc(sizeof(wchar_t));*/
     a,b,c=0.0;
     fscanf(FieldPoint,"%f,%f,%f,\n",&a,&b,&c);
     Np=Np+1;
     }
     
     float Champs[Np][3];  
     
     /* On revient au début du fichiers */
     rewind(FieldPoint);
     
     /* On rempli le tableau champs[][] */
     i=0;
     while(!feof(FieldPoint) && i<Np)
     {
     fscanf(FieldPoint,"%f,%f,%f,\n", &(Champs[i][0]),&(Champs[i][1]),&(Champs[i][2]));
     i=i+1;
     }
     
     /* Plus besoin du fichier, on le ferme */
     fclose(FieldPoint);   

     printf("Le fichier Champ[][] est rempli\n");
     for(i=0; i<Np;i++)
     {
     printf("%f,%f,%f\n",Champs[i][0],Champs[i][1],Champs[i][2]);
     }
     
/*Declaration de variable pour le calcul de la dimension de la grille*/
     double    typicalSizec;
     double    dimPhys;
     typicalSizec = fabs(ptsCentresMin);
     if(ptsCentresMax > typicalSizec) typicalSizec = ptsCentresMax;
     typicalSizec *= 1.5;
     typicalSize = typicalSizec;
     dimPhys = typicalSize * 3.0;
     float delta=((dimPhys) / (grilleSize));
     
/*Déclaration de constantes*/
     int j, k ,l ,finalsize;
     double  minfield, zi, ri;
     int  M;
     double x,y,z= 0.0;
     double xyz[4];
     double xyzco[4];

/* Conversion des angles en Radian */
     double thetaRad = theta * 0.017453293;
     double phiRad = phi * 0.017453293;  
     
/* Trie par ordre croissant des lignes de champs */
     int mini;
     double temp[3];
     for(i=0; i<Np-1; i++)
     {
     mini = i;
              for(j=i+1; j<Np; j++)
              {
              if(Champs[j][0] < Champs[mini][0]){mini = j;}
              } 
                  temp[0] = Champs[i][0];
                  temp[1] = Champs[i][1];
                  temp[2] = Champs[i][2];
                  
                  Champs[i][0] = Champs[mini][0];
                  Champs[i][1] = Champs[mini][1];
                  Champs[i][2] = Champs[mini][2];
                  
                  Champs[mini][0] = temp[0];
                  Champs[mini][1] = temp[1];
                  Champs[mini][2] = temp[2];
     }
     
/* On commence l'impression du fichier contenant les valeurs des surfaces intersectées*/
     g = fopen(InterSurfaceValue,"w");
     fPrintF(g, "Ligne de Champs; Surface Intersection; Surface In; Surface Out\n");
     fPrintF(g, " \n");

     int p= 0;
     while(p < Np)
     {
     
/* Nombre de lignes de champs inclue dans le calcul de la surface*/
     M = numbermin(Champs, Champs[p][0],Np);
     printf("le nombre de ligne de champs %f  est %i\n",Champs[p][0],M);
     
/* Newfield est un tableau ne contenant que la ligne de champs choisi ainsi que son z et son rayon*/
     double newfield[M][3];
     int e=0;
     for(i=0; i<Np; i++)
     {
             if(Champs[i][0] == Champs[p][0])
             {
              
              newfield[e][0] = Champs[i][0];
              newfield[e][1] = Champs[i][1];
              newfield[e][2] = Champs[i][2];
              e = e + 1;
             
             }
     }
     
/* Trie du tableau newfield de sorte que les z  soit ordonné par ordre croissant*/
     size_t Np1 = (sizeof newfield) / (sizeof newfield[2]);
     double temp[3];
     int min;
     printf("NP1 vaut %i\n",Np1);
     
     for(i=0; i<Np1-1; i++)
     {
     min = i;
              for(j=i+1; j<Np1; j++)
              {
              if(newfield[j][1] < newfield[min][1]){min = j;}
              } 
                  temp[0] = newfield[i][0];
                  temp[1] = newfield[i][1];
                  temp[2] = newfield[i][2];
                  
                  newfield[i][0] = newfield[min][0];
                  newfield[i][1] = newfield[min][1];
                  newfield[i][2] = newfield[min][2];
                  
                  newfield[min][0] = temp[0];
                  newfield[min][1] = temp[1];
                  newfield[min][2] = temp[2];
     }
/* Nous devons ensuite supprimer les doublons qui ne servent à rien du tout*/
   double finalfield[Np1][3];
   finalsize=0;
   l=0;

   for(i=0;i<Np1;i++)
   {
    l = 0;
    for(j=i+1; j<Np1;j++)
    {
      if(newfield[i][1] == newfield[j][1] && newfield[i][2] == newfield[j][2])
      {
      l = l +1;
      }
    }
    if(l==0)
    {  
    finalfield[finalsize][0] = newfield[i][0];
    finalfield[finalsize][1] = newfield[i][1];
    finalfield[finalsize][2] = newfield[i][2];
    printf(" %f, %f, %f,\n",finalfield[finalsize][0],finalfield[finalsize][1],finalfield[finalsize][2]);
    finalsize=finalsize+1;
    }
   }
   printf("finalsize vaut %i\n",finalsize);
         
/* On commence le calcul de la surface pour la ligne de champs en cours (boucle while) */
   float Surface = 0.0;
   float SurfaceIn = 0.0;
   float SurfaceOut = 0.0;
      
/******************************************************************************/
/* Intersection pour x supérieur à 0*/
/******************************************************************************/
   printf("Intersection du calcul surfacique pour x supérieur à 0\n");
   double xyzmax[3];
   double xyzmin[3];

   for(i = 0; i < finalsize - 1;i++)
   {
      zi = finalfield[i][1] ;
      ri = finalfield[i][2];
   j=0;
   while(j*(delta/2)<= 2* ri)
   {
   y = ri - j*(delta/2) ;
   k=0;
      while(k*(delta/2) <= 2*sqrt(ri*ri - y*y))
      {
      x = sqrt(ri*ri - y*y) - k*(delta/2)  ;
      
      xyz[0] = x * cos(phiRad) + sin(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad))+ (targetp[0] - cg[0]);
      xyz[1] = y*cos(thetaRad) - zi*sin(thetaRad)+ (targetp[1] - cg[1]);
      xyz[2] = (-x*sin(phiRad) + cos(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad)))+ (targetp[2] - cg[2]);
      
      xyzmax[0] = (sqrt(ri*ri - y*y) - (k+1)*(delta/2)) * cos(phiRad) + sin(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad))+ (targetp[0] - cg[0]) ;
      xyzmax[1] = y*cos(thetaRad) - zi*sin(thetaRad)+ (targetp[1] - cg[1]);
      xyzmax[2] = (-(sqrt(ri*ri - y*y) - (k+1)*(delta/2))*sin(phiRad) + cos(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad)))+ (targetp[2] - cg[2]);
      
      xyzmin[0] = (sqrt(ri*ri - y*y) - (k-1)*(delta/2)) * cos(phiRad) + sin(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad))+ (targetp[0] - cg[0]) ;
      xyzmin[1] = y*cos(thetaRad) - zi*sin(thetaRad)+ (targetp[1] - cg[1]);
      xyzmin[2] = (-(sqrt(ri*ri - y*y) - (k-1)*(delta/2))*sin(phiRad) + cos(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad)))+ (targetp[2] - cg[2]);
      
      if(ftot(xyz)>=1)
      {
        if(ftot(xyzmax)>=1 && ftot(xyzmin)<=1)
        {
        Surface = Surface + (delta)*(finalfield[i+1][1] - finalfield[i][1]);
        SurfaceOut = SurfaceOut + (delta)*(finalfield[i+1][1] - finalfield[i][1]);
        printf(" La Surface ajouté pour la ligne de champ %f est %f \n",Champs[p][0],Surface);
        printf(" La Surface Out ajouté pour la ligne de champ %f est %f \n",Champs[p][0],SurfaceOut);
        break;
        }
        else{break;}
      }
      else{k=k+1;}
      
      }
      j=j+1;
      }
      }
/*******************************************************************************/
/* Pour x inférieur à 0 */
/*******************************************************************************/
   printf("Intersection du calcul surfacique pour x inférieur à 0\n");
   for(i = 0; i < finalsize - 1;i++)
   {
      zi = finalfield[i][1];
      ri = finalfield[i][2];

   j=1;
   while(j*(delta/2)< 2* ri)
   {
   y = ri - j*(delta/2);
   k=0;
      while(k*(delta/2) <= 2*sqrt(ri*ri - y*y))
      {
      x = - sqrt(ri*ri - y*y) + k*(delta/2);
                          
      xyz[0] = x * cos(phiRad) + sin(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad))+ (targetp[0] - cg[0]) ;
      xyz[1] = y*cos(thetaRad) - zi*sin(thetaRad)+ (targetp[1] - cg[1]);
      xyz[2] = (-x*sin(phiRad) + cos(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad)))+ (targetp[2] - cg[2]);
      
      xyzmax[0] = (-sqrt(ri*ri - y*y) + (k+1)*(delta/2)) * cos(phiRad) + sin(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad))+ (targetp[0] - cg[0]);
      xyzmax[1] = y*cos(thetaRad) - zi*sin(thetaRad)+ (targetp[1] - cg[1]);
      xyzmax[2] = (-(-sqrt(ri*ri - y*y) + (k+1)*(delta/2))*sin(phiRad) + cos(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad)))+ (targetp[2] - cg[2]);
      
      xyzmin[0] = (-sqrt(ri*ri - y*y) + (k-1)*(delta/2)) * cos(phiRad) + sin(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad))+ (targetp[0] - cg[0]) ;
      xyzmin[1] = y*cos(thetaRad) - zi*sin(thetaRad)+(targetp[1] - cg[1]);
      xyzmin[2] = (-(-sqrt(ri*ri - y*y) + (k-1)*(delta/2))*sin(phiRad) + cos(phiRad)*(y*sin(thetaRad) + zi*cos(thetaRad)))+ (targetp[2] - cg[2]);
      
      if(ftot(xyz)>=1)
      {
        if(ftot(xyzmax)>=1 && ftot(xyzmin)<=1)
        {
        Surface = Surface + (delta) * (finalfield[i+1][1] - finalfield[i][1]);          
        SurfaceIn = SurfaceIn + (delta) * (finalfield[i+1][1] - finalfield[i][1]);
        printf(" La Surface ajouté pour la ligne de champ %f est %f \n",Champs[p][0],Surface);
        printf(" La Surface In ajouté pour la ligne de champ %f est %f \n",Champs[p][0],SurfaceIn);
        break;
        }
        else{break;}
      }
      else{k=k+1;}
      }
      j=j+1;
      }
      }
      fPrintF(g, " %f, %f, %f, %f \n",Champs[p][0], Surface, SurfaceIn, SurfaceOut);
      p = p + M;
      }
      fflush(g);
      fclose(g);  
      return;
}



