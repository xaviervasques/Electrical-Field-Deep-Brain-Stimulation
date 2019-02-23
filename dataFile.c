/* X.V. Last Revised 22/08/2006 */

#include  <stddef.h>
#include  <stdio.h>
#include  <stdlib.h>

#include  "utistdErr.h"
#include  "utistdStr.h"
#include  "utiStr.h"
#include  "utiVecDbl.h"
#include  "dataPointsStereo.h"
#include  "dataPoints.glob.h"
#include  "dataFile.h"

static char srcfilenamp[] = "dataFile";

/******************************************************************************/ 
/* On imprime dans un fichier dataPointsStereo.h les coordonnées stéréotaxique 
   du contour Pallidal pour pouvoir le transferer par la suite dans la 
   parti visu                                                                 */
/******************************************************************************/
   
void StereoPoint()
{
      static char dataPointsStereo[] = "C:\\Users\\Laura Cif\\Desktop\\Projet Thèse\\Programme Pallidus\\Programme 3D Visu et Calcul\\dataPointsStereo.h";
      FILE *f;
      int i;
        
      f = fopen(dataPointsStereo,"w");
      
fPrintF(f, "#ifndef  DATAPOINTS_STEREO_H \n");
fPrintF(f, "#define  DATAPOINTS_STEREO_H \n");
fPrintF(f, " \n");
fPrintF(f,"#define NPTS_ %i\n",NPTS_);
fPrintF(f, " \n");
fPrintF(f, "double liste_points_stereo[NPTS_ * 3] = { \n");

for (i=0;  i < (NPTS_*3) ;  i=i+3)
{ 
  fPrintF(f, " %f, %f, %f,\n", liste_points_stereo[i],liste_points_stereo[i + 1],liste_points_stereo[i +2]);
  printf(" %f, %f, %f \n",liste_points_stereo[i],liste_points_stereo[i +1],liste_points_stereo[i +2]);
}
  
fPrintF(f, "}; \n");
fPrintF(f, "#endif \n");
fflush(f);
fclose(f);

return;
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/

int dataFileRead(char *fileNamp)
{ 
  int       status = 0;
  FILE     *streamp;

  static size_t  cbufz = 256;
  static char    cbufp[256];
  char          *cp, *cpi, *cVp;
  int            i = 0;                                       /** indice du point **/
  int            id = 0;                                    /** id = 0-2 ~ xX,Y,Z **/
  double         dcoorp[3];
  double        *dp;

  static char    coorNamep[3] = {'X', 'Y', 'Z'};
  static char    form0p[]   = "%s::%s  BEGIN  attempt to read file=%s \n";
  static char    form1p[]   = "  Cannot open Initial Data Points File = %s \n";
  static char    form2p[]   = "  Coma NOT found : cannot get coordinate %s \n";
  static char    form3p[]   = "%s::%s  END\n"
                              "  point number NPTS = %d, coor number = %d \n";
  static char    form4p[]   = "  first point = { %f, %f, %f}\n";
  static char    form5p[]   = "  last  point = { %f, %f, %f}\n\n";
  static char    prognamp[] = "dataFileRead";


  NPTS = NPTS_;
  ptsInitVec.p = liste_points_stereo;
  ptsInitVec.x = NPTS * 3;
  ptsInitVec.z = NPTS * 3;

  fPrintF(stdout, form3p, srcfilenamp, prognamp, NPTS, ptsInitVec.x);
  dp = ptsInitVec.p;
  fPrintF(stdout, form4p, *dp, *(dp+1), *(dp +2));
  dp += 3*(NPTS - 1);
  fPrintF(stdout, form5p, *dp, *(dp+1), *(dp +2));
  return status;
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
int       dataTranslate()
{ 
  double         dste, dCentre, nd;
  double        *dIp, *dCp, *dp;
  int            i, id;
  static char    form0p[] = "%s::%s  BEGIN \n";
  static char    form1p[] = "  Initial CentreGravite = { %f, %f, %f} \n";
  static char    form2p[] = "  Initial MinMax X = { %f, %f}" 
                                   "   MinMax Y = { %f, %f}"
                                   "   MinMax Z = { %f, %f} \n\n";
  static char    form4p[] = "  Centres first point = { %f, %f, %f}\n";
  static char    form5p[] = "  Centres last  point = { %f, %f, %f}\n";
  static char    form6p[] = "  Centres MinMax X = { %f, %f}" 
                                   "   MinMax Y = { %f, %f}"
                                   "   MinMax Z = { %f, %f} \n\n";
  static char    form7p[] = "  Centres XYZMinMax = { %f, %f} \n\n";
 
  static char    form8p[] = "%s::%s  END \n\n";
  static char    prognamp[] = "dataTranslate";


  fPrintF(stdout, form0p, srcfilenamp, prognamp);
  
  /** Calcul centre de gravite **/
  
  dIp = ptsInitVec.p;
  ptsInitMinp[0] = ptsInitMaxp[0] = *dIp++;
  ptsInitMinp[1] = ptsInitMaxp[1] = *dIp++;
  ptsInitMinp[2] = ptsInitMaxp[2] = *dIp;

  cg[0] = cg[1] = cg[2] = 0.0;
  
  dIp = ptsInitVec.p;
  
  for(i=0; i < NPTS; i++)
  { for(id = 0;  id < 3; id++, dIp++)
    { dste = *dIp;
      cg[id] += dste;
           if(dste > ptsInitMaxp[id]) ptsInitMaxp[id] = dste;
      else if(dste < ptsInitMinp[id]) ptsInitMinp[id] = dste;
    }
  }
  nd = (double)NPTS;
  cg[0] = cg[0]/nd;
  cg[1] = cg[1]/nd;
  cg[2] = cg[2]/nd;

  /*cg[0] = targetp[0];  cg[1] = targetp[1];  cg[2] = targetp[2];*/

  fprintf(stdout, form1p, cg[0], cg[1], cg[2]);
  fprintf(stdout, form2p, ptsInitMinp[0], ptsInitMaxp[0],
                          ptsInitMinp[1], ptsInitMaxp[1],
                          ptsInitMinp[2], ptsInitMaxp[2]);
                    
 

                                                                     /** Centrage **/
  dblPVecAlloc(&ptsCentresVec, 3*NPTS);
  ptsCentresp = ptsCentresVec.p;

  dIp = ptsInitVec.p;
  dCp = ptsCentresp;
  ptsCentresMinp[0] = ptsCentresMaxp[0] = *dCp++;
  ptsCentresMinp[1] = ptsCentresMaxp[1] = *dCp++;
  ptsCentresMinp[2] = ptsCentresMaxp[2] = *dCp;

  dCp = ptsCentresp;
  for(i=0; i < NPTS; i++)
  { for(id = 0;  id < 3; id++, dIp++, dCp++)
    { dste = *dIp;
      dCentre = dste - cg[id];
           if(dCentre > ptsCentresMaxp[id]) ptsCentresMaxp[id] = dCentre;
      else if(dCentre < ptsCentresMinp[id]) ptsCentresMinp[id] = dCentre;
      *dCp = dCentre;
    }
  }

  dp = ptsCentresp;
  fPrintF(stdout, form4p, *dp, *(dp+1), *(dp +2));
  dp += 3*(NPTS - 1);
  fPrintF(stdout, form5p, *dp, *(dp+1), *(dp +2));
  fprintf(stdout, form6p, ptsCentresMinp[0], ptsCentresMaxp[0],
                          ptsCentresMinp[1], ptsCentresMaxp[1],
                          ptsCentresMinp[2], ptsCentresMaxp[2]);

  ptsCentresMin = ptsCentresMinp[0];
  if(ptsCentresMinp[1] < ptsCentresMin) ptsCentresMin = ptsCentresMinp[1];
  if(ptsCentresMinp[2] < ptsCentresMin) ptsCentresMin = ptsCentresMinp[2];
  ptsCentresMax = ptsCentresMaxp[0];
  if(ptsCentresMaxp[1] > ptsCentresMax) ptsCentresMax = ptsCentresMaxp[1];
  if(ptsCentresMaxp[2] > ptsCentresMax) ptsCentresMax = ptsCentresMaxp[2];
  fprintf(stdout, form7p, ptsCentresMin, ptsCentresMax);

  fprintf(stdout, form8p, srcfilenamp, prognamp);

/*
  dblPVecFree(&ptsInitVec);
*/
 
  return 0;
}
/******************************************************************************/
/******************************************************************************/
