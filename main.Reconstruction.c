#include  'stddef.h'
#include  'stdio.h'
#include  'stdlib.h'
#include  'math.h'
#include  'GL/glut.h'

#include  'utistdIO.h'
#include  'utistdErr.h'
#include  'utistdStr.h'

#include  'main.linear4.gl.V.def.h'


#include  'gl.linear4.h'
#include  'gm.linear4.initsolve.h'
#include  'gm.drawstate.h'
#include  'gm.drawstate.threshold.h'
#include  'gm.drawstate.E.h'
#include  'gm.drawstate.ST.h'
#include  'pallidus.draw.h'
#include  'pallidus.geom.h'

#include  'dataPoints.globalloc.h'
#include  'dataFile.h'
#include  'funcBasis.h'
#include  'syslin.h'
#include  'grille.h'
#include  'grille.glob.h'
#include  'glDraw.h'
#include  'main.Reconstruction.h'
#include  'funcTot.h'

#include "main.linear4.gl.V.def.h"

#define   VERSION_  "Décembre 2006"

void      choiceGridParameters();

/******************************************************************************
Programme principal
Remarque : Changement de la dimension de la grille : main.Reconstruction.h
******************************************************************************/

int main(int argx, char** argvpp)
{ 
  int       status;
  int       i;
  double    rayonc, rayon;
  
  static char    form1p[] = "arg[%d] = %s \n";
  static char    VERSION[] = VERSION_ ;
  static char    titleFormatp[] = "***** %s : Version %s *****\n" ;
  static char    prognamp[] = "main.reconstruction";
  
    
  
  fPrintF(stdout, titleFormatp, prognamp, VERSION);
  for(i = 0;  i < argx; i++) fPrintF(stderr, form1p, i, argvpp[i]);

  status = dataFileRead("bidon");
  dataTranslate();
  
  /** CHOIX  RAYON  FONCTIONS DE BASE **/
  rayonc = 110 * (ptsCentresMax - ptsCentresMin)/sqrt( (double)NPTS );
  rayon = rayonc;
  printf("******************************************************************\n");
  
  fprintf(stdout, "  rayon R calcule = %f \n  rayon R choisi = %f ptsCentresMax = %f et min %f\n",
                 rayonc, rayon,ptsCentresMax,ptsCentresMin);
                 
  
  
  fbasisSetR(rayon);

  /** CALCUL DES POIDS **/
  dblPVecRealloc(&weightVec, NPTS, 10);
  syslin(NPTS, ptsCentresVec.p, weightVec.p);

  /** CHOIX des PARAMETRES de la GRILLE **/
  choiceGridParameters();
  
  /** CALCUL de la fonction sur tous les points de la grille **/
  grilleValues();

  
  mainpara(argx, argvpp);
  glutInit(&argx, argvpp);
  
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize (300, 300);
  glutInitWindowPosition (400, 400);
  glutCreateWindow("X.V.");

  glInit();

  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutSpecialFunc(move);
  glutDisplayFunc(display);
  /*glutIdleFunc(display);*/
  
  glutMainLoop();
  
   
  exit(0);
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/

/* CHOIX des PARAMETRES de la GRILLE */
/*Il faut fixer la taille de la grille et la dimension physique de la grille : dans main.Reconstruction.h */

void      choiceGridParameters()
{ double    typicalSizec;
  double    dimPhys;
  /*int       grilleSize = 51;*/
  double    madim=0;
  
  /*CHOIX typicalSize */
  typicalSizec = fabs(ptsCentresMin);
  if(ptsCentresMax > typicalSizec) typicalSizec = ptsCentresMax;
  typicalSizec *= 1.5;

  typicalSize = typicalSizec;
/*
  typicalSize = 30.0;
  typicalSize = 10.0;
*/
  fprintf(stdout, "  typicalSize calculee = %f \n  typicalSize choisie  = %f \n",
                 typicalSizec, typicalSize);

  /* CHOIX  DIMENSIONS PHYSIQUE DE GRILLE */
  dimPhys = typicalSize * 3.0;
  grilleSetDim(dimPhys);
  fprintf(stdout, "  Dimensions Physiques Communes pour Grille = %f \n", dimPhys);

  /** CHOIX  Nombre de points pour GRILLE **/
  grilleSetSizes( grilleSize, grilleSize, grilleSize);
  
  printf("la dimphy est %f\n",dimPhys);
  
}
  
/*******************************************************************************/
/*******************************************************************************/

