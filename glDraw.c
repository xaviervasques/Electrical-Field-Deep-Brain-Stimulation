/*******************************************************************************
Xavier Vasques : X.V.
Last Revised : 29/03/2006
*******************************************************************************/

#include  <stdio.h>
#include  <string.h>
#include  <stdlib.h>
#include  <math.h>
#include  <GL/glut.h>

#include  "glDraw.globalloc.h"
#include  "dataPoints.glob.h"
#include  "grille.glob.h"
#include  "glDraw.h"
#include  "glDraw.Pallidus.h"
#include  "glDraw.Electrode.h"
#include  "Field.h"
#include  "StereoFrame.h"

static int     rotx = 0, roty = 0, rotz=0;                       /* rotations */
static double  transx = 0.0, transy = 0.0, transz = 0.0;      /* translations */
static double  scal = 1.0;                                  /* gestion du zoom*/

int       font = (int)GLUT_BITMAP_8_BY_13;
int       bitmapHeight = 13;
int       frame = 0, time = 0, timebase = 0;
char      s[80];

static char    srcfilenamp[] = "glDraw";

/******************************************************************************/
void      glInit() 
{ 
  GLfloat mat_specular[]     = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat mat_shininess[]    = { 100 };
  GLfloat light_position[]   = { 0.0, 0.0, 1.0, 0.0 };
  GLfloat diffuseMaterial[4] = { 0.5, 0.5, 0.5, 1.0 };
  static char    prognamp[] = "glInit";

  diffuseMaterial[0] = diffuseMaterial[1] = diffuseMaterial[2] = 0.8;
  NumList = glGenLists( 1 );

  fprintf(stderr, "%s::%s  BEGIN glDrawTriangles CALL \n", srcfilenamp, prognamp);
  glDrawTriangles();
  fprintf(stderr, "%s::%s        glDrawTriangles DONE \n", srcfilenamp, prognamp);

  glPixelStorei (GL_UNPACK_ALIGNMENT, 1);
  glClearColor (0.0, 0.0, 0.0, 0.0);
  glShadeModel (GL_FLAT);

  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, diffuseMaterial);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);

  fprintf(stderr, "%s::%s  GL light choices \n", srcfilenamp, prognamp);
  glLightfv(GL_LIGHT0, GL_POSITION, light_position);
  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  glEnable(GL_DEPTH_TEST);
  glColorMaterial(GL_FRONT_AND_BACK, GL_DIFFUSE);
  glEnable(GL_COLOR_MATERIAL);
  
  glEnable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  glBlendFunc(GL_ONE_MINUS_SRC_ALPHA, GL_SRC_ALPHA);
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
  glutSwapBuffers();
                
  return;
}
/******************************************************************************/


/******************************************************************************/
void      renderBitmapString(float x, float y, void *fontchoice, char *string)
{ char     *cp;

  /*  set position to start drawing fonts */                                        
  glRasterPos2f(x, y);                
  
  /* loop all the characters in the string */
  for (cp = string;  *cp != '\0';  cp++)
  { 
  glutBitmapCharacter(fontchoice, *cp);
  }
  return;
}
/******************************************************************************/

/******************************************************************************/
void      display(void)
{ int       i;
  double   *currentptp;
 
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  
  glLoadIdentity();

  glTranslatef( (GLfloat)transx, 0.0f, 0.0f);
  glTranslatef( 0.0f, (GLfloat)transy, 0.0f);
  glTranslatef( 0.0f, 0.0f, (GLfloat)transz);
  /* Rotation selon Touche X */
  glRotatef ( (GLfloat)rotx, 1.0, 0.0, 0.0);
  /* Rotation selon Touche Y */
  glRotatef ( (GLfloat)roty, 0.0, 1.0, 0.0);
  /* Rotation selon Touche Z */
  glRotatef ( (GLfloat)rotz, 0.0, 0.0, 1.0);
  /* Gestion du zoom (Touche: P et M) */
  glScalef( (GLfloat)scal, (GLfloat)scal, (GLfloat)scal);
  
  glPushMatrix();
  
  /* Tracer les triangles de la Liste */
  glColor4f(0.8f, 0.8f, 0.0f, 0.5f);
  glCallList( NumList );

  /* Tracer les points */
  glPointSize(3.0f);
  glBegin(GL_POINTS);
  glColor3f(0.0, 0.0, 0.0);
  glColor3f(0.8, 0.8, 0.8);
  for (i=0;  i < NPTS;  i++)
  { 
  currentptp = ptsCentresp + 3*i;
  glVertex3f( (float)currentptp[0], (float)currentptp[1], (float)currentptp[2] );
  }
  glEnd();
  
  /* Repère cartésien dans la structure*/
   
   /* Axe des X: rouge */
   /*glColor3f(1.0,0.0,0.0);   
   glBegin(GL_LINES);
     glVertex3f(0.0,  0.0, 0.0);
     glVertex3f(10.0,  0.0, 0.0);
   glEnd();
   glPushMatrix();   
   glTranslatef(10.0,  0.0, 0.0);
   glRotatef(90,0.0,1.0,0.0);
   glColor3f(1.0, 0.0, 0.0);
   glutSolidCone(1,4,50,50);
   glPopMatrix();*/
   
   /*Axe des Y: vert*/
   /*glColor3f(0.0,1.0,0.0);
   glBegin(GL_LINES);
     glVertex3f(0.0,  0.0, 0.0);
     glVertex3f(0.0,  10.0, 0.0);
   glEnd();
   glPushMatrix();
   glTranslatef(0.0,  10.0, 0.0);
   glRotatef(-90,1.0,0.0,0.0);
   glColor3f(0.0, 1.0, 0.0);
   glutSolidCone(1,4,50,50);
   glPopMatrix();*/
   
   /*AXE des Z: bleu*/
   /*glColor3f(0.0, 0.0, 1.0);
   glBegin(GL_LINES);
   glVertex3f(0.0,  0.0, 0.0);
   glVertex3f(0.0,  0.0, 10.0);
   glEnd();
   glPushMatrix();
   glTranslatef(0.0,  0.0, 10.0);
   glColor3f(0.0, 0.0, 1.0);
   glutSolidCone(1,4,50,50);
   glPopMatrix(); */

  glTranslatef(targetp[0] - cg[0], targetp[1] - cg[1], targetp[2] - cg[2]);
  
  glPushMatrix();
  
  /* Display of the isofield lines */
  glDisable(GL_LIGHTING);
   IsoField(thetaEl, phiEl);
  glEnable(GL_LIGHTING);
  
  /* Display of the electrode*/
  glDisable(GL_LIGHTING);
  glDrawElectrode1(thetaEl, phiEl);
  glEnable(GL_LIGHTING);
  
    
  /* Display of the stereotactic frame */
  glDisable(GL_LIGHTING);
  StereoFrame();
  glEnable(GL_LIGHTING);               
  
  glPopMatrix();
  
  glPopMatrix();

  /* Time Frame */
  frame++;
  time = glutGet(GLUT_ELAPSED_TIME);
  if(time - timebase > 1000) 
  { sprintf(s, "frame=%d/sec  time=%d sec  timebase=%d sec", frame, time/1000, timebase/1000);
    timebase = time;		
    frame = 0;
  }

  glLoadIdentity();
  
  /* Comments Display */
  glDisable(GL_LIGHTING);
  glColor3f(0.0f, 1.0f, 1.0f);
  renderBitmapString(-grilleDimensiondp[0] * 0.48, grilleDimensiondp[1] * 0.46,
                                   (void *)font, " Deep brain stimulation ");
  renderBitmapString(-grilleDimensiondp[0] * 0.48, grilleDimensiondp[1] * 0.42,
                                   (void *)font, s);

  glEnable(GL_LIGHTING);
  
  glutSwapBuffers();
  return;
}
/******************************************************************************/

/******************************************************************************/

void      reshape(int w, int h)
{ float     orthox, orthoy, orthoz;

  orthox = 0.5 * grilleDimensiondp[0] * (GLfloat)w/(GLfloat)h;
  orthoy = 0.5 * grilleDimensiondp[1];
  orthoz = 0.5 * grilleDimensiondp[2];

  glViewport (0, 0, (GLsizei) w, (GLsizei) h); 
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity ();

  glOrtho(-orthox, orthox, -orthoy, orthoy, -orthoz, orthoz);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  return;
}
/******************************************************************************/

/******************************************************************************/
/* Rotation: X,x,Y,y,Z,z. Zoom: p,m. Retour Initial: i,j,l */
void      keyboard(unsigned char key, int x, int y)
{

fprintf(stderr, "you typed character %2X = %c \n", key, key);

  switch (key)
  {
    case 'x':
         rotx = (rotx + 5) % 360;
         glutPostRedisplay();
         break;

    case 'X':
         rotx = (rotx - 5) % 360;
         glutPostRedisplay();
         break;

    case 'y':
         roty = (roty + 5) % 360;
         glutPostRedisplay();
         break;

    case 'Y':
         roty = (roty - 5) % 360;
         glutPostRedisplay();
         break;

    case 'z':
         rotz = (rotz + 5) % 360;
         glutPostRedisplay();
         break;

    case 'Z':
         rotz = (rotz - 5) % 360;
         glutPostRedisplay();
         break;

    case 'i':
    case 'I':
         rotx = roty = rotz =0;
         glutPostRedisplay();
         break;

    case 'j':
    case 'J':
         transy = transx = transz =0;
         glutPostRedisplay();
         break;
    case 'p':
    case 'P':
         scal = (scal + 0.05) ;
         glutPostRedisplay();
         break;
    case 'm':
    case 'M':
         scal = (scal - 0.05);
         glutPostRedisplay();
         break;
    case 'l':
    case 'L':
         scal=1.0;
         glutPostRedisplay();
         break;
  
    case 27:                                                   /** X1B = ESC char **/
    case 'q':
         exit(0);
         break;

    default:
         break;
  }
  return;
}
/******************************************************************************/

/******************************************************************************/
void      move(int key, int x, int y)
{

  switch (key) 
  {
    case GLUT_KEY_UP:
        transy = transy + 0.02 * grilleDimensiondp[1] ;
        break;

    case GLUT_KEY_DOWN:
        transy = transy - 0.02 * grilleDimensiondp[1] ;
        break;

    case GLUT_KEY_LEFT:
        transx = transx - 0.02 * grilleDimensiondp[0] ;
        break;

    case GLUT_KEY_RIGHT:
        transx = transx + 0.02 * grilleDimensiondp[0] ;
        break;
  }
  glutPostRedisplay();
  return;
}
/******************************************************************************/

