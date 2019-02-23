/* G.M.                           */
/* XV : Last Revised : 10/04/2006 */

#include  <stdio.h>
#include  <string.h>
#include  <stdlib.h>
#include  <math.h>
#include  <GL/glut.h>

#include  "glDraw.Electrode.h"

static char    srcfilenamp[] = "glDraw.Electrode";
static int     DEBUG = 0;

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      glDrawElectrode1(double theta, double phi)
{ float     electrodeRadius = 0.635, insulatorZ = 0.5, conductorZ = 1.5;
  float     insulatorWireZ = 26.25, coneZ = 1;
  float     zIC = 0.25 + 1.5 + 0.5 +1.5 + 0.5;                   /* altitude Cone **/
  static GLUquadricObj    *isol1p;
  static char    prognamp[] = "glDrawElectrode";


  isol1p = gluNewQuadric();
  gluQuadricDrawStyle(isol1p, GLU_LINE);
  glPushMatrix();

  glRotatef ( (GLfloat)phi,   0.0, 1.0, 0.0);
  glRotatef ( (GLfloat)theta, 1.0, 0.0, 0.0);
  
if(DEBUG)fprintf(stderr,"%s::%s BEGIN\n", srcfilenamp, prognamp);
  
  /*glLineWidth( 3.0);
  glColor3f(0.0, 0.0, 1.0);
  glBegin(GL_LINE_STRIP);
    glVertex3f(20.0,  0.0, 0.0);
    glVertex3f( 0.0,  0.0, 0.0);
    glVertex3f( 0.0, 10.0, 0.0);
  glEnd();
*/

          /** cone insulator **/
  glTranslatef(0.0, 0.0, zIC);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, 0.0, coneZ, 100,100);

         /** insulator 0 **/
  glTranslatef(0.0, 0.0, -insulatorZ);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorZ, 100,100);

         /** contact 0 **/
  glTranslatef(0.0, 0.0, -conductorZ);
  glColor3f(1.0, 0.0, 0.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, conductorZ, 100,100);

         /** insulator 1 **/
  glTranslatef(0.0, 0.0, -insulatorZ);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorZ, 100,100);

         /** contact 1 **/
  glTranslatef(0.0, 0.0, -conductorZ);
  glColor3f(1.0, 0.0, 0.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, conductorZ, 100,100);

         /** insulator 2 **/
  glTranslatef(0.0, 0.0, -insulatorZ);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorZ, 100,100);

         /** contact 2 **/
  glTranslatef(0.0, 0.0, -conductorZ);
  glColor3f(1.0, 0.0, 0.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, conductorZ, 100,100);

         /** insulator 3 **/
  glTranslatef(0.0, 0.0, -insulatorZ);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorZ, 100,100);

         /** contact 3 **/
  glTranslatef(0.0, 0.0, -conductorZ);
  glColor3f(1.0, 0.0, 0.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, conductorZ, 100,100);

        /** fil isolant **/
  glTranslatef(0.0, 0.0, -insulatorWireZ);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorWireZ, 100,100);
  
  glPopMatrix();
  return;
}
/******************************************************************************/
/******************************************************************************/

void      glDrawElectrode2(double theta, double phi)
{ float     electrodeRadius = 0.635, insulatorZ = 0.5, conductorZ = 1.5;
  float     insulatorWireZ = 26.25, coneZ = 1;
  float     zIC = 0.25 + 1.5 + 0.5 +1.5 + 0.5;                   /* altitude Cone **/
  static GLUquadricObj    *isol1p;
  static char    prognamp[] = "glDrawElectrode";


  isol1p = gluNewQuadric();
  gluQuadricDrawStyle(isol1p, GLU_LINE);
  glPushMatrix();

  glRotatef ( (GLfloat)phi,   0.0, 1.0, 0.0);
  glRotatef ( (GLfloat)theta, 1.0, 0.0, 0.0);


          /** cone insulator **/
  glTranslatef(0.0, 0.0, zIC);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, 0.0, coneZ, 100,100);

         /** insulator 0 **/
  glTranslatef(0.0, 0.0, -insulatorZ);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorZ, 100,100);

         /** contact 0 **/
  glTranslatef(0.0, 0.0, -conductorZ);
  glColor3f(1.0, 0.0, 0.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, conductorZ, 100,100);

         /** insulator 1 **/
  glTranslatef(0.0, 0.0, -insulatorZ);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorZ, 100,100);

         /** contact 1 **/
  glTranslatef(0.0, 0.0, -conductorZ);
  glColor3f(1.0, 0.0, 0.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, conductorZ, 100,100);

         /** insulator 2 **/
  glTranslatef(0.0, 0.0, -insulatorZ);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorZ, 100,100);

         /** contact 2 **/
  glTranslatef(0.0, 0.0, -conductorZ);
  glColor3f(1.0, 0.0, 0.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, conductorZ, 100,100);

         /** insulator 3 **/
  glTranslatef(0.0, 0.0, -insulatorZ);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorZ, 100,100);

         /** contact 3 **/
  glTranslatef(0.0, 0.0, -conductorZ);
  glColor3f(1.0, 0.0, 0.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, conductorZ, 100,100);

        /** fil isolant **/
  glTranslatef(0.0, 0.0, -insulatorWireZ);
  glColor3f(0.0, 0.0, 1.0);
  gluCylinder(isol1p, electrodeRadius, electrodeRadius, insulatorWireZ, 100,100);
  
  glPopMatrix();
  return;
}
/******************************************************************************/
/******************************************************************************/
