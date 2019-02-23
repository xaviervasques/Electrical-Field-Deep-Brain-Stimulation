/* Unité de Recherche sur les Mouvements Anormaux*/
/*         Xavier Vasques : 11042006             */
/*          Last Revised : 12042006              */                                
/*      Carde de stéreotaxie et l'origine        */

#include  <stdio.h>
#include  <GL/glut.h>
#include  "StereoFrame.h"
#include  "glDraw.globalloc.h"
#include  "dataPoints.glob.h"
#include  "grille.glob.h"

#include  "glDraw.h"
#include  "glDraw.Pallidus.h"
#include  "glDraw.Electrode.h"

void StereoFrame()
{

glRasterPos3f(5.0-targetp[0],  0.0-targetp[1], 36.0-targetp[2]);
glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, 'Z');

glColor3f(0.0,1.0,0.0);
glRasterPos3f(-5.0-targetp[0],  36.0-targetp[1], 0.0-targetp[2]);
glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, 'Y');

glColor3f(1.0,0.0,0.0);
glRasterPos3f(36.0-targetp[0],  -5.0-targetp[1], 0.0-targetp[2]);
glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, 'X');

     glColor3f(1.0,1.0,1.0);
     glPushMatrix();
     
     /* Plaque avant*/
     glBegin(GL_LINE_STRIP);
          glVertex3f(40-targetp[0],215-targetp[1],160-targetp[2]);
          glVertex3f(40-targetp[0],215-targetp[1],40-targetp[2]);
          glVertex3f(160-targetp[0],215-targetp[1],160-targetp[2]);
          glVertex3f(160-targetp[0],215-targetp[1],40-targetp[2]);
     glEnd();
     
     /*Plaque arriere*/
     glBegin(GL_LINE_STRIP);
          glVertex3f(40-targetp[0],-15-targetp[1],160-targetp[2]);
          glVertex3f(40-targetp[0],-15-targetp[1],40-targetp[2]);
          glVertex3f(160-targetp[0],-15-targetp[1],160-targetp[2]);
          glVertex3f(160-targetp[0],-15-targetp[1],40-targetp[2]);
     glEnd();
     
     /* Plaque gauche */
     glBegin(GL_LINE_STRIP);
          glVertex3f(5-targetp[0],40-targetp[1],160-targetp[2]);
          glVertex3f(5-targetp[0],40-targetp[1],40-targetp[2]);
          glVertex3f(5-targetp[0],160-targetp[1],160-targetp[2]);
          glVertex3f(5-targetp[0],160-targetp[1],40-targetp[2]);
     glEnd();
     
     /* Plaque droite */
     glBegin(GL_LINE_STRIP);
          glVertex3f(195-targetp[0],40-targetp[1],40-targetp[2]);
          glVertex3f(195-targetp[0],40-targetp[1],160-targetp[2]);
          glVertex3f(195-targetp[0],160-targetp[1],40-targetp[2]);
          glVertex3f(195-targetp[0],160-targetp[1],160-targetp[2]);
     glEnd();
     glPopMatrix();
  
   
   /* AXES d'origine*/
   
   /* Axe des X*/
   /* rouge*/
   
   glColor3f(1.0,0.0,0.0);    
   glBegin(GL_LINES);
     glVertex3f(0.0-targetp[0],  0.0-targetp[1], 0.0-targetp[2]);
     glVertex3f(30.0-targetp[0],  0.0-targetp[1], 0.0-targetp[2]);
   glEnd();
   glPushMatrix();   
   glTranslatef(30.0-targetp[0],  0.0-targetp[1], 0.0-targetp[2]);
   glRotatef(90,0.0,1.0,0.0);
   glColor3f(1.0, 0.0, 0.0);
   glutSolidCone(2,8,50,50);
   glPopMatrix();
   
   /*Axe des y*/
   /* vert*/
   glColor3f(0.0,1.0,0.0);
   glBegin(GL_LINES);
     glVertex3f(0.0-targetp[0],  0.0-targetp[1], 0.0-targetp[2]);
     glVertex3f(0.0-targetp[0],  30.0-targetp[1], 0.0-targetp[2]);
   glEnd();
   glPushMatrix();
   glTranslatef(0.0-targetp[0],  30.0-targetp[1], 0.0-targetp[2]);
   glRotatef(-90,1.0,0.0,0.0);
   glColor3f(0.0, 1.0, 0.0);
   glutSolidCone(2,8,50,50);
   glPopMatrix();
   
   /*AXE des z*/
   /* bleu*/
   glColor3f(0.0, 0.0, 1.0);
   glBegin(GL_LINES);
   glVertex3f(0.0-targetp[0],  0.0-targetp[1], 0.0-targetp[2]);
   glVertex3f(0.0-targetp[0],  0.0-targetp[1], 30.0-targetp[2]);
   glEnd();
   glTranslatef(0.0-targetp[0],  0.0-targetp[1], 30.0-targetp[2]);
   glPushMatrix();
   glColor3f(0.0, 0.0, 1.0);
   glutSolidCone(2,8,50,50);
   glPopMatrix(); 
   

 
}


