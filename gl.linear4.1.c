/*  URMAE/orientHaut/linear4.GL.V2/gl.linear4.1.c                             */
/*  Mennessier Gerard                 20010613                                */
/*  Last Revised : G.M.               20030513                                */

#include  <stddef.h>
#include  <math.h>
#include  <GL/glut.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiCurve.set.h"
#include  "utiCurve.2GL.h"

#include  "gl.linear4.1.globalloc.h"
#include  "gl.linear4.glob.h"

#include  "gl.linear4.1.h"

#include  "gm.drawstate.E.h"
#include  "pallidus.draw.h"
/******************************************************************************/
/*                                                                            */
/*             Drawing the Rho-Z plane of the Electrode Frame                 */
/*                                                                            */
/*                                WINDOW TITLE                                */
/*                        Electrode Frame RHO-Z window                        */
/*                                                                            */
/*  When type 'z' key, then use the (x,y) mouse position                      */
/*                     to set section "z E altitude"                          */
/******************************************************************************/

/******************************************************************************/
/*                                                                            */
/*  gmGl1InitWin(void)                                                        */
/*                                                                            */
/******************************************************************************/
void      gmGl1InitWin(void)
{ 
  gmDrawStateESetLimfromRhoZScaleRatio(1.0);

  glClearColor(1.0, 1.0, 1.0, 1.0);
  glShadeModel(GL_FLAT);
  glutSetWindowTitle("Electrode Frame RHO-Z window");

  win1ew = 0.5 * fullw;
  win1eh = 0.3 * fullh;
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmGl1Display(void)                                                        */
/******************************************************************************/
void      gmGl1Display(void)
{ int       myWIN = myWIN_1, win, wini, w, h;
  cSetVec  *csetvp;
  lcSegVec *pallidusLcSgVp;
  static    char    prognamp[] = "gl.linear4.1::gmGl1Display";

  myWIN = myWIN_1;
  wini = glutGetWindow();
  win = winp[myWIN];
  glutSetWindow(win);
  w = glutGet(GLUT_WINDOW_WIDTH);  h = glutGet(GLUT_WINDOW_HEIGHT);
fPrintF(stderr,"%s  myWIN=%d, win=%d, wini=%d, w=%d, h=%d\n",
                                                   prognamp, myWIN, win, wini, w, h);
  glClearColor(1.0, 1.0, 1.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);

  glLoadIdentity();
  glColor3f(0.0, 0.0, 0.0);
  glLineWidth(1.0);
                                                                 /** level curves **/
  gmDrawStateESetRZallCSet();
  csetvp = gmDrawStateEGetRZallCSetp();
  GLgraphSet(win, csetvp);
                                                          /** pallidus boundaries **/
  pallidusLcSgVp = pallidusERZGetlcSeg();
/*
  GLgraphLcSeg(win, pallidusLcSgVp);
*/

  glFlush();
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmGl1Reshape(int w, int h)                                                */
/******************************************************************************/
void      gmGl1Reshape(int w, int h)
{ int       myWIN = myWIN_1, win, wini;
  double    clipZp[2], clipRhop[2];
  double   *limRhop, *limZp;
  static    char      prognamp[] = "gl.linear4.1::gmGl1Reshape";

  myWIN = myWIN_1;
  wini = glutGetWindow();
  winwp[myWIN] = w;  winhp[myWIN] = h;
  win = winp[myWIN];
  glutSetWindow(win);
fprintf(stderr,"%s  myWIN=%d, win=%d, wini=%d, w=%d, h=%d\n",
                                                   prognamp, myWIN, win, wini, w, h);
  glViewport (0, 0, (GLsizei) w, (GLsizei) h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  limRhop = gmDrawStateEGetLimRho();
  limZp   = gmDrawStateEGetLimZ();
  clipRhop[0] = limRhop[0];  clipRhop[1] = limRhop[1];
  clipZp[0] = limZp[0];      clipZp[1] = limZp[1];
  glOrtho(clipRhop[0], clipRhop[1], clipZp[0], clipZp[1], -1.0, 1.0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  /*glutPostRedisplay();*/
  printf("************ Reshape 1 calling glFlush*******\n");
  glFlush();
  return;
}
/******************************************************************************/
/*                                                                            */
/*  gmGl1KB(unsigned char key, int x, int y)                                  */
/*                                                                            */
/*  When type 'z' key, then use the (x,y) mouse position                      */
/*                     to set section "z E altitude"                          */
/******************************************************************************/
void      gmGl1KB(unsigned char key, int x, int y)
{ double    z, rho;
  double    deltaz, deltarho;
  double   *limRhop, *limZp;
  int       myWIN = myWIN_1, win;
  int       w, h;
  static    char      prognamp[] = "gl.linear4.1::gmGl1KB";

  if(key != 'z') return;

  win = winp[myWIN];
  w = winwp[myWIN];  h = winhp[myWIN];

  limRhop = gmDrawStateEGetLimRho();
  limZp   = gmDrawStateEGetLimZ();
  deltarho = limRhop[1] - limRhop[0]; 
  deltaz   = limZp[1]   - limZp[0];  
  rho = limRhop[0] + (deltarho * x)/w;
  z   = limZp[1]   - (deltaz * y)/h ;
fPrintF(stderr,"%s  key=%d, x=%d, y=%d, rho=%f, z=%f \n", prognamp, key,x,y,rho,z);
  gmDrawStateESetzh(z);
  gmGl1Display();
  
  return;
}
/******************************************************************************/
/******************************************************************************/
