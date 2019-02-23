/*  URMAE/orientHaut/linear4.GL.V2/main.linear4.gl.c                          */
/*  Mennessier Gerard                 20010509                                */
/*  Last Revised : G.M.               20030509                                */

#include  <stddef.h>
#include  <math.h>
#include  <time.h>                                         /** for time and ctime **/
#include  <stdlib.h>                                               /** for system **/
#include  <GL/glut.h>

#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utistdStr.h"

#include  "main.linear4.gl.V.def.h"


#include  "gl.linear4.h"
#include  "gm.linear4.initsolve.h"
#include  "gm.drawstate.h"
#include  "gm.drawstate.threshold.h"
#include  "gm.drawstate.E.h"
#include  "gm.drawstate.ST.h"
#include  "pallidus.draw.h"
#include  "pallidus.geom.h"
#include  "main.Reconstruction.h"

/******************************************************************************/
/*                                                                            */
/* main.linear.gl                                                             */
/*                                                                            */
/******************************************************************************/
void       mainpara(int argxi, char *argvpp[])
{ static    char      prognamp[] = "main.linear4.gl" ;
  static    char      VERSION[] = "Xav Version" ;
  static    char      titleFormatp[] = "***** %s : Version %s *****\n" ;
/**  static    char    prognamp[] = "linearGL4/main.linear4.gl::MAIN";  **/
  time_t    timesec;

  static    char    sys_ps_cmd0[128] = "ps -efl ";
/*    "ps -o user,pid,ppid,pmem,rss,vsz,stime,etime,time,pri,comm -A "; */

  static    char    sys_ps_cmd1[128] = "ps -efl |grep \"menes\"";
/*    "ps -o user,pid,ppid,pmem,rss,vsz,stime,etime,time,pri,comm -A |grep \"menes\" "; */

  static    char    sys_ps_cmd2[128] = "ps -efl |grep  ";
/*    "ps -o user,pid,ppid,pmem,rss,vsz,stime,etime,time,pri,comm -A |grep    "; */

  static    char    timeformp[] = "     %s \n";
  int       argx, boolPrint1 = 1, boolPrint2 = 1;


  fPrintF(stderr,titleFormatp,prognamp,VERSION);
/* */
  /*system("uname -a");

  if(boolPrint1)
  { timesec = time(&timesec);
    fPrintF(stderr, timeformp, ctime(&timesec));
  }
  printf("etape 1 reussi\n");
  fPrintF(stderr, "  system 0 call :%s;\n",sys_ps_cmd0);
  system(sys_ps_cmd0);  fPrintF(stderr,"\n");
  fPrintF(stderr, "  system 0 was called\n");
  fPrintF(stdout, "  system 0 was called\n");


  fPrintF(stderr, "  system 1 call :%s;\n",sys_ps_cmd1);
  system(sys_ps_cmd1);  fPrintF(stderr,"\n");
  fPrintF(stderr, "  system 1 was called\n");
  fPrintF(stdout, "  system 1 was called\n");
*/
  if(boolPrint2)
  {
/*  strCpy(sys_ps_cmd2 + 68, argvpp[0]); */
    strCpy(sys_ps_cmd2 + 15, argvpp[0]);
    fPrintF(stderr, "  system 2 call :%s;\n",sys_ps_cmd2);
    system(sys_ps_cmd2);  fPrintF(stderr,"\n");
    fPrintF(stderr, "  system 2 was called\n");
    fPrintF(stdout, "  system 2 was called\n");
  }
goto JUMP1;
JUMP1:
/* */

  argx = argxi;
  gmInitBase(argx, argvpp);
  fPrintF(stderr, "\n");
  gmDrawStateInit();  
  gmDrawStateEInit();
  fPrintF(stderr, "\n");
  pallidusInitGeom(NULL);    
  pallidusInitDraw();
  gmDrawStateSTInit();
  fPrintF(stderr, "\n");
  gmDrawStateThresholdInit();
  fPrintF(stderr, "\n");
  gmGlInitGlut(&argx, argvpp);
  
}
/******************************************************************************/
/******************************************************************************/
