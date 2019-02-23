/*  ../libmy/utiCurve.points.c                                                */
/*  Mennessier Gerard                   20050207                              */
/*  Last revised M.G.                   20050208                              */
/*                                                                            */
 
#include  <stddef.h>
#include  <math.h>                                       /** for function  fabs() **/
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"
 
#include  "utiCurve.points.h"
#include  "utiCurve.GC.h"

static char    srcfilenamp[] = "utiCurve.points";
/******************************************************************************/
/*        cPointsPAlloc(cPoints *p, size_t nz)                                */
/******************************************************************************/
void      cPointsPAlloc(cPoints *p, size_t nz, char *progcallp)
{ static char   form1p[]= "called from %s.  malloc failed; desired double size=%d\n";
  static char   prognamp[] = "utiCurve.points::cPointsPAlloc";

  p->xp = (double *)malloc( 2*nz *sizeof(double) );
  if(p->xp == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, 2*nz);
  p->z = nz;  p->x = 0;
  *(p->xp) = 0.0;
  return;
}
/******************************************************************************/
/*        cPointsPRealloc(cPoints *p, size_t nz)                              */
/******************************************************************************/
void      cPointsPRealloc(cPoints *p, size_t neednz, size_t incrnz, char *progcallp)
{ double   *xpi;
  size_t    nzf;
  static char   form1p[]= "called from %s.  malloc failed; desired double size=%d\n";
  static char   form2p[]= "called from %s. realloc failed; desired double size=%d\n";
  static char   prognamp[] = "utiCurve.points::cPointsPRealloc";

  nzf = neednz + incrnz;
  xpi = p->xp;
  if(xpi == NULL)
  { p->xp = (double *)malloc( 2*nzf *sizeof(double) );
    if(p->xp == NULL) myErr1(-1, stderr, prognamp, form1p, progcallp, 2*nzf);
    p->z = nzf;  p->x = 0;  *(p->xp) = 0.0;
  }
  else
  { if(neednz <= p->z) return;
    p->xp = (double *)realloc(xpi, 2*nzf *sizeof(double) );
    if(p->xp == NULL) myErr1(-1, stderr, prognamp, form2p, progcallp, 2*nzf);
    p->z = nzf;
  }
  return;
}
/******************************************************************************/
/*        cPointsPrint(FILE  *bufp, cPoints *p)                               */
/******************************************************************************/
void      cPointsPrint(FILE  *bufp, cPoints *p)
{ int       i;
  double   *dp;
  static    char form1p[] = "%s::%s  cPointsp=%p, p->z=%d, p->x=%d ";
  static    char form2p[] = "{%f, %f} ";
  static    char prognamp[] = "cPointsPrint";

  fPrintF(bufp, form1p, srcfilenamp, prognamp, p, p->z, p->x);
  dp = p->xp;
  for(i = 0;  i < p->x;  i++,  dp++, dp++)
  { if(i%5 == 0) fPrintF(bufp, "\n");
    fPrintF(bufp, form2p, *dp, *(dp+1));
  }
  fPrintF(bufp, "\n");
  return;
}
/******************************************************************************/
/*        cPointsInc1(cPoints *p, double x, double y)                         */
/******************************************************************************/
void      cPointsInc1(cPoints *p, double x, double y)
{ size_t    newz, oldx;
  double   *dp;
  static    char prognamp[] = "cPointsInc1";

  oldx = p->x;  newz = oldx +1;
  cPointsPRealloc(p, newz, newz/8, prognamp);
  dp = p->xp + 2*(p->x);
  *dp++ = x;
  *dp = y;
  p->x = newz;
  return;
}
/******************************************************************************/
/*        cPointsIncNp(cPoints *p, size_t n, double *xyp)                     */
/*                                                                            */
/* Expecting n points, stored in 2*n double {x0,y0,x1,y1,x2,y2,...}           */
/******************************************************************************/
void      cPointsIncNp(cPoints *p, size_t n, double *xyp)
{ size_t    newz, oldx;
  int       i;
  double   *dp;
  static    char prognamp[] = "cPointsInc1";

  oldx = p->x;  newz = oldx +n;
  cPointsPRealloc(p, newz, newz/8, prognamp);
  dp = p->xp + 2*(p->x);
  for(i = 0;  i < n;  i++)
  { *dp++ = *xyp++;  *dp++ = *xyp++;
  }
  p->x = newz;
  return;
}
/******************************************************************************/
/******************************************************************************/
