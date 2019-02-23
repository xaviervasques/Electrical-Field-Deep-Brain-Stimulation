/*  ../libmy/minmaxd.c                                                        */
/*   Mennessier Gerard             960315                                     */
/*   Last revised M.G.             990415                                     */

#include  <stddef.h>
#include  "minmaxd.h"

/******************************************************************************/
/* minmaxd(n,&x,&xmin,&xmax)   x array de n real*8 {x(0),...,x(n-1)}          */
/*                                xmin,xmax=min,max  des x(i)                 */
/******************************************************************************/
void      minmaxd(size_t nx, double *xp, double *xminp, double* xmaxp)
{ size_t    i;
  double    xmin,xmax,x;
  
  xmax = xmin = *xp;  xp++;
  for(i = 1;  i < nx;  i++, xp++)
  { x = *xp;
    if(x < xmin)     { xmin = x;}
    else if(x > xmax){ xmax = x;}
  }
  *xminp = xmin;  *xmaxp = xmax;
  return;
}
/******************************************************************************/
/******************************************************************************/
