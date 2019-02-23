/*  libc/recipe2double/recipe.balanc.c                                        */
/** MATRIX INDICES FROM 0 to DIM-1                                           **/
/** MATRIX : DOUBLE version                                                  **/

#include  <math.h>
#define   RADIX 2.0

void      balanc(double **a, int n)
{ int       last,j,i;
  double    s,r,g,f,c,sqrdx;

  sqrdx = RADIX*RADIX;
  last=0;
  while(last == 0)
  { last=1;
    for(i=0;  i<n;  i++)
    { r =  c= 0.0;
      for(j=0;  j<n;  j++)
      { if(j != i) { c += fabs(a[j][i]);  r += fabs(a[i][j]);}
        if (c && r)
        { g=r/RADIX;
          f=1.0;
          s=c+r;
          while(c < g){ f *= RADIX;  c *= sqrdx;}
          g=r*RADIX;
          while(c > g){ f /= RADIX;  c /= sqrdx;}
          if ((c+r)/f < 0.95*s)
          { last = 0;
            g = 1.0/f;
            for(j=0;  j<n; j++) a[i][j] *= g;
            for(j=0;  j<n; j++) a[j][i] *= f;
          }
        }
      }
    }
  }
}
#undef RADIX
/* (C) Copr. 1986-92 Numerical Recipes Software 'i<. */
