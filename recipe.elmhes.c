/*  libc/recipe2double/recipe.elmhes.c                                        */
/** MATRIX INDICES FROM 0 to DIM-1                                           **/
/** MATRIX : DOUBLE version                                                  **/

#include  <math.h>
#define   SWAP(g,h) {y=(g);  (g)=(h);  (h)=y;}

void      elmhes(double **a, int n)
{ int       m,j,i,nm1;
  double    y,x;

  nm1 = n-1;
  for(m=1;  m<nm1;  m++)
  { x=0.0;
    i=m;
    for(j=m;  j<n;  j++)
    { if( fabs(a[j][m-1]) > fabs(x) ){ x=a[j][m-1];  i=j;}
    }
    if(i != m)
    { for(j=m-1;  j<n;  j++) SWAP(a[i][j],a[m][j])
      for(j=0;    j<n;  j++) SWAP(a[j][i],a[j][m])
    }
    if(x)
    { for(i=m+1;  i<n;  i++)
      { if((y=a[i][m-1]) != 0.0)
        { y /= x;  a[i][m-1]=y;
          for(j=m;  j<n;  j++) a[i][j] -= y*a[m][j];
          for(j=0;  j<n;  j++) a[j][m] += y*a[j][i];
        }
      }
    }
  }
}
#undef SWAP
/* (C) Copr. 1986-92 Numerical Recipes Software 'i<. */
