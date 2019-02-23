/*  libc/recipe2double/recipe.lubksb.c                                        */
/** VECTOR and MATRIX INDICES FROM 0 to N-1                                  **/
/** VECTOR and MATRIX : DOUBLE version                                       **/

/******************************************************************************/
/*                                                                            */
/* To solve linear equation  a[i][j] x[j] = b[i]                              */
/* USE: double **a, *b, d (for determinant)                                   */
/*      int      n, *indx                                                     */
/*      ...                                                                   */
/*      ludcmp(a,n,indx, &d);                                                 */
/*      lubksb(a,n,indx,b);                                                   */
/******************************************************************************/
void      lubksb(double **a, int n, int *indx, double *b)
{ int       i,ii=-1,ip,j;
  double    sum;

  for (i=0; i<n; i++) 
  { ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii >= 0){ for (j=ii; j<=i-1; j++) sum -= a[i][j]*b[j];}
    else if (sum) ii = i;
    b[i] = sum;
  }
  for (i=n-1; i>=0; i--) 
  { sum = b[i];
    for (j=i+1; j<n; j++) sum -= a[i][j]*b[j];
    b[i] = sum/a[i][i];
  }
}
/* (C) Copr. 1986-92 Numerical Recipes Software 'i<. */
/******************************************************************************/
/******************************************************************************/
