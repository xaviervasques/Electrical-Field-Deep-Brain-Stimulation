/*   ../libmy/echel.c                                                         */
/*   Mennessier Gerard             960320                                     */
/*   Last revised M.G.             990415                                     */

#include  <stddef.h>
#include  <math.h>
#include  "echel.h"

/******************************************************************************/
double    echel(double x)
{ double    roundto[8] = {1.0, 1.5, 2.0, 2.5, 4.0, 5.0, 8.0, 10.0 };
  double    dban,scale;
  int       iban,i;

  dban = log10(x);
  iban = (int)dban;  if(dban < 0.0){ iban--;}
  scale = pow(10.0, (double)iban);
  dban = x/scale;
  iban=7;   dban *= 2.0;
  for(i=0; i<7; i++){ if(dban < roundto[i] + roundto[i+1]){ iban=i; break;} }
  return  roundto[iban] * scale ;
}
/******************************************************************************/
/******************************************************************************/
