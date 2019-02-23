/*  ../libmy/utiByteOrderENDIAN.c                                             */
/*  Mennessier Gerard                   20020528                              */
/*  Last Revised : G.M.                 20041011                              */

#include  <stddef.h>
#include  <stdlib.h>                                                 /** for exit **/
#include "utistdIO.h"
#include "utistdErr.h"
#include "utiByteOrderENDIAN.h"

static char    srcfilenamp[] = "utiByteOrderENDIAN";
/******************************************************************************/
/*                                                                            */
/*          utiByteOrderGet()                                                 */
/*                                                                            */
/* return pointer on 2 short (boolean) {myLSBfirst, myMSBfirst}               */
/* which caracterize the Byte Order of the RUNNING MACHINE                    */
/*                                                                            */
/** for ENDIAN byte order : little/big <=> least/most first                  **/
/**  (-)_LITTLE_ENDIAN  or  (-)_BIG_ENDIAN                                   **/
/*                                                                            */
/** if boolLSBfirst = TRUE = 1  least significant First,                     **/
/**                                                  then boolMSBfirst = 0   **/
/** if boolMSBfirst = TRUE = 1  most  significant First,                     **/
/**                                                  then boolLSBfirst = 0   **/
/*                                                                            */
/** from   From X.h                                                          **/
/**        define  LSBFirst  0                                               **/
/**        define  MSBFirst  1                                               **/
/** from   utiByteOrderENDIAN.def.h                                          **/
/**        define  myLSBFirst  0                                             **/
/**        define  myMSBFirst  1                                             **/
/******************************************************************************/
short    *utiByteOrderGet()
{ static short   LM_SBfirstp[2] = {0,0};         /** {boolLSBfirst, boolMSBfirst} **/
  unionByte4     bo;
  static char    form1p[]   = "UNEXPECTED CASE STOP ui=X%8X\n";
  static char    prognamp[] = "utiByteOrderGet";

  bo.uch[0] = '\001';  bo.uch[1] = '\002';  bo.uch[2] = '\003';  bo.uch[3] = '\004';


  if(bo.ui == 0X01020304)     { LM_SBfirstp[0] = 0;  LM_SBfirstp[1] = 1;}
  else if(bo.ui == 0X04030201){ LM_SBfirstp[0] = 1;  LM_SBfirstp[1] = 0;}
  else                        
  { myErr2(myILLEGAL_CASE, stderr, srcfilenamp, prognamp, form1p, bo.ui);
  }
  return LM_SBfirstp;
}
/******************************************************************************/
/*                                                                            */
/*  given a pointer ucp of ix bytes_PAIRS (i.e. 2*ix bytes)                   */
/*        swap the bytes inside each pair                                     */
/******************************************************************************/
void    utiSwapB2(unsigned char *ucp, size_t ix)
{ unsigned char  uc;
  int     i;

  for(i = 0;  i < ix;  i++){ uc = *ucp; *ucp = *(ucp+1); ucp++; *ucp = uc; ucp++;}
  return;
}
/******************************************************************************/
/*                                                                            */
/*  given a pointer ucip of ix bytes_PAIRS (i.e. 2*ix bytes)                  */
/*        copy the swapped pairs into ucfp                                    */
/*  byte arrays ucip and ucfp are assumed NOT to OVERLAP                      */
/******************************************************************************/
void    utiCpySwapB2(unsigned char *ucfp, unsigned char *ucip, size_t ix)
{ int     i;

  for(i = 0;  i < ix;  i++){ *(ucfp+1) = *ucip++; *ucfp = *ucip++; ucfp += 2;}
  return;
}
/******************************************************************************/
/*                                                                            */
/*  given a pointer ucp of ix bytes_TETRADS (i.e. 4*ix bytes)                 */
/*        swap the bytes inside each tetrad                                   */
/******************************************************************************/
void    utiSwapB4(unsigned char *ucp, size_t ix)
{ unsigned char  uc;
  int     i;

  for(i = 0;  i < ix;  i++)
  { uc = *ucp;  *ucp = *(ucp+3);  *(ucp+3) = uc;
    ucp++;
    uc = *ucp;  *ucp = *(ucp+1);  *(ucp+1) = uc;
    ucp += 3;
  }
  return;
}
/******************************************************************************/
/*                                                                            */
/*  given a pointer ucip of ix bytes_TETRADS (i.e. 4*ix bytes)                */
/*        copy the swapped tetrads into ucfp                                  */
/*  byte arrays ucip and ucfp are assumed NOT to OVERLAP                      */
/******************************************************************************/
void    utiCpySwapB4(unsigned char *ucfp, unsigned char *ucip, size_t ix)
{ int     i;

  for(i = 0;  i < ix;  i++)
  { *(ucfp+3) = *ucip;  *ucfp = *(ucip+3);
    ucip++;  ucfp++;
    *(ucfp+1) = *ucip;  *ucfp = *(ucip+1);
    ucip += 3;  ucfp += 3;
  }
  return;
}
/******************************************************************************/
/*                                                                            */
/*  given a byte PAIR                                                         */
/*  swap the bytes and return unsigned short int value                        */
/******************************************************************************/
unsigned short utiSwapB2ushort(unsigned char *ucp)
{ unionByte4     ob;

  ob.uch[1] = *ucp++;  ob.uch[0] = *ucp;
  return  ob.ush[0];
}
/******************************************************************************/
/*                                                                            */
/*  given a byte TETRAD                                                       */
/*  swap the bytes and return unsigned int value                              */
/******************************************************************************/
unsigned int   utiSwapB4uint(unsigned char *ucp)
{ unionByte4     ob;

  ob.uch[3] = *ucp++; ob.uch[2] = *ucp++; ob.uch[1] = *ucp++;  ob.uch[0] = *ucp;
  return  ob.ui;
}
/******************************************************************************/
/******************************************************************************/
