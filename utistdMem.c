/*   ../libmy/utistdMem.c                                                     */
/*   Mennessier Gerard             940822                                     */
/*   Last revised M.G.             950824                                     */

#include  <stddef.h>                                          /** stddef.h (ANSI) **/
#include  "utistdMem.h"

/******************************************************************************/
/* mem  FAMILY : STANDARD                                                     */
/* see  GNU glibc source/sysdeps/generic/... for better algorithms            */
/******************************************************************************/

#ifndef  DEF_MAPMEM

/******************************************************************************/
/*        memChr(mp,c,n)                                                      */
/* ~= memchr  string.h (ANSI)                                                 */
/* search for the first occurence of a byte (in the n bytes mp)               */
/* equal to the value obtained by converting c to a char                      */
/* return a pointer to this byte if found, or NULL                            */
/******************************************************************************/
void     *memChr(const void *mp,int c,size_t n)
{ char     *cp, *cps, cc;

  cp=(char *)mp; cps=cp+n; cc=(char)c;
  for( ; cp<cps; cp++)
  { if(*cp == cc){return (void *)cp;}
  }
  return NULL;
}
/******************************************************************************/
/*        memCmp(m1p,m2p,n)                                                   */
/* ~= memcmp  string.h (ANSI)                                                 */
/* return <0 ,0, >0  if  m1p <, =, > m2p                                      */
/* where comparison is done on the first n characters ONLY                    */
/******************************************************************************/
int memCmp(const void *m1p, const void  *m2p, size_t n)
{ size_t    i;
  int       cd;
  char     *c1p, *c2p ;

  c1p=(char *)m1p; c2p=(char *)m2p;
  for(i=0; i<n; i++)
  { cd= (unsigned char)*c1p - (unsigned char)*c2p ;
    if(cd==0){ c1p++;c2p++; }
    else
    { if(cd>0){return cd;}
      return -1;
    }
  }
  return 0;
}
/******************************************************************************/
/*        memCpy(mfp,mip,n)                                                   */
/* ~= memcpy  string.h (ANSI)                                                 */
/* copy n bytes from mip to mfp                                               */
/* return mfp                                                                 */
/* WARNING : NO CARE in OVERLAP case !!!  memMove is safer!                   */
/******************************************************************************/
void     *memCpy(void *mfp, const void  *mip, size_t n)
{ size_t    i;
  char     *cfp, *cip ;

  cfp=(char *)mfp; cip=(char *)mip;
  for( i=0; i<n ; i++)
  { *cfp= *cip ; cfp++; cip++; }
  return mfp;
}
/******************************************************************************/
/*        memSet(mp,c,n)                                                      */
/* ~= memset  string.h (ANSI)                                                 */
/* sets the first n bytes in memory area mp                                   */
/* to the value obtained by converting c to a char                            */
/* return mp                                                                  */
/******************************************************************************/
void     *memSet(void *mp,int c,size_t n)
{ char     *cp, *cps, cc;

  cp=(char *)mp; cps=cp+n; cc=(char)c;
  for( ; cp<cps; cp++){ *cp=cc; }
  cp=mp; return (void *)cp;
}
/******************************************************************************/

#endif                                                      /** ifndef DEF_MAPMEM **/
#if  !defined(DEF_MAPMEM) || defined(SPARCM)

/******************************************************************************/
/*        memMove(mfp,mip,n)                                                  */
/* ~= memmove  string.h (ANSI)                                                */
/* copy n bytes from mip to mfp                                               */
/* return mfp                                                                 */
/* correct behaviour even if OVERLAP                                          */
/******************************************************************************/
void     *memMove(void *mfp, const void  *mip, size_t n)
{ size_t    i;
  char     *cfp, *cip ;

  cfp=(char *)mfp; cip=(char *)mip;
  if(cfp <= cip)
  { for( i=0; i<n ; i++){ *cfp= *cip ; cfp++; cip++; }
  }
  else
  { cfp=cfp+n-1; cip=cip+n-1;
    for( i=0; i<n ; i++){ *cfp= *cip ; cfp--; cip--; }
  }
  return mfp;
}
/******************************************************************************/

#endif                            /** if  !defined(DEF_MAPMEM) || defined(SPARCM) **/

/******************************************************************************/
/******************************************************************************/
