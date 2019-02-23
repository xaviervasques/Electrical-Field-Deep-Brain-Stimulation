/*   ../libmy/utistdStr.c                                                     */
/*   Mennessier Gerard             940822                                     */
/*   Last revised M.G.             950824                                     */

#include  <stddef.h>                                          /** stddef.h (ANSI) **/
#include  "utistdStr.h"

/******************************************************************************/
/* str  FAMILY : STANDARD                                                     */
/* see  GNU glibc source/sysdeps/generic/... for better algorithms            */
/******************************************************************************/

#ifndef  DEF_MAPSTR

/******************************************************************************/
/*        strChr(sp,c)                                                        */
/* ~= strchr  string.h (ANSI)                                                 */
/* search for the first occurence of a byte                                   */
/* equal to the value obtained by converting c to a char                      */
/* in the string sp (Zero terminated)                                         */
/* return a pointer to this byte if found, or NULL                            */
/******************************************************************************/
char     *strChr(const char *sp,int c)
{ char     *cp, cc;

  cp=(char *)sp; cc=(char)c;
  while( *cp )
  { if(*cp == cc){return cp;}
    cp++;
  }
  return NULL;
}

/******************************************************************************/
/*        strCmp(s1p,s2p)                                                     */
/* ~= strcmp  string.h (ANSI)                                                 */
/* return <0 ,0, >0  if  s1p <, =, > s2p                                      */
/******************************************************************************/
int       strCmp(const char *s1p, const char *s2p)
{ int     cd;

  while(1)
  { cd= (unsigned char)*s1p - (unsigned char)*s2p ;
    if(cd==0)
    { if(*s1p==0){return 0;}
      s1p++;s2p++;
    }
    else
    { if(cd>0){return cd;}
      return -1;
    }
  }
}

/******************************************************************************/
/*        strCpy(sfp,sip)                                                     */
/* ~= strcpy  string.h (ANSI)                                                 */
/* copy string sip (Zero terminated)  to sfp                                  */
/* sip and sfp strings should be NON-overlaping                               */
/* or  sip should be >= sfp                                                   */
/* return initial sfp                                                         */
/* WARNING : NO CARE in OVERLAP case !!!  strMove is safer!                   */
/******************************************************************************/
char     *strCpy(char *sfp,const char *sip)
{ char     *sp;

    sp=sfp;
    while( *sip ){ *sp = *sip; sip++;sp++;}
    *sp=0; 
    return sfp;
}

/******************************************************************************/
/*        strLen(sp)                                                          */
/* ~= strlen  string.h (ANSI)                                                 */
/* return the length of the Zero terminated string sp  (Zero not included)    */
/******************************************************************************/
size_t    strLen(const char  *sp)
{ size_t  nx=0;
  while( *sp ){nx++;sp++;}
  return nx;
}

#endif                                                      /** ifndef DEF_MAPSTR **/

/******************************************************************************/
/******************************************************************************/
