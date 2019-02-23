/*   ../libmy/utiMem.c                                                        */
/*   Mennessier Gerard             940822                                     */
/*   Last revised M.G.             990628                                     */

#include  <stddef.h>                                       /* stddef.h (ANSI) */
#include  "utistdErr.h"
#include  "utiMem.h"

/******************************************************************************/
/* see  GNU glibc source/sysdeps/generic/... for better algorithms            */
/******************************************************************************/

/******************************************************************************/
/* mem  FAMILY : OTHERS                                                       */
/******************************************************************************/

/******************************************************************************/
/*        memChrBack(mp,c,n)                                                  */
/* ~= memchr  string.h (ANSI)                                                 */
/* search for the first occurence of a byte (in the n BACKWARD bytes mp)      */
/* equal to the value obtained by converting c to a char                      */
/* return a pointer to this byte if found, or NULL                            */
/******************************************************************************/
void     *memChrBack(const void *mp,int c,size_t n)
{ char     *cp, *cps, cc;

  cp=(char *)mp; cps=cp- (ptrdiff_t)n; cc=(char)c;
  for( ; cp>cps; cp--)
  { if(*cp == cc){return (void *)cp;}
  }
  return NULL;
}

/******************************************************************************/
/*        memChr2(sp,c1,c2,valuep)                                            */
/* search for the first occurence of byte c1 or c2                            */
/* equal to the value obtained by converting c1 or c2 to a char               */
/* in the string sp (Zero terminated)                                         */
/* return a pointer to this byte if found, or NULL                            */
/*        and value +1 if ==c1, -1 if ==c2, 0 if not found, at adress valuep  */
/******************************************************************************/
void     *memChr2(const void *mp,int c1,int c2,int *valuep,size_t n)
{ char     *cp, *cps,  cc1,cc2;

  *valuep=0;
  cp=(char *)mp; cps=cp+n; cc1=(char)c1; cc2=(char)c2;
  for( ; cp<cps; cp++)
  { if(*cp==cc1){ *valuep= 1; return (void *)cp;}
    if(*cp==cc2){ *valuep=-1; return (void *)cp;}
  }
  return NULL;
}

/******************************************************************************/
/*        memDiffChr(mp,c,n)                                                  */
/* search for the first occurence (in the n bytes mp)                         */
/* of byte Different from c = byte                                            */
/* equal to the value obtained by converting c to a char                      */
/* return a pointer to this byte if found, or NULL                            */
/******************************************************************************/
void     *memDiffChr(const void *mp,int c,size_t n)
{ char     *cp, *cps, cc;

  cp=(char *)mp; cps=cp+n; cc=(char)c;
  for( ; cp<cps; cp++)
  { if(*cp != cc){return (void *)cp;}
  }
  return NULL;
}

/******************************************************************************/
/*        memCmp2(m1p,n1,m2p,n2)                                              */
/* compare the 2 sets of bytes  m1p (length n1) and m2p (length n2)           */
/* as beeing unsigned char                                                    */
/* return -1,0,1  if  m1p <, =, > m2p                                         */
/******************************************************************************/
int       memCmp2(const void *m1p, size_t n1, const void *m2p, size_t n2)
{ size_t    i;
  int       cd;
  char     *c1p, *c2p ;

  c1p=(char *)m1p; c2p=(char *)m2p;
  if(n1<=n2)
  { for(i=0; i<n1; i++)
    { cd= (unsigned char)*c1p - (unsigned char)*c2p ;
      if(cd==0){ c1p++; c2p++; }
      else
      { if(cd>0){return 1;}
        return -1;
      }
    }
    if(n1==n2){return 0;}
    else {return -1;}
  }
  else
  { for(i=0; i<n2; i++)
    { cd= (unsigned char)*c1p - (unsigned char)*c2p ;
      if(cd==0){ c1p++;c2p++; }
      else
      { if(cd>0){return 1;}
        return -1;
      }
    }
    return 1;
  }
}

/******************************************************************************/
/*        memCpyP(mfp,mip,n)                                                  */
/* copy n bytes from mip to mfp                                               */
/* return pointer on next to last copied byte                                 */
/* WARNING : NO CARE in OVERLAP case !!!  memMove is safer!                   */
/******************************************************************************/
void     *memCpyP(void *mfp, const void  *mip, size_t n)
{ size_t    i;
  char     *cfp, *cip ;

  cfp=(char *)mfp; cip=(char *)mip;
  for( i=0; i<n ; i++)
  { *cfp= *cip ; cfp++; cip++; }
  return cfp;
}

/******************************************************************************/
/*        memEq(m1p,n1,m2p,n2)                                                */
/* return 1 (TRUE)  if n1==n2 and                                             */
/*                     all bytes of m1p == corresponding bytes of m2p         */
/*     or 0 (FALSE) if m1p and m2p differ                                     */
/******************************************************************************/
int       memEq(const void *m1p,size_t n1,const void *m2p,size_t n2)
{ size_t    i;
  char     *c1p, *c2p ;

  if(n1!=n2){return 0;}

  c1p=(char *)m1p; c2p=(char *)m2p;
  for(i=0;i<n1;i++)
  { if(*c1p==*c2p){ c1p++; c2p++; }
    else { return 0;}
  }
  return 1;
}
/******************************************************************************/
/*        memGetNextSpe(mp,n,c)                                               */
/* search for the first occurence of the value obtained                       */
/* by converting c to a char                                                  */
/* return i, its relative address if found                                    */
/*     or -1  if not found                                                    */
/* ~~=memChr(mp,c,n) but return index, not absolute address                   */
/******************************************************************************/
int       memGetNextSpe(const void *mp,size_t n,int c)
{ int       i;
  char     *cp, cc;

  cp=(char *)mp; cc=c;
  for(i=0; i<n; i++){ if( *cp==cc ){return i;}
                       cp++;
                    }
  return (-1);
}
/******************************************************************************/
/*        memGetNextSpes(mp,n,sp)                                             */
/* search for the first occurence of one of the characters                    */
/* in the string sp (Zero terminated)                                         */
/* return i, its relative address if found                                    */
/*     or -1  if not found                                                    */
/******************************************************************************/
int       memGetNextSpes(const void *mp,size_t n,char *sp)
{ char      type[256] = 
            {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,

             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
            };
  size_t    i;
  char      *cp;

  while(*sp) {type[(int)*sp]='\01' ; sp++;}
  cp=(char *)mp;
  for(i=0; i<n; i++){ if( type[(int)*cp] ){return i;}
                       cp++; 
                    }
  return (-1);
}
/******************************************************************************/
/*        memAtoL(mp,n)                                                       */
/* ~~~= atol  or  strtol  stdlib.h                                            */
/* ascii in n bytes mp -> return long int  for base 10                        */
/* space (blanks) are ignored (NOT separators)                                */
/* check only expected characters and structure  in the entire string         */
/* expected structure : (sign) (ascii integer number)                         */
/* example:  "  -  6 4  3  "  will give the long int (-643)                   */
/******************************************************************************/
long      memAtoL(const void *mp,size_t n)
{ short     boolchsign = 0;
  int       digit, zero = '0' ;
  long      sum = 0, base = 10;
  char     *cp;
  char      form1[] = " character %c in %s is not a digit \n";

  cp = (char *)mp;
  while(n  &&  *cp == ' '){ cp++;  n--;}
  if(n)
  { if(*cp == '-'){ boolchsign = 1;  cp++;  n--;}
    else          { if(*cp == '+') { cp++;  n--;} }
  }
  while(n)
  {
    if(*cp == ' '){ cp++;  n--;  continue;}
    digit = (int)(*cp) - zero;
    if(digit < 0  ||  digit > 9){ myErr1(-1,stderr,"memAtoL",form1, *cp,mp); } 
    sum = base * sum + digit;
    cp++;  n--;
  }
  if(boolchsign) sum = -sum;
  return sum;
}
/******************************************************************************/
/*        memAtoLpart(sum,mp,n)                                               */
/*                                                                            */
/* expected structure : [spaces][sign](ascii integer number WITHOUT space)    */
/* return : pointer to next to last digit                                     */
/*          value into *sump                                                  */
/******************************************************************************/
char     *memAtoLpart(long *sump, const void *mp, size_t n)
{ short     boolchsign = 0;
  int       digit, zero = '0' ;
  long      sum = 0, base = 10;
  char     *cp;

  cp = (char *)mp;
  while(n  &&  *cp == ' '){ cp++;  n--;}                  /** skip leading spaces **/
  if(n)
  { if(*cp == '-'){ boolchsign = 1;  cp++;  n--;}
    else          { if(*cp == '+') { cp++;  n--;} }
  }
  while(n)
  { digit = (int)(*cp++) - zero;
    if(digit < 0  ||  digit > 9) break;        /** first NON digit => End of data **/
    sum = base * sum + digit;
    n--;
  }
  if(boolchsign) sum = -sum;
  *sump = sum;
  return  cp;
}
/******************************************************************************/
/*        memCpySkipChr(mfp,mip,n,c)                                          */
/* copy from mip to mfp , skipping if byte equal c                            */
/* return the new number of bytes                                             */
/* WARNING : NO CARE in OVERLAP case !!!                                      */
/******************************************************************************/
size_t    memCpySkipChr(void *mfp, const void  *mip, size_t n, int c)
{ size_t    i;
  char     *cfp, *cip, cc;

  cc = (char)c;  cfp=(char *)mfp;  cip=(char *)mip;
  for(i=0; i<n ; i++)
  { if(*cip != cc){ *cfp= *cip;  cfp++; }
    cip++;
  }
  return ( cfp- (char *)mfp );
}
/******************************************************************************/
/*        memCpySkipStr(mfp,mip,n,sp)                                         */
/* copy from mip to mfp , skipping if byte equal one char of string sp        */
/* return the new number of bytes                                             */
/* WARNING : NO CARE in OVERLAP case !!!                                      */
/******************************************************************************/
size_t    memCpySkipStr(void *mfp, const void  *mip, size_t n, char *sp)
{ size_t    i;
  char     *cfp, *cip ;
  char      type[256] = 
            {0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,

             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0,
             0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0
            };

  while(*sp) {type[(int)*sp]='\01' ; sp++;}

  cfp=(char *)mfp; cip=(char *)mip;
  for( i=0; i<n ; i++)
  { if( type[(int)(*cip)] ){ cip++; continue; }
    *cfp= *cip ; cfp++; cip++;
  }
  return ( cfp- (char *)mfp );
}

/******************************************************************************/
/*        memNormalizeSpace(mfp,mip,n)                                        */
/* copy from sip to sfp ,                                                     */
/* suppress leading and ending Spaces (Blanks)                                */
/* keep only 1 Space when several ones in initial string.                     */
/* return the new number of bytes                                             */
/* WARNING : NO CARE in OVERLAP case !!!                                      */
/******************************************************************************/
size_t    memNormalizeSpace(void *mfp, const void  *mip, size_t n)
{ char     *cfp, *cip ,*cps, *cfpi;
  char      space=' ' ;
  short     previousIsSpace;

  cfpi=cfp=(char *)mfp; cip=(char *)mip; cps=cip+n;
  previousIsSpace=1;
  for( ; cip<cps; cip++)
  { if(*cip == space)
    { if(previousIsSpace==0){ *cfp= space; cfp++; previousIsSpace=1; }
    }
    else{ *cfp= *cip ; cfp++; previousIsSpace=0; }
  }
                        /** cfp-1 allowed only if at least 1 char i.e. cfp>=mfp+1 **/
  if(cfp!=cfpi){ if( *(cfp-1) == space){ cfp--; } } 

  return  (size_t)(cfp - cfpi);
}
/******************************************************************************/
/*        memMem(m1p,n1,m2p,n2)                                               */
/* search the n2 bytes m2p  inside the n1 bytes m1p                           */
/* return : if match, pointer to beginning of match in m1p                    */
/*          else NULL                                                         */
/******************************************************************************/
void     *memMem(const void *m1p, size_t n1, const void *m2p, size_t n2)
{ size_t    i;
  int       match = 0;
  char     *c1p, *c2p, *c1fp, *cp;

  if(n1 < n2) return NULL;
  c1p = (char *)m1p;  c1fp = (c1p + n1) - (ptrdiff_t)n2;
  while(c1p <= c1fp)
  { match = 1;
    for(i = 0, cp = c1p, c2p = (char *)m2p;  i < n2;  i++)
    { if(*cp++ != *c2p++){ match = 0;  break;}
    }
    if(match){ return  (void*)c1p;}
    else     { c1p++;}
  }
  return NULL;
}
/******************************************************************************/
/******************************************************************************/
