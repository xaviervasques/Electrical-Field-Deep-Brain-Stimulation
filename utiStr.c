/*   ../libmy/utiStr.c                                                        */
/*   Mennessier Gerard             940822                                     */
/*   Last revised M.G.             970701                                     */

#include  <stddef.h>                                       /* stddef.h (ANSI) */
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utistdStr.h"
#include  "utiStr.h"

/******************************************************************************/
/* see  GNU glibc source/sysdeps/generic/... for better algorithms            */
/******************************************************************************/

/******************************************************************************/
/* str  FAMILY : OTHERS                                                       */
/******************************************************************************/

/******************************************************************************/
/*        strChr2(sp,c1,c2,valuep)                                            */
/* search for the first occurence of byte c1 or c2                            */
/* equal to the value obtained by converting c1 or c2 to a char               */
/* in the string sp (Zero terminated)                                         */
/* return a pointer to this byte if found, or NULL                            */
/*        and value +1 if ==c1, -1 if ==c2, 0 if not found, at adress valuep  */
/******************************************************************************/
char     *strChr2(const char *sp,int c1,int c2,int *valuep)
{ char     *cp, cc1,cc2;

  *valuep=0;
  cp=(char *)sp; cc1=(char)c1; cc2=(char)c2;
  while( *cp )
  { if(*cp==cc1){ *valuep= 1; return cp;}
    if(*cp==cc2){ *valuep=-1; return cp;}
    cp++;
  }
  return NULL;
}
/******************************************************************************/
/*        strDiffChr(sp,c)                                                    */
/* search for the first occurence of byte Different from c = byte             */
/* equal to the value obtained by converting c to a char                      */
/* in the string sp (Zero terminated)                                         */
/* return a pointer to this byte if found, or NULL                            */
/******************************************************************************/
char     *strDiffChr(const char *sp,int c)
{ char     *cp, cc;

  cp=(char *)sp; cc=(char)c;
  while( *cp )
  { if(*cp != cc){ return cp;}
    cp++;
  }
  return NULL;
}
/******************************************************************************/
/*        strUcmp(s1p,s2p)                                                    */
/* ~= strcmp  string.h (ANSI)                                                 */
/* return <0 ,0, >0  if  s1p <, =, > s2p                                      */
/******************************************************************************/
int       strUcmp(unsigned char *c1p,unsigned char *c2p)
{ unsigned  char  c1,c2;

  while(1)
  { c1=*c1p; c2=*c2p;
 
    if(c1==c2){ if(c1==0){return 0;}  c1p++;c2p++; }
    else 
    { if(c1 > c2){ return 1; }
      return -1;
    }
  }
}
/******************************************************************************/
/*        strCpyP(sfp,sip)                                                    */
/* copy string sip (Zero terminated)  to sfp                                  */
/* sip and sfp strings should be NON-overlaping                               */
/* or  sip should be >= sfp                                                   */
/* return pointer on new ZERO terminating char                                */
/* WARNING : NO CARE in OVERLAP case !!!  strMove is safer!                   */
/* same function as GCC stpcpy                                                */
/******************************************************************************/
char     *strCpyP(char *sfp,const char *sip)
{
    while( *sip ){ *sfp = *sip;  sip++;  sfp++;}
    *sfp = 0; 
    return sfp;
}
/******************************************************************************/
/*        strMove(sfp,sip)                                                    */
/* copy string sip (Zero terminated)  to sfp                                  */
/* sip and sfp strings may be overlaping                                      */
/* return pointer on new ZERO terminating char                                */
/******************************************************************************/
char     *strMove(char *sfp,const char *sip)
{ size_t    n;
  char     *cfp, *cip;

  if(sip >= sfp){ return strCpyP(sfp,sip);}
  else
  { n = strLen(sip);  cfp = sfp + n;   cip = (char*)sip + n;
    while(cip >= sip){ *cfp = *cip;  cfp--; cip--;}
    return  (sfp+n);
  }
}
/******************************************************************************/
/*        strEq(s1p,s2p)                                                      */
/* return 1 (TRUE)  if s1p == s2p                                             */
/*     or 0 (FALSE) if s1p and s2p differs                                    */
/******************************************************************************/
int       strEq(char *s1p,char *s2p)
{ 
  while(1)
  { if(*s1p == *s2p)
    { if(*s1p == 0){return 1;}
      s1p++; s2p++;
    }
    else
    { return 0;}
  }
}
/******************************************************************************/
/*        strGetNextSpe(sp,c)                                                 */
/* search for the first occurence of the value obtained                       */
/* by converting c to a char                                                  */
/* return its address if found                                                */
/*     or NULL  if not found                                                  */
/******************************************************************************/
char     *strGetNextSpe(char  *sp,int c)
{ char      cc;

  cc=(char)c ;
  while( *sp ){ if(*sp== cc){return sp;}
                sp++;
              }
  return NULL;
}
/******************************************************************************/
/*        strGetNextSpes(sp,spsp)                                             */
/* search for the first occurence of one of the characters                    */
/* in the string spsp (Zero terminated)                                       */
/* return its address if found                                                */
/*     or NULL  if not found                                                  */
/******************************************************************************/
char     *strGetNextSpes(char  *sp,char *spsp)
{ char  type[256] =
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

  while(*spsp){ type[(int)*spsp]='\01' ; spsp++;}
  while( *sp ){ if( type[(int)*sp] ){return sp;}
                 sp++;
              }
  return NULL;
}
/******************************************************************************/
/*        strGetNextDiff(sp,c)                                                */
/* search for the first occurence DIFFERENT from the value obtained           */
/* by converting c to a char                                                  */
/* return its address if found                                                */
/*     or NULL  if not found                                                  */
/******************************************************************************/
char     *strGetNextDiff(char  *sp,int c)
{ char      cc;

  cc=(char)c ;
  while( *sp ){ if(*sp != cc){return sp;}
                sp++;
              }
  return NULL;
}

/******************************************************************************/
/*        strAtoL(sp)                                                         */
/* ~~= atol  or  strtol  stdlib.h                                             */
/* ascii string sp -> return long int  for base 10                            */
/* space (blanks) are ignored (not separators)                                */
/* check only expected characters and structure  in the entire string         */
/* expected structure : (sign) (ascii integer number)                         */
/* example:  "  -  6 4  3  "  will give the long int (-643)                   */
/******************************************************************************/
long      int       strAtoL(const char  *sp)
{ char    chsign=0;
  int     digit, zero='0' ;
  long    sum=0, base=10;
  char   *cp;
  char    form1[]=" character %c in %s is not a digit \n" ;

  cp=(char *)sp;
  while( *cp==' ' ){cp++;}         
  if(*cp=='-'){ chsign='\01'; cp++;}
  else        { if(*cp=='+'){ cp++;} }
  while( *cp )
  { if(*cp==' '){ cp++; continue;}
    digit=(int)(*cp)-zero;
    if(digit<0 || digit >9){ myErr1(-1,stderr,"strAtoL",form1, *cp,sp); }  
    sum=base*sum+digit;
    cp++;
  }
  if(chsign){sum= -sum;}
  return sum;
}
/******************************************************************************/
/*        strCpySkipChr(sfp,sip,c)                                            */
/* copy from sip to sfp , skipping if char equal c                            */
/* return pointer on the terminating 0  of out string                         */
/* WARNING : NO CARE in OVERLAP case !!!                                      */
/******************************************************************************/
char     *strCpySkipChr(char *sfp, const char *sip, int c)
{ char     *cfp, *cip ;

  cfp=sfp; cip=(char*)sip;
  while(*cip)
  { if(*cip != c){ *cfp= *cip ; cfp++; }
    cip++;
  }
  *cfp=0 ;
  return  cfp;
}
/******************************************************************************/
/*        strCpySkipStr(sfp,sip,spsp)                                         */
/* copy from mip to mfp , skipping if byte equal one char of string sp        */
/* return pointer on the terminating 0  of out string                         */
/* WARNING : NO CARE in OVERLAP case !!!                                      */
/******************************************************************************/
char     *strCpySkipStr(char *sfp, const char *sip, char *spsp)
{
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

  while(*spsp) {type[(int)*spsp]='\01' ; spsp++;}

  cfp=sfp; cip=(char*)sip;
  while(*cip)
  { if( type[(int)(*cip)] ){ cip++; continue; }
    *cfp = *cip;  cfp++;  cip++;
  }
  *cfp=0 ;
  return cfp;
}
/******************************************************************************/
/*        strNormalizeSpace(sfp,sip)                                          */
/* copy from sip to sfp ,                                                     */
/* suppress leading and ending Spaces (Blanks)                                */
/* keep only 1 Space when several ones in initial string.                     */
/* return pointer on the terminating 0  of out string                         */
/* WARNING : NO CARE in OVERLAP case !!!                                      */
/******************************************************************************/
char     *strNormalizeSpace(char *sfp, char *sip)
{ char     *cfp, *cip, c;
  char      space = ' ' ;
  short     previousIsSpace = 1;

  cip = sip;  cfp = sfp; 
  while( (c = *cip) )
  { if(c == space)
    { if(previousIsSpace == 0){ *cfp = space;  cfp++;  previousIsSpace = 1;}
    }
    else{ *cfp = c;  cfp++;  previousIsSpace = 0;}
    cip++;
  }
                        /** cfp-1 allowed only if at least 1 char i.e. cfp>=sfp+1 **/
  if(cfp != sfp){ if( *(cfp-1) == space){ cfp--; } }  
  *cfp = 0;
  return  cfp;
}

/******************************************************************************/
/*        strRank(spp,strnb,newsp)                                            */
/* given a table of strnb strings (in spp) in  any order,                     */
/* return the rank of a new string newsp                                      */
/* rank=-1 if new string does not appear  in the table                        */
/* rank= 0 if new string equal the first  of the table                        */
/* rank= 1 if new string equal the second of the table                        */
/* ...                                                                        */
/******************************************************************************/
int       strRank(char  *spp[],int strnb,char *newsp)
{ int  i;

  for(i=0; i<strnb; i++){ if( strEq(newsp,spp[i]) ){ return i; } }
  return (-1);
}
/******************************************************************************/
/*        strClas(spp,nx,newsp,indexp)                                        */
/* given a table of nx strings (in spp) ordered in ascending order,           */
/* class a new string newsp                                                   */
/* store value of index in indexp                                             */
/* and return cmp                                                             */
/*                                                                            */
/*          spp[*indexp] <= newsp <  spp[*indexp +1]                          */
/*      or                  newsp < spp[0]          (cmp= 1 *indexp= -1)      */
/*      or spp[nx -1]    <  newsp                   (cmp=-1 *indexp=nx -1)    */
/*                                                                            */
/* cmp = -1  if newsp < spp[0]                                                */
/* cmp =  0  if newsp = spp[*indexp]                                          */
/* cmp = +1  if newsp > spp[nx -1]                                            */
/* cmp  not significant in other cases  -> ONLY cmp=0 is always significant   */
/******************************************************************************/
int       strClas(char *spp[],int nx,char *newsp,int *indexp)
{ int  indx,indc,indw;
  int  cmp;

  indx=nx-1; indc=indx/2; indw=nx;
  while(indw)
  { cmp=strUcmp( (unsigned char *)newsp, (unsigned char *)spp[indc] ); 
    if(cmp==0){ *indexp=indc; return 0;}
    if(cmp<0)
    { if(indc==0){ *indexp= -1;    return (-1);}
      if(indw==1){ *indexp=indc-1; return (-1);}
      indw=(indw+1)/2; indc=( (indc-indw)<0 )? 0:indc-indw ; continue;
    }
    else
    { if(indc==indx){ *indexp=indx; return 1;}
      if(indw==1)   { *indexp=indc; return 1;}
      indw=(indw+1)/2; indc=( (indc+indw)>indx )? indx:indc+indw ; continue;
    }  
  }
  myErr1(-1,stderr,"classStrp : ERROR ",
                          "Should never reach this line STOP, unless void list? \n");
  return  (-99);
} 
/******************************************************************************/
/*        strHalfSum(char *strfp,int *ip,char *str2p,char *str1p)             */
/* given strings str1p,str2p, fill string strfp by the  'half-sum'            */
/* ip points on a set of working memory.                                      */
/* allocation of n = max(length of str1p,str2p)  must have been done for ip   */
/* allocation of n+1  must have been done for strfp (+1 for ending zero)      */
/******************************************************************************/
void      strHalfSum(unsigned char*strfp,int *sp,unsigned char*c2p,unsigned char*c1p)
{ ptrdiff_t  ix;
  size_t     i;
  int       *lp, *lastp, odd;

  lp=sp;
  odd=0;
  while( *c1p|| *c2p)
  { if(odd){ *lp=256;}
    else   { *lp=0;}
    if(*c1p){ *lp += *c1p; c1p++; }
    if(*c2p){ *lp += *c2p; c2p++; }
    odd= (*lp) %2;
    *lp = (*lp) >> 1 ;  lp++;
  }
  lastp=lp-1;  ix=lp-sp;
  
  lp=lastp;
  for(i=ix; i>0; i--)
  { if( (*lp) >> 8 ){ *lp -=256; *(lp-1) +=1; }
    lp--;
  }
  for(i=ix; i>0; i--,sp++,strfp++) { *strfp = (unsigned char)(*sp) ; }
  *strfp=0;                                                /** ending string zero **/
  return;
}
/******************************************************************************/
/******************************************************************************/
