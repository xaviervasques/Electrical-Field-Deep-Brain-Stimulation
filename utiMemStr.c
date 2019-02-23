/*   ../libmy/utiMemStr.c                                                     */
/*   Mennessier Gerard             960123                                     */
/*   Last revised M.G.             960506                                     */

#include  <stddef.h>                                       /* stddef.h (ANSI) */
#include  "utiMemStr.h"
#include  "utistdErr.h"

/******************************************************************************/
/* see  GNU glibc source/sysdeps/generic/... for better algorithms            */
/******************************************************************************/

/******************************************************************************/
/*        memStrUcmp(mp,ix,sp)                                                */
/* compare a mem mp (length ix)  to a string sp                               */
/* return <0 ,0, >0  if  mp <, =, > sp                                        */
/******************************************************************************/
int       memStrUcmp(void *mp, size_t ix, unsigned char *sp)
{ unsigned  char   *c1p, *c2p ;
  unsigned  char    c1,c2;
  size_t    i;
  
  c1p=(unsigned char *)mp;  c2p=(unsigned char *)sp;
  for(i=ix; i>0; i--, c1p++,c2p++)
  { c1=*c1p; c2=*c2p;
    if(c1==c2){ continue; }
    else 
    { if(c1 > c2){ return 1; }
      return -1;
    }
  }
  if(*c2p){ return -1;}
  else    { return  0;}
}

/******************************************************************************/
/*        memStrClas(spp,strnb, mp,ix, indexp)                                */
/* given a table of strnb strings (in spp) ordered in ascending order,        */
/* class the mem mp, length ix                                                */
/* store value of index in indexp                                             */
/* and return cmp                                                             */
/*                                                                            */
/*          spp[*indexp] <= mp <  spp[*indexp +1]                             */
/*      or                  mp < spp[0]          (cmp= 1 *indexp= -1)         */
/*      or spp[strnb -1] <  mp                   (cmp=-1 *indexp=strnb -1)    */
/*                                                                            */
/* cmp = -1  if mp < spp[0]                                                   */
/* cmp =  0  if mp = spp[*indexp]                                             */
/* cmp = +1  if mp > spp[strnb -1]                                            */
/* cmp  not significant in other cases  -> ONLY cmp=0 is always significant   */
/******************************************************************************/
int       memStrClas(char  *spp[],int strnb,void *mp,size_t ix,int *indexp)
{ int       indx,indc,indw;
  int       cmp;
  unsigned  char   *cp;

  cp = (unsigned char *)mp;
  indx=strnb-1; indc=indx/2; indw=strnb;
  while(indw)
  { cmp= memStrUcmp( cp,ix, (unsigned char *)spp[indc] ); 
    if(cmp==0){ *indexp=indc; return 0; }
    if(cmp<0)
    { if(indc==0){ *indexp= -1;    return (-1) ;}
      if(indw==1){ *indexp=indc-1; return (-1) ;}
      indw=(indw+1)/2; indc=( (indc-indw)<0 )? 0:indc-indw ; continue;
    }
    else
    { if(indc==indx){ *indexp=indx; return 1 ;}
      if(indw==1)   { *indexp=indc; return 1 ;}
      indw=(indw+1)/2; indc=( (indc+indw)>indx )? indx:indc+indw ; continue;
    }  
  }
  myErr1(-1,stderr,"strClas ", 
                         "Should never reach this line STOP, unless void list? \n") ;
  return (-99);
} 
 
/******************************************************************************/
/*        memStrEq(mp,ix,sp)                                                  */
/* return 1 (TRUE)  if  mp == sp  (with length sp = ix)                       */
/*     or 0 (FALSE) if  mp and sp differ                                      */
/******************************************************************************/
int       memStrEq(void *mp,size_t ix,char *sp)
{ char     *cp;
  int       i;

  cp=(char *)mp;
  for(i=ix; i>0; i--, cp++,sp++)
  { if(*cp != *sp){ return 0;}
  }
  if( *sp){ return 0;}
  else    { return 1;}
}

/******************************************************************************/
/*        memStrRank(spp,strnb,mp,ix)                                         */
/* given a table of strnb strings (in spp) in  any order,                     */
/* return the rank of mp (length ix)                                          */
/* rank=-1 if new string does not appear  in the table                        */
/* rank= 0 if new string equal the first  of the table                        */
/* rank= 1 if new string equal the second of the table                        */
/* ...                                                                        */
/******************************************************************************/
int       memStrRank(char  *spp[],int strnb,void *mp,size_t ix)
{ int  i;

  for(i=0; i<strnb; i++){ if( memStrEq(mp,ix,spp[i]) ){ return i; } }
  return (-1);
}

/******************************************************************************/
/******************************************************************************/
