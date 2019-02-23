/*  ../libmy/utiClassRangeStrp.c                                              */
/*  Mennessier Gerard                   940104                                */
/*  Last Revised : M.G.                 970524                                */

#include  <stddef.h>
#include  "utistdErr.h"
#include  "utiAlloc.h"
#include  "utistdStr.h"
#include  "utiStr.h"

#define   MYSWAP(in,ix)  ( swap=indp[in],  indp[in]=indp[ix], indp[ix]=swap )
#define   MY_MAXSTRSIZE_     127
/******************************************************************************/
/** REMARQUE: utilise 1/2 somme et NON moyenne ==> no "malheureux problems"  **/
/** --> version 12                                                           **/
/** longueur maximale des String : MY_MAXSTRSIZE_                            **/
/******************************************************************************/

/******************************************************************************/
/*            ranindStrp(&x,&ind,nx)                                          */
/*                                                                            */
/*  x array de nx string pointers {x[0],..x[n-1]}                             */
/*  calcule les  int ind[i] tels que                                          */
/*  x[ind[0]] <= x[ind[1]] <= x[ind[2]] ....<=x[ind[n-1]]                     */
/******************************************************************************/
void      ranindStrp(unsigned char **xp,int *indp,size_t nz)
{
  int      *nxp;                                    /** pointeur sur array des nx **/
  short     initnxsize = 50;                     /** Taille 50 suffit d ordinaire **/
  size_t    nxsize;
  short     kx;                 /** indices des nx  =nombre d intervalle en cours **/
  int       n,i,in,ix,j1,j2,jmin,jmax,jnx;                    /** indices des ind **/
  int       jnmin=0, jnmax=0, j1ex=0, j2ex=0;                 /** indices des ind **/
  int       swap;                                         /** pour swaper 2 index **/
  short     j1l,j2l;              /** logical yes=1 si existe, no=0 si non existe **/

                   /**** WARNING : NEED length of all strings <= MY_MAXSTRSIZE_ ****/

  unsigned  char      xm[MY_MAXSTRSIZE_ +1] ;
  int       ip[MY_MAXSTRSIZE_];

  n=nz;  for(i=0; i<n; i++){ indp[i]=i;}
  if(n < 2 ){ return ;}

  nxsize = initnxsize;
  nxp= intAlloc(nxsize, "ranindStrp" );
  *nxp=n-1;  *(nxp+1)=0;  kx=1;

  while(kx >= 1)
  { in = *(nxp+kx);  ix = *(nxp+kx-1);
    if(in >= ix-1)                  
    {                                                     /** cas longueur 2 ou 1 **/
      if(in == ix-1)                                           /** cas longueur 2 **/
      { if( strUcmp(*(xp+indp[in]), *(xp+indp[ix]) )  > 0){ MYSWAP(in,ix); }
      }
      kx=kx-1;  *(nxp+kx)=  *(nxp+kx)+1 ;
    }
    else                              /** cas generique intervalle de longueur >2 **/
    { strHalfSum(xm,ip, *(xp+indp[in]), *(xp+indp[ix]) ) ;
                                              /** moyenne des extremitees calculee**/

/**recherche d un premier element >= moyenne xm; 
   doit toujours exister sauf arrondis malheureux lors de calcul de xm;
   si n existe pas alors choisir xm=valeur du premier element et declarer j1ex=in!**/
      j1l=0;   
      for(j1=in; j1<=ix; j1++)
      { if( strUcmp( *(xp+indp[j1]), xm) >= 0){ j1l=1; j1ex=j1; jnmin=j1+1; break;} 
      }
      if(j1l==0) 
      { j1l=1; j1ex=in; jnmin=in+1; strCpy((char*)xm,(char*) *(xp+indp[in]) ); }
                                                        /** donc j1l=1 toujours ! **/
                      /* recherche d un premier element en reculant < moyenne xm; **/
      j2l=0;
      for(j2=ix; j2>=jnmin; j2--)
      { if( strUcmp(*(xp+indp[j2]), xm) <= 0){ j2l=1; j2ex=j2; jnmax=j2-1; break;}
      }
      if(j2l){ jnx=j1ex; MYSWAP(j1ex,j2ex); }
      else
      { jnmax=jnmin-1;
        if(j1ex==in){ jnx=in;}
        else        { jnx=j1ex-1;}
      }
      /** fin recherche premier element et de premiere iteration de swap eventuel **/
                                                 /** suite des iterations de swap **/
      jmin=jnmin; jmax=jnmax;           
      while(jmin <= jmax)
      {    /** debut des 2 boucles de recherche de swap dans intervalle jmin,jmax **/
        j1l=0;   j2l=0;           /** mise a 0=false donc n'existe pas par defaut **/
        jnmin=jmax+1;       /** valeur par defaut  si boucle croissante tj fausse **/
        for(j1=jmin; j1<=jmax; j1++)
        { if( strUcmp( *(xp+indp[j1]), xm) >= 0){ j1l=1; j1ex=j1; jnmin=j1+1; break;}
        }      
        jnmax=jmin-1;     /** valeur par defaut  si boucle decroissante tj fausse **/
        for(j2=jmax; j2>=jmin; j2--)
        { if( strUcmp( *(xp+indp[j2]), xm) <= 0){ j2l=1; j2ex=j2; jnmax=j2-1; break;}
        }
                                                /** debut resultats des 2 boucles **/
        if(j1l)
        {                                               /** debut cas j1ex existe **/
          if(j2l)
          {                        /** debut cas j1ex et j2ex existent tous les 2 **/
            if(j1ex<j2ex){ jnx=j1ex;  MYSWAP(j1ex,j2ex); }
            else { jnx=j2ex; break;}   /** cas ou deja classe dans cet intervalle **/
          }                        /** fin   cas j1ex et j2ex existent tous les 2 **/
          else { jnx=j1ex-1; break;}       /** cas j1ex existe, j2ex n'existe pas **/
        }                                                 /** fin cas j1ex existe **/
        else
        {                                         /** debut cas j1ex n'existe pas **/
          if(j2l){ jnx=j2ex; break;}       /** cas j1ex n'existe pas, j2ex existe **/
          else                                /** cas j1ex et j2ex n'existent pas **/
          { myErr0(-1,stderr,"ranindStrp erreur : ne doit pas arriver \n" );
          }
        }                                          /** fin  cas j1ex n'existe pas **/
                                                /** fin   resultats des 2 boucles **/
        jmax=jnmax; jmin=jnmin;        
      }     /** fin des iterations de recherche de swap dans intervalle jmin,jmax **/
                                                    /** i.e. fin while(jmin<=jmax **/
                 /** creation nouvel intervalle sauf si de longueur 1 car trivial **/
      if(in >= jnx){ *(nxp+kx)=*(nxp+kx) +1 ;}
      else
      { kx++ ;
        if(kx >= nxsize)
        { nxp=intChkRealloc(nxp,&nxsize,(size_t)kx+1,10,"ranindStrp");
        }
        *(nxp+kx)=*(nxp+kx-1);  *(nxp+kx-1)=jnx;
      }
    }                                                       /** fin cas generique **/
  }                                                       /** fin  while(kx >= 1) **/
  free(nxp); return;
}
/******************************************************************************/
/*        classStrp(indexp,newsp,spp,nx)                                      */
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
/* cmp = +1  if newsp > spp[strnb -1]                                         */
/* cmp  not significant in other cases  -> ONLY cmp=0 is always significant   */
/******************************************************************************/
int       classStrp(int *indexp,char *newsp,char  *spp[],size_t nx)
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
/******************************************************************************/
