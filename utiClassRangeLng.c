/*  ../libmy/utiClassRangeLng.c                                               */
/*  Mennessier Gerard                   940104                                */
/*  Last Revised              M.G.      970421                                */

#include  <stddef.h>

#include  "utiClassRangeLng.h"
#include  "utistdErr.h"
#include  "utiAlloc.h"

#define   MYSWAP(in,ix)  ( swap=indp[in],  indp[in]=indp[ix], indp[ix]=swap )
/******************************************************************************/
/*            ranindLng(&x,&ind,nz)                                           */
/* entier -> gap -> no "malheureux pb"                                        */
/* version 11                                                                 */
/*                                                                            */
/*  x array de n long {x[0],..x[n-1]}                                         */
/*  calcule les int ind[i] tels que                                           */
/*  x[ind[0]] <= x[ind[1]] <= x[ind[2]] ....<=x[ind[n-1]]                     */
/******************************************************************************/
void      ranindLng(long *xp,int *indp,size_t nz)
{ int   *nxp;                                       /** pointeur sur array des nx **/
  size_t nxsize = 50;                            /** Taille 50 suffit d ordinaire **/
  short  kx;                    /** indices des nx  =nombre d intervalle en cours **/
  int    n,i,in,ix,j1,j2,jmin,jmax;                           /** indices des ind **/
  int    jnmin=0,jnmax=0,jnx,j1ex=0,j2ex=0;                   /** indices des ind **/
  int    swap;                                            /** pour swaper 2 index **/
  short  j1l,j2l;                 /** logical yes=1 si existe, no=0 si non existe **/
  double xmd;
  long   xm;

  n=nz;  for(i=0; i<n; i++){ indp[i]=i;}
  if(n < 2 ){ return ;}

  nxp = intAlloc(nxsize, "ranindLng");
  *nxp=n-1;  *(nxp+1)=0;  kx=1;

  while(kx >= 1)
  { in = *(nxp+kx);  ix = *(nxp+kx-1);
    if(in >= ix-1)                  
    {                                                     /** cas longueur 2 ou 1 **/
      if(in == ix-1){ if(*(xp+indp[in]) > *(xp+indp[ix])){ MYSWAP(in,ix);} }
      kx=kx-1;  *(nxp+kx) = *(nxp+kx)+1;
    }
    else                              /** cas generique intervalle de longueur >2 **/
    { xmd=0. ;
      for(i=in; i<=ix; i++) { xmd = xmd + *(xp+indp[i]) ;}
      xmd=xmd/( (double)(ix-in+1) );
      xm= (xmd >= 0)?  (xmd+1.e0) : (xmd-1.e0) ;             /** moyenne calculee **/

/**recherche d un premier element >= moyenne xm; 
   doit toujours exister sauf arrondis malheureux lors de calcul de xm;
   si n existe pas alors choisir xm=valeur du premier element et declarer j1ex=in!**/
      j1l=0;   
      for(j1=in; j1<=ix; j1++)
      { if( *(xp+indp[j1]) >= xm){ j1l=1; j1ex=j1; jnmin=j1+1; break;}
      }
      if(j1l==0){ j1l=1; j1ex=in; jnmin=in+1; xm= *(xp+ *(indp+in));}        
                                                        /** donc j1l=1 toujours ! **/
                      /* recherche d un premier element en reculant < moyenne xm; **/
      j2l=0;
      for(j2=ix; j2>=jnmin; j2--)
      { if( *(xp+indp[j2]) < xm){ j2l=1; j2ex=j2; jnmax=j2-1; break;}
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
        { if( *(xp+indp[j1]) >= xm){ j1l=1; j1ex=j1; jnmin=j1+1; break;}
        }      
        jnmax=jmin-1;     /** valeur par defaut  si boucle decroissante tj fausse **/
        for(j2=jmax; j2>=jmin; j2--)
        { if( *(xp+indp[j2]) <= xm){ j2l=1; j2ex=j2; jnmax=j2-1; break;}
        }
                                                /** debut resultats des 2 boucles **/
        if(j1l)
        {                                               /** debut cas j1ex existe **/
          if(j2l)
          {                        /** debut cas j1ex et j2ex existent tous les 2 **/
            if(j1ex<j2ex){ jnx=j1ex;  MYSWAP(j1ex,j2ex);}
            else { jnx=j2ex; break;}   /** cas ou deja classe dans cet intervalle **/
          }                        /** fin   cas j1ex et j2ex existent tous les 2 **/
          else { jnx=j1ex-1; break;}       /** cas j1ex existe, j2ex n'existe pas **/
        }                                                 /** fin cas j1ex existe **/
        else
        {                                         /** debut cas j1ex n'existe pas **/
          if(j2l){jnx=j2ex; break;}        /** cas j1ex n'existe pas, j2ex existe **/
          else                                /** cas j1ex et j2ex n'existent pas **/
          { myErr0(-1,stderr,"ranindLng erreur : ne doit pas arriver \n" );
          }
        }                                          /** fin  cas j1ex n'existe pas **/
                                                /** fin   resultats des 2 boucles **/
        jmax=jnmax; jmin=jnmin;
      }     /** fin des iterations de recherche de swap dans intervalle jmin,jmax **/
                                                    /** i.e. fin while(jmin<=jmax **/
                 /** creation nouvel intervalle sauf si de longueur 1 car trivial **/
/*      fprintf(stderr," kx=%6d ,in=%6d  ,ix=%6d ,jnx=%6d \n",kx,in,ix,jnx);  */
      if(in >= jnx){ *(nxp+kx)=*(nxp+kx) +1;}
      else
      { kx++;
        if(kx >= nxsize)
        { nxp=intChkRealloc(nxp,&nxsize,(size_t)kx+1,10,"ranindLng");
        }
        *(nxp+kx)=*(nxp+kx-1);  *(nxp+kx-1)=jnx;
      } 
    }                                                       /** fin cas generique **/
  }                                                       /** fin  while(kx >= 1) **/
  free(nxp); return;
}
/******************************************************************************/
/*        classLng(indexp,x,xp,nx)                                            */
/* given a table of nx long (in xp) ordered in ascending order,               */
/* class a new long x.                                                        */
/* store value of index in indexp                                             */
/* and return cmp                                                             */
/*                                                                            */
/*              xp[*indexp] <= x < xp[*indexp +1]                             */
/*           or                x < xp[0]          (cmp= 1 *indexp= -1)        */
/*           or xp[nx -1]   <  x                  (cmp=-1 *indexp=nx -1)      */
/*                                                                            */
/* cmp = -1  if  x < xp[0]                                                    */
/* cmp =  0  if  x = xp[*indexp]                                              */
/* cmp = +1  if  x > xp[nx -1]                                                */
/* cmp  not significant in other cases -> ONLY cmp=0 is always significant    */
/******************************************************************************/
int       classLng(int *indexp,long x,long *xp,size_t nx)
{ int       indx,indc,indw;
  long      diff;

  indx=nx-1; indc=indx/2; indw=nx;
  while(indw)
  { diff= x - xp[indc];
    if(diff==0){ *indexp=indc; return 0;}
    if(diff<0)
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
  myErr1(-1,stderr,"classLng : ERROR ", 
                           "Should never reach this line STOP, unless void list?\n");
  return (-99) ;
}
/******************************************************************************/
/******************************************************************************/
