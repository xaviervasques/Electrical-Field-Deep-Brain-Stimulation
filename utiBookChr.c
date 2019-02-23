/*  ../libmy/utiBookChr.c                                                     */
/*  Mennessier Gerard                   940822                                */
/*  Last revised : M.G.                 000801                                */

#include  "utiBookChr.h"

#include  "utistdErr.h"
#include  "utistdMem.h"
#include  "utistdStr.h"
#include  "utiAlloc.h"
#include  "utiVecChr.h"
#include  "utiVecPtr.h"

/******************************************************************************/
/*        chrBBookAlloc(bkp,nz)                                               */
/******************************************************************************/
void      chrBBookAlloc(chrBook *bkp,size_t  nz)
{  
  bkp->bp = chrAlloc(nz, "chrBBookAlloc");  bkp->bz = nz;  bkp->bx = 0;
  if(nz > 0){ *(bkp->bp) = 0;}
  return ;
}
/******************************************************************************/
/*        chrBBookRealloc(bkp,neednz,incrnz)                                  */
/******************************************************************************/
void      chrBBookRealloc(chrBook *bkp,size_t neednz,size_t incrnz)
{  
  bkp->bp = chrChkRealloc(bkp->bp, &(bkp->bz),neednz,incrnz,"chrBBookRealloc");
  return ;
}
/******************************************************************************/
/*        chrIBookAlloc(bkp,nz)                                               */
/******************************************************************************/
void      chrIBookAlloc(chrBook *bkp,size_t  nz)
{ 
  bkp->ip = ptrdfAlloc(nz, "chrIBookAlloc");  bkp->iz = nz;  bkp->ix = 0;
  if(nz > 0){ *(bkp->ip) = 0;}
  return ;
}
/******************************************************************************/
/*        chrIBookRealloc(bkp,needz,incrnz)                                   */
/******************************************************************************/
void      chrIBookRealloc(chrBook *bkp,size_t neednz,size_t incrnz)
{ 
  bkp->ip = ptrdfChkRealloc(bkp->ip, &(bkp->iz),neednz,incrnz,"chrIBookRealloc");
  return ;
}
/******************************************************************************/
/*        chrIBookFree(bkp)                                                   */
/******************************************************************************/
void      chrIBookFree(chrBook *bkp)
{ 
  free(bkp->ip); bkp->ip=NULL;  bkp->iz=0; bkp->ix=0; 
  return ;
}
/******************************************************************************/
/*        chrBBookFree(bkp)                                                   */
/******************************************************************************/
void      chrBBookFree(chrBook *bkp)
{ 
  free(bkp->bp); bkp->bp=NULL;  bkp->bz=0; bkp->bx=0; 
  return ;
}
/******************************************************************************/
/*        chrBookFree(bkp)                                                    */
/******************************************************************************/
void      chrBookFree(chrBook *bkp)
{ if(bkp == NULL){ return;}
  if( bkp->bp != NULL){ chrBBookFree(bkp);}
  if( bkp->ip != NULL){ chrIBookFree(bkp);}
  free(bkp);  bkp = NULL;
  return ;
}
/******************************************************************************/
/*        chrBookCPrint(bufp,bkp)                                              */
/******************************************************************************/
void      chrBookCPrint(FILE  *bufp, chrBook *bkp)
{ char     *cp;
  ptrdiff_t     *ip;
  size_t    ibc, iic, len;
  char      form1[]="book.bp=%p, book.bz=%d, book.bx=%d \n" ;
  char      form2[]="book.ip=%p, book.iz=%d, book.ix=%d \n" ;
  char      form3[]="char set number i=%d, ibeg=%d, iend=%d : \n" ;

  fPrintF(bufp,form1, bkp->bp,bkp->bz,bkp->bx);
  fPrintF(bufp,form2, bkp->ip,bkp->iz,bkp->ix);
  if(bkp->bp == 0){ return;}
  if(bkp->ip == 0){ return;}
  cp=bkp->bp;   ip=bkp->ip; 
  for(iic=0; iic<bkp->ix; iic++)
  { fPrintF(bufp,form3, iic,*ip,*(ip+1)-1 );
    len = *(ip+1) - *ip;  cp = bkp->bp + *ip;
    for(ibc=0; ibc<len; ibc++){ fPrintF(bufp, "%c", *(cp + ibc) );}
    fPrintF(bufp, "\n");
    ip++;
  }
  return ;
}
/******************************************************************************/
/*        chrBookXPrint(bkp)                                                  */
/******************************************************************************/
void      chrBookXPrint(FILE  *bufp, chrBook *bkp)
{ char     *cp;
  ptrdiff_t     *ip;
  size_t    ibc, iic, len;
  char      form1[]="book.bp=%p, book.bz=%d, book.bx=%d \n" ;
  char      form2[]="book.ip=%p, book.iz=%d, book.ix=%d \n" ;
  char      form3[]="char set number i=%d, ibeg=%d, iend=%d : \n" ;

  fPrintF(bufp,form1, bkp->bp,bkp->bz,bkp->bx);
  fPrintF(bufp,form2, bkp->ip,bkp->iz,bkp->ix);
  if(bkp->bp == 0){ return;}
  if(bkp->ip == 0){ return;}
  cp = bkp->bp;   ip = bkp->ip; 
  for(iic=0; iic<bkp->ix; iic++)
  { fPrintF(bufp,form3, iic,*ip,*(ip+1)-1 );
    len= *(ip+1) - *ip; cp=bkp->bp + *ip;
    for(ibc=0; ibc<len; ibc++){ fPrintF(bufp, "%2x", *(cp + ibc) ) ; }
    fPrintF(bufp, "\n");
    ip++;
  }
  return ;
}
/******************************************************************************/
/*        chrBookStrPrint(bkp)                                                */
/* print chrBook assuming that sets of char are all strings (null terminated) */
/******************************************************************************/
void      chrBookStrPrint(FILE  *bufp, chrBook *bkp)
{ char     *cp;
  ptrdiff_t     *ip;
  size_t    iic;
  char      form1[]="book.bp=%p, book.bz=%d, book.bx=%d \n" ;
  char      form2[]="book.ip=%p, book.iz=%d, book.ix=%d \n" ;
  char      form3[]="string number %d : \n" ;

  fPrintF(bufp,form1, bkp->bp,bkp->bz,bkp->bx);
  fPrintF(bufp,form2, bkp->ip,bkp->iz,bkp->ix);
  if(bkp->bp == 0){ return;}
  if(bkp->ip == 0){ return;}
  cp = bkp->bp;   ip = bkp->ip; 
  for(iic=0; iic<bkp->ix; iic++)
  { fPrintF(bufp,form3, iic,*ip ); cp = bkp->bp + *ip;
    fPrintF(bufp, "%s\n", cp ) ; 
    ip++;
  }
  return ;
}
/******************************************************************************/
/*        chrBookIPrint(bkp)                                                  */
/******************************************************************************/
void      chrBookIPrint(FILE  *bufp, chrBook *bkp)
{ ptrdiff_t     *ip;
  size_t    iic;
  char      form1[]= "book.bp=%p, book.bz=%d, book.bx=%d \n" ;
  char      form2[]= "book.ip=%p, book.iz=%d, book.ix=%d \n" ;
  char      form3[]= "indices list : ";
  char      form4[]= "%d, ";

  fPrintF(bufp,form1, bkp->bp,bkp->bz,bkp->bx);
  fPrintF(bufp,form2, bkp->ip,bkp->iz,bkp->ix);
  if(bkp->ip == 0){ return;}
  ip = bkp->ip; 
  fPrintF(bufp,form3);
  for(iic = 0;  iic <= bkp->ix;  iic++, ip++)
  { fPrintF(bufp,form4, (int)*ip );
  }
  fPrintF(bufp, "\n");
  return ;
}
/******************************************************************************/
/*        chrBookInc1Mem(bkp,cp,n)                                            */
/* add 1 set of char into a chrBook structure                                 */
/*  REMARK/WARNING : effective number of indices is number of string ix +1    */
/******************************************************************************/
void      chrBookInc1Mem(chrBook *bkp,char *cp,size_t n)
{ size_t    newix,newbx;
  char     *cbp;
  ptrdiff_t      ind;

  newix = bkp->ix + 1;  newbx = bkp->bx + n;
  chrIBookRealloc(bkp, newix+1 , bkp->iz /5);
  chrBBookRealloc(bkp, newbx   , bkp->bz /5);
  cbp = bkp->bp;  ind = *(bkp->ip + bkp->ix);  cbp += ind;
  memMove(cbp,cp,n);             bkp->bx = newbx;
  *(bkp->ip + newix) = ind + n;  bkp->ix = newix;
  return ;
}
/******************************************************************************/
/*        chrBookEnd(bkp,c)                                                   */
/* write int c (converted to char) in the last used memory                    */
/******************************************************************************/
void      chrBookEnd(chrBook *bkp,int c)
{ ptrdiff_t      ind;

  ind = *(bkp->ip + bkp->ix);
  *(bkp->bp + ind-1) = (char)c;
  return ;
}
/******************************************************************************/
/*        chrBookInc1Str(bkp,cp)                                              */
/* add 1 string (set of char 0 terminated) into a chrBook structure           */
/******************************************************************************/
void      chrBookInc1Str(chrBook *bkp,char *cp)
{ size_t    n;
  
  n = strLen(cp) +1;
  chrBookInc1Mem(bkp,cp,n);
  chrBookEnd(bkp,0);
  return ;
}
/******************************************************************************/
/*        chrBookEnlargeLast(bkp,cp,n)                                        */
/* enlarge the last set of char of the book by the n char pointed by cp       */
/******************************************************************************/
void      chrBookEnlargeLast(chrBook *bkp,char *cp,size_t n)
{ size_t    newbx;
  char     *cbp;
  ptrdiff_t      ind;

  newbx = bkp->bx + n;
  chrBBookRealloc(bkp, newbx, bkp->bz /5);
  cbp = bkp->bp;  ind = *(bkp->ip + bkp->ix);  cbp += ind;
  memMove(cbp,cp,n);  bkp->bx = newbx;  *(bkp->ip + bkp->ix) = ind + n;
  return;
}
/******************************************************************************/
/*        chrBookDecreaseLast(bkp,n)                                          */
/* decrease the last set of char of the book by  n  char                      */
/* It is user's responsability to choose n less or equal length of last set   */
/******************************************************************************/
void      chrBookDecreaseLast(chrBook *bkp,size_t n)
{ 
  bkp->bx -= n ;  *(bkp->ip + bkp->ix ) -= n ;
  return;
}
/******************************************************************************/
/*        chrBookRemoveLast(bkp)                                              */
/* remove the last set of char of the book                                    */
/******************************************************************************/
void      chrBookRemoveLast(chrBook *bkp)
{ if(bkp->ix == 0) return;
  bkp->ix -= 1 ;
  bkp->bx = *(bkp->ip + bkp->ix);
  return;
}


/******************************************************************************/
/*        chrBook2ptrVec                                                      */
/*                                                                            */
/* write in a ptrVec, the pointers to the strings of a chrBook.               */
/* assumes that  ptrVecp pVp  exists already.    ONLY realloc pVp->p.         */
/*                                                                            */
/* WARNING: if the base cBp->bp of the chrBook *cBp has changed,              */
/*          all pointers need to be recomputed                                */
/******************************************************************************/
ptrVecp   chrBook2ptrVec(ptrVecp pVp,chrBook *cBp)
{ char      *cp, *icp;
  ptrdiff_t *ip;
  size_t     ix;
  int        i;

  ix=cBp->ix;  ip=cBp->ip; icp=cBp->bp;
  pVp->p = ptrChkRealloc(pVp->p, &(pVp->z),ix,2, "chrBook2ptrVec" );
  pVp->x=0;
  for(i=ix; i>0; i--,ip++){ cp = icp+*ip;  ptrVecInc1(pVp,(void*)cp); }
  return pVp;
}
/******************************************************************************/
/*        chrBook2chrpp(chrBookp cBp)                                         */
/*                                                                            */
/* allocate and fill in  an array of pointers  to the strings of a chrBook    */
/* return a (char**)cpp  such that *cpp is the first pointer                  */
/*                                 *(cpp+1) is the second pointer ...         */
/*                                                                            */
/* WARNING: if the base cBp->bp of the chrBook *cBp has changed,              */
/*          all pointers need to be recomputed                                */
/******************************************************************************/
char    **chrBook2chrpp(chrBook *cBp)
{ char     *icp, **cpp, **cfpp;
  ptrdiff_t     *ip;
  size_t    ix;
  int       i;

  ix=cBp->ix;  ip=cBp->ip;  icp=cBp->bp;
  cpp = cfpp = (char**)ptrAlloc(ix, "chrBook2chrpp");
  for(i=ix; i>0; i--,ip++,cpp++){ *cpp = icp + *ip; }
  return   cfpp;
}
/******************************************************************************/
/******************************************************************************/
