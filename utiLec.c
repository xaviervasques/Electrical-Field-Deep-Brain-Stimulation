/*  ../libmy/utiLec.c                                                         */
/*  Mennessier Gerard                   941206                                */
/*  Last Revised : M.G.                 991025                                */

#include  <stddef.h>
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utistdMem.h"
#include  "utistdStr.h"
#include  "utiStr.h"
#include  "utiAlloc.h"
#include  "utiBookChr.h"
#include  "utiLec.h"

/******************************************************************************/
/*                              lec1all                                       */
/*                                                                            */
/* read the complete File *stream exactly as it is                            */
/*         (keeping NewLine's,...,and without adding ending 0 char)           */
/*         Can be used for binary files                                       */
/*                                                                            */
/* if reading stdin , the end of file is        Cntrl-D                       */
/*                                                                            */
/* output: lec1all pointer on the buffer,                                     */
/*         nb_bytep pointer on the effective size                             */
/******************************************************************************/
char     *lec1all(FILE *stream, size_t *nb_bytep)
{ char     *bufp;
  size_t    buftz;
  size_t    bufx,buftx;
  size_t    nb;
  static    char      prognam[]="lec1all" ;

  nb = 1024;  buftx = 0;  buftz = nb;
  bufp = chrAlloc(nb,prognam);
  while(  (bufx = fread(bufp+buftx, sizeof(char), nb, stream) )  )
  { buftx += bufx;
    if( buftx > 8*nb ) nb = 256 * (buftx / 1024);
    bufp = chrChkRealloc(bufp,&buftz,nb+buftx,buftx/8,prognam);
  }
  *nb_bytep = buftx;  return bufp;
}
/******************************************************************************/
/*                              lec2all                                       */
/*                                                                            */
/* read the complete File *stream                                             */
/* input : a line terminated by a Backslash \ ,means that it                  */
/*         continues on the next line.                                        */
/*         In this case the Baclslash and the immediately following NewLine   */
/*         are skipped.                                                       */
/*         Otherwise each NewLine is OverWrited by a ZERO char                */
/*         and is the last character of the string.                           */
/* output: a list of strings (ZERO terminated),                               */
/*         more precisely a pointer on a chrBook                              */
/* end   : either  2 succesives NewLine i.e 1 EMPTY LINE, or end of the file  */
/******************************************************************************/
chrBook  *lec2all(FILE *stream)
{ int       new,newstr;        /** newstr=1 if new str,  newstr=o if continuation **/
  size_t    bufz = 1024;
  chrBook  *bookp;
  char     *bufp, *c1p;
  ptrdiff_t      bufx;
  int       NLi = '\n' ;
  char      BSL = '\\' ;
  static    char      prognam[] ="lec2all" ;

  bufp = chrAlloc(bufz,prognam);
  bookp = chrBookAlloc(1,prognam); 
  chrBBookAlloc(bookp,bufz);  chrIBookAlloc(bookp,10);
  newstr = 1; 
  while( fgets(bufp, (int)bufz, stream) != NULL)
  { c1p = memChr(bufp,NLi,bufz);
    if(c1p == NULL){ bufx = bufz -1;  new =0;}
    else
    { if( *(c1p-1) == BSL){ bufx = c1p - bufp -1;  new =0;}
      else                { bufx = c1p - bufp +1;  new =1;  *c1p=0;}
    }
    if(new && (bufx == 1) ){ break;}                   /** STOP si 2 NLi de suite **/
    if(newstr) chrBookInc1Mem(bookp,bufp,bufx);
    else       chrBookEnlargeLast(bookp,bufp,bufx);
    newstr = new;
  }
  free(bufp);
  return bookp;
}
/******************************************************************************/
/*                              lec3all                                       */
/* read the complete File *stream                                             */
/* input : a line terminated by a Backslash \ ,means that it                  */
/*         continues on the next line.                                        */
/*         In this case the Baclslash and the immediately following NewLine   */
/*         are skipped.                                                       */
/*         Otherwise each NewLine is OverWrited by a ZERO char                */
/*         and is the last character of the string.                           */
/* output: a list of strings (ZERO terminated),                               */
/*         more precisely a pointer on a chrBook                              */
/* end   : end of the file                                                    */
/*         or   Cntrl-D       for stdin                                       */
/******************************************************************************/
chrBook  *lec3all(FILE *stream)
{ int       new,newstr;        /** newstr=1 if new str,  newstr=0 if continuation **/
  size_t    bufz = 1023;
  chrBook  *bookp;
  char     *bufp, *c1p;
  ptrdiff_t    bufx;
  int       NLi = '\n' ; 
  char      BSL = '\\' ;
  static    char      prognamp[] = "utiLec::lec3all";   /** static added 20010918 **/

/*
fPrintF(stderr,"%s : BEGIN with streamp=%p\n", prognamp, stream); 
*/
  bufp  = chrAlloc(bufz +1, prognamp);
  *bufp = 0;  bufp++;               /** Thus *(c1p-1) always exist if c1p != NULL **/
  bookp = chrBookAlloc(1, prognamp); 
  chrBBookAlloc(bookp,bufz);  chrIBookAlloc(bookp,10);
/*
fPrintF(stderr,"%s : bufp et bookp allocated\n", prognamp); 
*/
  newstr = 1; 
  while( fgets(bufp, (int)bufz, stream) != NULL)
  { 
/*
fPrintF(stderr,"lec3all :%s",bufp); 
*/
    c1p = memChr(bufp,NLi,bufz);
    if(c1p == NULL){ bufx = bufz-1;  new = 0;}
    else
    { if( *(c1p-1) == BSL){ bufx = c1p - bufp -1;  new = 0;}
      else                { bufx = c1p - bufp +1;  new = 1;  *c1p =0;}
    }
/* fPrintF(stderr,"lec3all, bufx=%d,bufz=%d,newstr=%d,new=%d\n", 
                                             (int)bufx,(int)bufz,newstr,new); */
    if(newstr) chrBookInc1Mem(bookp,bufp,bufx);
    else       chrBookEnlargeLast(bookp,bufp,bufx);
    newstr = new;
/*
fPrintF(stderr,"lec3all, bufx=%d,bk.bz=%d,bk.bx=%d,bk.iz=%d,bk.ix=%d,value:%s;\n",
                                         (int)bufx,(int)(bookp->bz),(int)(bookp->bx),
                                            (int)(bookp->iz),(int)(bookp->ix),bufp );
*/
  }
  bufp--;  free(bufp);
  return bookp;
}
/******************************************************************************/
/*                              leckbd                                        */
/*                                                                            */
/* read input from stream KEYBOARD                                            */
/* input : a line terminated by a Backslash \ ,means that it                  */
/*         continues on the next line.                                        */
/*         In this case the Baclslash and the immediately following NewLine   */
/*         are skipped.                                                       */
/*         Otherwise each NewLine is OverWrited by a ZERO char                */
/*         and is the last character of the string.                           */
/* output: a list of strings (ZERO terminated),                               */
/*         more precisely a pointer on a chrBook                              */
/* end   : either  2  EMPTY LINE, or given END-line                           */
/******************************************************************************/
chrBook  *leckbd(FILE *stream)
{ int       new,newstr;        /** newstr=1 if new str,  newstr=o if continuation **/
  int       existEND = 0,  count =2;
  char      endlinep[64] ="";
  int       endlinez = 64, endlinex;
  size_t    bufz = 1024;
  chrBook  *bookp;
  char     *bufp, *cp;
  ptrdiff_t    bufx;
  int       NLi ='\n' ; 
  char      BSL ='\\' ;
  char      messp[] = "give the END-line, ELSE type the ENTER/RETURN key" 
                                                     " (default is 2 empty lines)." ;
  char      prognam[]="leckbd" ;

  bufp = chrAlloc(bufz,prognam);
  bookp = chrBookAlloc(1,prognam); 
  chrBBookAlloc(bookp,bufz);  chrIBookAlloc(bookp,10);
  fPrintF(stdout,"%s\n",messp);
  cp = fgets(endlinep, endlinez, stdin);
  endlinex = strLen(endlinep) -1;
  if(endlinex){ existEND = 1;  endlinep[endlinex] = 0;}
  fPrintF(stdout, "          TYPE your input on the keyboard.\n");

  newstr = 1; 
  while(fgets(bufp,(int)bufz,stream) != NULL)
  { cp = memChr(bufp,NLi,bufz);
    if(cp == NULL){ bufx = bufz -1;  new =0;}
    else
    { if( *(cp-1) == BSL){ bufx = cp-bufp-1; new =0;}
      else               { bufx = cp-bufp+1; new =1;  *cp =0;}
    }

    if(newstr)
    { if(existEND){ if(strEq(bufp,endlinep)){ break;} }
      else        { if(endlinep[0] == NLi){ count--;  if(count == 0){ break;}}}
      chrBookInc1Mem(bookp,bufp,bufx);
      if(bufx != 1){ count = 2;}
    }
    else{ chrBookEnlargeLast(bookp,bufp,bufx); }
    newstr = new;
  }
  free(bufp);
  return bookp;
}

/******************************************************************************/
/*        getFile(char *form1)                                                */
/******************************************************************************/
FILE     *getFileStream(char *form1,char *form2,char *form4)
{ char      namep[128], *cp;
  int       namez=128;
  int       namex,max=10;
  FILE     *bufIN;
  FILE     *bufOUT = stdout ;  
  char      formap[] = "The file is :%s\n" ;
  char      formbp[] = "getFileStream" ;
  char      formcp[] = " cannot open file STOP\n" ;

  fPrintF(bufOUT,form2,form1); *namep=0;
  cp = NULL ;
  while(cp == NULL)
  { cp = fgets(namep,namez,stdin); 
/*    fPrintF(stderr,"cp=%p,namep=%s; \n",cp,namep);  */
    max--;
    if(max == 0) exit(-1);
  }
  namex = (int)strLen(namep);
  if(namex){ namex--;  *(namep+namex) = 0; }
  if(namex){ fPrintF(bufOUT,formap,namep);  bufIN = fopen(namep,"r");
             if(bufIN == NULL){ myErr1(-1,stderr,formbp,formcp); }
           }
  else     { fPrintF(bufOUT,form4); bufIN = stdin; 
           }
  fflush(bufOUT);
  return bufIN;
}

/******************************************************************************/
/*        getFileNameAndStream(char *mess1,char *mess2)                       */
/*                                                                            */
/* asks  for the File name, then open the corresponding FILE                  */
/*                          (stdin if NO name)                                */
/* ex : getFileNameAndStream("give the name of the file",                     */
/*                           " or type RETURN for direct keyboard input")     */
/* input  : 2 message strings.                                                */
/* return : the FILE pointer                                                  */
/******************************************************************************/
FILE     *getFileNameAndStream(char *mess1,char *mess2)
{ char      namep[128],  *cp;
  int       namez = 128;
  int       namex,  max = 10;
  FILE     *bufIN;
  char      formap[] = "The file is :%s\n" ;
  char      formbp[] = "utiLec::getFileNameAndStream" ;
  char      formcp[] = " cannot open file STOP\n" ;

  fPrintF(stdout,"%s%s\n", mess1,mess2);  *namep = 0;
  cp = NULL;
  while(cp == NULL)
  { cp = fgets(namep,namez,stdin);
    max--;  if(max == 0) exit (-1);
  }
  namex = (int)strLen(namep);
  if(namex){ namex--;  *(namep + namex) = 0;}
  if(namex)
  { fPrintF(stdout,formap,namep);  bufIN = fopen(namep,"r");
    if(bufIN == NULL) myErr1(-1,stderr,formbp,formcp);
  }
  else     
  { bufIN = stdin; 
  }
  fflush(stdout);
  return bufIN;
}

/******************************************************************************/
/*        lec2(char *form1)                                                   */
/******************************************************************************/
chrBook  *lec2(char *form1)
{ FILE     *bufIN;
  chrBook  *bookp;
  char      form2[]="%s          "
                  "ELSE type the RETURN key, then type your input on the keyboard\n";
  char      form3[]="          TYPE your input on the keyboard."
                     "  Twice RETURN to stop.\n";

  bufIN = getFileStream(form1,form2,form3);  bookp = lec2all(bufIN); 
  if(bufIN != stdin){ fclose(bufIN);}
  return bookp;
}
/******************************************************************************/
/*        lec3(char *form1)                                                   */
/*                                                                            */
/*  Cntrl-D  must be the character defined by eof for the current tty         */
/*  it can be checked by  stty -a     which should display   eof = ^d;        */
/*            set by      stty eof  ^d                                        */
/*  it puts an EOF on stdin which is reinitialized here by fseek              */
/******************************************************************************/
chrBook  *lec3(char *form1)
{ FILE     *bufIN;
  chrBook  *bookp;
  char      form2[]="%s          "
                  "ELSE type the RETURN key, then type your input on the keyboard\n";
  char      form3[]="          TYPE your input on the keyboard."
                     "  RETURN Cntrl-D to stop.\n";

  bufIN = getFileStream(form1,form2,form3);
  bookp = lec3all(bufIN); 
  if(bufIN != stdin){ fclose(bufIN);}
  else              { fseek(stdin,0,SEEK_SET);}
  return bookp;
}
/******************************************************************************/
/******************************************************************************/
