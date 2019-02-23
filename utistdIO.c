/*   ../libmy/utistdIO.c                                                      */
/*   Mennessier Gerard             940822                                     */
/*   Last revised M.G.             961108                                     */

/******************************************************************************/
/*  sprintf: SUN   return value char *   (initial memory address always)      */
/*           ANSI               int      (number of transmitted char)         */
/***** -->  AVOID  sprintf  USE sPrintF                                   *****/
/*  fputc  : SUN   has first arg of type  char                                */
/*           ANSI                         int                                 */
/***** -->  better to use  fPutC                                          *****/
/*  fileno: SUN   return value char                                           */
/*          POSIX              int                                            */
/***** -->  better to use  fileNo                                         *****/
/******************************************************************************/

#include  "utistdIO.h"

#if  defined(SPARCM)|| !defined(DEF_MAPIO)
#ifdef  SPARCM
#include  <varargs.h>    /* Standard include for variable number of arguments */
extern  int vfprintf(FILE *stream, const char *format, va_list ap) ;
extern  int vprintf (const char *, va_list);
extern  int vsprintf(char *membufp, const char *format, va_list ap) ;

#else
#include  <stdarg.h>     /* Standard include for variable number of arguments */
#endif                                                        /* ifdef SPARCM */

/******************************************************************************/
/*                            fPutc                                           */
/* ~= fputc    stdio.h (ANSI)                                                 */
/* (int,FILE*) in ANSI ; (char,FILE*) in SPARC                                */
/******************************************************************************/
int      fPutC(int ci, FILE *stream)
{ char     c;
  c=ci;
  return  fputc(c,stream);
}

/******************************************************************************/
/*                            sPrintF                                         */
/* ~= sprintf  stdio.h (ANSI)                                                 */
/* print into memory under control of a Format                                */
/******************************************************************************/
int   sPrintF(char *membufp, const char *formatp, ...)
{ int      st;
  va_list  vargp;
 
#ifdef  SPARCM
  va_start(vargp) ;
#else  
  va_start(vargp, formatp);
#endif
  st= vsprintf(membufp,formatp,vargp);                      /* stdio.h (ANSI) */
  va_end(vargp);
  return  st;
}
#endif                            /* if  defined(SPARCM)|| !defined(DEF_MAPIO */

#ifndef  DEF_MAPIO
/******************************************************************************/
/*                            fPrintF                                         */
/* ~= fprintf  stdio.h (ANSI)                                                 */
/* print into a File under control of a Format                                */
/******************************************************************************/
int   fPrintF(FILE *stream, const char *formatp, ...)
{ int      st;
  va_list  vargp;
 
#ifdef  SPARCM
  va_start(vargp) ;
#else
  va_start(vargp, formatp);
#endif
  st= vfprintf(stream,formatp,vargp);                       /* stdio.h (ANSI) */
  va_end(vargp);
  return  st;
}

#endif                                                   /* ifndef  DEF_MAPIO */


/******************************************************************************/
/*                            fileNo                                          */
/* ~= fileno  stdio.h                                                         */
/* SPARCM :           it is just a define with result type char               */
/* POSIX (not ANSI) : it is a function    with result type int                */
/******************************************************************************/
int       fileNo(FILE *stream)
{ return (int)( fileno(stream) );
}

/******************************************************************************/
/******************************************************************************/
