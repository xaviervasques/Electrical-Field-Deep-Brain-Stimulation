/*   ../libmy/utistdErr.c                                                     */
/*   Mennessier Gerard                  19940822                              */
/*   Last revised M.G.                  20040503                              */

#include  "utistdErr.h"
#include  "utistdIO.h"

#ifdef  SPARCM
#include  <varargs.h>    /* Standard include for variable number of arguments */
extern  int  vfprintf(FILE *stream, const char *format, va_list ap) ;
extern  int  vprintf (const char *, va_list);
extern  int  vsprintf(char *membufp, const char *format, va_list ap) ;

#else
#include  <stdarg.h>     /* Standard include for variable number of arguments */
#endif

extern    void      exit(int status);   /* to avoid including stdlib.h (ANSI) */
                                                    /* status=0 means Success */

/******************************************************************************/
/*                           myErr0                                           */
/******************************************************************************/
void      myErr0(int exitcod, FILE *stream, char *formatp, ...)
{ va_list   vargp;
 
#ifdef  SPARCM
  va_start(vargp) ;
#else  
  va_start(vargp, formatp);
#endif
  (void) vfprintf(stream, formatp, vargp);
  va_end(vargp);
  exit(exitcod);
}
/******************************************************************************/
/*                           myErr1                                           */
/******************************************************************************/
void      myErr1(int exitcod, FILE *stream, char *prognamp, char *formatp, ...)
{ va_list   vargp;
  static char    form1p[] = "ERROR in %s: ";
  
  fPrintF(stream, form1p, prognamp);
 
#ifdef  SPARCM
  va_start(vargp) ;
#else
  va_start(vargp, formatp);
#endif
  (void) vfprintf(stream, formatp, vargp);
  va_end(vargp);
  exit(exitcod);
}
/******************************************************************************/
/*                           myErr2                                           */
/******************************************************************************/
void      myErr2(int exitcod, FILE *stream, char *srcfilenamp, char *prognamp,
                                                                  char *formatp, ...)
{ va_list   vargp;
  static char    form1p[] = "ERROR in %s::%s ";

  fPrintF(stream, form1p, srcfilenamp, prognamp);
 
#ifdef  SPARCM
  va_start(vargp) ;
#else
  va_start(vargp, formatp);
#endif
  (void) vfprintf(stream, formatp, vargp);
  va_end(vargp);
  exit(exitcod);
}
/******************************************************************************/
/******************************************************************************/
