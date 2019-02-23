/*  URMAE/orientHaut/linear4.GL.V1/command.linear4.gl.                        */
/*  Mennessier Gerard                 20010507                                */
/*  Last Revised : G.M.               20020430                                */

#include  <stddef.h>
#include  <math.h>                                /** for absolute value : fabs() **/
#include  "utistdIO.h"
#include  "utistdErr.h"
#include  "utistdStr.h"
#include  "utiStr.h"
#include  "utiVecChr.h"

#include  "command.linear4.gl.h"

/******************************************************************************/
/*                                                                            */
/* main.linear   -newp newpinmods  -newm newneckmod  -ps PSFileName           */
/*               -rb BasisFilesRootName                                       */
/*               -o OutFileName                                               */
/*                                                                            */
/* where                                                                      */
/*   newpinmods = (0|m|p)^4                                                   */
/*   newneckmod = (MONO | BI)                                                 */
/*                                                                            */
/* Ex: main.linear  -newp 0m0p  -newm BI                                      */
/*                  -ps zzz.0p0p.ps                                           */
/*                  -o zzz.0p0p.out                                           */
/******************************************************************************/
/******************************************************************************/
/**  (filenamepp|bufpp)[10] =                                                **/
/**  Basisp[0], ...,Basisp[3], Outp, WritePSp, GridZ, GridR,                 **/
/******************************************************************************/
void      commandLine_linear_gl(FILE **bufpp, char **fileNamepp,
                              char *newpinmodstrp, char *newneckmodstrp, int *bip,
                                                            int argx, char *argvpp[])
{ int       i, is, strx;
  FILE     *bufOutp, *bufGridZp, *bufGridRp, *bufReadBasep;
  char     *nameOutp, *nameWritePSp;
  char     *zVerstrp = NULL, *rVerstrp = NULL;
  char     *rBaseDirp = NULL, *nameBaseFilep = NULL;
  char     *fullNameBaseFilep = NULL;
  char     *argp;
  int       bi = 0;
  chrVec   *chrVp = NULL;
  chrVec   *nameBaseFileVChrp = NULL;

  static    char      neckModeMONOp[] = "MONO",  neckModeBIp[] = "BI";
  static    char     *neckModespp[2] = {neckModeMONOp, neckModeBIp};
  static    char      nameOutDefp[] = "stdout";
  static    char      nameWritePSDefp[] = "zZ0000.ps";

  static    char      ZRVerstrpDefp[] = "0";
  static    char      ReadBaseDirDefp[] = ".";
  static    char      nameBaseFileDefp[] = "0000.NECK.MONO.00000000.0.0.txt";

  static    char      nameGridZp[] = "gridZ.linear4.INI.0";
  static    char      nameGridRp[] = "gridR.linear4.INI.NECK.0";
/*
  static    char      nameCylp[]  = "../linear4.Common/cylPot.linear4.INI.NECK.MONO";
  static    char      nameEltrdp[] =
                              "../linear4.Common/eltrdPot.linear4.INI.pppp.00000000";
*/

/*  static    chrVec    nameBaseFileVChrp[4]; */
  static    char      form1p[] = "  Cannot open %s = %s file\n";
  static    char      prognamp[] = "command.linear.gl::linear_command_line";


  fPrintF(stderr,"%s  BEGIN, argx=%d\n", prognamp, argx);
  if(argx < 1)
  { fPrintF(stderr,"argx=%d, MISSING arguments \n", argx);
    fPrintF(stderr,"SYNOPSIS: main.linear"
                              "  -newp newpinmods  -newm newneckmod  -ps PSFileName"
                                      "  -rb BasisFilesRootName  -o OutFileName \n");
    exit (1);
  }

  nameOutp = nameOutDefp;  bufOutp = NULL;
fPrintF(stderr," Default nameOut=%s \n", nameOutp);
  zVerstrp = rVerstrp = ZRVerstrpDefp;
fPrintF(stderr," Default zVersion=%s, rVersion=%s;\n", zVerstrp, rVerstrp);
  rBaseDirp = ReadBaseDirDefp;
  nameBaseFilep = nameBaseFileDefp;
fPrintF(stderr," Default rBaseDir=%s, nameBaseFile=%s;\n", rBaseDirp,nameBaseFilep);
  strcpy(newneckmodstrp, neckModespp[0]);
fPrintF(stderr," Default neckMode=%s \n", newneckmodstrp);
  strcpy(newpinmodstrp, "0000");
fPrintF(stderr," Default pinModes=%s \n", newpinmodstrp);

  nameWritePSp = nameWritePSDefp;
fPrintF(stderr," Default nameWritePS=%s;\n", nameWritePSp);

fPrintF(stderr," Default nameOut=%s \n", nameOutp);
fPrintF(stderr," Default zVersion=%s, rVersion=%s;\n", zVerstrp, rVerstrp);
fPrintF(stderr," Default rBaseDir=%s, nameBaseFile=%s;\n", rBaseDirp,nameBaseFilep);

  i = 1 ;
  while(i < argx)
  { argp = argvpp[i];
    fPrintF(stderr,"argvpp[%d]=%s ", i,argp) ;
    if(*argp == '-')
    { argp++;
      if(*argp == 'o')
      { i++;  nameOutp = argvpp[i];
        fPrintF(stderr," %s \n", nameOutp);
        bufOutp = fopen(nameOutp, "w");
      }
      else if(strEq(argp,"newp"))
      { i++;
        strcpy(newpinmodstrp, argvpp[i]);
        fPrintF(stderr," pinmodes=%s \n", newpinmodstrp);
      }
      else if(strEq(argp,"newm"))
      { i++;
        strcpy(newneckmodstrp, argvpp[i]);
        fPrintF(stderr," neck mode=%s \n", newneckmodstrp);
             if(strEq(newneckmodstrp, neckModespp[0]) ) bi = 0;
        else if(strEq(newneckmodstrp, neckModespp[1]) ) bi = 1;
      }
      else if(strEq(argp,"ps"))
      { i++;  nameWritePSp = argvpp[i];
        fPrintF(stderr," nameWritePSp = %s \n", nameWritePSp);
      }
    }
    i++ ;
  }
/*
goto JUMP1;
JUMP1:
*/

fPrintF(stderr,"nameOut=%s \n", nameOutp);
fPrintF(stderr,"zVersion=%s, rVersion=%s;\n", zVerstrp, rVerstrp);
fPrintF(stderr,"rBaseDir=%s, nameBaseFile=%s;\n", rBaseDirp,nameBaseFilep);
  fileNamepp[4] = nameOutp;
  if(bufOutp == NULL){ bufOutp = stdout;}
  bufpp[4] = bufOutp;

fPrintF(stderr,"nameWritePSp=%s \n", nameWritePSp);
  fileNamepp[5] = nameWritePSp;

  nameBaseFileVChrp = chrVecAlloc(4, prognamp);
fPrintF(stderr,"  nameBaseFileVChrp allocated at %p\n", (void*)nameBaseFileVChrp);
  chrVp = nameBaseFileVChrp;
fPrintF(stderr,"chrVp allocated at %p\n", (void*)chrVp);

/*
chrVecXPrint(stderr, chrVp);
chrVecPrint(stderr, chrVp);
*/

  for(is = 0;  is < 4;  is++, chrVp++)
  { fPrintF(stderr,"  is=%d; ", is);
    fPrintF(stderr,"  BaseFile is=%d  %s;\n", is, rBaseDirp);
    chrPVecAlloc(chrVp, 64);
    chrVecIncStr(chrVp, rBaseDirp);
    nameBaseFileVChrp[is].x--;
    chrVecIncStr(nameBaseFileVChrp + is, "/");
    nameBaseFileVChrp[is].x--;
chrVecXPrint(stderr, nameBaseFileVChrp + is);
chrVecPrint(stderr, nameBaseFileVChrp + is);
    strx = nameBaseFileVChrp[is].x;
    chrVecIncStr(nameBaseFileVChrp + is, nameBaseFilep);
    fullNameBaseFilep = nameBaseFileVChrp[is].p;
    *(fullNameBaseFilep + strx + is) = 'p' ;
fPrintF(stderr,"fileBaseName index %d = %s \n", is, fullNameBaseFilep);
    fileNamepp[is] = fullNameBaseFilep;
    bufReadBasep = fopen(fullNameBaseFilep, "r");
    if(bufReadBasep == NULL)
    { myErr1(-1,stderr,prognamp,form1p,"BaseFile", fullNameBaseFilep);}
    bufpp[is] = bufReadBasep;
  }

fPrintF(stderr,"nameGridZ=%s \n", nameGridZp);
  fileNamepp[6] = nameGridZp;
  bufGridZp = fopen(nameGridZp, "r");
  if(bufGridZp == NULL)
  { myErr1(-1,stderr,prognamp,form1p,"gridZ INIT", nameGridZp);}
  bufpp[6] = bufGridZp;

fPrintF(stderr,"nameGridR=%s \n", nameGridRp);
  fileNamepp[7] = nameGridRp;
  bufGridRp = fopen(nameGridRp, "r");
  if(bufGridRp == NULL)
  { myErr1(-1,stderr,prognamp,form1p,"gridR INIT", nameGridRp);}
  bufpp[7] = bufGridRp;

  *bip = bi;

fPrintF(stderr," Effective neckMode=%s \n", newneckmodstrp);
fPrintF(stderr," Effective pinModes=%s \n", newpinmodstrp);
fPrintF(stderr," Effective nameOut=%s \n", nameOutp);

  fPrintF(stderr,"%s  END\n", prognamp);
  return;
}
/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      nameWritePSsetPinmod(char *nameWritePSp, char *newPinmodesp)
{ char     *cp;

  cp = nameWritePSp + strLen(nameWritePSp) - 7;
  *cp++ = *newPinmodesp++;    *cp++ = *newPinmodesp++;    *cp++ = *newPinmodesp++;
  *cp   = *newPinmodesp;  
  return;
}

/******************************************************************************/
/*                                                                            */
/******************************************************************************/
void      nameWritePSsetObs(char *nameWritePSp, char *newObsp)
{ char     *cp;

  cp = nameWritePSp +  strLen(nameWritePSp) - 8;
  *cp = *newObsp;  
  return;
}
/******************************************************************************/
/******************************************************************************/
