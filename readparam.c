#include "helix.h"
#include <string.h>
#include <stdlib.h>

int ReadParam (struct sequence *Sequence)
{
/**************************************************************
 * Description: this subroutine reads the parameter (.dat) file
 * it returns 1 for success, 0 for failure
 **************************************************************/

  int nIndx = 0;
  char szLine[MAXSIZE];
  char* pszLine = 0;
  FILE *pDatFile = 0;
  struct residue *Seq = Sequence->residue;

  pDatFile = fopen(Sequence->szDatFile, "r");
  while (nIndx < Sequence->nLength)
  {
    fseek(pDatFile, 0, SEEK_SET);
    fgets(szLine, MAXSIZE, pDatFile);
    while (strcmp(Seq[nIndx].szName, strtok(szLine, "=")))
      if (!feof(pDatFile)) 
        fgets(szLine, MAXSIZE, pDatFile);
      else
      {
        printf("\aSequence error: %s undefined in %s\n",
               Seq[nIndx].szName, Sequence->szDatFile);
        return 0;
      }
    pszLine = strtok(0, "{");  /* a bug in strtok? */
    Seq[nIndx].w1 = atof(strtok(pszLine, ","));
    Seq[nIndx].w  = atof(strtok(0, ","));
    Seq[nIndx].vN = atof(strtok(0, ","));
    Seq[nIndx].vC = atof(strtok(0, ","));
    Seq[nIndx].N  = atof(strtok(0, ","));
    Seq[nIndx].C  = atof(strtok(0, "}"));
    ++nIndx;
  }
  fclose(pDatFile);

  return 1;
}
