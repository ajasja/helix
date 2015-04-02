#include "helix.h"
#include <string.h>

int FindRH (struct sequence *Sequence)
{
/************************************************************************
 * Description: this subroutine searches for R+ and H+ at C-terminus from 
 * n-4 to n-1 (end effects) and makes proper adjustments on w, C, vC
 * it returns 1 for successful search and acception, 0 for rejection
 ************************************************************************/

  int nIndx = 0;
  int flag = OFF, aFlag[MAXLEN];
  char szLine[MAXSIZE];
  struct residue *Seq = Sequence->residue;
  int nLength = Sequence->nLength;

  printf("\nSearching for Arg & His at C-terminus? Y/N (def=Y): ");
  gets(szLine);

  if (YES)
  {
    for (nIndx = nLength - 5; nIndx < nLength; ++nIndx)
      aFlag[nIndx] = OFF; /* initialize flags */

    if ('R' == Seq[nLength-5].szName[0] && '+' == Seq[nLength-5].cCharge)
    {
      printf("Arg is found at: %s(%d)\n", Seq[nLength-5].szName, nLength - 4);
      Seq[nLength-5].w += 0.1;
      flag = aFlag[nLength-5] = ON;
    }

    for (nIndx = nLength - 4; nIndx < nLength - 1; ++nIndx)
      if (strchr("RH", Seq[nIndx].szName[0]) && '+' == Seq[nIndx].cCharge)
      {
        printf("Arg or His is found at: %s(%d)\n", Seq[nIndx].szName, nIndx + 1);
        Seq[nIndx].w  += 0.2; 
        Seq[nIndx].vC += 0.004;
        flag = aFlag[nIndx] = ON;
        if ('R' == Seq[nIndx].szName[0] && nIndx != nLength - 4)
          Seq[nIndx].C *= 1.15;
        if ('H' == Seq[nLength-2].szName[0])
        {
          Seq[nLength-1].C *= 1.15; 
          aFlag[nLength-1] = ON;
        }
        if ('H' == Seq[nLength-3].szName[0])
        {
          Seq[nLength-2].vC *= 1.2; 
          Seq[nLength-3].vC *= 1.15; 
          Seq[nLength-2].C  *= 1.3;
          aFlag[nLength-2] = ON;
        }
        if ('H' == Seq[nLength-4].szName[0])
        {
          Seq[nLength-3].vC *= 1.2; 
          Seq[nLength-4].vC *= 1.1;
          aFlag[nLength-3] = ON;
        }
      }
    if (ON == flag)
    {
      PrintParam(Seq, nLength - 5, nLength - 2, aFlag);
      return 1;
    }
    else 
      printf(" No Arg or His are found at C-terminus!\n");
  }

  return 0;
}
