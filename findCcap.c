#include "helix.h"
#include <string.h>

int FindCcap (struct sequence *Sequence)
{
/****************************************************************
 * Description: this subroutine searches for C-caping boxes (Gly)
 * at C terminus (from n-4 to n), 
 * it returns 1 for successful search and acception 0 for rejection
 ****************************************************************/

  int nIndx = 0;
  int flag = OFF, aFlag[MAXLEN];
  char szLine[MAXSIZE];
  char szStr[] = "Adjust parameters to reflect identified interaction? Y/N (def=Y): ";
  struct residue *Seq = Sequence->residue;

  for (nIndx = Sequence->nLength - 5; nIndx < Sequence->nLength; ++nIndx)
  {
    aFlag[nIndx-1] = OFF; /* initialize flags */
    if ('G' == Seq[nIndx].szName[0])
    {
      printf("A C-cap is found at: %s(%d)\n", Seq[nIndx-1].szName, nIndx);
      printf(szStr); 
      gets(szLine);
      if (YES)
      {
        Seq[nIndx-1].C *= 1.09;
        if (Seq[nIndx+1].cCharge != '+')   
          Seq[nIndx-1].C *= 1.06;
        if ('\0' == Seq[nIndx+1].szName[0]) 
          Seq[nIndx-1].C *= 0.94;
        else if (strchr("GNQSTDCY", Seq[nIndx+1].szName[0]) ||
                !strcmp("NH2", Seq[nIndx+1].szName))
          Seq[nIndx-1].C *= 1.15;
        if (Seq[nIndx-5].cCharge != '-' && 
            Seq[nIndx-4].cCharge != '-' &&
            Seq[nIndx-3].cCharge != '-')
          Seq[nIndx-1].C *= 1.1;
        flag = aFlag[nIndx-1] = ON;
      }
    }
  }
  if (ON == flag)
  {
    PrintParam(Seq, Sequence->nLength - 6, Sequence->nLength - 1, aFlag);
    return 1;
  }
  else
    printf(" No C-caps are found and accepted!\n");

  return 0;
}
