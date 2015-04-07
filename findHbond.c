#include "helix.h"
#include <string.h>

int FindHbond (struct sequence *Sequence, bool Hbond_flag)
{
/*******************************************************************
 * Description: this subroutine searches for hydrogen-bonding effect
 * (side-chain/backbone or side-chain/side-chain H bonding) according
 * to the following patterns:
 * Q/E/D/K+ XXX E/D/K+/R+
 * E0 XXX E0
 * E0 XX E0
 * it returns 1 for successful search and acception, 0 for rejection
 *******************************************************************/

  int nIndx = 0, nLen = 0;
  float g_w[MAXLEN];
  int flag = OFF, aFlag[MAXLEN];
  char szLine[MAXSIZE];
  char szStr[] = "Hydrogen-bonding effect is found between:\n";
  struct residue *Seq = Sequence->residue;
  int nLength = Sequence->nLength;

  if(Hbond_flag)
    strcpy(szLine, "y");
  else
  {
    printf("\nSearch for hydrogen-bonding effects? Y/N (def=Y): "); 
    gets(szLine);
  }
  if (YES)
  {
    for (nIndx = 1; nIndx < nLength; ++nIndx)
    {
      g_w[nIndx] = 1.0;
      aFlag[nIndx] = OFF;
    }
    for (nIndx = 1; nIndx < nLength - 3; ++nIndx)
    {
      if(nIndx < nLength - 4 && 
         ('Q' == Seq[nIndx].szName[0]          &&
          strchr("ED", Seq[nIndx+4].szName[0]) ||
          strchr("ED", Seq[nIndx].szName[0])   &&
          strchr("KR", Seq[nIndx+4].szName[0]) && 
          '+' == Seq[nIndx+4].cCharge          ||
          'K' == Seq[nIndx].szName[0]          && 
          '+' == Seq[nIndx].cCharge            && 
          'E' == Seq[nIndx+4].szName[0]))
    /* Q/E/D/K+ XXX E/D/K+/R+ */
      {
        if(OFF == flag) 
          PRINT_HEADER(szStr, nIndx, nIndx + 4);
        g_w[nIndx+1] *= 1.2; 
        g_w[nIndx+3] *= 1.2;
        flag = ON; 
        nLen = MAX(nLen, nIndx + 3);
      } 
      else if (nIndx < nLength - 4           && 
               'E' == Seq[nIndx].szName[0]   &&
                0  == Seq[nIndx].cCharge     &&
               'E' == Seq[nIndx+4].szName[0] && 
                0  == Seq[nIndx+4].cCharge)
    /* E0 XXX E0 */
      {
        if (OFF == flag) 
          PRINT_HEADER(szStr, nIndx, nIndx+4);
        g_w[nIndx+1] *= 1.1; 
        g_w[nIndx+3] *= 1.1;
        flag = ON; 
        nLen = MAX(nLen, nIndx + 3);
      }
      if (nIndx < nLength - 3           && 
          'E' == Seq[nIndx].szName[0]   && 
           0  == Seq[nIndx].cCharge     &&
          'E' == Seq[nIndx+3].szName[0] && 
           0  == Seq[nIndx+3].cCharge)
    /* E0 XX E0 */
      {
        if (OFF == flag) 
          PRINT_HEADER(szStr, nIndx, nIndx + 3);
        g_w[nIndx+2] *= 1.26;
        flag = ON; 
        nLen = MAX(nLen, nIndx + 2);
      }
    }
    if (ON == flag)
    {
      printf("The parameters after adjustments are:\n");
      printf("\t   w'\t     w\t       vN\t vC\t   N\t     C\n");
      for (nIndx = 1; nIndx <= nLen; ++nIndx)
      {
        if (g_w[nIndx] > 1.0)
        {
          if (g_w[nIndx] > 1.2) 
            g_w[nIndx] = 1.2;
          Seq[nIndx].w *= g_w[nIndx];
          printf("%2d %s\t%f  %f  %f  %f  %f  %f\n",
                 nIndx + 1,
                 Seq[nIndx].szName,
                 Seq[nIndx].w1,
                 Seq[nIndx].w,
                 Seq[nIndx].vN,
                 Seq[nIndx].vC,
                 Seq[nIndx].N,
                 Seq[nIndx].C);
        } 
      }
      return 1;
    }
    else 
      printf(" No hydrogen-bond effects are found and accepted!\n");
  }
  return 0;
}
