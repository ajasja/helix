#include "helix.h"
#include <string.h>

int FindLipo (struct sequence *Sequence)
{
/******************************************************************
 * Description: this subroutine searches for lipophilic (dispersion) 
 * interacions according to the following patterns
 * L/F/I/M/V/W/Y/Nle XX L/F/I/M/V/W/Y/Nle
 * L/F/W/Y/M/V/I/Nle XXX L/F/W/Y/M/V/I/Nle
 * F/Y XXX H+
 * it returns 1 for successful search and acception, 0 for rejection
 ******************************************************************/

  int nIndx = 0;
  float g_w[MAXLEN], g_w1[MAXLEN]; /* multiplicative factor for w/w1 */
  int flag1 = OFF, flag2 = OFF,
      aFlag1[MAXLEN], aFlag2[MAXLEN], aFlag3[MAXLEN];
  char szLine[MAXSIZE];
  struct residue *Seq = Sequence->residue;
  int nLen = 0, nLength = Sequence->nLength;

  printf("\nSearch for lipophilic effects? Y/N (def=Y): "); 
  gets(szLine);
  if (YES)
  {
    for (nIndx = 1; nIndx < nLength; ++nIndx)
    { 
      g_w[nIndx] = g_w1[nIndx] = 1.0; /* initialize g_w/g_w1 */
      aFlag1[nIndx] = aFlag2[nIndx] = aFlag3[nIndx] = OFF;
    }
    for (nIndx = 1; nIndx < nLength - 3; ++nIndx)
    {
      if (strchr("LF", Seq[nIndx].szName[0]) &&
          strchr("LF", Seq[nIndx+3].szName[0]))
      /* L/F XX L/F */
      {
        g_w[nIndx]   *= 1.11; 
        g_w[nIndx+3] *= 1.11; 
        flag1 = aFlag1[nIndx] = ON;
      }
      else if ((strchr("LFIMVWY", Seq[nIndx].szName[0])   ||
               !strncmp("Nle",    Seq[nIndx].szName, 3))  &&
               (strchr("LFIMVWY", Seq[nIndx+3].szName[0]) ||
               !strncmp("Nle",    Seq[nIndx+3].szName, 3)))
      /* L/F/I/M/V/W/Y/Nle XX L/F/I/M/V/W/Y/Nle */
      {
        g_w[nIndx]   *= 1.06; 
        g_w[nIndx+3] *= 1.06; 
        flag1 = aFlag1[nIndx] = ON;
      }
      if (nIndx < nLength - 4 && 
          strchr("LF", Seq[nIndx].szName[0]) &&
          strchr("LF", Seq[nIndx+4].szName[0]))
      /* L/F XXX L/F */
      {
        g_w[nIndx]   *= 1.27; 
        g_w[nIndx+4] *= 1.27;
        flag1 = aFlag2[nIndx] = ON;
      }
      else if (nIndx < nLength - 4 && 
               strchr("WY", Seq[nIndx].szName[0]) &&
               strchr("WY", Seq[nIndx+4].szName[0]))
      /* W/Y XXX W/Y */
      {
        g_w[nIndx]   *= 1.05; 
        g_w[nIndx+4] *= 1.05;
        flag1 = aFlag2[nIndx] = ON;
      }
      else if (nIndx < nLength - 4 && 
               (strchr("IMV",  Seq[nIndx].szName[0])   ||
               !strncmp("Nle", Seq[nIndx].szName, 3))  &&
               (strchr("IMV",  Seq[nIndx+4].szName[0]) ||
               !strncmp("Nle", Seq[nIndx+4].szName, 3)))
      /* I/M/V/Nle XXX I/M/V/Nle */
      {
        g_w[nIndx]   *= 1.18; 
        g_w[nIndx+4] *= 1.18;
        flag1 = aFlag2[nIndx] = ON;
      }
      else if (nIndx < nLength - 4 && 
              (strchr("WYMVI", Seq[nIndx].szName[0])   ||
               !strncmp("Nle", Seq[nIndx].szName, 3))  &&
              (strchr("WYMVI", Seq[nIndx+4].szName[0]) ||
               !strncmp("Nle", Seq[nIndx+4].szName, 3)))
      /* W/Y/M/V/I/Nle XXX W/Y/M/V/I/Nle */
      {
        g_w[nIndx]   *= 1.15;
        g_w[nIndx+4] *= 1.15;
        flag1 = aFlag2[nIndx] = ON;
      }
      else if (nIndx < nLength - 4 && 
              (strchr("LFWYMVI", Seq[nIndx].szName[0])   ||
              !strncmp("Nle",    Seq[nIndx].szName, 3))  && 
              (strchr("LFWYMVI", Seq[nIndx+4].szName[0]) ||
              !strncmp("Nle",    Seq[nIndx+4].szName, 3)))
      /* L/F/W/Y/M/V/I/Nle XXX L/F/W/Y/M/V/I/Nle */
      {
        g_w[nIndx]   *= 1.22;
        g_w[nIndx+4] *= 1.22;
        flag1 = aFlag2[nIndx] = ON;
      }
    }

    for (nIndx = 1; nIndx < nLength - 1; ++nIndx)
    { 
      g_w1[nIndx] = g_w[nIndx];
      if (nIndx < nLength-4 && 
          strchr("FY", Seq[nIndx].szName[0]) &&
          'H' == Seq[nIndx+4].szName[0] && 
          '+' == Seq[nIndx+4].cCharge)
      /* F/Y XXX H+, aryl effect */
      {
        g_w[nIndx]    *= 1.4; 
        g_w[nIndx+2]  *= 1.28; 
        g_w1[nIndx+2] *= 1.28;
        flag1 = aFlag2[nIndx] = ON;
      }
      if (g_w[nIndx] > 2.2 && strchr("FY", Seq[nIndx].szName[0]))
        g_w[nIndx] = 2.2;
      else if (g_w[nIndx] > 1.8 && !strchr("FY", Seq[nIndx].szName[0]))
        g_w[nIndx] = 1.8;

      if(g_w1[nIndx] > 1.62) 
        g_w1[nIndx] = 1.62;
      if (ON == flag1 && OFF == flag2)
      {
        printf("Lipophillic interaction is found between:\n");
        flag2 = ON;
      }
      if (ON == aFlag1[nIndx]) /* ZXXZ */
      {
        printf(" %s(%d)\tand  %s(%d)\n", 
               Seq[nIndx].szName, 
               nIndx + 1, 
               Seq[nIndx+3].szName,
               nIndx + 4);
        aFlag1[nIndx] = OFF;
      }
      if (ON == aFlag2[nIndx]) /* ZXXXZ*/
      {
        if (strchr("LFWY", Seq[nIndx].szName[0]))
        {
          Seq[1].vN *= 1.2; 
          g_w[1] = g_w1[1] = 1.0;
          aFlag3[1] = ON; 
          nLen = MAX(nLen, 1);
        }
        g_w[nIndx+2] *= 1.16; 
        aFlag2[nIndx] = OFF;
        printf(" %s(%d)\tand  %s(%d)\n", 
               Seq[nIndx].szName,
               nIndx + 1,
               Seq[nIndx+4].szName,
               nIndx + 5);
      }
      if (g_w[nIndx] > 1.0)
      {
        Seq[nIndx].w1 *= g_w1[nIndx]; 
        Seq[nIndx].w  *= g_w[nIndx];
        aFlag3[nIndx] = ON; 
        nLen = MAX(nLen, nIndx);
      }
    }
    if (ON == flag1) 
    {
      PrintParam(Seq, 1, nLen, aFlag3);
      return 1;
    }
    else 
      printf(" No lipophilic interactions are found and accepted!\n");
  }

  return 0;
}
