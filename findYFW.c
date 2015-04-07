#include "helix.h"
#include <string.h>

int FindYFW(struct sequence *Sequence, bool YFW_flag)
{
/*************************************************************************
 * Description: this subroutine search for Tyr/Phe/Trp interactions
 * (lipophilic interactions) by matching the sequence against the following 
 * patterns then makes proper adjustments on w and w1
 * it returns 1 for sucessful search and acceptance, 0 for rejection
 *
 * Y XXX V/I/Nle (X for intermediate residues)
 * Y XXX L
 * L XXX Y/W
 * V/I XXX Y/W
 * Y XXX F/W/M/Nle/H0
 * F/W XXX L/H0
 * M/Nle XXX F/Y
 * F XXX M
 *************************************************************************/

  int nIndx = 0, nLen = 0;
  int flag = OFF, flag1 = OFF, aFlag[MAXLEN];
  char szLine[MAXSIZE];
  struct residue* Seq = Sequence->residue;

  if(YFW_flag)
    strcpy(szLine, "y");
  else
  {
    printf("\nSearching for Tyr/Phe/Trp interactions? Y/N (def=Y): ");
    gets(szLine);
  }
  if (YES)
  {
    for (nIndx = 0; nIndx < Sequence->nLength; ++nIndx)
      aFlag[nIndx] = OFF; /* initialize flags */

    for (nIndx = 0; nIndx < Sequence->nLength - 4; ++nIndx)
    {
      if ('Y' == Seq[nIndx].szName[0]           && 
          (strchr("VI", Seq[nIndx+4].szName[0]) ||
          !strcmp("Nle", Seq[nIndx+4].szName)))
      /* Y XXX V/I/Nle */ 
      {
        Seq[nIndx+2].w1 *= 1.25; 
        Seq[nIndx+2].w  *= 1.25;
        Seq[nIndx+4].w1 *= 1.10;
        Seq[nIndx+4].w  *= 1.10;
        flag = flag1 = aFlag[nIndx+2] = aFlag[nIndx+4] = ON; 
        nLen = MAX(nLen, nIndx+4);
      }
      else if ('Y' == Seq[nIndx].szName[0] &&
               'L' == Seq[nIndx+4].szName[0])
      /* Y XXX L */
      {
        Seq[nIndx+2].w1 *= 1.46;
        Seq[nIndx+2].w  *= 1.46;
        flag = flag1 = aFlag[nIndx+2] = ON; 
        nLen = MAX(nLen, nIndx+2);
      }
      else if ('L' == Seq[nIndx].szName[0] &&
               strchr("YW", Seq[nIndx+4].szName[0]))
      /* L XXX Y/W */
      {
        Seq[nIndx].w1   *= 1.15;  
        Seq[nIndx].w    *= 1.15;
        Seq[nIndx+2].w1 *= 1.3; 
        Seq[nIndx+2].w  *= 1.3;
        flag = flag1 = aFlag[nIndx] = aFlag[nIndx+2] = ON;
        nLen = MAX(nLen, nIndx+2);
      }
      else if (strchr("VI", Seq[nIndx].szName[0]) &&
               strchr("YW", Seq[nIndx+4].szName[0]))
      /* V/I XXX Y/W */
      {
        Seq[nIndx+2].w1 *= 1.16; 
        Seq[nIndx+2].w  *= 1.16;
        flag = flag1 = aFlag[nIndx+2] = ON; 
        nLen = MAX(nLen, nIndx+2);
      }
      else if ('Y' == Seq[nIndx].szName[0]            &&
               (strchr("FWM", Seq[nIndx+4].szName[0]) ||
               !strcmp("Nle", Seq[nIndx+4].szName)    ||
               'H'  == Seq[nIndx+4].szName[0]         && 
               '\0' == Seq[nIndx+4].cCharge))
      /* Y XXX F/W/M/Nle/H0 */
      {
        Seq[nIndx+2].w1 *= 1.18;
        Seq[nIndx+2].w  *= 1.18;
        flag = flag1 = aFlag[nIndx+2] = ON;
        nLen = MAX(nLen, nIndx+2);
      }
      else if (strchr("FW", Seq[nIndx].szName[0]) &&
               ('L'  == Seq[nIndx+4].szName[0] ||
                'H'  == Seq[nIndx+4].szName[0] &&
                '\0' == Seq[nIndx+4].cCharge))
      /* F/W XXX L/H0 */
      {
        Seq[nIndx+2].w1 *= 1.2; 
        Seq[nIndx+2].w  *= 1.2;
        flag = flag1 = aFlag[nIndx+2] = ON;
        nLen = MAX(nLen, nIndx+2);
      }
      else if (('M' == Seq[nIndx].szName[0]       ||
               !strcmp("Nle", Seq[nIndx].szName)) &&
               strchr("FY", Seq[nIndx+4].szName[0]))
      /* M/Nle XXX F/Y */
      {
        Seq[nIndx+2].w1 *= 1.04; 
        Seq[nIndx+2].w  *= 1.04;
        flag = flag1 = aFlag[nIndx+2] = ON;
        nLen = MAX(nLen, nIndx+2);
      }
      else if ('F' == Seq[nIndx].szName[0] &&
               'M' == Seq[nIndx+4].szName[0])
      {
      /* F XXX M */
        Seq[nIndx+2].w1 *= 1.4;
        Seq[nIndx+2].w  *= 1.4;
        flag = flag1 = aFlag[nIndx+2] = ON;
        nLen = MAX(nLen, nIndx+2);
      } 
      if (ON == flag1)
      {
        PRINT_HEADER("A Try interaction is found between:", nIndx, nIndx+4);
        flag1 = OFF;
      }
    }
    if (ON == flag)
    {
      PrintParam(Seq, 0, nLen, aFlag);
      return 1;
    }
    else printf(" No Tyr interactions are found and accepted!\n");
   }

  return 0;
}
