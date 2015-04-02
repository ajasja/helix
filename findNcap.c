#include "helix.h" 
#include <string.h>

int FindNcap (struct sequence *Sequence)
{
/************************************************************************
 * Description: this subroutine searches for N-capping boxes according to 
 * the following patterns then it makes proper adjustments on
 * S/T XX E-/D-/Q
 * N/D- XX E-/D-/Q
 * E- XX E-/D-/Q
 * G XX E-/D-/Q
 * it returns 2 for sucessful search, 1 for unsuccessful search,
 * 0 for aborted search
 ************************************************************************/

  int nIndx = 0, nKcap = 0, nLen = 0;
  int aFlag[MAXLEN];
  char szLine[MAXSIZE];
  char szStr1[] = " An N-capping box is found between: "; 
  char szStr2[] = "Adjust parameters to reflect identified interaction? Y/N (def=Y): "; 
  struct residue *Seq=Sequence->residue;

  printf("\nSearching for capping boxes? Y/N (def=Y): "); 
  gets(szLine);
  if (YES)
  {
    for (nIndx = 0; nIndx < Sequence->nLength; ++nIndx)
      aFlag[nIndx] = OFF; /* initialize flags */

    nIndx = 1;
    while (nIndx < Sequence->nLength - 4)
    {
      if (strchr("ST", Seq[nIndx].szName[0]) && strncmp(Seq[nIndx].szName, "Suc", 3))
      {
        if (strchr("ED", Seq[nIndx+3].szName[0]) && '-' == Seq[nIndx+3].cCharge)
        /* S/T XX E-/D- */
        {
          PRINT_HEADER(szStr1, nIndx, nIndx + 3);
          printf(szStr2);
          gets(szLine);
          if (YES)
          {
            Seq[nIndx-1].w1 *= 0.6;   
            Seq[nIndx-1].vN *= 0.7;
            Seq[nIndx-1].vC *= 1.1;  
            Seq[nIndx-1].N  *= 0.8;
            Seq[nIndx].w1 *= 0.4;
            Seq[nIndx].w  *= 0.5;
            Seq[nIndx].vC   *= 1.1;    
            Seq[nIndx].N    *= 1.2; 
            Seq[nIndx+1].vN *= 1.15; 
            Seq[nIndx+2].w1 *= 1.15;
            Seq[nIndx+3].w1 *= 1.2;
            aFlag[nIndx-1] = aFlag[nIndx] = aFlag[nIndx+1] = 
            aFlag[nIndx+2] = aFlag[nIndx+3] = ON;
            if (nIndx >= 2) 
            {
              Seq[nIndx-2].w1 *= .7; 
              aFlag[nIndx-2] = ON;
            }
            nKcap = nIndx; 
            nLen = MAX(nLen, nIndx + 3);
          }
        }

        else if ('Q' == Seq[nIndx+3].szName[0])
        /* S/T XX Q */
        {
          PRINT_HEADER(szStr1, nIndx, nIndx + 3);
          printf(szStr2); 
          gets(szLine);
          if (YES)
          {
            Seq[nIndx-1].w1 *= 0.7;   
            Seq[nIndx-1].vN *= 0.8;  
            Seq[nIndx-1].N  *= 0.8;
            Seq[nIndx].w1   *= 0.5;     
            Seq[nIndx].w    *= 0.7;     
            Seq[nIndx].N    *= 1.06;
            Seq[nIndx+2].w1 *= 1.06; 
            Seq[nIndx+3].w1 *= 1.05;
            aFlag[nIndx-1] = aFlag[nIndx] = 
            aFlag[nIndx+2] = aFlag[nIndx+3] = ON;
            if (nIndx - 2 >= 0) 
            {
              Seq[nIndx-2].w1 *= .8; 
              aFlag[nIndx-2] = ON;
            }
            nKcap = nIndx; 
            nLen  = MAX(nLen, nIndx + 3);
          }
        }
      }

      else if ('N' == Seq[nIndx].szName[0] || 
               ('D' == Seq[nIndx].szName[0] && '-' == Seq[nIndx].cCharge))
      {
        if (strchr("ED", Seq[nIndx+3].szName[0]) && '-' == Seq[nIndx+3].cCharge)
        /* N/D- XX E-/D- */
        {
          PRINT_HEADER(szStr1, nIndx, nIndx + 3);
          printf(szStr2); 
          gets(szLine);
          if (YES)
          {
            Seq[nIndx-1].w1 *= 0.7;   
            Seq[nIndx-1].vN *= 0.7;   
            Seq[nIndx-1].N  *= 0.7;
            Seq[nIndx].w1   *= 0.5;     
            Seq[nIndx].w    *= 0.6;      
            Seq[nIndx].N    *= 1.15;
            Seq[nIndx+1].vN *= 1.05; 
            Seq[nIndx+2].w1 *= 1.06; 
            Seq[nIndx+3].w1 *= 1.12;
            aFlag[nIndx-1] = aFlag[nIndx] = aFlag[nIndx+1] = 
            aFlag[nIndx+2] = aFlag[nIndx+3] = ON;
            if (nIndx - 2 >= 0) 
            {
              Seq[nIndx-2].w1 *= .8; 
              aFlag[nIndx-2] = ON;
            }
            nKcap = nIndx; 
            nLen  = MAX(nLen, nIndx + 3);
          }
        }

        else if ('Q' == Seq[nIndx+3].szName[0])
        /* N/D- XX Q */
        {
          PRINT_HEADER(szStr1, nIndx, nIndx + 3);
          printf(szStr2); 
          gets(szLine);
          if (YES)
          {
            Seq[nIndx-1].w1 *= 0.7; 
            Seq[nIndx-1].vN *= 0.8; 
            Seq[nIndx-1].N  *= 0.8;
            Seq[nIndx].w1   *= 0.7;   
            Seq[nIndx].w    *= 0.8;    
            Seq[nIndx+3].w1 *= 1.05;
            aFlag[nIndx-1] = aFlag[nIndx] = aFlag[nIndx+3] = ON;
            if (nIndx - 2 >= 0) 
            {
              Seq[nIndx-2].w1 *= 0.8; 
              aFlag[nIndx-2] = ON;
            }
            nKcap = nIndx; 
            nLen  = MAX(nLen, nIndx + 3);
          }
        }
      }
      else if ('E' == Seq[nIndx].szName[0] && '-' == Seq[nIndx].cCharge)
      {
        if (strchr("ED", Seq[nIndx+3].szName[0]) && '-' == Seq[nIndx+3].cCharge)
        /* E- XX E-/D- */
        {
          PRINT_HEADER(szStr1, nIndx, nIndx + 3);
          printf(szStr2); 
          gets(szLine);
          if (YES)
          {
            Seq[nIndx-1].w1 *= 0.7;  
            Seq[nIndx-1].N  *= 0.8;
            Seq[nIndx].w1   *= 0.5;    
            Seq[nIndx].w    *= 0.6;     
            Seq[nIndx].N    *= 1.3;    
            Seq[nIndx+1].vN *= 1.1; 
            Seq[nIndx+2].w1 *= 1.1; 
            Seq[nIndx+3].w1 *= 1.1;
            aFlag[nIndx-1] = aFlag[nIndx] = aFlag[nIndx+1] = 
            aFlag[nIndx+2] = aFlag[nIndx+3] = ON;
            if (nIndx - 2 >= 0) 
            {
              Seq[nIndx-2].w1 *= .8; 
              aFlag[nIndx-2] = ON;
            }
            nKcap = nIndx; 
            nLen = MAX(nLen, nIndx + 3);
          }
        }

        else if ('Q' == Seq[nIndx+3].szName[0])
        /* E- XX Q */
        {
          PRINT_HEADER(szStr1, nIndx, nIndx + 3);
          printf(szStr2); 
          gets(szLine);
          if (YES)
          {
            Seq[nIndx-1].w1 *= 0.7; 
            Seq[nIndx-1].w  *= 0.8; 
            Seq[nIndx-1].N  *= 0.8;
            Seq[nIndx].w1   *= 0.5; 
            Seq[nIndx].w    *= 0.8;     
            Seq[nIndx].N    *= 1.2; 
            Seq[nIndx+2].w1 *= 1.08;
            aFlag[nIndx-1] = aFlag[nIndx] = aFlag[nIndx+2] = ON;
            if (nIndx - 2 >= 0) 
            {
              Seq[nIndx-2].w1 *= 0.8; 
              aFlag[nIndx-2] = ON;
            } 
            nKcap = nIndx;
            nLen = MAX(nLen, nIndx + 2);
          }
        }
      }

      else if ('G' == Seq[nIndx].szName[0] && 
               strchr("ED", Seq[nIndx+3].szName[0]) &&
               '-' == Seq[nIndx+3].cCharge)
      /* G XX E-/D- */
      {
        PRINT_HEADER(szStr1, nIndx, nIndx + 3);
        printf(szStr2); 
        gets(szLine);
        if (YES)
        {
          Seq[nIndx-1].w1 *= 0.7;
          Seq[nIndx].vN   *= 1.05; 
          Seq[nIndx].N    *= 1.05;
          Seq[nIndx+3].w1 *= 1.05;
          aFlag[nIndx-1] = aFlag[nIndx] = aFlag[nIndx+3] = ON;
          nKcap = nIndx; 
          nLen = MAX(nLen, nIndx + 3);
        }
      }

      else if ('G' == Seq[nIndx].szName[0] && 'Q' == Seq[nIndx+3].szName[0])
      /* G XX Q */
      {
        PRINT_HEADER(szStr1, nIndx, nIndx + 3);
        printf(szStr2); 
        gets(szLine);
        if (YES)
        {
          Seq[nIndx+2].w1 *= .5;
          aFlag[nIndx] = aFlag[nIndx+2] = ON;
          nKcap = nIndx; 
          nLen = MAX(nLen, nIndx + 2);
        }
      }
      ++nIndx;
    }

    if (nKcap > 0) /* k is the initial cap position */
    {
      Seq[nKcap+2].w1 *= 3.0;
      if ('-' == Seq[nKcap+1].cCharge || '-' == Seq[nKcap+2].cCharge) 
        Seq[nKcap+1].vN *= 2.0; /* 1.44 * 1.4 */
      else 
        Seq[nKcap+1].vN *= 1.44;
      if ('-' == Seq[nKcap].cCharge   ||
          '-' == Seq[nKcap+1].cCharge || 
          '-' == Seq[nKcap+2].cCharge || 
          '-' == Seq[nKcap+3].cCharge)
      {
        Seq[nKcap+3].w1 *= 1.5; 
        aFlag[nKcap+3] = ON;
      }
      if (strchr("YFLVIM", Seq[nKcap-1].szName[0]) ||
         ('H' == Seq[nKcap-1].szName[0] && 0 == Seq[nKcap-1].cCharge))
      {
        Seq[nKcap].N *= 1.08; 
        aFlag[nKcap] = ON;
      }
      aFlag[nKcap+1] = aFlag[nKcap+2] = ON;

      PrintParam(Seq, 0, nLen, aFlag);
      return 2;
    }
    else 
    {
      printf(" No N-capping boxes are found and accepted!\n");
      return 1;
    }
  }
  return 0;
}
