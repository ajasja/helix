#include "helix.h" 
#include <string.h>

int FindCend (struct sequence *Sequence)
{
/********************************************************************
 * Description: this subroutine searches for end effects at C terminus
 * charged residues between n-4 and n
 * it returns 1 for sucessful search and acceptance, 0 for rejection
 ********************************************************************/

  int flag = OFF, aFlag[MAXLEN];
  char szLine[MAXSIZE];
  char szStr1[] = "Adjust parameters to reflect identified interaction? Y/N (def=Y): ";
  char szStr2[] = "End effect is found between: ";
  struct residue *Seq = Sequence->residue;
  int nLength = Sequence->nLength;

  printf("\nSearching for C-terminal end effects? Y/N (def=Y): "); 
  gets(szLine); 
  if (YES)
  {
    aFlag[nLength-5] = aFlag[nLength-4] = aFlag[nLength-3] = 
    aFlag[nLength-2] = aFlag[nLength-1] = OFF;

    if (!strcmp(Seq[nLength-1].szName, "NH2"))
    {
      if ('+' == Seq[nLength-5].cCharge)
      {
        PRINT_HEADER(szStr2, nLength-5, nLength-1);
        printf("\n%s", szStr1);
        gets(szLine);
        if (YES)
        {
          Seq[nLength-5].vC *= 1.06;
          Seq[nLength-5].w  *= 1.14;
          flag = aFlag[nLength-5] = ON;
        }
      }
      else if ('+' == Seq[nLength-4].cCharge)
      {
        PRINT_HEADER(szStr2, nLength-4, nLength-1);
        printf("\n%s", szStr1);
        gets(szLine);
        if (YES)
        {
          Seq[nLength-4].vC *= 1.08; 
          Seq[nLength-4].w  *= 1.18;
          flag = aFlag[nLength-4] = ON;
        }
      }
      else if('+' == Seq[nLength-3].cCharge)
      {
        PRINT_HEADER(szStr2, nLength-3, nLength-1);
        printf("\n%s", szStr1);
        gets(szLine);
        if (YES)
        {
          Seq[nLength-3].vC *= 1.1; 
          Seq[nLength-3].w  *= 1.18; 
          Seq[nLength-3].C  *= 1.12;
          flag = aFlag[nLength-3] = ON;
        }
      }
      else if ('+' == Seq[nLength-2].cCharge)
      {
        PRINT_HEADER(szStr2, nLength-2, nLength-1);
        printf("\n%s", szStr1); 
        gets(szLine);
        if (YES)
        {
          Seq[nLength-2].vC *= 1.1;
          Seq[nLength-2].C  *= 1.12;
          flag = aFlag[nLength-2] = ON;
        }
      }
    }
    else if (LOW == Sequence->salt && Sequence->pH >= 4.6)
    {
      Seq[nLength-1].C  *= .8;
      Seq[nLength-2].vC *= .8;
      Seq[nLength-2].C  *= .92;
      flag = aFlag[nLength-1] = aFlag[nLength-2] = ON;

      if (Seq[nLength-1].cCharge != '\0')
      {
        printf("%s%s(%d)\tand ", szStr2, Seq[nLength-1].szName, nLength);
        printf("%s(%d)\n%s", "COO-", nLength, szStr1); 
        gets(szLine);
        if (YES)
        {
          if ('+' == Seq[nLength-1].cCharge) 
            Seq[nLength-1].C *= 1.2;
          else                             
            Seq[nLength-1].C *= 0.84;
          flag = aFlag[nLength-1] = ON;
        }
      }
      else if (Seq[nLength-2].cCharge != '\0')
      {
        printf("%s%s(%d)\tand ", szStr2, Seq[nLength-2].szName, nLength - 1);
        printf("%s(%d)\n%s", "COO-", nLength, szStr1); 
        gets(szLine);
        if (YES)   
        {
          if ('+' == Seq[nLength-2].cCharge) 
          {
            Seq[nLength-2].vC *= 1.1; 
            Seq[nLength-2].C  *= 1.1;
          }
          else
          {
            Seq[nLength-2].vC *= .84;
            Seq[nLength-2].C  *= .9;
          }
          flag = aFlag[nLength-2] = ON;
        }
      }
      else if (Seq[nLength-3].cCharge != '\0')
      {
        printf("%s%s(%d)\tand ", szStr2, Seq[nLength-3].szName, nLength - 2);
        printf("%s(%d)\n%s", "COO-", nLength, szStr1);
        gets(szLine);
        if (YES)
        {
          if ('+' == Seq[nLength-3].cCharge)
          {
            Seq[nLength-3].vC *= 1.1; 
            Seq[nLength-3].w  *= 1.18;
          }
          else
          {
            Seq[nLength-3].vC *= .9;
            Seq[nLength-3].w  *= .84;
          }
          flag = aFlag[nLength-3] = ON;
        }
      }
      else if (Seq[nLength-4].cCharge != '\0')
      {
        printf("%s%s(%d)\tand ", szStr2, Seq[nLength-4].szName, nLength - 3);
        printf("%s(%d)\n%s", "COO-", nLength, szStr1); 
        gets(szLine);
        if (YES)
        {
          if ('+' == Seq[nLength-4].cCharge) 
            Seq[nLength-4].w *= 1.15;
          else                     
            Seq[nLength-4].w *= .87;
          flag = aFlag[nLength-4] = ON;
        }
      }
      else if (Seq[nLength-5].cCharge != 0)
      {
        printf("%s%s(%d)\tand ", szStr2, Seq[nLength-5].szName, nLength - 4);
        printf("%s(%d)\n%s", "COO-", nLength, szStr1);
        gets(szLine);
        if (YES)
        {
          if ('+' == Seq[nLength-5].cCharge) 
            Seq[nLength-5].w *= 1.11;
          else                            
            Seq[nLength-5].w *= 0.92;
          flag = aFlag[nLength-5] = ON;
        }
      }
    }
    if (ON == flag)
    {
      PrintParam(Seq, nLength - 5, nLength - 1, aFlag);
      return 1;
    }
    else 
      printf(" No C-terminal end effects are found and accepted!\n");
  }

  return 0;
}
