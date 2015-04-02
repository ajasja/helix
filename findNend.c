#include "helix.h"
#include <string.h>
#include <math.h>

int FindNend (struct sequence *Sequence)
{
/********************************************************************
 * Description: this subroutine searches for end effects at N terminus
 * charged residues between 1 and 5 
 * it returns 1 for successful search and acception, 0 for rejection
 ********************************************************************/

  int nIndx = 0, nLen = 0;
  int flag1 = OFF, flag2 = OFF, aFlag[MAXLEN];
  char szLine[MAXSIZE];
  char szStr1[] = "\nAdjust parameters to reflect identified interaction? Y/N (def=Y): ";
  char szStr2[] = "\nEnd effect is found between: ";
  struct residue *Seq = Sequence->residue;

  printf("\nSearching for N-terminal end effects? Y/N (def=Y): "); 
  gets(szLine);
  if (YES)
  {
    for (nIndx = 0; nIndx < 5; ++nIndx)
      aFlag[nIndx] = OFF; /* initialize flags */

    if( !strncmp(Seq[0].szName, "Suc", 3) && '-' == Seq[0].cCharge)
    {
      printf("%s%s(%d)\n", "End effect is found at ", Seq[0].szName, 1);
      printf(szStr1); 
      gets(szLine);
      if (YES)
      {
        Seq[1].N  *= 1.22; 
        Seq[2].N  *= 1.1;
        Seq[2].w1 *= 1.2; 
        Seq[3].w1 *= 1.1;
        flag1 = aFlag[1] = aFlag[2] = aFlag[3] = ON; 
        nLen = MAX(nLen, 3);
      }
      if ('+' == Seq[1].cCharge && '+' == Seq[2].cCharge)
      {
        printf("%s%s(%d)\tand %s(%d),%s(%d)\n",
                szStr2,
                Seq[0].szName, 1,
                Seq[1].szName, 2,
                Seq[2].szName, 3);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if (LOW == Sequence->salt) 
          {
            Seq[2].vN *= .75;
            Seq[2].w1 *= .75;
          }
          else                    
          {
            Seq[2].vN *= sqrt(.75); 
            Seq[2].w1 *= sqrt(.75);
          }
          flag1 = aFlag[2]=ON; 
          nLen = MAX(nLen, 2);
        }
      }
      else if ('+' == Seq[1].cCharge)
      {
        PRINT_HEADER(szStr2, 0, 1);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if (LOW == Sequence->salt)
          {
            Seq[0].N  *= .7; 
            Seq[1].N  *= .85; 
            Seq[1].vN *= .9;
          }
          else
          {
            Seq[0].N  *= sqrt(.7); 
            Seq[1].N  *= sqrt(.85);
            Seq[1].vN *= sqrt(.9);
          }
          flag1 = aFlag[0] = aFlag[1] = ON; 
          nLen = MAX(nLen, 1);
        } 
      }
      else if ('+' == Seq[2].cCharge)
      {
        PRINT_HEADER(szStr2, 0, 2);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if (LOW == Sequence->salt) 
            Seq[2].vN *= .85;
          else                    
            Seq[2].vN *= sqrt(.85);
          flag1 = aFlag[2] = ON; 
          nLen = MAX(nLen, 2);
        }
      }
    }
    else if (!strcmp(Seq[0].szName, "Ac")      || 
             !strcmp(Seq[0].szName, "sN")      ||
             !strncmp(Seq[0].szName, "Suc", 3) &&
             0 == Seq[0].cCharge)
    {
      if ('+' == Seq[1].cCharge)
      {
        PRINT_HEADER(szStr2, 0, 1);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if (LOW == Sequence->salt) 
            Seq[1].vN *= .65;
          else                    
            Seq[1].vN *= sqrt(.65);
          flag1 = aFlag[1] = ON; 
          nLen = MAX(nLen, 1);
        }
      }
      if ('+' == Seq[2].cCharge)
      {
        PRINT_HEADER(szStr2, 0, 2);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if (LOW == Sequence->salt) 
          {
            Seq[2].w1 *= .8; 
            Seq[2].vN *= .7;
          }
          else if(OFF == flag1)  
          {
            Seq[2].w1 *= sqrt(.8);
            Seq[2].vN *= sqrt(.7);
          }
          flag1 = aFlag[2] = ON; 
          nLen = MAX(nLen, 2);
        }
      }
      if ('+' == Seq[3].cCharge)
      {
        PRINT_HEADER(szStr2, 0, 3);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if (LOW == Sequence->salt) 
          {
            Seq[3].w1 *= .7; 
            Seq[3].vN *= .8;
          }
          else if (OFF == flag1)      
          {
            Seq[3].w1 *= sqrt(.7);
            Seq[3].vN *= sqrt(.8);
          }
          flag1 = aFlag[3] = ON; 
          nLen = MAX(nLen, 3);
        }
      }
      nIndx = 1;
      while (nIndx < 5 && OFF == flag2) /* only first occurance */
      {
        if ('D' == Seq[nIndx].szName[0] && '-' == Seq[nIndx].cCharge)
        {
          PRINT_HEADER(szStr2, 0, nIndx);
          printf(szStr1); 
          gets(szLine);
          if (YES)
          {
            if (LOW == Sequence->salt) 
              Seq[nIndx].N *= 1.14;
            else                    
              Seq[nIndx].N *= sqrt(1.14);
            Seq[nIndx].vN *= 1.12; 
            Seq[nIndx].w1 *= 1.25;
            flag1 = flag2 = aFlag[nIndx] = ON; 
            nLen = MAX(nLen, nIndx);
          }
        }
        else if ('E' == Seq[nIndx].szName[0] && '-' == Seq[nIndx].cCharge)
        {
          PRINT_HEADER(szStr2, 0, nIndx);
          printf(szStr1); 
          gets(szLine);
          if (YES)
          {
            if (LOW == Sequence->salt) 
              Seq[nIndx].N *= 1.26;
            else                    
              Seq[nIndx].N *= sqrt(1.26);
            Seq[nIndx].vN *= 1.08; 
            Seq[nIndx].w1 *= 1.1;
            flag1 = flag2 = aFlag[nIndx] = ON; 
            nLen = MAX(nLen, nIndx);
          }
        }
        ++nIndx;
      }
      flag2 = OFF;
    }
    else if (Sequence->pH <= 7.7)
    {
      if (Seq[0].cCharge != '-' && 
          Seq[1].cCharge != '-' && 
          Seq[2].cCharge != '-')
      {
        printf(" NH3+(1) w/o D-/E- at 1->3 recognized, ");
        printf("accept? Y/N (def=Y): "); 
        gets(szLine);
        if (YES)
        {
          Seq[2].w1 *= .96; 
          Seq[1].vN *= .78;
          Seq[0].N  *= .9;   
          Seq[2].vN *= .96;
          flag1 = aFlag[0] = aFlag[1] = aFlag[2] = ON; 
          nLen = MAX(nLen, nIndx);
        }
      }
      if ('+' == Seq[0].cCharge)
      {
        printf("%s%s\tand %s(%d)", szStr2, "NH3+(1)", Seq[0].szName, 1);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if (LOW == Sequence->salt==LOW) 
            Seq[0].N *= .68;
          else                    
            Seq[0].N *= sqrt(.68);
          flag1 = aFlag[0] = ON;
        }
      }
      if ('+' == Seq[1].cCharge)
      {
        printf("%s%s\tand %s(%d)", szStr2, "NH3+(1)", Seq[1].szName, 2);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if (LOW == Sequence->salt==LOW) 
          {
            Seq[1].N  *= .75; 
            Seq[1].vN *= .8;
          }
          else                    
          {
            Seq[1].N  *= sqrt(.75); 
            Seq[1].vN *= sqrt(.8);
          }
          flag1 = aFlag[1] = ON; 
          nLen = MAX(nLen, 1);
        }
      }
      if ('+' == Seq[2].cCharge)
      {
        printf("%s%s\tand %s(%d)", szStr2, "NH3+(1)", Seq[2].szName, 3);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if (LOW == Sequence->salt) 
          {
            Seq[2].vN *= .85;
            Seq[2].w1 *= .86;
          }
          else                    
          {
            Seq[2].vN *= sqrt(.85); 
            Seq[2].w1 *= sqrt(.86);
          }
          flag1 = aFlag[2] = ON; 
          nLen = MAX(nLen, 2);
        }
      }
      nIndx = 0; /* nIndx is the number of D-/E- found */
      if ('D' == Seq[0].szName[0] && '-' == Seq[0].cCharge)
      {
        printf("%s%s\tand %s(%d)", szStr2, "NH3+(1)", Seq[0].szName, 1);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          Seq[0].N *= 1.1;
          flag1 = aFlag[0] = ON; 
          ++nIndx;
        }
      }
      if ('-' == Seq[1].cCharge)
      {
        printf("%s%s\tand %s(%d)", szStr2, "NH3+(1)", Seq[1].szName, 2);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if ('D' == Seq[1].szName[0] && '-' == Seq[1].cCharge)
          {
            Seq[1].N  *= 1.1; 
            Seq[1].vN *= 1.08;
          }
          else if ('E' == Seq[1].szName[0] && '-' == Seq[1].cCharge)
          {
            Seq[1].N  *= 1.2;
            Seq[1].vN *= 1.1;
          }
          flag1 = aFlag[1] = ON; 
          ++nIndx; 
          nLen = MAX(nLen, 1);
        }
      }
      if ('-' == Seq[2].cCharge && nIndx < 2)
      {
        printf("%s%s\tand %s(%d)", szStr2, "NH3+(1)", Seq[2].szName, 3);
        printf(szStr1);
        gets(szLine);
        if (YES)
        {
          if ('D' == Seq[2].szName[0] && '-' == Seq[2].cCharge)
          {
            Seq[2].N  *= 1.1; 
            Seq[2].vN *= 1.08; 
            Seq[2].w1 *= 1.25;
          }
          else if ('E' == Seq[2].szName[0] && '-' == Seq[2].cCharge)
          {
            Seq[2].vN *= 1.18; 
            Seq[2].w1 *= 1.1;
          }
          flag1 = aFlag[2] = ON; 
          ++nIndx; 
          nLen = MAX(nLen, 2);
        }
      }
      if ('-' == Seq[3].cCharge && nIndx < 2)
      {
        printf("%s%s\tand %s(%d)", szStr2, "NH3+(1)", Seq[3].szName, 4);
        printf(szStr1); 
        gets(szLine);
        if (YES)
        {
          if ('D' == Seq[3].szName[0] && '-' == Seq[3].cCharge)
          {
            Seq[3].vN *= 1.08; 
            Seq[3].w1 *= 1.12;
          }
          else if ('E' == Seq[3].szName[0] && '-' == Seq[3].cCharge)
            Seq[3].w1 *= 1.12;
          flag1 = aFlag[3] = ON; 
          nLen = MAX(nLen, 3);
        }
      }
    }
    if (ON == flag1)
    {
      PrintParam (Seq, 0, nLen, aFlag);
      return 1;
    }
    else 
      printf(" No N-terminal end effects are found and accepted!\n");
  }

  return 0;
}
