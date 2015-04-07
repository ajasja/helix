#include "helix.h" 
#include <string.h>

#define PRINT(nCharge1, nCharge2) (printf(" %s(%d)\tand %s(%d)\n",\
                                   Seq[nCharge1].szName, nCharge1 + 1,\
                                   Seq[nCharge2].szName, nCharge2 + 1))

int FindCoul (struct sequence *Sequence, bool Coul_flag)
{
/********************************************************************
 * Description: this subroutine searches for Coulombic interactions 
 * according to the following patterns:
 * + - or - +
 * + X - or - X +
 * + XX +
 * - XX -
 * + XXX +
 * + X N/G X -
 * + X P X - 
 * + XX -
 * - XX +
 * - X N/G X +
 * - X P X +
 * + XXXXX +
 * - XXXXX -
 * R+/D- XXXXXX D-/R+
 * E- XXXXXXX R+
 * it returns 1 for successful search and acception, 0 for rejection
 ********************************************************************/

  int nIndx = 0, n2ndCharge = 0;
  int nLen = 0, nLength = Sequence->nLength;
  float g_w[MAXLEN], g_w1[MAXLEN];
  int flag = OFF, 
      aFlag1[MAXLEN], aFlag2[MAXLEN], aFlag3[MAXLEN], aFlag4[MAXLEN], 
      aFlag5[MAXLEN], aFlag6[MAXLEN], aFlag7[MAXLEN];
  char szLine[MAXSIZE];
  struct residue* Seq = Sequence->residue;

  if(Coul_flag)
    strcpy(szLine, "y");
  else
  {
    printf("\nSearch for coulombic interactions? Y/N (def=Y): "); 
    gets(szLine);
  }
  if (YES)
  {
    for (nIndx = 1; nIndx < nLength; ++nIndx) 
    {
      g_w[nIndx] = g_w1[nIndx] = 1.0;
      aFlag1[nIndx] = aFlag2[nIndx] = aFlag3[nIndx] = aFlag4[nIndx] =
      aFlag5[nIndx] = aFlag6[nIndx] = aFlag7[nIndx] = OFF;
    }
    for (nIndx = 1; nIndx < nLength-2; ++nIndx)
    {
      if ('+' * '-' == Seq[nIndx].cCharge * Seq[nIndx+1].cCharge)
    /* + - or - + */
      {
        g_w1[nIndx]   *= .95;
        g_w1[nIndx+1] *= .95;
        g_w[nIndx]    *= .95;
        g_w[nIndx+1]  *= .95;
        flag = aFlag1[nIndx] = ON; 
        nLen = MAX(nLen, nIndx+1);
      }
      else if ('+' * '-' == Seq[nIndx].cCharge * Seq[nIndx+2].cCharge)
    /* + X - or - X + */
      {
        g_w1[nIndx+1] *= .92;
        g_w[nIndx+1]  *= .92;
        flag = aFlag2[nIndx] = ON; 
        nLen = MAX(nLen, nIndx + 1);
      }
      if (nIndx < nLength - 3 && 
          '+' == Seq[nIndx].cCharge &&
          '+' == Seq[nIndx+3].cCharge)
    /* + XX + */
      { 
        g_w1[nIndx+1] *= .85; 
        g_w1[nIndx+2] *= .85;
        g_w[nIndx+1]  *= .95; 
        g_w[nIndx+2]  *= .95;
        flag = aFlag3[nIndx] = ON; 
        nLen = MAX(nLen, nIndx + 2);
      }
      else if (nIndx < nLength - 3 && 
               '-' == Seq[nIndx].cCharge &&
               '-' == Seq[nIndx+3].cCharge) 
    /* - XX - */
      {
        g_w1[nIndx+1] *= .95; 
        g_w1[nIndx+2] *= .95;
        g_w[nIndx+1]  *= .85;  
        g_w[nIndx+2]  *= .85;
        flag = aFlag3[nIndx] = ON; 
        nLen = MAX(nLen, nIndx + 2);
      }
      if (nIndx < nLength - 4 && 
          '+' == Seq[nIndx].cCharge &&
          '+' == Seq[nIndx+4].cCharge)
    /* + XXX + */
      {
        g_w1[nIndx+2] *= .6; 
        g_w[nIndx+2]  *= .78;
        flag = aFlag4[nIndx] = ON; 
        nLen = MAX(nLen, nIndx+2);
      }
      else if (nIndx < nLength - 4 && 
               '-' == Seq[nIndx].cCharge &&
               '-' == Seq[nIndx+4].cCharge)
    /* - XXX - */
      {
        g_w1[nIndx+2] *= .82; 
        g_w[nIndx+2]  *= .66;
        flag = aFlag4[nIndx] = ON; 
        nLen = MAX(nLen, nIndx + 2);
      }
      else if (nIndx < nLength - 4 && 
               '+' == Seq[nIndx].cCharge &&
               '-' == Seq[nIndx+4].cCharge)
      {
        if (strchr("NG", Seq[nIndx+2].szName[0])) /* + X N/G X - */
        {
          g_w1[nIndx+2] += .02; 
          g_w[nIndx+2]  += .1;
        }
        else if ('P' == Seq[nIndx+2].szName[0]) /* + X P X - */
          g_w[nIndx+2] += .03;

        g_w1[nIndx]   *= 1.06;  
        g_w1[nIndx+1] *= 1.08; 
        g_w1[nIndx+2] *= 1.2;
        g_w[nIndx]    *= 1.08;   
        g_w[nIndx+1]  *= 1.08;  
        g_w[nIndx+2]  *= 1.2;
        g_w[nIndx+3]  *= 1.09; 
        g_w[nIndx+4]  *= 1.07;
        flag = aFlag4[nIndx] = ON; 
        nLen = MAX(nLen, nIndx + 4);
      }
      else if (nIndx < nLength - 4 && 
               '+' == Seq[nIndx].cCharge &&
               '-' == Seq[nIndx+3].cCharge)
    /* + XX - */
      {
        g_w1[nIndx+1] *= 1.06; 
        g_w1[nIndx+2] *= 1.1; 
        g_w[nIndx]    *= 1.08;
        g_w[nIndx+3]  *= 1.08; 
        g_w[nIndx+4]  *= 1.06;
        flag = aFlag3[nIndx] = ON; 
        nLen = MAX(nLen, nIndx + 4);
      }
      else if (nIndx < nLength - 4 && 
               '-' == Seq[nIndx].cCharge &&
               '+' == Seq[nIndx+3].cCharge)
    /* - XX + */
      {
        g_w1[nIndx]   *= 1.08;  
        g_w1[nIndx+1] *= 1.1; 
        g_w1[nIndx+2] *= 1.08;
        g_w[nIndx+2]  *= 1.04; 
        g_w[nIndx+3]  *= 1.08; 
        g_w[nIndx+4]  *= 1.06;
        flag = aFlag3[nIndx] = ON; 
        nLen = MAX (nLen, nIndx + 4);
      }
      if (nIndx < nLength - 6 && 
          '-' == Seq[nIndx+1].cCharge &&
          '+' == Seq[nIndx+5].cCharge)
      {
        if (strchr("NG", Seq[nIndx+3].szName[0])) /* - X N/G X + */
        {
          g_w1[nIndx+3] += .1; 
          g_w[nIndx+3]  += .08;
        }
        else if ('P' == Seq[nIndx+3].szName[0])   /* - X P X + */
        {
          g_w1[nIndx+3] += .02; 
          g_w[nIndx+3]  += .03;
        }
        g_w1[nIndx]   *= 1.08;  
        g_w1[nIndx+1] *= 1.08; 
        g_w1[nIndx+2] *= 1.09;
        g_w1[nIndx+3] *= 1.21; 
        g_w1[nIndx+4] *= 1.08;
        g_w[nIndx+2]  *= 1.08;  
        g_w[nIndx+3]  *= 1.21;
        g_w[nIndx+4]  *= 1.09;
        g_w[nIndx+5]  *= 1.08;
        g_w[nIndx+6]  *= 1.08;
        flag = aFlag4[nIndx+1] = ON; 
        nLen = MAX(nLen, nIndx + 6);
      }
      else if (nIndx < nLength - 6 && 
               '+' == Seq[nIndx].cCharge &&
               '+' == Seq[nIndx+6].cCharge)
    /* + XXXXX + */
      {
        g_w1[nIndx]   *= .9;
        g_w1[nIndx+3] *= .9;
        g_w[nIndx+3]  *= .9; 
        g_w[nIndx+6]  *= .94;
        flag = aFlag5[nIndx] = ON; 
        nLen = MAX(nLen, nIndx + 6);
      }
      else if (nIndx < nLength - 6 && 
               '-' == Seq[nIndx].cCharge &&
               '-' == Seq[nIndx+6].cCharge)
    /* - XXXXX - */
      {
        g_w1[nIndx]   *= .94; 
        g_w1[nIndx+3] *= .9;
        g_w[nIndx+3]  *= .9;
        g_w[nIndx+6]  *= .9;
        flag = aFlag5[nIndx] = ON;
        nLen = MAX(nLen, nIndx + 6);
      }
      if (nIndx < nLength - 7           && 
          'R' == Seq[nIndx].szName[0]   && 
          '+' == Seq[nIndx].cCharge     && 
          'D' == Seq[nIndx+7].szName[0] && 
          '-' == Seq[nIndx+7].cCharge)
    /* R+ XXXXXX D- */
      {
        g_w[nIndx+2] *= 1.1; 
        g_w[nIndx+5] *= 1.2;
        flag = aFlag6[nIndx] = ON; 
        nLen = MAX(nLen, nIndx + 5);
      }
      else if (nIndx < nLength - 7           && 
               'D' == Seq[nIndx].szName[0]   &&
               '-' == Seq[nIndx].cCharge     &&
               'R' == Seq[nIndx+7].szName[0] && 
               '+' == Seq[nIndx+7].cCharge)
    /* D- XXXXXX R+ */
      {
        g_w1[nIndx+1] *= 1.1; 
        g_w1[nIndx+2] *= 1.1;
        g_w[nIndx+3]  *= 1.1;
        g_w[nIndx+5]  *= 1.2;
        flag = aFlag6[nIndx] = ON; 
        nLen = MAX(nLen, nIndx + 5);
      }
      if (nIndx < nLength - 8           && 
          'E' == Seq[nIndx].szName[0]   &&
          '-' == Seq[nIndx].cCharge     && 
          'R' == Seq[nIndx+8].szName[0] && 
          '+' == Seq[nIndx+8].cCharge)
    /* E- XXXXXXX R+ */
      {
        g_w1[nIndx+2] *= 1.3; 
        g_w[nIndx+2]  *= 1.1;
        g_w[nIndx+4]  *= 1.2; 
        g_w[nIndx+6]  *= 1.26;
        flag = aFlag7[nIndx] = ON;
        nLen = MAX(nLen, nIndx + 6);
      }
    }

    if (ON == flag)
    {
      printf("Columbic interaction is found between:\n");
      for (nIndx = 1; nIndx <= nLen; ++nIndx)
      {
        if (ON == aFlag1[nIndx])
          PRINT(nIndx, nIndx + 1);
        if (ON == aFlag2[nIndx])
          PRINT(nIndx, nIndx + 2);
        if (ON == aFlag3[nIndx])
          PRINT(nIndx, nIndx + 3);
        if (ON == aFlag4[nIndx])
          PRINT(nIndx, nIndx + 4);
        if (ON == aFlag5[nIndx])
          PRINT(nIndx, nIndx + 6);
        if (ON == aFlag6[nIndx])
          PRINT(nIndx, nIndx + 7);
        if (ON == aFlag7[nIndx])
          PRINT(nIndx, nIndx + 8);
      }
      printf("The parameters after adjustment are:\n");
      printf("\t   w'\t     w\t       vN\t vC\t   N\t     C\n");
      for (nIndx = 1; nIndx <= nLen; ++nIndx)
      {
        g_w[nIndx]  = MIN(g_w[nIndx], 1.35);  
        g_w1[nIndx] = MIN(g_w1[nIndx], 1.35);
        if (g_w1[nIndx] != 1.0 || g_w[nIndx] != 1.0)
        {
          Seq[nIndx].w1 *= g_w1[nIndx]; 
          Seq[nIndx].w  *= g_w[nIndx];
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
    else printf(" No coulombic interactions are found!\n\n");
  }

  return 0;
}
