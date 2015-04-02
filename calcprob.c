#include "matrix.h"

int CalcProb (struct sequence *Sequence)
{
/*******************************************************************
 * Description: this subroutine calculates the configurational 
 * probabilities of each residue being in w1, w, vN or vC state
 * for example, Pw = d(lnZ)/d(lnw) = (w/Z0)*(dZw/dw)
 * where Z = W[0].W[1]...W[12]
 * indices 1 - 8 correspond to matrices W[0] - W[7]
 * indices 9 - n-5 correspond to matrices W[7]
 * indices n-4 - n correspond to matrices W[8] - W[12]
 * it returns 1 for success, 0 for failure
 *******************************************************************/

  int nIndx = 0, nMatrix = 0, nResidue = 0;
  float Z0, aZu[MAXLEN], aZw1[MAXLEN], aZw[MAXLEN], aZvN[MAXLEN], aZvC[MAXLEN];
  struct residue *Seq = Sequence->residue;
  int nLength = Sequence->nLength;

  printf("Starting matrix calculation: ");
  fflush(stdout);

  for (nIndx = 1, nMatrix = 0; nIndx < nLength; ++nIndx)
  {
    nMatrix = (nIndx < 7 || nIndx > nLength - 6) ? nMatrix + 1 : 7; 
    Z0 = MtplyMatrix(W[nMatrix], nIndx, 0, 0.0, Sequence);
  }

  for (nResidue = 1; nResidue < nLength - 1; ++nResidue)
  {
    nIndx = (1 == nResidue) ? 1 : nResidue - 1;
    while (nIndx < nLength)
    {
      if (nIndx <= 7)
        nMatrix = nIndx;
      else if (nIndx > 7 && nIndx < nLength - 5)
        nMatrix = 7;
      else
        nMatrix = 12 - (nLength - 1 - nIndx);
      if (nIndx == nResidue)
      {
        MtplyMatrix(ClearMatrix(W[nMatrix], u),    nIndx, nResidue, u,    Sequence);
        MtplyMatrix(ClearMatrix(W[nMatrix], w1_m), nIndx, nResidue, w1_m, Sequence);
        MtplyMatrix(ClearMatrix(W[nMatrix], w_m),  nIndx, nResidue, w_m,  Sequence);
        MtplyMatrix(ClearMatrix(W[nMatrix], vN_m), nIndx, nResidue, vN_m, Sequence);
        MtplyMatrix(ClearMatrix(W[nMatrix], vC_m), nIndx, nResidue, vC_m, Sequence);
      }
      else
      {
        aZu[nResidue]  = MtplyMatrix(W[nMatrix], nIndx, nResidue, u,    Sequence);
        aZw1[nResidue] = MtplyMatrix(W[nMatrix], nIndx, nResidue, w1_m, Sequence);
        aZw[nResidue]  = MtplyMatrix(W[nMatrix], nIndx, nResidue, w_m,  Sequence);
        aZvN[nResidue] = MtplyMatrix(W[nMatrix], nIndx, nResidue, vN_m, Sequence);
        aZvC[nResidue] = MtplyMatrix(W[nMatrix], nIndx, nResidue, vC_m, Sequence);
      }
      ++nIndx;
    }
    /* calculate fractional helicity and probabilities of each */
    /* residue being in different states                       */
    Seq[nResidue].fH  = 1.0 - aZu[nResidue] / Z0;
    Seq[nResidue].Pw1 = aZw1[nResidue] / Z0;
    Seq[nResidue].Pw  = aZw[nResidue]  / Z0;
    Seq[nResidue].PvN = aZvN[nResidue] / Z0;
    Seq[nResidue].PvC = aZvC[nResidue] / Z0;
    printf(". ");
    fflush(stdout);
  }
  printf("done!\n\n");

  return 1;
}
