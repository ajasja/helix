#include "matrix.h"

float MtplyMatrix(struct matrix W, /* matrix of the next residue */
                  int nIndx,       /* index along the sequence   */
                  int nResidue,    /* index of the residue       */
                  float type,      /* type of probability        */
                  struct sequence *Sequence)
{
/*************************************************************
 * Description: this subroutine performs matrix multiplication
 * (1 by m * m by n = 1 by n), the result may be used to set the
 * initial matrix during the next calculation
 * it returns the 1 by 1 matrix (a float) in the end
 *************************************************************/

  int nRow = 0, nCol = 0;
  float aSum[MAXDIM] = {0.0}; /* array to hold the 1 by n matrix */
  struct matrix* pW = 0;
  /* static structures to save intermediate       */
  /* results for different types of probabilities */
  static struct matrix W0, Ww1, Ww, WvN, WvC, Wu;
  static struct matrix Winit = {1, 3, {u, u, u}};
  struct residue *Seq = Sequence->residue;

  if (type == 0.0)
    pW = &W0;
  else if (type == w1_m) 
    pW = &Ww1;
  else if (type == w_m)  
    pW = &Ww;
  else if (type == vN_m) 
    pW = &WvN;
  else if (type == vC_m) 
    pW = &WvC;
  else if (type == u)
    pW = &Wu;

  if (1 == nIndx || nIndx == nResidue - 1)
    *pW = Winit;         /* set the matrix to Winit */
  if (nIndx <= nResidue - 1 && type != u)
    return pW->Wi[0][0]; /* same result as Wu */

  for (nCol = 0; nCol < W.nCol; ++nCol)
    for (nRow = 0; nRow < W.nRow; ++nRow)
      if (W.Wi[nRow][nCol] == w1_m)
        aSum[nCol] += pW->Wi[0][nRow] * Seq[nIndx].w1;
      else if(W.Wi[nRow][nCol] == w_m)
        aSum[nCol] += pW->Wi[0][nRow] * Seq[nIndx].w;
      else if (W.Wi[nRow][nCol] == vN_m)
        aSum[nCol] += pW->Wi[0][nRow] * Seq[nIndx-1].N * Seq[nIndx].vN;
      else if (W.Wi[nRow][nCol] == vC_m)
        aSum[nCol] += pW->Wi[0][nRow] * Seq[nIndx+1].C * Seq[nIndx].vC;
      else if (W.Wi[nRow][nCol] == u)
        aSum[nCol] += pW->Wi[0][nRow];

  for (nCol = 0; nCol < W.nCol; ++nCol)
    pW->Wi[0][nCol] = aSum[nCol];
  pW->nCol = W.nCol;

  if (type == u && nIndx == nResidue - 1)
    Winit = Wu;  /* set Winit for the next residue */

  return pW->Wi[0][0];
}
