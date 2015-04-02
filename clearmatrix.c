#include "matrix.h"

struct matrix ClearMatrix(struct matrix W, float type)
{
/*****************************************************************
 * Description: this subroutine sets every element in the matrix
 * except the designated type to 0, which is equivalent to 
 * (type)*dW/d(type)
 * it returns a copy of the cleared matrix
 *****************************************************************/

  int nRow = 0, nCol = 0;

  for (nCol = 0; nCol < W.nCol; ++nCol)
    for (nRow = 0; nRow < W.nRow; ++nRow)
      if (W.Wi[nRow][nCol] != 0.0 && W.Wi[nRow][nCol] != type) 
        W.Wi[nRow][nCol] = 0.0;

  return W;
}
