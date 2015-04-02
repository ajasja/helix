#include "helix.h"

void PrintParam (struct residue *Seq, int nFirst, int nLast, int *aFlag)
{
/*******************************************************************
 * Description: this subroutine prints out the parameters of a residue
 * when its flag is set.
 *******************************************************************/

  int nIndx = 0;

  printf("The parameters after adjustments are:\n");
  printf("\t   w'\t     w\t       vN\t vC\t   N\t     C\n");
  for (nIndx = nFirst; nIndx <= nLast; ++nIndx)
  {
    if (ON == aFlag[nIndx])
      printf("%2d %s\t%f  %f  %f  %f  %f  %f\n", 
             nIndx+1, 
             Seq[nIndx].szName,
             Seq[nIndx].w1,
             Seq[nIndx].w, 
             Seq[nIndx].vN,
             Seq[nIndx].vC,
             Seq[nIndx].N,
             Seq[nIndx].C);
  }
}
