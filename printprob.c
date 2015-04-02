#include "helix.h"
#define PRINT(format, value) fprintf(pOutFile, format, value);\
                              printf(format, value)

void PrintProb(struct sequence *Sequence)
{
/*******************************************************************
 * Description: this subroutine prints out the list of probabilities
 * of each residue at the end of the program 
 *******************************************************************/

  int nIndx = 0;
  FILE *pOutFile = 0;
  float avg = 0.0;
  struct residue *Seq=Sequence->residue;
  int nLength = Sequence->nLength;

  if (!(pOutFile = fopen("helix.out", "w")))
    return;

  PRINT("%s\n", "****** Helix v.1.5 (C) 1995-96 Andersen Group ******");
  PRINT("Sequence file: %s\n", Sequence->szSeqFile);
  PRINT("Database file: %s\n", Sequence->szDatFile);
  PRINT("pH = %.1f  ", Sequence->pH);
  if (LOW == Sequence->salt)
  {
    PRINT("%s\n", "Low Salt Condition");
  }
  else
  {
    PRINT("%s\n", "High Salt Condition"); 
  }
  PRINT("\n\t%s", "P(w')  P(w)   P(vN)  P(vC) P(w'+w)  fH"); 
  fprintf(pOutFile, "\n 1 %s\t%.3f  %.3f  %.3f  %.3f  %.3f  %.3f",
	  Seq[0].szName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  printf("\n 1 %s\t%.3f  %.3f  %.3f  %.3f  %.3f  %.3f",
	  Seq[0].szName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

  for (nIndx = 1; nIndx < nLength - 1; ++nIndx)
  {
    fprintf(pOutFile, 
            "\n%2d %s\t%.3f  %.3f  %.3f  %.3f  %.3f %.3f", 
            nIndx + 1, 
            Seq[nIndx].szName,
            Seq[nIndx].Pw1, 
            Seq[nIndx].Pw,
            Seq[nIndx].PvN,
            Seq[nIndx].PvC,
            Seq[nIndx].Pw1 + Seq[nIndx].Pw,
            Seq[nIndx].fH);
    printf("\n%2d %s\t%.3f  %.3f  %.3f  %.3f  %.3f  %.3f", 
           nIndx + 1,
           Seq[nIndx].szName,
           Seq[nIndx].Pw1, 
           Seq[nIndx].Pw, 
           Seq[nIndx].PvN, 
           Seq[nIndx].PvC,
           Seq[nIndx].Pw1 + Seq[nIndx].Pw,
           Seq[nIndx].fH);
    avg += Seq[nIndx].Pw1 + Seq[nIndx].Pw;
  }
  fprintf(pOutFile, "\n%2d %s\t%.3f  %.3f  %.3f  %.3f  %.3f %.3f",
          nLength, Seq[nLength-1].szName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  printf("\n%2d %s\t%.3f  %.3f  %.3f  %.3f  %.3f  %.3f", 
          nLength, Seq[nLength-1].szName, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  fprintf(pOutFile, "\n <P(w'+w)>(n-4) = %.4f  <P(w'+w)>(n-3) = %.4f",
          avg/(nLength-4), avg/(nLength-3));
  printf("\n <P(w'+w)>(n-4) = %.4f  <P(w'+w)>(n-3) = %.4f",
         avg/(nLength-4), avg/(nLength-3));
  PRINT("\n\n%s\n", "*************** End of the Output ****************");
  fclose(pOutFile);
}
