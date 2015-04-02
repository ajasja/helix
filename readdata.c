#include "helix.h"
#include <string.h>

int ReadData(struct sequence *Sequence)
{
/*******************************************************************
 * Description: this subroutine gets the initial input from the
 * sequence file, parameter file and the previous output file
 * it returns 1 for successful reading, 0 for failure
 *******************************************************************/

  int nIndex = 0;
  float pH = 7.0;
  FILE *pOutFile = 0, *pSeqFile = 0, *pDatFile = 0;
  char szLine[MAXSIZE];
  struct residue *Seq = Sequence->residue;

/* getting default sequence and parameter file names as well as pH value */
  if (pOutFile = fopen("helix.out", "r"))
  {
    fgets(szLine, 80, pOutFile);
    fscanf(pOutFile, "%*s%*s%s", Sequence->szSeqFile);
    fscanf(pOutFile, "%*s%*s%s", Sequence->szDatFile);
    if (1 != fscanf(pOutFile, "%*s%*s%f", &pH))
      pH = 7.0;
    fclose(pOutFile);
  }

  printf("-->> Maximum Sequence Length: %d <<--\n", MAXLEN);
  printf("Enter the sequence file (def=%s)-> ", Sequence->szSeqFile);
  gets(szLine);
  if (szLine[0] != '\0')
    strncpy(Sequence->szSeqFile, szLine, MAXSIZE);
  if (!(pSeqFile = fopen(Sequence->szSeqFile, "r")))
  {
    printf("\aSequence file not found: %s\n", Sequence->szSeqFile);
    return 0;
  }
  fclose(pSeqFile);

  printf("Enter the database file (def=%s)-> ", Sequence->szDatFile); 
  gets(szLine);
  if (szLine[0] != '\0') 
    strncpy(Sequence->szDatFile, szLine, MAXSIZE);
  if (!(pDatFile = fopen(Sequence->szDatFile, "r")))
  {
    printf("\aDatabase file not found: %s\n", Sequence->szDatFile);
    return 0;
  }
  fclose(pDatFile);

  printf("Want to change pH (%.1f)? Y/N (def=N): ", pH);
  gets(szLine);
  if ('y' == szLine[0] || 'Y' == szLine[0])
  {
    printf("Enter the new pH value: ");
    while (scanf("%f", &pH) != 1 || pH < 0.0 || pH > 13.0)
    { 
      while (getchar() != '\n');
      printf("\aInvalid pH value, try it again: ");
    }
    while (getchar() != '\n');
  }
  Sequence->pH = pH;

  printf("High salt or low salt concentration? H/L (def=L): "); 
  gets(szLine);
  if ('h' == szLine[0] || 'H' == szLine[0]) 
    Sequence->salt = HIGH;
  else
    Sequence->salt = LOW;

  if (!ReadSeq(Sequence) || !ReadParam(Sequence))
    return 0;

  printf("\nThe sequence and initial parameters are:\n");
  printf("\t  w'\t  w\t  vN\t  vC\t  N\t  C\n");
  for (nIndex = 0; nIndex < Sequence->nLength; ++nIndex)
    printf("%2d %s\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
           nIndex + 1, 
           Seq[nIndex].szName, 
           Seq[nIndex].w1, 
           Seq[nIndex].w, 
           Seq[nIndex].vN, 
           Seq[nIndex].vC, 
           Seq[nIndex].N,
           Seq[nIndex].C);

  return 1;
}
