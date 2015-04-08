#include "helix.h"
#include <string.h>

int ReadData(struct sequence *Sequence, char* cmdseq, bool seq_flag, char* salt, float pH)
{
/*******************************************************************
 * Description: this subroutine gets the initial input from the
 * sequence file, parameter file and the previous output file,
 * OR from the command line
 * it returns 1 for successful reading, 0 for failure
 *******************************************************************/

  int nIndex = 0;
  FILE *pOutFile = 0, *pSeqFile = 0, *pDatFile = 0;
  char szLine[MAXSIZE];
  struct residue *Seq = Sequence->residue;

  /* get sequence filename at runtime */
  if (!seq_flag)
  {
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
  }

  /* set database file (always use default helix.db) */
  strcpy(szLine,"helix.db");  //always use the default database file
  if (szLine[0] != '\0') 
    strncpy(Sequence->szDatFile, szLine, MAXSIZE);
  if (!(pDatFile = fopen(Sequence->szDatFile, "r")))
  {
    printf("\aDatabase file not found: %s\n", Sequence->szDatFile);
    return 0;
  }
  fclose(pDatFile);

  /* set pH */
  Sequence->pH = pH;

  /* set salt concentration */
  strcpy(szLine, salt);
  if ('h' == szLine[0] || 'H' == szLine[0]) 
    Sequence->salt = HIGH;
  else
    Sequence->salt = LOW;

  /* get sequence from command line or specified file */
  if (seq_flag)
  {
    if (!ReadCmdSeq(Sequence, cmdseq) || !ReadParam(Sequence))
      return 0;
  }
  else
  {
    if (!ReadSeq(Sequence) || !ReadParam(Sequence))
      return 0;
  }

  /* print out summary */
  dprint("\nThe sequence and initial parameters are:\n");
  dprint("\t  w'\t  w\t  vN\t  vC\t  N\t  C\n");
  for (nIndex = 0; nIndex < Sequence->nLength; ++nIndex)
    dprint("%2d %s\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\t%5.3f\n", 
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
