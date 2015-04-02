#include "helix.h"
#include <string.h>

int ReadSeq (struct sequence *Sequence)
{
/*********************************************************
 * Description: this subroutine reads sequence (.seq) file 
 * it returns 1 for success, 0 for failure
 *********************************************************/
  int nLength = 0;
  char szLine[MAXSIZE], szFirstLine[MAXSIZE];
  char *pszName = 0;
  FILE *pSeqFile = 0;
  struct residue *Seq = Sequence->residue;
  int pH = Sequence->pH;

  pSeqFile = fopen(Sequence->szSeqFile, "r");
  fgets(szFirstLine, MAXSIZE, pSeqFile);

  /* check to see if there are any comments on top */
  while (!strchr(szFirstLine, '{'))
    if (!feof(pSeqFile))
    {
      strncpy(Sequence->szComment, szFirstLine, MAXSIZE);
      fgets(szFirstLine, MAXSIZE, pSeqFile);
    }
    else
    {
      printf("\aFile input error: %s\n", Sequence->szSeqFile);
      return 0;
    }

  /* go through the sequence to get individual residue name */
  strcpy(szLine, " ");
  strcat(szLine, szFirstLine);
  while (!feof(pSeqFile))
  {
    if (!nLength) 
      strtok(szLine,"{");
    else
    {
      pszName = strtok(szLine, ",}\n");
      strncpy(Seq[nLength++].szName, pszName, MAXSIZE);
    }
    while (pszName = strtok(0, ",}\n"))
      strncpy(Seq[nLength++].szName, pszName, MAXSIZE);
    fgets(szLine, MAXSIZE, pSeqFile);

    if (nLength >= MAXLEN)
    {
      printf("File input error: exceeding maximum length!\n");
      return 0;
    }
  }

  Sequence->nLength = nLength;
  SetCharge(Sequence);

  fclose(pSeqFile);
  return 1;
}
