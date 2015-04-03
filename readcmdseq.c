#include "helix.h"
#include <string.h>

int ReadCmdSeq (struct sequence *Sequence, char* cmdseq)
{
/****************************************************************
 * Description: this subroutine reads sequence from string cmdseq 
 * it returns 1 for success, 0 for failure
 ***************************************************************/
  int nLength = 0;
  char szLine[MAXSIZE], szFirstLine[MAXSIZE];
  char *pszName = 0;
  struct residue *Seq = Sequence->residue;
  int pH = Sequence->pH;

  strcpy(szFirstLine, cmdseq);

  strcpy(szLine, " ");
  strcat(szLine, szFirstLine);
  if (!nLength) 
    strtok(szLine,"{");
  else
  {
    pszName = strtok(szLine, ",}\n");
    strncpy(Seq[nLength++].szName, pszName, MAXSIZE);
  }
  while (pszName = strtok(0, ",}\n"))
    strncpy(Seq[nLength++].szName, pszName, MAXSIZE);

  if (nLength >= MAXLEN)
  {
    printf("Sequence input error: exceeding maximum length!\n");
    return 0;
  }

  Sequence->nLength = nLength;
  SetCharge(Sequence);

  return 1;
}
