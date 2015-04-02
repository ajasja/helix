#include "helix.h"
#include <string.h>

int SetCharge(struct sequence *Sequence)
/********************************************************************* 
 * Description: this subroutine sets the charge of a residue according 
 * to the pH
 * it returns 1 for success, 0 for failure
 *********************************************************************/
{
  int nIndx = 0;
  struct residue *Seq = Sequence->residue;
  int pH = Sequence->pH;

  for (nIndx = 0; nIndx < Sequence->nLength; ++nIndx)
  {
    Seq[nIndx].cCharge = 0; /* set the default charge to neutral */
    switch(Seq[nIndx].szName[0])
    {
      case 'D': if (pH >= 3.8)  
                  Seq[nIndx].cCharge = '-';
                if ('\0' == Seq[nIndx].szName[1])
                  strcat(Seq[nIndx].szName, (pH >= 3.8 ? "-" : "0"));
                break;
      case 'E': if (pH >= 4.6)  
                  Seq[nIndx].cCharge = '-';
                if ('\0' == Seq[nIndx].szName[1])
                  strcat(Seq[nIndx].szName, (pH >= 4.6 ? "-" : "0"));
                break;
      case 'H': if (pH <= 6.8)   
                  Seq[nIndx].cCharge='+';
                if ('\0' == Seq[nIndx].szName[1])
                  strcat(Seq[nIndx].szName, (pH <= 6.8 ? "+" : "0"));
                break;
      case 'K': if (pH < 11.5)
                  Seq[nIndx].cCharge = '+';
                if ('\0' == Seq[nIndx].szName[1])
                  strcat(Seq[nIndx].szName, (pH < 11.5 ? "+" : "0"));
                break;
      case 'R': Seq[nIndx].cCharge = '+';
                if ('\0' == Seq[nIndx].szName[1])
                  strcat(Seq[nIndx].szName, "+");
                break;
      default:  if (!strcmp(Seq[nIndx].szName, "Suc"))
                {
                  if (pH >= 4.6)  
                    Seq[nIndx].cCharge = '-';
                  if ('\0' == Seq[nIndx].szName[3])
                    strcat(Seq[nIndx].szName, (pH >= 4.6 ? "-" : "0"));
                }
    }
  }

  return 1;
}
