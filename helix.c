/**************************************************
 * Project: HELIX
 * Version: 1.5
 * File: helix.c,
 * Author: Hui Tong
 * (C) 1995-1997 University of Washington
 * Last Updated: 1/8/97
 * Description: this program simluates 2-state
 * helix-coil transitions in peptides and proteins
 **************************************************/

#include "helix.h"

void main(int argc,char *argv[]) 
{
  int fNcap = 0;
  struct sequence Sequence;

/* input data from sequence and parameter files */
  if (!ReadData(&Sequence))
    return;
/* search for Tyr interactions */
  FindYFW(&Sequence);

/* search for Arg/His at C-terminus */
  FindRH(&Sequence);

/* searching for N-capping boxes */
  fNcap = FindNcap(&Sequence);

/* searching for C-capping boxes */
  if (fNcap)
    FindCcap(&Sequence);

/* searching for end effects at N terminus */
  if (0 == fNcap || 1 == fNcap)
    FindNend(&Sequence);

/* searching for end effects at C terminus */
  FindCend(&Sequence);

/* searching for lipophilic (hydrophobic) interacions */
  FindLipo(&Sequence);

/* searching for hydrogen-bonding effect*/
  FindHbond(&Sequence);

/* searching for Coulombic interactions */
  FindCoul(&Sequence);

/* calculate conformation probabilities at each residue */
  CalcProb(&Sequence);

/* print them out */
  PrintProb(&Sequence);
}
