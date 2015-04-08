/**************************************************
 * Project: HELIX
 * Version: 1.6
 * File: helix.c,
 * Author: Hui Tong
 * Modified by Igor Drobnak
 * (C) 1995-1997 University of Washington
 * (C) 2015 National institute of Chemistry, Slovenia
 * Last Updated: 3/4/2015
 * Description: this program simulates 2-state
 * helix-coil transitions in peptides and proteins
 **************************************************/

#include "helix.h"
#include <string.h>
#include <stdio.h>

void main(int argc,char *argv[]) 
{
  int fNcap = 0;
  struct sequence Sequence;
  int i=0;  //iterator
  
  /* parameters, read in from the command line */
  char cmdseq[MAXSIZE];  //sequence to be read from the command line
  bool seq_flag=false;
  bool YFW_flag, RH_flag, cap_flag, end_flag, lipo_flag, Hbond_flag, Coul_flag, all_flag=false;
  bool min_flag=false;
  char salt[MAXSIZE];  //string indicating salt concentration (H=hi / L=low)
  strcpy(salt, "L");  //default salt concentration is low
  float pH=7.0;  //pH at which to calculate charges
  char out[MAXSIZE];  //name of output file
  strcpy(out, "helix.out");

  /* read command-line arguments */
  printf("\n****************** Helix v1.6 *******************\n");
  dprint("No. of command line parameters: %d\n", argc-1);
  for (i=1; i<argc; i++)
  {
    dprint("argument %d: %s\n", i, argv[i]);
    
    // read sequence
    if (strstr(argv[i], "--seq="))
    {
      strcpy(cmdseq, strpbrk(argv[i], "=")+1);
      printf("Input sequence: %s\n", cmdseq);
      seq_flag=true;
    }
    
    // read pH
    else if (strstr(argv[i], "--pH="))
    {
      char* temp=strpbrk(argv[i], "=")+1;
      if(sscanf(temp, "%f", &pH)!=1 || pH<0.0 || pH>13.0)
      {
        pH=7.0;
        printf("Invalid pH given: %s. Falling back to default (7.0).\n",temp);
      }
      printf("pH: %f\n", pH);
    }
    
    // read salt concentration
    else if (strstr(argv[i], "--salt="))
    {
      strcpy(salt,strpbrk(argv[i], "=")+1);
      printf("Salt concentration: %s\n", salt);
    }
    
    // silently include all features in the model
    else if (strcmp(argv[i], "--findAll")==0)
    {
      printf("Looking up all features\n");
      all_flag=true;
    }
    
    // print only minimal output
    else if (strcmp(argv[i], "--min")==0)
    {
      printf("Writing minimal output only\n");
      min_flag=true;
    }
    
    // read name of output file
    else if (strstr(argv[i], "--out="))
    {
      strcpy(out,strpbrk(argv[i], "=")+1);
      printf("Output to file: %s\n", out);
    }
    
    // print out help information
    else if (strcmp(argv[i], "-h") || strcmp(argv[i], "--help"))
    {
      puts("\n**********************");
      puts("Helix program, v 1.6");
      puts("Original code by Hui Tong, modified by Igor Drobnak");
      puts("**********************");
      puts("Usage: helix [options]\n");
      puts("--seq=\"{S,E,Q}\"   performs calculation on the specified sequence (default asks for filename at runtime)");
      puts("                  Note: capital letters, commas, quotation marks and braces \"{}\" are all obligatory\n");
      puts("--pH=<x>          sets pH value to <x> (default = 7.0)");
      puts("                  Note: <x> must be a number between 0.0 and 13.0\n");
      puts("--salt=[H|h]      uses high salt concentration in calculations (default = low salt)\n");
      puts("--findAll         use all available features in the model calculations (default asks about each feature at runtime)\n");
      puts("--min             produces a minimal output file only (single column containing % helicity)\n");
      puts("--out=<filename>  writes output to file <filename> (default = helix.out)\n");
      puts("--help, -h        print this help information\n");
      return;
    }
  }
  
  if (all_flag)  //set all individual flags to true
    YFW_flag=RH_flag=cap_flag=end_flag=lipo_flag=Hbond_flag=Coul_flag= true;
  
  printf("----------------\n");
  
/* input data from sequence and parameter files */
  if (!ReadData(&Sequence, cmdseq, seq_flag, salt, pH))
    return;
  
/* search for Tyr interactions */
  FindYFW(&Sequence, YFW_flag);

/* search for Arg/His at C-terminus */
  FindRH(&Sequence, RH_flag);

/* searching for N-capping boxes */
  fNcap = FindNcap(&Sequence, cap_flag);

/* searching for C-capping boxes */
  if (fNcap)
    FindCcap(&Sequence, cap_flag);

/* searching for end effects at N terminus */
  if (0 == fNcap || 1 == fNcap)
    FindNend(&Sequence, end_flag);

/* searching for end effects at C terminus */
  FindCend(&Sequence, end_flag);

/* searching for lipophilic (hydrophobic) interactions */
  FindLipo(&Sequence, lipo_flag);

/* searching for hydrogen-bonding effect*/
  FindHbond(&Sequence, Hbond_flag);

/* searching for Coulombic interactions */
  FindCoul(&Sequence, Coul_flag);

/* calculate conformation probabilities at each residue */
  CalcProb(&Sequence);

/* print them out */
  if(min_flag)
    PrintProbMin(&Sequence, out);
  else
    PrintProb(&Sequence, out);
  
  dprint("Done!\n");
}
