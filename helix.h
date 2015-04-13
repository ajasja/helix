#ifndef HELIX_H
#define HELIX_H

#include <stdio.h>
#include <stdbool.h>

/* maximum string size */
#define MAXSIZE 200
/* maximum number of residues */
#define MAXLEN 60
/* salt concentraion */
#define LOW 1
/* salt concentraion */
#define HIGH 2
/* parameter changed flag */
#define ON 3
/* parameter unchanged flag */
#define OFF 4
#define YES ('y' == szLine[0] || 'Y' == szLine[0] || '\0' == szLine[0])

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) > (y) ? (y) : (x))

#define PRINT_HEADER(szComment,nFirst,nLast) printf(szComment);\
printf(" %s(%d) and %s(%d)\n",Seq[nFirst].szName,(nFirst)+1,\
Seq[nLast].szName,(nLast)+1)

#ifdef DEBUG
    #define dprint printf
#else
    #define dprint
#endif

struct residue {
  char  szName[MAXSIZE];
  float w1; 
  float w; 
  float vN; 
  float vC; 
  float N; 
  float C;
  char  cCharge;
  float Pw1;
  float Pw;
  float PvN;
  float PvC;
  float fH;
};

struct sequence {
  struct residue residue[MAXLEN];
  int            nLength;
  float          pH;
  int            salt;
  char           szDatFile[MAXSIZE];
  char           szSeqFile[MAXSIZE];
  char           szComment[MAXSIZE];
};

int ReadData  (struct sequence *Sequence, char* cmdseq, bool seq_flag, char* salt, float pH);
int ReadParam (struct sequence *Sequence);
int ReadSeq   (struct sequence *Sequence);
int ReadCmdSeq(struct sequence *Sequence, char* cmdseq);
int SetCharge (struct sequence *Sequence);

int FindYFW   (struct sequence *Sequence, bool YFW_flag);
int FindRH    (struct sequence *Sequence, bool RH_flag);
int FindNcap  (struct sequence *Sequence, bool cap_flag);
int FindCcap  (struct sequence *Sequence, bool cap_flag);
int FindNend  (struct sequence *Sequence, bool end_flag);
int FindCend  (struct sequence *Sequence, bool end_flag);
int FindLipo  (struct sequence *Sequence, bool lipo_flag);
int FindHbond (struct sequence *Sequence, bool Hbond_flag);
int FindCoul  (struct sequence *Sequence, bool Coul_flag);

int CalcProb (struct sequence *Sequence);

void PrintParam (struct residue *Seq, int nFirst, int nLast, int *aFlag);
void PrintProb (struct sequence *Sequence, char* out);
void PrintProbMin (struct sequence *Sequence, char* out);

#endif
