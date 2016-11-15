/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
#ifndef BWT_h
#define BWT_h

#include "common.h"
#include "fmi.h"
#include "suffixArray.h"

typedef struct {
  IndexType len;      // Length of bwt (not counting initial zeros)
  int nseq;
  uchar *bwt;

  // Alphabet
  int alen;
  char *alphabet;

  FMI *f;
  suffixArray *s;

} BWT;


typedef struct _SI_ {
  IndexType start;  // Start of suffix interval
  int len;          // Interval length
  int qi;           // Position in query
  int ql;           // Length in query (if relevant)
  int count;        // Used to count matches below current
  int score;
  struct _SI_ *next;
  struct _SI_ *samelen;
} SI;



/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
void write_BWT_header(BWT *b, FILE *bwtfile);
BWT *read_BWT(FILE *bwtfile);
BWT *readIndexes(FILE *fp);
void get_suffix(FMI *fmi, suffixArray *s, IndexType i, int *iseq, IndexType *pos);
uchar *retrieve_seq(int snum, BWT *b);
IndexType InitialSI(FMI *f, uchar ct, IndexType *si);
IndexType UpdateSI(FMI *f, uchar ct, IndexType *si, IndexType *newsi);
void recursive_free_SI(SI *si);
SI *maxMatches(FMI *f, char *str, int len, int L, int max_matches);
SI *maxMatches_withStart(FMI *f, char *str, int len, int L, int max_matches, IndexType si0, IndexType si1, int offset);
SI *greedyExact(FMI *f, char *str, int len, int L, int jump);
/* FUNCTION PROTOTYPES END */

#endif
