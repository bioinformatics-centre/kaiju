/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef SUFFIXARRAY_h
#define SUFFIXARRAY_h

#include "common.h"
#include "fmi.h"
#include <pthread.h>

typedef struct {
  IndexType len;      // Length of (actual) SA ( = bwtlen)

  // Parallel processing
  FMI *f;             // Pointer to suffix array
  int interval;       // Interval used for parallel processing of sequences when building suffix array
  int next_start;
  pthread_mutex_t lock;

  // Suffix array checkpoints
  IndexType ncheck;   // Number of checkpoints
  uchar *sa;          // Actual array holding SA checkpoints
  int chpt_exp;       // Exponent of checkpoint distance
  int nbytes;         // Number of bytes used per entry
  int sbits;          // Number of bits used for encoding sequence number
  int pbits;          // Number of bits used to encode position
  long mask;          // Mask for lowest pbits bits
  long check;         // Used to check if we are at a checkpoint

  // Alphabet
  int alen;
  char *alphabet;

  // Sequence information
  int nseq;           // Number of sequences
  char **ids;         // IDs
  int *sortorder;     // Sort order for sequences
  IndexType *seqlengths;  // lengths of sequences
  IndexType maxlength;    // Maximum length of sequences
  int *seq;           // Sequence number corresponding to termination symbol number i

} BWT;



typedef struct _SI_ {
  IndexType start;  // Start of suffix interval
  int len;          // Interval length
  int qi;           // Position in query
  int ql;           // Length in query (if relevant)
  int count;        // Used to count matches below current
  struct _SI_ *next;
  struct _SI_ *samelen;
} SI;



/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
void translate2numbers(uchar *s, const IndexType slen, char *alphabet, int alen);
void translate2numbers_ptr(uchar *s, const IndexType slen, uchar * trans); 
BWT *read_BWT_header(FILE *bwtfile);
uchar *read_BWT(IndexType len, int alen, char *alphabet, FILE *bwtfile);
void write_BWT_bin(const BWT *s, FILE *fp);
BWT *read_BWT_bin(FILE *fp);
void init_suffixArray(BWT *s, FMI *f, int chpt_exp);
void fill_suffixArray(int from, int to, BWT *s);
void *fill_suffixArray_parallel(BWT *s);
void get_suffix(FMI *fmi, BWT *s, IndexType i, int *iseq, IndexType *pos);
uchar *retrieve_seq(int snum, BWT *s, FMI *f);
IndexType UpdateSI(FMI *f, uchar ct, IndexType *si, IndexType *newsi);
void recursive_free_SI(SI *si);
SI *maxMatches(FMI *f, char *str, int len, int L, int max_matches);
SI *maxMatches_fromlastposonly(FMI *f, char *str, int len, int L, int max_matches);
SI *maxMatches_withStart(FMI *f, char *str, int len, int L, int max_matches, IndexType si0, IndexType si1, int offset);
SI *allMatches(FMI *f, char *str, int len, int L);
void SIprintQmatch(SI *si, char *query, char *alphabet, FILE *fp);
void printMatches(FMI *f, char *str, int len, int L, char *alphabet);
/* FUNCTION PROTOTYPES END */

#endif
