/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
#ifndef SUFFIXARRAY_h
#define SUFFIXARRAY_h

#include "common.h"
#include "fmi.h"
#include "sequence.h"

typedef struct {
  IndexType len;      // Length of (actual) SA ( = bwtlen)

  // Suffix array checkpoints
  IndexType ncheck;   // Number of checkpoints
  uchar *sa;          // Actual array holding SA checkpoints
  int chpt_exp;       // Exponent of checkpoint distance
  int nbytes;         // Number of bytes used per entry
  int sbits;          // Number of bits used for encoding sequence number
  int pbits;          // Number of bits used to encode position
  long mask;          // Mask for lowest pbits bits
  long check;         // Used to check if we are at a checkpoint

  // Sequence information
  int nseq;              // Number of sequences
  char **ids;            // IDs (in order of forward sorted seqs)
  int *seqTermOrder;     // Order of sequence termination
  IndexType *seqlengths; // lengths of sequences
  IndexType maxlength;   // Maximum length of sequences
  int hash_step;         // distance in hash table
  SEQstruct **hash;
  char *seqstart;   // Start of sequence

} suffixArray;


/* Decode long from n bytes */
static inline long uchar2long(uchar *c, int n) {
  long val=*c++;
  while ( --n >0 ) val = (val<<8) + *c++;
  return val;
}

/*
  For SA entry k, return seq no. (in *nseq) and position within (*pos)
  Entry consists of nbytes bytes starting at position sa+k*nbytes.
*/
static inline void suffixArray_decode_number(int *nseq, long *pos, long k, suffixArray *s) {
  long val = uchar2long( (s->sa + k * s->nbytes), s->nbytes);
  *nseq = (int)(val>>s->pbits);
  *pos = val & s->mask;
}





/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
void suffixArray_make_hash(SEQstruct *base, suffixArray *s, int Hstep);
suffixArray *init_suffixArray(SEQstruct *ss, int chpt_exp);
void write_suffixArray_checkpoints(char **sa, IndexType start, IndexType length,
				   suffixArray *s, FILE *sa_file);
void write_suffixArray_header(suffixArray *s, FILE *fp);
suffixArray *read_suffixArray_header(FILE *fp);
void read_suffixArray_body(suffixArray *s, FILE *fp);
void write_suffixArray(suffixArray *s, FILE *fp);
/* FUNCTION PROTOTYPES END */

#endif
