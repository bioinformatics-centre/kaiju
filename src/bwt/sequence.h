/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
#ifndef SEQUENCE_h
#define SEQUENCE_h

#include "common.h"

typedef struct __SEQstruct__ {
  char *id;
  char *descr;     // Description (stuff following id). Allocated with id.
  char rc;         // Reverse complement. 0=forward 1=complement
  long len;
  long pos;        // Index of sequence in long allocation
  char *start;
  long id_filepos;
  long seq_filepos;
  int sort_order;
  struct __SEQstruct__ *next;
} SEQstruct;

typedef struct {
  int len;        // Alphabet length
  int caseSens;
  char *a;        // Alphabet sequence (0 terminated, first char is terminator, last may be wildcard)
  char *trans;    // Translate char c to int i: trans[c]=i
  char *comp;     // DNA complement comp[a]=t, etc.
} AlphabetStruct;



static inline int letter2number(char c, AlphabetStruct *a) { return a->trans[(int)c];}
static inline char number2letter(int i, AlphabetStruct *a) { return a->a[i];}


/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
SEQstruct *alloc_SEQstruct();
void free_SEQstruct(SEQstruct *ss);
AlphabetStruct *alloc_AlphabetStruct(char *a, int caseSens, int revcomp);
void free_AlphabetStruct(AlphabetStruct *astruct);
void translate2numbers(uchar *s, const IndexType slen, AlphabetStruct *astruct);
/* FUNCTION PROTOTYPES END */


#endif
