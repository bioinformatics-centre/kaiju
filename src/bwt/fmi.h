/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
#ifndef COMPACTFMI_h
#define COMPACTFMI_h

#include "common.h"

/* Simple FM index with a checkpoint for every 256 letters */
typedef struct {
  int alen;           // Length of alphabet
  IndexType bwtlen;   // Total length of BWT
  uchar *bwt;         // BWT string
  int N1;             // Total number of entries in index 1 (bwtlen>>ex1 +1);
  int N2;             // Total number of entries in index 2 (bwtlen>>ex2 +1);
  IndexType **index1; // FM index1 (one array per letter)
  ushort **index2;    // Counts relative to index1 checkpoints (assuming 16 bit int)
  int *startLcode;    // start numbers for byte encoding of letter and number
} FMI;




/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
FMI *alloc_FMI(uchar *bwt, IndexType bwtlen, int alen);
FMI *read_fmi(FILE *fp);
void write_fmi(const FMI *f, FILE *fp);
IndexType FMindex(FMI *f, uchar ct, IndexType k);
IndexType FMindexCurrent(FMI *f, uchar *c, IndexType k);
void FMindexAll(FMI *f, IndexType k, IndexType *fmia);
void FMIrecode(FMI *fmi);
FMI *makeIndex(uchar *bwt, long bwtlen, int alen);
FMI *makeIndex_OLD(uchar *bwt, long bwtlen, int alen);
/* FUNCTION PROTOTYPES END */


#endif
