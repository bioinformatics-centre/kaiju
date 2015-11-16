/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */


/*
  The FM index (fmi) is used to identify the number of times the letter c occurred BEFORE
  the current location k. (This is how it is in this implementation, others may include
  current location in count).

  We supply these functions:
  IndexType FMindex(FMI *fmi, uchar c, IndexType k):
           returns the FMI value at position k for letter c
  IndexType FMindexCurrent(FMI *fmi, IndexType k):
           returns the FMI value at position k for the letter in bwt[k]
           (for efficiency - could have used FMindex with the proper c)
  void FMindexAll(FMI *f, IndexType k, IndexType *fmia) {
           returns the FMI value at position k for all letters in alphabet
  int UpdateSI(FMI *fmi, uchar c, IndexType *si):
           returns next suffix interval for char c and interval si.
	   OVER-WRITES si!!
           Returns length of new interval
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "simplefmi.h"

//#define TESTING


#define ex2 8    // Exponent for checkpoints index 2 (ex2<ex1!!)

#include "fmicommon.h"


FMI *alloc_FMI(uchar *bwt, IndexType bwtlen, int alen) {
  FMI *f = alloc_FMI_common(bwt, bwtlen, alen, sizeof(ushort));
  return f;
}



FMI *read_fmi(FILE *fp) {
  FMI *f = read_fmi_common(sizeof(ushort),fp);
  return f;
}



void write_fmi(const FMI *f, FILE *fp) {
  write_fmi_common(f, sizeof(ushort), fp);
  // fwrite(f->startLcode,sizeof(int),f->alen+1,fp);
}



/***********************************************
 *
 * Querying FMI
 *
 *
 ***********************************************/




static int lettercount(const uchar letter, const uchar *str, const int len) {
  int i, c=0;
  for (i=0; i<len; ++i) if (str[i]==letter) ++c;
  return c;
}



/* Get the checkpointed FMI value for k  */
static inline IndexType fmi_chpt_value(const FMI *f, const IndexType k, const uchar c) {
  return f->index1[k>>ex1][c] + f->index2[k>>ex2][c];
}





/* Return the FMI value for target letter ct at position k */
static inline IndexType FMindex_undirectional(FMI *f, uchar ct, IndexType k) {
  IndexType k2 = k&round2;

  return fmi_chpt_value(f,k,ct)+lettercount(ct,f->bwt+k2,k-k2);
}



/* Return the FMI value for target letter ct at position k */
static inline IndexType FMindex_directional(FMI *f, uchar ct, IndexType k) {
  int dir = fmi_direction(k);
  IndexType k2 = k&round2;

  if (dir<0) return fmi_chpt_value_with_dir(f,k,ct,dir)+lettercount(ct,f->bwt+k2,k-k2);
  else return fmi_chpt_value_with_dir(f,k,ct,dir)-lettercount(ct,f->bwt+k,fmi_end_length(k, f->bwtlen));
}



IndexType FMindex(FMI *f, uchar ct, IndexType k) {
  return FMindex_directional(f, ct, k);
}


/* Return the letter (in *c) and the FMI value for the BWT letter at
   position k
*/
IndexType FMindexCurrent(FMI *f, uchar *c, IndexType k) {
  *c=f->bwt[k];
  return FMindex(f, *c, k);
}


/*
  Return the FMI value for all letters at position k
  A result (fmia) array of length alen must be supplied (not checked!)
  Not much optimization can be done, so it is just using the FMindex function.
*/
void FMindexAll(FMI *f, IndexType k, IndexType *fmia) {
  uchar c;
  for (c=0; c<f->alen; ++c) fmia[c] = FMindex(f, c, k);
}





/***********************************************
 *
 * Building FMI
 *
 *
 ***********************************************/


FMI *makeIndex(uchar *bwt, long bwtlen, int alen) {
  return makeIndex_common(bwt, bwtlen, alen);
}


#ifdef TESTING
#include "testing.c"
#endif
