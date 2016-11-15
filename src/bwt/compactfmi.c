/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

/*
  The FM index (fmi) is used to identify the number of times the letter c occurred BEFORE
  the current location k. (This is how it is in this implementation, others may include
  current location in count).

  In this implementation, the FMI is stored only with the letter in the BWT, so at
  location k, we can access fmi ONLY for letter bwt[k]. For example:

  BWT      A    C    A    G    G    T    G    A    T    C
  fmi    107   32  108   77   78   99   79  109  100   33

  The numbers in the FMI are further coded and checkpointed.
  Additionally, we use the actual suffix array (SA)
  position of the corresponding suffix rather than the number from the
  first suffix starting with that letter (which is how the FMI is
  normally described).

  So we supply these functions:
  IndexType FMindex(FMI *fmi, uchar c, IndexType k):
           returns the FMI value at position k for letter c
  IndexType FMindexCurrent(FMI *fmi, IndexType k):
           returns the FMI value at position k for the letter in bwt[k]
           (for efficiency - could have used FMindex with the proper c)
  void FMindexAll(FMI *f, IndexType k, IndexType *fmia) {
           returns the FMI value at position k for all letters in alphabet
*/


#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "compactfmi.h"

// #define TESTING

#define ex2 8    // Exponent for checkpoints index 2 (ex2<ex1!!)

#include "fmicommon.h"


/*
  maxNcode[a] contains the max number encoded with letter a
  numbers 0..256 correponds to letter and number of letters 0..maxNcode[a]-1

  0 letter  0 number  0
  1         0         1
  2         0         2
  :         :         :
  m0-2      0         m0-2
  m0-1      0         255 (m0-1)

  m0        1         0
  m0+1      1         1
  :	    :         :
  m0+m1-1   1         255 (m1-1)

  m0+m1     2         0
  m0+m1+1   2         1
  :         :         :

*/


static uchar lcode[256];
static uchar ncode[256];

static inline uchar fmi_decode_letter(uchar code) { return lcode[code]; }
static inline uchar fmi_decode_number(uchar code) { return ncode[code]; }


static void fmi_fill_codes(int alen, int *startLcode) {
  int a, n, k;

  // for (a=0;a<alen+1;++a) fprintf(stderr,"fmi_fill_codes %d %d\n",a,startLcode[a]);

  for (a=0;a<alen;++a) {
    n=0;
    for (k=startLcode[a]; k<startLcode[a+1]-1; ++k) {
      lcode[k]=a;
      ncode[k]=n++;
    }
    lcode[k]=a;
    ncode[k]=255;
  }
}



/* Encode letter c and number n as byte code */
static inline uchar encode_letter_number(uchar c, int n, int *startLcode) {
  int max;

  max = startLcode[c+1] - startLcode[c]-1;
  if (n>max) n=max;
  return startLcode[c]+n;
}



/*
  Count the number of each letter in the BWT.
  Use the fraction to set the maxNcode.
  maxNcode has to be at least 2 for any letter, to encode 0 or more than 0.
 */
static int *find_startLcode(const int alen, const uchar *bwt, const IndexType blen) {
  IndexType i, tot=0, *count = (IndexType *)calloc(alen,sizeof(IndexType));
  int *maxNcode = (int *)calloc(alen,sizeof(int));
  int *startLcode = (int *)calloc(alen+1,sizeof(int));
  int a, sum=0, min, max=0;

  for (i=0;i<blen;++i) count[bwt[i]]+=1;

  for (a=0; a<alen;++a) tot+=count[a];
  for (a=0; a<alen;++a) {
    maxNcode[a] = 256*((double)count[a]/tot);
    if (maxNcode[a]<2) maxNcode[a] = 2;
    if (maxNcode[a]>maxNcode[max]) max = a;
    sum += maxNcode[a];
  }

  /* If numbers add up to less than 256, add to most common letter */
  if (sum<256) maxNcode[max] += 256-sum;

  /* If sum greater than 256, subtract from the least abundant letters */
  while (sum>256) {
    min=0;
    for (a=1; a<alen;++a) {
      if (maxNcode[min]<=2) min=a;
      if (maxNcode[a]>2 && maxNcode[a]<maxNcode[min]) min = a;
    }
    maxNcode[min] -= 1;
    sum -= 1;
  }

  // fprintf(stderr,"find_maxNcode ");
  // for (a=0; a<alen;++a) fprintf(stderr," %d:%d",a,maxNcode[a]);
  // fprintf(stderr,"\n");

  /* Now find start for each letter */
  startLcode[0]=0;
  startLcode[alen]=256;
  for (a=0;a<alen;++a) startLcode[a+1]=startLcode[a]+maxNcode[a];

  free(count);
  free(maxNcode);
  return startLcode;
}



FMI *alloc_FMI(uchar *bwt, IndexType bwtlen, int alen) {
  FMI *f = alloc_FMI_common(bwt, bwtlen, alen, sizeof(ushort));
  f->startLcode = find_startLcode(alen, bwt, bwtlen);
  fmi_fill_codes(alen,f->startLcode);
  return f;
}




FMI *read_fmi(FILE *fp) {
  FMI *f = read_fmi_common(sizeof(ushort),fp);
  f->startLcode = (int *)malloc((f->alen+1)*sizeof(int));
  fread(f->startLcode,sizeof(int),f->alen+1,fp);
  fmi_fill_codes(f->alen,f->startLcode);
  return f;
}




void write_fmi(const FMI *f, FILE *fp) {
  write_fmi_common(f, sizeof(ushort), fp);
  fwrite(f->startLcode,sizeof(int),f->alen+1,fp);
}




/***********************************************
 *
 * Querying FMI
 *
 *
 ***********************************************/



/* Return the number of letter c from the last checkpoint

   IT ASSUMES THAT YOU ARE AT LETTER C IN THE BWT!!!

   If the max value is encountered in the FMI value coded in BWT, this function
   will count in order to return the correct FMI value
   Note that it return the number*direction (negative for forward search)
*/
static inline int fmi_bwt2number(const uchar c, uchar *bwt, const int direction) {
  int n, k=0;

  /* Search if n==255  */
  while ( ( n = fmi_decode_number(*bwt) ) ==255 ) {
    k += 1;
    /* Find next letter equal to c */
    bwt+=direction;
    while ( fmi_decode_letter(*bwt) != c) bwt+=direction;
  }

  if (direction <0) return n+k;
  else return -n-k-1;
}




/* Get difference between checkpoints surrounding k
   For rare letters, this can tell if there is a letter of type c in the interval
   at all (not actually used at the moment)
*/
static inline int fmi_chpt_difference(FMI *f, IndexType k, uchar c, int direction) {
  IndexType chpt1, chpt1up, chpt2;
  int diff;

  chpt2 = k>>ex2;
  diff = (int)(f->index2[chpt2+1][c] - f->index2[chpt2][c]);

  /* Special case: k is just below chpt1 */
  if ( direction > 0 ) {
    chpt1 = k>>ex1;
    chpt1up = (chpt2+1)>>(ex1-ex2);
    if (chpt1 != chpt1up) diff += (int)(f->index1[c][chpt1up] - f->index1[c][chpt1]);
  }

  return diff;
}





/*
  Search for nearest letter in direction dir.
  Stop if bound is reached
  Returns NULL if letter is NOT found
*/
static inline uchar *find_closest_letter_with_bound(const uchar ct, uchar *bwt,
				       const int dir, const uchar *bound) {
  while ( ct != fmi_decode_letter(*bwt) ) {
    if (bwt == bound) { return NULL; }
    bwt += dir;
  }
  return bwt;
}



/* Return the FMI value for target letter ct at position k
   Search for closest from (and including) current position k
   and return the fmi value

   It is possible to optimize this function by using nleft as in FMindexAll

*/
IndexType FMindex(FMI *f, uchar ct, IndexType k) {
  uchar c, *bwt, *bwtstop;
  int direction;
  IndexType fmi, delta=0;

  bwt = f->bwt+k;
  if (k<f->bwtlen) c = fmi_decode_letter(*bwt);
  else c=255;
  direction=fmi_direction(k);

  fmi = fmi_chpt_value_with_dir(f, k, ct, direction);
  /* If we don't have target letter at pos k, find nearest */
  if (c != ct) {
    /* BWT position at lower checkpoint */
    bwtstop = f->bwt + (k&round2);

    /* Letter not found if we are at a check point already */
    if (bwt==bwtstop) bwt=NULL;
    else {
      /* Find the boundary for the search when direction = +1 */
      if ( direction>0 ) {
	bwtstop += size2;
	if (bwtstop>=f->bwt + f->bwtlen) {
	  bwtstop=f->bwt + f->bwtlen;
	  if (k>=f->bwtlen) bwt=bwtstop-1;
	}
	bwtstop -=1;
      }
      /* You have to add one to the fmi value if direction < 0 AND ct!=bwt[k] AND a letter is found */
      else delta=1;

      if (bwt==bwtstop) bwt=NULL;
      else bwt = find_closest_letter_with_bound(ct, bwt+direction, direction, bwtstop);
    }
  }

  /* If letter is encountered, add proper value */
  if (bwt) fmi += delta + fmi_bwt2number(ct, bwt, direction);

  return fmi;
}



/* Return the FMI value for the BWT letter at position k */
static IndexType FMindexHere(const FMI *f, uchar *bwt, const uchar c, const IndexType k) {
  int n, direction;

  // Is k above or below midpoint of index2?
  direction=fmi_direction(k);

  // Get number
  n = fmi_bwt2number(c, bwt, direction);

  return n + fmi_chpt_value_with_dir(f, k, c, direction);
}



/* Return the FMI value for the BWT letter at position k */
IndexType FMindexCurrent(FMI *f, uchar *c, IndexType k) {
  uchar *bwt;
  int n, direction;

  // Read letter
  bwt = f->bwt + k;
  *c = fmi_decode_letter(*bwt);

  return FMindexHere(f,bwt,*c,k);
}



/* Return the FMI value for all letters at position k
   Search for closest from (and including) current position k
   and return the fmi value.
   A result (fmia) array of length alen must be supplied (not checked!)
*/
void FMindexAll(FMI *f, IndexType k, IndexType *fmia) {
  uchar c, *bwt;
  int i, n, nleft, direction;
  IndexType fmi;

  bwt = f->bwt+k;
  direction=fmi_direction(k);
  for (i=0;i<f->alen;++i) fmia[i] = size2; // If letter has not been found, count = size2

  /* The total count of letters from lower checkpoint */
  nleft = k-(k&round2);

  // We do not count the current letter anyway, when dir = -1
  if (direction<0) --bwt;
  else nleft = size2 - nleft;  // boundary for the search when direction = +1

  if ( direction>0 && nleft > f->bwtlen-k) nleft = f->bwtlen-k;

  while ( nleft>0 ) {
    c = fmi_decode_letter(*bwt);
    if (c==0) break;
    // k+=direction;  // only for debug
    // DPRINT("nleft=%d k=%ld dist chkpt=%d c=%d fmi=%ld ",nleft,k,(int)(k-(k&round2)),c,fmia[c] );
    if ( fmia[c]>=size2 ) {
      n = fmi_decode_number(*bwt);
      if (n<255) {           // Letter count found
	n += fmia[c]-size2+1;
	if (direction<0) fmia[c] = n;
	else fmia[c] = -n;
	nleft -= n;
      }
      else ++fmia[c];  // Counting the number of times the letter returns 255
    }
    // DPRINT("fmi=%ld nleft=%d\n",fmia[c],nleft);
    bwt += direction;
  }

  for (i=0;i<f->alen;++i) {
    if ( fmia[i]>=size2 ) fmia[i]=0;
    fmia[i] += fmi_chpt_value_with_dir(f, k, (uchar)i, direction);
  }

}



/***********************************************
 *
 * Building FMI
 *
 *
 ***********************************************/


/*
  Assume that checkpoints are done
*/
void FMIrecode(FMI *fmi) {
  IndexType i, j, ii, R1, R2, *total;
  int a, *current, *deltaFmi;
  uchar *sbwt;

  current = (int *)calloc(fmi->alen,sizeof(int));
  deltaFmi = (int *)calloc(size2,sizeof(int));

  // Note that current char is NOT counted
  i=0;
  R2=0;
  sbwt=fmi->bwt;
  for (ii=0; ii<fmi->bwtlen; ++ii) {
    /* Check if we are at a checkpoint 2 */
    if ( ii>0 && !(ii&check2) ) {
      R2 = ii>>ex2;
      /* Insert values in checkpoint 2 */
      for (j=0; j<size2>>1;++j) sbwt[j] = encode_letter_number(sbwt[j],deltaFmi[j],fmi->startLcode);
      for (   ; j<size2;   ++j)	sbwt[j] = encode_letter_number(sbwt[j],(current[sbwt[j]] - deltaFmi[j])-1,fmi->startLcode);
      // for (   ; j<size2;   ++j)	sbwt[j] = encode_letter_number(sbwt[j],(current[sbwt[j]] - deltaFmi[j]),fmi->startLcode);
      for (a=0;a<fmi->alen;++a) current[a]=0;
      i=0;
      sbwt+=size2;
    }
    a = sbwt[i];
    deltaFmi[i] = current[a];
    current[a] += 1;
    ++i;
  }

  /* Code differences in BWT */
  for (j=0; j<(size2>>1) && j<i; ++j) sbwt[j] = encode_letter_number(sbwt[j],(uchar)deltaFmi[j],fmi->startLcode);
  for (   ; j<size2      && j<i; ++j) sbwt[j] = encode_letter_number(sbwt[j],(uchar)(current[sbwt[j]] - deltaFmi[j])-1,fmi->startLcode);


  free(current);
  free(deltaFmi);
}




/* 
*/
FMI *makeIndex(uchar *bwt, long bwtlen, int alen) {
  FMI *fmi;

  fmi = makeIndex_common(bwt, bwtlen, alen);
  FMIrecode(fmi);
  return fmi;
}






#ifdef TESTING
#include "testing.c"
#endif
