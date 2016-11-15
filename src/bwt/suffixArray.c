/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
/**********************************************************************

   Functions for check-pointed suffix array


 **********************************************************************/



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <time.h>

#include "sequence.h"
// #include "bwt.h"
#include "suffixArray.h"



/***********************************************
 *
 * Bit magic
 *
 ***********************************************/



static inline void suffixArray_set_masks(suffixArray *s) {
  s->mask = (1<<(s->pbits))-1;
  s->check = (1<<(s->chpt_exp))-1;
}


/* Encode long in n bytes (assumes that number can fit in n bytes) */
static inline void long2uchar(long k, uchar *c, int n) {
  while ( n-- >0 ) { c[n] = k; k = k>>8; }
}



/*
  Encode seq no. (nseq) and position (pos) in nbytes bytes in *code
*/
static inline void suffixArray_encode_number(int nseq, long pos, uchar *code, suffixArray *s) {
  long val = nseq;
  val = (val << s->pbits)  + pos;
  long2uchar(val, code, s->nbytes);
}



/* Returns the bits needed to encode a number (without using log2) */
static int bitsNeeded(long k) { int i=0; while (k>>i) ++i; return i; }




/***********************************************
 *
 * Sequence look-up
 *
 * Used only during construction of SA
 *
 ***********************************************/



/*
  To get from a position in the concatenated sequence to an actual sequence
  a hashing scheme is introduced. For every Hstep there is a
  pointer to the SEQstruct.
*/
void suffixArray_make_hash(SEQstruct *base, suffixArray *s, int Hstep) {
  SEQstruct *ss, **spt;
  long i, limit;

  s->hash_step = Hstep;
  // The number of entries in hash table
  i = 2+base->len/s->hash_step;
  // alloc table
  spt = s->hash = (SEQstruct **)malloc(i*sizeof(SEQstruct *));
  s->seqstart = base->start;

  ss=base;
  i=0;
  // Walk through all sequences
  while (ss->next) {
    ss = ss->next;
    limit = (long)(ss->start - s->seqstart) + ss->len;
    while (i*s->hash_step <= limit) spt[i++] = ss;
  }
}



/*
  Get pointer to sequence for a given position in concatenated string.
*/
static SEQstruct *hash_lookupSeq(char *suffix, suffixArray *s) {
  int i;
  //long pos = (long)(suffix - s->seqstart);
  SEQstruct *ss;

  //i = pos/s->hash_step;
  i = (long)(suffix - s->seqstart)/s->hash_step;
  ss = s->hash[i];
  // while ( pos > ss->pos + ss->len ) ss=ss->next;
  while ( suffix > ss->start + ss->len ) ss=ss->next;
  return ss;
}





/***********************************************
 *
 * Build SA, output SA
 *
 ***********************************************/






/* Allocate and initialize suffix array structure to hold SA checkpoint values
   Does NOT allocate space for array (->sa, ->ids, etc)
*/
suffixArray *init_suffixArray(SEQstruct *ss, int chpt_exp) {
  IndexType bwtlen, maxseqlen;
  int nseq;
  SEQstruct *cur;
  suffixArray *s = (suffixArray *)malloc(sizeof(suffixArray));

  /* Calculate the real length of the BWT and find max sequence length */
  cur = ss->next;
  bwtlen=nseq=0;
  maxseqlen=0;
  while (cur) {
    bwtlen += cur->len+1;
    if (cur->len>maxseqlen) maxseqlen=cur->len;
    ++nseq;
    cur = cur->next;
  }

  // fprintf(stderr,"bwtlen %ld nseq %d\n",bwtlen,nseq);

  s->len = bwtlen;
  s->nseq = nseq;
  s->chpt_exp = chpt_exp;
  s->maxlength = maxseqlen;
  s->hash_step = 0;

  s->ncheck = (s->len>>chpt_exp) - (s->nseq>>chpt_exp);
  s->sbits = bitsNeeded(s->nseq);
  s->pbits = bitsNeeded(s->maxlength);
  s->nbytes = (7+s->sbits+s->pbits)/8;
  suffixArray_set_masks(s);

  // s->sa = (uchar*)malloc(s->ncheck*s->nbytes*sizeof(uchar));
  // s->seq = (int*)malloc(s->nseq*sizeof(int));

  // Make seq hash. For step length use (arbitrarily) half the av. seq length
  suffixArray_make_hash(ss, s, (int)(0.5*bwtlen/nseq));

  s->sa=NULL;
  s->seqTermOrder=NULL;
  s->seqlengths=NULL;

  return s;
}




static void writeseq(char *s, int l, char *term, FILE *fp) {
  int i, n=0;
  char *a = "*ACDEFGHIKLMNPQRSTVWYX";
  for (i=0; i<l; ++i) {
    if (n || s[i]==0) { n=1; fprintf(fp,"%c",a[0]); }
    else fprintf(fp,"%c",a[s[i]]);
  }
  if (term) fprintf(fp,"%s",term);
}


/* Go through a suffix array and look up SA checkpoints, and write in files
 */
void write_suffixArray_checkpoints(char **sa, IndexType start, IndexType length,
				   suffixArray *s, FILE *sa_file) {
  IndexType i, k;
  uchar code[32];
  SEQstruct *seq;

  //  fprintf(stderr,"chkpt start=%ld end=%ld\n",start,start+length);


  /* NB: The suffixes being sorted do not include sequence ends (0), so
     the actual start is start+nseq
  */
  k = start + s->nseq;
  // k = start;

  for (i=0; i<length; ++i, ++k) {
    if ( !(k&s->check) ) {
      // Use position in long concatenated sequence
      seq=hash_lookupSeq(sa[i], s);
      suffixArray_encode_number(seq->sort_order,(long)(sa[i]-seq->start), code, s);
      fwrite(code,1,s->nbytes,sa_file);
      --(s->ncheck);

      /*
      if ( seq->start > sa[i] || seq->start+seq->len<sa[i]) {
	fprintf(stderr,"ERROR: k=%ld pos=%ld start=%ld len=%ld iseq=%d id=%s\n",
		k,(long)(sa[i]-seq->start),seq->pos,seq->len,seq->sort_order,seq->id);
      }
      */
    }
  }
}




/* Assumes that length of string is <256 - otherwise truncate */
static inline void writeString(char *str, FILE *fp) {
  uchar ul;
  int l;
  l = strlen(str);
  if (l>255) l=255;
  ul = l;
  fwrite(&ul,sizeof(uchar),1,fp);
  fwrite(str,sizeof(char),l,fp);
}



static inline char *readString(FILE *fp) {
  uchar ul;
  int l;
  char *str;
  fread(&ul,sizeof(uchar),1,fp);
  l = ul;
  str = malloc((l+1)*sizeof(char));
  fread(str,sizeof(char),l,fp);
  str[l]='\0';
  return str;
}





/* Write SA header */
void write_suffixArray_header(suffixArray *s, FILE *fp) {
  int i;
  fwrite(&(s->len),sizeof(IndexType),1,fp);
  fwrite(&(s->ncheck),sizeof(IndexType),1,fp);
  fwrite(&(s->chpt_exp),sizeof(int),1,fp);
  fwrite(&(s->nbytes),sizeof(int),1,fp);
  fwrite(&(s->sbits),sizeof(int),1,fp);
  fwrite(&(s->pbits),sizeof(int),1,fp);
  fwrite(&(s->mask),sizeof(long),1,fp);
  fwrite(&(s->check),sizeof(long),1,fp);

  fwrite(&(s->nseq),sizeof(int),1,fp);
  for (i=0; i<s->nseq;++i) writeString(s->ids[i], fp);
  fwrite(s->seqTermOrder,sizeof(int),s->nseq,fp);  
  fwrite(s->seqlengths,sizeof(IndexType),s->nseq,fp);
  // fwrite(&(s->maxlength),sizeof(IndexType),1,fp);
}



/* Read SA  */
suffixArray *read_suffixArray_header(FILE *fp) {
  int i;
  suffixArray *s = (suffixArray*)malloc(sizeof(suffixArray));

  fread(&(s->len),sizeof(IndexType),1,fp);
  fread(&(s->ncheck),sizeof(IndexType),1,fp);
  fread(&(s->chpt_exp),sizeof(int),1,fp);
  fread(&(s->nbytes),sizeof(int),1,fp);
  fread(&(s->sbits),sizeof(int),1,fp);
  fread(&(s->pbits),sizeof(int),1,fp);
  fread(&(s->mask),sizeof(long),1,fp);
  fread(&(s->check),sizeof(long),1,fp);

  fread(&(s->nseq),sizeof(int),1,fp);
  s->ids = (char **)malloc(s->nseq*sizeof(char*));
  for (i=0; i<s->nseq;++i) s->ids[i] = readString(fp);
  s->seqTermOrder = (int *)malloc(s->nseq*sizeof(int));
  fread(s->seqTermOrder,sizeof(int),s->nseq,fp);  
  s->seqlengths = (IndexType *)malloc(s->nseq*sizeof(IndexType));
  fread(s->seqlengths,sizeof(IndexType),s->nseq,fp);

  //  s->sa = (uchar *)malloc(s->ncheck*s->nbytes*sizeof(uchar));
  //  fread(s->sa,sizeof(uchar),s->ncheck*s->nbytes,fp);

  s->sa = NULL;
  s->maxlength=0;
  s->hash=NULL;
  s->hash_step=0;

  return s;
}




/* Read SA  */
void read_suffixArray_body(suffixArray *s, FILE *fp) {
  s->sa = (uchar *)malloc(s->ncheck*s->nbytes*sizeof(uchar));
  fread(s->sa,sizeof(uchar),s->ncheck*s->nbytes,fp);
}



void write_suffixArray(suffixArray *s, FILE *fp) {
  write_suffixArray_header(s,fp);
  fwrite(s->sa,sizeof(uchar),s->ncheck*s->nbytes,fp);
}

