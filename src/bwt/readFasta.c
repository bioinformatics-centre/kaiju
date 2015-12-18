/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
/*
  Functions to read a fasta file into one long sequence
  (including reverse complement if specified). Using a linked list of
  SEQstruct (defined in sequence.h)

  The first struct points to the long sequence (->start), has the
  alphabet in the ->id and the total length of the whole sequence
  (->len). The following contains the individual sequences.

  Letters in in the sequence are translated according to an array
  (translation) with an entry for each of the 128 ASCII codes.
  Symbols that are not letters (spaces, newlines etc) should be
  negative.

*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "sequence.h"
#include "readFasta.h"
#include "common.h"

/* Return the size of a file in bytes */
static long file_size(FILE *fp) {
  long size;
  fseek(fp, 0L, SEEK_END);
  size = ftell(fp);
  rewind(fp);
  return size;
}



static inline void test_seq_alloc(char *seq, char *seqend) {
  if (seq>=seqend) ERROR("Not enough memory allocated for sequence",1);
}



/* Assumes that filepointer is at beginning of id (after ">")
   Allocates SEQstruct and returns
   Truncates ID line if longer than MaxIDlen, but always reads to end of line
*/
static SEQstruct *read_id(FILE *fp, char *idline) {
  SEQstruct *ret;
  int c=1;
  char *cptr;
  const int MaxIDlen = 1000;

  if ( !fgets(idline,MaxIDlen,fp) ) {
    fprintf(stderr,"read_id: Should read an identifier, but didn't\n");
    free(idline);
    return NULL;
  }

  /* Absorb until end of line if longer than MaxIDlen */
  cptr=idline+strlen(idline)-1;
  if (*cptr!='\n') { do c=getc(fp); while (c!='\n' && c!=EOF); }
  if (c==EOF) {
    fprintf(stderr,"read_fasta_forward: EOF while reading ID line\n");
    free(idline);
    return NULL;
  }
  *cptr = '\0';

  /* Alloc sequence structure and copy ID */
  ret = alloc_SEQstruct();
  ret->id = strdup(idline);
  /* Separate id and description */
  cptr=ret->id;
  while (*cptr && !isblank(*cptr)) ++cptr;
  if (*cptr) ret->descr = cptr+1;
  *cptr=0;

  return ret;
}



/*
  Read all sequences into one long array and ids into a linked list of IDs.
  string translate holds the char used for each ascii symbol

  Sequences are each terminated by padding term symbols

  Note: seq must point to enough allocated memory!!

  Returns SEQstruct that starts a linked list. In first element:
     ->start        points to beginning of sequence allocation  
     ->len          (long) holds the total length of concattenated sequence
     ->sort_order   (int) holds the total number of sequences read

*/
static SEQstruct *read_fasta_forward(FILE *fp, char *seq, long length, char *translate, char term, int padding) {
  const int MaxIDlen = 10000;
  long filepos;
  int c, i;
  SEQstruct *ret=NULL, *ss;
  char *seqend = seq+length;
  char *idline = (char *)malloc((MaxIDlen+2)*sizeof(char));


  /* Absorb until "\n>" (in beginning of file) */
  do c=getc(fp); while (c!='>' && c!=EOF);
  if (c==EOF) return NULL;

  /* Alloc sequence structure to hold anchor and hold total length of
     allocation. */
  ret=ss=alloc_SEQstruct();
  ret->start = seq;
  ret->sort_order = 0;

  /* Add padding in beginning */
  for (i=0; i<padding; ++i) *seq++ = term;

  while (c!=EOF) {
    /* Now at ">" in beginning of line. Read id line */
    filepos=ftell(fp);
    ss->next = read_id(fp,idline);
    if (ss->next==NULL) ERROR("Failed to read ID",3);
    ss=ss->next;

    ss->start = seq;
    ss->pos = (long)(seq - ret->start);
    ss->id_filepos = filepos;
    ss->seq_filepos = ftell(fp);
    ret->sort_order += 1;

    /* Read sequence */
    int lastc=0;
    c=getc(fp);
    while ( c!=EOF ) {
      /* New fasta entry is starting */
      if (c=='>' && lastc=='\n') {
	ss->len = seq - ss->start;
	// Add padding
	for (i=0; i<padding; ++i) { *seq++ = term; test_seq_alloc(seq,seqend); }
	break;
      }
      if (translate[c]>=0) {
        *seq++ = translate[c];
	test_seq_alloc(seq,seqend);
	ss->len += 1;
      }
      lastc=c;
      c=getc(fp);
    }
  }

  if (seq != ret->start) for (i=0; i<padding; ++i) { *seq++ = term; test_seq_alloc(seq,seqend); }
  ret->len = seq - ret->start;

  free(idline);

  //fprintf(stderr,"read_fasta_forward returning\n");

  return ret;
}



SEQstruct *revcompSEQstruct(SEQstruct *ss, char *s, char *translate) {
  SEQstruct *rr = alloc_SEQstruct();
  long i;
  char *r, *f;

  rr->id = ss->id;
  rr->descr = ss->descr;
  rr->rc = 1;
  rr->len = ss->len;
  rr->start = s;
  f = ss->start;

  r=rr->start+rr->len-1;
  for (i=0; i<ss->len; ++i) *r-- = translate[(int)(*f++)];

  return rr;
}



static void revcomp_all(SEQstruct *fbase, char *translate, char term, int padding) {
  int i;
  SEQstruct *rbase=NULL, *ss, *rr;
  char *s=fbase->start+fbase->len;

  /* Reverse complement each */
  ss = fbase;
  while (ss->next) {
    ss=ss->next;
    if (!rbase) rr=rbase=revcompSEQstruct(ss,s, translate);
    else { rr->next=revcompSEQstruct(ss,s,translate); rr=rr->next; }
    rr->pos = (long)(s-fbase->start);
    s+=rr->len;
    for (i=0; i<padding; ++i) *s++ = term;
    fbase->sort_order += 1;  // Holds number of seqs!
  }
  ss->next = rbase;
  fbase->len = (long)(s-fbase->start);
}




/*
  Note that complement is the complement alphabet (if alph is ACGTN, complement is TGCAN)
  If complement==NULL, it is not reverse transcribed

  If length is not given, the file size is used. Otherwise length MUST be at least as long
  as the concatenated string!!

  Sequence is padded with "padding" term chars
*/

SEQstruct *readFasta(FILE *fp, long length, char *transtab, char *complement, char term, int padding) {
  char *seq, *tmp;
  char *revtrans;
  int i, al;
  long totlen;
  SEQstruct *ss, *sst;

  /* alloc space for long sequence */
  if (length) totlen = length;
  else totlen = file_size(fp);
  if (complement) totlen *= 2;
  seq = (char *)malloc(totlen*sizeof(char));

  if (!seq) {
    fprintf(stderr,"readFasta: Failed to alloc seq of length %ld\n",totlen);
    return NULL;
  }

  ss = read_fasta_forward(fp, seq, totlen, transtab, term, padding);

  if (!ss) ERROR("readFasta: No sequences read",5);

  if (complement) {
    if (totlen<2*ss->len) ERROR("readFasta: Not enough space allocated for reverse complement",2);

    al = strlen(complement);
    revtrans = (char*)malloc(128*sizeof(char));
    /* OK. This is not trivial
       Original alphabet is already translated to integers,
               e.g. for "*ACGT" A->1, C->2, G->3, T->4
       We want "*ACGT" -> "*TGCA", so 1->4 (A->T), 2->3 (C->G), etc
       All other chars stay the same
     */
    for (i=0;i<128;++i) revtrans[i] = i;
    for (i=0; i<al;++i) revtrans[i] = transtab[(int)complement[i]];

    revcomp_all(ss, revtrans, term, padding);
    free(revtrans);
  }

  /* Down-size to actual size */
  tmp = ss->start;
  ss->start = (char *)realloc((void*)tmp,ss->len);

  /* If down-sizing by realloc actually copies to new array: */
  if (tmp != ss->start) {
    sst=ss;
    while (sst->next) { sst=sst->next; sst->start=ss->start+sst->pos; }
  }

  // fprintf(stderr,"readFasta returning\n");

  return ss;
}









