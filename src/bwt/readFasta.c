/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

/*
  Functions to read a fasta file into one long sequence
  (including reverse complement). Using a linked list of SEQstruct
  (defined in sequence.h)

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

// #define TEST

/* Return the size of a file in bytes */
static long file_size(FILE *fp) {
  long size;
  fseek(fp, 0L, SEEK_END);
  size = ftell(fp);
  rewind(fp);
  return size;
}



static inline void test_seq_alloc(char *seq, char *seqend) {
  if (seq>=seqend) {fprintf(stderr,"Not enough memory allocated for sequence\n"); exit(1); }
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
  char *idline = (char *)malloc((MaxIDlen+2)*sizeof(char));
  char *seqend = seq+length;


  /* Absorb until "\n>" (in beginning of file) */
  do c=getc(fp); while (c!='>' && c!=EOF);
  if (c==EOF) return NULL;

  /* Alloc sequence structure to hold anchor and hold total length of
     allocation. */
  ret=ss=alloc_SEQstruct();
  ret->start = seq;
  ret->sort_order = 0;

  /* First char is set to term */
  *seq++=term;

  while (c!=EOF) {
    /* Now at ">" in beginning of line. Read id line */
    filepos=ftell(fp);
    if ( !fgets(idline,MaxIDlen,fp) ) {
      fprintf(stderr,"read_fasta_forward: Should read an identifier, but didn't\n");
      break;
      /* free(idline); return NULL; */
    }
    char *cptr=idline+strlen(idline)-1;
    if (*cptr!='\n') { do c=getc(fp); while (c!='\n' && c!=EOF); } /* Read to EOL */
    if (c==EOF) { break; }
    *cptr = '\0';

    /* Alloc sequence structure and copy ID */
    ss->next = alloc_SEQstruct();
    ss=ss->next;
    ss->id = strdup(idline);
    cptr=ss->id;
    while (*cptr && !isblank(*cptr)) ++cptr;
    if (*cptr) ss->descr = cptr+1;
    *cptr=0;
    ss->start = seq;
    ss->id_filepos = filepos;
    ss->seq_filepos = ftell(fp);
    ret->sort_order += 1;

    /* Read sequence */
    int lastc=0;
    c=getc(fp);
    while ( c!=EOF ) {
      // fprintf(stderr,"%c",c);
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


/* Space must be allocated and pointed to by ss->seq with length ss->len */
/* static */ void revcomp_seq(SEQstruct *ss, char *translate) {
  char *f;
  char *r;
  f = ss->start;
  r = ss->start+2*ss->len-1;
  while (f<r) { *r-- = translate[(int)(*f++)]; }
  ss->len *= 2;
}



/* static */ void set_positions(SEQstruct *ss, int revcomp) {
  long L;
  char *base = ss->start;

  if (revcomp) {
    L = 2*ss->len-1;
    while (ss->next) {
      ss=ss->next;
      ss->fpos = ss->start - base;
      ss->rpos = L-ss->fpos;
    }
  }
  else {
    while (ss->next) {
      ss=ss->next;
      ss->fpos = ss->start - base;
    }
  }
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
  seq = (char *)malloc((totlen+padding)*sizeof(char));

  if (!seq) {
    fprintf(stderr,"readFasta: Failed to alloc seq of length %ld\n",totlen);
    return NULL;
  }
  // fprintf(stderr,"readFasta: Alloc seq of length %ld\n",totlen);

  ss = read_fasta_forward(fp, seq, totlen, transtab, term, padding);

  /* ss->start now is = seq */

  if (complement) {
    al = strlen(complement);
    revtrans = (char*)malloc(128*sizeof(char));
    /* OK. This is not trivial
       Original alphabet is already translated to integers,
               e.g. for "ACGT" A->0, C->1, G->2, T->3
       We want "ACGT" -> "TGCA", so 0->3 (A->T), 1->2 (C->G), etc
       All other chars stay the same
     */
    for (i=0;i<128;++i) revtrans[i]=i;
    for (i=0; i<al;++i) revtrans[i] = transtab[(int)complement[i]];
    revcomp_seq(ss,revtrans);
    free(revtrans);
  }

  /* set fpos and rpos */
  /* Down-size to actual size */
  tmp = ss->start;
  if (complement) {
    set_positions(ss, 1);
    ss->start = (char *)realloc((void*)(ss->start),2*ss->len+padding+1);
  }
  else {
    set_positions(ss, 0);
    ss->start = (char *)realloc((void*)(ss->start),ss->len+1);
  }

  /* If down-sizing by realloc actually copies array: */
  if (tmp != ss->start) {
    sst=ss;
    while (sst->next) { sst=sst->next; sst->start=ss->start+sst->fpos; }
  }

  // for (i=0;i<padding;++i) *(ss->start+ss->len+i) = term;

  // fprintf(stderr,"readFasta returning\n");

  return ss;
}







/* Makes a translation table from an alphabet to a translation, so
        table[alphabet[i]] = translation[i]

   If translation==NULL, the letter alphabet[i] is translated to i.

   Sequences are changed to upper case ACGT and all non-ACGT are changed to N
   0 terminate strings

   Letters not in alphabet are translated to dummy.
   If dummy==0, dummy is set to the translation of the last char in alphabet
   (assumed to be a "wildcard" character)

   Translation for non-characters is -1.

   casesens !=0, means case sensitive, otherwise case insentitive.

   Returns an array (char *table) of length 128

   The length of the translation table (which may contain zeros)
   HAS to be as long (or longer) than alphabet.

   If 0 is not in the alphabet, it is translated to 0
*/
char *translation_table(char *alphabet, char *translation, char dummy, int casesens) {
  int i, l, freetrans=0;
  char *table = (char*)malloc(128*sizeof(char));

  l = strlen(alphabet);

  if (translation==0) {
    translation = (char *)malloc(l*sizeof(char));
    for (i=0; i<l; ++i) translation[i]=i;
    freetrans=1;
  }
  if (dummy==0) dummy = translation[l-1];

  table[0]=0;
  for (i=1; i<128; ++i) { table[i]=-1; if (isalpha((char)i)) table[i]=dummy; }

  if (!casesens) {
    for (i=0; i<l; ++i) {
      table[toupper(alphabet[i])]=translation[i];
      table[tolower(alphabet[i])]=translation[i];
    }
  }
  else {
    for (i=0; i<l; ++i) table[alphabet[i]] = translation[i];
  }

  if (freetrans) free(translation);

  return table;
}



/*
  First long is the length of the concatenated sequence.
 */
void writeSequenceFile(SEQstruct *base, char *alph, FILE *fp) {
  SEQstruct *ss;
  int n, l;
  long x[3];

  /* Write sequence (incl length) */
  fwrite_n_bytes(base->start,base->len,fp);

  /* Write alphabet  */
  fwrite_n_bytes(alph,-1,fp);

  /* Count number of seqs and write */
  n=0; ss = base;
  while (ss->next) { n++; ss = ss->next; }
  fwrite((void *)&n,sizeof(int),1,fp);

  ss=base;
  while (ss->next) {
    ss = ss->next;
    fwrite_n_bytes(ss->id,-1,fp);
    /* write length, fpos, rpos */
    x[0]=ss->len;
    x[1]=ss->fpos;
    x[2]=ss->rpos;
    fwrite((void *)x,sizeof(long),3,fp);
  }
}



/*
  NB: Alphabet is returned in base->id!!!

  Reads from current position in file (after reading the sequence)
  Use readIndex, if you want to seek to index.
*/
static SEQstruct *readIndex_noSeek(FILE *fp) {
  SEQstruct *ss, *base;
  int k;
  long n,l;
  long x[3];

  /* Read alphabet */
  base = alloc_SEQstruct();
  base->id = fread_n_bytes(&n,fp);
  fprintf(stderr,"readIndex_noseek: Alph %ld >>%s<<\n",n,base->id);

  /* Read number of seqs: */
  fread((void *)&k,sizeof(int),1,fp);

  ss=base;
  while (k>0) {
    ss->next=alloc_SEQstruct();
    ss = ss->next;
    /* read id */
    ss->id = fread_n_bytes(&n,fp);
    /* read length, fpos, rpos */
    fread((void *)x,sizeof(long),3,fp);
    ss->len=x[0];
    ss->fpos=x[1];
    ss->rpos=x[2];
    ss->start = NULL;
    --k;
  }

  return base;
}


/*
  Assumes that file is rewound and seeks to beginning of index
*/
SEQstruct *readIndex(FILE *fp) {
  SEQstruct *ss, *base;
  long l;

  /* Read seq length */
  fread((void *)&l,sizeof(long),1,fp);
  fprintf(stderr,"readIndex: Length of seq %ld\n",l);
  /* Go to end of seq */
  fseek(fp,l,SEEK_CUR);

  fprintf(stderr,"readIndex: At start of index\n");

  base = readIndex_noSeek(fp);

  fprintf(stderr,"readIndex: Done\n");

  return base;
}




uchar *readSeqOnly(long *len, FILE *fp, int padding, char term) {
  uchar *c;
  int i;
  c = fread_n_bytes_plus(len,fp, padding);
  for (i=0;i<padding;++i) *(c+*len+i) = term;
  return c;
}



SEQstruct *readSeq(FILE *fp, int padding, char term) {
  SEQstruct *base, *ss;
  long seqlen;
  uchar *seq;

  seq = readSeqOnly(&seqlen, fp, padding, term);

  base = readIndex_noSeek(fp);

  /* Set the ->start pointer to beginning of sequence (is it necessary??) */
  ss=base;
  while (ss->next) {
    ss = ss->next;
    ss->start = base->start+ss->fpos;
  }

  return base;
}




/*
  To get from a position in the concatenated sequence to an actual sequence position,
  a hashing scheme is introduced. For K = (int)(i/2^Estep) = i>>Estep, there is a
  pointer to the SEQstruct.
  Estep is the 2-exponent, and it should be around log2(length of shortest seq).
*/
SEQstruct **make_seq_hash(SEQstruct *base, int Estep) {
  SEQstruct *ss,**spt;
  int i,j, limit;

  i = (base->len)>>Estep;
  spt = (SEQstruct **)malloc(i*sizeof(SEQstruct *));

  ss=base;
  i=0;
  j=base->len -1;
  while (ss->next) {
    ss = ss->next;
    limit = (ss->fpos+ss->len-1)>>Estep;
    while (i<limit) spt[i++] = ss;
    limit = (ss->rpos-ss->len+1)>>Estep;
    while (j>=limit) spt[j--] = ss;
  }
  return spt;
}



/*
  Get pointer to sequence for a given position in concatenated string.
*/
SEQstruct *lookupSeq(IndexType n, SEQstruct **hash, int Estep) {
  int k;

  k = n>>Estep;
  while ( n >= hash[k]->fpos + hash[k]->len ) ++k;
  return hash[k];
}


