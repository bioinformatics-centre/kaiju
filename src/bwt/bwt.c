/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>

#include "common.h"
#include "bwt.h"
#include "fmi.h"



/***********************************************************************
 *
 * Suffix array functions and functions to search FMI
 *
 ***********************************************************************/



/***********************************************
 *
 * Bit magic
 *
 ***********************************************/



static inline void suffixArray_set_masks(BWT *s) {
  s->mask = (1<<(s->pbits))-1;
  s->check = (1<<(s->chpt_exp))-1;
}


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
static inline void suffixArray_decode_number(int *nseq, long *pos, int k, BWT *s) {
  long val = uchar2long( (s->sa + k * s->nbytes), s->nbytes);
  *nseq = val>>s->pbits;
  *pos = val & s->mask;
}


/* Encode long in n bytes (assumes that number can fit in n bytes) */
static inline void long2uchar(long k, uchar *c, int n) {
  while ( n-- >0 ) { c[n] = k; k = k>>8; }
}


/*
  For SA entry k, encode seq no. (nseq) and position within (pos) in
  s->nbytes bytes
*/
static inline void suffixArray_encode_number(int nseq, long pos, int k, BWT *s) {
  long val = nseq;
  val = (val << s->pbits)  + pos;
  long2uchar(val, (s->sa + k * s->nbytes), s->nbytes);
}


/* Returns the bits needed to encode a number (without using log2) */
static int bitsNeeded(long k) { int i=0; while (k>>i) ++i; return i; }



/***********************************************
 *
 * I/O stuff
 *
 ***********************************************/



void translate2numbers_ptr(uchar *s, const IndexType slen, uchar * trans) {
	IndexType k;
	for (k=0;k<slen;++k) s[k]=trans[s[k]];
}


void translate2numbers(uchar *s, const IndexType slen, char *alphabet, int alen) {
uchar *trans;
int i;
IndexType k;
trans = (uchar *)calloc(128,sizeof(uchar));
for (i=0; i<alen;++i) trans[(int)alphabet[i]]=i;
for (k=0;k<slen;++k) s[k]=trans[s[k]];
free(trans);
}



/*
  Read BWT header from file (made by mkbwt) and resturn in BWT struct
*/
BWT *read_BWT_header(FILE *bwtfile) {
  BWT *s=(BWT *)calloc(sizeof(BWT),1);
  const int linelen=1000;
  char line[linelen], dummy[linelen];
  int i, k, l;

  // Read header stuff
  fprintf(stderr,"Reading BWT header ... ");

  fgets(line,linelen-1,bwtfile);
  sscanf(line,"%s %ld",dummy,&(s->len));
  fgets(line,linelen-1,bwtfile);
  sscanf(line,"%s %d",dummy,&(s->alen));
  s->alphabet = (char *)calloc(s->alen+1,sizeof(char));
  fgets(line,linelen-1,bwtfile);
  sscanf(line,"%s %s",dummy,s->alphabet);
  fgets(line,linelen-1,bwtfile);
  sscanf(line,"%s %d",dummy,&(s->nseq));

  fprintf(stderr,"DONE\n");

  // Read ids
  fprintf(stderr,"Reading BWT seq info ... ");
  s->ids = (char **)calloc(s->nseq,sizeof(char*));
  s->sortorder = (int *)malloc(s->nseq*sizeof(int));
  s->seqlengths = (IndexType *)malloc(s->nseq*sizeof(IndexType));
  s->maxlength = 0;
  for (i=0; i<s->nseq; ++i) {
    fgets(line,linelen-1,bwtfile);
    sscanf(line,"%s %d %d",dummy,&l,&k);
    s->ids[i] = strdup(dummy);
    s->seqlengths[i]=l;
    if (l>s->maxlength) s->maxlength=l;
    s->sortorder[i]=k;
  }
  fprintf(stderr,"DONE\n");

  return s;
}

/*
  Read BWT itself from file (made by mkbwt)
  if alphabet != NULL, translate the sequence (otherwise assume it is translated)
*/
uchar *read_BWT(IndexType len, int alen, char *alphabet, FILE *bwtfile) {
  uchar *bwt;
  const int linelen=1000;
  char line[linelen];
  int translate;

  // First line is "#BWT translate=yes" or "#BWT translate=no" or just "#BWT"
  fgets(line,linelen-1,bwtfile);
  if (strcmp("#BWT translate=no",line)==0) translate = 0;
  else translate=1;

  bwt = (uchar *)malloc(len*sizeof(uchar));
  fprintf(stderr,"Reading BWT itself ... ");
  fread(bwt, sizeof(char), len, bwtfile);
  fprintf(stderr,"DONE\n");

  // Translate BWT from letters to numbers
  if (translate) {
    fprintf(stderr,"Translating BWT ... ");
    translate2numbers(bwt, len, alphabet, alen);
    fprintf(stderr,"DONE\n");
  }

  return bwt;
}


/* Assumes that length of string is <256 - otherwise truncate */
static inline void writeString(char *str, int l, FILE *fp) {
  uchar ul;
  if (!l) {
    l = strlen(str);
    if (l>255) l=255;
  }
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




/* Write the suffix array in file (binary) */
void write_BWT_bin(const BWT *s, FILE *fp) {
  int i;
  fwrite(&(s->len),sizeof(IndexType),1,fp);
  fwrite(&(s->ncheck),sizeof(IndexType),1,fp);
  fwrite(&(s->chpt_exp),sizeof(int),1,fp);
  fwrite(&(s->nbytes),sizeof(int),1,fp);
  fwrite(&(s->sbits),sizeof(int),1,fp);
  fwrite(&(s->pbits),sizeof(int),1,fp);
  fwrite(&(s->mask),sizeof(long),1,fp);
  fwrite(&(s->check),sizeof(long),1,fp);
  fwrite(s->sa,sizeof(uchar),s->ncheck*s->nbytes,fp);

  fwrite(&(s->alen),sizeof(int),1,fp);
  fwrite(s->alphabet,sizeof(char),s->alen,fp);

  fwrite(&(s->nseq),sizeof(int),1,fp);
  for (i=0; i<s->nseq;++i) writeString(s->ids[i], 0, fp);
  fwrite(s->sortorder,sizeof(int),s->nseq,fp);
  fwrite(s->seqlengths,sizeof(IndexType),s->nseq,fp);
  fwrite(&(s->maxlength),sizeof(IndexType),1,fp);
  fwrite(s->seq,sizeof(int),s->nseq,fp);
}



BWT *read_BWT_bin(FILE *fp) {
  int i;
  BWT *s = (BWT *)malloc(sizeof(BWT));

  fread(&(s->len),sizeof(IndexType),1,fp);
  fread(&(s->ncheck),sizeof(IndexType),1,fp);
  fread(&(s->chpt_exp),sizeof(int),1,fp);
  fread(&(s->nbytes),sizeof(int),1,fp);
  fread(&(s->sbits),sizeof(int),1,fp);
  fread(&(s->pbits),sizeof(int),1,fp);
  fread(&(s->mask),sizeof(long),1,fp);
  fread(&(s->check),sizeof(long),1,fp);
  s->sa = (uchar *)malloc(s->ncheck*s->nbytes*sizeof(uchar));
  fread(s->sa,sizeof(uchar),s->ncheck*s->nbytes,fp);

  fread(&(s->alen),sizeof(int),1,fp);
  s->alphabet = (char *)malloc((s->alen+1)*sizeof(char));
  fread(s->alphabet,sizeof(char),s->alen,fp);
  s->alphabet[s->alen] = '\0';

  fread(&(s->nseq),sizeof(int),1,fp);
  s->ids = (char **)malloc(s->nseq*sizeof(char*));
  for (i=0; i<s->nseq;++i) s->ids[i] = readString(fp);
  s->sortorder=(int*)malloc(sizeof(int)*s->nseq);
  fread(s->sortorder,sizeof(int),s->nseq,fp);
  s->seqlengths = (IndexType*)malloc(sizeof(IndexType)*s->nseq);
  fread(s->seqlengths,sizeof(IndexType),s->nseq,fp);
  fread(&(s->maxlength),sizeof(IndexType),1,fp);
  s->seq=(int*)malloc(sizeof(int)*s->nseq);
  fread(s->seq,sizeof(int),s->nseq,fp);

  return s;
}




/***********************************************
 *
 * Build SA
 *
 ***********************************************/





/* Initialize suffix array structure to hold SA checkpoint values
   Call after reading BWT
*/
void init_suffixArray(BWT *s, FMI *f, int chpt_exp) {
  s->f = f;
  s->chpt_exp = chpt_exp;
  s->ncheck = (s->len>>chpt_exp)+1;
  s->sbits = bitsNeeded(s->nseq);
  s->pbits = bitsNeeded(s->maxlength);
  s->nbytes = (7+s->sbits+s->pbits)/8;
  suffixArray_set_masks(s);

  s->sa = (uchar*)malloc(s->ncheck*s->nbytes*sizeof(uchar));
  s->seq = (int*)malloc(s->nseq*sizeof(int));

  // For parallel processing
  s->interval = 1+s->nseq/100;
  s->next_start=0;
  pthread_mutex_init(&(s->lock), NULL);
}



/*
  Fill suffix array from BWT in checkpoints
  Assumes that suffix array already contains sequence information
  Can be done in parallel, and called with from/to variables set
  Call after init_suffixArray
*/
void fill_suffixArray(int from, int to, BWT *s) {
  IndexType k, l;
  uchar c;
  int i, all=0;
  FMI *fmi = s->f;

  if (to==0) {from=0; to=s->nseq; all=1; }

  // Run through all sequences
  int ten_percent = (int)(s->nseq/10);
  for (i=from; i<to; ++i) {
    if ( all && (i+1)%ten_percent==0 ) fprintf(stderr,"%d%% ... ",10*(int)((i+1)/ten_percent));
    l = s->seqlengths[i];
    k = (IndexType)(s->sortorder[i]);
    while (l-->0) {
      k = FMindexCurrent(fmi, &c, k);
      if ( !(k&s->check) ) {
	suffixArray_encode_number(i, l, k>>s->chpt_exp, s);
      }
    }
    k = FMindexCurrent(fmi, &c, k);
    if (c!=0) fprintf(stderr,"fill_suffixArray: SOMETHING WRONG sequence %d %s %ld %c\n",i, s->ids[i], l, s->alphabet[c]);
    s->seq[k]=i;   // This is the k'th termination in BWT, which then corresponds to sequence i
  }
}


/* Parallel filling of SA */
void *fill_suffixArray_parallel(BWT *s) {
  int from,to;

  while (1) {
    pthread_mutex_lock(&(s->lock));
    if (s->next_start>=s->nseq) { pthread_mutex_unlock(&(s->lock)); break; }
    from=s->next_start;
    s->next_start = to = from+s->interval;
    if ( to>s->nseq) s->next_start = to = s->nseq;
    pthread_mutex_unlock(&(s->lock));
    fill_suffixArray(from, to, s);
  }
  return NULL;
}




/***********************************************
 *
 * Query SA && FMI
 *
 ***********************************************/



/* Find suffix for suffix number i
   Return sequence number in *iseq and position in *pos
*/
void get_suffix(FMI *fmi, BWT *s, IndexType i, int *iseq, IndexType *pos) {
  IndexType k=0;
  uchar c=1;

  while ( c && (i & s->check) ) {
    i = FMindexCurrent(fmi,&c,i);
    ++k;
  }

  if (c) { suffixArray_decode_number(iseq, pos, i>>s->chpt_exp, s); *pos += k; }
  else { *iseq = s->seq[i]; *pos=k-1; }

}



/*
  Reconstruct sequence number snum (according to the original order of the sequence file)
*/
uchar *retrieve_seq(int snum, BWT *s, FMI *f) {
  IndexType l, k;
  uchar *seq, *cur;

  l = s->seqlengths[snum];
  k = (IndexType)(s->sortorder[snum]);
  seq = (uchar *)malloc((l+2)*sizeof(uchar));

  cur = seq+l+1;
  *cur-- = 0; *cur-- = 0;
  while (cur>=seq) {
    k = FMindexCurrent(f, cur, k);
    --cur;
  }
  return seq;
}




/* Find whole (initial) suffix interval for letter ct */
static inline IndexType InitialSI(FMI *f, uchar ct, IndexType *si) {
  IndexType r = f->N1-1;
  si[0]=f->index1[r][ct];
  if (ct<f->alen-1) si[1]=f->index1[r][ct+1];
  else si[1]=f->bwtlen;
  return si[1]-si[0]+1;
}



/* Finds the suffix interval for letter ct
   If newsi==NULL and an SI is found, it OVERWRITES values in si
   Note that the actual SI is from si[0] to si[1]-1
*/
IndexType UpdateSI(FMI *f, uchar ct, IndexType *si, IndexType *newsi) {
  IndexType nsi[2];

  if (si[0]>si[1]) return 0;

  nsi[0] = FMindex(f, ct, si[0]);
  nsi[1] = FMindex(f, ct, si[1]);

  if (nsi[0]>=nsi[1]) return 0;

  if (!newsi) newsi=si;
  newsi[0] = nsi[0];
  newsi[1] = nsi[1];

  return nsi[1] - nsi[0];
}



static SI *alloc_SI(IndexType *si, int query_pos, int query_len){
  SI *r = (SI *)malloc(sizeof(SI));
  r->start = si[0];
  r->len=(int)(si[1]-si[0]);
  r->qi = query_pos;
  r->ql = query_len;
  r->count = 0;
  r->next = NULL;
  r->samelen = NULL;
  return r;
}



void recursive_free_SI(SI *si) {
  if (!si) return;
  if (si->next) recursive_free_SI(si->next);
  if (si->samelen) recursive_free_SI(si->samelen);
  free(si);
}


/*
  Free matches of each length until there are at least max (we would get <max
  if more were freed).
  Returns min length of retained matches
*/
static inline int free_until_max_SI(SI *si, int max) {
  int n;
  SI *cur;
  if (!si || si->count<=max ) return 0;
  n = si->count;
  cur = si;
  while (cur->next && n - cur->next->count < max) cur=cur->next;
  // Now cur->next==NULL OR totalcount>=max for cur
  if (cur->next) {
    n = cur->next->count;
    recursive_free_SI(cur->next);
    cur->next=NULL;
    while (si) {si->count -= n; si=si->next; }
  }
  return cur->ql;
}



// Brute force stupid sorting (assuming short lists)
static SI *insert_SI_sorted(SI *base, SI *new) {
  SI *tmp;
  new->count=new->len;
  if (base==NULL) { return new; }
  if (base->ql<new->ql) {
    new->next=base;
    new->count += base->count;
    return new;
  }
  tmp=base;
  while (tmp->next && tmp->next->ql >= new->ql) {
    tmp->count += new->len;
    tmp=tmp->next;
  }
  tmp->count +=new->len;
  // Now tmp is >= new AND (tmp->next<new OR tmp->next==NULL)
  if (tmp->ql==new->ql) {
    new->samelen=tmp->samelen;
    if (tmp->samelen) new->count += tmp->samelen->count;
    tmp->samelen=new;
  }
  else {
    new->next=tmp->next;
    if (tmp->next) new->count += tmp->next->count;
    tmp->next=new;
  }
  return base;
}




/* Find maximal matches of str of length L in a sorted linked list
   Returns null if there are no matches
   If max_matches==0, not limit imposed
*/
SI *maxMatches(FMI *f, char *str, int len, int L, int max_matches) {
  SI *first=NULL, *cur=NULL;
  IndexType si[2], l;
  int i, j, k;

  // Go through the sequence from the back
  for (j=len-1; j>=L-1; --j) {
    i=j;
    InitialSI(f, str[i], si);
    // Extend backward
    while ( i-- > 0 ) {
      if ( UpdateSI(f, str[i], si, NULL) == 0) break;
    }
    i+=1;
    l = j-i+1;
    if (l>=L) {
      // If the begin of the match (i) equals the the previous, it is within previous match
      if ( !cur || i < cur->qi ) {
	cur = alloc_SI(si, i, l);
	first = insert_SI_sorted(first, cur);
	// If max_matches is set, check to see if max is reached and reset L
	if (max_matches>0) {
	  k = free_until_max_SI(first, max_matches);
	  if (k>L) L=k;
    // The latest si may be freed if too short - then set it =NULL
    if (l<k) cur=NULL;
	}
      }
    }
    // If the last match reached beginning of sequence, no need to continue
    if (i<=1) break;
  }

  return first;
}

//PETER:
SI *maxMatches_fromlastposonly(FMI *f, char *str, int len, int L, int max_matches) {
  SI *first=NULL, *cur=NULL;
  IndexType si[2], l;
  int i, j, k;

  // Go through the sequence from the back
  j=len-1;
  //for (j=len-1; j>=L-1; --j) {
    i=j;
    InitialSI(f, str[i], si);
    // Extend backward
    while ( i-- > 0 ) {
      if ( UpdateSI(f, str[i], si, NULL) == 0) break;
    }
    i+=1;
    l = j-i+1;
    if (l>=L) {
      // If the begin of the match (i) equals the the previous, it is within previous match
      if ( !cur || i < cur->qi ) {
	cur = alloc_SI(si, i, l);
	first = insert_SI_sorted(first, cur);
	// If max_matches is set, check to see if max is reached and reset L
	if (max_matches>0) {
	  k = free_until_max_SI(first, max_matches);
	  if (k>L) L=k;
    // The latest si may be freed if too short - then set it =NULL
    if (l<k) cur=NULL;
	}
      }
    }
    // If the last match reached beginning of sequence, no need to continue
    //if (i<=1) break;
  //}

  return first;
}

SI *maxMatches_withStart(FMI *f, char *str, int len, int L, int max_matches, IndexType si0, IndexType si1, int offset) {
  SI *first=NULL, *cur=NULL;
  IndexType si[2], l;
  int i, j, k;
  si[0] = si0;
  si[1] = si1;

  // Go through the sequence from the back
  j=len-1;
  //for (j=len-1; j>=L-1; --j) {
    i=j-offset+1;
    //InitialSI(f, str[i], si);
    // Extend backward
    while ( i-- > 0 ) {
      if ( UpdateSI(f, str[i], si, NULL) == 0) break;
    }
    i+=1;
    l = j-i+1;
    if (l>=L) {
      // If the begin of the match (i) equals the the previous, it is within previous match
      if ( !cur || i < cur->qi ) {
	cur = alloc_SI(si, i, l);
	first = insert_SI_sorted(first, cur);
	// If max_matches is set, check to see if max is reached and reset L
	if (max_matches>0) {
	  k = free_until_max_SI(first, max_matches);
	  if (k>L) L=k;
    // The latest si may be freed if too short - then set it =NULL
    if (l<k) cur=NULL;
	}
      }
    }
    // If the last match reached beginning of sequence, no need to continue
    //if (i<=1) break;
  //}

  return first;
}



static SI *maxMatches_old (FMI *f, char *str, int len, int L) {
  SI *first=NULL, *cur=NULL;
  IndexType si[2], l;
  int i, j;

  // Go through the sequence from the back
  for (j=len-1; j>=L-1; --j) {
    i=j;
    InitialSI(f, str[i], si);
    // Extend backward
    while ( i-- > 0 ) {
      if ( UpdateSI(f, str[i], si, NULL) == 0) break;
    }
    i+=1;
    l = j-i+1;
    if (l>=L) {
      // If the begin of the match (i) equals the the previous, it is within previous match
      if ( !cur || i < cur->qi ) {
	// fprintf(stderr,"MATCH %d %d-%d %ld\n",(int)(si[1]-si[0]), i,j, l);
	cur = alloc_SI(si, i, l);
	first = insert_SI_sorted(first, cur);
      }
    }
  }

  return first;
}


/* Find the all matches of str of length L. Return in a linked list
   Returns null if there are no matches
*/
SI *allMatches(FMI *f, char *str, int len, int L) {
  SI *first=NULL, *cur=NULL, *tmp;
  IndexType si[2], l;
  int i, j;

  for (i=0; i<=len-L; ++i) {
    j = i+L-1;
    l = InitialSI(f, str[j], si);
    while ( j-- > i && l>0 ) l = UpdateSI(f, str[j], si, NULL);
    if (l>0) {
      tmp=cur;
      cur = alloc_SI(si, i,L);
      if (first) tmp->next = cur;
      else first = cur;
    }
  }

  return first;
}




void SIprintQmatch(SI *si, char *query, char *alphabet, FILE *fp) {
  int i;
  for (i=0; i<si->ql; ++i) fprintf(fp,"%c",alphabet[query[si->qi+i]]);
}


/* Find matches of min length L and print SI length for each */
void printMatches(FMI *f, char *str, int len, int L, char *alphabet) {
  IndexType si[2], l1, l2;
  double x;
  int i, j;

  for (i=L-1; i<len; ++i) {
    j = i;
    l1 = InitialSI(f, str[j], si);
    while ( j-- > 0 && l1>0 ) {
      l2=l1;
      l1 = UpdateSI(f, str[j], si, NULL);
      if (i-j >= L) printf("%d %d %c %d %d\n",j,i,alphabet[str[j]],(int)l1,(int)l2);
    }
  }
}

