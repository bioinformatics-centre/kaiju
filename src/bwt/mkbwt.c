/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

/*
  Finds the BWT of a set of sequences (concatenated)

  Input: fasta file

  ALPHABET
  ========

  Characters not in the alphabet are set to the last char - so if you want to
  include a wildcard, the alphabet should have e.g. X as the last charater.

  For proteins the alphabet might look like this
       ACDEFGHIKLMNPQRSTVWYX
  
  Internally the seq termination character is 0 and this is sorted first.

  OUTPUT
  ======
  The BWT is output in the same alphabet as the sequence, but contains
  termination characters


gcc -g -o sufSort -l pthread sufSort.c

*/

// Use multikey quick sort
#define MKQS

// Use repeat sorting (not necessary with MKQS -?)
#define REPSORT

// DEBUG1 prints info on function calls etc by each worker
// #define DEBUG1

// DEBUG2 prints the sorted seqs to stdout (use only for small sets!)
// Number is the seq length printed
// #define DEBUG2 50


#ifdef DEBUG1
#define DEBUG1LINE(x) x; fflush(stderr);
#else
#define DEBUG1LINE(x) 
#endif

#ifdef DEBUG2
#define DEBUG2LINE(x) x
#else
#define DEBUG2LINE(x) 
#endif




#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <unistd.h>
#include <pthread.h>
#include <time.h>

#include "multikeyqsort.h"
#include "common.h"
#include "sequence.h"
#include "readFasta.h"
#include "mkbwt_vars.h"
#include "suffixArray.h"



/* Global var for setting of number for worker */
int workerNum=0;


/* Status flags */
#define UNINIT 0
#define FILLED 1
#define SORTED 2
#define BWT 3
#define DONE 4
#define SAFILE 6
#define PROCESS 5


typedef struct {
  int wlen;     /* prefix length */
  IndexType start; /* Number of the first suffix */
  int wn;       /* prefix (word) number */
  int jump;     /* =wlen except for words with terminator for which it is the dist to term */
  char *word;   /* prefix (word) */
  char **sa;    /* Suffixes which needs to be sorted (pointers to *seq) */
  char **cur;   /* Current suffix */
  long len;     /* Length of sa to be sorted */
  int alen;     /* Length of alphabet */
  char *alphabet; /* Alphabet */
  long slen;    /* Length of the whole sequence */
  char *seq;    /* Sequence to be sorted */
  char *bwt;    /* BWT for this bucket */
  int status;
} Bucket;




/* Stack of Buckets to be processed by the threads + info needed for threads */
typedef struct {
  int alen;
  char *alphabet;
  int wlen;       /* Word length for sorting */
  long slen;      /* Length of the whole sequence */
  char *seq;      /* Sequence to be sorted */

  int nbuckets;   /* Number of buckets */
  long *bucket_size;  /* Bucket sizes */
  Bucket **b;     /* Pointers to Buckets that needs sorting */

  int sort;       /* Next bucket for sorting */
  int sawrite;    /* Next bucket to write SA checkpoints */
  int bwtwrite;   /* Next bucket to write bwt */
  int nthreads;   /* Number of threads */
  int nfill;      /* Number of buckets to fill in each filling round */
  int nextfill;   /* Number of buckets to fill next */

  suffixArray *sa_struct;

  pthread_mutex_t lock;  /* Lock to be set by the thread while grabbing the next job */
  pthread_mutex_t lock_fill;  /* Lock to be set by the thread while grabbing the next job */

  FILE *bwtfile;  /* Pointer to an open file for writing of BWT */
  FILE *safile;   /* Pointer to an open file for writing SA */

} BucketStack;





/* Functions for calculating word numbers

   Numbering is the standard way - e.g. for DNA:
   0 ***    6 *A*    12 *C*
   1 **A    7 *AA       :
   2 **C    8 *AC
     :        :
   5 **N   11 *AN    17 *CN

   etc...

   Formula: w[wlen-1]*alen^0 + w[wlen-2]*alen^1 + ... + w[0]*alen^(wlen-1)

   Table is made with this, so _letter_numbers_[n][a] = a*alen^n

*/
static int **_letter_numbers_;
static void fill_letter_numbers(int alen, int wlen) {
  int i,power=1, n=wlen;
  _letter_numbers_ = (int **)malloc(n*sizeof(int *));
  while ( n-- >0 ) {
    _letter_numbers_[n] = (int *)malloc(alen*sizeof(int));
    for ( i=0; i<alen; ++i) _letter_numbers_[n][i] = i*power;
    power *= alen;
  }
}


/* NOTE: We are using negative chars for sorting terminal symbols, so we have s[n]>0,
   which means that we always add 0 for a terminating symbol
*/
static inline int word_number(char *s, int n, int alen) {
  int w = 0;
  while (n-- > 0) {
    if (s[n]>0) w += _letter_numbers_[n][ (int)s[n] ];
  }
  return w;
}



/* Word corresponding to number
   l=word length
   al=alphabet length
*/
static int *powers=NULL;
static void free_powers() { if (powers) free(powers); powers=NULL; }
static void fill_powers(int l, int al) {
  int i;
  free_powers();
  powers = (int *)malloc(l*sizeof(int));
  powers[0]=1;
  for (i=1 ; i<l ; ++i) powers[i]=al*powers[i-1];
}


/*
   *w must point to allocated space
   If alphabet is not given, letters are encoded as numbers from 0 to alen-1
   word is terminated by 0
*/
static char *number2word(int n, int alen, char *w, int wlen, char *alphabet) {
  int i=0, nw, l=wlen;
  while (l>0) { nw = n/powers[--l]; w[i++]=nw; n -= nw*powers[l]; }
  w[i]=0;
  if (alphabet) for (i=0; i<wlen; ++i) w[i]=alphabet[w[i]];
  return w;
}


/* Alloc Bucket, but don't alloc suffix array */
Bucket *allocBucket(long slen, char *seq, int wlen, int wn, int alen,
		    char *alphabet, long bucket_size, IndexType sa_start) {
  int i;
  Bucket *bucket = (Bucket *)malloc(sizeof(Bucket));

  bucket->wlen=wlen;
  bucket->start = sa_start;
  bucket->wn = wn;
  bucket->word = (char *)malloc(wlen*sizeof(char)+1);
  number2word(wn,alen,bucket->word,wlen,alphabet);
  bucket->len=bucket_size;
  bucket->sa=NULL;
  bucket->cur=NULL;
  bucket->alen=alen;
  bucket->alphabet=alphabet;
  bucket->slen=slen;
  bucket->seq=seq;
  bucket->bwt=NULL;
  bucket->status=UNINIT;

  /* Check if word contains 0 */
  for (i=0; i<bucket->wlen; ++i) if (bucket->word[i]==alphabet[0]) break;
  bucket->jump = i;

  return bucket;
}





/* Calculate bucket sizes
   Returns an array with bucket sizes
*/
long *bucketSizes(char *seq, long len, int alen, int n, int *nbuckets) {
  char *end;
  int i,l;
  long *b;

  fill_letter_numbers(alen, n);

  /* Number of buckets */
  l = 1;
  for (i=0; i<n; ++i) l *= alen;
  if (l>1e6) {
    fprintf(stderr,"bucketSizes: n=%d is too big for alphabet len=%d .. n^(alen+1)=%d\n",n,alen,l);
    exit(7676);
  }
  *nbuckets=l;

  b=(long *)calloc(l,sizeof(long));

  end = seq+len-1;
  while (seq<end) {
    if ( *seq>0 && *seq<alen) b[word_number(seq,n,alen)] += 1;
    ++seq;
  }

  return b;
}



BucketStack *initBucketStack(int alen, char *alphabet, long slen, char *seq, FILE *bwtfile,
			     FILE *safile, suffixArray *sa_struct, int wlen) {
  int i,j;
  IndexType sa_pos=0;
  BucketStack *bs=(BucketStack *)malloc(sizeof(BucketStack));

  bs->alen = alen;
  bs->alphabet = alphabet;
  bs->slen = slen;
  bs->seq = seq;
  bs->wlen = wlen;
  bs->bwtfile = bwtfile;
  bs->sa_struct=sa_struct;
  bs->safile = safile;

  /* Calculate bucket sizes */
  bs->bucket_size = bucketSizes(seq, slen, alen, wlen, &(bs->nbuckets));

  /* Allocate buckets */
  bs->b = (Bucket **)calloc(bs->nbuckets,sizeof(Bucket *));

  fill_powers(wlen, alen);
  for (i=0; i<bs->nbuckets; ++i) {
    bs->b[i] = allocBucket(slen, seq, wlen, i, alen, alphabet, bs->bucket_size[i], sa_pos);
    sa_pos += bs->bucket_size[i];
  }

  bs->sort = 0;
  bs->sawrite = -1;
  bs->bwtwrite = -1;
  bs->nextfill = 0;

  pthread_mutex_init(&(bs->lock), NULL);
  pthread_mutex_init(&(bs->lock_fill), NULL);

  return bs;
}



#ifdef DEBUG2
/* This is for debugging
   Prints sequence from pointer current to end (or max)
   if current==NULL, icurrent is used instead
*/
void print_seq(char *seq, char *current, long icurrent, long slen, int max, char *alphabet, int alen, FILE *fp) {
  long i;
  if (current) icurrent = current-seq;
  /* Print bwt */
  i = icurrent-1; if (i<0) i=slen-1;
  fprintf(fp,"%c ",alphabet[seq[i]]);
  if (slen>icurrent+max) max=icurrent+max;
  else max=slen;
  for (i=icurrent; i<max; ++i) fprintf(fp,"%c",(seq[i]>=alen?'a'+seq[i]-alen:alphabet[seq[i]]));
  fprintf(fp,"\n");
}
#endif




/* Fill m next buckets */
void fillBuckets(BucketStack *bs, int istart, int m) {
  int i, iend;
  Bucket **bucket=bs->b;
  char *seq, *end;

  /* The interval of word numbers to fill */
  iend = istart+m;
  if (iend > bs->nbuckets) iend = bs->nbuckets;

  /* alloc buckets  */
  for (i=istart; i<iend; ++i) {
    if (bucket[i]->len)
      bucket[i]->sa = bucket[i]->cur = (char **)malloc(bucket[i]->len*sizeof(char *));
    else bucket[i]->sa = bucket[i]->cur = NULL;
  }

  /* Fill buckets */
  seq = bs->seq;
  end = seq+bs->slen-1;
  while (seq<end) {
    if (*seq>0 ) {
      i=word_number(seq,bs->wlen,bs->alen);
      if (i>=istart && i<iend) {
	*(bucket[i]->cur) = seq+bucket[i]->jump;
	bucket[i]->cur += 1;
      }
    }
    ++seq;
  }

  /* Flag that these buckets are now ready for sorting */
  for (i=istart; i<iend; ++i) { bucket[i]->status = FILLED; }

}



/* Lib function strcmp (in string.h) uses UNSIGNED char for comparison */
static inline int compare_strings(const void *str1, const void *str2) {
  char *s1, *s2;
  s1 = *(char **)str1;
  s2 = *(char **)str2;
  while (*s1 && *s2 && *s1==*s2) { ++s1; ++s2; }
  return (int)( *s1 - *s2 );
}




#ifdef REPSORT

/*
  long repeats of the same letter (such as Ns in genomes) are extremely
  slow to sort. This function takes care of such repeats
 */
void repeatSuffixSort(char **s, int l, char a, int jump) {
  int clow, chigh;
  int count;
  char *tmp, **test, **low, **high, *limit;

  // fprintf(stderr,"repeatSuffixSort %d\n",a);

  /* Identify all suffixes starting with letter different from a
     and move to appropriate location in array
     Note that if a suffix starts with a it is ignored, because it'll be handled later by
     "back-propagation" of those suffixes starting with a letter != a
     (here a suffix is the suffix coming after the word that defines the bucket)
  */
  test = low = s;
  high = s+l-1;
  limit = *s;           // Limit is the lowest pointer value encluntered
  while (test<=high) {
    if (*test<limit) limit = *test;
    if (**test<a) { *low = *test; ++low; ++test; }
    else {
      if (**test>a) {
	tmp = *high;
	*high = *test;
	--high;
	*test = tmp;
      }
      else ++test;
    }
  }
  clow = (int)(low-s);
  chigh = (int)(s+l-high-1);
  limit -= jump;

  /* Now sort the high and low intervals */
#ifdef MKQS
  multikeyqsort(s,clow);
  multikeyqsort(s+l-chigh,chigh);
#else
  qsort(s,clow,sizeof(char*),compare_strings);
  qsort(s+l-chigh,chigh,sizeof(char*),compare_strings);
#endif

  if (clow+chigh<l) {
    /* Now "backpropagate" low results */
    count = clow;
    low = s+clow;
    test = s;
    while (count) {
      tmp = *test-jump-1;
      if (tmp>=limit && *tmp==a) { *low = *test-1; ++low; }
      else --count;
      ++test;
    }

    /* Now "backpropagate" high results */
    count = chigh;
    high = s+l-chigh-1;
    test = s+l-1;
    while (count) {
      tmp = *test-jump-1;
      if (tmp>=limit && *tmp==a) { *high = *test-1; --high; }
      else --count;
      --test;
    }
  }

}



/* Check for homopolymers (repeats)
   Returns the character number if it is a hp. Otherwise return -1
*/
static inline int checkHomoPol(char *s, int wlen) {
  int i, j;
  // The suffix comes after the word of the bucket
  s -= wlen;
  i=s[0];
  for (j=1;j<wlen;++j) if (s[j]!=i) { i=-1; break; }
  return i;
}

#endif

/* Sorts a single bucket and return pointers to right location (-wlen) */
void sortBucket(Bucket *b) {
  int i, h;

  if (b->len==0) {b->status=SORTED; return; }

#ifdef REPSORT

  h = checkHomoPol(b->sa[0],b->wlen);
  if (h>0) {
    // fprintf(stderr,"sortBucket: sorting homopolymer %d\n",h);
    repeatSuffixSort(b->sa, b->len, h, b->wlen);
  }
  else
#endif
#ifdef MKQS
  multikeyqsort(b->sa,b->len);
#else
  qsort(b->sa,b->len,sizeof(char*),compare_strings);
#endif

  for (i=0; i<b->len; ++i) b->sa[i]-=b->jump;
  b->status=SORTED;
}



/* Get the BWT for bucket */
void bwtBucket(Bucket *b) {
  int k;

  if (b->len) {
    b->bwt = malloc(b->len);
    for (k=0; k<b->len; ++k) b->bwt[k]=*(b->sa[k]-1);
  }
  b->status = BWT;
}


/* free the SA for bucket */
static inline void free_sa(Bucket *b) {
  int k;

  if (b->len) {
    free(b->sa);
    b->sa=NULL;
  }
}




void bwtWriteBucket(Bucket *b, FILE *bwtfile) {
  int k;

#ifdef DEBUG2
  if (b->len) {
    int i;
    for (i=0; i<b->len; ++i) {
      fprintf(stdout,"DEBUG2%8d %4d ", (int)(b->sa[i]-b->seq), b->wn );
      print_seq(b->seq, b->sa[i], 0, b->slen, DEBUG2, b->alphabet, b->alen, stdout);
    }
  }
  free(b->sa);
  b->sa=NULL;
#endif
  if (b->len) {
    // Write actual bwt
    fwrite(b->bwt,sizeof(char),b->len, bwtfile);
    // Free BWT
    free(b->bwt);
    b->bwt=NULL;
  }
  b->status = DONE;
}





/* This is the worker function that sorts buckets */
/*
  A bucket needs to either be
  - sorted
  - back-propagated (to get BWT)
  - written

  The first two can be done in parallel
  WRITING MUST BE SERIAL
*/
void *BucketSorter(void *x) {
  BucketStack *bs = (BucketStack *)x;
  Bucket *b;
  int i;
  int wn=++workerNum;

  DEBUG1LINE(fprintf(stderr,"Worker %d STARTS\n",wn));

  while( bs->bwtwrite < bs->nbuckets ) {

    pthread_mutex_lock(&(bs->lock));
    /* Skip empty buckets - does not work -? */
    // while (bs->b[bs->write]->len ==0) ++bs->write;

    // While waiting for the unlocking, bs->sort might have been incremented
    if (bs->sort>=bs->nbuckets) {
      pthread_mutex_unlock(&(bs->lock));
      break;
    }


    /* Write SA checkpoints */
    if ( bs->sawrite >=0 && bs->sawrite < bs->nbuckets && bs->b[bs->sawrite]->status == BWT ) {
      b=bs->b[bs->sawrite];
      b->status=PROCESS;
      pthread_mutex_unlock(&(bs->lock));
      DEBUG1LINE(fprintf(stderr,"Worker %d: write SA in file for word %d %s\n",wn,b->wn,b->word));
      write_suffixArray_checkpoints(b->sa, b->start, b->len, bs->sa_struct, bs->safile);
      DEBUG1LINE(fprintf(stderr,"Worker %d: write SA in file for word %d %s DONE\n",wn,b->wn,b->word));
      /* Free suffix array */
      free_sa(b);
      b->status=SAFILE;
      ++(bs->sawrite);
      continue;
    }



    /* Write BWT */
    if ( bs->bwtwrite >=0 && bs->bwtwrite < bs->nbuckets && bs->b[bs->bwtwrite]->status == SAFILE ) {
      b=bs->b[bs->bwtwrite];
      b->status=PROCESS;
      pthread_mutex_unlock(&(bs->lock));
      DEBUG1LINE(fprintf(stderr,"Worker %d: Write BWT for word %d %s\n",wn,b->wn,b->word));
      bwtWriteBucket(b,bs->bwtfile);
      DEBUG1LINE(fprintf(stderr,"Worker %d: Write BWT for word %d %s DONE\n",wn,b->wn,b->word));
      ++(bs->bwtwrite);
      continue;
    }



    /* If you made it this far and sorting is finished, exit worker */
    if (bs->sort>=bs->nbuckets) {
      pthread_mutex_unlock(&(bs->lock_fill));
      break;
    }

    /* Look ahead and start filling buckets */
    if ( bs->sort+bs->nfill > bs->nextfill ) {
      int fs = bs->nextfill;
      int fl = bs->nfill;
      if ( bs->nbuckets - bs->nextfill < 3*bs->nfill/2) {
	fl = bs->nbuckets - bs->nextfill;
	bs->nextfill = bs->nbuckets+bs->nfill;
      }
      else bs->nextfill += fl;

      pthread_mutex_lock(&(bs->lock_fill));
      pthread_mutex_unlock(&(bs->lock));

      DEBUG1LINE(fprintf(stderr,"Worker %d: fill buckets for %d..%d\n",wn,fs,fs+fl-1));
      fillBuckets(bs,fs,fl);
      DEBUG1LINE(fprintf(stderr,"Worker %d: fill buckets for %d..%d DONE\n",wn,fs,fs+fl-1));
      pthread_mutex_unlock(&(bs->lock_fill));
      continue;
    }

    /* Stop and wait for the filling to stop, if bucket is not filled */
    b = bs->b[bs->sort];
    if (b->status==UNINIT) {
      DEBUG1LINE(fprintf(stderr,"Worker %d: WAITING for fill\n",wn));
      pthread_mutex_lock(&(bs->lock_fill));
      pthread_mutex_unlock(&(bs->lock_fill));
    }


    /* Start sorting */
    bs->sort += 1;
    pthread_mutex_unlock(&(bs->lock));
    DEBUG1LINE(fprintf(stderr,"Worker %d: sort for word %d %s\n",wn,b->wn,b->word));
    sortBucket(b);
    DEBUG1LINE(fprintf(stderr,"Worker %d: sort for word %d %s DONE\n",wn,b->wn,b->word));
    /* Calculate BWT */
    DEBUG1LINE(fprintf(stderr,"Worker %d: Calc BWT for word %d %s\n",wn,b->wn,b->word));
    bwtBucket(b);
    DEBUG1LINE(fprintf(stderr,"Worker %d: Calc BWT for word %d %s DONE\n",wn,b->wn,b->word));

  }

  DEBUG1LINE(fprintf(stderr,"Worker %d RETURNS\n",wn));
  return NULL;
}




static int compare_SEQstruct(const void *x1, const void *x2) {
  SEQstruct *ss1, *ss2;
  ss1 = *(SEQstruct **)x1;
  ss2 = *(SEQstruct **)x2;
  return compare_strings( (void*)&(ss1->start), (void*)&(ss2->start) );
}




/*
  When looking up a sequence in the suffix array, the correct sort order of
  the sequences are needed.

  Take a linked list of sequences and sort them
  Record the sort rank in ->sort_order of the sequence

  This function also creates ->seqlengths, ->ids, and ->seqTermOrder

  Assumes that ->sort_order holds the terminator sort order

  Initial SEQstruct points to first sequence and total

 */
static void SortSeqs(SEQstruct *ss, suffixArray *s) {
  int i, nseq = s->nseq;
  SEQstruct **ssarray, *cur;

  ssarray = (SEQstruct **)malloc(nseq*sizeof(SEQstruct*));
  s->ids = (char **)malloc(nseq*sizeof(char*));
  s->seqlengths = (IndexType *)malloc(nseq*sizeof(IndexType));
  s->seqTermOrder = (int *)malloc(nseq*sizeof(int));

  cur = ss->next;
  for (i=0; i<nseq; ++i) { ssarray[i]=cur; cur=cur->next; }

  qsort(ssarray,nseq,sizeof(SEQstruct*),compare_SEQstruct);

  // for (i=0; i<nseq; ++i) fprintf(stderr,"%s\n",ssarray[i]->id);

  for (i=0; i<nseq; ++i) {
    cur = ssarray[i];
    /* Save ID and length in sorted order */
    s->ids[i] = cur->id;
    s->seqlengths[i] = cur->len;
    /* Save the seqTermOrder */
    s->seqTermOrder[i] = cur->sort_order;
    /* Record the sort order */
    cur->sort_order = i;
  }

  free(ssarray);
}




/*
  Compare two suffixes str1 & str2.
  Use 0 termination of sequences (concatenated into one long seq).

  If they are equal, compare backwards!

  backwards comparison is also 0 terminated (so the concatenated string
  should start with a 0).

  Note that the backward comparison is relatively rare, so aggressive
  optimization is not so important.

  Backward comparison happens only when two suffixes are equal distance from
  the 0 termination (otherwise the shortest sorts first).

  Note also that there is zero padding, so the sequence may start with
  a 0, because the suffix starts after the word. This happens only at ends
  of sequenecs, and therefore we immediately start back-sorting.

*/
static int compare_strings_reverse(const void *str1, const void *str2) {
  int i;
  char *s1, *s2;
  i = compare_strings(str1,str2);
  // i = strcmp(*(char **)str1,*(char **)str2);
  if (i==0) {
    /* If word size was known, more could be subtracted */
    s1 = *(char **)str1-1;
    s2 = *(char **)str2-1;
    /* Skip zeros backward (this can only happen at ends of sequences) */
    while (*s1=='\0' && *s2=='\0') { s1--; s2--; }
    while ( i==0 ) {
      i=(int)*s1 -(int)*s2;
      if ( *s1==0 || *s2==0 ) break;
      s1--; s2--;
    }
  }

  return i;
}




static int compare_SEQstruct_reverse(const void *x1, const void *x2) {
  SEQstruct *ss1, *ss2;
  char *s1, *s2;

  ss1 = *(SEQstruct **)x1;
  ss2 = *(SEQstruct **)x2;

  s1=ss1->start+ss1->len;
  s2=ss2->start+ss2->len;

  //  return compare_reverse(ss1->start+ss1->len, ss2->start+ss2->len);
  return compare_strings_reverse( (void*)&s1, (void*)&s2);
}



/*
  The suffixes starting with 0 are reverse sorted.

  In order to be able to recover the sequence, it is therefore necesary to
  know the sort order of the zeros.

  Take a linked list of sequences and sort them reversely.

  Record in sort_order
 */
void revSortSeqs(SEQstruct *ss) {
  int i, nseq = ss->sort_order;
  SEQstruct **ssarray, *cur;

  ssarray = (SEQstruct **)malloc(nseq*sizeof(SEQstruct*));

  cur = ss->next;
  for (i=0; i<nseq; ++i) { ssarray[i]=cur; cur=cur->next; }
  qsort(ssarray,nseq,sizeof(SEQstruct*),compare_SEQstruct_reverse);

  /* Record the sort order */
  for (i=0; i<nseq; ++i) ssarray[i]->sort_order = i;

  free(ssarray);
}


/* Set the sort_order equal to order in linked list */
void readOrder(SEQstruct *ss) {
  int i, nseq = ss->sort_order;
  SEQstruct *cur;

  cur = ss->next;
  for (i=0; i<nseq; ++i) { cur->sort_order = i; cur=cur->next; }
}

/*
  Encode the sort_order in the end of each sequence
  We're assuming that there is enough padding to hold the encoded number
  Numbers are encoded with char alen to 127
*/
void encodeOrder(SEQstruct *ss, int alen) {
  int i, wl, nseq = ss->sort_order;
  SEQstruct *cur;
  char alphabet[128];
  //int bigAlen = 128-alen;
  int bigAlen = 100;

  /* Word len */
  wl = 1+(int)(log10(nseq)/log10(bigAlen));
  /* Alphabet */
  // for (i=0; i<bigAlen; ++i) alphabet[i]=i+alen;
  for (i=0; i<bigAlen; ++i) alphabet[i]=i-bigAlen-1;
  alphabet[i]=0;

  fill_powers(wl, bigAlen);

  cur = ss->next;
  while (cur) {
    number2word(cur->sort_order, bigAlen, cur->start+cur->len, wl, alphabet);
    // fprintf(stderr,"%s\n",cur->id);
    cur=cur->next;
  }
}


/* Write the initial part of the BWT corresponding to the term symbols */
void write_term(SEQstruct *ss, suffixArray *sa, FILE *fp) {
  SEQstruct *cur;
  int i;
  char *bwt = (char *)malloc(sa->nseq*sizeof(char));

  cur = ss->next;
  while (cur) {
    i = sa->seqTermOrder[cur->sort_order];
    bwt[i] = *(cur->start+cur->len-1);
    cur = cur->next;
  }
  fwrite(bwt,1,sa->nseq,fp);
  free(bwt);
}



/*
  If the string equals DNA, RNA or protein, the alphabet is set
  to appropriate string.

 */
char *read_alphabet(char *alphabet, char term) {
  char *a;
  int l;

  if (!term) term='*';

  if (strcmp(alphabet,"DNA")==0) { l=6; a = malloc(l+2); strcpy(a+1,"ACGTN"); }
  else {
    if (strcmp(alphabet,"RNA")==0) { l=6; a = malloc(l+2); strcpy(a+1,"ACGUN"); }
    else {
      if (strcmp(alphabet,"protein")==0) { l=22; a = malloc(l+2); strcpy(a+1,"ACDEFGHIKLMNPQRSTVWYX"); }
      else {
	l = 1+strlen(alphabet);
	a = malloc(l+2);
	strcpy(a+1,alphabet);
      }
    }
  }
  a[0]=term;

  return a;
}




// call with NULL to initialize
void print_time(char *text) {
  static clock_t tic;
  clock_t tac;
  if (text) {
    tac=clock();
    fprintf(stderr,"%s time = %fs\n",text,(double)(tac-tic)/CLOCKS_PER_SEC);
    tic = tac;
  }
  else tic=clock();
}



int main(int argc, char **argv) {
  long approx_seqlen, bwtlen;
  int alen, i, j, k;
  char *alphabet, *filename;
  AlphabetStruct *astruct;
  FILE *fp;
  FILE *bwtfile, *sa_file;
  int wlen, nseq;
  SEQstruct *ss;
  suffixArray *sa_struct;
  int padding=10;          // The zero padding between sequences. Needs to be able to hold seq order
                           // and word length. 10 is more than enough for the first.

  print_time(NULL);

  /* Parsing options and arguments */
  OPT_read_cmdline(opt_struct, argc, argv);
  if (help) { OPT_help(opt_struct); exit(0); }
  OPT_print_vars(stderr, opt_struct, "# ", 0);

  /* Input file and approximate sequence length */
  if (length>=1.e-6) approx_seqlen = (long)( (length+0.1)*1.e6 );
  else approx_seqlen = 0;
  if (infilename) {
    fp = fopen(infilename,"r");
    if (!fp) ERRORs("mkbwt: Can't open file %s\n",infilename, 1);
  }
  else {
    if (!approx_seqlen) ERROR("mkbwt: You need to specify length (-l) in MB, when sequences are read on stdin",2);
    fp = stdin;
  }

  /* Check if output file name is given */
  if (!outfilename) {
    if (!infilename) ERROR("mkbwt: You need to specify name for output file if sequences are read on stdin",2);
    outfilename = infilename;
  }

  /* First set alphabet (allocated in the option parsing code) */
  alphabet = read_alphabet(Alphabet,term[0]);
  astruct = alloc_AlphabetStruct(alphabet, caseSens, revComp);
  alen = astruct->len;

  /* Set word length */
  wlen = 2;
  if (alen<=16) wlen = 3;
  if (alen<=8)  wlen = 4;
  if (padding<wlen) padding=wlen;

  /* Read sequences from fasta file */
  ss = readFasta(fp, approx_seqlen, astruct->trans, astruct->comp, 0, padding);
  if (infilename) fclose(fp);

  print_time("Sequences read");

  fprintf(stderr,"SLEN %ld\nNSEQ %d\nALPH %s",ss->len,ss->sort_order,alphabet);
  if (astruct->comp) fprintf(stderr," (%s)",astruct->comp);
  fprintf(stderr,"\n");

  /* ss->next points to the first seq. Other fields of ss contains other info */
#ifdef DEBUG2
  print_seq(ss->start, ss->start, 0, ss->len, 50, alphabet, alen, stderr);
#endif

  /* Open files */
  int l=strlen(outfilename);
  filename = (char *)malloc((l+10)*sizeof(char));
  strcpy(filename,outfilename);
  strcpy(filename+l,".bwt");
  bwtfile = fopen(filename,"w");
  strcpy(filename+l,".sa");
  sa_file = fopen(filename,"w");
  free(filename);

  // Alloc and init suffix array
  sa_struct = init_suffixArray(ss, checkpoint);

  DEBUG1LINE(fprintf(stderr,"SA initiated\n"));

  bwtlen=sa_struct->len;
  nseq=sa_struct->nseq;

  // Write header for bwtfile
  fwrite(&bwtlen,sizeof(IndexType),1,bwtfile);
  fwrite(&nseq,sizeof(int),1,bwtfile);
  fwrite(&alen,sizeof(int),1,bwtfile);
  fwrite(alphabet,sizeof(char),alen,bwtfile);

  DEBUG1LINE(fprintf(stderr,"BWT header written\n"));

  /* Sorting sequence terminations

     Seq terminations can in principle be sorted anyway you like.
     Two methods are implemented: sorted as read or sorted reverse.
     If they are reverse sorted, the BWT will be more compressible.

     Sorting is implemented by extending the sequence with words that
     sort in the required sort order.
   */
  if (revsort) revSortSeqs(ss);
  else readOrder(ss);
  /* Encode ordering in sequence */
  encodeOrder(ss,alen);

  DEBUG1LINE(fprintf(stderr,"Order encoded\n"));

  /* Alloc bucket stack */
  BucketStack *wbs = initBucketStack(alen,alphabet,ss->len,ss->start,bwtfile, sa_file, sa_struct, wlen);

  DEBUG1LINE(fprintf(stderr,"Bucket stack initiated\n"));

  /* Number of buckets to fill. Ad hoc at the moment... */
  wbs->nfill = wbs->nbuckets/20+nThreads;

  /* Start nThreads-1 to begin with */
  pthread_t *worker = (pthread_t *)malloc(nThreads*sizeof(pthread_t));
  for (i=0; i<nThreads-1; ++i) {
    pthread_create(&(worker[i]), NULL, BucketSorter, wbs);
  }

  /* Do other work in main thread independent of SA sorting */

  /* Do the sorting of seqs */
  SortSeqs(ss, sa_struct);
  DEBUG1LINE(fprintf(stderr,"Sequences sorted\n"));

  /* Write first part of BWT */
  write_term(ss, sa_struct, bwtfile);
  DEBUG1LINE(fprintf(stderr,"BWT for term chars written\n"));

  //sa_struct->seqTermOrder = revSortSeqs(ss, bwtfile);

  /* Write SA header */
  write_suffixArray_header(sa_struct, sa_file);
  free(sa_struct->seqTermOrder); sa_struct->seqTermOrder=NULL;
  free(sa_struct->ids); sa_struct->ids = NULL;
  DEBUG1LINE(fprintf(stderr,"SA header written\n"));

  /* Allow for writing of other buckets */
  wbs->sawrite = wbs->bwtwrite = 0;

  /* Start last thread */
  pthread_create(&(worker[nThreads-1]), NULL, BucketSorter, wbs);

  /* Join workers */
  for (i=0; i<nThreads; ++i) {
    DEBUG1LINE(fprintf(stderr,"Join worker %d\n",i));
    pthread_join(worker[i], NULL);
  }
  free(worker);

  /* Print last */
  while ( wbs->bwtwrite < wbs->nbuckets ) {
    Bucket *b = wbs->b[wbs->bwtwrite];
    DEBUG1LINE(fprintf(stderr,"Finish bucket for word %s\n",b->word));
    if (b->status == BWT ) {
      /* Write SA checkpoints */
      write_suffixArray_checkpoints(b->sa, b->start, b->len, wbs->sa_struct, wbs->safile);
      /* free SA */
      free_sa(b);
    }
    /* Write BWT */
    bwtWriteBucket(b,wbs->bwtfile);
    ++(wbs->bwtwrite);
  }

  fprintf(stderr,"SA NCHECK=%ld\n",sa_struct->ncheck);

  fclose(bwtfile);
  fclose(sa_file);

  print_time("Sorting done, ");

  /* Free a lot of stuf.... */

}


