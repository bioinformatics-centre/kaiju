/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
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

// #define DEBUG1
// #define DEBUG2

#ifdef DEBUG1
#define DEBUG1LINE(x) x
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
#include <unistd.h>
#include <pthread.h>
#include <time.h>

#include "sequence.h"
#include "readFasta.h"
#include "mkbwt_vars.h"


/* Global var for setting of number for worker */
int workerNum=0;


/* Status flags */
#define UNINIT 0
#define FILLED 1
#define SORTED 2
#define BWT 3
#define DONE 4
#define PROCESS 5


typedef struct {
  int wlen;     /* prefix length */
  int wn;       /* prefix (word) number */
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
  int write;      /* Next bucket to be written */
  int nthreads;   /* Number of threads */
  int nfill;      /* Number of buckets to fill in each filling round */
  int nextfill;   /* Number of buckets to fill next */

  pthread_mutex_t lock;  /* Lock to be set by the thread while grabbing the next job */
  pthread_mutex_t lock_fill;  /* Lock to be set by the thread while grabbing the next job */

  FILE *bwtfile;  /* Pointer to an open file for writing of BWT */

} BucketStack;




/* Functions for calculating word numbers
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



static inline int word_number(char *s, int n) {
  int w = 0;
  while (n-- > 0) w += _letter_numbers_[n][ (int)s[n] ];
  return w;
}


/* Word corresponding to number
   Mostly for debugging
*/
static int *powers=NULL;
static void fill_powers(int l, int al) {
  int i;
  powers = (int *)malloc(l*sizeof(int));
  powers[0]=1;
  for (i=1 ; i<l ; ++i) powers[i]=al*powers[i-1];
}

static char *number2word(int n, int alen, char *w, int wlen, char *alphabet) {
  int i=0, nw, l=wlen;
  if (!powers) fill_powers(wlen, alen);
  while (l>0) { nw = n/powers[--l]; w[i++]=nw; n -= nw*powers[l]; }
  w[i]=0;
  if (alphabet) for (i=0; i<wlen; ++i) w[i]=alphabet[w[i]];
  return w;
}


/* Alloc Bucket, but don't alloc suffix array */
Bucket *allocBucket(long slen, char *seq, int wlen, int wn, int alen,
		    char *alphabet, long bucket_size) {
  int i;
  Bucket *bucket = (Bucket *)malloc(sizeof(Bucket));

  bucket->wlen=wlen;
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
    if (*seq) b[word_number(seq,n)] += 1;
    ++seq;
  }

  return b;
}



BucketStack *initBucketStack(int alen, char *alphabet, long slen, char *seq, FILE *bwtfile, int wlen) {
  int i,j;
  BucketStack *bs=(BucketStack *)malloc(sizeof(BucketStack));

  bs->alen = alen;
  bs->alphabet = alphabet;
  bs->slen = slen;
  bs->seq = seq;
  bs->wlen = wlen;
  bs->bwtfile = bwtfile;

  /* Calculate bucket sizes */
  bs->bucket_size = bucketSizes(seq, slen, alen, wlen, &(bs->nbuckets));

  /* Allocate buckets */
  bs->b = (Bucket **)calloc(bs->nbuckets,sizeof(Bucket *));

  for (i=0; i<bs->nbuckets; ++i) bs->b[i] = allocBucket(slen, seq, wlen, i,
	      alen, alphabet, bs->bucket_size[i]);

  bs->sort = 0;
  bs->write = -1;
  bs->nextfill = 0;

  pthread_mutex_init(&(bs->lock), NULL);
  pthread_mutex_init(&(bs->lock_fill), NULL);

  return bs;
}



#ifdef DEBUG2
/* This is for debugging
   Prints sequence from pointer current to end
   if current==NULL, icurrent is used instead
*/
void print_seq_toend(char *seq, char *current, long icurrent, long slen, char *alphabet, FILE *fp) {
  long i;
  int max = 30;
  if (current) icurrent = current-seq;
  /* Print bwt */
  i = icurrent-1; if (i<0) i=slen-1;
  fprintf(fp,"%c ",alphabet[seq[i]]);
  if (slen>icurrent+max) max=icurrent+max;
  else max=slen;
  for (i=icurrent; i<max; ++i) fprintf(fp,"%c",alphabet[seq[i]]);
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
    if (*seq) {
      i=word_number(seq,bs->wlen);
      if (i>=istart && i<iend) {
	*(bucket[i]->cur) = seq+bs->wlen;
	bucket[i]->cur += 1;
      }
    }
    ++seq;
  }

  /* Flag that these buckets are now ready for sorting */
  for (i=istart; i<iend; ++i) { bucket[i]->status = FILLED; }

}



/*
  Compare two suffixes str1 & str2
  Use 0 termination
*/
static int compare_strings(const void *str1, const void *str2) {
  return strcmp(*(char **)str1,*(char **)str2);
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
  i = strcmp(*(char **)str1,*(char **)str2);
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


/*
  long repeats of the same letter (such as Ns in genomes) are extremely
  slow to sort. This function takes care of such repeats
 */
void repeatSuffixSort(char **s, int l, char a, int wlen) {
  int clow, chigh;
  char *tmp, **test, **low, **high, *limit;

  /* Identify all suffixes starting with letter different from a
   and move to appropriate location in array */
  test = low = s;
  high = s+l-1;
  limit = *s;
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
  limit -= wlen;

  /* Now sort the high and low intervals */
  qsort(s,clow,sizeof(char*),compare_strings_reverse);
  qsort(s+l-chigh,chigh,sizeof(char*),compare_strings_reverse);

  /* Now "backpropagate" low results */
  int count = clow;
  low = s+clow;
  test = s;
  while (count) {
    tmp = *test-wlen-1;
    if (tmp>=limit && *tmp==a) { *low = *test-1; ++low; }
    else --count;
    ++test;
  }

  /* Now "backpropagate" high results */
  count = chigh;
  high = s+l-chigh-1;
  test = s+l-1;
  while (count) {
    tmp = *test-wlen-1;
    if (tmp>=limit && *tmp==a) { *high = *test-1; --high; }
    else --count;
    --test;
  }

}



/* Check for homolymers (repeats)
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



/* Sorts a single bucket and return pointers to right location (-wlen) */
void sortBucket(Bucket *b) {
  int i, h;

  if (b->len==0) {b->status=SORTED; return; }

  h = checkHomoPol(b->sa[0],b->wlen);
  if (h>0) {
    // fprintf(stderr,"sortBucket: sorting homopolymer %d\n",h);
    repeatSuffixSort(b->sa, b->len, h, b->wlen);
  }
  else qsort(b->sa,b->len,sizeof(char*),compare_strings_reverse);

  for (i=0; i<b->len; ++i) b->sa[i]-=b->wlen;
  b->status=SORTED;
}



/* Get the BWT for bucket */
void bwtBucket(Bucket *b) {
  int k;

  if (b->len) {
    b->bwt = malloc(b->len+1);
    for (k=0; k<b->len; ++k) b->bwt[k]=b->alphabet[*(b->sa[k]-1)];
    b->bwt[k] = '\0';
#ifndef DEBUG2
    free(b->sa);
    b->sa=NULL;
#endif
  }
  b->status = BWT;
}



void bwtWriteBucket(Bucket *b, FILE *bwtfile) {
  int k;
#ifdef DEBUG2
  if (b->len) {
    int i;
    for (i=0; i<b->len; ++i) {
      fprintf(bwtfile,"%10d %4d ", (int)(b->sa[i]-b->seq), b->wn );
      print_seq_toend(b->seq, b->sa[i], 0, b->slen, b->alphabet, bwtfile);
    }
  }
#else
  if (b->len) {
    fwrite(b->bwt,sizeof(char),b->len, bwtfile);
    free(b->bwt);
    b->bwt=NULL;
  }
#endif
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

  while( bs->sort < bs->nbuckets ) {

    pthread_mutex_lock(&(bs->lock));

    // While waiting for the unlocking, bs->sort might have been incremented
    if (bs->sort>=bs->nbuckets) {
      pthread_mutex_unlock(&(bs->lock));
      break;
    }

    if ( bs->write >=0 && bs->write < bs->nbuckets && bs->b[bs->write]->status == BWT ) {
      b=bs->b[bs->write];
      b->status=PROCESS;
      pthread_mutex_unlock(&(bs->lock));
      DEBUG1LINE(fprintf(stderr,"Worker %d: writes BWT to file for word %d %s\n",wn,b->wn,b->word));
      bwtWriteBucket(b,bs->bwtfile);
      DEBUG1LINE(fprintf(stderr,"Worker %d: writes BWT to file for word %d %s DONE\n",wn,b->wn,b->word));
      ++(bs->write);
      continue;
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
    if (b->status==0) {
      DEBUG1LINE(fprintf(stderr,"Worker %d: WAITING for fill\n",wn));
      pthread_mutex_lock(&(bs->lock_fill));
      pthread_mutex_unlock(&(bs->lock_fill));
    }

    /* Start sorting */
    bs->sort += 1;
    pthread_mutex_unlock(&(bs->lock));
    DEBUG1LINE(fprintf(stderr,"Worker %d: sort + bwt for word %d %s\n",wn,b->wn,b->word));
    sortBucket(b);

    /* Get the bwt immediately for the newly sorted */
    bwtBucket(b);
    DEBUG1LINE(fprintf(stderr,"Worker %d: sort + bwt for word %d %s DONE\n",wn,b->wn,b->word));
  }

  DEBUG1LINE(fprintf(stderr,"Worker %d RETURNS\n",wn));
  return NULL;
}



// NOT USED. DELETE.
static int compare_reverse(const char *s1, const char *s2) {
  int i=0;
  /* Skip zeros backward (this can only happen at ends of sequences) */
  while (*s1==0 && *s2==0) { s1--; s2--; }
  while ( *s1 && *s2 && i==0 ) { i=(int)*s1-- -(int)*s2--; }

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



void write_ids_bin(SEQstruct *ss, FILE *fp) {
  int l=0;
  SEQstruct *cur;
  // Total length of IDs
  cur = ss->next;
  while (cur) {
    l += strlen(cur->id)+1;
    cur = cur->next;
  }

  cur = ss->next;
  while (cur) {
    fprintf(fp,"%s %ld %d\n",cur->id,cur->len, cur->sort_order);
    cur = cur->next;
  }

}


/*
  The suffixes starting with 0 (terminator) are not sorted in the suffix array
  and are therefore not included in the BWT.

  In order to be able to recover the sequence, it is therefore necesary to
  know the sort order of the zeros.

  Take a linked list of sequences and sort them reversely.
  Record the sort rank in ->sort_order
  Initial SEQstruct points to first sequence and total number of seqs in
  ->sort_order
 */
void revSortSeqs(SEQstruct *ss, char *alphabet, FILE *fp) {
  int i, nseq = ss->sort_order;
  SEQstruct **ssarray, *cur;

  ssarray = (SEQstruct **)malloc(nseq*sizeof(SEQstruct*));

  cur = ss->next;
  for (i=0; i<nseq; ++i) { ssarray[i]=cur; cur=cur->next; }

  qsort(ssarray,nseq,sizeof(SEQstruct*),compare_SEQstruct_reverse);

  /* Record the sort order */
  for (i=0; i<nseq; ++i) ssarray[i]->sort_order = i;

  /* Write result of this sorting */
  cur = ss->next;
  while (cur) {
    fprintf(fp,"%s %ld %d\n",cur->id,cur->len, cur->sort_order);
    cur = cur->next;
  }

  /* Write the BWT for the termination char */
  fprintf(fp,"#BWT\n");
  for (i=0; i<nseq; ++i) {
    cur = ssarray[i];
    fprintf(fp,"%c",alphabet[cur->start[cur->len-1]]);
  }

  free(ssarray);
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
  char *alphabet;
  FILE *fp;
  FILE *bwtfile;
  int wlen;
  SEQstruct *ss, *cur;

  print_time(NULL);

  /* Parsing options and arguments */
  OPT_read_cmdline(opt_struct, argc, argv);
  if (help) { OPT_help(opt_struct); exit(0); }
  OPT_print_vars(stderr, opt_struct, "# ", 0);

  /* First set alphabet (allocated in the option parsing code) */
  alphabet = read_alphabet(Alphabet,term[0]);
  alen = strlen(alphabet);

  /* Set word length */
  wlen = 2;
  if (alen<=16) wlen = 3;
  if (alen<=8)  wlen = 4;

  /* Input file and approximate sequence length */
  if (length>=1.e-6) approx_seqlen = (long)( (length+0.1)*1.e6 );
  else approx_seqlen = 0;
  if (infilename) {
    fp = fopen(infilename,"r");
    if (!fp) { fprintf(stderr,"mkbwt: Can't open file %s\n",infilename); exit(1); }
  }
  else {
    if (!approx_seqlen) {
      fprintf(stderr,"mkbwt: You need to specify length (-l) in MB, when sequences are read on stdin\n");
      exit(1);
    }
    fp = stdin;
  }

  /* Read sequences from fasta file */
  char *tt = translation_table(alphabet, NULL, alen-1, caseSens);
  ss = readFasta(fp, approx_seqlen, tt, NULL, 0, wlen);
  free(tt);
  if (infilename) fclose(fp);

  print_time("Sequences read");

  // fprintf(stderr,"SLEN %ld\nALPH %s\n",ss->len,alen,alphabet);

#ifdef DEBUG2
  print_seq_toend(ss->start, ss->start, 0, ss->len, alphabet, stderr);
#endif

  /* Open file for BWT */
  if (outfilename) bwtfile = fopen(outfilename,"w");
  else bwtfile = stdout;

  /* Calculate the real length of the BWT */
  cur = ss->next;
  bwtlen=0;
  while (cur) {
    bwtlen += cur->len+1;
    cur = cur->next;
  }
  fprintf(stderr,"SLEN %ld\nALPH %s\n",ss->len,alphabet);


  fprintf(bwtfile,"#LEN %ld\n#ALEN %d\n#ALPH %s\n",bwtlen, alen, alphabet);
  fprintf(bwtfile,"#NSEQ %d\n",ss->sort_order);

  /* Alloc all bucket stack */
  BucketStack *wbs = initBucketStack(alen,alphabet,ss->len,ss->start,bwtfile, wlen);

  /* Number of buckets to fill. Ad hoc at the moment... */
  wbs->nfill = wbs->nbuckets/20+nThreads;

  /* Start nThreads-1 to begin with */
  pthread_t *worker = (pthread_t *)malloc(nThreads*sizeof(pthread_t));
  for (i=0; i<nThreads-1; ++i) {
    // fprintf(stderr,"start worker\n");
    pthread_create(&(worker[i]), NULL, BucketSorter, wbs);
  }

  /* Do other work in main thread independent of SA sorting */

  /* Do the sorting and writing of BWT for terminating symbols */
  revSortSeqs(ss, alphabet,bwtfile);

  /* Allow for writing of other buckets */
  wbs->write = 0;

  /* Start last thread */
  pthread_create(&(worker[nThreads-1]), NULL, BucketSorter, wbs);

  /* Join workers */
  for (i=0; i<nThreads; ++i) pthread_join(worker[i], NULL);
  free(worker);

  /* Print last */
  while ( wbs->write < wbs->nbuckets ) {
    Bucket *b = wbs->b[wbs->write];
    bwtWriteBucket(b,wbs->bwtfile);
    ++(wbs->write);
  }

  fclose(bwtfile);

  print_time("Sorting done");

  /* Free a lot of stuf.... */

}


