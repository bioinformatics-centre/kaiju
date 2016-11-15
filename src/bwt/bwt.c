/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>

#include "common.h"
#include "bwt.h"
#include "fmi.h"
#include "suffixArray.h"



/***********************************************************************
 *
 * The BWT struct collects everything and points to both the suffix array
 * and the FMI
 *
 ***********************************************************************/





/***********************************************
 *
 * I/O stuff
 *
 ***********************************************/





/*
	 Write BWT header
	 */
void write_BWT_header(BWT *b, FILE *bwtfile) {
	fwrite(&(b->len),sizeof(IndexType),1,bwtfile);
	fwrite(&(b->nseq),sizeof(int),1,bwtfile);
	fwrite(&(b->alen),sizeof(int),1,bwtfile);
	fwrite(b->alphabet,sizeof(char),b->alen,bwtfile);
}


/*
	 Read BWT from file (made by mkbwt) and resturn in BWT struct
	 */
static BWT *read_BWT_header(FILE *bwtfile) {
	BWT *b=(BWT *)malloc(sizeof(BWT));

	fread(&(b->len),sizeof(IndexType),1,bwtfile);
	fread(&(b->nseq),sizeof(int),1,bwtfile);
	fread(&(b->alen),sizeof(int),1,bwtfile);
	b->alphabet=(char *)calloc(sizeof(char),b->alen+1);
	fread(b->alphabet,sizeof(char),b->alen,bwtfile);

	return b;
}


/*
	 Read BWT from file (made by mkbwt) and return in BWT struct
	 */
BWT *read_BWT(FILE *bwtfile) {
	BWT *b=read_BWT_header(bwtfile);
	b->bwt = (uchar *)malloc(b->len*sizeof(uchar));
	fread(b->bwt,sizeof(uchar),b->len,bwtfile);
	return b;
}


/*
	 Read indexes from one file (made by mkfmi) and resturn in BWT struct
	 */
BWT *readIndexes(FILE *fp) {
	BWT *b=read_BWT_header(fp);

	b->bwt=NULL;

	b->s = read_suffixArray_header(fp);
	read_suffixArray_body(b->s, fp);
	b->f = read_fmi(fp);

	return b;
}





/***********************************************
 *
 * Query SA && FMI
 *
 ***********************************************/



/* Find suffix for suffix number i
	 Return sequence number in *iseq and position in *pos
	 */
void get_suffix(FMI *fmi, suffixArray *s, IndexType i, int *iseq, IndexType *pos) {
	IndexType k=0;
	uchar c=1;

	while ( c && (i & s->check) ) {
		i = FMindexCurrent(fmi,&c,i);
		++k;
	}

	if (c) {
		suffixArray_decode_number(iseq, pos,
				(i>>s->chpt_exp)-((s->nseq-1)>>s->chpt_exp)-1, s);
		*pos += k;
	}
	else { *iseq = i; *pos=k-1; }

}


/*
	 Reconstruct sequence number snum (according to the original order of the sequence file)
	 */
uchar *retrieve_seq(int snum, BWT *b) {
	IndexType l, k;
	uchar *seq, *cur;

	l = b->s->seqlengths[snum];
	k = (IndexType)(b->s->seqTermOrder[snum]);
	seq = (uchar *)malloc((l+2)*sizeof(uchar));

	cur = seq+l+1;
	*cur-- = 0; *cur-- = 0;
	while (cur>=seq) {
		k = FMindexCurrent(b->f, cur, k);
		--cur;
	}
	return seq;
}


/* Find whole (initial) suffix interval for letter ct */
IndexType InitialSI(FMI *f, uchar ct, IndexType *si) {
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
	new->count=new->len;               // ->len is the si length = number of matches
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
	 Returns null if there are no matchesBWT *b
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
				// fprintf(stderr,"MATCH %d %d-%d %ld\n",(int)(si[1]-si[0]), i,j, l);
				cur = alloc_SI(si, i, l);
				first = insert_SI_sorted(first, cur);
				// If max_matches is set, check to see if max is reached and reset L
				if (max_matches>0) {
					k = free_until_max_SI(first, max_matches);
					if (k>L) L=k;
					// The latest si may be freed if too short - then set it to NULL
					if (l<k) cur=NULL;
				}
			}
		}
		// If the last match reached beginning of sequence, no need to continue                                                                                                                                    
		if (i<=1) break;   
	}

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
	//	if ( !cur || i < cur->qi ) { // cur is never != NULL here
			//cur = alloc_SI(si, i, l);
			first = alloc_SI(si, i, l);
			//first = insert_SI_sorted(first, cur);
			// If max_matches is set, check to see if max is reached and reset L
			//if (max_matches>0) { this function is always called with max_matches=1, so save the if
			//	k = free_until_max_SI(first, max_matches);
			//	if (k>L) L=k;
				// The latest si may be freed if too short - then set it =NULL
			//	if (l<k) cur=NULL;
			//}
	//	}
	}
	// If the last match reached beginning of sequence, no need to continue
	//if (i<=1) break;
	//}

	return first;
}



/* Find maximal matches (longer than L) of str of length len in a linked list.
	 Returns all matches of maximal length.
	 Returns null if there are no matches

	 if jump is positive, it jumps by L-jump after having a match of length L
	 Note that L is dynamic (max length found)
	 */
SI *greedyExact(FMI *f, char *str, int len, int L, int jump) {
	SI *tmp, *first=NULL, *cur=NULL;
	IndexType si[2], l;
	int i, j, delta=1;

	if (jump>=0) delta=L-jump;

	// Go through the sequence from the back
	for (j=len-1; j>=L-1; j-=delta) {
		i=j;
		InitialSI(f, str[i], si);
		// Extend backward
		while ( i-- > 0 ) {
			if ( UpdateSI(f, str[i], si, NULL) == 0) break;
		}
		i+=1; // Start of match
		l = j-i+1;

		if (l>=L) {
			if (l>L) {
				recursive_free_SI(first);  // Free the shorter ones
				first = NULL;
				L=l;
				if (jump>=0) delta=L-jump;
			}
			cur = first;
			first = alloc_SI(si, i, l);
			first->samelen=cur;
		}
		if (i<=1) break;
	}

	return first;
}



