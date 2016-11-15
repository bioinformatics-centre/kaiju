/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
/*
  Copied from 
  http://www.drdobbs.com/database/sorting-strings-with-three-way-radix-qui/184410724
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#ifndef min 
#define min(a, b) ((a)<=(b) ? (a) : (b)) 
#endif


static inline void swap(char *a[], int i, int j) 
{     char *t = a[i];
      a[i] = a[j];
      a[j] = t; 
}
static inline void vecswap(char *a[], int i, int j, int n) 
{     while (n-- > 0)
         swap(a, i++, j++); 
}


#define ch(i) a[i][depth] 


/* Faster version */

/* Std strcmp use unsigned comparison */
static inline int my_strcmp(const char *s1, const char *s2) {
  while (*s1 && *s2 && *s1==*s2) { ++s1; ++s2; }
  return (int)( *s1 - *s2 );
}


static inline int med3func(char *a[], int ia, int ib, int ic, int depth) 
{   int va, vb, vc;
    if ((va=ch(ia)) == (vb=ch(ib)))
         return ia;
    if ((vc=ch(ic)) == va || vc == vb)
         return ic;
    return va < vb ?
          (vb < vc ? ib : (va < vc ? ic : ia ) )
        : (vb > vc ? ib : (va < vc ? ia : ic ) ); 
} 


void inssort(char *a[], int n, int depth) 
{   int i, j;
    for (i = 1; i < n; i++)
      for (j = i; j > 0; j--) {
         if (my_strcmp(a[j-1]+depth, a[j]+depth) <= 0)
             break;
         swap(a, j, j-1);
      } 
}  


void ssort2(char *a[], int n, int depth) 
{    int le, lt, gt, ge, r, v;
     int pl, pm, pn, d;

     if (n <= 10) {
        inssort(a, n, depth);
        return;
     }

     pl = 0;
     pm = n/2;
     pn = n-1;
     if (n > 50) {
        d = n/8;
        pl = med3func(a, pl, pl+d, pl+2*d,depth);
        pm = med3func(a, pm-d, pm, pm+d,depth);
        pn = med3func(a, pn-2*d, pn-d, pn,depth);
     }
     pm = med3func(a, pl, pm, pn,depth);
     swap(a, 0, pm);
     v = ch(0);
     for (le = 1; le < n && ch(le) == v; le++)
       ;  
     if (le == n) {
         if (v != 0) ssort2(a, n, depth+1);
         return;
     }
     lt = le;
     gt = ge = n-1;
     for (;;) {
         for ( ; lt <= gt && ch(lt) <= v; lt++)
             if (ch(lt) == v) swap(a, le++, lt);
         for ( ; lt <= gt && ch(gt) >= v; gt--) {
	   if (ch(gt) == v) swap(a, gt, ge--);
	 }
         if (lt > gt)
             break;
         swap(a, lt++, gt--);
     }
     r = min(le, lt-le);
     vecswap(a, 0, lt-r, r);
     r = min(ge-gt, n-ge-1);
     vecswap(a, lt, n-r, r);
     ssort2(a, lt-le, depth);
     if (v != 0)
       ssort2(a + lt-le, le + n-ge-1, depth+1);
     ssort2(a + n-(ge-gt), ge-gt, depth); 
}


//void ssort2main(char *a[], int n) 
void multikeyqsort(char *a[], int n)
{ ssort2(a, n, 0); }
