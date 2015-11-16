/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "common.h"

/* Writes the number of bytes (long) followed by the string.
   If n<0, strlen((char*)c)+1 is used - so the terminating zero is included
*/
void fwrite_n_bytes(void *c, long n, FILE *fp) {
  if (n<0) n = strlen((char*)c)+1;
  fwrite((void *)&n,sizeof(long),1,fp); 
  fwrite(c,1,n,fp); 
}


void *fread_n_bytes(long *n, FILE *fp) {
  void *c = NULL;
  fread((void *)n,sizeof(long),1,fp);
  c = (void *)malloc(*n);
  fread(c,1,*n,fp);
  return c;
}

/* Allocates extra bytes */
void *fread_n_bytes_plus(long *n, FILE *fp, int extra) {
  void *c = NULL;
  fread((void *)n,sizeof(long),1,fp);
  c = (void *)malloc(*n+extra);
  fread(c,1,*n,fp);
  return c;
}

