/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef COMMON_h
#define COMMON_h

typedef unsigned char uchar;
typedef unsigned short int ushort;
typedef unsigned int uint;
typedef long int IndexType;


/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
void fwrite_n_bytes(void *c, long n, FILE *fp);
void *fread_n_bytes(long *n, FILE *fp);
void *fread_n_bytes_plus(long *n, FILE *fp, int extra);
/* FUNCTION PROTOTYPES END */

#endif
