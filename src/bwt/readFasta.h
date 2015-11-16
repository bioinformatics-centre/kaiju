/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef READFASTA_h
#define READFASTA_h

#include "common.h"

/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
SEQstruct *readFasta(FILE *fp, long length, char *transtab, char *complement, char term, int padding);
char *translation_table(char *alphabet, char *translation, char dummy, int casesens);
void writeSequenceFile(SEQstruct *base, char *alph, FILE *fp);
SEQstruct *readIndex(FILE *fp);
uchar *readSeqOnly(long *len, FILE *fp, int padding, char term);
SEQstruct *readSeq(FILE *fp, int padding, char term);
SEQstruct **make_seq_hash(SEQstruct *base, int Estep);
SEQstruct *lookupSeq(IndexType n, SEQstruct **hash, int Estep);
/* FUNCTION PROTOTYPES END */

#endif
