/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
#ifndef READFASTA_h
#define READFASTA_h

#include "common.h"

/* FUNCTION PROTOTYPES BEGIN  ( by funcprototypes.pl ) */
SEQstruct *revcompSEQstruct(SEQstruct *ss, char *s, char *translate);
SEQstruct *readFasta(FILE *fp, long length, char *transtab, char *complement, char term, int padding);
char *translation_table(char *alphabet, char *translation, char dummy, int casesens);
/* FUNCTION PROTOTYPES END */

#endif
