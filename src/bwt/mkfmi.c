/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>

#include "common.h"
#include "fmi.h"
#include "bwt.h"
#include "mkfmi_vars.h"



int main (int argc, char **argv) {
  int i;
  uchar *bwt;
  FILE *fp=NULL;
  FMI *f;
  BWT *s;
  char filename[1000];

  /* Parsing options and arguments */
  OPT_read_cmdline(opt_struct, argc, argv);
  if (help) { OPT_help(opt_struct); exit(0); }
  OPT_print_vars(stderr, opt_struct, "# ", 0);

  if (!outfilenm) {
    fprintf(stderr,"You have to specify name for output files (first argumetn)\n");
    exit(5);
  }

  if (!bwtfile) fp=stdin;
  else {
    fp = fopen(bwtfile,"r");
    if (!fp) { fprintf(stderr,"Could not open BWT file, %s\n",bwtfile); exit(2); }
  }
  s = read_BWT_header(fp);
  bwt = read_BWT(s->len,s->alen,s->alphabet, fp);
  fclose(fp);

  fprintf(stderr,"BWT of length %ld has been read with %d sequencs, alphabet=%s\n",
	  s->len,s->nseq, s->alphabet); 

  fprintf(stderr,"Constructing FM index ... ");
  f = makeIndex(bwt, s->len, s->alen);
  fprintf(stderr,"DONE\n");

  fprintf(stderr,"Constructing suffix array checkpoints, START\n");
  init_suffixArray(s, f, checkpoint);
  // Filling intervals of suffix checkpoints in parallel
  if (nThreads<2) nThreads=2;
  pthread_t *worker = (pthread_t *)malloc(nThreads*sizeof(pthread_t));
  for (i=0; i<nThreads-1; ++i) pthread_create(&(worker[i]), NULL, fill_suffixArray_parallel, s);

  fprintf(stderr,"Writing FM index to file ... ");
  sprintf(filename,"%s.fmi",outfilenm);
  fp = fopen(filename,"w");
  write_fmi(f, fp);
  fclose(fp);
  fprintf(stderr,"DONE\n");

  // Starting last worker
  if (nThreads>1) pthread_create(&(worker[nThreads-1]), NULL, fill_suffixArray_parallel, s);

  for (i=0; i<nThreads; ++i) pthread_join(worker[i], NULL);
  free(worker);

  fprintf(stderr,"Constructing suffix array checkpoints, DONE\n");

  fprintf(stderr,"Writing suffix array checkpoints and sequence info to file ... ");
  sprintf(filename,"%s.sa",outfilenm);
  fp = fopen(filename,"w");
  write_BWT_bin(s, fp);
  fclose(fp);
  fprintf(stderr,"DONE\n");

}
