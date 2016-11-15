/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "common.h"
#include "fmi.h"
#include "bwt.h"
#include "suffixArray.h"
#include "mkfmi_vars.h"

void error(char *format, char *arg) {
  fprintf(stderr,"ERROR: ");
  fprintf(stderr,format,arg);
  exit(1);
}


int main (int argc, char **argv) {
  int l;
  FILE *fp=NULL;
  BWT *b;
  char *filename;

  /* Parsing options and arguments */
  OPT_read_cmdline(opt_struct, argc, argv);
  if (help) { OPT_help(opt_struct); exit(0); }
  OPT_print_vars(stderr, opt_struct, "# ", 0);

  if (!filenm) {
    fprintf(stderr,"You have to specify name for index files (first argument)\n");
    exit(5);
  }

  l=strlen(filenm);
  filename = (char *)malloc((l+10)*sizeof(char));
  strcpy(filename,filenm);

  /* Read BWT */
  strcpy(filename+l,".bwt");
  fp = fopen(filename,"r");
  if (!fp) error("File %s containing BWT could not be opened for reading\n",filename);
  fprintf(stderr,"Reading BWT from file %s ... ",filename);
  b = read_BWT(fp);
  fclose(fp);
  fprintf(stderr,"DONE\n");
  fprintf(stderr,"BWT of length %ld has been read with %d sequencs, alphabet=%s\n",
	  b->len, b->nseq, b->alphabet); 

  /* Read SA */
  strcpy(filename+l,".sa");
  fp = fopen(filename,"r");
  if (!fp) error("File %s containing SA could not be opened for reading\n",filename);
  fprintf(stderr,"Reading suffix array from file %s ... ",filename);
  b->s = read_suffixArray_header(fp);
  /* If the whole SA is saved, don't read it! */
  if (b->s->chpt_exp > 0) read_suffixArray_body(b->s,fp);
  fclose(fp);
  fprintf(stderr,"DONE\n");

  /* Concatenate stuff in fmi file */
  strcpy(filename+l,".fmi");
  fp = fopen(filename,"w");
  if (!fp) error("File %s for FMI could not be opened for reading\n",filename);
  fprintf(stderr,"Writing BWT header and SA to file  %s ... ",filename);
  write_BWT_header(b, fp);
  write_suffixArray(b->s,fp);
  fprintf(stderr,"DONE\n");

  fprintf(stderr,"Constructing FM index\n");
  b->f = makeIndex(b->bwt, b->len, b->alen);
  fprintf(stderr,"\nDONE\n");

  fprintf(stderr,"Writing FM index to file ... ");
  write_fmi(b->f,fp);
  fclose(fp);
  fprintf(stderr,"DONE\n");

  if (removecmd) {
    int cl = strlen(removecmd);
    char *command = malloc(cl+2*l+20);
    sprintf(command,"%s %s.sa %s.bwt",removecmd,filenm,filenm);
    fprintf(stderr,"Removing files with this command: %s\n",command);
    system(command);
    free(command);
  }
  else {
    strcpy(filename+l,".bwt");
    fprintf(stderr,"\n  !!  You can now delete files %s and ",filename);
    strcpy(filename+l,".sa");
    fprintf(stderr,"%s  !!\n\n",filename);
    free(filename);
  }

}
