/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */



typedef struct __SEQstruct__ {
  char *id;
  char *descr;     // Description (stuff following id). Allocated with id.
  long len;
  long fpos;       // Sequence position for sequence
  long rpos;       // Sequence position for start of rev. complement
  char *start;
  long id_filepos;
  long seq_filepos;
  int sort_order;
  struct __SEQstruct__ *next;
} SEQstruct;


static SEQstruct *alloc_SEQstruct() {
  SEQstruct *ss=(SEQstruct *)malloc(sizeof(SEQstruct));
  ss->len = 0;
  ss->fpos = 0;
  ss->rpos = 0;
  ss->id = NULL;
  ss->descr = NULL;
  ss->start = NULL;
  ss->id_filepos = 0;
  ss->seq_filepos = 0;
  ss->sort_order = 0;
  ss->next=NULL;
  return ss;
}

static void free_SEQstruct(SEQstruct *ss) {
  if (ss) {
    if (ss->id) free(ss->id);
    free(ss);
  }
}

/* Assumes that base->start points to the whole sequence */
static void recursive_free_SEQstruct(SEQstruct *base) {
  SEQstruct *ss, *next;
  if (base->start) free(base->start);
  ss=base;
  next=ss->next;
  while (ss) {
    free_SEQstruct(ss);
    ss=next;
    if (next) next=ss->next;
  }
}

