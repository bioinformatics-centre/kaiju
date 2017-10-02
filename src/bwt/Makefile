CC = gcc
#CFLAGS  = -g
CFLAGS = -O3 -g -Wno-unused-result
LDLIBS = -lpthread -lm

ifeq ($(uname -s), "Darwin")
LD_LIBS_STATIC = -Wl,-all_load -lpthread -Wl,-noall_load -lm
else
LD_LIBS_STATIC = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lm
endif

all: mkbwt mkfmi Makefile

mkbwt: mkbwt.o readFasta.o suffixArray.o multikeyqsort.o sequence.o

mkfmi: mkfmi.o bwt.o suffixArray.o compactfmi.o

mkbwt.o: mkbwt_vars.h mkbwt.c common.h multikeyqsort.h sequence.h

mkfmi.o: mkfmi_vars.h mkfmi.c fmi.h common.h

sequence.o: sequence.h common.h

readFasta.o: readFasta.c readFasta.h sequence.h common.h

compactfmi.o: compactfmi.c compactfmi.h common.h fmicommon.h

suffixArray.o: suffixArray.c suffixArray.h common.h sequence.h

bwt.o: bwt.c bwt.h fmi.h common.h

multikeyqsort.o: multikeyqsort.c multikeyqsort.h

clean:
	rm -f mkfmi mkbwt

static: LDFLAGS = -static
static: LDLIBS = $(LD_LIBS_STATIC)
static: all

debug: all

.PHONY: clean static debug
