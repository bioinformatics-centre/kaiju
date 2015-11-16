CC      = gcc # -static -lm
#CFLAGS  = -g
CFLAGS  = -w -O3 -g
LDLIBS = -lpthread


ALL: mkbwt mkfmi

mkbwt: mkbwt.o readFasta.o common.o

mkfmi: mkfmi.o bwt.o simplefmi.o common.o

mkbwt.o: mkbwt_vars.h mkbwt.c

mkfmi.o: mkfmi_vars.h mkfmi.c fmi.h

readFasta.o: readFasta.c readFasta.h sequence.h common.h

common.o: common.c common.h

simplefmi.o: simplefmi.c simplefmi.h common.h fmicommon.h

bwt.o: bwt.c bwt.h fmi.h common.h

suffixArray.o: suffixArray.c suffixArray.h bwt.h fmi.h common.h


clean:
	rm -f *.o *.h.old *~ mkbwt mkfmi

.PHONY: clean 
