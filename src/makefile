CC = gcc
CXX = g++
CFLAGS = -O3 -DNDEBUG
CXXFLAGS = -O3 -pthread -std=c++11 -DNDEBUG
LDLIBS = -lpthread -lz
INCLUDES	= -I./include -I./include/ncbi-blast+

BLASTOBJS = include/ncbi-blast+/algo/blast/core/pattern.o \
			 include/ncbi-blast+/algo/blast/core/blast_posit.o \
			 include/ncbi-blast+/algo/blast/composition_adjustment/matrix_frequency_data.o \
			 include/ncbi-blast+/algo/blast/core/blast_dynarray.o \
			 include/ncbi-blast+/algo/blast/core/matrix_freq_ratios.o \
			 include/ncbi-blast+/algo/blast/core/blast_encoding.o \
			 include/ncbi-blast+/algo/blast/core/blast_stat.o \
			 include/ncbi-blast+/algo/blast/core/blast_filter.o \
			 include/ncbi-blast+/algo/blast/core/blast_util.o \
			 include/ncbi-blast+/algo/blast/core/blast_message.o \
			 include/ncbi-blast+/algo/blast/core/ncbi_erf.o \
			 include/ncbi-blast+/algo/blast/core/blast_options.o \
			 include/ncbi-blast+/algo/blast/core/ncbi_math.o \
			 include/ncbi-blast+/algo/blast/core/blast_program.o \
			 include/ncbi-blast+/algo/blast/core/ncbi_std.o \
			 include/ncbi-blast+/algo/blast/core/blast_psi_priv.o \
			 include/ncbi-blast+/util/tables/raw_scoremat.o \
			 include/ncbi-blast+/algo/blast/core/blast_query_info.o \
			 include/ncbi-blast+/algo/blast/core/blast_seg.o

BWTOBJS = bwt/bwt.o bwt/compactfmi.o bwt/sequence.o bwt/suffixArray.o

ifeq ($(uname -s), "Darwin")
LD_LIBS_STATIC = -Wl,-all_load -lpthread -lz -Wl,-noall_load
else
LD_LIBS_STATIC = -Wl,--whole-archive -lpthread -lz -Wl,--no-whole-archive
endif


all: makefile kaiju kaiju2krona mergeOutputs kaijuReport kaijux kaijup convertNR addTaxonNames bwt/mkbwt
	mkdir -p ../bin && cp kaiju kaijux kaijup kaiju2krona mergeOutputs kaijuReport convertNR addTaxonNames bwt/mkbwt bwt/mkfmi ../util/gbk2faa.pl ../util/makeDB.sh ../util/taxonlist.tsv ../util/convert_mar_to_kaiju.py ../bin/

# use bwt/mkbwt as target for compiling everything in the bwt folder first
bwt/mkbwt:
	$(MAKE) -C bwt/ $(MAKECMDGOALS)

kaiju: makefile bwt/mkbwt kaiju.o ReadItem.o Config.o ConsumerThread.o util.o $(BLASTOBJS)
	$(CXX) $(LDFLAGS) -o kaiju kaiju.o ReadItem.o Config.o ConsumerThread.o util.o $(BWTOBJS) $(BLASTOBJS) $(LDLIBS)

kaijux: makefile bwt/mkbwt kaijux.o ReadItem.o Config.o ConsumerThread.o ConsumerThreadx.o util.o $(BLASTOBJS)
	$(CXX) $(LDFLAGS) -o kaijux kaijux.o ReadItem.o Config.o ConsumerThread.o ConsumerThreadx.o util.o $(BWTOBJS) $(BLASTOBJS) $(LDLIBS)

kaijup: makefile bwt/mkbwt kaijup.o ReadItem.o Config.o ConsumerThread.o ConsumerThreadx.o ConsumerThreadp.o util.o $(BLASTOBJS)
	$(CXX) $(LDFLAGS) -o kaijup kaijup.o ReadItem.o Config.o ConsumerThread.o ConsumerThreadx.o ConsumerThreadp.o util.o $(BWTOBJS) $(BLASTOBJS) $(LDLIBS)

kaiju2krona: makefile bwt/mkbwt kaiju2krona.o util.o
	$(CXX) $(LDFLAGS) -o kaiju2krona kaiju2krona.o util.o $(BWTOBJS)

mergeOutputs: makefile bwt/mkbwt mergeOutputs.o util.o
	$(CXX) $(LDFLAGS) -o mergeOutputs mergeOutputs.o util.o $(BWTOBJS)

kaijuReport: makefile bwt/mkbwt kaijuReport.o util.o
	$(CXX) $(LDFLAGS) -o kaijuReport kaijuReport.o util.o $(BWTOBJS)

addTaxonNames: makefile bwt/mkbwt addTaxonNames.o util.o
	$(CXX) $(LDFLAGS) -o addTaxonNames addTaxonNames.o util.o $(BWTOBJS)

convertNR: makefile bwt/mkbwt Config.o convertNR.o util.o $(BLASTOBJS)
	$(CXX) $(LDFLAGS) -o convertNR convertNR.o Config.o util.o $(BWTOBJS) $(BLASTOBJS)


%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
%.o : %.cpp version.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<


clean:
	rm -f -v bwt/mkbwt bwt/mkfmi kaiju kaijux kaijup kaiju2krona mergeOutputs kaijuReport convertNR addTaxonNames ../bin/*
	find . -name "*.o" -delete
	$(MAKE) -C bwt/ clean

static: LDFLAGS = -static
static: LDLIBS = $(LD_LIBS_STATIC)
static: all

debug: CXXFLAGS = -O3 -pthread -std=c++11 -g -Wall -Wpedantic -Wextra -Wconversion -fno-omit-frame-pointer
debug: CFLAGS = -g -O3 -Wall -Wno-uninitialized
debug: all

.PHONY: clean debug static
