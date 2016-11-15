CC = gcc
CXX = g++
CFLAGS = -g -O3 -DNDEBUG -Wall -Wno-uninitialized
CXXFLAGS = -ansi -pedantic -O3 -pthread -std=c++11 -DNDEBUG -g -Wall -Wconversion -Wno-unused-function
LDLIBS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
INCLUDES	= -I./include/ProducerConsumerQueue/src -I./include/ncbi-blast+

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


all: makefile kaiju kaiju2krona mergeOutputs kaijuReport kaijux kaijup convertNR addTaxonNames bwt/mkbwt
	mkdir -p ../bin && cp kaiju kaijux kaijup kaiju2krona mergeOutputs kaijuReport convertNR addTaxonNames bwt/mkbwt bwt/mkfmi ../util/gbk2faa.pl ../util/makeDB.sh ../util/taxonlist.tsv ../bin/

# use bwt/mkbwt as target for compiling everything in the bwt folder first
bwt/mkbwt:
	$(MAKE) -C bwt/ $(MAKECMDGOALS)

kaiju: makefile bwt/mkbwt kaiju.o ReadItem.o Config.o ConsumerThread.o util.o $(BLASTOBJS)
	$(CXX) $(LDFLAGS) -o kaiju kaiju.o ReadItem.o Config.o ConsumerThread.o util.o bwt/bwt.o bwt/compactfmi.o bwt/sequence.o bwt/suffixArray.o $(BLASTOBJS) $(LDLIBS)

kaijux: makefile bwt/mkbwt kaijux.o ReadItem.o Config.o ConsumerThread.o ConsumerThreadx.o util.o $(BLASTOBJS)
	$(CXX) $(LDFLAGS) -o kaijux kaijux.o ReadItem.o Config.o ConsumerThread.o ConsumerThreadx.o util.o bwt/bwt.o bwt/compactfmi.o bwt/sequence.o bwt/suffixArray.o $(BLASTOBJS) $(LDLIBS)

kaijup: makefile bwt/mkbwt kaijup.o ReadItem.o Config.o ConsumerThread.o ConsumerThreadx.o ConsumerThreadp.o util.o $(BLASTOBJS)
	$(CXX) $(LDFLAGS) -o kaijup kaijup.o ReadItem.o Config.o ConsumerThread.o ConsumerThreadx.o ConsumerThreadp.o util.o bwt/bwt.o bwt/compactfmi.o bwt/sequence.o bwt/suffixArray.o $(BLASTOBJS) $(LDLIBS)

kaiju2krona: makefile kaiju2krona.o util.o
	$(CXX) $(LDFLAGS) -o kaiju2krona kaiju2krona.o util.o

mergeOutputs: makefile mergeOutputs.o util.o
	$(CXX) $(LDFLAGS) -o mergeOutputs mergeOutputs.o util.o

kaijuReport: makefile kaijuReport.o util.o
	$(CXX) $(LDFLAGS) -o kaijuReport kaijuReport.o util.o

addTaxonNames: makefile addTaxonNames.o util.o
	$(CXX) $(LDFLAGS) -o addTaxonNames addTaxonNames.o util.o

convertNR: makefile Config.o convertNR.o bwt/mkbwt $(BLASTOBJS)
	$(CXX) $(LDFLAGS) -o convertNR convertNR.o Config.o util.o ./bwt/sequence.o $(BLASTOBJS)


%.o : %.c
	$(CC) $(CFLAGS) $(INCLUDES) -c -o $@ $<
%.o : %.cpp version.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<


clean:
	rm -f -v bwt/mkbwt bwt/mkfmi kaiju kaijux kaijup kaiju2krona mergeOutputs kaijuReport convertNR addTaxonNames ../bin/*
	find . -name "*.o" -delete
	$(MAKE) -C bwt/ clean

static: LDFLAGS = -static
static: all

.PHONY: clean funstatic
