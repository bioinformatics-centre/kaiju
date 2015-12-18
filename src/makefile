CXX = g++
CXXFLAGS = -ansi -pedantic -O3 -pthread -std=c++11 -g
LDLIBS = -lpthread
INCLUDES	= -I./ProducerConsumerQueue/src

all: makefile bwt/mkbwt bwt/mkfmi kaiju kaiju2krona mergeOutputs kaijuReport kaijux
	mkdir -p ../bin && cp kaiju kaijux kaiju2krona mergeOutputs kaijuReport bwt/mkbwt bwt/mkfmi ../util/gbk2faa.pl ../bin/

bwt/mkbwt bwt/mkfmi bwt/bwt.o bwt/compactfmi.o bwt/sequence.o bwt/suffixArray.o:
	$(MAKE) -C bwt/

kaiju: makefile ./bwt/bwt.o ./bwt/compactfmi.o ./bwt/sequence.o ./bwt/suffixArray.o kaiju.o ReadItem.o Config.o ConsumerThread.o  
	$(CXX) -o kaiju kaiju.o ReadItem.o Config.o ConsumerThread.o ./bwt/bwt.o ./bwt/compactfmi.o ./bwt/sequence.o ./bwt/suffixArray.o $(LDLIBS) 

kaijux: makefile ./bwt/bwt.o ./bwt/compactfmi.o ./bwt/sequence.o ./bwt/suffixArray.o kaijux.o ReadItem.o Config.o ConsumerThread.o  ConsumerThreadx.o 
	$(CXX) -o kaijux kaijux.o ReadItem.o Config.o ConsumerThread.o ConsumerThreadx.o ./bwt/bwt.o ./bwt/compactfmi.o ./bwt/sequence.o ./bwt/suffixArray.o $(LDLIBS) 

kaiju2krona: makefile kaiju2krona.o
	$(CXX) -o kaiju2krona kaiju2krona.o

mergeOutputs: makefile mergeOutputs.o
	$(CXX) -o mergeOutputs mergeOutputs.o

kaijuReport: makefile kaijuReport.o
	$(CXX) -o kaijuReport kaijuReport.o

Config.o: Config.cpp Config.hpp
	$(CXX) $(CXXFLAGS) -c Config.cpp

ConsumerThread.o: ConsumerThread.cpp ConsumerThread.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ConsumerThread.cpp

ConsumerThreadx.o: ConsumerThreadx.cpp ConsumerThreadx.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c ConsumerThreadx.cpp

ReadItem.o: ReadItem.cpp makefile
	$(CXX) $(CXXFLAGS) -c ReadItem.cpp

kaiju.o: kaiju.cpp makefile ProducerConsumerQueue/src/ProducerConsumerQueue.hpp ProducerConsumerQueue/src/ProducerConsumerQueue.tpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c kaiju.cpp

kaijux.o: kaijux.cpp makefile ProducerConsumerQueue/src/ProducerConsumerQueue.hpp ProducerConsumerQueue/src/ProducerConsumerQueue.tpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c kaijux.cpp

kaiju2krona.o: kaiju2krona.cpp makefile
	$(CXX) $(CXXFLAGS) -c kaiju2krona.cpp

mergeOutputs.o: mergeOutputs.cpp makefile
	$(CXX) $(CXXFLAGS) -c mergeOutputs.cpp

kaijuReport.o: kaijuReport.cpp makefile
	$(CXX) $(CXXFLAGS) -c kaijuReport.cpp

clean: 
	rm -f -v *.o kaiju kaijux kaiju2krona mergeOutputs kaijuReport ../bin/*
	$(MAKE) -C bwt/ clean

.PHONY: clean 
