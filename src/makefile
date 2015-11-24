CXX = g++
CXXFLAGS = -Wall -Wextra -ansi -pedantic -O3 -pthread -std=c++11 -g
LDFLAGS = -pthread
INCLUDES	= -I./ProducerConsumerQueue/src

all: makefile bwt/mkbwt bwt/mkfmi kaiju kaiju2krona mergeOutputs kaijuReport
	mkdir -p ../bin && cp kaiju kaiju2krona mergeOutputs kaijuReport bwt/mkbwt bwt/mkfmi ../util/gbk2faa.pl ../bin/

bwt/mkbwt bwt/mkfmi bwt/bwt.o bwt/simplefmi.o:
	$(MAKE) -C bwt/

kaiju: makefile ./bwt/bwt.o ./bwt/simplefmi.o kaiju.o Config.o ReadItem.o Config.o ConsumerThread.o  
	$(CXX) -o kaiju $(LDFLAGS) kaiju.o ReadItem.o Config.o ConsumerThread.o ./bwt/bwt.o ./bwt/simplefmi.o 

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

ReadItem.o: ReadItem.cpp makefile
	$(CXX) $(CXXFLAGS) -c ReadItem.cpp

kaiju.o: kaiju.cpp makefile ProducerConsumerQueue/src/ProducerConsumerQueue.hpp ProducerConsumerQueue/src/ProducerConsumerQueue.tpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c kaiju.cpp

kaiju2krona.o: kaiju2krona.cpp makefile
	$(CXX) $(CXXFLAGS) -c kaiju2krona.cpp

mergeOutputs.o: mergeOutputs.cpp makefile
	$(CXX) $(CXXFLAGS) -c mergeOutputs.cpp

kaijuReport.o: kaijuReport.cpp makefile
	$(CXX) $(CXXFLAGS) -c kaijuReport.cpp

clean: 
	rm -f -v *.o kaiju kaiju2krona mergeOutputs kaijuReport ../bin/*
	$(MAKE) -C bwt/ clean

.PHONY: clean 
