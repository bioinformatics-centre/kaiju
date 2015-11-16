/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef CONSUMER_THREAD_H
#define CONSUMER_THREAD_H

#define NDEBUG

#include <unordered_map>
#include <unordered_set>
#include <list>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <mutex>
#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>
#include <string>
#include <cstring>
#include <climits>
#include <map>

#include "./ProducerConsumerQueue/src/ProducerConsumerQueue.hpp"
#include "ReadItem.hpp"
#include "Config.hpp"

using namespace std;

class Fragment {
	public:
	string seq = "";
	uint num_mm = 0;
	int diff = 0;
	uint pos_lastmm = 0;
	IndexType si0, si1;
	int matchlen;
	Fragment(string s) { seq = s; }
	Fragment(string s, uint n, uint p, int d) {
		seq = s; 
		num_mm=n; 
		diff=d; 
		pos_lastmm=p;
	}
	Fragment(string s, uint n, uint p, int d, SI * si) { 
		seq = s;
		num_mm=n;
		diff=d; 
		pos_lastmm=p;
		si0 = si->start;
		si1 = si->start+(IndexType)si->len;
		matchlen=si->ql;
	}
	Fragment(string s, uint n, uint p, SI * si) { 
		seq = s; 
		num_mm=n;
		pos_lastmm=p;
		si0 = si->start;
		si1 = si->start+(IndexType)si->len;
		matchlen=si->ql;
	}
	Fragment(string s, uint n, uint p) { 
		seq = s; 
		num_mm=n; 
		pos_lastmm=p;
	}
}; 

class ConsumerThread {
	ProducerConsumerQueue<ReadItem*> * myWorkQueue;        

	unordered_map<uint64_t,unsigned int> node2depth;

	uint8_t codon_to_int(const char* codon);
	uint8_t revcomp_codon_to_int(const char* codon);

	uint8_t nuc2int[256];   
	uint8_t compnuc2int[256];   
	char codon2aa[256];
	uint8_t aa2int[256];   

	unordered_map<char, vector<char>> blosum_subst;
	int8_t blosum62diag[20];
	int8_t b62[20][20];

	string translations[6];
	multimap<int,Fragment *,std::greater<uint>> fragments;
	vector<SI *> best_matches_SI; 
	vector<SI *> longest_matches_SI; 
	vector<string> best_matches; 
	vector<string> longest_fragments; 
	set<uint64_t> match_ids; 

	int best_match_score = 0;
	string extraoutput = "";

	Config * config;
	ostringstream output;
	uint64_t count = 0;
	uint64_t classify_length();
	uint64_t classify_greedyblosum();

	void clearFragments();
	int calcScore(const string &, int);
	int calcScore(const char *, uint, uint, int);

	void addAllMismatchVariantsAtPosSI(Fragment *,uint,size_t,SI *); // to be used by greedyblosum
	Fragment * getNextFragment(int);

	void eval_match_scores(SI *si, Fragment *);
	void ids_from_SI_recursive(SI *si);
	void ids_from_SI(SI *si);
	void getAllFragmentsBits(string & line);
	void flush_output();

	public:        
	ConsumerThread(ProducerConsumerQueue<ReadItem*>* workQueue, Config * config);        
	void doWork(); 

	
};
#endif

