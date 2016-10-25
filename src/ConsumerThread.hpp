/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef CONSUMERTHREAD_H
#define CONSUMERTHREAD_H

#include <stdint.h>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <list>
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
#include <utility>
#include <functional>
#include <locale>

#include "include/ProducerConsumerQueue/src/ProducerConsumerQueue.hpp"
#include "ReadItem.hpp"
#include "Config.hpp"

#include "include/ncbi-blast+/algo/blast/core/blast_seg.h"
#include "include/ncbi-blast+/algo/blast/core/blast_filter.h"
#include "include/ncbi-blast+/algo/blast/core/blast_encoding.h"

extern "C" {
#include "./bwt/bwt.h"
}

using namespace std;

class Fragment {
	public:
	string seq = "";
	unsigned int num_mm = 0;
	int diff = 0;
	bool SEGchecked = false;
	unsigned int pos_lastmm = 0;
	IndexType si0, si1;
	int matchlen;
	Fragment(string s) { seq = s; }
	Fragment(string s, bool b) { seq = s; SEGchecked = true; }
	Fragment(string s, unsigned int n, unsigned int p, int d) {
		seq = s;
		num_mm=n;
		diff=d;
		pos_lastmm=p;
	}
	Fragment(string s, unsigned int n, unsigned int p, int d, IndexType arg_si0,IndexType arg_si1,int len) {
		seq = s;
		num_mm=n;
		diff=d;
		pos_lastmm=p;
		si0 = arg_si0;
		si1 = arg_si1;
		matchlen=len;
		SEGchecked = true; // fragments with substitutions have been checked before
	}
	Fragment(string s, unsigned int n, unsigned int p, int d, SI * si) {
		seq = s;
		num_mm=n;
		diff=d;
		pos_lastmm=p;
		si0 = si->start;
		si1 = si->start+(IndexType)si->len;
		matchlen=si->ql;
	}
	Fragment(string s, unsigned int n, unsigned int p, SI * si) {
		seq = s;
		num_mm=n;
		pos_lastmm=p;
		si0 = si->start;
		si1 = si->start+(IndexType)si->len;
		matchlen=si->ql;
	}
	Fragment(string s, unsigned int n, unsigned int p) {
		seq = s; 
		num_mm=n; 
		pos_lastmm=p;
	}
}; 

class ConsumerThread {
	protected:
	ProducerConsumerQueue<ReadItem*> * myWorkQueue;        

	unordered_map<uint64_t,unsigned int> node2depth;

	uint8_t codon_to_int(const char* codon);
	uint8_t revcomp_codon_to_int(const char* codon);

	uint8_t nuc2int[256];   
	uint8_t compnuc2int[256];   
	char codon2aa[256];
	uint8_t aa2int[256];   

	map<char, vector<char>> blosum_subst;
	int8_t blosum62diag[20];
	int8_t b62[20][20];

	string translations[6];
	multimap<unsigned int,Fragment *,std::greater<unsigned int>> fragments;
	vector<SI *> best_matches_SI; 
	vector<SI *> longest_matches_SI; 
	vector<string> best_matches; 
	vector<string> longest_fragments; 
	set<uint64_t> match_ids; 

	unsigned int best_match_score = 0;
	string extraoutput = "";

	Config * config;
	ostringstream output;
	uint32_t read_count = 0;
	uint64_t classify_length();
	uint64_t classify_greedyblosum();

	void clearFragments();
	unsigned int calcScore(const string &);
	unsigned int calcScore(const string &, int);
	unsigned int calcScore(const string &, size_t, size_t, int);

	void addAllMismatchVariantsAtPosSI(const Fragment *,unsigned int, size_t, SI *); // used in Greedy mode
	Fragment * getNextFragment(unsigned int);

	void eval_match_scores(SI *si, Fragment *);
	void ids_from_SI_recursive(SI *si);
	void ids_from_SI(SI *si);
	void getAllFragmentsBits(const string & line);
	void flush_output();

	public:        
	ConsumerThread(ProducerConsumerQueue<ReadItem*>* workQueue, Config * config);        
	void doWork(); 

	
};
#endif

