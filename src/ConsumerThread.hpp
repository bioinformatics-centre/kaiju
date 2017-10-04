/* This file is part of Kaiju, Copyright 2015-2017 Peter Menzel and Anders Krogh,
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

#include "ReadItem.hpp"
#include "Config.hpp"
#include "util.hpp"

#include "ProducerConsumerQueue/src/ProducerConsumerQueue.hpp"
#include "algo/blast/core/blast_seg.h"
#include "algo/blast/core/blast_filter.h"
#include "algo/blast/core/blast_encoding.h"

extern "C" {
#include "bwt/bwt.h"
}

const double LN_2 = 0.6931471805;
/* using values for ungapped BLOSUM62 matrix from ncbi-blast+/algo/blast/core/blast_stat.c:263 */
const double LAMBDA = 0.3176;
const double LN_K = -2.009915479;  // K = 0.134

class Fragment {
	public:
	std::string seq = "";
	unsigned int num_mm = 0;
	int diff = 0;
	bool SEGchecked = false;
	unsigned int pos_lastmm = 0;
	IndexType si0, si1;
	int matchlen;
	Fragment(std::string s) { seq = s; }
	Fragment(std::string s, bool b) { seq = s; SEGchecked = true; }
	Fragment(std::string s, unsigned int n, unsigned int p, int d) {
		seq = s;
		num_mm=n;
		diff=d;
		pos_lastmm=p;
	}
	Fragment(std::string s, unsigned int n, unsigned int p, int d, IndexType arg_si0,IndexType arg_si1,int len) {
		seq = s;
		num_mm=n;
		diff=d;
		pos_lastmm=p;
		si0 = arg_si0;
		si1 = arg_si1;
		matchlen=len;
		SEGchecked = true; // fragments with substitutions have been checked before
	}
	Fragment(std::string s, unsigned int n, unsigned int p, int d, SI * si) {
		seq = s;
		num_mm=n;
		diff=d;
		pos_lastmm=p;
		si0 = si->start;
		si1 = si->start+(IndexType)si->len;
		matchlen=si->ql;
	}
	Fragment(std::string s, unsigned int n, unsigned int p, SI * si) {
		seq = s;
		num_mm=n;
		pos_lastmm=p;
		si0 = si->start;
		si1 = si->start+(IndexType)si->len;
		matchlen=si->ql;
	}
	Fragment(std::string s, unsigned int n, unsigned int p) {
		seq = s;
		num_mm=n;
		pos_lastmm=p;
	}
};

class ConsumerThread {
	protected:
	ProducerConsumerQueue<ReadItem*> * myWorkQueue;

	std::unordered_map<uint64_t,unsigned int> node2depth;

	uint8_t codon_to_int(const char* codon);
	uint8_t revcomp_codon_to_int(const char* codon);

	uint8_t nuc2int[256];
	uint8_t compnuc2int[256];
	char codon2aa[256];
	uint8_t aa2int[256];

	std::map<char, std::vector<char>> blosum_subst;
	int8_t blosum62diag[20];
	int8_t b62[20][20];

	std::string translations[6];
	std::multimap<unsigned int,Fragment *,std::greater<unsigned int>> fragments;
	std::vector<SI *> best_matches_SI;
	std::vector<SI *> longest_matches_SI;
	std::vector<std::string> best_matches;
	std::vector<std::string> longest_fragments;
	std::set<uint64_t> match_ids;
	std::set<std::string> match_dbnames;

	unsigned int best_match_score = 0;
	std::string extraoutput = "";

	double query_len;

	Config * config;
	std::ostringstream output;
	uint32_t read_count = 0;
	uint64_t classify_length();
	uint64_t classify_greedyblosum();

	void clearFragments();
	unsigned int calcScore(const std::string &);
	unsigned int calcScore(const std::string &, int);
	unsigned int calcScore(const std::string &, size_t, size_t, int);

	void addAllMismatchVariantsAtPosSI(const Fragment *,unsigned int, size_t, SI *); // used in Greedy mode
	Fragment * getNextFragment(unsigned int);

	void eval_match_scores(SI *si, Fragment *);
	void ids_from_SI_recursive(SI *si);
	void ids_from_SI(SI *si);
	void getAllFragmentsBits(const std::string & line);
	void flush_output();

	public:
	ConsumerThread(ProducerConsumerQueue<ReadItem*>* workQueue, Config * config);
	void doWork();


};
#endif

