/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef CONFIG_H
#define CONFIG_H

#include <string.h>
#include <iostream>
#include <list>
#include <unordered_map>
#include <map>
#include <set>
#include <fstream>
#include <iterator>
#include <stdint.h>
#include <mutex>

extern "C" {
#include "./bwt/fmi.h"
#include "./bwt/bwt.h"
}

using namespace std;

enum Mode { MEM, GREEDYBLOSUM };

class Config {
	public:	
		Mode mode;
		bool debug;
		bool verbose;
		unsigned int mismatches;
		unsigned int min_fragment_length;
		int min_score; 
		unsigned int seed_length;
		ostream * out_stream;
		unordered_map<uint64_t,uint64_t> * nodes;

		std::mutex out_mutex;

		FMI * fmi;
		BWT * bwt;

		uchar * trans;

		Config();
		~Config();

		void init();
		uint64_t lca_from_ids(unordered_map<uint64_t,unsigned int> & node2depth, set<uint64_t> & ids);

};


#endif
