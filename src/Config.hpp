/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
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

#include "include/ncbi-blast+/algo/blast/core/blast_seg.h"
#include "include/ncbi-blast+/algo/blast/core/blast_filter.h"
#include "include/ncbi-blast+/algo/blast/core/blast_encoding.h"

extern "C" {
#include "./bwt/fmi.h"
#include "./bwt/bwt.h"
#include "./bwt/sequence.h"
}


enum Mode { MEM, GREEDYBLOSUM };

class Config {
	public:
		Mode mode;
		bool debug;
		bool SEG = false;
		SegParameters * blast_seg_params;
		bool verbose;
		bool input_is_protein;
		unsigned int mismatches;
		unsigned int min_fragment_length;
		unsigned int min_score;
		unsigned int seed_length;
		std::ostream * out_stream;
		std::unordered_map<uint64_t,uint64_t> * nodes;

		std::mutex out_mutex;

		FMI * fmi;
		BWT * bwt;

		AlphabetStruct * astruct;

		Config();
		~Config();

		void init();
		uint64_t lca_from_ids(std::unordered_map<uint64_t,unsigned int> & node2depth, std::set<uint64_t> & ids);

};


#endif
