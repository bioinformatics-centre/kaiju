/* This file is part of Kaiju, Copyright 2015-2019 Peter Menzel and Anders Krogh,
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


enum Mode { MEM, GREEDY };

class Config {
	public:
		Mode mode = GREEDY;
		size_t max_matches_SI = 20; // maximum number of best matches with same score, used for LCA and output
		size_t max_match_ids = 20; // maximum number of ids to print, used for LCA and output
		size_t max_match_acc = 20; // maximum number of accession numbers to print, used for LCA and output

		bool debug = false;
		bool verbose = false;

		bool SEG = true;
		bool input_is_protein = false;
		unsigned int min_fragment_length = 11; // in MEM and Greedy modes
		unsigned int mismatches = 3; // in Greedy mode
		unsigned int min_score = 65; // in Greedy mode
		unsigned int seed_length = 7; // in Greedy mode
		bool use_Evalue = true; // can only be used in Greedy mode
		double min_Evalue = 0.01; // can only be used in Greedy mode
		double db_length;

		SegParameters * blast_seg_params;

		std::ostream * out_stream;
		std::unordered_map<uint64_t,uint64_t> * nodes;

		FMI * fmi;
		BWT * bwt;

		AlphabetStruct * astruct;

		Config();
		~Config();

		void init();

};


#endif
