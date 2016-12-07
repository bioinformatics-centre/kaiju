/*************************************************
  Kaiju

  Authors: Peter Menzel <pmenzel@gmail.com> and
           Anders Krogh <krogh@binf.ku.dk>

  Copyright (C) 2015,2016 Peter Menzel and Anders Krogh

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  You should have received a copy of the GNU General Public License
  along with this program, see file LICENSE.
  If not, see <http://www.gnu.org/licenses/>.

  See the file README.md for documentation.
**************************************************/

#include <stdint.h>
#include <getopt.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <deque>
#include <stdexcept>

#include "include/ProducerConsumerQueue/src/ProducerConsumerQueue.hpp"
#include "ReadItem.hpp"
#include "ConsumerThread.hpp"
#include "Config.hpp"
#include "version.hpp"
#include "util.hpp"

extern "C" {
#include "./bwt/bwt.h"
}


void usage(char *progname);

int main(int argc, char** argv) {


	Config * config = new Config();

	std::unordered_map<uint64_t,uint64_t> * nodes = new std::unordered_map<uint64_t,uint64_t>();

	std::string nodes_filename = "";
	std::string fmi_filename = "";
	std::string sa_filename = "";
	std::string in1_filename = "";
	std::string in2_filename = "";
	std::string output_filename;

	Mode mode = MEM;

	int min_score = 65;
	int min_fragment_length = 11;
	int seed_length = 7;
	int mismatches  = 0;

	int num_threads = 1;
	bool verbose = false;
	bool debug = false;
	bool paired  = false;
	bool input_is_protein = false;
	bool SEG_check = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "a:hdpxvn:m:e:l:t:f:i:j:s:z:o:")) != -1) {
		switch (c)  {
			case 'a': {
									if("mem" == std::string(optarg)) mode = MEM;
									else if("greedy" == std::string(optarg)) mode = GREEDYBLOSUM;
									else { std::cerr << "-a must be a valid mode.\n"; usage(argv[0]); }
									break;
								}
			case 'h':
				usage(argv[0]);
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'p':
				input_is_protein = true; break;
			case 'x':
				SEG_check = true; break;
			case 'o':
				output_filename = optarg; break;
			case 'f':
				fmi_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'i':
				in1_filename = optarg; break;
			case 'j': {
									in2_filename = optarg;
									paired = true;
									break;
								}
			case 'l': {
									try {
										seed_length = std::stoi(optarg);
									}
									catch(const std::invalid_argument& ia) {
										std::cerr << "Invalid argument in -l " << optarg << std::endl;
									}
									catch (const std::out_of_range& oor) {
										std::cerr << "Invalid argument in -l " << optarg << std::endl;
									}
									break;
								}
			case 's': {
									try {
										min_score = std::stoi(optarg);
									}
									catch(const std::invalid_argument& ia) {
										std::cerr << "Invalid argument in -s " << optarg << std::endl;
									}
									catch (const std::out_of_range& oor) {
										std::cerr << "Invalid argument in -s " << optarg << std::endl;
									}
									break;
								}
			case 'm': {
									try {
										min_fragment_length = std::stoi(optarg);
									}
									catch(const std::invalid_argument& ia) {
										std::cerr << "Invalid argument in -m " << optarg << std::endl;
									}
									catch (const std::out_of_range& oor) {
										std::cerr << "Invalid argument in -m " << optarg << std::endl;
									}
									break;
								}
			case 'e': {
									try {
										mismatches = std::stoi(optarg);
									}
									catch(const std::invalid_argument& ia) {
										std::cerr << "Invalid numerical argument in -e " << optarg << std::endl;
									}
									catch (const std::out_of_range& oor) {
										std::cerr << "Invalid numerical argument in -e " << optarg << std::endl;
									}
									break;
								}
			case 'z': {
									try {
										num_threads = std::stoi(optarg);
									}
									catch(const std::invalid_argument& ia) {
										std::cerr << "Invalid argument in -z " << optarg << std::endl;
									}
									catch (const std::out_of_range& oor) {
										std::cerr << "Invalid argument in -z " << optarg << std::endl;
									}
									break;
								}
			default:
								usage(argv[0]);
		}
	}
	if(min_score <= 0) { error("Min Score (-s) must be greater than 0."); usage(argv[0]); }
	if(num_threads <= 0) {  error("Number of threads (-z) must be greater than 0."); usage(argv[0]); }
	if(min_fragment_length <= 0) { error("Min fragment length (-m) must be greater than 0."); usage(argv[0]); }
	if(mismatches < 0) { error("Number of mismatches must be >= 0."); usage(argv[0]); }
	if(seed_length < 7) { error("Seed length must be >= 7."); usage(argv[0]); }
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file, using the -t option."); usage(argv[0]); }
	if(fmi_filename.length() == 0) { error("Please specify the location of the FMI file, using the -f option."); usage(argv[0]); }
	if(in1_filename.length() == 0) { error("Please specify the location of the input file, using the -i option."); usage(argv[0]); }
	if(paired && input_is_protein) { error("Protein input only supports one input file."); usage(argv[0]); }

	if(debug) {
		std::cerr << "Parameters: \n";
		std::cerr << "  minimum fragment length for matches: " << min_fragment_length << "\n";
		std::cerr << "  minimum blosum score for matches: " << min_score << "\n";
		std::cerr << "  max number of mismatches within a match: "  << mismatches << "\n";
		std::cerr << "  run mode: "  << mode << "\n";
		std::cerr << "  input file 1: " << in1_filename << "\n";
		if(in2_filename.length() > 0)
			std::cerr << "  input file 2: " << in2_filename << "\n";
	}

	config->mode = mode;
	config->nodes = nodes;
	config->debug = debug;
	config->verbose = verbose;
	config->input_is_protein = input_is_protein;
	config->min_score = min_score;
	config->min_fragment_length = min_fragment_length;
	config->seed_length = seed_length;
	config->mismatches = mismatches;
	config->SEG = SEG_check;

	if(verbose) std::cerr << getCurrentTime() << " Reading database" << std::endl;

	std::ifstream nodes_file;
	nodes_file.open(nodes_filename.c_str());
	if(!nodes_file.is_open()) { std::cerr << "Error: Could not open file " << nodes_filename << std::endl; usage(argv[0]); }
	if(verbose) std::cerr << " Reading taxonomic tree from file " << nodes_filename << std::endl;
	parseNodesDmp(*nodes,nodes_file);
	nodes_file.close();

	{
	if(verbose) std::cerr << " Reading index from file " << fmi_filename << std::endl;
	FILE * fp = fopen(fmi_filename.c_str(),"r");
	if (!fp) { std::cerr << "Could not open file " << fmi_filename << std::endl; usage(argv[0]); }
	BWT * b = readIndexes(fp);
	fclose(fp);
	if(debug) fprintf(stderr,"BWT of length %ld has been read with %d sequences, alphabet=%s\n", b->len,b->nseq, b->alphabet);
	config->bwt = b;
	config->fmi = b->f;
	}

	config->init();

	if(output_filename.length()>0) {
		if(verbose) std::cerr << "Output file: " << output_filename << std::endl;
		std::ofstream * read2id_file = new std::ofstream();
		read2id_file->open(output_filename);
		if(!read2id_file->is_open()) {  std::cerr << "Could not open file " << output_filename << " for writing" << std::endl; exit(EXIT_FAILURE); }
		config->out_stream = read2id_file;
	}
	else {
		config->out_stream = &std::cout;
	}

	ProducerConsumerQueue<ReadItem*>* myWorkQueue = new ProducerConsumerQueue<ReadItem*>(500);
	std::deque<std::thread> threads;
	std::deque<ConsumerThread *> threadpointers;
	for(int i=0; i < num_threads; i++) {
		ConsumerThread * p = new ConsumerThread(myWorkQueue, config);
		threadpointers.push_back(p);
		threads.push_back(std::thread(&ConsumerThread::doWork,p));
	}

	std::ifstream in1_file, in2_file;
	in1_file.open(in1_filename);

	if(!in1_file.is_open()) {  std::cerr << "Could not open file " << in1_filename << std::endl; exit(EXIT_FAILURE); }
	if(in2_filename.length() > 0) {
		in2_file.open(in2_filename);
		if(!in2_file.is_open()) {  std::cerr << "Could not open file " << in2_filename << std::endl; exit(EXIT_FAILURE); }
	}

	bool firstline_file1 = true;
	bool firstline_file2 = true;
	bool isFastQ_file1 = false;
	bool isFastQ_file2 = false;
	std::string line_from_file;
	line_from_file.reserve(2000);
	std::string suffixStartCharacters = " /\t";
	std::string name;
	std::string sequence1;
	std::string sequence2;
	sequence1.reserve(2000);
	if(paired) sequence2.reserve(2000);

	if(verbose) std::cerr << getCurrentTime() << " Start classification using " << num_threads << " threads." << std::endl;

	while(getline(in1_file,line_from_file)) {
		if(line_from_file.length() == 0) { continue; }
		if(firstline_file1) {
			char fileTypeIdentifier = line_from_file[0];
			if(fileTypeIdentifier == '@') {
				isFastQ_file1 = true;
			}
			else if(fileTypeIdentifier != '>') {
				std::cerr << "Auto-detection of file type for file " << in1_filename << " failed."  << std::endl;
				exit(EXIT_FAILURE);
			}
			firstline_file1 = false;
		}
		if(isFastQ_file1) {
			// remove '@' from beginning of line
			line_from_file.erase(line_from_file.begin());
			// delete suffixes like '/1' or ' 1:N:0:TAAGGCGA' from end of read name
			size_t n = line_from_file.find_first_of(suffixStartCharacters);
			if(n != std::string::npos) { line_from_file.erase(n); }
			name = line_from_file;
			// read sequence line
			getline(in1_file,line_from_file);
			sequence1 = line_from_file;
			// skip + lin
			in1_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			// skip quality score line
			in1_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		else { //FASTA
			// remove '>' from beginning of line
			line_from_file.erase(line_from_file.begin());
			// delete suffixes like '/1' or ' 1:N:0:TAAGGCGA' from end of read name
			size_t n = line_from_file.find_first_of(suffixStartCharacters);
			if(n != std::string::npos) { line_from_file.erase(n); }
			name = line_from_file;
			// read lines until next entry starts or file terminates
			sequence1.clear();
			while(!(in1_file.peek()=='>' || in1_file.peek()==EOF)) {
				getline(in1_file,line_from_file);
				sequence1.append(line_from_file);
			}
		} // end FASTA

		strip(sequence1); // remove non-alphabet chars

		if(paired) {
			line_from_file = "";
			while(line_from_file.length() == 0) {
				if(!getline(in2_file,line_from_file)) {
					//that's the border case where file1 has more entries than file2
					std::cerr << "Error: File " << in1_filename <<" contains more reads then file " << in2_filename  <<std::endl;
					in1_file.close();
					in2_file.close();
					exit(EXIT_FAILURE);
				}
			}
			if(firstline_file2) {
				char fileTypeIdentifier = line_from_file[0];
				if(fileTypeIdentifier == '@') {
					isFastQ_file2 = true;
				}
				else if(fileTypeIdentifier != '>') {
					std::cerr << "Auto-detection of file type for file " << in2_filename << " failed."  << std::endl;
					exit(EXIT_FAILURE);
				}
				firstline_file2 = false;
			}
			if(isFastQ_file2) {
				// remove '@' from beginning of line
				line_from_file.erase(line_from_file.begin());
				// delete suffixes like '/2' or ' 2:N:0:TAAGGCGA' from end of read name
				size_t n = line_from_file.find_first_of(suffixStartCharacters);
				if(n != std::string::npos) { line_from_file.erase(n); }
				if(name != line_from_file) {
					std::cerr << "Error: Read names are not identical between the two input files. Probably reads are not in the same order in both files." << std::endl;
					in1_file.close();
					in2_file.close();
					exit(EXIT_FAILURE);
				}
				// read sequence line
				getline(in2_file,line_from_file);
				sequence2 = line_from_file;
				// skip + line
				in2_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
				// skip quality score line
				in2_file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			}
			else { // FASTA
				// remove '>' from beginning of line
				line_from_file.erase(line_from_file.begin());
				// delete suffixes like '/2' or ' 2:N:0:TAAGGCGA' from end of read name
				size_t n = line_from_file.find_first_of(suffixStartCharacters);
				if(n != std::string::npos) { line_from_file.erase(n); }
				if(name != line_from_file) {
					std::cerr << "Error: Read names are not identical between the two input files" << std::endl;
					in1_file.close();
					in2_file.close();
					exit(EXIT_FAILURE);
				}
				sequence2.clear();
				while(!(in2_file.peek()=='>' || in2_file.peek()==EOF)) {
					getline(in2_file,line_from_file);
					sequence2.append(line_from_file);
				}
			}
			strip(sequence2); // remove non-alphabet chars
			myWorkQueue->push(new ReadItem(name, sequence1, sequence2));
		} // not paired
		else {
			myWorkQueue->push(new ReadItem(name, sequence1));
		}

	} // end main loop around file1

	myWorkQueue->pushedLast();

	if(in1_file.is_open()) in1_file.close();

	if(paired && in2_file.is_open()) {
		if(getline(in2_file,line_from_file) && line_from_file.length()>0) {
			std::cerr << "Warning: File " << in2_filename <<" has more reads then file " << in1_filename  <<std::endl;
		}
		in2_file.close();
	}

	while(!threads.empty()) {
		threads.front().join();
		threads.pop_front();
		delete threadpointers.front();
		threadpointers.pop_front();
	}
	if(verbose) std::cerr << getCurrentTime() << " Finished." << std::endl;

	config->out_stream->flush();
	if(output_filename.length()>0) {
		((std::ofstream*)config->out_stream)->close();
		delete ((std::ofstream*)config->out_stream);
	}

	delete myWorkQueue;
	delete config;
	delete nodes;
	return EXIT_SUCCESS;
}

void usage(char *progname) {
	fprintf(stderr, "Kaiju %s\n",KAIJUVERSION);
	fprintf(stderr, "Copyright 2015,2016 Peter Menzel, Anders Krogh\n");
	fprintf(stderr, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -f kaiju_db.fmi -i reads.fastq [-j reads2.fastq]\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file\n");
	fprintf(stderr, "   -f FILENAME   Name of database (.fmi) file\n");
	fprintf(stderr, "   -i FILENAME   Name of input file containing reads in FASTA or FASTQ format\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -j FILENAME   Name of second input file for paired-end reads\n");
	fprintf(stderr, "   -o FILENAME   Name of output file. If not specified, output will be printed to STDOUT\n");
	fprintf(stderr, "   -z INT        Number of parallel threads (default: 1)\n");
	fprintf(stderr, "   -a STRING     Run mode, either \"mem\"  or \"greedy\" (default: mem)\n");
	fprintf(stderr, "   -e INT        Number of mismatches allowed in Greedy mode (default: 0)\n");
	fprintf(stderr, "   -m INT        Minimum match length (default: 11)\n");
	fprintf(stderr, "   -s INT        Minimum match score in Greedy mode (default: 65)\n");
	fprintf(stderr, "   -x            Enable SEG low complexity filter\n");
	fprintf(stderr, "   -p            Input sequences are protein sequences\n");
	fprintf(stderr, "   -v            Enable verbose output\n");
	//fprintf(stderr, "   -d            Enable debug output.\n");
	exit(EXIT_FAILURE);
}

