/* This file is part of Kaiju, Copyright 2015-2018 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <deque>
#include <stdexcept>

#include "zstr/zstr.hpp"
#include "ProducerConsumerQueue/src/ProducerConsumerQueue.hpp"

#include "ReadItem.hpp"
#include "ConsumerThreadp.hpp"
#include "Config.hpp"
#include "version.hpp"
#include "util.hpp"

extern "C" {
#include "bwt/bwt.h"
}


void usage(char *progname);

int main(int argc, char** argv) {


	Config * config = new Config();

	std::string fmi_filename;
	std::string in1_filename;
	std::string output_filename;

	int num_threads = 1;
	bool verbose = false;
	bool debug = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "a:hdxvn:m:e:E:l:f:i:s:z:o:")) != -1) {
		switch (c)  {
			case 'a': {
									if("mem" == std::string(optarg)) config->mode = MEM;
									else if("greedy" == std::string(optarg)) config->mode = GREEDY;
									else { std::cerr << "-a must be a valid mode.\n"; usage(argv[0]); }
									break;
								}
			case 'h':
				usage(argv[0]);
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'x':
				config->SEG = true; break;
			case 'o':
				output_filename = optarg; break;
			case 'f':
				fmi_filename = optarg; break;
			case 'i':
				in1_filename = optarg; break;
			case 'l': {
									try {
										int seed_length = std::stoi(optarg);
										if(seed_length < 7) { error("Seed length must be >= 7."); usage(argv[0]); }
										config->seed_length = (unsigned int)seed_length;
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
										int min_score = std::stoi(optarg);
										if(min_score <= 0) { error("Min Score (-s) must be greater than 0."); usage(argv[0]); }
										config->min_score = (unsigned int)min_score;
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
										int min_fragment_length = std::stoi(optarg);
										if(min_fragment_length <= 0) { error("Min fragment length (-m) must be greater than 0."); usage(argv[0]); }
										config->min_fragment_length = (unsigned int)min_fragment_length;
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
										int mismatches = std::stoi(optarg);
										if(mismatches < 0) { error("Number of mismatches must be >= 0."); usage(argv[0]); }
										config->mismatches = (unsigned int)mismatches;
									}
									catch(const std::invalid_argument& ia) {
										std::cerr << "Invalid numerical argument in -e " << optarg << std::endl;
									}
									catch (const std::out_of_range& oor) {
										std::cerr << "Invalid numerical argument in -e " << optarg << std::endl;
									}
									break;
								}
			case 'E': {
									try {
										config->min_Evalue = std::stod(optarg);
										if(config->min_Evalue <= 0.0) { error("E-value threshold must be greater than 0."); usage(argv[0]); }
										config->use_Evalue = true;
									}
									catch(const std::invalid_argument& ia) {
										std::cerr << "Invalid numerical argument in -E " << optarg << std::endl;
									}
									catch (const std::out_of_range& oor) {
										std::cerr << "Invalid numerical argument in -E " << optarg << std::endl;
									}
									break;
								}
			case 'z': {
									try {
										num_threads = std::stoi(optarg);
										if(num_threads <= 0) {  error("Number of threads (-z) must be greater than 0."); usage(argv[0]); }
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
	if(fmi_filename.length() == 0) { error("Please specify the location of the FMI file, using the -f option."); usage(argv[0]); }
	if(in1_filename.length() == 0) { error("Please specify the location of the input file, using the -i option."); usage(argv[0]); }
	if(config->use_Evalue && config->mode != GREEDY) { error("E-value calculation is only available in Greedy mode. Use option: -a greedy"); usage(argv[0]); }

	if(debug) {
		std::cerr << "Parameters: \n";
		std::cerr << "  minimum match length: " << config->min_fragment_length << "\n";
		std::cerr << "  minimum blosum62 score for matches: " << config->min_score << "\n";
		std::cerr << "  seed length for greedy matches: " << config->seed_length << "\n";
		if(config->use_Evalue)
			std::cerr << "  minimum E-value: " << config->min_Evalue << "\n";
		std::cerr << "  max number of mismatches within a match: "  << config->mismatches << "\n";
		std::cerr << "  run mode: "  << ((config->mode==MEM) ? "MEM" : "Greedy") << "\n";
		std::cerr << "  input file 1: " << in1_filename << "\n";
	}

	config->debug = debug;
	config->verbose = verbose;

	if(verbose) std::cerr << getCurrentTime() << " Reading database" << std::endl;

	readFMI(fmi_filename,config);

	config->init();

	if(output_filename.length()>0) {
		std::cerr << "Output file: " << output_filename << std::endl;
		std::ofstream * read2id_file = new std::ofstream();
		read2id_file->open(output_filename);
		if(!read2id_file->is_open()) {  error("Could not open file " + output_filename + " for writing"); exit(EXIT_FAILURE); }
		config->out_stream = read2id_file;
	}
	else {
		config->out_stream = &std::cout;
	}

	ProducerConsumerQueue<ReadItem*>* myWorkQueue = new ProducerConsumerQueue<ReadItem*>(500);
	std::deque<std::thread> threads;
	std::deque<ConsumerThreadp *> threadpointers;
	for(int i=0; i < num_threads; i++) {
		ConsumerThreadp* p = new ConsumerThreadp(myWorkQueue, config);
		threadpointers.push_back(p);
		threads.push_back(std::thread(&ConsumerThreadp::doWork,p));
	}

	zstr::ifstream* in1_file = NULL;
	try {
		in1_file = new zstr::ifstream(in1_filename);
		if(!in1_file->good()) {  error("Could not open file " + in1_filename); exit(EXIT_FAILURE); }
	}
	catch(std::exception e) {
		error("Could not open file " + in1_filename); exit(EXIT_FAILURE);
	}

	bool isFastQ = false;
	bool firstline = true;
	std::string line_from_file;
	line_from_file.reserve(2000);
	std::string name;
	std::string sequence;
	sequence.reserve(2000);

	if(verbose) std::cerr << getCurrentTime() << " Start search using " << num_threads << " threads." << std::endl;

	while(getline(*in1_file,line_from_file)) {
		if(firstline) {
			char fileTypeIdentifier = line_from_file[0];
			if(fileTypeIdentifier == '@') {
				isFastQ = true;
			}
			else if(fileTypeIdentifier != '>') {
				error("Auto-detection of file type for file " + in1_filename + " failed.");
				exit(EXIT_FAILURE);
			}
			firstline = false;
		}
		if(isFastQ) {
			// remove '@' from beginning of line
			line_from_file.erase(line_from_file.begin());
			name = line_from_file;
			// read sequence line
			getline(*in1_file,line_from_file);
			sequence = line_from_file;
			// skip + line
			in1_file->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
			// skip quality score line
			in1_file->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}
		else { //FASTA
			// remove '>' from beginning of line
			line_from_file.erase(line_from_file.begin());
			name = line_from_file;
			// read lines until next entry starts or file terminates
			sequence.clear();
			while(!(in1_file->peek()=='>' || in1_file->peek()==EOF)) {
				getline(*in1_file,line_from_file);
				sequence.append(line_from_file);
			}
		} // end FASTA

		strip(sequence); // remove non-alphabet chars

		myWorkQueue->push(new ReadItem(name, sequence));

	} // end main loop around file1

	myWorkQueue->pushedLast();

	delete in1_file;

	while(!threads.empty()) {
		threads.front().join();
		threads.pop_front();
		delete  threadpointers.front();
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
	return EXIT_SUCCESS;
}


void usage(char *progname) {
	fprintf(stderr, "Kaiju %s\n",KAIJUVERSION);
	fprintf(stderr, "Copyright 2015-2018 Peter Menzel, Anders Krogh\n");
	fprintf(stderr, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:\n   %s -f proteins.fmi -i reads.fastq\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -f FILENAME   Name of database file (.fmi) file\n");
	fprintf(stderr, "   -i FILENAME   Name of input file containing reads in FASTA or FASTQ format\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -o FILENAME   Name of output file. If not specified, output will be printed to STDOUT\n");
	fprintf(stderr, "   -z INT        Number of parallel threads for classification (default: 1)\n");
	fprintf(stderr, "   -a STRING     Run mode, either \"mem\"  or \"greedy\" (default: greedy)\n");
	fprintf(stderr, "   -e INT        Number of mismatches allowed in Greedy mode (default: 3)\n");
	fprintf(stderr, "   -m INT        Minimum match length (default: 11)\n");
	fprintf(stderr, "   -s INT        Minimum match score in Greedy mode (default: 65)\n");
	fprintf(stderr, "   -E FLOAT      Minimum E-value in Greedy mode\n");
	fprintf(stderr, "   -x            Enable SEG low complexity filter (enabled by default)\n");
	fprintf(stderr, "   -X            Disable SEG low complexity filter\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	//fprintf(stderr, "   -d            Enable debug output.\n");
	exit(EXIT_FAILURE);
}

