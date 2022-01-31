/* This file is part of Kaiju, Copyright 2015-2017 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <climits>

#include "zstr/zstr.hpp"
#include "Config.hpp"
#include "util.hpp"

void usage(char *progname);

int main(int argc, char **argv) {

	Config * config = new Config();

	std::unordered_map<uint64_t,uint64_t> * nodes = new std::unordered_map<uint64_t,uint64_t>();
	std::unordered_map<uint64_t,uint64_t> * merged = new std::unordered_map<uint64_t,uint64_t>();

	std::unordered_map<uint64_t,unsigned int> node2depth;
	std::unordered_map<std::string,uint64_t> acc2taxid;
	std::unordered_set<std::string> excluded_accessions;

	std::string nodes_filename;
	std::string merged_filename;
	std::string list_filename;
	std::string acc_taxid_filename;
	std::string excluded_accession_filename;
	std::string nr_filename;
	std::string out_filename;

	//uint64_t include_ids [3] = { 2, 2157, 10239 };  //Bacteria, Archaea, Viruses
	std::unordered_set<uint64_t> include_ids;
	bool verbose = false;
	bool debug = false;
	bool addAcc = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "ahdvrm:l:g:t:i:o:e:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'a':
				addAcc = true; break;
			case 'e':
				excluded_accession_filename = optarg; break;
			case 'l':
				list_filename = optarg; break;
			case 'm':
				merged_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'g':
				acc_taxid_filename = optarg; break;
			case 'i':
				nr_filename = optarg; break;
			case 'o':
				out_filename = optarg; break;
			default:
				usage(argv[0]);
		}
	}
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file, using the -t option."); usage(argv[0]); }
	if(merged_filename.length() == 0) { error("Please specify the location of the merged.dmp file, using the -m option."); usage(argv[0]); }
	if(acc_taxid_filename.length() == 0) { error("Please specify the location of the prot.accession2taxid file, using the -g option."); usage(argv[0]); }
	if(out_filename.length() == 0) { error("Please specify the name of the output file, using the -o option."); usage(argv[0]); }

	config->nodes = nodes;
	config->debug = debug;
	config->verbose = verbose;

	std::ifstream nodes_file;
	nodes_file.open(nodes_filename.c_str());
	if(!nodes_file.is_open()) { error("Could not open file " + nodes_filename); exit(EXIT_FAILURE); }
	std::cerr << getCurrentTime() << " Reading taxonomic tree from file " << nodes_filename << std::endl;
	parseNodesDmp(*nodes, nodes_file);
	nodes_file.close();

	std::ifstream merged_file;
	merged_file.open(merged_filename.c_str());
	if(!merged_file.is_open()) { error("Could not open file " + merged_filename); exit(EXIT_FAILURE); }
	std::cerr << getCurrentTime() << " Reading file " << merged_filename << std::endl;
	parseMergedDmp(*merged, merged_file);
	merged_file.close();

	if(list_filename.length()==0) {
		std::cerr << "No taxa list specified, using Archaea, Bacteria, and Viruses." << std::endl;
		include_ids.insert((uint64_t)2);
		include_ids.insert((uint64_t)2157);
		include_ids.insert((uint64_t)10239);
	}
	else {
		std::ifstream list_file;
		list_file.open(list_filename);
		if(!list_file.is_open()) { error("Could not open file " + list_filename); exit(EXIT_FAILURE); }
		std::cerr << getCurrentTime() << " Reading taxa from file " << list_filename << std::endl;
		std::string line;
		while(getline(list_file, line)) {
			if(line.length() == 0) { continue; }
			size_t start = line.find_first_of("0123456789");
			if(start == std::string::npos) {
				continue;
			}
			size_t end = line.find_first_not_of("0123456789",start+1);
			if(end == std::string::npos) {
				end = start + line.length()-start;
			}
			try {
				uint64_t taxid = stoul(line.substr(start,end-start));
				if(debug)	std::cerr << "Found taxon id " << taxid << ", start=" << start << " end = " <<end<< std::endl;
				if(nodes->count(taxid) > 0) {
					include_ids.insert(taxid);
				}
				else {
					std::cerr << "Warning: Taxon ID " << taxid << " was not found in taxonomic tree. Skipping." << std::endl;
				}
			}
			catch(const std::invalid_argument& ia) {
				std::cerr << "Error: Found bad taxon id in line: " << line << std::endl;
			}
			catch (const std::out_of_range& oor) {
				std::cerr << "Error: Found bad number (out of range error) in line: " << line << std::endl;
			}
		}
		list_file.close();
	}

	acc2taxid.reserve(400e6);

	zstr::ifstream* acc_taxid_file = nullptr;
	try {
		acc_taxid_file = new zstr::ifstream(acc_taxid_filename);
		if(!acc_taxid_file->good()) {  error("Could not open file " + acc_taxid_filename); exit(EXIT_FAILURE); }
	} catch(std::exception e) { error("Could not open file " + acc_taxid_filename); exit(EXIT_FAILURE); }
	std::cerr << getCurrentTime() << " Reading accession to taxon id map from file " << acc_taxid_filename << std::endl;
	std::string line;
	getline(*acc_taxid_file, line); // skip header line
	while(getline(*acc_taxid_file, line)) {
		if(line.length() == 0) { continue; }
		try {
			size_t start = line.find('\t',0);
			size_t end = line.find("\t",start+1);
			//std::string acc = line.substr(start+1,end-start-1);
			uint64_t taxid = strtoul(line.c_str() + end + 1,NULL,10);
			if(taxid == ULONG_MAX) {
				std::cerr << "Found bad taxid number (out of range error) in line: " << line << std::endl;
				continue;
			}
			// check if taxon id is in nodes.dmp, otherwise check merged.dmp
			if(nodes->count(taxid) == 0) {
				if(merged->count(taxid) > 0) {
					if(verbose) std::cerr << "Taxon ID " << taxid << " for accession " << line.substr(start+1,end-start-1) << " was replaced by " << merged->at(taxid) << "\n";
					taxid = merged->at(taxid);
					if(nodes->count(taxid) == 0) {
						if(verbose) std::cerr << "Taxon ID " << taxid << " was not found in nodes.dmp\n";
					}
					else {
						acc2taxid.emplace(line.substr(start+1,end-start-1), taxid);
					}
				}
				else {
					if(verbose) std::cerr << "Taxon ID " << taxid << " for accession " << line.substr(start+1,end-start-1) << " was not found in nodes.dmp and merged.dmp\n";
				}
			}
			else { // taxid was found in nodes.dmp, so add accession to acc2taxid
				//std::cerr << acc << "\t" << taxid << "\n";
				acc2taxid.emplace(line.substr(start+1,end-start-1), taxid);
			}
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Found bad identifier in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
		}
	}
	delete acc_taxid_file;

	if(excluded_accession_filename.length()>0) {
		std::ifstream list_file;
		list_file.open(excluded_accession_filename);
		if(!list_file.is_open()) { error("Could not open file " + excluded_accession_filename); exit(EXIT_FAILURE); }
		std::cerr << getCurrentTime() << " Reading excluded accessions from file " << excluded_accession_filename << std::endl;
		std::string line;
		while(getline(list_file, line)) {
			if(line.length() == 0) { continue; }
			excluded_accessions.emplace(line);
		}
	}

	std::ifstream inputfile;
	if(nr_filename.length()>0) {
		inputfile.open(nr_filename);
		if(!inputfile.is_open()) { error("Could not open file " + nr_filename); exit(EXIT_FAILURE); }
	}

	if(verbose) std::cerr << "Writing to file " << out_filename << std::endl;
	std::ofstream out_file;
	out_file.open(out_filename);
	if(!out_file.is_open()) {  error("Could not open file " + out_filename + " for writing."); exit(EXIT_FAILURE); }

	std::cerr << getCurrentTime() << " Processing NR file " << nr_filename << std::endl;

	bool skip = true;
	bool first = true;
	std::ostringstream output;
	uint64_t outlinecount = 0;
	std::set<uint64_t> ids;
	while(getline(inputfile.is_open() ? inputfile : std::cin, line)){
		if(line.length() == 0) { continue; }
		if(line[0]=='>') {
			std::string first_acc;
			ids.clear();
			if(debug) std::cerr << "processing line " << line << std::endl;
			skip = false;
			size_t start = 1, end = 0;
			while((end = line.find(' ',start)) != std::string::npos) {
				// acc is between start and end
				std::string acc = line.substr(start, end - start);
				if(excluded_accessions.count(acc) > 0) {
					skip = true;
					break;
				}
				auto pos = acc2taxid.find(acc);
				if(pos != acc2taxid.end() && pos->second > 0) {
					if(addAcc && first_acc.empty()) { first_acc = acc; } // use first Accession that has taxon id as first part of DB identifier
					if(debug) std::cerr << "Accession " << acc << " belongs to taxon id " << pos->second  << std::endl;
					ids.insert(pos->second);
				}
				else if(verbose) { std::cerr << "Accession " << acc <<" was either not found in " << acc_taxid_filename << " or in " << nodes_filename << "\n"; }
				//look for next ID
				if((start = line.find("\x01",end+1)) == std::string::npos) { // no more entries
					break;
				}
				else {
					start++;
				}
			}
			if(skip) {
				if(verbose) std::cerr << "Found excluded accession nr. on line " << line << "\n";
				continue;
			}

			skip = true;
			if(!ids.empty()) {
				bool keep = false;
				uint64_t lca = (ids.size()==1) ?  *(ids.begin()) : lca_from_ids(config, node2depth, ids);
				if(debug) std::cerr << "LCA=" << lca << std::endl;
				if(nodes->count(lca)==0) {
					if(verbose) std::cerr << "Taxon ID " << lca << " not found in taxonomy!" << std::endl;
					continue;
				}
				uint64_t id = lca;
				while(nodes->count(id)>0 && id != 1) {
					if(include_ids.count(id) > 0) {
						keep = true;
						break;
					}
					id = nodes->at(id);
				}
				if(keep) {
					if(!first) { output << "\n";  } else { first = false; }
					output << ">";
					if(addAcc) output << first_acc << "_";
					output << lca << "\n";
					skip = false;
					outlinecount++;
				}
				else if(debug) { std::cerr << "Skipping sequence starting with taxon id " << *(ids.begin()) << std::endl; }
			}
			else if(verbose) std::cerr << "Could not find any taxonomy id for line " << line << "\n";
		}
		else {
			if(!skip) {
				size_t p = 0;
				while((p = line.find_first_of("ARNDCQEGHILKMFPSTWYV",p)) != std::string::npos) {
					output << line[p];
					p++;
				}
				outlinecount++;
			}
		}
		if(outlinecount%40000==0) {
			out_file << output.str();
			output.str("");
		}
	}
	output << std::endl;
	out_file << output.str();
	if(inputfile.is_open())
		inputfile.close();
	out_file.close();

	std::cerr << getCurrentTime() << " Finished." << std::endl;
	return EXIT_SUCCESS;
}

void usage(char *progname) {
	print_usage_header();
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -g prot.accession2taxid -i nr\n", progname);
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file.\n");
	fprintf(stderr, "   -g FILENAME   Name of prot.accession2taxid file.\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -a            Prefix taxon ID in database names with the first accession number per record.\n");
	fprintf(stderr, "   -i FILENAME   Name of NR file. If this option is not used, then the program will read from STDIN.\n");
	fprintf(stderr, "   -l FILENAME   Name of file with taxon IDs. Only records having one of these IDs as ancestor in the taxonomy will be used.\n");
	fprintf(stderr, "   -e FILENAME   Name of file with accession numbers that will be excluded.\n");
	exit(EXIT_FAILURE);
}

