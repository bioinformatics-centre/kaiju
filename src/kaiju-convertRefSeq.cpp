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

	std::string nodes_filename;
	std::string merged_filename;
	std::string list_filename;
	std::string acc_taxid_filename;
	std::string out_filename;

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
			case 'l':
				list_filename = optarg; break;
			case 'm':
				merged_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'g':
				acc_taxid_filename = optarg; break;
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

	acc2taxid.reserve(300e6);

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
			if((start = line.find('\t', 0)) == std::string::npos) {
				std::cerr << "Error parsing line: " << line << std::endl;
				continue;
			}
			//std::string acc = line.substr(start+1,end-start-1);
			if(line.substr(0, 3) != "WP_") {
				continue;
			}
			uint64_t taxid = strtoul(line.c_str() + start + 1,NULL,10);
			if(taxid == ULONG_MAX) {
				std::cerr << "Found bad taxid number (out of range error) in line: " << line << std::endl;
				continue;
			}
			if(taxid == 0) {
				if(debug) std::cerr << "Taxon ID is zero in line: " << line << std::endl;
				continue;
			}
			// check if taxon id is in nodes.dmp, otherwise check merged.dmp
			if(nodes->count(taxid) == 0) {
				if(merged->count(taxid) > 0) {
					if(debug) std::cerr << "Taxon ID " << taxid << " for accession " << line.substr(0, start) << " was replaced by " << merged->at(taxid) << "\n";
					taxid = merged->at(taxid);
					if(nodes->count(taxid) == 0) {
						if(verbose) std::cerr << "Taxon ID " << taxid << " was not found in nodes.dmp\n";
					}
					else {
						acc2taxid.emplace(line.substr(0, start-1), taxid);
					}
				}
				else {
					if(verbose) std::cerr << "Taxon ID " << taxid << " for accession " << line.substr(0, start) << " was not found in nodes.dmp and merged.dmp\n";
				}
			}
			else { // taxid was found in nodes.dmp, so add accession to acc2taxid
				//std::cerr << acc << "\t" << taxid << "\n";
				acc2taxid.emplace(line.substr(0, start), taxid);
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

	if(verbose) std::cerr << "Writing to file " << out_filename << std::endl;
	std::ofstream out_file;
	out_file.open(out_filename);
	if(!out_file.is_open()) {  error("Could not open file " + out_filename + " for writing."); exit(EXIT_FAILURE); }

	std::cerr << getCurrentTime() << " Processing STDIN " << std::endl;

	bool skip = true;
	bool first = true;
	std::ostringstream output;
	uint64_t outlinecount = 0;
	while(getline(std::cin, line)){
		if(line.length() == 0) { continue; }
		if(line[0]=='>') {
			if(debug) std::cerr << "processing line " << line << std::endl;
			uint64_t tax_id = 0;
			std::string acc = "";
			skip = true;
			size_t start = 1, end = 0;
			// lookg for space
			if((end = line.find(' ',start)) != std::string::npos) {
				// acc is between start and end
				acc = line.substr(start, end - start);
				auto pos = acc2taxid.find(acc);
				if(pos != acc2taxid.end() && pos->second > 0) {
					if(debug) std::cerr << "Accession " << acc << " belongs to taxon id " << pos->second  << std::endl;
					tax_id =  pos->second;
					uint64_t temp_id = tax_id;
					while(nodes->count(temp_id) > 0 && temp_id != 1) {
						if(include_ids.count(temp_id) > 0) {
							skip = false;
							break;
						}
						temp_id = nodes->at(temp_id);
					}
				}
				else {
					if(verbose) { std::cerr << "Accession " << acc <<" was not found in " << acc_taxid_filename << "\n"; }
				}
			}

			if(!skip) {
					if(!first) { output << "\n";  } else { first = false; }
					output << ">";
					if(addAcc) output << acc << "_";
					output << tax_id << "\n";
					outlinecount++;
				}
				else if(debug) { std::cerr << "Skipping sequence with header: " << line << std::endl; }
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
	out_file.close();

	std::cerr << getCurrentTime() << " Finished." << std::endl;
	return EXIT_SUCCESS;
}

void usage(char *progname) {
	print_usage_header();
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -m merged.dmp -g prot.accession2taxid \n", progname);
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file.\n");
	fprintf(stderr, "   -m FILENAME   Name of merged.dmp file.\n");
	fprintf(stderr, "   -g FILENAME   Name of prot.accession2taxid.FULL.gz file.\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -a            Prefix taxon ID with the accession number.\n");
	fprintf(stderr, "   -v            Verbose mode\n");
	fprintf(stderr, "   -d            Debug mode\n");
	exit(EXIT_FAILURE);
}

