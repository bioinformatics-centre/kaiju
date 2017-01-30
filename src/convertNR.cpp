/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <vector>
#include <unordered_set>
#include <iterator>

#include "Config.hpp"
#include "version.hpp"
#include "util.hpp"

void usage(char *progname);
bool has_parent(uint64_t id, uint64_t parent, std::unordered_map<uint64_t,uint64_t> & nodes);

int main(int argc, char **argv) {

	std::unordered_map<uint64_t,uint64_t> nodes;
	std::unordered_map<uint64_t,unsigned int> node2depth;
	std::unordered_map<std::string,uint64_t> acc2taxid;

	std::string nodes_filename;
	std::string list_filename;
	std::string acc_taxid_filename;
	std::string nr_filename;
	std::string out_filename;

	//uint64_t include_ids [3] = { 2, 2157, 10239 };  //Bacteria, Archaea, Viruses
	std::unordered_set<uint64_t> include_ids;
	bool verbose = false;
	bool debug = false;
	bool addcount = false;
	Config * config = new Config();

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hcdvrl:g:t:i:o:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'c':
				addcount = true; break;
			case 'l':
				list_filename = optarg; break;
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
	if(acc_taxid_filename.length() == 0) { error("Please specify the location of the prot.accession2taxid file, using the -g option."); usage(argv[0]); }
	if(out_filename.length() == 0) { error("Please specify the name of the output file, using the -o option."); usage(argv[0]); }

	config->nodes = &nodes;
	config->debug = debug;
	config->verbose = verbose;

	std::ifstream mapfile;
	mapfile.open(nodes_filename);
	if(!mapfile.is_open()) { std::cerr << "Error: Could not open file " << nodes_filename << std::endl; usage(argv[0]); }
	std::cerr << "Reading taxonomic tree from file " << nodes_filename << std::endl;
	std::string line;
	while(getline(mapfile, line)) {
		if(line.length() == 0) { continue; }
		try {
			size_t end = line.find_first_not_of("0123456789");
			uint64_t node = stoul(line.substr(0,end));
			size_t start = line.find_first_of("0123456789",end);
			end = line.find_first_not_of("0123456789",start+1);
			uint64_t parent = stoul(line.substr(start,end-start));
			nodes.emplace(node,parent);
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Found bad number in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
		}
	}
	mapfile.close();

	if(list_filename.length()==0) {
		std::cerr << "No taxa list specified, using Archaea, Bacteria, and Viruses." << std::endl;
		include_ids.insert((uint64_t)2);
		include_ids.insert((uint64_t)2157);
		include_ids.insert((uint64_t)10239);
	}
	else {
		std::ifstream list_file;
		list_file.open(list_filename);
		if(!list_file.is_open()) { std::cerr << "Error: Could not open file " << list_filename << std::endl; usage(argv[0]); }
		std::cerr << "Reading taxa from file " << list_filename << std::endl;
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
				if(nodes.count(taxid) > 0) {
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



	std::ifstream acc_taxid_file;
	acc_taxid_file.open(acc_taxid_filename);
	if(!acc_taxid_file.is_open()) { std::cerr << "Error: Could not open file " << acc_taxid_filename << std::endl; usage(argv[0]); }
	std::cerr << "Reading accession to taxon id map from file " << acc_taxid_filename << std::endl;
	getline(acc_taxid_file, line); // skip header line
	while(getline(acc_taxid_file, line)) {
		if(line.length() == 0) { continue; }
		try {
			size_t start = line.find('\t',0);
			size_t end = line.find("\t",start+1);
			std::string acc = line.substr(start+1,end-start-1);
			start = end+1;
			end = line.find_first_of("\t\n",start);
			uint64_t taxid = stoul(line.substr(start,end-start));
//			std::cerr << acc << "\t" << taxid << "\n";
			acc2taxid.emplace(acc,taxid);
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Found bad identifier in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
		}
	}
	acc_taxid_file.close();

	std::ifstream inputfile;
	if(nr_filename.length()>0) {
		inputfile.open(nr_filename);
		if(!inputfile.is_open()) { std::cerr << "Error: Could not open file " << nr_filename << std::endl; usage(argv[0]); }
	}

	if(verbose) std::cerr << "Writing to file " << out_filename << std::endl;
	std::ofstream out_file;
	out_file.open(out_filename);
	if(!out_file.is_open()) {  std::cerr << "Error: Could not open file " << out_filename << " for writing!" << std::endl; usage(argv[0]); }

	std::cerr << "Processing NR file " << nr_filename << std::endl;

	bool skip = true;
	bool first = true;
	std::ostringstream output;
	uint64_t fa_counter = 1;
	uint64_t outlinecount = 0;
	std::set<uint64_t> ids;
	while(getline(inputfile.is_open() ? inputfile : std::cin, line)){
		if(line.length() == 0) { continue; }
		if(line[0]=='>') {
			ids.clear();
			if(debug) std::cerr << "processing line " << line << std::endl;
			skip = true;
			size_t start = 1, end = 0;
			while((end = line.find(' ',start)) != std::string::npos) {
				// acc is between start and end
				std::string acc = line.substr(start, end - start);
				if(acc2taxid.count(acc)>0 && acc2taxid.at(acc)>0 && nodes.count(acc2taxid.at(acc))>0) {
					if(debug) std::cerr << "Accession " << acc << " belongs to taxon id " << acc2taxid.at(acc)  << std::endl;
					ids.insert(acc2taxid.at(acc));
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

			if(ids.size()>0) {
				bool keep = false;
				uint64_t lca = (ids.size()==1) ?  *(ids.begin()) : config->lca_from_ids(node2depth, ids);
				if(debug) std::cerr << "LCA=" << lca << std::endl;
				if(nodes.count(lca)==0) { std::cerr << "Taxon ID " << lca << " not found in taxonomy!" << std::endl; continue; }
				uint64_t id = lca;
				while(nodes.count(id)>0 && id != 1) {
					if(include_ids.count(id) > 0) {
						keep = true;
						break;
					}
					id = nodes.at(id);
				}
				if(keep) {
					if(!first) { output << "\n";  } else { first = false; }
					output << ">";
					if(addcount) output << fa_counter++ << "_";
					output << lca << "\n"; //endl;
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

	return EXIT_SUCCESS;
}

bool has_parent(uint64_t id, uint64_t parent, std::unordered_map<uint64_t,uint64_t> & nodes) {

	if(nodes.count(id)==0) { std::cerr << "Taxon ID " << id << " not found in taxonomy!" << std::endl; return false; }
	if(nodes.count(parent)==0) { std::cerr << "Taxon ID " << parent << " not found in taxonomy!" << std::endl; return false; }
	while(nodes.count(id)>0 && id != 1) {
		if(id == parent) { return true; }
		id = nodes.at(id);
	}
	return false;
}


void usage(char *progname) {
	fprintf(stderr, "Kaiju %s\n",KAIJUVERSION);
	fprintf(stderr, "Copyright 2015,2016 Peter Menzel, Anders Krogh\n");
	fprintf(stderr, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -g prot.accession2taxid -i nr\n", progname);
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file.\n");
	fprintf(stderr, "   -g FILENAME   Name of prot.accession2taxid file.\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of NR file. If this option is not used, then the program will read from STDIN.\n");
	fprintf(stderr, "   -l FILENAME   Name of file with taxon IDs. These IDs must be contained in nodes.dmp and denote the extracted clades from the NR file.\n");
	exit(EXIT_FAILURE);
}

