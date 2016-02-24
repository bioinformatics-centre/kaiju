/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <getopt.h>
#include <vector>
#include <iterator>

#include "Config.hpp"

void usage(char *progname);
bool has_parent(uint64_t id, uint64_t parent, std::unordered_map<uint64_t,uint64_t> & nodes);

using namespace std;

int main(int argc, char **argv) {

	unordered_map<uint64_t,uint64_t> nodes;
	unordered_map<uint64_t,unsigned int> node2depth;
	unordered_map<uint64_t,uint64_t> gi2taxid;

	string nodes_filename;
	string gi_taxid_filename;
	string nr_filename;
	string out_filename;

	uint64_t include_ids [3] = { 2, 2157, 10239 };  //Bacteria, Archaea, Viruses
	bool verbose = false;
	bool debug = false;
	bool addcount = false;
	Config * config = new Config();

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	char c;
	while ((c = getopt (argc, argv, "hcdvrg:t:i:o:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'd':
				debug = true; break;
			case 'v':
				verbose = true; break;
			case 'c':
				addcount = true; break;
			case 't':
				nodes_filename = optarg; break;
			case 'g':
				gi_taxid_filename = optarg; break;
			case 'i':
				nr_filename = optarg; break;
			case 'o':
				out_filename = optarg; break;
			default:
				usage(argv[0]);
		}
	}
	if(nodes_filename.length() == 0) { cerr << "Error: Please specify the location of the nodes.dmp file, using the -t option."  << endl; usage(argv[0]); }
	if(gi_taxid_filename.length() == 0) { cerr << "Error: Please specify the location of the gi_taxid_prot.dmp file, using the -g option."  << endl; usage(argv[0]); }
	if(out_filename.length() == 0) { cerr << "Error: Please specify the name of the output file, using the -o option."  << endl; usage(argv[0]); }

	config->nodes = &nodes;
	config->debug = debug;
	config->verbose = verbose;

	ifstream mapfile;
	mapfile.open(nodes_filename.c_str());
	if(!mapfile.is_open()) { cerr << "Error: Could not open file " << nodes_filename << endl; usage(argv[0]); }
	cerr << "1/3 Reading taxonomic tree from file " << nodes_filename << endl;

	string line;
	while(getline(mapfile, line)) {
		if(line.length() == 0) { continue; } 
		try {
			size_t end = line.find_first_not_of("0123456789");
			uint64_t node = stoul(line.substr(0,end));
			size_t start = line.find_first_of("0123456789",end);
			end = line.find_first_not_of("0123456789",start+1);
			uint64_t parent = stoul(line.substr(start,end-start));
			nodes.insert(make_pair(node,parent));  //maybe the nodes->at(node) = parent;  would be faster?!
		}
		catch(const std::invalid_argument& ia) {
			cerr << "Found bad number in line: " << line << endl; 
		}
		catch (const std::out_of_range& oor) {
			cerr << "Found bad number (out of range error) in line: " << line << endl; 
		}
	}
	mapfile.close();

	ifstream gi_taxid_file;
	gi_taxid_file.open(gi_taxid_filename.c_str());
	if(!gi_taxid_file.is_open()) { cerr << "Error: Could not open file " << gi_taxid_filename << endl; usage(argv[0]); }
	cerr << "2/3 Reading gi to taxon id map from file " << gi_taxid_filename << endl;
	while(getline(gi_taxid_file, line)) {
		if(line.length() == 0) { continue; }
		size_t end = line.find_first_not_of("0123456789");
		try {
			uint64_t gi = stoul(line.substr(0,end));
			size_t start = line.find_first_of("0123456789",end);
			end = line.find_first_not_of("0123456789",start+1);
			uint64_t taxid = stoul(line.substr(start,end-start));
			gi2taxid.insert(make_pair(gi,taxid));  //maybe the nodes->at(node) = parent;  would be faster?!
		}
		catch(const std::invalid_argument& ia) {
			cerr << "Found bad identifier in line: " << line << endl; 
		}
		catch (const std::out_of_range& oor) {
			cerr << "Found bad number (out of range error) in line: " << line << endl; 
		}
	}
	gi_taxid_file.close();

	ifstream inputfile;
	if(nr_filename.length()>0) {
		inputfile.open(nr_filename);
		if(!inputfile.is_open()) { cerr << "Error: Could not open file " << nr_filename << endl; usage(argv[0]); }
	}

	if(verbose) cerr << "Writing to file " << out_filename << endl;
	ofstream out_file;
	out_file.open(out_filename);    
	if(!out_file.is_open()) {  cerr << "Could not open file " << out_filename << " for writing" << endl; exit(EXIT_FAILURE); }

	cerr << "3/3 Processing NR file " << nr_filename << endl;

	bool skip = true;
	bool first = true;
	ostringstream output;
	uint64_t fa_counter = 1;
	uint64_t outlinecount = 0;
	set<uint64_t> ids;
	while(getline(inputfile.is_open() ? inputfile : cin, line)){
		if(line.length() == 0) { continue; }
		if(line[0]=='>') {
			ids.clear();
			if(debug) cerr << "processing line " << line << endl;
			skip = true;
			size_t start = 0, stop = 0;
			while((start = line.find("gi|",start)) != string::npos) {
				if((stop = line.find("|",start+3)) != string::npos && stop - start -3 > 0) {
					if(debug) cerr << line.substr(start+3,stop-start-3) << endl;
					if(line.substr(start+3,stop-start-3).find_first_not_of("0123456789")==string::npos) {
						if(debug) cerr << "Found id at " << start << "-" << stop << "=" << line.substr(start + 3,stop - start -3) << endl;
						try {
							uint64_t id = stoul(line.substr(start + 3,stop - start -3));
							if(id > 0 && gi2taxid.count(id)>0 && gi2taxid.at(id)>0 && nodes.count(gi2taxid.at(id))>0) {
								if(debug) cerr << "gi " << id << " belongs to taxon id " << gi2taxid.at(id)  << endl;
								ids.insert(gi2taxid.at(id));
							}
							else if(verbose) { cerr << "ID " << id <<" was either not found in gi_taxid_prot.dmg or in nodes.dmp." << endl; }
						}
						catch(const std::invalid_argument& ia) {
							if(verbose) { cerr << "Found bad GI at pos " << start << " in line " << line << endl; }
						}
						catch (const std::out_of_range& oor) {
							if(verbose) { std::cerr << "Found bad GI (out of range error) at position " << start << " in line " << line << endl; }
						}
					}
				}
				else { 
					if(verbose) { cerr << "Found bad GI at pos " << start << " in line " << line << endl; }
				}
				start=start+3;
			}
			if(ids.size()>0) { 
				uint64_t lca = (ids.size()==1) ?  *(ids.begin()) : config->lca_from_ids(node2depth, ids);
				bool keep = false;
				for(auto include_id : include_ids) {
					if(has_parent(lca,include_id,nodes)) {
						keep = true;
						break;
					}
				}
				if(keep) {
					if(!first) { output << "\n";  } else { first = false; }
					output << ">";
					if(addcount) output << fa_counter++ << "_";
					output << lca << "\n"; //endl;
					skip = false;
					outlinecount++;
				}
				else if(debug) { cerr << "Skipping sequence starting with taxon id " << *(ids.begin()) << endl; }
			} 
			else if(verbose) cerr << "Could not find taxonomy id for line " << line << endl;
		}
		else {
			if(!skip) {
				size_t p = 0;
				while((p = line.find_first_of("ARNDCQEGHILKMFPSTWYV",p)) != string::npos) {
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
	output << endl;
	out_file << output.str();
	if(inputfile.is_open())
		inputfile.close();
	out_file.close();

	return EXIT_SUCCESS;
}

bool has_parent(uint64_t id, uint64_t parent, unordered_map<uint64_t,uint64_t> & nodes) {

	if(nodes.count(id)==0) { cerr << "Taxon ID " << id << " not found in taxonomy!" << endl; return false; }
	if(nodes.count(parent)==0) { cerr << "Taxon ID " << parent << " not found in taxonomy!" << endl; return false; }
	while(nodes.count(id)>0 && id != 1) {
		if(id == parent) { return true; }
		id = nodes.at(id);
	}
	return false;
}


void usage(char *progname) { 
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -g gi_taxid_prot.dmp -i nr\n", progname);
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file.\n");
	fprintf(stderr, "   -g FILENAME   Name of gi_taxid_prot.dmp file.\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of NR file. If this option is not used, then the program will read from STDIN.\n");
	exit(EXIT_FAILURE);
}

