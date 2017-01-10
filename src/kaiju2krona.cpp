/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <string>
#include <vector>
#include <stdexcept>

#include "version.hpp"
#include "util.hpp"

void usage(char *progname);

int main(int argc, char** argv) {


	std::unordered_map<uint64_t,uint64_t> nodes;
	std::unordered_map<uint64_t, std::string> node2name;

	std::string nodes_filename = "";
	std::string names_filename = "";
	std::string in1_filename = "";
	std::string out_filename;

	bool verbose = false;
	bool count_unclassified = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "huvn:t:i:o:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'v':
				verbose = true; break;
			case 'u':
				count_unclassified = true; break;
			case 'o':
				out_filename = optarg; break;
			case 'n':
				names_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'i':
				in1_filename = optarg; break;
			default:
				usage(argv[0]);
		}
	}

	if(out_filename.length() == 0) { error("Error: Please specify the name of the output file, using the -o option."); usage(argv[0]); }
	if(names_filename.length() == 0) { error("Please specify the location of the names.dmp file with the -n option."); usage(argv[0]); }
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file, using the -t option."); usage(argv[0]); }
	if(in1_filename.length() == 0) { error("Please specify the location of the input file, using the -i option."); usage(argv[0]); }

	std::ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file.is_open()) { std::cerr << "Error: Could not open file " << nodes_filename << std::endl; usage(argv[0]); }
	if(verbose) std::cerr << "Reading taxonomic tree from file " << nodes_filename << std::endl;
	parseNodesDmp(nodes,nodes_file);
	nodes_file.close();

	std::ifstream names_file;
	names_file.open(names_filename);
	if(!names_file.is_open()) { std::cerr << "Error: Could not open file " << names_filename << std::endl; usage(argv[0]); }
	if(verbose) std::cerr << "Reading taxon names from file " << names_filename << std::endl;
	parseNamesDmp(node2name,names_file);
	names_file.close();

	if(verbose) std::cerr << "Processing " << in1_filename <<"..." << "\n";

	std::ifstream in1_file;
	in1_file.open(in1_filename);
	if(!in1_file.is_open()) {  std::cerr << "Could not open file " << in1_filename << std::endl; exit(EXIT_FAILURE); }

	std::unordered_map<uint64_t, uint64_t> node2hitcount;
	long num_unclassified = 0;

	std::string line;
	while(getline(in1_file,line)) {
		if(line.length() == 0) { continue; }
		if(line[0] != 'C') {
			if(count_unclassified) num_unclassified++;
			continue;
		}

		size_t found = line.find('\t');
		found = line.find('\t',found+1);
		size_t end = line.find_first_not_of("0123456789",found+1);
		try {
			uint64_t taxonid = stoul(line.substr(found,end-found));
			if(node2hitcount.count(taxonid)>0)
				node2hitcount[taxonid]++;
			else
				node2hitcount[taxonid] = 1;
			}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Found bad taxon id in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Found bad taxon id (out of range error) in line: " << line << std::endl;
		}

	} // end main loop around file1


	if(in1_file.is_open()) in1_file.close();


	if(verbose) std::cerr << "Writing to file " << out_filename << std::endl;
	std::ofstream krona_file;
	krona_file.open(out_filename);
	if(!krona_file.is_open()) {  std::cerr << "Could not open file " << out_filename << " for writing" << std::endl; exit(EXIT_FAILURE); }
	for(auto  it : node2hitcount) {
		uint64_t id = it.first;
		if(nodes.count(id)==0) {
			std::cerr << "Warning: Taxon ID " << id << " found in input file is not contained in taxonomic tree file "<< nodes_filename << ".\n";
			continue;
		}
		if(node2name.count(id)==0) {
			std::cerr << "Warning: Taxon ID " << id << " found in input file is not contained in names.dmp file "<< names_filename << ".\n";
			continue;
		}
		std::vector<std::string> lineage;
		lineage.push_back(node2name.at(id));
		while(nodes.count(id)>0 && id != nodes.at(id)) {
			if(node2name.count(nodes.at(id))==0) {
				std::cerr << "Warning: Taxon ID " << nodes.at(id) << " found in input file is not contained in names file "<< names_filename << ".\n";
			}
			else {
				lineage.insert(lineage.begin(),node2name.at(nodes.at(id)));
			}
			id = nodes.at(id);
		}
		krona_file << it.second ;
		for(auto  itl : lineage) krona_file << "\t" << itl;
		krona_file << "\n";
	}
	if(count_unclassified && num_unclassified>0) {
		krona_file << num_unclassified << "\tUnclassified" << std::endl;
	}
	krona_file.close();

	return EXIT_SUCCESS;
}

void usage(char *progname) {
	fprintf(stderr, "Kaiju %s\n",KAIJUVERSION);
	fprintf(stderr, "Copyright 2015,2016 Peter Menzel, Anders Krogh\n");
	fprintf(stderr, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju2krona.out\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of input file\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file\n");
	fprintf(stderr, "   -n FILENAME   Name of names.dmp file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	fprintf(stderr, "   -u            Include count for unclassified reads in output.\n");
	exit(EXIT_FAILURE);
}

