/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <algorithm>
#include <string>
#include <list>
#include <utility>
#include <stdexcept>
#include <deque>

#include "util.hpp"

void usage(char *progname);

int main(int argc, char** argv) {

	std::unordered_map<uint64_t,uint64_t> nodes;
	std::unordered_map<uint64_t, std::string> node2name;
	std::unordered_map<uint64_t, std::string> node2rank;
	std::unordered_map<uint64_t, std::string> node2path;

	std::string nodes_filename = "";
	std::string names_filename = "";
	std::string in_filename = "";
	std::string out_filename;

	bool filter_unclassified = false;

	bool verbose = false;

	bool full_path = false;
	bool specified_ranks = false;
	std::string ranks_arg;
	std::list<std::string> ranks_list;
	std::set<std::string> ranks_set;


	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hvpur:n:t:i:o:")) != -1) {
		switch(c)  {
			case 'h':
				usage(argv[0]);
			case 'v':
				verbose = true; break;
			case 'u':
				filter_unclassified = true; break;
			case 'p':
				full_path = true; break;
			case 'o':
				out_filename = optarg; break;
			case 'n':
				names_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'i':
				in_filename = optarg; break;
			case 'r': {
				specified_ranks = true;
				ranks_arg = optarg; break; }
			default:
				usage(argv[0]);
		}
	}
	if(names_filename.length() == 0) { error("Please specify the location of the names.dmp file with the -n option."); usage(argv[0]); }
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file, using the -t option."); usage(argv[0]); }
	if(in_filename.length() == 0) { error("Please specify the location of the input file, using the -i option."); usage(argv[0]); }
	if(ranks_arg.length() > 0 && full_path) { error("Please use either option -r or -p, but not both of them."); usage(argv[0]); }

	/* parse user-supplied rank list into list and set */
	if(ranks_arg.length() > 0) {
		size_t begin = 0;
		size_t pos = -1;
		std::string rankname;
		while((pos = ranks_arg.find(",",pos+1)) != std::string::npos) {
			rankname = ranks_arg.substr(begin,(pos - begin));
			if(rankname.length()==0 || rankname==",") { begin=pos+1; continue; }
			ranks_list.emplace_back(rankname);
			ranks_set.emplace(rankname);
			begin = pos+1;
		}
		rankname = ranks_arg.substr(begin);
		if(!(rankname.length()==0 || rankname==",")) {
			ranks_set.emplace(rankname);
			ranks_list.emplace_back(rankname);
		}
	}

	/* read nodes.dmp */
	std::ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file.is_open()) { std::cerr << "Error: Could not open file " << nodes_filename << std::endl; usage(argv[0]); }
	if(verbose) std::cerr << "Reading taxonomic tree from file " << nodes_filename << std::endl;
	parseNodesDmpWithRank(nodes,node2rank,nodes_file);
	nodes_file.close();

	/* read names.dmp */
	std::ifstream names_file;
	names_file.open(names_filename);
	if(!names_file.is_open()) { std::cerr << "Error: Could not open file " << names_filename << std::endl; usage(argv[0]); }
	if(verbose) std::cerr << "Reading taxon names from file " << names_filename << std::endl;
	parseNamesDmp(node2name,names_file);
	names_file.close();

	std::ifstream in_file;
	in_file.open(in_filename);
	if(!in_file.is_open()) {  std::cerr << "Could not open file " << in_filename << std::endl; exit(EXIT_FAILURE); }

	std::ostream * out_stream;
	if(out_filename.length()>0) {
		if(verbose) std::cerr << "Output file: " << out_filename << std::endl;
		std::ofstream * output_file = new std::ofstream();
		output_file->open(out_filename);
		if(!output_file->is_open()) {  std::cerr << "Could not open file " << out_filename << " for writing" << std::endl; exit(EXIT_FAILURE); }
		out_stream = output_file;
	}
	else {
		out_stream = &std::cout;
	}

	if(verbose) std::cerr << "Processing " << in_filename <<"..." << "\n";

	std::string line;
	while(getline(in_file,line)) {
		if(line.length() == 0) { continue; }
		if(line[0] != 'C') {
			if(!filter_unclassified) {
				*out_stream << line << "\n";
			}
			continue;
		}

		size_t found = line.find('\t');
		found = line.find('\t',found+1);
		size_t end = line.find_first_not_of("0123456789",found+1);
		uint64_t taxonid;
		try {
			taxonid = stoul(line.substr(found,end-found));
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Found bad taxon id in line: " << line << std::endl;
			*out_stream << line << "\n";
			continue;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Error: Found bad taxon id (out of range error) in line: " << line << std::endl;
			*out_stream << line << "\n";
			continue;
		}

		if(nodes.count(taxonid)==0) {
			std::cerr << "Warning: Taxon ID " << taxonid << " in output file is not contained in taxonomic tree file "<< nodes_filename << ".\n";
			*out_stream << line << "\n";
			continue;
		}
		if(node2name.count(taxonid)==0) {
			std::cerr << "Warning: Taxon ID " << taxonid << " in output file is not found in file "<< names_filename << ".\n";
			*out_stream << line << "\n";
			continue;
		}
		if(full_path || specified_ranks) {
			if(node2path.count(taxonid)>0) { // look if path is already saved
				*out_stream << line << '\t' << node2path.at(taxonid) << "\n";
				continue;
			}
			std::deque<std::string> lineage; // for full_path
			std::map<std::string,std::string> curr_rank_values;
			if(specified_ranks) { //set the values for all specified ranks to NA, which will be overwritten by the actual values if they are found
				for(auto it : ranks_list) {
					curr_rank_values.emplace(it,"NA");
				}
			}
			//  go from leaf to root starting at taxonid and gather values for ranks
			uint64_t id = taxonid;
			while(nodes.count(id)>0 && id != nodes.at(id)) {
				std::string taxon_name;
				if(specified_ranks) {
					if(node2rank.count(id)==0 || node2rank.at(id)=="no rank") {  // no rank name
						id = nodes.at(id);
						continue;
					}
					std::string rank_name = node2rank.at(id);
					if(ranks_set.count(rank_name)==0) { // rank name is not in specified list of ranks
						id = nodes.at(id);
						continue;
					}
					taxon_name = getTaxonNameFromId(node2name, id, names_filename);
					curr_rank_values[rank_name] = taxon_name;
				}
				else { //full path
					taxon_name = getTaxonNameFromId(node2name, id, names_filename);
					lineage.emplace_front(taxon_name);
				}
				id = nodes.at(id);
			} // end while

			// assemble lineage into one string
			std::string lineage_text;
			if(specified_ranks) {
				for(auto it : ranks_list) {
					lineage_text += curr_rank_values[it];
					lineage_text += "; ";
				}
			}
			else { // full path
				for(auto  itl : lineage) {
					lineage_text += itl;
					lineage_text += "; ";
				}
			}
			// now lineage_text contains the final lineage for the taxon
			node2path.emplace(taxonid,lineage_text);
			*out_stream << line << '\t' << lineage_text << "\n";
		}
		else {
			*out_stream << line << '\t' << getTaxonNameFromId(node2name, taxonid, names_filename) << "\n";
		}
	}  // end while getline

	if(in_file.is_open()) {
		in_file.close();
	}
	out_stream->flush();
	if(out_filename.length()>0) {
		((std::ofstream*)out_stream)->close();
		delete ((std::ofstream*)out_stream);
	}

	return 0;

}

void usage(char *progname) {
	print_usage_header();
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju-names.out\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of input file\n");
	fprintf(stderr, "   -o FILENAME   Name of output file. If not specified, output will be printed to STDOUT.\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file\n");
	fprintf(stderr, "   -n FILENAME   Name of names.dmp file.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -u            Unclassified reads are not contained in the output.\n");
	fprintf(stderr, "   -p            Print full taxon path.\n");
	fprintf(stderr, "   -r            Print taxon path containing only ranks specified by a comma-separated list,\n");
	fprintf(stderr, "                 for example: superkingdom,phylum,class,order,family,genus,species\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

