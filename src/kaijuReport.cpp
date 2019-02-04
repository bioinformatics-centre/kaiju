/* This file is part of Kaiju, Copyright 2015-2019 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <set>
#include <list>
#include <deque>
#include <algorithm>
#include <string>
#include <functional>
#include <utility>
#include <stdexcept>
#define __STDC_FORMAT_MACROS
#include <inttypes.h>

#include "util.hpp"

void usage(char *progname);

int main(int argc, char** argv) {


	std::unordered_map<uint64_t,uint64_t> nodes;
	std::unordered_map<uint64_t, std::string> node2name;
	std::unordered_map<uint64_t, std::string> node2rank;

	uint64_t taxonid_viruses = 10239;

	std::string nodes_filename = "";
	std::string names_filename = "";
	std::string in_filename = "";
	std::string out_filename;
	float min_percent = 0.0;
	int min_read_count = 0;

	bool filter_unclassified = false;
	bool full_path = false;
	bool verbose = false;
	std::string rank;

	bool specified_ranks = false;
	std::string ranks_arg;
	std::list<std::string> ranks_list;
	std::set<std::string> ranks_set;

	// ------------------------------------------- START -------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt (argc, argv, "hpvur:n:t:i:o:m:c:l:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'v':
				verbose = true; break;
			case 'u':
				filter_unclassified = true; break;
			case 'r':
				rank = optarg; break;
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
			case 'c': {
									try {
										min_read_count = std::stoi(optarg);
									}
									catch(const std::invalid_argument& ia) {
										std::cerr << "Invalid number in -c " << optarg << std::endl;
									}
									catch (const std::out_of_range& oor) {
										std::cerr << "Invalid number in -c " << optarg << std::endl;
									}
									break;
								}
			case 'm': {
									try {
										min_percent = std::stof(optarg);
									}
									catch(const std::invalid_argument& ia) {
										std::cerr << "Invalid number in -m " << optarg << std::endl;
									}
									catch (const std::out_of_range& oor) {
										std::cerr << "Invalid number in -m " << optarg << std::endl;
									}
									break;
								}
			case 'l': {
				specified_ranks = true;
				ranks_arg = optarg; break; }
			default:
				usage(argv[0]);
		}
	}
	if(names_filename.length() == 0) { error("Please specify the location of the names.dmp file with the -n option."); usage(argv[0]); }
	if(nodes_filename.length() == 0) { error("Please specify the location of the nodes.dmp file with the -t option."); usage(argv[0]); }
	if(out_filename.length() == 0) { error("Please specify the name of the output file with the -o option."); usage(argv[0]); }
	if(in_filename.length() == 0) { error("Please specify the location of the input file with the -i option."); usage(argv[0]); }
	if(rank.length() == 0) { error("Please specify the rank (phylum, class, order, family, genus, or species) with the -r option."); usage(argv[0]); }
	if(ranks_arg.length() > 0 && full_path) { error("Please use either option -r or -l, but not both of them."); usage(argv[0]); }
	if(!(rank.compare("phylum")==0 || rank.compare("class")==0 || rank.compare("order")==0 || rank.compare("family")==0 || rank.compare("genus")==0 || rank.compare("species")==0)) {
		error("Rank must be one of: phylum, class, order, family, genus, species."); usage(argv[0]);
	}
	if(min_read_count < 0) {
		error("Min required read count (-c) must be >= 0"); usage(argv[0]);
	}
	if(min_percent < 0.0 || min_percent > 100.0) {
		error("Min required percent (-m) must be between 0.0 and 100.0"); usage(argv[0]);
	}
	if(min_percent > 0.0 && min_read_count > 0) {
		error("Either specify minimum percent with -m or minimum read count with -c."); usage(argv[0]);
	}

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

		if(ranks_set.count(rank)==0) {
			error("Specified rank " + rank + " is not contained in rank list supplied with option -l"); usage(argv[0]);
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

	if(verbose) std::cerr << "Processing " << in_filename <<"..." << "\n";

	std::map<uint64_t, uint64_t> node2hitcount;
	uint64_t unclassified = 0;
	uint64_t totalreads = 0;
	uint64_t total_virus_reads = 0;
	std::string line;
	while(std::getline(in_file,line)) {

		if(line.length() == 0) { continue; }
		totalreads++;
		if(line[0] != 'C') { unclassified++; continue; }

		size_t found = line.find('\t');
		found = line.find('\t',found+1);
		size_t end = line.find_first_not_of("0123456789",found+1);
		try {
			uint64_t taxonid = stoul(line.substr(found,end-found));
			if(nodes.count(taxonid)==0) {
				std::cerr << "Warning: Taxon ID " << taxonid << " in output file is not contained in taxonomic tree file "<< nodes_filename << ".\n";
				continue;
			}
			if(is_ancestor(nodes,taxonid_viruses,taxonid)) {
				total_virus_reads++;
			}
			if(node2hitcount.count(taxonid)>0)
				node2hitcount[taxonid]++;
			else
				node2hitcount[taxonid] = 1;
			}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Found bad taxon id in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Error: Found bad taxon id (out of range error) in line: " << line << std::endl;
		}
	}
	if(in_file.is_open()) in_file.close();

	std::map<uint64_t, uint64_t> node2summarizedhits;

	for(auto it : node2hitcount) {
		uint64_t id = it.first;
		uint64_t reads = it.second;
		while(nodes.count(id)>0 && id != nodes.at(id)) {
			(node2summarizedhits.count(id) > 0) ?  node2summarizedhits[id] += reads : node2summarizedhits[id]  = reads;
			id = nodes.at(id);
		}
	}

	if(filter_unclassified)
		totalreads -= unclassified;


	// Go through node2summarizedhits and check each node at the specified rank
	// if it is above threshold, then add it to a sorted map for later printing
	// Viruses are omitted here.
	uint64_t reads_at_rank_sum = 0;
	uint64_t reads_at_rank_below_percent_threshold = 0;
	uint64_t reads_at_rank_below_count_threshold = 0;

	std::multimap<uint64_t,uint64_t ,std::greater<uint64_t>> sorted_count2ids;
	for(auto it : node2summarizedhits) {
		uint64_t id = it.first;
		uint64_t count = it.second;
		if(is_ancestor(nodes,taxonid_viruses,id)) {
			continue;
		}
		if(node2rank.count(id)==0) { std::cerr << "Error: No rank specified for taxonid " << id << std::endl; continue; }
		if(rank == node2rank[id]) {
			if((int)count >= min_read_count) {
				float percent = (float)it.second/(float)totalreads*100;
				if(percent >= min_percent)
					sorted_count2ids.emplace(count,id);
				else
					reads_at_rank_below_percent_threshold += count;
			} else {
				reads_at_rank_below_count_threshold += count;
			}
			reads_at_rank_sum += count;
		}
	}

	if(filter_unclassified)
		assert(totalreads >= reads_at_rank_sum);
	else
		assert(totalreads >= unclassified + reads_at_rank_sum);

	uint64_t viruses = (node2summarizedhits.count(taxonid_viruses) > 0) ?  node2summarizedhits[taxonid_viruses] : 0;
	assert(total_virus_reads==viruses);

	uint64_t above = (filter_unclassified) ? totalreads - reads_at_rank_sum  : totalreads - unclassified - reads_at_rank_sum;
	above -= viruses;

	/* ---------- print output ---------------- */
	if(verbose) std::cerr << "Writing to file " << out_filename << std::endl;
	FILE * report_file = fopen(out_filename.c_str(),"w");
	if(report_file==NULL) {  std::cerr << "Could not open file " << out_filename << " for writing" << std::endl; exit(EXIT_FAILURE); }
	fprintf(report_file,"        %%\t    reads\t%s\n",rank.c_str());
	fprintf(report_file,"-------------------------------------------\n");

	for(auto it : sorted_count2ids) {
		std::string name;
		if(full_path || specified_ranks) {
			uint64_t id = it.second;
			std::deque<std::string> lineage; // for full_path
			std::map<std::string,std::string> curr_rank_values;
			if(specified_ranks) { //set the values for all specified ranks to NA, which will be overwritten by the actual values if they are found
				for(auto it : ranks_list) {
					curr_rank_values.emplace(it,"NA");
				}
			}
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
					lineage_text += curr_rank_values[it] + ";";
				}
			}
			else { // full path
				for(auto  itl : lineage) {
					lineage_text += itl + "; ";
				}
			}
			// now lineage_text contains the final lineage for the taxon
			name = lineage_text;
		}
		else {
			name = getTaxonNameFromId(node2name, it.second, names_filename);
		}
		float percent = (float)it.first/(float)totalreads*100.0f;
		fprintf(report_file,"%9.6f\t%9" PRIu64 "\t%s\n", percent, it.first, name.c_str() );
	}

	fprintf(report_file,"-------------------------------------------\n");
	fprintf(report_file,"%9.6f\t%9" PRIu64 "\tViruses\n", (float)viruses/(float)totalreads*100.0, viruses);
	fprintf(report_file,"%9.6f\t%9" PRIu64 "\tcannot be assigned to a %s \n", (float)above/(float)totalreads*100.0, above, rank.c_str());
	if(min_read_count > 0)
		fprintf(report_file,"%9.6f\t%9" PRIu64 "\tbelong to a %s having less than %i reads\n", (float)reads_at_rank_below_count_threshold/(float)totalreads*100.0, reads_at_rank_below_count_threshold, rank.c_str(), min_read_count);
	if(min_percent > 0.0)
		fprintf(report_file,"%9.6f\t%9" PRIu64 "\tbelong to a %s with less than %g%% of all reads\n", (float)reads_at_rank_below_percent_threshold/(float)totalreads*100.0, reads_at_rank_below_percent_threshold, rank.c_str(), min_percent);
	fprintf(report_file,"-------------------------------------------\n");
	if(filter_unclassified)
		fprintf(report_file,"%9.6f\t%9" PRIu64 "\tunclassified\n", (float)unclassified/(float)(totalreads+unclassified)*100.0, unclassified );
	else
		fprintf(report_file,"%9.6f\t%9" PRIu64 "\tunclassified\n", (float)unclassified/(float)(totalreads)*100.0, unclassified );

	fclose(report_file);

	return EXIT_SUCCESS;
}


void usage(char *progname) {
	print_usage_header();
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju.report\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of input file\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file\n");
	fprintf(stderr, "   -n FILENAME   Name of names.dmp file.\n");
	fprintf(stderr, "   -r STRING     Taxonomic rank, must be one of: phylum, class, order, family, genus, species\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -m FLOAT      Number in [0, 100], denoting the minimum required percentage for the taxon to be reported (default: 0.0)\n");
	fprintf(stderr, "   -c INT        Integer number > 0, denoting the minimum required number of reads for the taxon to be reported (default: 0)\n");
	fprintf(stderr, "   -u            Unclassified reads are not counted for the total reads when calculating percentages for classified reads.\n");
	fprintf(stderr, "   -p            Print full taxon path.\n");
	fprintf(stderr, "   -l            Print taxon path containing only ranks specified by a comma-separated list,\n");
	fprintf(stderr, "                 for example: superkingdom,phylum,class,order,family,genus,species\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Only one of the options -m and -c may be used at a time.\n");
	exit(EXIT_FAILURE);
}

