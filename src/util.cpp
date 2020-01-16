/* This file is part of Kaiju, Copyright 2015-2017 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include "util.hpp"

extern "C" {
#include "./bwt/bwt.h"
}

void print_usage_header() {
	fprintf(stderr, "Kaiju %s\n",KAIJUVERSION);
	fprintf(stderr, "Copyright 2015-2020 Peter Menzel, Anders Krogh\n");
	fprintf(stderr, "License GPLv3+: GNU GPL version 3 or later <http://gnu.org/licenses/gpl.html>\n");
	fprintf(stderr, "\n");
}

void error(const std::string e) {
	std::cerr << "Error: " << e << std::endl << std::endl;
}

inline bool isalpha(const char & c) {
	return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

void strip(std::string & s) {
		for(auto it = s.begin(); it!=s.end(); ++it) {
			if(!isalpha(*it)) {
				s.erase(it);
				it--;
			}
		}
}

std::string getCurrentTime() {
	time_t t = time(0);
	char buffer[9] = {0};
	strftime(buffer, 9, "%H:%M:%S", localtime(&t));
	return std::string(buffer);
}

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(const std::unordered_map<uint64_t,uint64_t> & nodes, const std::string & id1, const std::string & id2) {

		uint64_t node1;
		uint64_t node2;
		try {
			node1 = stoul(id1);
			node2 = stoul(id2);
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Error: Bad number in taxon id" << std::endl;
			return false;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Error: Bad number (out of range error) in taxon id" << std::endl;
			return false;
		}

		return is_ancestor(nodes,node1,node2);
}

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(const std::unordered_map<uint64_t,uint64_t> & nodes, uint64_t node1, uint64_t node2) {
		if(nodes.count(node1)==0) { std::cerr << "Taxon ID " << node1 << " not found in taxonomy!" << std::endl; return false; }
		if(nodes.count(node2)==0) { std::cerr << "Taxon ID " << node2 << " not found in taxonomy!" << std::endl; return false; }
		if(node2==node1) {
			return true;
		}
		/* climb up from node 2 and return true if encountering node 1 */
		while(nodes.count(node2)>0 && node2 != nodes.at(node2)) {
			node2 = nodes.at(node2);
			if(node2==node1) {
				return true;
			}
		}
		return false;
}

void parseNodesDmp(std::unordered_map<uint64_t,uint64_t> & nodes, std::ifstream & nodes_file) {
		nodes.reserve(2e6);
		std::string line;
		while(std::getline(nodes_file, line)) {
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
			catch(const std::out_of_range& oor) {
				std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
			}
		}
}

void parseNodesDmpWithRank(std::unordered_map<uint64_t,uint64_t> & nodes, std::unordered_map<uint64_t,std::string> & node2rank, std::ifstream & nodes_file) {
	nodes.reserve(2e6);
	std::string line;
	while(std::getline(nodes_file, line)) {
		if(line.length() == 0) { continue; }
		try {
			size_t end = line.find_first_not_of("0123456789");
			//cerr << "end=" << end << "\t";
			uint64_t node = stoul(line.substr(0,end));
			size_t start = line.find_first_of("0123456789",end);
			//cerr << "start=" << start <<"\t";
			end = line.find_first_not_of("0123456789",start+1);
			//cerr << "end=" << end <<"\t";
			uint64_t parent = stoul(line.substr(start,end-start));
			start = line.find_first_of("abcdefghijklmnopqrstuvwxyz",end);
			//cerr << "start=" << start <<","<< line[start] << "\t";
			end = line.find_first_not_of("abcdefghijklmnopqrstuvwxyz ",start);
			//cerr << "end=" << end << "\t";
			std::string rank = line.substr(start,end-start);
			nodes.emplace(node,parent);
			node2rank.emplace(node,rank);
			//cerr << node << "\t" << parent << "\t" <<rank << "\n";

		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Found bad number in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
		}
	}
}

void parseNamesDmp(std::unordered_map<uint64_t,std::string> & names, std::ifstream & names_file) {
	std::string line;
	while(std::getline(names_file, line)) {
		if(line.length() == 0) { continue; }
		try {
			if(line.find("scientific name")==std::string::npos) continue;
			size_t start = line.find_first_of("0123456789");
			size_t end = line.find_first_not_of("0123456789",start);
			uint64_t node_id = stoul(line.substr(start,end-start));
			start = line.find_first_not_of("\t|",end);
			end = line.find_first_of("\t|",start+1);
			std::string name = line.substr(start,end-start);
			names.emplace(node_id,name);
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Found bad number in line: " << line << std::endl;
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Found bad number (out of range error) in line: " << line << std::endl;
		}
	}


}


std::string getTaxonNameFromId(const std::unordered_map<uint64_t,std::string> & node2name, uint64_t id, const std::string & names_filename) {
	std::string taxon_name;
	if(node2name.count(id)==0) {
		std::cerr << "Warning: Taxon ID " << id << " is not found in file "<< names_filename << "." << std::endl;
		taxon_name = "taxonid:"; taxon_name += std::to_string(id);
	}
	else {
		taxon_name = node2name.at(id);
	}
	return taxon_name;
}

uint64_t lca_from_ids(Config * config, std::unordered_map<uint64_t,unsigned int> & node2depth, const std::set<uint64_t> & ids) {

	size_t num_ids = ids.size();
	if(num_ids == 1) {
		return *(ids.begin());
	}
	uint64_t * leafs = (uint64_t *) calloc(num_ids,sizeof(uint64_t));
	unsigned int shallowest_depth = 100000;
	unsigned int index = 0;
	for(auto it : ids) {

		if(config->nodes->count(it)==0) {
			if(config->verbose) std::cerr << "Warning: Taxon ID " << it << " in database is not contained in taxonomic tree.\n";
			num_ids--;
			continue;
		}

		// check if this id was already seen, then skip it
		leafs[index++] = it;

		//if id is alrady in the depth map then do not add it.
		auto pos = node2depth.find(it);
		if(pos == node2depth.end()) {
			unsigned int depth = 1;
			uint64_t id = it;
			while(config->nodes->count(id)>0 && id != config->nodes->at(id)) {
				depth++;
				id = config->nodes->at(id);
			}
			node2depth.emplace(it,depth);
			//cerr << "Inserting to depth map: " << *it <<" -> " << depth << endl;
			if(depth < shallowest_depth) { shallowest_depth = depth; }
		}
		else if(pos->second < shallowest_depth) {
			shallowest_depth = pos->second;
		}
	}

	if(num_ids<=0) {
		free(leafs);
		return 0;
	}

	//cerr << "shallowest depth = " << shallowest_depth << endl;

	for(int index = 0; index < num_ids; ++index) {
		for(int i = node2depth.at(leafs[index]) - shallowest_depth; i > 0; i--) {
			leafs[index] = config->nodes->at(leafs[index]);
		}
	}

	while(true) {
		//foreach element in the list, check if id is the same, otherwise go one level up in tree, i.e. one more iteration
		uint64_t first = leafs[0];
		bool found = true;
		//for (auto it=leafs.begin(); it != leafs.end(); ++it) {
		for(size_t index = 0; index < num_ids; ++index) {
			if(first != leafs[index]) {
				found = false;
			}
			leafs[index] = config->nodes->at(leafs[index]);
		}
		if(found) {
			free(leafs);
			return first;
		}
	}
	free(leafs);

}

void readFMI(std::string fmi_filename, Config * config) {

	if(config->verbose) std::cerr << " Reading index from file " << fmi_filename << std::endl;
	FILE * fp = fopen(fmi_filename.c_str(),"r");
	if(!fp) { error("Could not open file " + fmi_filename); exit(EXIT_FAILURE); }
	BWT * b = readIndexes(fp);
	fclose(fp);
	if(config->debug) fprintf(stderr,"BWT of length %ld has been read with %d sequences, alphabet=%s\n", b->len, b->nseq, b->alphabet);
	config->bwt = b;
	config->fmi = b->f;

}
