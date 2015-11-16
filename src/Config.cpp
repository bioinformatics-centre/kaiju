/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include "Config.hpp"

using namespace std;

Config::Config() { // constructor
}

Config::~Config() { // Destructor
		free(trans);
}

void Config::init() {
	trans = (uchar *)calloc(128,sizeof(uchar));
	for(int i=0; i<bwt->alen;++i) trans[(int)bwt->alphabet[i]]=i;
}


uint64_t lca_from_ids(unordered_map<uint64_t,unsigned int> & node2depth, set<uint64_t> & ids, unordered_map<uint64_t,uint64_t> * nodes) {

	if(ids.size() == 1) {
		return *(ids.begin());
	}
	uint num_ids = ids.size();
	uint64_t * leafs = (uint64_t *) calloc(num_ids,sizeof(uint64_t));
	unsigned int shallowest_depth = 100000;
	uint index = 0;
	for(auto it = ids.begin() ; it != ids.end(); ++it) {
		uint64_t id = *it;	

		if(nodes->count(id)==0) {
			cerr << "Warning: Taxon ID " << id << " in database is not contained in taxonomic tree.\n"; 
			num_ids--;
			continue;
		}

		// check if this id was already seen, then skip it
		leafs[index++] = id;

		//if id is alrady in the depth map then do not add it.
		if(node2depth.count(id)==0) {
			unsigned int depth = 1;
			while(nodes->count(id)>0 && id != nodes->at(id)) {
				depth++;
				id = nodes->at(id);	
			}
			node2depth.insert(pair<uint64_t,unsigned int>(*it,depth));
			//cerr << "Inserting to depth map: " << *it <<" -> " << depth << endl;
			if(depth < shallowest_depth) { shallowest_depth = depth; }
		}
		else if(node2depth.at(*it) < shallowest_depth) { shallowest_depth = node2depth.at(*it); }
	}

	if(num_ids<=0) {
		free(leafs);
		return 0;
	}

	//cerr << "shallowest depth = " << shallowest_depth << endl;

	// bring all IDs up to the same depth
	/*for (auto it=leafs.begin(); it != leafs.end(); ++it) {
		//cerr << "Bringing leaf " << *it << " to depth " << shallowest_depth <<" from depth " << node2depth->at(*it) << endl;
		for(int i = node2depth.at(*it) - shallowest_depth; i > 0; i--) {
			*it = nodes->at(*it);
		}
	}*/
	for(uint index = 0; index < num_ids; ++index) {
		for(int i = node2depth.at(leafs[index]) - shallowest_depth; i > 0; i--) {
			leafs[index]	= nodes->at(leafs[index]);
		}
	}

	while(true) {
		//foreach element in the list, check if id is the same, otherwise go one level up in tree, i.e. one more iteration 
		uint64_t first = leafs[0];
		bool found = true;
		//for (auto it=leafs.begin(); it != leafs.end(); ++it) {
		for(uint index = 0; index < num_ids; ++index) {
			if(first != leafs[index]) {
				found = false;
			}
			leafs[index] = nodes->at(leafs[index]);
		}
		if(found) {
			free(leafs);
			return first;
		}
	}
	free(leafs);

}

