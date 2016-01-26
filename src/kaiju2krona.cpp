/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <locale>
#include <string>
#include <functional>


void usage(char *progname);

using namespace std;

int main(int argc, char** argv) {


	unordered_map<uint64_t,uint64_t> nodes;
	unordered_map<uint64_t, string> node2name;

	string nodes_filename = "";
	string names_filename = "";
	string in1_filename = "";
	string out_filename;

	bool verbose = false;
	bool count_unclassified = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	char c;
	while ((c = getopt (argc, argv, "huvn:t:i:o:")) != -1) {
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
	if(names_filename.length() == 0) { cerr << "Error: Please specify the location of the names.dmp file, using the -n option."  << endl; usage(argv[0]); }
	if(nodes_filename.length() == 0) { cerr << "Error: Please specify the location of the nodes.dmp file, using the -t option."  << endl; usage(argv[0]); }
	if(out_filename.length() == 0) { cerr << "Error: Please specify the name of the output file, using the -o option."  << endl; usage(argv[0]); }
	if(in1_filename.length() == 0) { cerr << "Error: Please specify the location of the input file, using the -i option."  << endl; usage(argv[0]); }
	

	ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file.is_open()) { cerr << "Error: Could not open file " << nodes_filename << endl; usage(argv[0]); }
	if(verbose) cerr << "Reading taxonomic tree from file " << nodes_filename << endl;
	string line;
	while(getline(nodes_file, line)) {
		if(line.length() == 0) { continue; }
		try {
			size_t end = line.find_first_not_of("0123456789");
			uint64_t node = stoul(line.substr(0,end));
			size_t start = line.find_first_of("0123456789",end);
			end = line.find_first_not_of("0123456789",start+1);
			uint64_t parent = stoul(line.substr(start,end-start));
			nodes.insert(make_pair(node,parent));
		}
		catch(const std::invalid_argument& ia) {
			cerr << "Found bad number in line: " << line << endl; 
		}
		catch (const std::out_of_range& oor) {
			cerr << "Found bad number (out of range error) in line: " << line << endl; 
		}
	}
	nodes_file.close();


	ifstream names_file;
	names_file.open(names_filename);
	if(!names_file.is_open()) { cerr << "Error: Could not open file " << names_filename << endl; usage(argv[0]); }
	if(verbose) cerr << "Reading taxon names from file " << names_filename << endl;
	while(getline(names_file, line)) {
		if(line.length() == 0) { continue; }
		try {
			if(line.find("scientific name")==string::npos) continue;			
			size_t start = line.find_first_of("0123456789");
			size_t end = line.find_first_not_of("0123456789",start);
			uint64_t node = stoul(line.substr(start,end-start));
			start = line.find_first_not_of("\t|",end);
			end = line.find_first_of("\t|",start+1);
			string name = line.substr(start,end-start);
			node2name.insert(make_pair(node,name)); 
		}
		catch(const std::invalid_argument& ia) {
			cerr << "Found bad number in line: " << line << endl; 
		}
		catch (const std::out_of_range& oor) {
			cerr << "Found bad number (out of range error) in line: " << line << endl; 
		}
	}
	names_file.close();

	if(verbose) cerr << "Processing " << in1_filename <<"..." << "\n";

	ifstream in1_file;
	in1_file.open(in1_filename);    
	
	if(!in1_file.is_open()) {  cerr << "Could not open file " << in1_filename << endl; exit(EXIT_FAILURE); }
	
	map<uint64_t, uint64_t> node2hitcount;
	long num_unclassified = 0;
	
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
			cerr << "Found bad taxon id in line: " << line << endl; 
		}
		catch (const std::out_of_range& oor) {
			cerr << "Found bad taxon id (out of range error) in line: " << line << endl; 
		}
		

	} // end main loop around file1


	if(in1_file.is_open()) in1_file.close();
	
	
	if(verbose) cerr << "Writing to file " << out_filename << endl;
	ofstream krona_file;
	krona_file.open(out_filename);    
	if(!krona_file.is_open()) {  cerr << "Could not open file " << out_filename << " for writing" << endl; exit(EXIT_FAILURE); }
	for(auto  it : node2hitcount) {
		uint64_t id = it.first;
		if(nodes.count(id)==0) {
			cerr << "Warning: Taxon ID " << id << " found in input file is not contained in taxonomic tree file "<< nodes_filename << ".\n";
			continue;
		}
		if(node2name.count(id)==0) {
			cerr << "Warning: Taxon ID " << id << " found in input file is not contained in names.dmp file "<< names_filename << ".\n";
			continue;
		}
		vector<string> lineage;
		lineage.push_back(node2name.at(id));
		while(nodes.count(id)>0 && id != nodes.at(id)) {
			if(node2name.count(nodes.at(id))==0) {
				cerr << "Warning: Taxon ID " << nodes.at(id) << " found in input file is not contained in names file "<< names_filename << ".\n";
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
		krona_file << num_unclassified << "\tUnclassified" << endl;
	}
	krona_file.close();

	return EXIT_SUCCESS;    
}


void usage(char *progname) { 
	fprintf(stderr, "Usage:\n   %s -t nodes.dmp -n names.dmp -i kaiju.out -o kaiju2krona.out\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of input file\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file, only required of -c is set to lca\n");
	fprintf(stderr, "   -n FILENAME   Name of names.dmp file.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -v            Enable verbose output.\n");
	fprintf(stderr, "   -u            Include count for unclassified reads in output.\n");
	exit(EXIT_FAILURE);
}

