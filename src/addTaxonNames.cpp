/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <getopt.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <map>
#include <algorithm>
#include <string>
#include <utility>
#include <stdexcept>
#include <deque>


void usage(char *progname);

using namespace std;

int main(int argc, char** argv) {


	unordered_map<uint64_t,uint64_t> nodes;
	unordered_map<uint64_t, string> node2name;
	unordered_map<uint64_t, string> node2rank;
	unordered_map<uint64_t, string> node2path;

	string nodes_filename = "";
	string names_filename = "";
	string in_filename = "";
	string out_filename;

	bool filter_unclassified = false;
	bool full_path = false;
	bool verbose = false;
	string rank;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hvpun:t:i:o:")) != -1) {
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
			default:
				usage(argv[0]);
		}
	}
	if(names_filename.length() == 0) { cerr << "Error: Please specify the location of the names.dmp file with the -n option."  << endl; usage(argv[0]); }
	if(nodes_filename.length() == 0) { cerr << "Error: Please specify the location of the nodes.dmp file with the -t option."  << endl; usage(argv[0]); }
	if(in_filename.length() == 0) { cerr << "Error: Please specify the location of the input file with the -i option."  << endl; usage(argv[0]); }

	ifstream nodes_file;
	nodes_file.open(nodes_filename);
	if(!nodes_file.is_open()) { cerr << "Error: Could not open file " << nodes_filename << endl; usage(argv[0]); }
	if(verbose) cerr << "Reading taxonomic tree from file " << nodes_filename << endl;
	string line;
	while(getline(nodes_file, line)) {
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
			end = line.find_first_not_of("abcdefghijklmnopqrstuvwxyz",start);
			//cerr << "end=" << end << "\t";
			string rank = line.substr(start,end-start);
			nodes.insert(make_pair(node,parent));
			node2rank.insert(make_pair(node,rank));
			//cerr << node << "\t" << parent << "\t" <<rank << "\n";

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


	ifstream in_file;
	in_file.open(in_filename);    
	if(!in_file.is_open()) {  cerr << "Could not open file " << in_filename << endl; exit(EXIT_FAILURE); }

	ostream * out_stream;
	if(out_filename.length()>0) {
		if(verbose) cerr << "Output file: " << out_filename << endl;
		ofstream * output_file = new ofstream();
		output_file->open(out_filename);    
		if(!output_file->is_open()) {  cerr << "Could not open file " << out_filename << " for writing" << endl; exit(EXIT_FAILURE); }
		out_stream = output_file;
	}
	else {
		out_stream = &cout;
	}

	if(verbose) cerr << "Processing " << in_filename <<"..." << "\n";
	
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
		try {
			uint64_t taxonid = stoul(line.substr(found,end-found));
			if(nodes.count(taxonid)==0) {
				cerr << "Warning: Taxon ID " << taxonid << " in output file is not contained in taxonomic tree file "<< nodes_filename << ".\n"; 
				continue;
			}
			if(node2name.count(taxonid)==0) {
				cerr << "Warning: Taxon ID " << taxonid << " in output file is not found in file "<< names_filename << ".\n";
				continue;
			}
			if(full_path) {
				if(node2path.count(taxonid)>0) { // look if path is already saved
					*out_stream << line << '\t' << node2path.at(taxonid) << "\n";
					continue;
				}
				deque<string> lineage;
				uint64_t id = taxonid;
				lineage.push_front(node2name.at(id));
				while(nodes.count(id)>0 && id != nodes.at(id)) {
					string name;
					if(node2name.count(nodes.at(id))==0) {
						cerr << "Warning: Taxon ID " << nodes.at(id) << " is not found in file "<< names_filename << ".\n";
						name = "n/a";
					}
					else {
						name = node2name.at(id);
					}
					id = nodes.at(id);
					lineage.push_front(name);
				}
				string lineage_text;
				for(auto  itl : lineage) {
					lineage_text += itl;
					lineage_text += "; ";
				}
				node2path.insert(make_pair(taxonid,lineage_text));
				*out_stream << line << '\t' << lineage_text << "\n";
			}
			else {
				*out_stream << line << '\t' << node2name.at(taxonid) << "\n";
			}
		}
		catch(const std::invalid_argument& ia) {
			cerr << "Found bad taxon id in line: " << line << endl; 
		}
		catch (const std::out_of_range& oor) {
			cerr << "Found bad taxon id (out of range error) in line: " << line << endl; 
		}
	}  // end while getline

	if(in_file.is_open()) {
		in_file.close();
	}
	out_stream->flush();
	if(out_filename.length()>0) {
		((ofstream*)out_stream)->close();
		delete ((ofstream*)out_stream);
	}

}


void usage(char *progname) { 
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
	fprintf(stderr, "   -v            Enable verbose output.\n");
	fprintf(stderr, "\n");
	exit(EXIT_FAILURE);
}

