/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <locale>
#include <string>
#include <functional>

#include "ConsumerThread.hpp"

void usage(char *progname);
string strip(const string &s);
bool isalpha(char & c);
string calc_lca(unordered_map<uint64_t,uint64_t> *, string, string);

using namespace std;

int main(int argc, char** argv) {

	unordered_map<uint64_t,uint64_t> * nodes = new unordered_map<uint64_t,uint64_t>();

	string nodes_filename = "";
	string in1_filename = "";
	string in2_filename = "";
	string out_filename;
	string conflict = "1";

	ostream * out_stream;

	bool verbose = false;
	bool debug = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	char c;
	while ((c = getopt (argc, argv, "hdvc:n:t:i:j:o:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'v':
				verbose = true; break;
			case 'd':
				debug = true; break;
			case 'c':
				conflict = optarg; break;
			case 'o':
				out_filename = optarg; break;
			case 't':
				nodes_filename = optarg; break;
			case 'i':
				in1_filename = optarg; break;
			case 'j':
				in2_filename = optarg; break;
			default:
								usage(argv[0]);
		}
	}
	if(!(conflict=="1" || conflict=="2" || conflict=="lca")) { cerr << "Error: Value of argument -c must be either 1, 2, or lca." <<endl; usage(argv[0]); }
	if(conflict=="lca" && nodes_filename.length() == 0) { cerr << "Error: LCA mode requires the name of the nodes.dmp file, using the -t option."  << endl; usage(argv[0]); }
	if(in1_filename.length() == 0) { cerr << "Error: Specify the name of the first input file, using the -i option."  << endl; usage(argv[0]); }
	if(in2_filename.length() == 0) { cerr << "Error: Specify the name of the second input file, using the -j option."  << endl; usage(argv[0]); }
	
	if(nodes_filename.length() > 0) {
		ifstream nodes_file;
		nodes_file.open(nodes_filename.c_str());
		if(!nodes_file.is_open()) { cerr << "Error: Could not open file " << nodes_filename << endl; usage(argv[0]); }
		string line;
		while(getline(nodes_file, line)) {
			if(line.length() == 0) { continue; }
			try {
				size_t end = line.find_first_not_of("0123456789");
				uint64_t node = stoul(line.substr(0,end));
				size_t start = line.find_first_of("0123456789",end);
				end = line.find_first_not_of("0123456789",start+1);
				uint64_t parent = stoul(line.substr(start,end-start));
				nodes->insert(make_pair(node,parent));  //maybe the nodes->at(node) = parent;  would be faster?!
			}
			catch(const std::invalid_argument& ia) {
				cerr << "Found bad number in line: " << line << endl; 
			}
			catch (const std::out_of_range& oor) {
				cerr << "Found bad number (out of range error) in line: " << line << endl; 
			}
		}
		nodes_file.close();
	}


	if(out_filename.length()>0) {
		ofstream * filestream = new ofstream();
		filestream->open(out_filename);    
		if(!filestream->is_open()) {  cerr << "Could not open file " << out_filename << " for writing" << endl; exit(EXIT_FAILURE); }
		out_stream = filestream;
	}
	else {
		out_stream = &cout;
	}

	ifstream in1_file, in2_file;
	in1_file.open(in1_filename);    
	if(!in1_file.is_open()) {  cerr << "Could not open file " << in1_filename << endl; exit(EXIT_FAILURE); }
	in2_file.open(in2_filename);    
	if(!in2_file.is_open()) {  cerr << "Could not open file " << in2_filename << endl; exit(EXIT_FAILURE); }

	
	string line;
	line.reserve(500);

	uint count = 0;

	uint countC1 = 0;
	uint countC2 = 0;
	uint countC12 = 0;
	uint countC3 = 0;
	uint countC1notC2 = 0;
	uint countC2notC1 = 0;


	while(getline(in1_file,line)) {                		
		count++;
	
		//if(debug) cerr << "Count=" << count << endl;
		//if(debug) cerr << line << endl;
		// get the three values
		char classified1 = line[0];
		size_t index1 = line.find('\t');
		if(index1 == string::npos) { cerr << "Error Could not parse line " << count << " in file " << in1_filename << endl; break; }
		size_t index2 = line.find('\t',index1+1);
		if(index2 == string::npos) { cerr << "Error Could not parse line " << count << " in file " << in1_filename << endl; break; }
		string name1 = line.substr(index1+1,index2-index1-1);
		size_t end = line.find_first_not_of("0123456789",index2+1);
		if(end == string::npos) { if(index2 < line.length()-1) end=line.length()-1; else { cerr << "Error Could not parse line " << count << " in file " << in1_filename << endl; break; }}
		string taxon_id1 = line.substr(index2+1,end-index2);

		//if(debug) cerr << "Name1=" << name1 <<" ID1=" << taxon_id1 << endl;

		if(!getline(in2_file,line)) {
			//that's the border case where file1 has more entries than file2
			cerr << "Error: File " << in1_filename <<" has more lines then file " << in2_filename  <<endl;
			break; 
		}
		//if(debug) cerr << line << endl;
		
		// get the three values for second file
		char classified2 = line[0];
		index1 = line.find('\t');
		if(index1 == string::npos) { cerr << "Error Could not parse line " << count << " in file " << in2_filename << endl; break; }
		index2 = line.find('\t',index1+1);
		if(index2 == string::npos) { cerr << "Error Could not parse line " << count << " in file " << in2_filename << endl; break; }
		string name2 = line.substr(index1+1,index2-index1-1);
		end = line.find_first_not_of("0123456789",index2+1);
		if(end == string::npos) { if(index2 < line.length()-1) end=line.length()-1; else { cerr << "Error Could not parse line " << count << " in file " << in2_filename << endl; break; }}
		string taxon_id2 = line.substr(index2+1,end-index2);

		//if(debug) cerr << "Name2=" << name1 <<" ID2=" << taxon_id2 << endl;
		
		if(name1 != name2) {
			cerr << "Error: Read names are not identical between the two input files" << endl;
			break;
		}
		if(!(classified1=='C' || classified1 =='U')) {
			cerr << "Error: Line " << count << " in file "<< in1_filename << " does not start with C or U. "<< endl;
			break;
		}
		if(!(classified2=='C' || classified2 =='U')) {
			cerr << "Error: Line " << count << " in file "<< in2_filename << " does not start with C or U. "<< endl;
			break;
		}

		if(classified1=='C' && classified2=='C')  {
			string lca;
			if(taxon_id1==taxon_id2) {
				lca = taxon_id1;
			}
			else { // different taxon ids in both files
				if(conflict=="1")
					lca = taxon_id1;
				else if(conflict=="2")
					lca = taxon_id2;
				else {
					assert(conflict=="lca");
					lca = calc_lca(nodes, taxon_id1, taxon_id2); 
					if(lca=="0") { cerr << "Error while calculating LCA of " << taxon_id1 << " and " << taxon_id2 << " in line " << count << endl; break; }
					if(debug) cerr << "LCA of "<< taxon_id1 << " and " << taxon_id2 << " is "  << lca << endl;
				}
			}
			countC1++; countC2++; countC12++; countC3++;
			(*out_stream) << "C" << "\t" << name1 << "\t" << lca <<"\n";
		}
		else if(classified1=='C') {
			assert(classified2=='U');
			countC1++; countC1notC2++; countC3++;
			(*out_stream) << "C" << "\t" << name1 << "\t" << taxon_id1 <<"\n";

		}
		else if(classified2=='C') {
			assert(classified1=='U');
			countC2++; countC2notC1++; countC3++;
			(*out_stream) << "C" << "\t" << name1 << "\t" << taxon_id2 <<"\n";

		}
		else {
			assert(classified1=='U' && classified2=='U');
			(*out_stream) << "U" << "\t" << name1 << "\t0\n";
		}

		if(count%20000==0) {
			out_stream->flush();
		}

	} // end main loop around file1


	if(in1_file.is_open()) in1_file.close();
	if(in2_file.is_open()) {
		if(getline(in2_file,line) && line.length()>0) {
			cerr << "Warning: File " << in2_filename <<" has more lines then file " << in1_filename  <<endl;
		}
		in2_file.close();
	}

	out_stream->flush();
	if(out_filename.length()>0) {
		((ofstream*)out_stream)->close();
		delete ((ofstream*)out_stream);
	}

	if(verbose) {
		fprintf(stderr, "Number of all reads in input:\t%10u\n",count);
		fprintf(stderr, "         classified in file1:\t%10u  %6.2f%%\n",countC1,((double)countC1/(double)count*100.0));
		fprintf(stderr, "            but not in file2:\t%10u  %6.2f%%\n",countC1notC2,((double)countC1notC2/(double)count*100.0));
		fprintf(stderr, "         classified in file2:\t%10u  %6.2f%%\n",countC2,((double)countC2/(double)count*100.0));
		fprintf(stderr, "            but not in file1:\t%10u  %6.2f%%\n",countC2notC1,((double)countC2notC1/(double)count*100.0));
		fprintf(stderr, "          classified in both:\t%10u  %6.2f%%\n",countC12,((double)countC12/(double)count*100.0));
		fprintf(stderr, "         combined classified:\t%10u  %6.2f%%\n",countC3,((double)countC3/(double)count*100.0));
	}

	delete nodes;
	return EXIT_SUCCESS;    
}

inline bool isalpha(char & c) { 
	return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z');
}

string strip(const string &s) {
		string result;
		result.reserve(s.length());
		for(auto it=s.cbegin(); it!=s.cend(); ++it) {
			if(isalpha(*it))
				result += *it;
		}
		return result;
}

void usage(char *progname) { 
	fprintf(stderr, "Usage:\n   %s -i in1.tsv -j in2.tsv [-o outfile.tsv] [-c 1|2|lca] [-t nodes.dmp] [-v] [-d]\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of first input file\n");
	fprintf(stderr, "   -j FILENAME   Name of second input file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "   -c STRING     Conflict resolution mode, must be 1, 2 or lca (default: 1)\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file, only required when -c is set to lca\n");
	fprintf(stderr, "   -v            Enable verbose output, which will print a summary in the end.\n");
	fprintf(stderr, "   -d            Enable debug output.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "NOTE: Both input files need to be sorted by the read name in the second column.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "The option -c determines the method of resolving conflicts in the taxonomic assignment for a read.\n");
	fprintf(stderr, "If set to 1, then the taxon id from the first input file is used.\n");
	fprintf(stderr, "If set to 2, then the taxon id from the second input file is used.\n");
	fprintf(stderr, "If set to lca, then the LCA of the two taxon ids from both input files is used. This option requires to set the name of the nodes.dmp file.\n");
	exit(EXIT_FAILURE);
}

string calc_lca(unordered_map<uint64_t,uint64_t> * nodes, string id1, string id2) {

		uint64_t node1;
		uint64_t node2; 
		try {
			node1 = stoul(id1);
			node2 = stoul(id2);
		}
		catch(const std::invalid_argument& ia) {
			cerr << "Bad number in taxon id" << endl; 
			return 0;
		}
		catch (const std::out_of_range& oor) {
			cerr << "Bad number (out of range error) in taxon id" << endl; 
			return 0;
		}

		unordered_set<uint64_t> lineage1;
		lineage1.insert(node1);
		while(nodes->count(node1)>0 && node1 != nodes->at(node1)) {
			lineage1.insert(nodes->at(node1));
			node1 = nodes->at(node1);	
		}
		
		uint64_t lca = node2;
		do {
			lca = nodes->at(lca);	
		} while(lineage1.count(lca)==0 && lca != nodes->at(lca));


		return to_string(lca);
}


