/* This file is part of Kaiju, Copyright 2015-2019 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <set>
#include <locale>
#include <string>
#include <stdexcept>
#include <cstdarg>

#include "util.hpp"

void usage(const char * progname);
std::string calc_lca(std::unordered_map<uint64_t,uint64_t> &, const std::string &, const std::string &);

int main(int argc, char** argv) {

	std::unordered_map<uint64_t,uint64_t> nodes;

	std::string nodes_filename = "";
	std::string in1_filename = "";
	std::string in2_filename = "";
	std::string out_filename;
	std::string conflict = "lca";

	std::ostream * out_stream;

	bool verbose = false;
	bool debug = false;
	bool use_score = false;

	// --------------------- START ------------------------------------------------------------------
	// Read command line params
	int c;
	while ((c = getopt(argc, argv, "hdvsc:n:t:i:j:o:")) != -1) {
		switch (c)  {
			case 'h':
				usage(argv[0]);
			case 'v':
				verbose = true; break;
			case 'd':
				debug = true; break;
			case 's':
				use_score = true; break;
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
	if(!(conflict=="1" || conflict=="2" || conflict=="lca" || conflict=="lowest")) { error("Value of argument -c must be either 1, 2, lca, or lowest."); usage(argv[0]); }
	if(conflict=="lca" && nodes_filename.length() == 0) { error("Error: LCA mode requires the name of the nodes.dmp file, using the -t option."); usage(argv[0]); }
	if(in1_filename.length() == 0) { error("Specify the name of the first input file, using the -i option."); usage(argv[0]); }
	if(in2_filename.length() == 0) { error("Specify the name of the second input file, using the -j option."); usage(argv[0]); }

	if(nodes_filename.length() > 0) {
		std::ifstream nodes_file;
		nodes_file.open(nodes_filename);
		if(!nodes_file.is_open()) { std::cerr << "Error: Could not open file " << nodes_filename << std::endl; usage(argv[0]); }
		parseNodesDmp(nodes,nodes_file);
		nodes_file.close();
	}


	if(out_filename.length()>0) {
		std::ofstream * filestream = new std::ofstream();
		filestream->open(out_filename);
		if(!filestream->is_open()) {  std::cerr << "Could not open file " << out_filename << " for writing" << std::endl; exit(EXIT_FAILURE); }
		out_stream = filestream;
	}
	else {
		out_stream = &std::cout;
	}

	std::ifstream in1_file, in2_file;
	in1_file.open(in1_filename);
	if(!in1_file.is_open()) {  std::cerr << "Could not open file " << in1_filename << std::endl; exit(EXIT_FAILURE); }
	in2_file.open(in2_filename);
	if(!in2_file.is_open()) {  std::cerr << "Could not open file " << in2_filename << std::endl; exit(EXIT_FAILURE); }


	std::string line;
	line.reserve(500);

	unsigned int count = 0;

	unsigned int countC1 = 0;
	unsigned int countC2 = 0;
	unsigned int countC12 = 0;
	unsigned int countC3 = 0;
	unsigned int countC1notC2 = 0;
	unsigned int countC2notC1 = 0;


	while(getline(in1_file,line)) {
		count++;

		//if(debug) std::cerr << "Count=" << count << std::endl;
		//if(debug) std::cerr << line << std::endl;
		// get the three values
		char classified1 = line[0];
		if(!(classified1=='C' || classified1 =='U')) { std::cerr << "Error: Line " << count << " in file "<< in1_filename << " does not start with C or U. "<< std::endl; break; }
		size_t index_tab1 = line.find('\t');
		if(index_tab1 == std::string::npos) { std::cerr << "Error: Could not parse line " << count << " in file " << in1_filename << std::endl; break; }
		size_t index_tab2 = line.find('\t', index_tab1 + 1);
		if(index_tab2 == std::string::npos) { std::cerr << "Error: Could not parse line " << count << " in file " << in1_filename << std::endl; break; }
		std::string name1 = line.substr(index_tab1 + 1, index_tab2 - index_tab1 - 1);
		std::string taxon_id1;
		std::string score1s = "0";
		if(use_score && classified1=='C') { // look for third tab
			size_t index_tab3 = line.find('\t', index_tab2 + 1);
			if(index_tab3 == std::string::npos) { std::cerr << "Error: No score column (4th col) found in line " << count << " in file " << in1_filename << std::endl; break; }
			// taxon id is between tab2 and tab3
			taxon_id1 = line.substr(index_tab2 + 1, index_tab3 - index_tab2 - 1);
			// score follows after tab3
			size_t end = line.find_first_not_of(".0123456789",index_tab3 + 1);
			if(end == std::string::npos) {
				if(index_tab3 < line.length()-1){
					end = line.length()-1;
				}
				else { std::cerr << "Error Could not parse line " << count << " in file " << in1_filename << std::endl; break; }
			}
			else {
				end -= 1;
			}
			score1s = line.substr(index_tab3 + 1, end - index_tab3);
		}
		else { // look to the end of the taxon id (which starts after second tab), it can either be the line end or some other columns following
			size_t end = line.find_first_not_of("0123456789",index_tab2 + 1);
			if(end == std::string::npos) {
				if(index_tab2 < line.length()-1){
					end = line.length()-1;
				}
				else { std::cerr << "Error Could not parse line " << count << " in file " << in1_filename << std::endl; break; }
			}
			else {
				end -= 1;
			}
			taxon_id1 = line.substr(index_tab2 + 1, end - index_tab2);
		}

		if(debug) std::cerr << "Name1=" << name1 <<" ID1=" << taxon_id1 << (use_score ? " Score="+score1s+"\n" : "\n");

		if(!getline(in2_file,line)) {
			//that's the border case where file1 has more entries than file2
			std::cerr << "Error: File " << in1_filename <<" has more lines then file " << in2_filename  <<std::endl;
			break;
		}
		//if(debug) std::cerr << line << std::endl;

		// get the three values for second file
		char classified2 = line[0];
		if(!(classified2=='C' || classified2 =='U')) { std::cerr << "Error: Line " << count << " in file "<< in2_filename << " does not start with C or U. "<< std::endl; break; }
		index_tab1 = line.find('\t');
		if(index_tab1 == std::string::npos) { std::cerr << "Error: Could not parse line " << count << " in file " << in1_filename << std::endl; break; }
		index_tab2 = line.find('\t', index_tab1 + 1);
		if(index_tab2 == std::string::npos) { std::cerr << "Error: Could not parse line " << count << " in file " << in1_filename << std::endl; break; }
		std::string name2 = line.substr(index_tab1 + 1, index_tab2 - index_tab1 - 1);
		std::string taxon_id2;
		std::string score2s = "0";
		if(use_score and classified2=='C') { // look for third tab
			size_t index_tab3 = line.find('\t', index_tab2 + 1);
			if(index_tab3 == std::string::npos) { std::cerr << "Error: No score column (4th col) found in line " << count << " in file " << in2_filename << std::endl; break; }
			// taxon id is between tab2 and tab3
			taxon_id2 = line.substr(index_tab2 + 1, index_tab3 - index_tab2 - 1);
			// score follows after tab3
			size_t end = line.find_first_not_of(".0123456789",index_tab3 + 1);
			if(end == std::string::npos) {
				if(index_tab3 < line.length()-1){
					end = line.length()-1;
				}
				else { std::cerr << "Error Could not parse line " << count << " in file " << in1_filename << std::endl; break; }
			}
			else {
				end -= 1;
			}
			score2s = line.substr(index_tab3 + 1, end - index_tab3);
		}
		else { // look to the end of the taxon id (which starts after second tab), it can either be the line end or some other columns following
			size_t end = line.find_first_not_of("0123456789",index_tab2 + 1);
			if(end == std::string::npos) {
				if(index_tab2 < line.length()-1){
					end = line.length()-1;
				}
				else { std::cerr << "Error Could not parse line " << count << " in file " << in1_filename << std::endl; break; }
			}
			else {
				end -= 1;
			}
			taxon_id2 = line.substr(index_tab2 + 1, end - index_tab2);
		}

		if(debug) std::cerr << "Name2=" << name2 <<" ID2=" << taxon_id2 << (use_score ? " Score="+score2s+"\n" : "\n");

		if(name1 != name2) {
			std::cerr << "Error: Read names are not identical between the two input files on line " << count << std::endl;
			break;
		}

		if(classified1=='C' && classified2=='C')  {
			std::string lca;

			std::string output_score;
			bool same_score = true;
			double score1d = 0.0;
			double score2d = 0.0;
			// parse score strings
			if(use_score) {
				try { score1d = std::stod(score1s); } catch (const std::exception&) { std::cerr << "Error while parsing score on line " << count << "in file " << in1_filename << std::endl; }
				try { score2d = std::stod(score2s); } catch (const std::exception&) { std::cerr << "Error while parsing score on line " << count << "in file " << in2_filename << std::endl; }
				if(score1d != score2d) { same_score = false; }
			}

			if(taxon_id1==taxon_id2) {
				lca = taxon_id1;
				if(use_score) { output_score = (score2d > score1d) ? score2s : score1s; }
			}
			else { // different taxon ids in both files
				if(same_score) {
					output_score = score1s;
					if(conflict=="1") {
						lca = taxon_id1;
					}
					else if(conflict=="2") {
						lca = taxon_id2;
					}
					else if(conflict=="lowest") {
						if(is_ancestor(nodes,taxon_id1,taxon_id2)) {
							lca = taxon_id2;
						}
						else if(is_ancestor(nodes,taxon_id2,taxon_id1)) {
							lca = taxon_id1;
						}
						else {
							lca = calc_lca(nodes, taxon_id1, taxon_id2);
						}
						if(lca=="0") { std::cerr << "Error while calculating lowest node of " << taxon_id1 << " and " << taxon_id2 << " in line " << count << ", setting taxon id to " << taxon_id1 << std::endl; lca=taxon_id1; }
					}
					else {
						assert(conflict=="lca");
						lca = calc_lca(nodes, taxon_id1, taxon_id2);
						if(lca=="0") { std::cerr << "Error while calculating lowest node of " << taxon_id1 << " and " << taxon_id2 << " in line " << count << ", setting taxon id to " << taxon_id1 << std::endl; lca=taxon_id1; }
						if(debug) std::cerr << "LCA of "<< taxon_id1 << " and " << taxon_id2 << " is "  << lca << std::endl;
					}
				}
				else { //different scores
					if(score1d > score2d) {
						lca = taxon_id1;
						output_score = score1s;
					}
					else {
						assert(score1d < score2d);
						lca = taxon_id2;
						output_score = score2s;
					}
				}
			}
			countC1++; countC2++; countC12++; countC3++;
			(*out_stream) << "C" << "\t" << name1 << "\t" << lca << (use_score ? "\t"+output_score+"\n" : "\n");
		}
		else if(classified1=='C') {
			assert(classified2=='U');
			countC1++; countC1notC2++; countC3++;
			(*out_stream) << "C" << "\t" << name1 << "\t" << taxon_id1 << (use_score ? "\t"+score1s+"\n" : "\n");
		}
		else if(classified2=='C') {
			assert(classified1=='U');
			countC2++; countC2notC1++; countC3++;
			(*out_stream) << "C" << "\t" << name1 << "\t" << taxon_id2 << (use_score ? "\t"+score2s+"\n" : "\n");
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
			std::cerr << "Warning: File " << in2_filename <<" has more lines then file " << in1_filename  <<std::endl;
		}
		in2_file.close();
	}

	out_stream->flush();
	if(out_filename.length()>0) {
		((std::ofstream*)out_stream)->close();
		delete ((std::ofstream*)out_stream);
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

	return EXIT_SUCCESS;
}

void usage(const char * progname) {
	print_usage_header();
	fprintf(stderr, "Usage:\n   %s -i in1.tsv -j in2.tsv [-o outfile.tsv] [-c 1|2|lca|lowest] [-s] [-t nodes.dmp] [-v] [-d]\n", progname);
	fprintf(stderr, "\n");
	fprintf(stderr, "Mandatory arguments:\n");
	fprintf(stderr, "   -i FILENAME   Name of first input file\n");
	fprintf(stderr, "   -j FILENAME   Name of second input file\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "Optional arguments:\n");
	fprintf(stderr, "   -o FILENAME   Name of output file.\n");
	fprintf(stderr, "   -c STRING     Conflict resolution mode, must be 1, 2,  lca, or lowest (default: lca)\n");
	fprintf(stderr, "   -t FILENAME   Name of nodes.dmp file, only required when -c is set to lca\n");
	fprintf(stderr, "   -s            Use 4th column with classification score to give precedence to taxon with better score.\n");
	fprintf(stderr, "   -v            Enable verbose output, which will print a summary in the end.\n");
	fprintf(stderr, "   -d            Enable debug output.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "NOTE: Both input files need to be sorted by the read name in the second column.\n");
	fprintf(stderr, "\n");
	fprintf(stderr, "The option -c determines the method of resolving conflicts in the taxonomic assignment for a read.\n");
	fprintf(stderr, "Possible values are '1', '2', 'lca', 'lowest':\n");
	fprintf(stderr, "  '1' -> the taxon id from the first input file is used.\n");
	fprintf(stderr, "  '2' -> the taxon id from the second input file is used.\n");
	fprintf(stderr, "  'lca' -> the least common ancestor of the two taxon ids from both input files is used.\n");
	fprintf(stderr, "  'lowest' -> the lower rank of the two taxa is used if they are within the same lineage. Otherwise the LCA is used.\n");
	fprintf(stderr, "When using values 'lca' or 'lowest', the path to the file nodes.dmp needs to be specified via option -t.\n");
	exit(EXIT_FAILURE);
}



std::string calc_lca(std::unordered_map<uint64_t,uint64_t> & nodes, const std::string & id1, const std::string & id2) {

		uint64_t node1;
		uint64_t node2;
		try {
			node1 = stoul(id1);
			node2 = stoul(id2);
		}
		catch(const std::invalid_argument& ia) {
			std::cerr << "Warning: Bad number in taxon id" << std::endl;
			return "0";
		}
		catch (const std::out_of_range& oor) {
			std::cerr << "Warning: Bad number (out of range error) in taxon id" << std::endl;
			return "0";
		}

		if(nodes.count(node1)==0 && nodes.count(node2)==0) {
			std::cerr << "Warning: Taxon IDs " << node1 << " and " << node2 << " from both input files are not contained in taxonomic tree.\n";
			return "0";
		}
		else if(nodes.count(node1)==0) {
			std::cerr << "Warning: Taxon ID " << node1 << " in first input file is not contained in taxonomic tree.\n";
			return std::to_string(node2);
		}
		else if(nodes.count(node2)==0) {
			std::cerr << "Warning: Taxon ID " << node2 << " in second input file is not contained in taxonomic tree.\n";
			return std::to_string(node1);
		}

		std::set<uint64_t> lineage1;
		lineage1.emplace(node1);
		while(nodes.count(node1)>0 && node1 != nodes.at(node1)) {
			lineage1.emplace(nodes.at(node1));
			node1 = nodes.at(node1);
		}

		uint64_t lca = node2;
		do {
			lca = nodes.at(lca);
		} while(lineage1.count(lca)==0 && lca != nodes.at(lca));


		return std::to_string(lca);
}


