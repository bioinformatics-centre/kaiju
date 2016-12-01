/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include "util.hpp"

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

void parseNodesDmp(std::unordered_map<uint64_t,uint64_t> & nodes, std::ifstream & nodes_file) {
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
