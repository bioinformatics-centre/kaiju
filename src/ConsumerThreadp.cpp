/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include "ConsumerThreadp.hpp"

void ConsumerThreadp::doWork() {
	ReadItem * item = NULL;
	while(myWorkQueue->pop(&item)) {
		assert(item != NULL);
		read_count++;

		if(item->sequence1.length() < config->min_fragment_length) {
			output << "U\t" << item->name << "\t0\n";
			delete item;
			continue;
		}

		for (auto & c: item->sequence1) {
			c = (char)toupper(c);
		}
		size_t start = 0;
		size_t pos = item->sequence1.find_first_not_of("ACDEFGHIKLMNPQRSTVWY");
		while(pos != std::string::npos) {
			if(pos-start >= config->min_fragment_length) {
				std::string subseq =  item->sequence1.substr(start,pos-start);
				//cerr << "subseq=" << subseq << endl;
				if(config->mode==GREEDYBLOSUM) {
					const unsigned int score = calcScore(subseq);
					if(score >= config->min_score) {
						fragments.insert(std::pair<unsigned int,Fragment *>(score,new Fragment(subseq)));
					}
				}
				else {
					fragments.insert(std::pair<unsigned int,Fragment *>(subseq.length(),new Fragment(subseq)));
				}
			}
			start = pos+1;
			pos = item->sequence1.find_first_not_of("ACDEFGHIKLMNPQRSTVWY", pos + 1);
		}
		//add remaining sequence, which corresponds to the whole sequence if no invalid char was found
		std::string subseq = item->sequence1.substr(start,item->sequence1.length()-start);
		if(subseq.length() >= config->min_fragment_length) {
			if(config->mode==GREEDYBLOSUM) {
				const unsigned int score = calcScore(subseq);
				if(score >= config->min_score) {
					fragments.insert(std::pair<unsigned int,Fragment *>(score,new Fragment(subseq)));
				}
			}
			else {
				fragments.insert(std::pair<unsigned int,Fragment *>(subseq.length(),new Fragment(subseq)));
			}
		}

		if(config->debug) std::cerr << fragments.size()  << " fragments found in the read."<< "\n";

		if(fragments.empty()) {
			output << "U\t" << item->name << "\t0\n";
			delete item;
			continue;
		}

		extraoutput = "";
		if(config->mode == MEM) {
			classify_length();
		}
		else if(config->mode == GREEDYBLOSUM) {
			classify_greedyblosum();
		}
		else { // this should not happen
			assert(false);
		}

		if(extraoutput.length() > 0) {
			output << "C\t" << item->name << "\t" << extraoutput << "\n";
		}
		else  {
			output << "U\t" << item->name << "\n";
		}

		delete item;

		clearFragments();

		if(read_count==40000) {
			flush_output();
			read_count = 0;
		}

	}

	flush_output();

}



