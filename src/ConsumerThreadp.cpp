/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include "ConsumerThreadp.hpp"

void ConsumerThreadp::doWork() {        
	ReadItem * item = NULL;
	while(myWorkQueue->pop(&item)) {    	
		assert(item != NULL);
		count++;

		if((!item->paired && item->sequence1.length() < config->min_fragment_length) || 
			(item->paired && item->sequence1.length() < config->min_fragment_length && item->sequence2.length() < config->min_fragment_length)) {
			output << "U\t" << item->name << "\t0\n";
			delete item;
			continue;
		}		


		unsigned int score = calcScore(item->sequence1,0);
		if(score < config->min_score) {
			output << "U\t" << item->name << "\t0\n";
			delete item;
			continue;
		}

		fragments.insert(std::pair<unsigned int,Fragment *>(score,new Fragment(item->sequence1)));

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

		if(count%40000==0) flush_output();

	}

	flush_output();

}



