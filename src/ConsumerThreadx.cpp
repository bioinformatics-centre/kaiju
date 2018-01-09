/* This file is part of Kaiju, Copyright 2015-2017 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include "ConsumerThreadx.hpp"



void ConsumerThreadx::classify_greedyblosum() {

		best_matches_SI.clear();
		best_matches.clear();
		best_match_score = 0;

		while(1) {
			Fragment * t = getNextFragment(best_match_score);
			if(!t) break;
			const std::string fragment = t->seq;
			const size_t length = fragment.length();
			const unsigned int num_mm = t->num_mm;

			if(config->debug) { std::cerr << "Searching fragment "<< fragment <<  " (" << length << ","<< num_mm << "," << t->diff << ")" << "\n"; }
			char * seq = new char[length+1];
			std::strcpy(seq, fragment.c_str());

			translate2numbers((uchar *)seq, (unsigned int)length, config->astruct);

			SI * si = NULL;
			if(num_mm > 0) {
				if(num_mm == config->mismatches) { //after last mm has been done, we need to have at least reached the min_length
					si = maxMatches_withStart(config->fmi, seq, (unsigned int)length, config->min_fragment_length, 1,t->si0,t->si1,t->matchlen);
				}
				else {
					si = maxMatches_withStart(config->fmi, seq, (unsigned int)length, t->matchlen, 1,t->si0,t->si1,t->matchlen);
				}
			}
			else {
				si = maxMatches(config->fmi, seq, (unsigned int)length, config->seed_length, 0); //initial matches
			}

			if(!si) {// no match for this fragment
				if(config->debug) std::cerr << "No match for this fragment." << "\n";
				delete[] seq;
				delete t;
				continue; // continue with the next fragment
			}
			if(config->debug) std::cerr << "Longest match has length " << (unsigned int)si->ql <<  "\n";

			if(config->mismatches > 0 && num_mm < config->mismatches) {
				SI * si_it = si;
				while(si_it) {
					unsigned int match_right_end = si_it->qi + si_it->ql - 1;
					if(num_mm > 0) assert(match_right_end  == length - 1); // greedy matches end always at the end
					if(config->debug) std::cerr << "Match from " << si_it->qi << " to " << match_right_end << ": " << fragment.substr(si_it->qi, match_right_end -  si_it->qi +1)  << " (" << si_it->ql << ")\n";
					if(si_it->qi > 0 && match_right_end + 1 >= config->min_fragment_length) {
						//1. match must end before beginning of fragment, i.e. it is extendable
						//2. remaining fragment, from zero to end of current match, must be longer than minimum length of accepted matches
						const size_t erase_pos = (match_right_end < length - 1) ? match_right_end + 1 : std::string::npos;
						addAllMismatchVariantsAtPosSI(t,(unsigned int)(si_it->qi - 1),erase_pos,si_it);
					}
					si_it = si_it->samelen ? si_it->samelen : si_it->next;
				}
			}

			if((unsigned int)si->ql < config->min_fragment_length) { // match was too short
				if(config->debug) { std::cerr << "Match of length " << si->ql << " is too short\n"; }
				delete[] seq;
				delete t;
				recursive_free_SI(si);
				continue; // continue with the next fragment
			}

			eval_match_scores(si, t);

			delete[] seq;
			delete t;

		} // end current fragment

		if(best_matches_SI.empty()) {
			return;
		}

		if(config->use_Evalue) {
			//calc e-value and only return match if > cutoff

			double bitscore = (LAMBDA * best_match_score - LN_K) / LN_2;
			double Evalue = config->db_length * query_len * pow(2, -1 * bitscore);
			if(config->debug) std::cerr << "E-value = " << Evalue << std::endl;

			if(Evalue > config->min_Evalue) {
				for(auto itm : best_matches_SI) {
					free(itm);
				}
				return;
			}
		}
		
		match_ids.clear();

		for(auto itm : best_matches_SI) {
			ids_from_SI(itm);
		}
		for(auto itm : best_matches_SI) {
			//recursive_free_SI(itm);
			free(itm);
		}

		std::stringstream ss;
		ss << best_match_score << "\t" ;
		for(auto it : match_ids) ss << it << ",";
		ss  << "\t";
		for(auto it : best_matches) ss << it << ",";
		extraoutput = ss.str();

}

void ConsumerThreadx::classify_length() {

		unsigned int longest_match_length = 0;
		longest_matches_SI.clear();
		longest_fragments.clear();

		while(1) {
			Fragment * t = getNextFragment(longest_match_length);
			if(!t) break;// searched all fragments that are longer than best match length
			const std::string fragment = t->seq;
			const unsigned int length = (unsigned int)fragment.length();

			if(config->debug) { std::cerr << "Searching fragment "<< fragment <<  " (" << length << ")" << "\n"; }
			char * seq = new char[length+1];
			std::strcpy(seq, fragment.c_str());

			translate2numbers((uchar *)seq, length, config->astruct);
			//use longest_match_length here too:
			SI * si = maxMatches(config->fmi, seq, length, std::max(config->min_fragment_length,longest_match_length),  1);

			if(!si) {// no match for this fragment
				if(config->debug) std::cerr << "No match for this fragment." << "\n";
				delete[] seq;
				delete t;
				continue; // continue with the next fragment
			}


			// just get length here and save si when it is longest
			if(config->debug) std::cerr << "Longest match is length " << (unsigned int)si->ql << "\n";
			if((unsigned int)si->ql > longest_match_length) {
				for(auto itm : longest_matches_SI)
					recursive_free_SI(itm);
				longest_matches_SI.clear();
				longest_matches_SI.push_back(si);
				longest_match_length = (unsigned int)si->ql;
				if(config->verbose) {
					longest_fragments.clear();
					longest_fragments.push_back(fragment.substr(si->qi,si->ql));
				}
			}
			else if((unsigned int)si->ql == longest_match_length) {
				longest_matches_SI.push_back(si);
				if(config->verbose)
					longest_fragments.push_back(fragment.substr(si->qi,si->ql));
			}
			else {
				recursive_free_SI(si);
				si = NULL;
			}
			delete[] seq;
			delete t;

		} // end current fragment

		if(longest_matches_SI.empty()) {
			return;
		}
		match_ids.clear();
		for(auto itm : longest_matches_SI) {
			ids_from_SI_recursive(itm);
		}
		for(auto itm : longest_matches_SI) {
			recursive_free_SI(itm);
		}

		std::stringstream ss;
		ss << longest_match_length << "\t";
		for(auto it : match_ids) ss << it << ",";
		ss  << "\t";
		for(auto it : longest_fragments) ss << it << ",";
		extraoutput = ss.str();

}

void ConsumerThreadx::doWork() {
	ReadItem * item = NULL;
	while(myWorkQueue->pop(&item)) {
		assert(item != NULL);
		read_count++;
		if(read_count > 20000) {
			flush_output();
			read_count = 0;
		}

		if((!item->paired && item->sequence1.length() < config->min_fragment_length*3) ||
			(item->paired && item->sequence1.length() < config->min_fragment_length*3 && item->sequence2.length() < config->min_fragment_length*3)) {
			output << "U\t" << item->name << "\t0\n";
			delete item;
			continue;
		}

		extraoutput = "";

		query_len = static_cast<double>(item->sequence1.length()) / 3.0;
		if(item->sequence1.length() >= config->min_fragment_length*3) {
			if(config->debug) std::cerr << "Getting fragments for read: "<< item->sequence1 << "\n";
			getAllFragmentsBits(item->sequence1);
		}
		if(item->paired) {
			query_len += static_cast<double>(item->sequence2.length()) / 3.0;
			if(item->sequence2.length() >= config->min_fragment_length*3) {
				if(config->debug) std::cerr << "Getting fragments for 2nd read: " << item->sequence2 << "\n";
				getAllFragmentsBits(item->sequence2);
			}
		}
		if(config->debug) std::cerr << fragments.size()  << " fragments found in the read."<< "\n";

		if(config->mode == MEM) {
			classify_length();
		}
		else if(config->mode == GREEDY) {
			classify_greedyblosum();
		}
		else { // this should not happen
			assert(false);
		}

		if(extraoutput.length() > 0) {
			output << "C\t" << item->name << "\t" << extraoutput << "\n";
			if(config->debug) {
				output << "C\t" << item->name << "\t" << extraoutput << "\n";
			}

		}
		else  {
			output << "U\t" << item->name << "\n";
			if(config->debug) {
				std::cerr << "U\t" << item->name << "\n";
			}

		}

		delete item;

		clearFragments();

	}

	flush_output();

}


void ConsumerThreadx::ids_from_SI(SI *si) {
	IndexType k, pos;
	int iseq;
	for (k=si->start; k<si->start+si->len; ++k) {
		if(match_ids.size() > config->max_match_ids) {
			break;
		}
		get_suffix(config->fmi, config->bwt->s, k, &iseq, &pos);
		match_ids.insert(config->bwt->s->ids[iseq]);
	}
}

void ConsumerThreadx::ids_from_SI_recursive(SI *si) {
	SI * si_it = si;
	while(si_it) {
		IndexType k, pos;
		int iseq;
		for (k=si_it->start; k<si_it->start+si_it->len; ++k) {
			if(match_ids.size() > config->max_match_ids) {
				break;
			}
			get_suffix(config->fmi, config->bwt->s, k, &iseq, &pos);
			match_ids.insert(config->bwt->s->ids[iseq]);
		} // end for
		si_it = si_it->samelen;
	} // end while all SI with same length
}
