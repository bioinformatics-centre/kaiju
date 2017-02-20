/* This file is part of Kaiju, Copyright 2015-2017 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include "Config.hpp"


Config::Config() { // constructor
}

Config::~Config() { // Destructor
	free(astruct->trans);
	free(astruct->a);
	free(astruct);
	if(SEG) {
		SegParametersFree(blast_seg_params);
	}
}

void Config::init() {
	astruct = alloc_AlphabetStruct(bwt->alphabet,0,0);
	if(SEG) {
		blast_seg_params = SegParametersNewAa();
		blast_seg_params->overlaps = TRUE;
	}
}


