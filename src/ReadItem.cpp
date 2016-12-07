/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include "ReadItem.hpp"

ReadItem::ReadItem(std::string name, std::string sequence1) {
	this->name = name;
	this->sequence1 = sequence1;
}

ReadItem::ReadItem(std::string name, std::string sequence1, std::string sequence2) {
	this->paired = true;
	this->name = name;
	this->sequence1 = sequence1;
	this->sequence2 = sequence2;
}

