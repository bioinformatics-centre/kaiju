/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#include "ReadItem.hpp"

ReadItem::ReadItem(const std::string & n, const std::string & s1) : name(n), sequence1(s1) { }

ReadItem::ReadItem(const std::string & n, const std::string & s1, const std::string & s2) : name(n), sequence1(s1), sequence2(s2), paired(true) { }


