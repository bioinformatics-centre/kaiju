/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef KAIJU_UTIL_H
#define KAIJU_UTIL_H

#include <string>
#include <iostream>
#include <time.h>
#include <unordered_map>
#include <fstream>

void error(const std::string e);
void strip(std::string &s);
bool isalpha(const char & c);
std::string getCurrentTime();

void parseNodesDmp(std::unordered_map<uint64_t,uint64_t> &, std::ifstream &);
void parseNodesDmpWithRank(std::unordered_map<uint64_t,uint64_t> &, std::unordered_map<uint64_t,std::string> &, std::ifstream &);
void parseNamesDmp(std::unordered_map<uint64_t,std::string> &, std::ifstream &);
std::string getTaxonNameFromId(std::unordered_map<uint64_t,std::string> &, uint64_t, std::string &);

#endif
