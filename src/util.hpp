/* This file is part of Kaiju, Copyright 2015-2017 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef KAIJU_UTIL_H
#define KAIJU_UTIL_H

#include <stdlib.h>
#include <string>
#include <iostream>
#include <set>
#include <time.h>
#include <unordered_map>
#include <fstream>

#include "Config.hpp"

void error(const std::string e);

void strip(std::string &s);

bool isalpha(const char & c);

std::string getCurrentTime();

void parseNodesDmp(std::unordered_map<uint64_t,uint64_t> &, std::ifstream &);

void parseNodesDmpWithRank(std::unordered_map<uint64_t,uint64_t> &, std::unordered_map<uint64_t,std::string> &, std::ifstream &);

void parseNamesDmp(std::unordered_map<uint64_t,std::string> &, std::ifstream &);

std::string getTaxonNameFromId(std::unordered_map<uint64_t,std::string> &, uint64_t, std::string &);

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(std::unordered_map<uint64_t,uint64_t> &, const std::string &, const std::string &);

/* returns true if node1 is ancestor of node2  or if node1==node2*/
bool is_ancestor(std::unordered_map<uint64_t,uint64_t> &, uint64_t, uint64_t);

uint64_t lca_from_ids(Config *,std::unordered_map<uint64_t,unsigned int> &, std::set<uint64_t> &);

void readFMI(std::string fmi_filename, Config * config);

#endif
