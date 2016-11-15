/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef KAIJU_UTIL_H
#define KAIJU_UTIL_H

#include <string>
#include <iostream>
#include <time.h>

void error(const std::string e);
void strip(std::string &s);
bool isalpha(const char & c);
std::string getCurrentTime();

#endif
