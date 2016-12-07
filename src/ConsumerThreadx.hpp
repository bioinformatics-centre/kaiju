/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef CONSUMERTHREADX_H
#define CONSUMERTHREADX_H

#include <stdint.h>
#include <assert.h>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <list>
#include <cmath>
#include <algorithm>
#include <mutex>
#include <iostream>
#include <sstream>
#include <vector>
#include <iterator>
#include <string>
#include <cstring>
#include <climits>
#include <map>

#include "include/ProducerConsumerQueue/src/ProducerConsumerQueue.hpp"
#include "ReadItem.hpp"
#include "Config.hpp"
#include "ConsumerThread.hpp"


class ConsumerThreadx: public ConsumerThread  {
	protected:
	void classify_length();
	void classify_greedyblosum();
	void ids_from_SI(SI *);
	void ids_from_SI_recursive(SI *);
	std::set<char *> match_ids;

	public:
	ConsumerThreadx(ProducerConsumerQueue<ReadItem*>* workQueue, Config * config) : ConsumerThread(workQueue, config) { };
	void doWork();

};
#endif

