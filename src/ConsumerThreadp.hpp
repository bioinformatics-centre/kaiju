/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef CONSUMERTHREADP_H
#define CONSUMERTHREADP_H

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
#include "ConsumerThreadx.hpp"


class ConsumerThreadp: public ConsumerThreadx  {

	public:
	ConsumerThreadp(ProducerConsumerQueue<ReadItem*>* workQueue, Config * config) : ConsumerThreadx(workQueue, config) { };
	void doWork();

};
#endif

