/* This file is part of Kaiju, Copyright 2015,2016 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef READ_ITEM_H
#define READ_ITEM_H

#include <string>

class ReadItem {
    public:
        std::string name;
        bool paired = false;
        std::string sequence1;
        std::string sequence2;
        ReadItem(std::string name, std::string sequence1);
        ReadItem(std::string name, std::string sequence1, std::string sequence2);
};
#endif
