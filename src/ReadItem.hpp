/* This file is part of Kaiju, Copyright 2015 Peter Menzel and Anders Krogh,
 * Kaiju is licensed under the GPLv3, see the file LICENSE. */

#ifndef READ_ITEM_H
#define READ_ITEM_H

#include <string>

using namespace std;
class ReadItem {
    public:
        string name;
        bool paired = false;
        string sequence1;
        string sequence2;    
        ReadItem(string name, string sequence1);
        ReadItem(string name, string sequence1, string sequence2);
};
#endif
