/*
 Copyright (C) 2016 Héctor Condori Alagón.

 This file is part of readmapper, the massive Smith-Waterman pairwise sequence aligner.

 readmapper is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 readmapper is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with readmapper.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "readmapper.hpp"

#include "buffer.hpp"

class FASTAReader
{
private:
    int16_t default_ = 0;
    
    int16_t* seqs_;
    int maxSeqLength_;
    int numSeqs_;
    
    Buffer<int16_t> *seqs;
    std::vector<std::string> *seqs_ids;
    
    int vector_size;
    
    std::vector<int> seq_sizes;
    int current_size;
    int current_seq;
    int seq_length;
    int current_pos = 0;
    FILE* file_;
    int current;
    bool inheader;
    
    void add_char(char ch);
    void add_char_to_id(char ch);
    bool accept(char ch);
    bool acceptSpace();
    bool acceptPrintable();
    bool acceptPrintableBut(char ch);
public:
    FASTAReader(char* filename);
    FASTAReader(std::string filename);
    void setDefault(int value);
    //int16_t* getSeqs() { return this->seqs_; }
    bool next(Buffer<int16_t> *seqs, int *seqs_len, std::vector<std::string> *ids, int factor);
    bool next(Buffer<int16_t> *seqs, int *seqs_len, std::vector<std::string> *ids);
};

