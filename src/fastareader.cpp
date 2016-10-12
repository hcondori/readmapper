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

#include "fastareader.hpp"

FASTAReader::FASTAReader(char* filename) :
    file_(fopen(filename, "r")), seq_sizes(VSIZE)
{
    current = fgetc_unlocked(file_);
}

FASTAReader::FASTAReader(std::string filename) :
    file_(fopen(filename.c_str(), "r")), seq_sizes(VSIZE)
{
    current = fgetc_unlocked(file_);
}

void FASTAReader::add_char(char ch)
{
    int seq_length = vector_size / VSIZE;
    
    if(this->seq_sizes[current_seq] == seq_length)
    {
        vector_size = 2 * vector_size;
        this->seqs->resize(vector_size);
    }
    (*this->seqs)[VSIZE * seq_sizes[current_seq] + current_seq] = current;
    seq_sizes[current_seq]++;
    current = fgetc_unlocked(file_);
}

void FASTAReader::add_char_to_id(char ch)
{
    (*this->seqs_ids)[current_seq].push_back(ch);
    current = fgetc_unlocked(file_);
}

bool FASTAReader::accept(char ch)
{
    if(current == ch)
    {
        current = fgetc_unlocked(file_);
        return true;
    }
    //current = fgetc_unlocked(file_);
    return false;
}

bool FASTAReader::acceptSpace()
{
    if(isspace(current))
    {
        current = fgetc_unlocked(file_);
        return true;
    }
    return false;
}


bool FASTAReader::acceptPrintable()
{
    if(isprint(current))
    {
        //std::cout << "Hello!!!" << std::endl;
        if(!inheader)
            add_char(current);
        else
            add_char_to_id(current);
            //current = fgetc_unlocked(file_);
        return true;
    }
    return false;
}

bool FASTAReader::acceptPrintableBut(char ch)
{
    if(isprint(current) && current != ch)
    {
        if(!inheader)
            add_char(current);
        else
            current = fgetc_unlocked(file_);
        return true;
    }
    return false;
}

bool FASTAReader::next(Buffer<int16_t> *seqs, int *seqs_len, 
                       std::vector<std::string> *ids, int factor)
{
    this->seqs_ids = ids;
    
    
    for(std::string &id : (*seqs_ids))
        id.clear();
        
    this->seqs = seqs;
    this->vector_size = seqs->size();
    bool status = false;
    //this->seq_length = VSIZE / seq_count;
    std::fill(seq_sizes.begin(), seq_sizes.end(), this->default_);
    seqs->clear(default_);
    for(current_seq = 0; current_seq < VSIZE; current_seq++)
    {
        if(accept('>'))
        {
            inheader = true;
            while(acceptPrintable());
            acceptSpace();
            inheader = false;
            while(acceptPrintableBut('>') || acceptSpace());
            status = true;
        }
    }
    
    std::copy_n(seq_sizes.data(), VSIZE, seqs_len);
    return status;
}

bool FASTAReader::next(Buffer<int16_t> *seqs, int *seqs_len,
                       std::vector<std::string> *ids)
{
    return this->next(seqs, seqs_len, ids, 1);
}

void FASTAReader::setDefault(int value)
{
    this->default_ = value;
}
