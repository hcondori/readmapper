#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "readmapper.hpp"

#include "buffer.hpp"

class FASTAReader
{
private:
    int16_t* seqs_;
    int maxSeqLength_;
    int numSeqs_;
    
    Buffer<int16_t> *seqs;
    
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
    bool accept(char ch);
    bool acceptSpace();
    bool acceptPrintable();
    bool acceptPrintableBut(char ch);
public:
    FASTAReader(char* filename);
    FASTAReader(std::string filename);
    //int16_t* getSeqs() { return this->seqs_; }
    bool next(Buffer<int16_t> *seqs, int factor);
    bool next(Buffer<int16_t> *seqs);
};
