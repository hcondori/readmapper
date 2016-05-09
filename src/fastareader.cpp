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
            current = fgetc_unlocked(file_);
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

bool FASTAReader::next(Buffer<int16_t> *seqs, int *seqs_len, int factor)
{
    this->seqs = seqs;
    this->vector_size = seqs->size();
    bool status = false;
    //this->seq_length = VSIZE / seq_count;
    std::fill(seq_sizes.begin(), seq_sizes.end(), 0);
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

bool FASTAReader::next(Buffer<int16_t> *seqs, int *seqs_len)
{
    return this->next(seqs, seqs_len, 1);
}
