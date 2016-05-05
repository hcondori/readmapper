#include <algorithm>
#include <iostream>
#include <vector>
#include <omp.h>

#include "readmapper.hpp"

#include "buffer.hpp"
#include "sequence.hpp"
#include "fastareader.hpp"
#include "sw.hpp"

#define DEFAULT_SEQ_LEN 256

int main(int argc, char* argv[])
{
    if(argc < 3)
    {
        std::cout
            << "Error: no se especificaron suficientes archivos de entrada." 
            << std::endl;
        return 1;
    }
    
    std::string filename1 = argv[1];
    std::string filename2 = argv[2];
    
    FASTAReader reader1(filename1);
    FASTAReader reader2(filename2);
    
    //matriz de sustitucion
    int smatrix[]{ 5, -4, -4, -4,
                -4,  5, -4, -4,
                -4, -4,  5, -4,
                -4, -4, -4,  5};
    int gap_open   = 10;
    int gap_extend =  1;
    int match = 5;
    int mismatch = -4;
    
    #pragma omp parallel
    {
        int seq_len = DEFAULT_SEQ_LEN;
        
        //container vectors for sequences
        Buffer<int16_t> seqs1(seq_len * VSIZE, 64);
        Buffer<int16_t> seqs2(seq_len * VSIZE, 64);
        
        //legths of sequences
        int seqs1_len[VSIZE];
        int seqs2_len[VSIZE];
        
        //containter for flags
        Buffer<int8_t> flags(seq_len * seq_len * VSIZE, 64);
        int16_t scores[VSIZE];
        int16_t ipos[VSIZE];
        int16_t jpos[VSIZE];
        bool more_seqs = true;
        
        //read one batch
        #pragma omp critical
        more_seqs = reader1.next(&seqs1) && reader2.next(&seqs2);
        
        do
        { 
            //int x = reader1.getSequence().size();
            //int y = reader2.getSequence().size();
            
            smith_waterman(seqs1.data(), seqs2.data(), match, mismatch, gap_open, gap_extend, 
                           flags.data(), scores, ipos, jpos);
            //alignment.print();
            
            for(int i =0; i< VSIZE; i++)
                std::cout << scores[i] << std::endl;
            
            #pragma omp critical
            more_seqs = reader1.next(&seqs1) && reader2.next(&seqs2);
        }
        while(more_seqs);
    }
    return 0;
}
