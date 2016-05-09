#include <algorithm>
#include <iostream>
#include <vector>
#include <omp.h>

#include "readmapper.hpp"

#include "backtrack.hpp"
#include "buffer.hpp"
#include "fastareader.hpp"
#include "sw.hpp"


/*
* Wrapper for FASTAReader
*/
bool read_seqs(FASTAReader &reader1, FASTAReader &reader2, 
                    Buffer<int16_t> *seqs1, Buffer<int16_t> *seqs2,
                    int *seqs1_len, int *seqs2_len)
{
    bool ans;
    #pragma omp critical
    ans = reader1.next(seqs1, seqs1_len) && reader2.next(seqs2, seqs2_len);
    return ans; 
}

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
        Buffer<int16_t> seqs1(seq_len * VSIZE, ALNSIZE);
        Buffer<int16_t> seqs2(seq_len * VSIZE, ALNSIZE);
        
        //legths of sequences
        int seqs1_len[VSIZE];
        int seqs2_len[VSIZE];
        
        //containter for flags
        Buffer<int16_t> flags(seq_len * seq_len * VSIZE, ALNSIZE);
        int16_t __attribute((aligned(ALNSIZE))) scores[VSIZE];
        int16_t __attribute((aligned(ALNSIZE))) ipos[VSIZE];
        int16_t __attribute((aligned(ALNSIZE))) jpos[VSIZE];
        bool more_seqs = true;
        
        //alignments
        char aln1[512];
        char aln2[512];
        
        //max sizes
        int x, y;
        
        //
        int x0, y0;
        
        //read one batch
        while(read_seqs(reader1, reader2, &seqs1, &seqs2, seqs1_len, seqs2_len))
        {
            std::cout << "Hola" << std::endl; 
            x = *std::max_element(seqs1_len, seqs1_len + VSIZE);
            y = *std::max_element(seqs2_len, seqs2_len + VSIZE);
            
            smith_waterman(seqs1.data(), seqs2.data(), match, mismatch, gap_open, gap_extend, 
                           flags.data(), scores, ipos, jpos, x, y);
            
            for(int i = 0; i < VSIZE; i++)
            {
                //sw_backtrack(i, flags.data(), seqs1.data(), seqs2.data(), x, y,
                    //aln1, aln2, ipos[i], jpos[i], x0, y0);
                
            }
            //std::cout << ipos[0] << std::endl;
        }
    }
    return 0;
}
