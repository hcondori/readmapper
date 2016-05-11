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
        
        //alignments
        char aln1[128];
        char aln2[128];
        
        //max sizes
        int max_x, max_y;
        
        //alignment start position
        int x0, y0;
        
        //read one batch
        while(read_seqs(reader1, reader2, &seqs1, &seqs2, seqs1_len, seqs2_len))
        {
            max_x = *std::max_element(seqs1_len, seqs1_len + VSIZE);
            max_y = *std::max_element(seqs2_len, seqs2_len + VSIZE);
            std::cout << "max x=" << max_x << " , max y=" << max_y << std::endl;

            smith_waterman(seqs1.data(), seqs2.data(), match, mismatch, gap_open, gap_extend, 
                           flags.data(), scores, ipos, jpos, max_x + 1, max_y + 1);
            
            for(int i = 0; i < VSIZE; i++)
            {
                sw_backtrack(i, flags.data(), seqs1.data(), seqs2.data(), max_x+1, max_y+1,
                    aln1, aln2, ipos[i], jpos[i], x0, y0);
                if(i==0)
                {    
                    std::cout << 
                        "flag=" << flags[((ipos[i] - 4)*(max_y+1) + jpos[i]-4)*VSIZE] << std::endl;
                    puts(aln1); puts(aln2);
                    std::cout << "x=" << ipos[i] << " , y=" << jpos[i] << std::endl;
                    std::cout << "x0=" << x0 << " , y0=" << y0 << std::endl;
                }
                
            }
            
        }
    }
    return 0;
}
