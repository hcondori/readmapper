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

#include <cstdio>

#include "readmapper.hpp"

#include "backtrack.hpp"
#include "buffer.hpp"
#include "fastareader.hpp"
#include "sw.hpp"
#include "printer.hpp"


/*
* Wrapper for FASTAReader
*/
bool read_seqs(FASTAReader &reader1, FASTAReader &reader2, 
                    Buffer<int16_t> *seqs1, Buffer<int16_t> *seqs2,
                    int *seqs1_len, int *seqs2_len, 
                    std::vector<std::string> *seqs1_ids,
                    std::vector<std::string> *seqs2_ids)
{
    bool ans;
    
    #pragma omp critical
    ans = reader1.next(seqs1, seqs1_len, seqs1_ids) && 
        reader2.next(seqs2, seqs2_len, seqs2_ids);
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
    
    reader1.setDefault(0);
    reader2.setDefault(1);
    
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
        
        //containers for ids
        std::vector<std::string> seqs1_ids(VSIZE);
        std::vector<std::string> seqs2_ids(VSIZE);
        
        
        
        //legths of sequences
        int seqs1_len[VSIZE];
        int seqs2_len[VSIZE];
        
        //containter for flags
        Buffer<int8_t> flags(seq_len * seq_len * VSIZE, ALNSIZE);
        int16_t __attribute((aligned(ALNSIZE))) scores[VSIZE];
        int16_t __attribute((aligned(ALNSIZE))) ipos[VSIZE];
        int16_t __attribute((aligned(ALNSIZE))) jpos[VSIZE];
        
        //containers for arrays
        int16_t inf = gap_open + gap_extend + 1;
        //int16_t aF[256 * VSIZE] __attribute((aligned(ALNSIZE))) = {(int16_t)(-inf)};
        //int16_t aH[256 * VSIZE] __attribute((aligned(ALNSIZE))) = {0};

        int bsize = 128 * VSIZE;
	
        //Buffer<int16_t> E(bsize, ALNSIZE);
        Buffer<int16_t> F(bsize, ALNSIZE);
        Buffer<int16_t> H(bsize, ALNSIZE);
        //int16_t __attribute((aligned(ALNSIZE))) H[128 * VSIZE];
        
        //alignments
        char aln1[256];
        char aln2[256];
        
        //max sizes
        int max_x, max_y;
        
        //alignment start position
        int x0, y0;
        
        while(read_seqs(reader1, reader2, &seqs1, &seqs2, seqs1_len, seqs2_len, 
              &seqs1_ids, &seqs2_ids))
        {
            max_x = *std::max_element(seqs1_len, seqs1_len + VSIZE) + 1;
            max_y = *std::max_element(seqs2_len, seqs2_len + VSIZE) + 1;
            //E.clear(-inf);
            F.clear(-inf);
            H.clear(0);
            //flags.clear(0);
	    
            smith_waterman(seqs1.data(), seqs2.data(), match, mismatch, gap_open, gap_extend, 
			   flags.data(), scores, ipos, jpos, max_x, max_y, F.data(), H.data());
            
            for(int i = 0; i < VSIZE; i++)
            {
                //std::cout << scores[i] << std::endl;
                //std::cout << ipos[i] << std::endl;
                //std::cout << jpos[i] << std::endl;
                sw_backtrack(i, flags.data(), seqs1.data(), seqs2.data(), max_x, max_y,
                    aln1, aln2, ipos[i], jpos[i], x0, y0);
                //puts(aln1);
                //puts(aln2);
                print_alignment (stdout, seqs1_ids, seqs2_ids, scores, 
                    aln1, aln2, strlen(aln1), i);
            }
        }
    }
    return 0;
}
