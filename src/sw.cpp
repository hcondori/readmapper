/*
 Copyright (C) 2016 Héctor Condori Alagón.

 This file is part of readmapper, the massive Smith-Waterman pairwise sequence aligner.

 ALN is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ALN is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ALN.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <algorithm>
#include <cstdint>
#include <omp.h>

#include "readmapper.hpp"

/*
 * Smith-Waterman with match/mismatch values
 */
void
fill_table_16_to_16_i16(int8_t*  __restrict__ flags, int16_t*  __restrict__ seqs1, 
                        int16_t*  __restrict__ seqs2, int x, int y,
                        int16_t match, int16_t mismatch, const int16_t gap_open,
                        const int16_t gap_extend, int16_t*  __restrict__ scores,
                        int16_t*  __restrict__ ipos, int16_t*  __restrict__ jpos)
{
    const int16_t mask1 = 1;
    const int16_t mask2 = 2;
    const int16_t mask4 = 4;
    const int16_t mask8 = 8;

    int16_t __attribute((aligned(ALNSIZE))) s1[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) s2[VSIZE];

    int16_t __attribute((aligned(ALNSIZE))) E[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) E_sub[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) F[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) F_sub[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) diag[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_diag[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_left[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) score[VSIZE];
    
    int16_t __attribute((aligned(ALNSIZE))) c_up[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) c_left[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) b_up[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) b_left[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_gt_0[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_eq_diag[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_eq_E[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_eq_F[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_ne_E[VSIZE];
        
    int16_t __attribute((aligned(ALNSIZE))) flag[VSIZE];

    std::fill_n(flags, VSIZE * y, 0);
    
    int16_t __attribute((aligned(ALNSIZE))) is_zero[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) are_equal[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_gt_scores[VSIZE];

    std::fill_n(scores, VSIZE, 0);
    
    int inf = gap_open + gap_extend + 1;
    //int16_t aF[y * VSIZE] __attribute((aligned(ALNSIZE))) = {(int16_t)(-inf)};
    //int16_t aH[y * VSIZE] __attribute((aligned(ALNSIZE))) = {0};
    
    int16_t aF[256 * VSIZE] __attribute((aligned(ALNSIZE))) = {(int16_t)(-inf)};
    int16_t aH[256 * VSIZE] __attribute((aligned(ALNSIZE))) = {0};
    
    int16_t *sp1 = seqs1;
    int16_t *sp2;

    for(int i = 1; i < x; i++, sp1 += VSIZE)
    {
        sp2 = seqs2;
        std::copy_n(sp1, VSIZE, s1);
        std::fill_n(E, VSIZE, -inf);
        std::fill_n(H_diag, VSIZE, 0);
        std::fill_n(H, VSIZE, 0);
        std::fill_n(flags + VSIZE * i * y, VSIZE, 0);
        for (int j = 1; j < y; j++, sp2 += VSIZE)
        {
            std::copy_n(sp2, VSIZE, s2);
            std::copy_n(aH + VSIZE * j, VSIZE, H_left);
            std::copy_n(aF + VSIZE * j, VSIZE, F_sub);
            #pragma omp simd safelen(16) 
            for(int k = 0; k < VSIZE; k++)
            {
                is_zero[k] = s2[k] == 0;
                is_zero[k] = s1[k] == 0 && is_zero[k];                
                are_equal[k] = s1[k] == s2[k] && !is_zero[k];                               
                score[k] = are_equal[k]? match : mismatch;
                
                diag[k] = H_diag[k] + score[k];
                H_diag[k] = H_left[k];
                
                E_sub[k] = E[k] - gap_extend;        //for now, E is E_up
                E[k] = H[k] - gap_open;              //for now, H is H_up
                E[k] = std::max(E[k], E_sub[k]);      //TODO: may need declare simd

                F_sub[k] = F_sub[k] - gap_extend;
                F[k] = H_left[k] - gap_open;
                F[k] = std::max(F[k], F_sub[k]);
                
                H[k] = std::max(E[k], F[k]);
                H[k] = std::max(H[k], diag[k]);
                H[k] = std::max(H[k], (int16_t)0);
                
                
                //logic tests
                H_gt_0[k] = H[k] > 0;
                H_eq_E[k] = H[k] == E[k];
                H_eq_F[k] = H[k] == F[k];
                H_eq_diag[k] = H == diag;
                H_gt_scores[k] = H[k] > scores[k];
                
                // ********FLAGS********
                
                //c_up
                c_up[k] = E[k] == E_sub[k]; // E[i,j] == E[i,j-1]-gap_extent ?
                flag[k] = c_up[k] & mask1;

                //c_left
                c_left[k] = F[k] == F_sub[k]; // F[i,j] == F[i-1,j]-gap_extent ?
                flag[k] |= c_left[k] & mask2;
                
                //b_up
                b_up[k] = (H_eq_E[k] || H_eq_diag[k]) && H_gt_0[k];
                flag[k] |=  b_up[k] & mask4;
                
                //b_left
                b_left[k] = ((H_eq_E[k] && !H_eq_F[k]) || H_eq_diag[k]) && H_gt_0[k];
                flag[k] = b_left[k] & mask8;
                
                ipos[k] = H_gt_scores[k]? 0 : i;
                jpos[k] = H_gt_scores[k]? 0 : j;
                scores[k] = std::max(H[k], scores[k]);
            }
            std::copy_n(F, VSIZE, aF + VSIZE * j);
            std::copy_n(H, VSIZE, aH + VSIZE * j);
            std::copy_n(flag, VSIZE, flags + VSIZE * (y * i + j));
        }
    }
}


void smith_waterman(int16_t __restrict__ *seqs1, int16_t __restrict__  *seqs2, 
                    int match, int mismatch, int gap_open, int gap_extend, 
                    int8_t __restrict__ *flags, int16_t __restrict__  *scores, int16_t* 
                    __restrict__ ipos, int16_t * __restrict__ jpos)
{
    fill_table_16_to_16_i16(flags, seqs1, seqs2, 100, 100,
                        match, mismatch, gap_open, gap_extend, scores,
                        ipos, jpos);
    
}