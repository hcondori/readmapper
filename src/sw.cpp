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
#include <cstdint>
#include <omp.h>

#include "readmapper.hpp"

/*
 * Smith-Waterman with match/mismatch values
 */
void
fill_table_i16(int16_t*  __restrict__ flags, int16_t*  __restrict__ seqs1, 
               int16_t*  __restrict__ seqs2, const int x, const int y,
               const int16_t match, const int16_t mismatch, const int16_t gap_open,
               const int16_t gap_extend, int16_t*  __restrict__ scores,
               int16_t*  __restrict__ ipos, int16_t*  __restrict__ jpos,
               int16_t __restrict__ *aF, int16_t __restrict__ *aH)
{
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
    
    
    int16_t __attribute((aligned(ALNSIZE))) flag[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) c_up[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) c_left[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) b_up[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) b_left[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_gt_0[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_eq_diag[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_eq_E[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_eq_F[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_ne_E[VSIZE];    
    int16_t __attribute((aligned(ALNSIZE))) is_zero[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) are_equal[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_gt_scores[VSIZE];

    std::fill_n(flags, VSIZE * y, 0);
    std::fill_n(scores, VSIZE, 0);
    std::fill_n(ipos, VSIZE, 0);
    std::fill_n(jpos, VSIZE, 0);
    
    int16_t inf = gap_open + gap_extend + 1;
    //int16_t aF[y * VSIZE] __attribute((aligned(ALNSIZE))) = {(int16_t)(-inf)};
    //int16_t aH[y * VSIZE] __attribute((aligned(ALNSIZE))) = {0};
    
    //int16_t aF[256 * VSIZE] __attribute((aligned(ALNSIZE))) = {(int16_t)(-inf)};
    //int16_t aH[256 * VSIZE] __attribute((aligned(ALNSIZE))) = {0};
    
    int16_t *sp1 = seqs1;
    int16_t *sp2;

    #pragma unroll(4)
    for(int i = 1; i < x; i++, sp1 += VSIZE)
    {
        sp2 = seqs2;
        std::copy_n(sp1, VSIZE, s1);
        std::fill_n(E, VSIZE, -inf);
        std::fill_n(H_diag, VSIZE, 0);
        std::fill_n(H, VSIZE, 0);
        std::fill_n(flags + VSIZE * i * y, VSIZE, 0);
        #pragma unroll(4)
        for (int j = 1; j < y; j++, sp2 += VSIZE)
        {
            std::copy_n(sp2, VSIZE, s2);
            std::copy_n(aH + VSIZE * j, VSIZE, H_left);
            std::copy_n(aF + VSIZE * j, VSIZE, F_sub);
            #pragma omp simd aligned(ipos, jpos, scores: ALNSIZE)
            for(int k = 0; k < VSIZE; k++)
            {
                //is_zero[k] = s1[k] == 0;
                //are_equal[k] = !is_zero[k] & (s1[k] == s2[k]);
                are_equal[k] = s1[k] - s2[k];
                
                score[k] = are_equal[k]? mismatch : match;
                //score[k] = are_equal[k]? match : mismatch;
                
                diag[k] = H_diag[k] + score[k];
                H_diag[k] = H_left[k];              //for the next iteration
                
                E_sub[k] = E[k] - gap_extend;       //for now, E is E_up
                E[k] = H[k] - gap_open;             //for now, H is H_up
                E[k] = std::max(E[k], E_sub[k]);

                F_sub[k] = F_sub[k] - gap_extend;
                F[k] = H_left[k] - gap_open;
                F[k] = std::max(F[k], F_sub[k]);
                
                H[k] = std::max(E[k], F[k]);
                H[k] = std::max(H[k], diag[k]);
                H[k] = std::max(H[k], (int16_t)0);
                
                //logic tests
		//H_gt_0[k] = H[k] > 0;
                H_gt_0[k] = (-H[k] & (int16_t)(0X80000000)) >> 15;
                H_eq_E[k] = H[k] == E[k];
                H_eq_F[k] = H[k] == F[k];
                H_eq_diag[k] = H[k] == diag[k];
		//H_gt_scores[k] = H[k] > socres[k];
                H_gt_scores[k] = (scores[k] - H[k]) & (int16_t)(0X80000000) >> 15;
                
                // ********FLAGS********
                
                //c_up
                c_up[k] = E[k] == E_sub[k]; // E[i,j] == E[i,j-1]-gap_extent ?
                flag[k] = c_up[k];
                
                //c_left
                c_left[k] = F[k] == F_sub[k]; // F[i,j] == F[i-1,j]-gap_extent ?
                flag[k] |= c_left[k] << 1;
                
                //b_up
                b_up[k] = (H_eq_E[k] | H_eq_diag[k]) & H_gt_0[k];
                flag[k] |=  b_up[k] << 2;
                
                //b_left
                b_left[k] = ((H_eq_F[k] & ~H_eq_E[k]) | H_eq_diag[k]) & H_gt_0[k];
                flag[k] |= b_left[k] << 3;
                
                ipos[k] = H_gt_scores[k]? i : ipos[k];
                jpos[k] = H_gt_scores[k]? j : jpos[k];
                scores[k] = std::max(H[k], scores[k]);
            }
            std::copy_n(F, VSIZE, aF + VSIZE * j);
            std::copy_n(H, VSIZE, aH + VSIZE * j);
            std::copy_n(flag, VSIZE, flags + VSIZE * (y * i + j));
        }
    }
}


/*
 * Smith-Waterman with match/mismatch values with E, F, H tables
 */
void
fill_table_i16_2(int16_t*  __restrict__ flags, int16_t*  __restrict__ seqs1, 
		 int16_t*  __restrict__ seqs2, const int x, const int y,
		 const int16_t match, const int16_t mismatch, const int16_t gap_open,
		 const int16_t gap_extend, int16_t*  __restrict__ scores,
		 int16_t*  __restrict__ ipos, int16_t*  __restrict__ jpos,
		 int16_t __restrict__ *E, int16_t __restrict__ *F,
		 int16_t __restrict__ *H)
{
    int16_t __attribute((aligned(ALNSIZE))) s1[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) s2[VSIZE];

    int16_t __attribute((aligned(ALNSIZE))) score[VSIZE];
 
    int16_t __attribute((aligned(ALNSIZE))) is_zero[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) are_equal[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_gt_scores[VSIZE];
    
    std::fill_n(flags, VSIZE * y, 0);
    std::fill_n(scores, VSIZE, 0);
    std::fill_n(ipos, VSIZE, 0);
    std::fill_n(jpos, VSIZE, 0);
    
    int16_t inf = gap_open + gap_extend + 1;
    
    int16_t *sp1 = seqs1;
    int16_t *sp2;

    int current = VSIZE * (y + 1);
    int diagonal = 0;
    int up = VSIZE * y;
    int left = VSIZE;
    
    #pragma unroll(4)
    for(int i = 1; i < x; i++, sp1 += VSIZE)
    {
	current += VSIZE;
	diagonal += VSIZE;
	up += VSIZE;
	left += VSIZE;
	
        sp2 = seqs2;
        std::copy_n(sp1, VSIZE, s1);
        #pragma unroll(4)
        for (int j = 1; j < y; j++, sp2 += VSIZE)
        {
            std::copy_n(sp2, VSIZE, s2);
	    #pragma omp simd aligned(ipos, jpos, scores, E, F, H: ALNSIZE)
            for(int k = 0; k < VSIZE; k++)
            {
                //is_zero[k] = s1[k] == 0;
                //are_equal[k] = !is_zero[k] & (s1[k] == s2[k]);
                are_equal[k] = s1[k] - s2[k];
                
                score[k] = are_equal[k]? mismatch : match;
                
                E[current + k] = std::max(E[left + k] - gap_extend,
					  H[left + k] - gap_open);

                F[current + k] = std::max(F[up + k] - gap_extend,
					  H[up + k] - gap_open);
                
                H[current + k] = std::max(E[current + k], F[current + k]);
                H[current + k] = std::max(H[current + k], (int16_t)(H[diagonal + k] + score[k]));
                H[current + k] = std::max(H[current + k], (int16_t)0);
                
                H_gt_scores[k] = H[current + k] > scores[k];
                
                
                ipos[k] = H_gt_scores[k]? i : ipos[k];
                jpos[k] = H_gt_scores[k]? j : jpos[k];
                scores[k] = std::max(H[k], scores[k]);
            }
            current += VSIZE;
	    diagonal += VSIZE;
	    up += VSIZE;
	    left += VSIZE;
        }
    }
}


void smith_waterman(int16_t __restrict__ *seqs1, int16_t __restrict__  *seqs2, 
                    const int16_t match, const int16_t mismatch, 
                    const int16_t gap_open, const int gap_extend, 
                    int16_t __restrict__ *flags, int16_t __restrict__  *scores, 
                    int16_t* __restrict__ ipos, int16_t * __restrict__ jpos, 
                    const int x, const int y, 
                    int16_t __restrict__ *aF, int16_t __restrict__ *aH)
{
    fill_table_i16(flags, seqs1, seqs2, x, y,
                   match, mismatch, gap_open, gap_extend, scores,
                   ipos, jpos, aF, aH);
    
}

void smith_waterman_2(int16_t __restrict__ *seqs1, int16_t __restrict__  *seqs2, 
		      const int16_t match, const int16_t mismatch, 
		      const int16_t gap_open, const int gap_extend, 
		      int16_t __restrict__ *flags, int16_t __restrict__  *scores, 
		      int16_t* __restrict__ ipos, int16_t * __restrict__ jpos, 
		      const int x, const int y, 
		      int16_t __restrict__ *E, int16_t __restrict__ *F,
		      int16_t __restrict__ *H)
{
    fill_table_i16_2(flags, seqs1, seqs2, x, y,
                   match, mismatch, gap_open, gap_extend, scores,
		     ipos, jpos, E, F, H);
    
}

