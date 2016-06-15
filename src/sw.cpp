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
fill_table_i16(int8_t*  __restrict__ flags, int16_t*  __restrict__ seqs1, 
               int16_t*  __restrict__ seqs2, const int x, const int y,
               const int16_t match, const int16_t mismatch, const int16_t gap_open,
               const int16_t gap_extend, int16_t*  __restrict__ scores,
               int16_t*  __restrict__ ipos, int16_t*  __restrict__ jpos,
               int16_t __restrict__ *aF, int16_t __restrict__ *aH)
{
    int16_t __attribute((aligned(ALNSIZE))) s1[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) s2[VSIZE];

    int16_t __attribute((aligned(ALNSIZE))) E[VSIZE];
    //int16_t __attribute((aligned(ALNSIZE))) E_sub[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) F[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) F_sub[VSIZE];
    //int16_t __attribute((aligned(ALNSIZE))) diag[VSIZE];
    //int16_t __attribute((aligned(ALNSIZE))) H[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_diag[VSIZE];
    int16_t __attribute((aligned(ALNSIZE))) H_left[VSIZE];
    //int16_t __attribute((aligned(ALNSIZE))) score[VSIZE];

    int16_t E_sub;
    //int16_t F_sub;
    int16_t diag;
    int16_t H;
    int16_t score;
    

    int8_t  flag;
    int16_t  c_up;
    int16_t  c_left;
    int16_t  b_up;
    int16_t  b_left;
    int16_t  H_gt_0;
    int16_t  H_eq_diag;
    int16_t  H_eq_E;
    int16_t  H_eq_F;
    int16_t  H_ne_E;    
    int16_t  is_zero;
    int16_t  are_equal;
    int16_t H_gt_scores;

    int8_t __attribute((aligned(ALNSIZE))) flag_tmp[VSIZE];


    std::fill_n(flags, VSIZE * y, 0);
    std::fill_n(scores, VSIZE, 0);
    std::fill_n(ipos, VSIZE, 0);
    std::fill_n(jpos, VSIZE, 0);
    
    int16_t inf = gap_open + gap_extend + 1;
    
    //int16_t aF[256 * VSIZE] __attribute((aligned(ALNSIZE))) = {(int16_t)(-inf)};
    //int16_t aH[256 * VSIZE] __attribute((aligned(ALNSIZE))) = {0};
    
    int16_t *sp1 = seqs1;
    int16_t *sp2;

    //#pragma unroll(4)
    for(int i = 1; i < x; i++, sp1 += VSIZE)
    {
        sp2 = seqs2;
        std::copy_n(sp1, VSIZE, s1);
        std::fill_n(E, VSIZE, -inf);
        std::fill_n(H_diag, VSIZE, 0);
        //std::fill_n(H, VSIZE, 0);
	H = 0;
        std::fill_n(flags + VSIZE * i * y, VSIZE, 0);
        //#pragma unroll(4)
        for (int j = 1; j < y; j++, sp2 += VSIZE)
        {
            std::copy_n(sp2, VSIZE, s2);
            std::copy_n(aH + VSIZE * j, VSIZE, H_left);
            std::copy_n(aF + VSIZE * j, VSIZE, F_sub);
	    //F_sub = 0;
            #pragma simd aligned(ipos, jpos, scores: ALNSIZE)    \
		private(diag, score, E_sub,				\
			flag, c_up, c_left, b_up, b_left, H_gt_0, H_eq_diag, \
			H_eq_E, H_eq_F, H_ne_E, is_zero, are_equal) firstprivate(H)
            for(int k = 0; k < VSIZE; k++)
            {
                //is_zero[k] = s1[k] == 0;
                //are_equal[k] = !is_zero[k] & (s1[k] == s2[k]);
                are_equal = s1[k] - s2[k];
                
                score = are_equal? mismatch : match;
                
                diag = H_diag[k] + score;
                H_diag[k] = H_left[k];              //for the next iteration
                
                E_sub = E[k] - gap_extend;       //for now, E is E_up
                E[k] = H - gap_open;             //for now, H is H_up
                E[k] = std::max(E[k], E_sub);

                F_sub[k] = F_sub[k] - gap_extend;
                F[k] = H_left[k] - gap_open;
                F[k] = std::max(F[k], F_sub[k]);
                 
                H = std::max(E[k], F[k]);
                H = std::max(H, diag);
                H = std::max(H, (int16_t)0);
                
                //logic tests
		H_gt_0 = H > 0;
                //H_gt_0 = (-H & (int16_t)(0X8000)) >> 15;
                H_eq_E = H == E[k];
                H_eq_F = H == F[k];
                H_eq_diag = H == diag;
		H_gt_scores = H > scores[k];
                //H_gt_scores = (scores[k] - H) & (int16_t)(0X8000);
                
                // ********FLAGS********
                
                //c_up
                c_up = E[k] == E_sub; // E[i,j] == E[i,j-1]-gap_extent ?
                flag = c_up;
                
                //c_left
                c_left = F[k] == F_sub[k]; // F[i,j] == F[i-1,j]-gap_extent ?
                flag |= c_left << 1;
                
                //b_up
                b_up = (H_eq_E | H_eq_diag) & H_gt_0;
                flag |=  b_up << 2;
                
                //b_left
                b_left = ((H_eq_F & ~H_eq_E) | H_eq_diag) & H_gt_0;
                flag |= b_left << 3;
                
                ipos[k] = H_gt_scores? i : ipos[k];
                jpos[k] = H_gt_scores? j : jpos[k];
                scores[k] = std::max(H, scores[k]);

		//flag_tmp[k] = flag;
		
		flags[VSIZE * (y * i + j) + k] = flag;
		aH[VSIZE * j + k] = H;
            }
	    std::copy_n(F, VSIZE, aF + VSIZE * j);
            //std::copy_n(H, VSIZE, aH + VSIZE * j);
            //std::copy_n(flag_tmp, VSIZE, flags + VSIZE * (y * i + j));
        }
    }
}

void
fill_table_i16_2(int8_t*  __restrict__ flags, const int16_t*  __restrict__ seqs1, 
               const int16_t*  __restrict__  seqs2, const int x, const int y,
               const int16_t match, const int16_t mismatch, const int16_t gap_open,
               const int16_t gap_extend, int16_t*  __restrict__ scores,
               int16_t*  __restrict__ ipos, int16_t*  __restrict__ jpos,
               int16_t __restrict__ *aF2, int16_t __restrict__  *aH2)
{
    int16_t s1, s2;
    int16_t diag;
    int16_t E;
    int16_t E_sub;
    int16_t F;
    int16_t F_sub;
    int16_t H_diag;
    int16_t H_left;
    int16_t H;
    int16_t score;

    int16_t flag;
    int16_t c_up;
    int16_t c_left;
    int16_t b_up;
    int16_t b_left;
    int16_t H_gt_0;
    int16_t H_eq_diag;
    int16_t H_eq_E;
    int16_t H_eq_F;
    int16_t H_ne_E;    
    int16_t is_zero;
    int16_t H_gt_scores;

    /*bool flag;
      bool c_up;
      bool c_left;
      bool b_up;
      bool b_left;
      bool H_gt_0;
      bool H_eq_diag;
      bool H_eq_E;
      bool H_eq_F;
      bool H_ne_E;    
      bool is_zero;
      bool H_gt_scores;*/

    const int16_t inf = gap_open + gap_extend + 1;

    int8_t __attribute((aligned(ALNSIZE))) flag_tmp[VSIZE];


    std::fill_n(flags, VSIZE * y, 0);
    //std::fill_n(scores, VSIZE, 0);
    std::fill_n(ipos, VSIZE, 0);
    std::fill_n(jpos, VSIZE, 0);
    
    int16_t aF[128 * VSIZE] __attribute((aligned(ALNSIZE))) = {(int16_t)(-inf)};
    int16_t aH[128 * VSIZE] __attribute((aligned(ALNSIZE))) = {0};
    std::fill_n(aF, VSIZE * 128, -inf);
    std::fill_n(aH, VSIZE * 128, 0);
    
    int16_t temp_i, temp_j, temp_score;

    for(int i = 1; i < x; i++)
    	std::fill_n(flags + VSIZE * i * y, VSIZE, 0);

    int i, j;

    #pragma omp simd aligned(ipos, jpos, scores, aF: ALNSIZE)		\
	private(s1, s2, diag, score, E_sub, E, F_sub, F, H_left, H_diag, H, \
		flag, c_up, c_left, b_up, b_left, H_gt_0, H_eq_diag,	\
		H_eq_E, H_eq_F, H_ne_E, is_zero, H_gt_scores,		\
		temp_i, temp_j, temp_score, i, j)
    for(int k = 0; k < VSIZE; k++)
    {
	temp_score = 0;
	temp_i = 0;
	temp_j = 0;
	for(i = 1; i < x; i++)
	{
	    for (j = 1, H = 0, H_diag = 0, E = -inf; j < y; j++)
	    {
		s1 = seqs1[VSIZE * (i - 1) + k];
		s2 = seqs2[VSIZE * (j - 1) + k];
		F_sub  = aF[VSIZE * j + k];
		H_left = aH[VSIZE * j + k];
		//H_left = aH[j];
		
		//is_zero[k] = s1[k] == 0;
		//are_equal[k] = !is_zero[k] & (s1[k] == s2[k]);
                
		score = (s1 == s2)? match : mismatch;
                
		diag = H_diag + score;
		H_diag = H_left;              //for the next iteration
                
		E_sub = E - gap_extend;       //for now, E is E_up
		E = H - gap_open;             //for now, H is H_up
		E = std::max(E, E_sub);

		F_sub = F_sub - gap_extend;
		F = H_left - gap_open;
		F = std::max(F, F_sub);
                 
		H = std::max(E, F);
		H = std::max(H, diag);
		H = std::max(H, (int16_t)0);
                
		//logic tests
		H_gt_0 = H > 0;
		H_eq_E = H == E;
		H_eq_F = H == F;
		H_eq_diag = H == diag;
		H_gt_scores = H > temp_score;
                
		// ********FLAGS********
                
		//c_up
		c_up = E == E_sub; // E[i,j] == E[i,j-1]-gap_extent ?
		flag = c_up;
                
		//c_left
		c_left = F == F_sub; // F[i,j] == F[i-1,j]-gap_extent ?
		flag |= c_left << 1;
			
		//b_up
		b_up = (H_eq_E | H_eq_diag) & H_gt_0;
		flag |=  b_up << 2;
                
		//b_left
		b_left = ((H_eq_F & ~H_eq_E) | H_eq_diag) & H_gt_0;
		flag |= b_left << 3;
                
		temp_i = H_gt_scores? i : temp_i;
		temp_j = H_gt_scores? j : temp_j;
		temp_score = std::max(H, temp_score);
		
		flags[VSIZE * (y * i + j) + k] = flag;
		aF[VSIZE * j + k] = F;
		aH[VSIZE * j + k] = H;
		
		//aF[k] = F;
		//aH[k] = H;
		
		//aH[j] = H;
	    }
	}
	scores[k] = temp_score;
	ipos[k] = temp_i;
	jpos[k] = temp_j;
    }
}



/*
 * Smith-Waterman with match/mismatch values with E, F, H tables
 */
void
fill_table_i16_3(int16_t*  __restrict__ flags, int16_t*  __restrict__ seqs1, 
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
                    int8_t __restrict__ *flags, int16_t __restrict__  *scores, 
                    int16_t __restrict__ *ipos, int16_t __restrict__ *jpos, 
                    const int x, const int y, 
                    int16_t __restrict__ *aF, int16_t __restrict__ *aH)
{
    fill_table_i16_2(flags, seqs1, seqs2, x, y,
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
    fill_table_i16_3(flags, seqs1, seqs2, x, y,
                   match, mismatch, gap_open, gap_extend, scores,
		     ipos, jpos, E, F, H);
    
}

