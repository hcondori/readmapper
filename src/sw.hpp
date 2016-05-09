#pragma once

#include "readmapper.hpp"

void smith_waterman(int16_t __restrict__ *seqs1, int16_t __restrict__  *seqs2, 
                    int match, int mismatch, int gap_open, int gap_extend, 
                    int16_t __restrict__ *flags, int16_t __restrict__  *scores, int16_t* 
                    __restrict__ ipos, int16_t * __restrict__ jpos, int x, int y);