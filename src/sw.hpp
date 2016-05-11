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

#pragma once

#include "readmapper.hpp"

void smith_waterman(int16_t __restrict__ *seqs1, int16_t __restrict__  *seqs2, 
                    const int16_t match, const int16_t mismatch, 
                    const int16_t gap_open, const int gap_extend, 
                    int16_t __restrict__ *flags, int16_t __restrict__  *scores, 
                    int16_t* __restrict__ ipos, int16_t * __restrict__ jpos, 
                    const int x, const int y);
