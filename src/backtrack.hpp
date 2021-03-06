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

#include <cstdio>
#include <cstdint>

void
sw_backtrack (int index, int8_t* flags, int16_t* seqs1, int16_t* seqs2, int w, int h,
              char* aln1, char* aln2, int x, int y, int& x0, int& y0);

