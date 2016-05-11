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

#include <iostream>

#include "readmapper.hpp"

#include "backtrack.hpp"

void
inplace_reverse (char* str, int len)
{
    char temp;
    for (int i = 0; i < len / 2; i++)
    {
        temp = str[i];
        str[i] = str[len - i - 1];
        str[len - i - 1] = temp;
    }
}

/*
 *
 */
void
sw_backtrack (int index, int16_t* flags, int16_t* seqs1, int16_t* seqs2, int w, int h,
              char* aln1, char* aln2, int x, int y, int& x0, int& y0)
{
    int16_t d_mask  = 0B0000000000001100; //diagonal
    int16_t bl_mask = 0B0000000000001000; //begin left
    int16_t bu_mask = 0B0000000000000100; //begin up
    int16_t cl_mask = 0B0000000000000010; //continue left
    int16_t cu_mask = 0B0000000000000001; //continue up

    int c = 0;

    while ((flags[(x * h + y) * VSIZE + index] & d_mask))
    {
        std::cout << "x=" << x << ", y=" << y 
            << ", flag=" << flags[(h*x + y)*VSIZE + index] << std::endl;
        if ((flags[(x * h + y) * VSIZE + index] & d_mask) == d_mask)
        {
            aln1[c] = seqs1[(--x) * VSIZE + index];
            aln2[c] = seqs2[(--y) * VSIZE + index];
            c++;
        }
        else if (flags[(x * h + y) * VSIZE + index] & bl_mask)
        {
            do
            {
                aln1[c] = seqs1[(x - 1) * VSIZE + index];
                aln2[c] = '-';
                c++;
            }
            while (flags[(x-- * h + y) * VSIZE + index] & cl_mask);
        }
        else
        {
            do
            {
                aln1[c] = '-';
                aln2[c] = seqs2[(y - 1) * VSIZE + index];
                c++;
            }
            while (flags[(x * h + y--) * VSIZE + index] & cu_mask);
        }
    }
    x0 = x + 1;
    y0 = y + 1;
    aln1[c] = '\0';
    aln2[c] = '\0';
    inplace_reverse(aln1, c);
    inplace_reverse(aln2, c);  
}

