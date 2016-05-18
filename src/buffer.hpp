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

#include <algorithm>

template <typename T>
class Buffer
{
private:
    T *data_;
    int size_;
    int alignment_size_;
public:
    
    Buffer(int size, int alignment_size)
    {
        this->size_ = size;
        this->alignment_size_ = alignment_size;
        this->data_ = (T*)aligned_alloc(alignment_size, size * sizeof(T));
    }

    T* data()
    {
        return this->data_;
    }

    T& operator[](const int index)
    {
        return this->data_[index];
    }

    int size()
    {
        return this->size_;
    }

    int resize(int new_size)
    {
        //TODO: what if new size < old size ?
        this->size_ = new_size;
        T *temp = (T*)aligned_alloc(this->alignment_size_, new_size * sizeof(T));
        std::copy_n(this->data_, this->size_, temp);
        free(this->data_);
        this->data_ = temp;
        return 0;
    }
    
    void clear(T pattern)
    {
        std::fill_n(this->data_, this->size_, pattern);
    }
    
    void clear()
    {
        clear(0);
    }

    ~Buffer()
    {
        if(this->data_ != nullptr)
            free(this->data_);
    }
};
