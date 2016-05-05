#pragma once

template <class T>
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
        //TODO
    }

    ~Buffer()
    {
        if(this->data_ != nullptr)
            free(this->data_);
    }
};
