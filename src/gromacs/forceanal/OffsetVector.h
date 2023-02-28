/*
    OffsetVector.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/07/28
    Description: Fixed length vector with non-zero begin.
*/

#ifndef SRC_GROMACS_FORCEANAL_OFFSETVECTOR_H_
#define SRC_GROMACS_FORCEANAL_OFFSETVECTOR_H_

// C++ STL
#include <cstdint>
#include <stdexcept>
#include <vector>

namespace ForceAnal
{

template<typename T>
class OffsetVector
{
using index = uint64_t;
using pointer = T*;
using const_pointer = const T*;
using reference = T&;
using const_reference = const T&;

public:
    OffsetVector() : offset(0), length(0) {}

    index size() { return length; }

    index id_begin() { return offset; }

    index id_end() { return offset + length; }

    pointer data() { return container.data(); }

    void clear()
    {
        container.assign(length, T());
    }

    void resize(index sloc, index eloc)
    {
        offset = sloc;
        length = eloc - sloc;
        container.resize(length, T());
    }

    reference operator[](const index idx)
    {
        if (idx < offset) throw std::out_of_range("Index out of bound of vector");
        index loc = idx - offset;
        if (loc >= length) throw std::out_of_range("Index out of bound of vector");
        return container[loc];
    }

private:
    index offset;

    index length;

    std::vector<T> container;
};

} // ForceAnal

#endif /* SRC_GROMACS_FORCEANAL_OFFSETVECTOR_H_ */