/*
    InteractionType.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/09/30
    Description: Force interaction type identifier.
*/

#include "InteractionType.h"

uint8_t ForceAnal::itype2index(ForceAnal::InteractionType itype)
{
    if (itype == 0) return 0;

    uint8_t n = 7;
    if (itype >> 4 == 0)
    {
        n -= 4;
        itype <<= 4;
    }
    if (itype >> 6 == 0)
    {
        n -= 2;
        itype <<= 2;
    }
    n += (itype >> 7);
    return n;
}