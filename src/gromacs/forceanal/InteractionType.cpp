/*
 * InteractionType.cpp
 *
 *  Created on: Sept 30, 2020
 *      Author: Pengbo Song
 */

#include "InteractionType.h"

uint8_t ForceAnal::index_itype(ForceAnal::InteractionType itype)
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