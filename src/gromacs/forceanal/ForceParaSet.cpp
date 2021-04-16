/*
 * ForceParaSet.cpp
 *
 *  Created on: Apr 14, 2021
 *      Author: Pengbo Song
 */

#include "ForceParaSet.h"

namespace ForceAnal {

ForceParaSet::ForceParaSet()
 : result_filename("result"),
   write_in_binary(true),
   Nevery(1),
   Nrepeat(1),
   Nfreq(1),
   group1_id(0),
   group2_id(0)
{
    reset_filename();
    check_average_steps();
}

ForceParaSet::~ForceParaSet()
{

}

void ForceParaSet::reset_filename()
{
    if (write_in_binary)
        result_full_filename = result_filename + ".fxt";
    else
        result_full_filename = result_filename + ".for";
}

void ForceParaSet::check_average_steps()
{
    uint64_t Nfreq_mult = uint64_t(Nevery) * uint64_t(Nrepeat);
    if (Nfreq < Nfreq_mult)
        Nfreq = Nfreq_mult;
}

}
