/*
 * ForceData.h
 *
 *  Created on: Sept 30, 2020
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_FORCEDATA_H_
#define SRC_GROMACS_FORCEANAL_FORCEDATA_H_

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"
#include "InteractionType.h"
#include "AtomForce.h"

namespace ForceAnal {

class SummedMode
{
    public:

    SummedMode(int64_t grp1_sid, int64_t grp1_eid, int64_t grp2_sid, int64_t grp2_eid)
    {
        forces.reset(grp1_sid, grp1_eid);
        for (uint64_t i = 0; i < forces.length; ++i)
        {
            forces[i].reset(grp2_sid, grp2_eid);
        }
    }

    SummedForce forces;
};

class DetailedMode
{
    public:

    DetailedMode(int64_t grp1_sid, int64_t grp1_eid, int64_t grp2_sid, int64_t grp2_eid)
    {
        forces.reset(grp1_sid, grp1_eid);
        for (uint64_t i = 0; i < forces.length; ++i)
        {
            forces[i].reset(grp2_sid, grp2_eid);
            for (uint64_t j = 0; j < forces[i].length; ++j)
                forces[i][j].init();
        }
    }

    DetailedForce forces;
};

template <class ForceMode>
class ForceData : ForceMode
{
public:
    ForceData(int64_t group1_sid, int64_t group1_eid, int64_t group2_sid, int64_t group2_eid, real write_threshold, double average_factor)
     : ForceMode(group1_sid, group1_eid, group2_sid, group2_eid),
       avg_factor(average_factor)
    {
         threshold = std::max<real>(write_threshold * write_threshold, NONZERO_LIMIT);
    }

    void add_detailed_force(int affected, int applied, InteractionType itype, rvec force);

    void clear();

    void average_forces();

    void write_forces_txt(const std::string &fname, int32_t frameid);

    void write_forces_bin(const std::string &fname, int32_t frameid);

private:
    friend class ForceAnalysis;

    real threshold;

    double avg_factor;

};

}

#endif /* SRC_GROMACS_FORCEANAL_FORCEDATA_H_ */
