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

    SummedMode()
    {}

    SummedForce forces;
};

class DetailedMode
{
    public:

    DetailedMode()
    {}

    DetailedForce forces;
};

template <class ForceMode>
class ForceData : ForceMode
{
public:
    ForceData(real write_threshold, double average_factor);

    void add_detailed_force(int affected, int applied, InteractionType itype, rvec force);

    void clear();

    void average_forces();

    void write_forces_txt(std::string fname, int32_t frameid);

    void write_forces_bin(std::string fname, int32_t frameid);

private:
    friend class ForceAnalysis;

    real threshold;

    double avg_factor;

    bool summed;

    std::map<int, std::map<int, DetailedForce>> detailed_forces;

    std::map<int, std::map<int, InteractionForce>> summed_forces;

};

}

#endif /* SRC_GROMACS_FORCEANAL_FORCEDATA_H_ */
