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
#include "DetailedForce.h"

namespace ForceAnal {

class ForceData
{
public:
    ForceData(bool summed_mode, real write_threshold, double average_factor);

    ~ForceData();

    void add_detailed_force(int affected, int applied, InteractionType itype, rvec force);

    void clear_detailed_forces();

    void clear_summed_forces();

    void clear();

    void accumulate_summed_forces();

    void average_summed_forces_laststep();

    uint32_t pairforce_count();

    void write_avg_forces(std::ofstream& res_stream, int frameid, InteractionType out_itype, bool write_bin);

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
