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
    ForceData(std::string const& result_filename);

    ~ForceData();

    void add_detailed_force(int affected, int applied, InteractionType itype, rvec force);

    void clear_detailed_forces();

    void accumulate_summed_forces();

    void average_summed_forces_laststep(double avg_factor = 1.0);

    uint32_t pairforce_count();
    
    void write_average_summed_forces(double avg_factor = 1.0);

private:
    friend class ForceAnalysis;

    std::ofstream result_file;

    std::map<int, std::map<int, DetailedForce>> detailed_forces;

    std::map<int, std::map<int, InteractionForce>> summed_forces;

    bool accumulated_mode = true;
};

}

#endif /* SRC_GROMACS_FORCEANAL_FORCEDATA_H_ */
