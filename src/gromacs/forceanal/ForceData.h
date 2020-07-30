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

    void add_detailed_force(int affected, int applied, InteractionType interact_type, rvec force);
    
    void write_detailed_forces();

private:
    friend class ForceAnalysis;

    std::ofstream result_file;

    std::vector<DetailedForce> detailed_forces;
};

}

#endif /* SRC_GROMACS_FORCEANAL_FORCEDATA_H_ */
