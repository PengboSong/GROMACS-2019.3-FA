/*
 * ForceData.cpp
 *
 *  Created on: Sept 30, 2020
 *      Author: Pengbo Song
 */

#include "ForceData.h"

namespace ForceAnal {

ForceData::ForceData(std::string const& result_filename)
{
    result_file.open(result_filename, std::ios::binary);
    // Raise runtime error if file handler cannot be opened
    if (!result_file.is_open())
    {
        
    }
}

ForceData::~ForceData()
{
    result_file.close();
}

void ForceData::add_detailed_force(int affected, int applied, InteractionType interact_type, rvec force)
{
    detailed_forces.push_back(
        DetailedForce(affected, applied, interact_type, force)
    );
}

void ForceData::write_detailed_forces()
{
    uint32_t forces_len = detailed_forces.size();
    if (result_file.is_open())
    {
        result_file.write((char*)&forces_len, sizeof(forces_len));
        real fx, fy, fz;
        for (auto it = detailed_forces.begin(); it != detailed_forces.end(); ++it)
        {
            fx = it->f[0];
            fy = it->f[1];
            fz = it->f[2];
            result_file.write((char*)&(it->i), sizeof(int));
            result_file.write((char*)&(it->j), sizeof(int));
            result_file.write((char*)&(it->type), sizeof(InteractionType));
            result_file.write((char*)&fx, sizeof(real));
            result_file.write((char*)&fy, sizeof(real));
            result_file.write((char*)&fz, sizeof(real));
        }
    }
    else
    {
        std::runtime_error("Cannot open ForceAnalysis result file.");
        return ;
    }
    
    detailed_forces.clear();
}

}
