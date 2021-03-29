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

void ForceData::add_detailed_force(int affected, int applied, InteractionType itype, rvec force)
{
    detailed_forces[affected][applied].append(itype, force);
}

void ForceData::clear_detailed_forces()
{
    for (auto it = detailed_forces.begin(); it != detailed_forces.end(); ++it)
        for (auto iit = it->second.begin(); iit != it->second.end(); ++iit)
            iit->second.clear();
}

void ForceData::accumulate_summed_forces()
{
    InteractionType itype;
    rvec f;
    for (auto it = detailed_forces.begin(); it != detailed_forces.end(); ++it)
    {
        for (auto iit = it->second.begin(); iit != it->second.end(); ++iit)
        {
            iit->second.accumulate(&itype, &f);
            summed_forces[it->first][iit->first] += InteractionForce(itype, f);
            iit->second.clear();
        }
    }
}

void ForceData::average_summed_forces_laststep(double avg_factor)
{
    for (auto it = detailed_forces.begin(); it != detailed_forces.end(); ++it)
        for (auto iit = it->second.begin(); iit != it->second.end(); ++iit)
            summed_forces[it->first][iit->first] *= avg_factor;
}

uint32_t ForceData::pairforce_count()
{
    uint32_t count_num = 0;
    for (auto it = detailed_forces.begin(); it != detailed_forces.end(); ++it)
        count_num += it->second.size();
    return count_num;
}

void ForceData::write_average_summed_forces(double avg_factor)
{
    uint32_t forces_len = pairforce_count();
    if (result_file.is_open())
    {
        result_file.write((char*)&forces_len, sizeof(forces_len));
        int i, j;
        InteractionType itype;
        real fx, fy, fz;
        for (auto it = summed_forces.begin(); it != summed_forces.end(); ++it)
        {
            i = it->first;
            for (auto iit = it->second.begin(); iit != it->second.end(); ++iit)
            {
                j = iit->first;
                itype = iit->second.itype;
                fx = iit->second.f[0] * avg_factor;
                fy = iit->second.f[1] * avg_factor;
                fz = iit->second.f[2] * avg_factor;
                result_file.write((char*)&i, sizeof(int));
                result_file.write((char*)&j, sizeof(int));
                result_file.write((char*)&itype, sizeof(InteractionType));
                result_file.write((char*)&fx, sizeof(real));
                result_file.write((char*)&fy, sizeof(real));
                result_file.write((char*)&fz, sizeof(real));
                iit->second.clear();
            }
        }
    }
    else
    {
        std::runtime_error("Cannot open ForceAnalysis result file.");
        return ;
    }
}

}
