/*
 * ForceData.cpp
 *
 *  Created on: Sept 30, 2020
 *      Author: Pengbo Song
 */

#include "ForceData.h"

#include "gromacs/math/vec.h"
#include "gromacs/utility/fatalerror.h"

namespace ForceAnal {

ForceData::ForceData(bool summed_mode, real write_threshold, double average_factor)
 : summed(summed_mode),
   threshold(write_threshold),
   avg_factor(average_factor)
{
}

ForceData::~ForceData()
{
}

void ForceData::add_detailed_force(int affected, int applied, InteractionType itype, rvec force)
{
    detailed_forces[affected][applied].append(itype, force);
}

void ForceData::clear_detailed_forces()
{
    for (auto it = detailed_forces.begin(); it != detailed_forces.end(); ++it)
    {
        for (auto iit = it->second.begin(); iit != it->second.end(); ++iit)
            iit->second.clear();
        it->second.clear();
    }
}

void ForceData::clear_summed_forces()
{
    for (auto it = summed.begin(); it != summed.end(); ++it)
    {
        for (auto iit = it->second.begin(); iit != it->second.end(); ++iit)
            iit->second.clear();
        it->second.clear();
    }
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

void ForceData::average_summed_forces_laststep()
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

void ForceData::write_avg_forces(std::ofstream &res_stream, int frameid, InteractionType out_itype, bool write_bin)
{
    if (!res_stream.is_open())
        gmx_fatal(FARGS, "GROMACS Force Analysis module can not write force data to binary file.\n");
    
    int i, j, k;
    InteractionType itype;

    std::vector<PairForce> out_forces;
    if (summed)
    {
        for (auto it = summed_forces.begin(); it != summed_forces.end(); ++it)
        {
            i = it->first;
            for (auto iit = it->second.begin(); iit != it->second.end(); ++iit)
            {
                j = iit->first;
                PairForce pf(i, j, iit->second.itype, iit->second.f);
                pf *= avg_factor;
                if (pf.f < threshold)
                    continue;
                out_forces.push_back(pf);
            }
        }
    }
    else
    {
        for (auto it = detailed_forces.begin(); it != detailed_forces.end(); ++it)
        {
            i = it->first;
            for (auto iit = it->second.begin(); iit != it->second.end(); ++iit)
            {
                j = iit->first;
                for (k = 0; k < Interact_COUNT; ++k)
                {
                    if ((out_itype & 1 << k) != 0)
                    {
                        itype = 1 << k;
                        PairForce pf(i, j, itype, iit->second.f[k]);
                        pf *= avg_factor;
                        if (pf.f < threshold)
                            continue;
                        out_forces.push_back(pf);
                    }
                }
            }
        }
    }

    if (write_bin)
    {
        size_t forces_count = out_forces.size();
        res_stream.write((char*)&frameid, sizeof(int));
        res_stream.write((char*)&forces_count, sizeof(uint32_t));
        for (const PairForce &pf : out_forces)
        {
            res_stream.write((char*)&pf.i, sizeof(int));
            res_stream.write((char*)&pf.j, sizeof(int));
            res_stream.write((char*)&pf.type, sizeof(InteractionType));
            res_stream.write((char*)&pf.fx, sizeof(real));
            res_stream.write((char*)&pf.fy, sizeof(real));
            res_stream.write((char*)&pf.fz, sizeof(real));
            res_stream.write((char*)&pf.f, sizeof(real));
        }
    }
    else
    {
        res_stream << "START FRAME " << frameid << std::endl;
        for (const PairForce &pf : out_forces)
        {
            res_stream << pf.i << ' ' << pf.j << ' ' << pf.fx << ' ' << pf.fy << ' ' << pf.fz << ' ' << pf.type << std::endl;
        }
        res_stream << "END FRAME " << frameid << std::endl;
    }
}

}
