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

template <class ForceMode>
ForceData<ForceMode>::ForceData(real write_threshold, double average_factor)
 : avg_factor(average_factor)
{
    threshold = std::max<real>(write_threshold * write_threshold, NONZERO_LIMIT);
}

template <class ForceMode>
ForceData<ForceMode>::~ForceData()
{
}

template <>
void ForceData<SummedMode>::add_detailed_force(int affected, int applied, InteractionType itype, rvec force)
{
    forces[affected].add(applied, itype, force);
}

template <>
void ForceData<DetailedMode>::add_detailed_force(int affected, int applied, InteractionType itype, rvec force)
{
    forces[affected][applied].add(itype, force);
}

template <>
void ForceData<SummedMode>::clear()
{
    for (uint64_t i = 0; i < forces.length; ++i)
        forces[i].init();
}

template <>
void ForceData<DetailedMode>::clear()
{
    for (uint64_t i = 0; i < forces.length; ++i)
        forces[i].init();
}

template <>
void ForceData<SummedMode>::average_forces()
{
    for (uint64_t i = 0; i < forces.length; ++i)
        forces[i] *= avg_factor;
}

template <>
void ForceData<DetailedMode>::average_forces()
{
    for (uint64_t i = 0; i < forces.length; ++i)
        for (uint64_t j = 0; j < forces[i].length; ++j)
            forces[i][j] *= avg_factor;
}

template <>
void ForceData<SummedMode>::write_forces_txt(std::string fname, int frameid)
{
    std::ofstream txtstream(fname, std::ios::app);

    if (!txtstream.is_open())
        gmx_fatal(FARGS, "GROMACS Force Analysis module can not write force data to file.\n");
    
    txtstream << "START FRAME " << frameid << std::endl;

    int ai, aj;
    int64_t idx;
    InteractionType itype;
    real fx, fy, fz, f;
    AtomForce force_ai;
    for (uint64_t i = 0; i < forces.length; ++i)
    {
        ai = i - forces.offset;
        force_ai = forces[i];
        for (uint64_t j = 0; j < force_ai.atomn; ++j)
        {
            aj = j - force_ai.offset;
            idx = 4 * j;
            itype = force_ai.itypes[j];
            fx = force_ai.forces[idx + XX];
            fy = force_ai.forces[idx + YY];
            fz = force_ai.forces[idx + ZZ];
            f = fx * fx + fy * fy + fz * fz;
            // Filter forces that are too small
            if (f > threshold)
            {
                txtstream << ai << ' ' << aj << ' ' << fx << ' ' << fy << ' ' << fz << ' ' << (int)itype << std::endl;
            }
        }
    }

    txtstream << "END FRAME " << frameid << std::endl;

    txtstream.close();
}

template <>
void ForceData<DetailedMode>::write_forces_txt(std::string fname, int frameid)
{
    std::ofstream txtstream(fname, std::ios::app);

    if (!txtstream.is_open())
        gmx_fatal(FARGS, "GROMACS Force Analysis module can not write force data to file.\n");
    
    txtstream << "START FRAME " << frameid << std::endl;

    int ai, aj;
    uint8_t idx;
    std::vector<real> force_aij;
    for (uint64_t i = 0; i < forces.length; ++i)
    {
        ai = i - forces.offset;
        for (uint64_t j = 0; j < forces[i].length; ++j)
        {
            aj = j - forces[i].offset;
            // Filter forces that are too small
            if (forces[i][j].sum() > threshold)
            {
                force_aij = forces[i][j].forces;
                txtstream << "Pair " << ai << ' ' << aj << std::endl;
                for (uint8_t k = 0; k < (Interact_COUNT + 1); ++k)
                {
                    idx = 4 * k;
                    txtstream << force_aij[idx + XX] << ' ' << force_aij[idx + YY] << ' ' << force_aij[idx + ZZ] << ' ' << force_aij[idx + XYZ] << std::endl; 
                }
            }
        }
    }

    txtstream << "END FRAME " << frameid << std::endl;

    txtstream.close();
}

template <>
void ForceData<SummedMode>::write_forces_bin(std::string fname, int32_t frameid)
{
    std::ofstream binstream(fname, std::ios::app | std::ios::binary);
    
    if (!binstream.is_open())
        gmx_fatal(FARGS, "GROMACS Force Analysis module can not write force data to file.\n");
    
    int32_t forces_count = 0;
    uint64_t initpos = binstream.tellp();
    binstream.write((char*)&frameid, sizeof(int32_t));
    binstream.write((char*)&forces_count, sizeof(int32_t));
    
    int32_t ai, aj;
    int64_t idx;
    InteractionType itype;
    real fx, fy, fz, f;
    AtomForce force_ai;
    for (uint64_t i = 0; i < forces.length; ++i)
    {
        ai = i - forces.offset;
        force_ai = forces[i];
        for (uint64_t j = 0; j < force_ai.atomn; ++j)
        {
            aj = j - force_ai.offset;
            idx = 4 * j;
            itype = force_ai.itypes[j];
            fx = force_ai.forces[idx + XX];
            fy = force_ai.forces[idx + YY];
            fz = force_ai.forces[idx + ZZ];
            f = fx * fx + fy * fy + fz * fz;
            // Filter forces that are too small
            if (f > threshold)
            {
                ++forces_count;
                binstream.write((char*)&ai, sizeof(int32_t));
                binstream.write((char*)&aj, sizeof(int32_t));
                binstream.write((char*)&itype, sizeof(InteractionType));
                binstream.write((char*)&fx, sizeof(real));
                binstream.write((char*)&fy, sizeof(real));
                binstream.write((char*)&fz, sizeof(real));
                binstream.write((char*)&f, sizeof(real));
            }
        }
    }

    // Get forces count and back to write at the start of the frame
    binstream.seekp(initpos + sizeof(int32_t));
    binstream.write((char*)&frameid, sizeof(int32_t));
    binstream.seekp(std::ios::end);

    binstream.close();
}

template <>
void ForceData<DetailedMode>::write_forces_bin(std::string fname, int frameid)
{
    std::ofstream binstream(fname, std::ios::app | std::ios::binary);

    if (!binstream.is_open())
        gmx_fatal(FARGS, "GROMACS Force Analysis module can not write force data to file.\n");
    
    binstream << "START FRAME " << frameid << std::endl;

    int32_t forces_count = 0;
    uint64_t initpos = binstream.tellp();
    binstream.write((char*)&frameid, sizeof(int32_t));
    binstream.write((char*)&forces_count, sizeof(int32_t));

    int ai, aj;
    std::vector<real> force_aij;
    for (uint64_t i = 0; i < forces.length; ++i)
    {
        ai = i - forces.offset;
        for (uint64_t j = 0; j < forces[i].length; ++j)
        {
            aj = j - forces[i].offset;
            // Filter forces that are too small
            if (forces[i][j].sum() > threshold)
            {
                ++forces_count;
                force_aij = forces[i][j].forces;
                binstream.write((char*)&ai, sizeof(int32_t));
                binstream.write((char*)&aj, sizeof(int32_t));
                binstream.write((char*)&force_aij, Interact_FORCEVEC_LEN * sizeof(real));
            }
        }
    }

    // Get forces count and back to write at the start of the frame
    binstream.seekp(initpos + sizeof(int32_t));
    binstream.write((char*)&frameid, sizeof(int32_t));
    binstream.seekp(std::ios::end);

    binstream.close();
}

}
