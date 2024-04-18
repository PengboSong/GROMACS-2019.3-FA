/*
    ForceData.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/09/30
    Description: Lowlevel data container and analysis tool for Force Analysis.
*/

#include "ForceData.h"

namespace ForceAnal {

SummedMode::SummedMode()
{
}

SummedMode::SummedMode(const GrpIdx& grp1idx, const GrpIdx& grp2idx)
{
    forces.resize(grp1idx.first, grp1idx.second);
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
        forces[i].resize(grp2idx.first, grp2idx.second);
}

void SummedMode::add_detailed_force(int affected, int applied, InteractionType itype, rvec force)
{
    forces[affected].add(applied, itype, force);
}

void SummedMode::clear()
{
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
        forces[i].clear();
}

DetailedMode::DetailedMode()
{
}

DetailedMode::DetailedMode(const GrpIdx& grp1idx, const GrpIdx& grp2idx)
{
    forces.resize(grp1idx.first, grp1idx.second);
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
        forces[i].resize(grp2idx.first, grp2idx.second);
}

void DetailedMode::add_detailed_force(int affected, int applied, InteractionType itype, rvec force)
{
    forces[affected][applied].add(itype, force);
}

void DetailedMode::clear()
{
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
        for (uint64_t j = forces[i].id_begin(); j < forces[i].id_end(); ++j)
            forces[i][j].clear();
}

ListMode::ListMode()
{
}

void ListMode::add_detailed_force(int affected, int applied, InteractionType itype, rvec force)
{
    forces.push_back(PairForce(affected, applied, itype, force));
}

void ListMode::clear()
{
    forces.clear();
    std::vector<PairForce>().swap(forces);
}

ForceData::ForceData()
 : lowlim(NONZERO_LIMIT)
{
}

ForceData::ForceData(real squared_threshold)
 : lowlim(squared_threshold)
{
}

SummedData::SummedData()
{
}

SummedData::SummedData(const GrpIdx& grp1idx, const GrpIdx& grp2idx, real squared_threshold)
 : SummedMode(grp1idx, grp2idx), ForceData(squared_threshold)
{
}

void SummedData::export_forces(SummedData &other)
{
    InteractionType *itypes;
    real *fs;
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
    {
        itypes = forces[i].itypes_data();
        fs = forces[i].forces_data();
        for (uint64_t j = 0; j < forces[i].size(); ++j)
        {
            fs = forces[i].forces_data();
            rvec f_ij;
            f_ij[XX] = fs[4 * j + XX];
            f_ij[YY] = fs[4 * j + YY];
            f_ij[ZZ] = fs[4 * j + ZZ];
            other.add_detailed_force(i, j + forces[i].id_begin(), itypes[j], f_ij);
        }
    }
}

void SummedData::scale_forces(real scaling_factor)
{
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
        forces[i] *= scaling_factor;
}

void SummedData::write_forces_txt(std::ofstream& txtstream)
{
    if (!txtstream.is_open()) return;

    int ai, aj;
    int64_t idx;
    InteractionType itype, *itypes;
    real fx, fy, fz, f, *fs;
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
    {
        ai = static_cast<int>(i);
        itypes = forces[i].itypes_data();
        fs = forces[i].forces_data();
        for (uint64_t j = 0; j < forces[i].size(); ++j)
        {
            aj = j + forces[i].id_begin();
            idx = 4 * j;
            itype = itypes[j];
            fx = fs[idx + XX];
            fy = fs[idx + YY];
            fz = fs[idx + ZZ];
            f = fx * fx + fy * fy + fz * fz;
            // Filter forces that are too small
            if (f > lowlim)
            {
                txtstream << ai << ' ' << aj << ' ' << fx << ' ' << fy << ' ' << fz << ' ' << (int)itype << std::endl;
            }
        }
    }
}

void SummedData::write_forces_bin(std::ofstream& binstream, uint32_t& forces_count, const int64_t& saddr, int64_t& eaddr)
{
    if (!binstream.is_open()) return;
    binstream.seekp(saddr);

    int32_t ai, aj;
    int64_t idx;
    InteractionType itype, *itypes;
    real fx, fy, fz, f, *fs;
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
    {
        ai = static_cast<int>(i);
        itypes = forces[i].itypes_data();
        fs = forces[i].forces_data();
        for (uint64_t j = 0; j < forces[i].size(); ++j)
        {
            aj = j + forces[i].id_begin();
            idx = 4 * j;
            itype = itypes[j];
            fx = fs[idx + XX];
            fy = fs[idx + YY];
            fz = fs[idx + ZZ];
            f = fx * fx + fy * fy + fz * fz;
            // Filter forces that are too small
            if (f > lowlim)
            {
                f = std::sqrt(f);
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
    eaddr = binstream.tellp();
}

DetailedData::DetailedData()
{
}

DetailedData::DetailedData(const GrpIdx& grp1idx, const GrpIdx& grp2idx, real squared_threshold)
 : DetailedMode(grp1idx, grp2idx), ForceData(squared_threshold)
{
}

void DetailedData::scale_forces(real scaling_factor)
{
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
        for (uint64_t j = forces[i].id_begin(); j < forces[i].id_end(); ++j)
            forces[i][j] *= scaling_factor;
}

void DetailedData::write_forces_txt(std::ofstream& txtstream)
{
    if (!txtstream.is_open()) return;

    int ai, aj;
    uint8_t idx;
    real* force_aij;
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
    {
        ai = static_cast<int>(i);
        for (uint64_t j = forces[i].id_begin(); j < forces[i].id_end(); ++j)
        {
            aj = static_cast<int>(j);
            // Filter forces that are too small
            if (forces[i][j].sum() > lowlim)
            {
                force_aij = forces[i][j].data();
                txtstream << "Pair " << ai << ' ' << aj << std::endl;
                for (uint8_t k = 0; k < (Interact_COUNT + 1); ++k)
                {
                    idx = 4 * k;
                    txtstream << force_aij[idx + XX] << ' ' << force_aij[idx + YY] << ' ' << force_aij[idx + ZZ] << ' ' << force_aij[idx + XYZ] << std::endl; 
                }
            }
        }
    }
}

void DetailedData::write_forces_bin(std::ofstream& binstream, uint32_t& forces_count, const int64_t& saddr, int64_t& eaddr)
{
    if (!binstream.is_open()) return;
    binstream.seekp(saddr);

    int ai, aj;
    real* force_aij;
    for (uint64_t i = forces.id_begin(); i < forces.id_end(); ++i)
    {
        ai = static_cast<int>(i);
        for (uint64_t j = forces[i].id_begin(); j < forces[i].id_end(); ++j)
        {
            aj = static_cast<int>(j);
            // Filter forces that are too small
            if (forces[i][j].sum() > lowlim)
            {
                ++forces_count;
                binstream.write((char*)&ai, sizeof(int32_t));
                binstream.write((char*)&aj, sizeof(int32_t));
                binstream.write((char*)forces[i][j].data(), Interact_FORCEVEC_LEN * sizeof(real));
            }
        }
    }
    eaddr = binstream.tellp();
}

ListData::ListData()
{
}

ListData::ListData(real squared_threshold)
 : ListMode(), ForceData(squared_threshold)
{
}

void ListData::write_forces_txt(std::ofstream& txtstream)
{
    if (!txtstream.is_open()) return;

    for (const PairForce &pf : forces)
    {
        // Filter forces that are too small
        if (pf.f[XYZ] > lowlim)
        {
            txtstream << pf.i << ' ' << pf.j << ' ' << pf.f[XX] << ' ' << pf.f[YY] << ' ' << pf.f[ZZ] << ' ' << (int)pf.type << std::endl;
        }
    }
}

void ListData::write_forces_bin(std::ofstream& binstream, uint32_t& forces_count, const int64_t& saddr, int64_t& eaddr)
{
    if (!binstream.is_open()) return;
    binstream.seekp(saddr);

    for (const PairForce &pf : forces)
    {
        // Filter forces that are too small
        if (pf.f[XYZ] > lowlim)
        {
            ++forces_count;            
            binstream.write((char*)&pf.i, sizeof(int32_t));
            binstream.write((char*)&pf.j, sizeof(int32_t));
            binstream.write((char*)&pf.type, sizeof(InteractionType));
            binstream.write((char*)&pf.f, sizeof(real) * 4);
        }
    }
    eaddr = binstream.tellp();
}

}
