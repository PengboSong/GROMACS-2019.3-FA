/*
    ForceData.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/09/30
    Description: Lowlevel data container and analysis tool for Force Analysis.
*/

#ifndef SRC_GROMACS_FORCEANAL_FORCEDATA_H_
#define SRC_GROMACS_FORCEANAL_FORCEDATA_H_

#include <iostream>
#include <fstream>

#include "gromacs/utility/fatalerror.h"

#include "ForceAnalDef.h"
#include "ForceDataContainer.h"

namespace ForceAnal {

class SummedMode
{
public:
    SummedMode();

    SummedMode(const GrpIdx& grp1idx, const GrpIdx& grp2idx);

    void add_detailed_force(int affected, int applied, InteractionType itype, rvec force);

    void clear();

private:
    friend class SummedData;

    SummedForce forces;
};

class DetailedMode
{
public:
    DetailedMode();

    DetailedMode(const GrpIdx& grp1idx, const GrpIdx& grp2idx);

    void add_detailed_force(int affected, int applied, InteractionType itype, rvec force);

    void clear();

private:
    friend class DetailedData;

    DetailedForce forces;
};

class ListMode
{
public:
    ListMode();

    void add_detailed_force(int affected, int applied, InteractionType itype, rvec force);

    void clear();

private:
    friend class ListData;

    PairForceList forces;
};

class ForceData
{
public:
    ForceData();

    ForceData(real squared_threshold, real average_factor);

    void write_forces_txt(std::ofstream& txtstream);

    void write_forces_bin(std::ofstream& binstream, uint32_t& forces_count, const uint64_t& saddr, uint64_t& eaddr);

    real lowlim;

    real avgfactor;
};

class SummedData : public SummedMode, public ForceData
{
public:
    SummedData();

    SummedData(const GrpIdx& grp1idx, const GrpIdx& grp2idx, real squared_threshold, uint64_t Naverage);

    void average_forces();

    void write_forces_txt(std::ofstream& txtstream);

    void write_forces_bin(std::ofstream& binstream, uint32_t& forces_count, const uint64_t& saddr, uint64_t& eaddr);
};

class DetailedData : public DetailedMode, public ForceData
{
public:
    DetailedData();

    DetailedData(const GrpIdx& grp1idx, const GrpIdx& grp2idx, real squared_threshold, uint64_t Naverage);

    void average_forces();

    void write_forces_txt(std::ofstream& txtstream);

    void write_forces_bin(std::ofstream& binstream, uint32_t& forces_count, const uint64_t& saddr, uint64_t& eaddr);
};

class ListData : public ListMode, public ForceData
{
public:
    ListData();

    ListData(real squared_threshold);

    void write_forces_txt(std::ofstream& txtstream);

    void write_forces_bin(std::ofstream& binstream, uint32_t& forces_count, const uint64_t& saddr, uint64_t& eaddr);
};

}

#endif /* SRC_GROMACS_FORCEANAL_FORCEDATA_H_ */
