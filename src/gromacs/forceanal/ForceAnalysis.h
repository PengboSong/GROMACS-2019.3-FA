/*
 * ForceAnalysis.h
 *
 *  Created on: Sept 19, 2020
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_FORCEANALYSIS_H_
#define SRC_GROMACS_FORCEANAL_FORCEANALYSIS_H_

#include <cstdint>
#include <cstdio>
#include <string>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"
#include "InteractionType.h"
#include "ForceData.h"


class ForceAnalysis
{
public:
    ForceAnalysis();

    ~ForceAnalysis();

    void add_pairforce(int i, int j, ForceAnal::InteractionType type, rvec fi);

    void add_nonbonded(int i, int j, real pf_coul, real pf_vdw, real dx, real dy, real dz);

    void add_nonbonded_coulomb(int i, int j, real pf_coul, real dx, real dy, real dz);

    void add_nonbonded_vdw(int i, int j, real pf_vdw, real dx, real dy, real dz);

    void add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k);

    void add_dihedral(int i, int j, int k, int l, rvec f_i, rvec f_j, rvec f_k, rvec f_l);

private:
    std::string result_filename = "result.bin";

    ForceAnal::ForceData forces;
};

#endif /* SRC_GROMACS_FORCEANAL_FORCEANALYSIS_H_ */
