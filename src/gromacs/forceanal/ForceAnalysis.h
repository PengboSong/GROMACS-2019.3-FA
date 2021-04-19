/*
 * ForceAnalysis.h
 *
 *  Created on: Sept 19, 2020
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_FORCEANALYSIS_H_
#define SRC_GROMACS_FORCEANAL_FORCEANALYSIS_H_

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
#include <cstdint>
#include <cstdio>
#include <string>

#include "InteractionType.h"
#include "ForceData.h"
#include "ForceParaSet.h"

class ForceAnalysis : public ForceAnal::ForceParaSet
{
public:
    ForceAnalysis();
    ForceAnalysis(int nfile, const t_filenm fnm[], gmx_mtop_t *mtop);
    ~ForceAnalysis();

    void add_pairforce(int i, int j, ForceAnal::InteractionType type, rvec f_ij);

    void add_nonbonded(int i, int j, real pf_coul, real pf_vdw, real dx, real dy, real dz);

    void add_nonbonded_coulomb(int i, int j, real pf_coul, real dx, real dy, real dz);

    void add_nonbonded_vdw(int i, int j, real pf_vdw, real dx, real dy, real dz);

    void add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k, rvec r_ij, rvec r_kj, rvec r_ik);

    void add_dihedral(int ai, int aj, int ak, int al, rvec f_i, rvec f_j, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl);

    void set_average_steps(int Nevery_, int Nrepeat_, int Nfreq_);

    void write_frame(bool write_last_frame = false);

private:
    uint64_t frame_count;

    ForceAnal::ForceData forces;
};

#else

struct ForceAnalysis {};

#endif

#ifdef __cplusplus
extern "C" {
#endif

void FA_add_nonbonded(class ForceAnalysis *FA, int i, int j, real pf_coul, real pf_vdw, real dx, real dy, real dz);

void FA_add_nonbonded_coulomb(class ForceAnalysis *FA, int i, int j, real pf_coul, real dx, real dy, real dz);

void FA_add_nonbonded_vdw(class ForceAnalysis *FA, int i, int j, real pf_vdw, real dx, real dy, real dz);

#ifdef __cplusplus
}
#endif

#endif /* SRC_GROMACS_FORCEANAL_FORCEANALYSIS_H_ */
