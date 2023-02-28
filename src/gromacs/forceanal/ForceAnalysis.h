/*
    ForceAnalysis.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/09/19
    Description: Core module for Force Analysis.
*/

#ifndef SRC_GROMACS_FORCEANAL_FORCEANALYSIS_H_
#define SRC_GROMACS_FORCEANAL_FORCEANALYSIS_H_

#define FORCEANAL_DEBUG

#include "gromacs/math/paddedvector.h"

#ifdef __cplusplus
#include <cstdint>
#include <cstdio>
#include <string>

#include "ForceAnalDef.h"
#include "ForceData.h"
#include "ForceParaSet.h"

class ForceAnalysis : public ForceAnal::ForceParaSet
{
public:
    ForceAnalysis();
    ForceAnalysis(int nfile, const t_filenm fnm[], gmx_mtop_t *mtop);
    ~ForceAnalysis();

    void init_outfiles();

    void next_frame() { ++frame_count; }

    bool in_grp1(const int idx);

    bool in_grp2(const int idx);

    bool in_grp(int& i, int& j, rvec f_ij);

    void add_pairforce(int i, int j, ForceAnal::InteractionType type, rvec f_ij);

    void add_nonbonded(int i, int j, real pf_coul, real pf_vdw, real dx, real dy, real dz);

    void add_nonbonded_coulomb(int i, int j, real pf_coul, real dx, real dy, real dz);

    void add_nonbonded_vdw(int i, int j, real pf_vdw, real dx, real dy, real dz);

    void tri_res(rvec f_i, rvec r_ij, rvec r_ik, rvec f_ij, rvec f_ik);

    void add_angle(int ai, int aj, int ak, rvec f_i, rvec gmx_unused f_j, rvec f_k, rvec r_ij, rvec r_kj, rvec r_ik);

    void add_dihedral(int ai, int aj, int ak, int al, rvec f_i, rvec f_j, rvec f_k, rvec gmx_unused f_l, rvec r_ij, rvec r_kj, rvec r_kl);

    void write_frame(bool write_last_frame = false);

    void write_forces();

    void write_atom_forces(const char* fnm, const rvec* f);

    void write_tot_forces(const rvec* f);

private:
    uint32_t frame_count;

    ForceAnal::SummedData summed_forces;

    ForceAnal::DetailedData detailed_forces;

    ForceAnal::ListData listed_forces;
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
