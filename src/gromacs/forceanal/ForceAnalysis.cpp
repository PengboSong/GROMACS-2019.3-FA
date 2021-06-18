/*
 * ForceAnalysis.cpp
 *
 *  Created on: Sept 19, 2020
 *      Author: Pengbo Song
 */

#include <algorithm>

#include "gromacs/math/vectypes.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/extended_vec.h"
#include "gromacs/utility/real.h"
#include "ForceAnalysis.h"
#include "InteractionType.h"

ForceAnalysis::ForceAnalysis()
 : frame_count(0),
   summed_forces(threshold, 1.0 / Naverage),
   detailed_forces(threshold, 1.0 / Naverage)
{
}

ForceAnalysis::ForceAnalysis(int nfile, const t_filenm fnm[], gmx_mtop_t *mtop)
 : ForceAnal::ForceParaSet(nfile, fnm, mtop),
   frame_count(0),
   summed_forces(threshold, 1.0 / Naverage),
   detailed_forces(threshold, 1.0 / Naverage)
{
}

ForceAnalysis::~ForceAnalysis()
{    
}

void ForceAnalysis::add_pairforce(int i, int j, ForceAnal::InteractionType type, rvec f_ij)
{
    if (i < j)
    {
        if (summed_mode)
            summed_forces.add_detailed_force(i, j, type, f_ij);
        else
            detailed_forces.add_detailed_force(i, j, type, f_ij);
    }
    else
    {
        rvec f_ji;
        rvec_opp(f_ij, f_ji);
        if (summed_mode)
            summed_forces.add_detailed_force(j, i, type, f_ji);
        else
            detailed_forces.add_detailed_force(j, i, type, f_ji);
    }
}

void ForceAnalysis::add_nonbonded(int i, int j, real pf_coul, real pf_vdw, real dx, real dy, real dz)
{
    add_nonbonded_coulomb(i, j, pf_coul, dx, dy, dz);
    add_nonbonded_vdw(i, j, pf_vdw, dx, dy, dz);
}

void ForceAnalysis::add_nonbonded_coulomb(int i, int j, real pf_coul, real dx, real dy, real dz)
{
    rvec coul_force;
    coul_force[0] = dx * pf_coul;
    coul_force[1] = dy * pf_coul;
    coul_force[2] = dz * pf_coul;
    add_pairforce(i, j, ForceAnal::Interact_COULOMB, coul_force);
}

void ForceAnalysis::add_nonbonded_vdw(int i, int j, real pf_vdw, real dx, real dy, real dz)
{
    rvec lj_force;
    lj_force[0] = dx * pf_vdw;
    lj_force[1] = dy * pf_vdw;
    lj_force[2] = dz * pf_vdw;
    add_pairforce(i, j, ForceAnal::Interact_VDW, lj_force);
}

void ForceAnalysis::add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k, rvec r_ij, rvec r_kj, rvec r_ik)
{    
    rvec f_ij, f_ik, f_kj;
    svmul(iprod(f_i, r_ij) / norm2(r_ij), r_ij, f_ij);
    svmul(iprod(f_i, r_ik) / norm2(r_ik), r_ik, f_ik);
    svmul(iprod(f_k, r_kj) / norm2(r_kj), r_kj, f_kj);
    add_pairforce(ai, aj, ForceAnal::Interact_ANGLE, f_ij);
    add_pairforce(ai, ak, ForceAnal::Interact_ANGLE, f_ik);
    add_pairforce(ak, aj, ForceAnal::Interact_ANGLE, f_kj);
}

void ForceAnalysis::add_dihedral(int ai, int aj, int ak, int al, rvec f_i, rvec f_j, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl)
{
    rvec r_ik, r_il, r_jl;
    rvec f_ij, f_ik, f_il, f_jl, f_kj, f_kl;
    rvec_sub(r_ij, r_kj, r_ik);
    rvec_add(r_ik, r_kl, r_il);
    rvec_sub(r_kj, r_kl, r_jl);
    svmul(iprod(f_i, r_ij) / norm2(r_ij), r_ij, f_ij);
    svmul(iprod(f_i, r_ik) / norm2(r_ik), r_ik, f_ik);
    svmul(iprod(f_i, r_il) / norm2(r_il), r_il, f_il);
    svmul(iprod(f_j, r_jl) / norm2(r_jl), r_jl, f_jl);
    svmul(iprod(f_k, r_kj) / norm2(r_kj), r_kj, f_kj);
    svmul(iprod(f_k, r_kl) / norm2(r_kl), r_kl, f_kl);
    add_pairforce(ai, aj, ForceAnal::Interact_ANGLE, f_ij);
    add_pairforce(ai, ak, ForceAnal::Interact_ANGLE, f_ik);
    add_pairforce(ai, al, ForceAnal::Interact_ANGLE, f_il);
    add_pairforce(aj, al, ForceAnal::Interact_ANGLE, f_jl);
    add_pairforce(ak, aj, ForceAnal::Interact_ANGLE, f_kj);
    add_pairforce(ak, al, ForceAnal::Interact_ANGLE, f_kl);
}

void ForceAnalysis::write_frame(bool write_last_frame)
{
    ++frame_count;

    if ((frame_count % Naverage) == 0)
    {
        if (summed_mode)
        {
            summed_forces.average_forces();
            if (!res_txt_fn.empty())
                summed_forces.write_forces_txt(res_txt_fn, frame_count);
            if (!res_bin_fn.empty())
                summed_forces.write_forces_bin(res_bin_fn, frame_count);
            summed_forces.clear();
        }
        else
        {
            detailed_forces.average_forces();
            if (!res_txt_fn.empty())
                detailed_forces.write_forces_txt(res_txt_fn, frame_count);
            if (!res_bin_fn.empty())
                detailed_forces.write_forces_bin(res_bin_fn, frame_count);
            detailed_forces.clear();
        }
    }
    
    if (write_last_frame)
        frame_count = 0;   // Reset frame counter
}

void FA_add_nonbonded(class ForceAnalysis *FA, int i, int j, real pf_coul, real pf_vdw, real dx, real dy, real dz)
{
    FA->add_nonbonded(i, j, pf_coul, pf_vdw, dx, dy, dz);
}

void FA_add_nonbonded_coulomb(class ForceAnalysis *FA, int i, int j, real pf_coul, real dx, real dy, real dz)
{
    FA->add_nonbonded_coulomb(i, j, pf_coul, dx, dy, dz);
}

void FA_add_nonbonded_vdw(class ForceAnalysis *FA, int i, int j, real pf_vdw, real dx, real dy, real dz)
{
    FA->add_nonbonded_vdw(i, j, pf_vdw, dx, dy, dz);
}
