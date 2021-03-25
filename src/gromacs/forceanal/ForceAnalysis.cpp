/*
 * ForceAnalysis.cpp
 *
 *  Created on: Sept 19, 2020
 *      Author: Pengbo Song
 */

#include "gromacs/math/vectypes.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/extended_vec.h"
#include "gromacs/utility/real.h"
#include "ForceAnalysis.h"
#include "InteractionType.h"

ForceAnalysis::ForceAnalysis()
 : forces(result_filename)
{
    
}

ForceAnalysis::~ForceAnalysis()
{
    
}

void ForceAnalysis::add_pairforce(int i, int j, ForceAnal::InteractionType type, rvec f_ij)
{
    rvec f_ji;
    rvec_opp(f_ij, f_ji);
    forces.add_detailed_force(i, j, type, f_ij);
    forces.add_detailed_force(j, i, type, f_ji);
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
    rvec r_ji, r_jk, r_ki, f_ij, f_ik, f_ji, f_jk, f_ki, f_kj;
    rvec_opp(r_ij, r_ji);
    rvec_opp(r_ik, r_ki);
    rvec_opp(r_kj, r_jk);
    svmul(iprod(f_i, r_ij) / norm2(r_ij), r_ij, f_ij);
    svmul(iprod(f_i, r_ik) / norm2(r_ik), r_ik, f_ik);
    svmul(iprod(f_j, r_ji) / norm2(r_ji), r_ji, f_ji);
    svmul(iprod(f_j, r_jk) / norm2(r_jk), r_jk, f_jk);
    svmul(iprod(f_k, r_ki) / norm2(r_ki), r_ki, f_ki);
    svmul(iprod(f_k, r_kj) / norm2(r_kj), r_kj, f_kj);
    forces.add_detailed_force(ai, aj, ForceAnal::Interact_ANGLE, f_ij);
    forces.add_detailed_force(ai, ak, ForceAnal::Interact_ANGLE, f_ik);
    forces.add_detailed_force(aj, ai, ForceAnal::Interact_ANGLE, f_ji);
    forces.add_detailed_force(aj, ak, ForceAnal::Interact_ANGLE, f_jk);
    forces.add_detailed_force(ak, ai, ForceAnal::Interact_ANGLE, f_ki);
    forces.add_detailed_force(ak, aj, ForceAnal::Interact_ANGLE, f_kj);
}

void ForceAnalysis::add_dihedral(int ai, int aj, int ak, int al, rvec f_i, rvec f_j, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl)
{
    rvec r_ik, r_il, r_ji, r_jk, r_jl, r_ki, r_li, r_lj, r_lk;
    rvec f_ij, f_ik, f_il, f_ji, f_jk, f_jl, f_ki, f_kj, f_kl, f_li, f_lj, f_lk;
    rvec_opp(r_ij, r_ji);
    rvec_opp(r_kj, r_jk);
    rvec_opp(r_kl, r_lk);
    rvec_sub(r_ij, r_kj, r_ik);
    rvec_sub(r_ik, r_lk, r_il);
    rvec_sub(r_jk, r_lk, r_jl);
    rvec_opp(r_ik, r_ki);
    rvec_opp(r_il, r_li);
    rvec_opp(r_jl, r_lj);
    svmul(iprod(f_i, r_ij) / norm2(r_ij), r_ij, f_ij);
    svmul(iprod(f_i, r_ik) / norm2(r_ik), r_ik, f_ik);
    svmul(iprod(f_i, r_il) / norm2(r_il), r_il, f_il);
    svmul(iprod(f_j, r_ji) / norm2(r_ji), r_ji, f_ji);
    svmul(iprod(f_j, r_jk) / norm2(r_jk), r_jk, f_jk);
    svmul(iprod(f_j, r_jl) / norm2(r_jl), r_jl, f_jl);
    svmul(iprod(f_k, r_ki) / norm2(r_ki), r_ki, f_ki);
    svmul(iprod(f_k, r_kj) / norm2(r_kj), r_kj, f_kj);
    svmul(iprod(f_k, r_kl) / norm2(r_kl), r_kl, f_kl);
    svmul(iprod(f_l, r_li) / norm2(r_li), r_li, f_li);
    svmul(iprod(f_l, r_lj) / norm2(r_lj), r_lj, f_lj);
    svmul(iprod(f_l, r_lk) / norm2(r_lk), r_lk, f_lk);
    forces.add_detailed_force(ai, aj, ForceAnal::Interact_ANGLE, f_ij);
    forces.add_detailed_force(ai, ak, ForceAnal::Interact_ANGLE, f_ik);
    forces.add_detailed_force(ai, al, ForceAnal::Interact_ANGLE, f_il);
    forces.add_detailed_force(aj, ai, ForceAnal::Interact_ANGLE, f_ji);
    forces.add_detailed_force(aj, ak, ForceAnal::Interact_ANGLE, f_jk);
    forces.add_detailed_force(aj, al, ForceAnal::Interact_ANGLE, f_jl);
    forces.add_detailed_force(ak, ai, ForceAnal::Interact_ANGLE, f_ki);
    forces.add_detailed_force(ak, aj, ForceAnal::Interact_ANGLE, f_kj);
    forces.add_detailed_force(ak, al, ForceAnal::Interact_ANGLE, f_kl);
    forces.add_detailed_force(al, ai, ForceAnal::Interact_ANGLE, f_li);
    forces.add_detailed_force(al, aj, ForceAnal::Interact_ANGLE, f_lj);
    forces.add_detailed_force(al, ak, ForceAnal::Interact_ANGLE, f_lk);
}

void ForceAnalysis::write_frame()
{
    forces.write_detailed_forces();
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
