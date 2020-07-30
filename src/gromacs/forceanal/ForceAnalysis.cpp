/*
 * ForceAnalysis.cpp
 *
 *  Created on: Sept 19, 2020
 *      Author: Pengbo Song
 */

#include "ForceAnalysis.h"

ForceAnalysis::ForceAnalysis()
 : forces(result_filename)
{
    
}

ForceAnalysis::~ForceAnalysis()
{
    
}

void ForceAnalysis::add_pairforce(int i, int j, ForceAnal::InteractionType type, rvec fi)
{
    rvec fj;
    fj[0] = -fi[0];
    fj[1] = -fi[1];
    fj[2] = -fi[2];
    forces.add_detailed_force(i, j, type, fi);
    forces.add_detailed_force(j, i, type, fj);
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

void ForceAnalysis::add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k)
{
    forces.add_detailed_force(ai, -1, ForceAnal::Interact_ANGLE, f_i);
    forces.add_detailed_force(aj, -1, ForceAnal::Interact_ANGLE, f_j);
    forces.add_detailed_force(ak, -1, ForceAnal::Interact_ANGLE, f_k);
}

void ForceAnalysis::add_dihedral(int i, int j, int k, int l, rvec f_i, rvec f_j, rvec f_k, rvec f_l)
{
    forces.add_detailed_force(i, -1, ForceAnal::Interact_DIHEDRAL, f_i);
    forces.add_detailed_force(j, -1, ForceAnal::Interact_DIHEDRAL, f_j);
    forces.add_detailed_force(k, -1, ForceAnal::Interact_DIHEDRAL, f_k);
    forces.add_detailed_force(l, -1, ForceAnal::Interact_DIHEDRAL, f_l);
}
