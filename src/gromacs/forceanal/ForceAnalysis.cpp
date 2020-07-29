/*
 * ForceAnalysis.cpp
 *
 *  Created on: Sept 19, 2020
 *      Author: Pengbo Song
 */

#include "ForceAnalysis.h"

ForceAnalysis::ForceAnalysis()
{
    result_file.open(result_filename, std::ios::binary);
    // Raise runtime error if file handler cannot be opened
    if (!result_file.is_open())
    {
        
    }
}

ForceAnalysis::~ForceAnalysis()
{
    result_file.close();
}

void ForceAnalysis::add_bond(int i, int j, rvec force)
{
    forces.add(i, InteractionType::Interact_BOND, force);
    forces.add(j, InteractionType::Interact_BOND, force);
}

void ForceAnalysis::add_polar_bond(int i, int j, rvec force)
{
    forces.add(i, InteractionType::Interact_BOND, force);
    forces.add(j, InteractionType::Interact_BOND, force);
}

void ForceAnalysis::add_14_interaction(int i, int j, rvec force)
{
    forces.add(i, InteractionType::Interact_BOND, force);
    forces.add(j, InteractionType::Interact_BOND, force);
}

void ForceAnalysis::add_nonbonded(int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz)
{
    rvec coul_force, lj_force;
    coul_force[0] = dx * pf_coul;
    coul_force[1] = dy * pf_coul;
    coul_force[2] = dz * pf_coul;
    lj_force[0] = dx * pf_coul;
    lj_force[1] = dy * pf_coul;
    lj_force[2] = dz * pf_coul;
    forces.add(i, InteractionType::Interact_COULOMB, coul_force);
    forces.add(j, InteractionType::Interact_LJ, lj_force);
}

void ForceAnalysis::add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k)
{
    forces.add(ai, InteractionType::Interact_ANGLE, f_i);
    forces.add(aj, InteractionType::Interact_ANGLE, f_j);
    forces.add(ak, InteractionType::Interact_ANGLE, f_k);
}

void ForceAnalysis::add_dihedral(int i, int j, int k, int l, rvec f_i, rvec f_j, rvec f_k, rvec f_l)
{
    forces.add(i, InteractionType::Interact_DIHEDRAL, f_i);
    forces.add(j, InteractionType::Interact_DIHEDRAL, f_j);
    forces.add(k, InteractionType::Interact_DIHEDRAL, f_k);
    forces.add(l, InteractionType::Interact_DIHEDRAL, f_l);
}

void ForceAnalysis::write_detailed_forces()
{
    int force_number = forces.detailed_forces.size();
    result_file.write((char*)&force_number, sizeof(force_number));
    int i, interact_type;
    real fx, fy, fz;
    for (auto it = forces.detailed_forces.begin(); it != forces.detailed_forces.end(); ++it)
    {
        i = it->ai;
        interact_type = (int)it->type;
        fx = it->f[0];
        fy = it->f[1];
        fz = it->f[2];
        result_file.write((char*)&i, sizeof(int));
        result_file.write((char*)&interact_type, sizeof(int));
        result_file.write((char*)&fx, sizeof(real));
        result_file.write((char*)&fy, sizeof(real));
        result_file.write((char*)&fz, sizeof(real));
    }
}
