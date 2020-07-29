/*
 * ForceAnalysis.h
 *
 *  Created on: Sept 19, 2020
 *      Author: Pengbo Song
 */

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

enum class InteractionType : std::int8_t {
    Interact_BOND,
    Interact_POLAR,
    Interact_ANGLE,
    Interact_DIHEDRAL,
    Interact_LJ,
    Interact_COULOMB,
    Interact_14,
    Interact_NONBONDED,
};

struct DetailedForce
{
    DetailedForce(int i, InteractionType interact_type, rvec force)
    {
        ai = i;
        type = interact_type;
        for (int j = 0; j < 3; ++j)
            f[j] = force[j];
    }

    ~DetailedForce()
    {}

    int ai;
    InteractionType type;
    rvec f;
};

struct AtomForce
{
    AtomForce()
    {}

    ~AtomForce()
    {}

    void add(int i, InteractionType interact_type, rvec force)
    {
        detailed_forces.push_back(DetailedForce(i, interact_type, force));
    }

    std::vector<DetailedForce> detailed_forces;
};

class ForceAnalysis
{
public:
    ForceAnalysis();

    ~ForceAnalysis();

    void add_bond(int i, int j, rvec force);

    void add_polar_bond(int i, int j, rvec force);

    void add_14_interaction(int i, int j, rvec force);

    void add_nonbonded(int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz);

    void add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k);

    void add_dihedral(int i, int j, int k, int l, rvec f_i, rvec f_j, rvec f_k, rvec f_l);

    void write_detailed_forces();

private:
    std::string result_filename = "result.bin";

    std::ofstream result_file;

    AtomForce forces;
}
