/*
    InteractionType.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/09/30
    Description: Force interaction type identifier.
*/

#ifndef SRC_GROMACS_FORCEANAL_INTERACTIONTYPE_H_
#define SRC_GROMACS_FORCEANAL_INTERACTIONTYPE_H_

// C++ STL
#include <cstdint>

namespace ForceAnal {
    using InteractionType = uint8_t;

    static const InteractionType Interact_NONE      =      0;
    static const InteractionType Interact_BOND      = 1 << 0;
    static const InteractionType Interact_POLAR     = 1 << 1;
    static const InteractionType Interact_ANGLE     = 1 << 2;
    static const InteractionType Interact_DIHEDRAL  = 1 << 3;
    static const InteractionType Interact_1_4       = 1 << 4;
    static const InteractionType Interact_COULOMB   = 1 << 5;
    static const InteractionType Interact_VDW       = 1 << 6;

    static const uint8_t Interact_COUNT = 7;

    static const InteractionType Interact_BONDED    = Interact_BOND + Interact_POLAR + Interact_ANGLE + Interact_DIHEDRAL;
    static const InteractionType Interact_NONBONDED = Interact_1_4 + Interact_COULOMB + Interact_VDW;
    static const InteractionType Interact_ALL       = Interact_BONDED + Interact_NONBONDED;

    uint8_t itype2index(InteractionType itype);
}

#endif /* SRC_GROMACS_FORCEANAL_INTERACTIONTYPE_H_ */
