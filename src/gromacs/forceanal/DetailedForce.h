/*
 * DetailedForce.h
 *
 *  Created on: Sept 30, 2020
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_DETAILEDFORCE_H_
#define SRC_GROMACS_FORCEANAL_DETAILEDFORCE_H_

#include "gromacs/math/vectypes.h"
#include "InteractionType.h"

namespace ForceAnal {    
    struct DetailedForce
    {
        DetailedForce(int affected, int applied, InteractionType interact_type, rvec force)
        {
            i = affected;
            j = applied;
            type = interact_type;
            for (int k = 0; k < DIM; ++k)
                f[k] = force[k];
        }

        ~DetailedForce()
        {}

        int i, j;
        InteractionType type;
        rvec f;
    };
}

#endif /* SRC_GROMACS_FORCEANAL_DETAILEDFORCE_H_ */
