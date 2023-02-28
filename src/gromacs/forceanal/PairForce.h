/*
    PairForce.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/09/30
    Description: Data container for a single force pair.
*/

#ifndef SRC_GROMACS_FORCEANAL_PAIRFORCE_H_
#define SRC_GROMACS_FORCEANAL_PAIRFORCE_H_

#include "ForceAnalDef.h"

namespace ForceAnal {
    struct PairForce
    {
        PairForce(int affected, int applied, InteractionType interact_type, rvec force)
        {
            i = affected;
            j = applied;
            type = interact_type;
            f[XX] = force[XX];
            f[YY] = force[YY];
            f[ZZ] = force[ZZ];
            f[XYZ] = norm2(force);
        }

        void operator*=(double factor)
        {
            for (int i = 0; i < 4; ++i)
                f[i] *= factor;
        }

        int i, j;
        InteractionType type;
        rvec4 f;
    };
}

#endif /* SRC_GROMACS_FORCEANAL_PAIRFORCE_H_ */
