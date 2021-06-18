/*
 * InteractionForce.h
 *
 *  Created on: June 15, 2021
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_INTERACTIONFORCE_H_
#define SRC_GROMACS_FORCEANAL_INTERACTIONFORCE_H_

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#include "InteractionType.h"

namespace ForceAnal {
    static const uint8_t XYZ = 3;
    static const uint8_t Interact_FORCEVEC_LEN = 4 * (Interact_COUNT + 1);

    struct InteractionForce
    {
        InteractionForce()
        {
            init();
        }

        ~InteractionForce()
        {}

        void init()
        {
            forces.assign(Interact_FORCEVEC_LEN, 0.);
        }

        void add(InteractionType itype, rvec f)
        {
            uint8_t idx = index_itype(itype);
            if (idx != 0)
            {
                idx = 4 * (idx - 1);
                forces[idx + XX] = f[XX];
                forces[idx + YY] = f[YY];
                forces[idx + ZZ] = f[ZZ];
            }
        }

        real sum()
        {
            uint8_t offset = 4 * Interact_COUNT;
            real fx, fy, fz;
            for (uint8_t i = 0; i < offset;)
            {
                fx = forces[i++];
                fy = forces[i++];
                fz = forces[i++];
                forces[i++] = fx * fx + fy * fy + fz * fz;
                forces[offset + XX] += fx;
                forces[offset + YY] += fy;
                forces[offset + ZZ] += fz;
            }
            fx = forces[offset + XX];
            fy = forces[offset + YY];
            fz = forces[offset + ZZ];
            forces[offset + XYZ] = fx * fx + fy * fy + fz * fz;
            return forces[offset + XYZ];
        }

        void operator*=(double factor)
        {
            for (uint8_t i = 0; i < Interact_FORCEVEC_LEN; ++i)
                forces[i] *= factor;
        }

        std::vector<real> forces;
    };
}

#endif /* SRC_GROMACS_FORCEANAL_INTERACTIONFORCE_H_ */
