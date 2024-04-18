/*
    InteractionForce.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/06/15
    Description: Container for force vectors with different interaction types.
*/

#ifndef SRC_GROMACS_FORCEANAL_INTERACTIONFORCE_H_
#define SRC_GROMACS_FORCEANAL_INTERACTIONFORCE_H_

#include "ForceAnalDef.h"

namespace ForceAnal {

static const uint8_t Interact_FORCEVEC_LEN = 4 * (Interact_COUNT + 1);

class InteractionForce
{
using index = uint8_t;

public:
    InteractionForce()
    {
        forces.resize(Interact_FORCEVEC_LEN, 0.);
    }

    real* data() { return forces.data(); }

    void clear()
    {
        forces.assign(Interact_FORCEVEC_LEN, 0.);
    }

    void add(InteractionType itype, rvec f)
    {
        index idx = itype2index(itype);
        if (idx != 0)
        {
            idx = 4 * (idx - 1);
            forces[idx + XX] += f[XX];
            forces[idx + YY] += f[YY];
            forces[idx + ZZ] += f[ZZ];
        }
    }

    real sum()
    {
        index offset = 4 * Interact_COUNT;
        real sumfx = 0., sumfy = 0., sumfz = 0., sumf = 0.;
        real fx, fy, fz;
        for (index i = 0; i < offset;)
        {
            fx = forces[i++];
            fy = forces[i++];
            fz = forces[i++];
            forces[i++] = real_norm2(fx, fy, fz);
            sumfx += fx;
            sumfy += fy;
            sumfz += fz;
        }
        forces[offset + XX] = sumfx;
        forces[offset + YY] = sumfy;
        forces[offset + ZZ] = sumfz;
        sumf = real_norm2(sumfx, sumfy, sumfz);
        forces[offset + XYZ] = sumf;
        return sumf;
    }

    void operator*=(double factor)
    {
        for (index i = 0; i < Interact_FORCEVEC_LEN; ++i)
            forces[i] *= factor;
    }

private:
    std::vector<real> forces;
};

}

#endif /* SRC_GROMACS_FORCEANAL_INTERACTIONFORCE_H_ */
