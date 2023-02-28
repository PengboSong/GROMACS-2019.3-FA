/*
    AtomForce.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/07/16
    Description: 
*/

#ifndef SRC_GROMACS_FORCEANAL_ATOMFORCE_H_
#define SRC_GROMACS_FORCEANAL_ATOMFORCE_H_

#include "ForceAnalDef.h"
#include "InteractionForce.h"

namespace ForceAnal {

class AtomForce
{
using index = uint64_t;

public:
    AtomForce() : offset(0), atomn(0), forcelen(0) {}

    void clear()
    {
        itypes.assign(atomn, 0);
        forces.assign(forcelen, 0.);
    }

    index size() { return atomn; }

    index id_begin() { return offset; }

    index id_end() { return offset + atomn; }

    InteractionType* itypes_data() { return itypes.data(); }

    real* forces_data() { return forces.data(); }

    void resize(index sloc, index eloc)
    {
        offset = sloc;
        atomn = eloc - sloc;
        forcelen = 4 * atomn;
        itypes.resize(atomn, 0);
        forces.resize(forcelen, 0.);
    }

    void add(index ai, InteractionType itype, rvec force)
    {
        if (ai < offset) return;            
        index loc = ai - offset;
        if (loc >= atomn) return;
        index floc = 4 * loc;
        forces[floc + XX] += force[XX];
        forces[floc + YY] += force[YY];
        forces[floc + ZZ] += force[ZZ];
        itypes[loc] |= itype;
    }

    void norm2()
    {
        real fx, fy, fz;
        for (index i = 0; i < forcelen;)
        {
            fx = forces[i++];
            fy = forces[i++];
            fz = forces[i++];
            forces[i++] = fx * fx + fy * fy + fz * fz;
        }
    }

    void operator*=(double factor)
    {
        for (index i = 0; i < forcelen; ++i)
            forces[i] *= factor;
    }

private:
    index offset;

    index atomn;

    index forcelen;

    std::vector<InteractionType> itypes;

    std::vector<real> forces;
};

}

#endif /* SRC_GROMACS_FORCEANAL_ATOMFORCE_H_ */
