/*
 * DetailedForce.h
 *
 *  Created on: Sept 30, 2020
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_DETAILEDFORCE_H_
#define SRC_GROMACS_FORCEANAL_DETAILEDFORCE_H_

#include "gromacs/math/extended_vec.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"
#include "ForceAnalConst.h"
#include "InteractionType.h"

namespace ForceAnal {
    struct PairForce
    {
        PairForce(int affected, int applied, InteractionType interact_type, rvec force)
        {
            i = affected;
            j = applied;
            type = interact_type;
            fx = force[XX];
            fy = force[YY];
            fz = force[ZZ];
            f = norm(force);
        }

        ~PairForce()
        {}

        void operator*=(double factor)
        {
            fx *= factor;
            fy *= factor;
            fz *= factor;
            f *= factor;
        }

        int i, j;
        InteractionType type;
        real fx, fy, fz, f;
    };

    struct InteractionForce
    {
        InteractionForce()
        {
            itype = 0;
            clear_rvec(f);
        }

        InteractionForce(InteractionType itype_, rvec f_)
        {
            itype = itype_;
            copy_rvec(f_, f);
        }

        ~InteractionForce()
        {}

        void operator+=(const InteractionForce &F)
        {
            itype = (itype | F.itype);
            rvec_inc(f, F.f);
        }

        void operator*=(double factor)
        {
            dsvimul(f, factor);
        }

        void clear()
        {
            itype = 0;
            clear_rvec(f);
        }

        InteractionType itype;
        rvec f;
    };

    struct DetailedForce
    {
        DetailedForce()
        {
            for (int i = 0; i < Interact_COUNT; ++i)
                clear_rvec(f[i]);
        }

        ~DetailedForce()
        {}

        void clear()
        {
            for (int i = 0; i < Interact_COUNT; ++i)
                clear_rvec(f[i]);
        }

        void append(InteractionType itype_, rvec f_)
        {
            int i = 0;
            while ((itype_ >> i & 1) == 0 && i < Interact_COUNT) ++i;
            rvec_inc(f[i], f_);                    
        }

        void accumulate(InteractionType *itype_, rvec *f_)
        {
            *itype_ = 0;
            clear_rvec(*f_);
            for (int i = 0; i < Interact_COUNT; ++i)
            {
                if (!rvec_small(f[i], NONZERO_LIMIT))
                {
                    itype_ += 1 << i;
                    rvec_inc(*f_, f[i]);
                }
            }
        }
        
        rvec f[Interact_COUNT];
    };
}

#endif /* SRC_GROMACS_FORCEANAL_DETAILEDFORCE_H_ */
