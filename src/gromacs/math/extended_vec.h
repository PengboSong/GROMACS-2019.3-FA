#ifndef SRC_GROMACS_MATH_EXTENDED_VEC_H_
#define SRC_GROMACS_MATH_EXTENDED_VEC_H_

#include <cmath>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

static inline void rvec_opp(const rvec a, rvec b)
{
    real x, y, z;

    x = -a[XX];
    y = -a[YY];
    z = -a[ZZ];

    b[XX] = x;
    b[YY] = y;
    b[ZZ] = z;
}

#endif /* SRC_GROMACS_MATH_EXTENDED_VEC_H_ */