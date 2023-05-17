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

static inline void rvec_opp(rvec a)
{
    real x, y, z;

    x = -a[XX];
    y = -a[YY];
    z = -a[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

static inline bool rvec_small(const rvec a, real b)
{
    if (a[XX] > b)
        return false;
    if (a[YY] > b)
        return false;
    if (a[ZZ] > b)
        return false;
    return true;
}

static inline bool rvec_abs_small(const rvec a, real b)
{
    if (std::abs(a[XX]) > b)
        return false;
    if (std::abs(a[YY]) > b)
        return false;
    if (std::abs(a[ZZ]) > b)
        return false;
    return true;
}

static inline void svimul(rvec a, real b)
{
    real x, y, z;

    x = b * a[XX];
    y = b * a[YY];
    z = b * a[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

static inline void dsvimul(rvec a, double b)
{
    double x, y, z;

    x = b * a[XX];
    y = b * a[YY];
    z = b * a[ZZ];

    a[XX] = x;
    a[YY] = y;
    a[ZZ] = z;
}

static inline void copy_rvec_to_2dvec(const rvec a, double *b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
}

static inline void copy_rvec_to_3dvec(const rvec a, double *b)
{
    b[XX] = a[XX];
    b[YY] = a[YY];
    b[ZZ] = a[ZZ];
}

#endif /* SRC_GROMACS_MATH_EXTENDED_VEC_H_ */