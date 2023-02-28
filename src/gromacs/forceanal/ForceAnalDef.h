/*
    ForceAnalDef.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/04/19
    Description: Collection of definitions in Force Analysis module.
*/

#ifndef SRC_GROMACS_FORCEANAL_FORCEANALDEF_H_
#define SRC_GROMACS_FORCEANAL_FORCEANALDEF_H_

// C++ STL
#include <cstdint>
#include <vector>
#include <map>
#include <string>
#include <utility>

// GROMACS include
#include "gromacs/math/extended_vec.h"   // Additional rvec operations
#include "gromacs/math/vec.h"            // GROMACS rvec operations
#include "gromacs/math/vectypes.h"       // For rvec
#include "gromacs/topology/ifunc.h"      // For rvec4
#include "gromacs/utility/real.h"        // For real

// ForceAnal module
#include "BiMap.h"
#include "InteractionType.h"
#include "OffsetVector.h"

namespace ForceAnal
{
    /*
    ForceAnal module uses 4-size force vectors. The first 3 real variables are
    X, Y, Z components of a force vector, the same as rvec or dvec in GROMACS,
    and the last one is used for the scalar value of the force. GROMACS provides
    3 aliases to fetch X, Y, Z components at the front, they are XX, YY and ZZ,
    respectively. Therefore, another alias to fetch the scalar value to the end
    is required.
    */
    static const uint8_t XYZ = 3;

    // Force with a smaller scalar value than that will be ignored
    static const real NONZERO_LIMIT = 1e-6;

    using atomindex = int32_t;

    // atom2res map
    using Atom2Res = std::map<atomindex, atomindex>;

    // res2atom map
    using Res2Atom = std::map<atomindex, std::vector<atomindex>>;

    // Group index IDs
    using GrpIdx = std::pair<atomindex, atomindex>;

    using OutputType = uint8_t;

    static const OutputType OUT_NOTHING =      0;
    static const OutputType OUT_VECTOR  = 1 << 0;
    static const OutputType OUT_SCALAR  = 1 << 1;

    static const std::string INP_YES = "yes";
    static const std::string INP_NO = "no";

    enum class DATA_MODE : uint8_t
    {
        None = 0U,
        SummedMode = 1U,
        DetailedMode = 2U,
        ListMode = 3U,
        AtomForceMode = 4U,
    };

    enum class FORCE_UNIT : uint8_t
    {
        Atom = 1U,
        Residue = 2U,
        Molecule = 3U,
    };

    static inline real real_norm2(const real fx, const real fy, const real fz)
    {
        return fx * fx + fy * fy + fz * fz;
    }

    static inline real real_norm(const real fx, const real fy, const real fz)
    {
        return std::sqrt(real_norm2(fx, fy, fz));
    }
} // ForceAnal

#endif /* SRC_GROMACS_FORCEANAL_FORCEANALDEF_H_ */
