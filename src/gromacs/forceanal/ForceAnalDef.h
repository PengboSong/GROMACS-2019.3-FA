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
#include <ostream>
#include <string>
#include <utility>
#include <unordered_set>
#include <unordered_map>

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

    // Paired index IDs
    using GrpIdx = std::pair<atomindex, atomindex>;

    // Group atom map
    using AtomMap = std::unordered_map<atomindex, atomindex>;

    static inline std::ostream &operator<<(std::ostream &os, const GrpIdx& idx)
    {
        os << '[' << idx.first << ',' << idx.second << ')';
        return os;
    }

    using OutputType = uint8_t;

    static const OutputType OUT_NOTHING =      0;
    static const OutputType OUT_VECTOR  = 1 << 0;
    static const OutputType OUT_SCALAR  = 1 << 1;

    static inline std::ostream &operator<<(std::ostream &os, OutputType otp)
    {
        switch (otp)
        {
            case OUT_NOTHING:
                os << "None";
                break;
            case OUT_VECTOR:
                os << "Vector";
                break;
            case OUT_SCALAR:
                os << "Scalar";
                break;
            case OUT_VECTOR + OUT_SCALAR:
                os << "Vector + Scalar";
                break;
        }
        return os;
    }

    static const std::string INP_YES = "yes";
    static const std::string INP_NO = "no";

    enum class DATA_MODE : uint8_t
    {
        None          = 0U,
        SummedMode    = 1U,
        DetailedMode  = 2U,
        ListMode      = 3U,
        AtomForceMode = 4U
    };
    
    static inline std::ostream &operator<<(std::ostream &os, DATA_MODE mode)
    {
        switch (mode)
        {
            case DATA_MODE::None:
                os << "None";
                break;
            case DATA_MODE::SummedMode:
                os << "SummedMode";
                break;
            case DATA_MODE::DetailedMode:
                os << "DetailedMode";
                break;
            case DATA_MODE::ListMode:
                os << "ListMode";
                break;
            case DATA_MODE::AtomForceMode:
                os << "AtomForceMode";
                break;
        }
        return os;
    }

    enum class FORCE_UNIT : uint8_t
    {
        None     = 0U,
        Atom     = 1U,
        Residue  = 2U,
        Molecule = 3U,
        Group    = 4U,
        System   = 5U
    };
    
    static inline std::ostream &operator<<(std::ostream &os, FORCE_UNIT fu)
    {
        switch (fu)
        {
            case FORCE_UNIT::None:
                os << "None";
                break;
            case FORCE_UNIT::Atom:
                os << "Atom";
                break;
            case FORCE_UNIT::Residue:
                os << "Residue";
                break;
            case FORCE_UNIT::Molecule:
                os << "Molecule";
                break;
            case FORCE_UNIT::Group:
                os << "Group";
                break;
            case FORCE_UNIT::System:
                os << "System";
                break;
        }
        return os;
    }

    enum class OUT_FORCE_TYPE : uint8_t
    {
        None                     = 0U,
        AtomForceTotal           = 1U,
        AtomForceNonbonded       = 2U,
        AtomForceNonbondedBonded = 3U,
    };
    
    static inline std::ostream &operator<<(std::ostream &os, OUT_FORCE_TYPE oftp)
    {
        switch (oftp)
        {
            case OUT_FORCE_TYPE::None:
                os << "None";
                break;
            case OUT_FORCE_TYPE::AtomForceTotal:
                os << "Total Atom Force";
                break;
            case OUT_FORCE_TYPE::AtomForceNonbonded:
                os << "Nonbonded Atom Force";
                break;
            case OUT_FORCE_TYPE::AtomForceNonbondedBonded:
                os << "Nonbonded + Bonded Atom Force";
                break;
        }
        return os;
    }

    static inline real real_norm2(const real fx, const real fy, const real fz)
    {
        return fx * fx + fy * fy + fz * fz;
    }

    static inline real real_norm(const real fx, const real fy, const real fz)
    {
        return std::sqrt(real_norm2(fx, fy, fz));
    }

    static inline std::string modfnm(const std::string fnm, const std::string prefix, const std::string suffix)
    {
        return prefix + (suffix.empty() ? fnm : (fnm.find_first_of('.') == std::string::npos ? fnm + suffix : fnm.substr(0, fnm.find_first_of('.')) + suffix + fnm.substr(fnm.find_first_of('.'))));
    }
} // ForceAnal

#endif /* SRC_GROMACS_FORCEANAL_FORCEANALDEF_H_ */
