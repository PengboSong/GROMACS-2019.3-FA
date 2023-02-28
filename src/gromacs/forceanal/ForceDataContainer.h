/*
    ForceDataConatiner.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/07/28
    Description: Lowlevel data container collections.
*/

// ForceAnal module
#include "ForceAnalDef.h"
#include "AtomForce.h"
#include "InteractionForce.h"
#include "PairForce.h"

namespace ForceAnal
{
    using DetailedAtomForce = OffsetVector<InteractionForce>;
    using DetailedForce = OffsetVector<DetailedAtomForce>;

    using SummedForce = OffsetVector<AtomForce>;

    using PairForceList = std::vector<PairForce>;

}