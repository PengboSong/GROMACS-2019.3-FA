/*
 * ForceIOtype.h
 *
 *  Created on: Apr 19, 2021
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_FORCEIOTYPE_H_
#define SRC_GROMACS_FORCEANAL_FORCEIOTYPE_H_

#include <cstdint>
#include <string>

namespace ForceAnal {

using OutputType = uint8_t;

static const OutputType OUT_NOTHING =      0;
static const OutputType OUT_VECTOR  = 1 << 0;
static const OutputType OUT_SCALAR  = 1 << 1;

static const std::string INP_YES = "yes";
static const std::string INP_NO = "no";

}

#endif /* SRC_GROMACS_FORCEANAL_FORCEIOTYPE_H_ */