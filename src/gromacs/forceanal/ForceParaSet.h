/*
 * ForceParaSet.h
 *
 *  Created on: Apr 14, 2021
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_FORCESETTINGS_H_
#define SRC_GROMACS_FORCEANAL_FORCESETTINGS_H_

#include "gromacs/commandline/filenm.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus

#include <string>
#include <cstdint>

namespace ForceAnal {

class ForceParaSet
{
public:
    ForceParaSet();
    ForceParaSet(int nfile, const t_filenm fnm[], gmx_mtop_t *mtop);
    ~ForceParaSet();

    void reset_filename();

    void check_average_steps();

protected:
    std::string result_filename;
    std::string result_full_filename;
    
    bool write_in_binary;

    uint32_t Nevery;
    uint32_t Nrepeat;
    uint64_t Nfreq;

    int group1_id, group2_id;
};

}

#else

struct ForceParaSet
{}

#endif

#endif /* SRC_GROMACS_FORCEANAL_FORCESETTINGS_H_ */