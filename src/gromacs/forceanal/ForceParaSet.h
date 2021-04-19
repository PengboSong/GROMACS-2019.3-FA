/*
 * ForceParaSet.h
 *
 *  Created on: Apr 14, 2021
 *      Author: Pengbo Song
 */

#ifndef SRC_GROMACS_FORCEANAL_FORCEPARASET_H_
#define SRC_GROMACS_FORCEANAL_FORCEPARASET_H_

#include "ForceIOtype.h"

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

    void check_average_steps();

    void set_parameters(int nfile, const t_filenm fnm[]);

    bool checkterm2bool(const char* term, bool def = true);

protected:
    std::string result_binary_filename;
    std::string result_text_filename;

    std::string summed_term;
    std::string vector_term;
    std::string scalar_term;

    bool summed_mode;

    OutputType output_type;

    real threshold;

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

#endif /* SRC_GROMACS_FORCEANAL_FORCEPARASET_H_ */