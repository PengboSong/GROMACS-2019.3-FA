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

    static bool checkterm2bool(const std::string term, bool def = true);

    static std::string handle_empty_string(const char *str);

protected:
    std::string outpara_fn;

    std::string res_bin_fn;
    std::string res_txt_fn;

    std::string summed_term;
    std::string vector_term;
    std::string scalar_term;

    bool summed_mode;

    OutputType output_type;

    real threshold;

    uint64_t Naverage;

    int group1_id, group2_id;
};

}

#else

struct ForceParaSet
{}

#endif

#endif /* SRC_GROMACS_FORCEANAL_FORCEPARASET_H_ */