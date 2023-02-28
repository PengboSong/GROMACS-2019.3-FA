/*
    ForceParaSet.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/04/14
    Description: Initialize parameters for Force Analysis.
*/

#ifndef SRC_GROMACS_FORCEANAL_FORCEPARASET_H_
#define SRC_GROMACS_FORCEANAL_FORCEPARASET_H_

#include "gromacs/commandline/filenm.h"
#include "gromacs/topology/topology.h"

#include "ForceAnalDef.h"

namespace ForceAnal {

class ForceParaSet
{
public:
    ForceParaSet();
    ForceParaSet(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global);
    ~ForceParaSet();

    void setParas(int nfile, const t_filenm fnm[]);

    void setGroup();

    void mapAtomRes(gmx_mtop_t *top_global);

    static bool checkterm2bool(const char* term, bool def = true);

    static std::string handle_empty_string(const char *str);

    static void handle_index(FORCE_UNIT forceunit, const int* block, const int blocknr, const std::vector<atomindex>& resmap, GrpIdx& grpidx, std::vector<atomindex>& excl);

protected:
    // Output ForceAnal parameter filename (*.par)
    std::string outpara_fn;

    // Output ForceAnal binary format force data filename (-fo *.for)
    std::string res_bin_fn;

    // Output ForceAnal text format force data filename (-ft *.fxt)
    std::string res_txt_fn;

    // Output ForceAnal binary format atom force data filename (-fa *.fat)
    std::string totf_bin_fn;

    // Input index file (Can not read from gmx_groups_t)
    std::string index_fn;

    // Output ForceAnal map binary data filename (-fmp *.fmp)
    std::string map_fn;

    // Selected mode to handle force data
    DATA_MODE datamode;

    // Selected force output type (scalar, vector or scalar+vector)
    OutputType output_type;

    // Output data in atom unit, residue or molecule
    FORCE_UNIT forceunit;

    // 
    real threshold;

    // Total atom number in system
    uint32_t atomn;

    // Total residue number in system
    uint32_t resn;

    // Total molecule number in system
    uint32_t moln;

    // 
    uint64_t Naverage;

    // Map atom ID to residue ID (both atom and residue IDs begin from 0)
    // e.g. For atom 1 in residue 1, resmap[0] = 0
    std::vector<atomindex> resmap;

    // Map atom ID to molecule ID (both atom and residue IDs begin from 0)
    // e.g. For atom 1 in molecule 1, molmap[0] = 0
    std::vector<atomindex> molmap;
    
    // Group for force analysis
    std::string grp1nm, grp2nm;

    // Packed group index for group1 & group2 (grp1_start, grp1_end, grp2_start, grp2_end)
    // For atom in group 1, its ID should be in range grp1_start <= aid < grp1_end
    GrpIdx grp1idx, grp2idx;

    // Index excluded from group1 & group2
    std::vector<atomindex> exclgrp1, exclgrp2;
};

}

#endif /* SRC_GROMACS_FORCEANAL_FORCEPARASET_H_ */