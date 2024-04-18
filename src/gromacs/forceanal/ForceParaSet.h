/*
    ForceParaSet.h
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/04/14
    Description: Initialize parameters for Force Analysis.
*/

#ifndef SRC_GROMACS_FORCEANAL_FORCEPARASET_H_
#define SRC_GROMACS_FORCEANAL_FORCEPARASET_H_

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/variant.h"

#include "ForceAnalDef.h"

namespace ForceAnal {

struct GroupPairParaSet {
    // Output ForceAnal binary format force data filename (-fo *.for)
    std::string res_bin_fn;

    // Output ForceAnal text format force data filename (-ft *.fxt)
    std::string res_txt_fn;

    // Selected mode to handle force data
    DATA_MODE datamode;

    // Selected force output type (scalar, vector or scalar+vector)
    OutputType output_type;

    // Output data in atom unit, residue or molecule
    FORCE_UNIT fu1, fu2;

    // 
    real threshold;
    
    // Group for force analysis
    std::string grp1nm, grp2nm;

    // Packed group index for group1 & group2 (grp1_start, grp1_end, grp2_start, grp2_end)
    // For atom in group 1, its ID should be in range grp1_start <= aid < grp1_end
    GrpIdx grp1idx, grp2idx;

    // Atoms list and mapping index for group1 & group2
    AtomMap grp1, grp2;

    // Force container
    gmx::Variant forces;
};

static constexpr int MAXGRPPAIRN = 32;

class ForceParaSet
{
public:
    ForceParaSet();
    ForceParaSet(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global, const t_inputrec *inputrec);
    ~ForceParaSet();

    void setParas(int nfile, const t_filenm fnm[]);

    void setGroupPairParas(std::vector<t_inpfile>* inp, warninp_t wi, int gi, GroupPairParaSet& paras);

    void mapAtomRes(gmx_mtop_t *top_global);

    static inline bool checkterm2bool(const char* term, bool def = true)
    {
        if (term == INP_YES)
            return true;
        else if (term == INP_NO)
            return false;
        else
            return def;
    }

    static inline std::string handle_empty_string(const char* str)
    {
        return str != nullptr ? std::string(str) : std::string();
    }

    static inline FORCE_UNIT check_forceunit(const char* str)
    {
        std::string kw(str);
        if (kw == "atom")          return FORCE_UNIT::Atom;
        else if (kw == "residue")  return FORCE_UNIT::Residue;
        else if (kw == "molecule") return FORCE_UNIT::Molecule;
        else if (kw == "group")    return FORCE_UNIT::Group;
        else if (kw == "system")   return FORCE_UNIT::System;
        return FORCE_UNIT::Residue;
    }

    static void handle_index(FORCE_UNIT forceunit, int gi, const int* block, const int blocknr, const std::vector<atomindex>& amap, GrpIdx& grpidx, AtomMap& grp);

    static void format_paraset(const uint32_t atomn, const GroupPairParaSet& para);

    static bool check_grpmap(uint32_t atomn, const GrpIdx& grpidx, const AtomMap& grp);

protected:
    // Output ForceAnal parameter filename (*.par)
    std::string outpara_fn;

    // Output ForceAnal binary format force data filename (-fo *.for)
    std::string res_bin_fn;

    // Output ForceAnal text format force data filename (-ft *.fxt)
    std::string res_txt_fn;

    // Output ForceAnal binary format summed force (atom force total, residue force total, nonbonded, nonbonded + bonded) data filename (-fa *.fat)
    std::string totf_bin_fn;
    std::string restotf_bin_fn;
    std::string atomf_nb_bin_fn;
    std::string atomf_nb_b_bin_fn;

    // Output ForceAnal binary format force deviation data filename (-fd *.fat)
    std::string fdev_bin_fn;

    // Input index file (Can not read from gmx_groups_t)
    std::string index_fn;

    // Output ForceAnal map binary data filename (-fmp *.fmp)
    std::string map_fn;

    // Used to check force resolution accuracy
    // Resolved pairwise forces should satisfy equations like
    // Fi = Fij + Fik and Fij = Fji. If numerical solution of expressions
    // including Fi - Fij - Fik and Fij - Fji is close enough to 0, then
    // those necessary equations is considered as satisfied. The upper bound
    // of deviation is limited by this value.
    real force_threshold;

    // Total atom number in system
    uint32_t atomn;

    // Total residue number in system
    uint32_t resn;

    // Total molecule number in system
    uint32_t moln;

    // Total index group number in system
    uint32_t grpn;

    // 
    uint64_t Naverage;

    // Map atom ID to residue ID (both atom and residue IDs begin from 0)
    // e.g. For atom 1 in residue 1, resmap[0] = 0
    std::vector<atomindex> resmap;

    // Map atom ID to molecule ID (both atom and molecule IDs begin from 0)
    // e.g. For atom 1 in molecule 1, molmap[0] = 0
    std::vector<atomindex> molmap;

    // Type of electrostatics
    int eeltype;

    // Type of Van der waals
    int vdwtype;

    // Number of group pairs
    int ngrppairs;

    // Group parameter settings
    std::vector<GroupPairParaSet> grppairparas;
};

}

#endif /* SRC_GROMACS_FORCEANAL_FORCEPARASET_H_ */