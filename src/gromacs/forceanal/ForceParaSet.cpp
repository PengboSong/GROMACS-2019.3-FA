/*
    ForceParaSet.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/04/14
    Description: Initialize parameters for Force Analysis.
*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <utility>

#include "gromacs/commandline/filenm.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/filestream.h"

#include "ForceParaSet.h"

namespace ForceAnal {

ForceParaSet::ForceParaSet()
 : outpara_fn("faout.par"),
   force_threshold(5.0E-3F),
   atomn(0),
   resn(0),
   moln(0),
   Naverage(1),
   eeltype(eelCUT),
   vdwtype(evdwCUT)
{
}

ForceParaSet::ForceParaSet(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global, const t_inputrec *inputrec)
 : outpara_fn("faout.par"),
   force_threshold(5.0E-3F),
   atomn(top_global->natoms),
   Naverage(1),
   eeltype(inputrec->coulombtype),
   vdwtype(inputrec->vdwtype)
{    
    // -fo is optional
    res_bin_fn = handle_empty_string(opt2fn_null("-fo", nfile, fnm));
    // -ft is optional
    res_txt_fn = handle_empty_string(opt2fn_null("-ft", nfile, fnm));
    // -fa is optional
    totf_bin_fn = handle_empty_string(opt2fn_null("-fa", nfile, fnm));
    // Derived other summed force data filenames from the given one
    restotf_bin_fn = modfnm(totf_bin_fn, "", "_res");
    atomf_nb_bin_fn = modfnm(totf_bin_fn, "", "_nb");
    atomf_nb_b_bin_fn = modfnm(totf_bin_fn, "", "_nb+b");
    // -fd is optional
    fdev_bin_fn = handle_empty_string(opt2fn_null("-fd", nfile, fnm));
    // -fn is optional
    index_fn = handle_empty_string(opt2fn_null("-fn", nfile, fnm));
    // -fmp is optional
    map_fn = handle_empty_string(opt2fn_null("-fmp", nfile, fnm));

    mapAtomRes(top_global);
    
    setParas(nfile, fnm);
}

ForceParaSet::~ForceParaSet()
{
}

void ForceParaSet::handle_index(FORCE_UNIT forceunit, int gi, const int* block, const int blocknr, const std::vector<atomindex>& amap, GrpIdx& grpidx, AtomMap& grp)
{
    if (forceunit == FORCE_UNIT::Atom)
    {
        grpidx.first  = *std::min_element(block, block + blocknr);
        grpidx.second = *std::max_element(block, block + blocknr) + 1;
        for (int ni = 0; ni < blocknr; ++ni) grp[block[ni]] = block[ni];
    }
    else if (forceunit == FORCE_UNIT::Residue || forceunit == FORCE_UNIT::Molecule)
    {
        // resmap/molmap: index 0 -> atom 1
        std::unordered_set<int> uniqres;
        for (int ni = 0; ni < blocknr; ++ni) uniqres.insert(amap[block[ni]]);
        grpidx.first  = *std::min_element(uniqres.begin(), uniqres.end());
        grpidx.second = *std::max_element(uniqres.begin(), uniqres.end()) + 1;

        // Instead of adding all atoms in residues, only selecting atoms in given block
        for (int ni = 0; ni < blocknr; ++ni) grp[block[ni]] = amap[block[ni]];
        // for (int ai = 0; ai < amap.size(); ++ai) if (uniqres.count(amap[ai])) grp[ai] = amap[ai];
    }
    else if (forceunit == FORCE_UNIT::Group)
    {
        grpidx.first  = gi;
        grpidx.second = gi + 1;
        for (int ni = 0; ni < blocknr; ++ni) grp[block[ni]] = gi;
    }
}

void ForceParaSet::format_paraset(const uint32_t atomn, const GroupPairParaSet &para)
{
    std::cout << "Write Text Format Force Data = " << (para.res_txt_fn.empty() ? "NO" : "YES") << std::endl;
    std::cout << "Write Binary Format Force Data = " << (para.res_bin_fn.empty() ? "NO" : "YES") << std::endl;
    std::cout << "Force Analysis Mode = " << para.datamode << std::endl;
    std::cout << "Output Data Type = " << para.output_type << std::endl;
    std::cout << "--*-- GROUP A --*--" << std::endl;
    std::cout << "Group Name = " << para.grp1nm << std::endl;
    std::cout << "Force Unit = " << para.fu1 << std::endl;
    std::cout << "Number of Atoms = " << para.grp1.size() << std::endl;
    std::cout << "Nodes Range = " << para.grp1idx << std::endl;
    std::cout << "Atom Map Verified = " << (check_grpmap(atomn, para.grp1idx, para.grp1) ? "Success" : "Failed") << std::endl << std::endl;
    std::cout << "--*-- GROUP B --*--" << std::endl;
    std::cout << "Group Name = " << para.grp2nm << std::endl;
    std::cout << "Force Unit = " << para.fu2 << std::endl;
    std::cout << "Number of Atoms = " << para.grp2.size() << std::endl;
    std::cout << "Nodes Range = " << para.grp2idx << std::endl;
    std::cout << "Atom Map Verified = " << (check_grpmap(atomn, para.grp2idx, para.grp2) ? "Success" : "Failed") << std::endl;
}

bool ForceParaSet::check_grpmap(uint32_t atomn, const GrpIdx &grpidx, const AtomMap &grp)
{
    for (uint32_t ai = 0; ai < atomn; ++ai)
        if (grp.count(ai) && (grp.at(ai) < grpidx.first || grp.at(ai) >= grpidx.second))
            return false;
    return true;
}

void ForceParaSet::mapAtomRes(gmx_mtop_t *top_global)
{
    // Clear resi./mol map
    atomn = top_global->natoms;

    resn = moln = 0;
    resmap.reserve(atomn);
    molmap.reserve(atomn);

    t_atom* patom = nullptr;
    t_atoms* patoms = nullptr;
    for (std::size_t moltypei = 0; moltypei < top_global->moltype.size(); ++moltypei)
    {
        patoms = &top_global->moltype[moltypei].atoms;
        patom = patoms->atom;
        int nmol = top_global->molblock[moltypei].nmol, nres = patoms->nres;
        for (int mi = 0; mi < nmol; ++mi)
        {
            for (int ai = 0; ai < patoms->nr; ++ai)
            {
                resmap.push_back(resn + patom[ai].resind);
                molmap.push_back(moln + mi);
            }
            resn += nres;
        }
        moln += nmol;
    }
    
    if (!map_fn.empty())
    {
        std::ofstream mapstream(map_fn, std::ios::binary | std::ios::trunc);
        mapstream.write((char*)&atomn, sizeof(uint32_t));
        mapstream.write((char*)resmap.data(), sizeof(atomindex) * resmap.size());
        mapstream.write((char*)molmap.data(), sizeof(atomindex) * molmap.size());
    }
}

void ForceParaSet::setParas(int nfile, const t_filenm fnm[])
{
    warninp_t wi = init_warning(FALSE, 0);
    std::string FA_paraset_fn;
    if (opt2bSet("-fp", nfile, fnm))
        FA_paraset_fn = handle_empty_string(opt2fn("-fp", nfile, fnm));
    
    std::vector<t_inpfile> inp;
    if (!FA_paraset_fn.empty())
    {
        gmx::TextInputFile inpara_stream(FA_paraset_fn);
        inp = read_inpfile(&inpara_stream, FA_paraset_fn.c_str(), wi);
    }
    else
        inp.clear();

    if (!checkterm2bool(get_estr(&inp, "summed-force-atom", "yes")))
        totf_bin_fn.clear();
    if (!checkterm2bool(get_estr(&inp, "summed-force-residue", "yes")))
        restotf_bin_fn.clear();
    if (!checkterm2bool(get_estr(&inp, "summed-force-nonbonded", "no")))
        atomf_nb_bin_fn.clear();
    if (!checkterm2bool(get_estr(&inp, "summed-force-nonbonded-bonded", "no")))
        atomf_nb_b_bin_fn.clear();

    Naverage = get_eint64(&inp, "naverage", 1, wi);

    ngrppairs = get_eint(&inp, "group-pairs", 0, wi);
    if (ngrppairs > MAXGRPPAIRN)
        gmx_warning("Force Analysis module only supports up to %d group pairs.", MAXGRPPAIRN);
    grppairparas.resize(ngrppairs);

    for (int gi = 0; gi < ngrppairs; ++gi) setGroupPairParas(&inp, wi, gi, grppairparas[gi]);

    // Format output parameters
    std::cout << std::endl << "Write Atom Summed Force = " << (totf_bin_fn.empty() ? "NO" : "YES") << std::endl;
    std::cout << "Write Residue Summed Force = " << (restotf_bin_fn.empty() ? "NO" : "YES") << std::endl;
    std::cout << "Write Nonbonded Summed Force = " << (atomf_nb_bin_fn.empty() ? "NO" : "YES") << std::endl;
    std::cout << "Write Nonbonded + Bonded Summed Force = " << (atomf_nb_b_bin_fn.empty() ? "NO" : "YES") << std::endl;
    std::cout << "Write Force Decomposition Residual Error = " << (fdev_bin_fn.empty() ? "NO" : "YES") << std::endl;
    std::cout << std::endl << "Average Frequency = " << Naverage << std::endl;
    std::cout << "Number of Group Pairs = " << ngrppairs << std::endl;
    for (int gi = 0; gi < ngrppairs; ++gi)
    {
        std::cout << "==*==*== GROUP " << gi << " ==*==*==" << std::endl;
        format_paraset(atomn, grppairparas[gi]);
        std::cout << "==*==*==*==*==*==*==*==*==" << std::endl << std::endl;
    }

    gmx::TextOutputFile outpara_stream(outpara_fn);
    write_inpfile(&outpara_stream, outpara_fn.c_str(), &inp, FALSE, WriteMdpHeader::yes, wi);
}

void ForceParaSet::setGroupPairParas(std::vector<t_inpfile>* inp, warninp_t wi, int gi, GroupPairParaSet &paras)
{
    std::string pfx = "grp" + std::to_string(gi) + "-", sfx = "_grp" + std::to_string(gi);

    std::string datamode_term = get_estr(inp, pfx + "analmode", "summed");
    if (datamode_term == "summed")
        paras.datamode = DATA_MODE::SummedMode;
    else if (datamode_term == "detailed")
        paras.datamode = DATA_MODE::DetailedMode;
    else if (datamode_term == "list")
        paras.datamode = DATA_MODE::ListMode;
    else
        paras.datamode = DATA_MODE::None;

    if (res_bin_fn.empty() && res_txt_fn.empty())
        paras.datamode = DATA_MODE::None;

    if (!res_bin_fn.empty())
        paras.res_bin_fn = modfnm(res_bin_fn, "", sfx);
    if (!res_txt_fn.empty())
        paras.res_txt_fn = modfnm(res_txt_fn, "", sfx);
    
    paras.output_type = OUT_NOTHING;
    if (checkterm2bool(get_estr(inp, pfx + "vector", "yes")))
        paras.output_type |= OUT_VECTOR;
    if (checkterm2bool(get_estr(inp, pfx + "scalar", "yes")))
        paras.output_type |= OUT_SCALAR;
    
    real t = get_ereal(inp, pfx + "threshold", 1.0E-3F, wi);
    paras.threshold = t * t;

    paras.fu1 = check_forceunit(get_estr(inp, pfx + "forceunitA", "residue"));
    paras.fu2 = check_forceunit(get_estr(inp, pfx + "forceunitB", "residue"));
    paras.grp1nm = get_estr(inp, pfx + "groupA", "System");
    paras.grp2nm = get_estr(inp, pfx + "groupB", "System");

    if (index_fn.empty())
        throw std::runtime_error("Index file is required for group pairs analysis.");
    
    char** grpnms;
    t_blocka* grps = init_index(index_fn.c_str(), &grpnms);
    for (int gi = 0; gi < grps->nr; ++gi)
    {
        int blocks = grps->index[gi];
        if (grpnms[gi] == paras.grp1nm)
            handle_index(paras.fu1, gi, &grps->a[blocks], grps->index[gi + 1] - blocks, paras.fu1 == ForceAnal::FORCE_UNIT::Molecule ? molmap : resmap, paras.grp1idx, paras.grp1);
        if (grpnms[gi] == paras.grp2nm)
            handle_index(paras.fu2, gi, &grps->a[blocks], grps->index[gi + 1] - blocks, paras.fu2 == ForceAnal::FORCE_UNIT::Molecule ? molmap : resmap, paras.grp2idx, paras.grp2);
    }
    if (paras.fu1 == FORCE_UNIT::System)
    {
        paras.grp1idx = std::make_pair<atomindex, atomindex>(0, 1);
        for (atomindex ai = 0; ai < atomn; ++ai) paras.grp1[ai] = 0;
    }
    if (paras.fu2 == FORCE_UNIT::System)
    {
        paras.grp2idx = std::make_pair<atomindex, atomindex>(0, 1);
        for (atomindex ai = 0; ai < atomn; ++ai) paras.grp2[ai] = 0;
    }        
}
}
