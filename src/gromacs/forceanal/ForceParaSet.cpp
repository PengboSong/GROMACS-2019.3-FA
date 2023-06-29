/*
    ForceParaSet.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2021/04/14
    Description: Initialize parameters for Force Analysis.
*/

#include <algorithm>
#include <set>
#include <fstream>
#include <utility>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/topology/index.h"
#include "gromacs/utility/filestream.h"

#include "ForceParaSet.h"

namespace ForceAnal {

ForceParaSet::ForceParaSet()
 : outpara_fn("faout.par"),
   datamode(DATA_MODE::None),
   output_type(OUT_NOTHING),
   forceunit(FORCE_UNIT::Atom),
   threshold(1.0E-3F),
   force_threshold(5.0E-3F),
   atomn(0),
   resn(0),
   moln(0),
   Naverage(1),
   grp1idx(0, 0),
   grp2idx(0, 0),
   eeltype(eelCUT),
   vdwtype(evdwCUT)
{
}

ForceParaSet::ForceParaSet(int nfile, const t_filenm fnm[], gmx_mtop_t *top_global, const t_inputrec *inputrec)
 : outpara_fn("faout.par"),
   datamode(DATA_MODE::None),
   output_type(OUT_NOTHING),
   forceunit(FORCE_UNIT::Atom),
   threshold(1.0E-3F),
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
    // -fd is optional
    fdev_bin_fn = handle_empty_string(opt2fn_null("-fd", nfile, fnm));
    // -fn is optional
    index_fn = handle_empty_string(opt2fn_null("-fn", nfile, fnm));
    // -fmp is optional
    map_fn = handle_empty_string(opt2fn_null("-fmp", nfile, fnm));

    // Derived other atom force data filenames from the specified one
    std::string::size_type sz = totf_bin_fn.find_first_of('.');
    if (sz == std::string::npos)
    {
        atomf_nb_bin_fn = totf_bin_fn + "_nb";
        atomf_nb_b_bin_fn = totf_bin_fn + "_nb+b";
    }
    else
    {
        std::string atomf_prefix = totf_bin_fn.substr(0, sz), atomf_suffix = totf_bin_fn.substr(sz);
        atomf_nb_bin_fn = atomf_prefix + "_nb" + atomf_suffix;
        atomf_nb_b_bin_fn = atomf_prefix + "_nb+b" + atomf_suffix;
    }
    
    setParas(nfile, fnm);

    mapAtomRes(top_global);

    setGroup();
}

ForceParaSet::~ForceParaSet()
{
}

bool ForceParaSet::checkterm2bool(const char* term, bool def)
{
    if (term == INP_YES)
        return true;
    else if (term == INP_NO)
        return false;
    else
        return def;
}

std::string ForceParaSet::handle_empty_string(const char *str)
{
    if (str != nullptr)
        return std::string(str);
    else
        return std::string();
}

void ForceParaSet::handle_index(FORCE_UNIT forceunit, const int* block, const int blocknr, const std::vector<atomindex>& resmap, GrpIdx& grpidx, GrpIdx& grpaid, std::vector<atomindex>& excl)
{
    int grpstart, grpend, grplen;
    if (forceunit == FORCE_UNIT::Atom)
    {
        grpstart = *std::min_element(block, block + blocknr);
        grpend = *std::max_element(block, block + blocknr) + 1;
        grplen = grpend - grpstart;
        grpidx.first = grpstart;
        grpidx.second = grpend;
        grpaid.first = grpstart;
        grpaid.second = grpend;
        if (blocknr < grplen)
        {
            std::vector<int> filled_range(grplen, 0);
            for (int ni = 0; ni < blocknr; ++ni)
                filled_range[block[ni] - grpstart] = 1;
            for (std::size_t idx = 0; idx < filled_range.size(); ++idx)
                if (filled_range[idx] == 0)
                    excl.push_back(grpstart + idx);
        }
    }
    else
    {
        grpidx.first = *std::min_element(block, block + blocknr);
        grpidx.second = *std::max_element(block, block + blocknr) + 1;

        const char* identifier = (forceunit == FORCE_UNIT::Residue) ? "resi." : "mol";
        
        std::set<int> uniqres;
        // resmap/molmap: index 0 -> atom 1
        for (int ni = 0; ni < blocknr; ++ni)
            uniqres.insert(resmap[block[ni]]);
        grpstart = *std::min_element(uniqres.begin(), uniqres.end());
        grpend = *std::max_element(uniqres.begin(), uniqres.end()) + 1;
        grplen = grpend - grpstart;
        grpidx.first = grpstart;
        grpidx.second = grpend;
        if (uniqres.size() < static_cast<std::size_t>(grplen))
        {
            std::vector<int> filled_range(grplen, 0);
            for (const int& resi : uniqres)
                filled_range[resi - grpstart] = 1;
            for (std::size_t idx = 0; idx < filled_range.size(); ++idx)
                if (filled_range[idx] == 0)
                    excl.push_back(grpstart + idx);
        }
    }
    const char* identifier;
    switch (forceunit)
    {
        case FORCE_UNIT::Atom:
            identifier = "atom";
            break;
        case FORCE_UNIT::Residue:
            identifier = "resi.";
            break;
        case FORCE_UNIT::Molecule:
            identifier = "mol";
            break;
    }
    printf("Force group starts from %s %d to %s %d\n", identifier, grpstart, identifier, grpend - 1);
}

void ForceParaSet::mapAtomRes(gmx_mtop_t *top_global)
{
    // Clear resi./mol map
    atomn = top_global->natoms;

    resn = 0;
    resmap.clear();
    resmap.reserve(atomn);

    moln = 0;
    molmap.clear();
    molmap.reserve(atomn);

    int nmol = 0, nres = 0;
    t_atom* patom = nullptr;
    t_atoms* patoms = nullptr;
    for (std::size_t moltypei = 0; moltypei < top_global->moltype.size(); ++moltypei)
    {
        patoms = &top_global->moltype[moltypei].atoms;
        patom = patoms->atom;
        nmol = top_global->molblock[moltypei].nmol;
        nres = patoms->nres;
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
    
    if (!map_fn.empty() && ((forceunit == FORCE_UNIT::Residue) || (forceunit == FORCE_UNIT::Molecule)))
    {
        std::ofstream mapstream(map_fn, std::ios::binary | std::ios::trunc);
        uint8_t filecode = static_cast<uint8_t>(forceunit);
        mapstream.write((char*)&filecode, sizeof(uint8_t));
        mapstream.write((char*)&atomn, sizeof(uint32_t));
        switch (forceunit)
        {
            case FORCE_UNIT::Residue:
                mapstream.write((char*)resmap.data(), sizeof(atomindex) * resmap.size());
                break;
            case FORCE_UNIT::Molecule:
                mapstream.write((char*)molmap.data(), sizeof(atomindex) * molmap.size());
                break;
        }
    }
}

void ForceParaSet::setGroup()
{
    if (index_fn.empty())
    {
        grp1idx = std::make_pair<atomindex, atomindex>(0, atomn);
        grp2idx = std::make_pair<atomindex, atomindex>(0, atomn);
        grp1aid = std::make_pair<atomindex, atomindex>(0, atomn);
        grp2aid = std::make_pair<atomindex, atomindex>(0, atomn);
    }
    else
    {
        char** grpnms;
        t_blocka* grps = init_index(index_fn.c_str(), &grpnms);
        for (int gi = 0; gi < grps->nr; ++gi)
        {
            int blocks = grps->index[gi];
            if (grpnms[gi] == grp1nm)
                handle_index(forceunit, &grps->a[blocks], grps->index[gi + 1] - blocks, resmap, grp1idx, grp1aid, exclgrp1);
            if (grpnms[gi] == grp2nm)
                handle_index(forceunit, &grps->a[blocks], grps->index[gi + 1] - blocks, resmap, grp2idx, grp2aid, exclgrp2);
        }
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

    std::string datamode_term = get_estr(&inp, "data-mode", "summed");

    if (datamode_term == "summed")
        datamode = DATA_MODE::SummedMode;
    else if (datamode_term == "detailed")
        datamode = DATA_MODE::DetailedMode;
    else if (datamode_term == "list")
        datamode = DATA_MODE::ListMode;
    else
        datamode = DATA_MODE::None;

    if (res_bin_fn.empty() && res_txt_fn.empty())
        datamode = DATA_MODE::None;
    
    if (checkterm2bool(get_estr(&inp, "vector", "yes")))
        output_type |= OUT_VECTOR;
    if (checkterm2bool(get_estr(&inp, "scalar", "yes")))
        output_type |= OUT_SCALAR;

    std::string forceunit_term = get_estr(&inp, "force-unit", "residue");
    if (forceunit_term == "atom")
        forceunit = FORCE_UNIT::Atom;
    else if (forceunit_term == "residue")
        forceunit = FORCE_UNIT::Residue;
    else if (forceunit_term == "molecule")
        forceunit = FORCE_UNIT::Molecule;
    
    threshold = get_ereal(&inp, "threshold", 1.0E-3F, wi);
    threshold *= threshold;
    Naverage = get_eint64(&inp, "naverage", 1, wi);
    grp1nm = get_estr(&inp, "group1", "System");
    grp2nm = get_estr(&inp, "group2", "System");

    gmx::TextOutputFile outpara_stream(outpara_fn);
    write_inpfile(&outpara_stream, outpara_fn.c_str(), &inp, FALSE, WriteMdpHeader::yes, wi);
}

}
