/*
 * ForceParaSet.cpp
 *
 *  Created on: Apr 14, 2021
 *      Author: Pengbo Song
 */

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/utility/filestream.h"

#include "ForceParaSet.h"

namespace ForceAnal {

ForceParaSet::ForceParaSet()
 : outpara_fn("faout.par"),
   summed_mode(true),
   output_type(OUT_NOTHING),
   threshold(1e-3),
   Naverage(1),
   group1_id(0),
   group2_id(0),
   group1_sid(1),
   group1_eid(1),
   group2_sid(1),
   group2_eid(1)
{
}

ForceParaSet::ForceParaSet(int nfile, const t_filenm fnm[], gmx_mtop_t *mtop)
 : outpara_fn("faout.par"),
   summed_mode(true),
   output_type(OUT_NOTHING),
   threshold(1e-3),
   Naverage(1),
   group1_id(0),
   group2_id(0)
{
    set_parameters(nfile, fnm);
    
    res_bin_fn = handle_empty_string(opt2fn("-fo", nfile, fnm));
    res_txt_fn = handle_empty_string(opt2fn("-ft", nfile, fnm));
}

ForceParaSet::~ForceParaSet()
{
}

bool ForceParaSet::checkterm2bool(const std::string term, bool def)
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

void ForceParaSet::set_parameters(int nfile, const t_filenm fnm[])
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

    summed_term = get_estr(&inp, "summed", "yes");
    vector_term = get_estr(&inp, "vector", "yes");
    scalar_term = get_estr(&inp, "scalar", "yes");

    summed_mode = checkterm2bool(summed_term);
    if (checkterm2bool(vector_term))
        output_type |= OUT_VECTOR;
    if (checkterm2bool(scalar_term))
        output_type |= OUT_SCALAR;
    
    threshold = get_ereal(&inp, "threshold", 1e-6, wi);
    Naverage = get_eint64(&inp, "naverage", 1, wi);
    group1_sid = get_eint64(&inp, "grp1_start", 1, wi);
    group1_eid = get_eint64(&inp, "grp1_end", 1, wi);
    group2_sid = get_eint64(&inp, "grp2_start", 1, wi);
    group2_eid = get_eint64(&inp, "grp2_end", 1, wi);

    gmx::TextOutputFile outpara_stream(outpara_fn);
    write_inpfile(&outpara_stream, outpara_fn.c_str(), &inp, FALSE, WriteMdpHeader::yes, wi);
}

}
