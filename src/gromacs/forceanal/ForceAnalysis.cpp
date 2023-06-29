/*
    ForceAnalysis.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/09/19
    Description: Core module for Force Analysis.
*/

// C++ STL
#include <algorithm>

// GROMACS
#include "gromacs/linearalgebra/gmx_lapack.h"
#include "gromacs/linearalgebra/matrix.h"

// ForceAnal module
#include "ForceAnalysis.h"

ForceAnalysis::ForceAnalysis()
 : frame_count(0)
{
}

ForceAnalysis::ForceAnalysis(int nfile, const t_filenm fnm[], gmx_mtop_t *mtop, const t_inputrec *inputrec)
 : ForceAnal::ForceParaSet(nfile, fnm, mtop, inputrec),
   frame_count(0)
{
    init_outfiles();

    force_deviation.resizeWithPadding(atomn);

    switch (datamode)
    {
        case ForceAnal::DATA_MODE::SummedMode:
            summed_forces = ForceAnal::SummedData(grp1idx, grp2idx, threshold);
            break;
        case ForceAnal::DATA_MODE::DetailedMode:
            detailed_forces = ForceAnal::DetailedData(grp1idx, grp2idx, threshold);
            break;
        case ForceAnal::DATA_MODE::ListMode:
            listed_forces = ForceAnal::ListData(threshold);
            break;
        default:
            break;
    }
}

ForceAnalysis::~ForceAnalysis()
{
}

void ForceAnalysis::init_outfiles()
{
    // If output result files already exist, clean file contents    
    if (!res_txt_fn.empty())
    {
        std::ofstream txtfile(res_txt_fn, std::ios::out | std::ios::trunc);
        if (txtfile.is_open())
            txtfile.close();
    }
    if (!res_bin_fn.empty())
    {
        std::ofstream binfile(res_bin_fn, std::ios::out | std::ios::trunc);
        if (binfile.is_open())
        {
            uint8_t filecode = static_cast<uint8_t>(datamode);
            binfile.write((char*)&filecode, sizeof(uint8_t));
            // Instead of writing group index range data, using standalone atom map
            /*
            binfile.write((char*)&grp1idx.first, sizeof(ForceAnal::atomindex));
            binfile.write((char*)&grp1idx.second, sizeof(ForceAnal::atomindex));
            binfile.write((char*)&grp2idx.first, sizeof(ForceAnal::atomindex));
            binfile.write((char*)&grp2idx.second, sizeof(ForceAnal::atomindex));
            */
            binfile.close();
        }
    }
    if (!totf_bin_fn.empty())
    {
        std::ofstream totfile(totf_bin_fn, std::ios::out | std::ios::trunc);
        if (totfile.is_open())
        {
            uint8_t filecode = static_cast<uint8_t>(ForceAnal::DATA_MODE::AtomForceMode);
            totfile.write((char*)&filecode, sizeof(uint8_t));
            // Instead of writing group index range data, using standalone atom map
            /*
            totfile.write((char*)&grp1idx.first, sizeof(ForceAnal::atomindex));
            totfile.write((char*)&grp1idx.second, sizeof(ForceAnal::atomindex));
            */
            totfile.close();
        }
    }
    if (!atomf_nb_bin_fn.empty())
    {
        std::ofstream nbatomfile(atomf_nb_bin_fn, std::ios::out | std::ios::trunc);
        if (nbatomfile.is_open())
        {
            uint8_t filecode = static_cast<uint8_t>(ForceAnal::DATA_MODE::AtomForceMode);
            nbatomfile.write((char*)&filecode, sizeof(uint8_t));
            nbatomfile.close();
        }
    }
    if (!atomf_nb_b_bin_fn.empty())
    {
        std::ofstream nbbatomfile(atomf_nb_bin_fn, std::ios::out | std::ios::trunc);
        if (nbbatomfile.is_open())
        {
            uint8_t filecode = static_cast<uint8_t>(ForceAnal::DATA_MODE::AtomForceMode);
            nbbatomfile.write((char*)&filecode, sizeof(uint8_t));
            nbbatomfile.close();
        }
    }
    if (!fdev_bin_fn.empty())
    {
        std::ofstream devfile(fdev_bin_fn, std::ios::out | std::ios::trunc);
        if (devfile.is_open())
        {
            uint8_t filecode = static_cast<uint8_t>(ForceAnal::DATA_MODE::AtomForceMode);
            devfile.write((char*)&filecode, sizeof(uint8_t));
            // Instead of writing group index range data, using standalone atom map
            /*
            devfile.write((char*)&grp1idx.first, sizeof(ForceAnal::atomindex));
            devfile.write((char*)&grp1idx.second, sizeof(ForceAnal::atomindex));
            */
            devfile.close();
        }
    }

}

bool ForceAnalysis::atom_in_grp1(const int idx)
{
    return (idx >= grp1aid.first) && (idx < grp1aid.second);
}

bool ForceAnalysis::atom_in_grp2(const int idx)
{
    return (idx >= grp2aid.first) && (idx < grp2aid.second);
}

bool ForceAnalysis::in_grp1(const int idx)
{
    // For atom in group 1, its ID should be in range grp1_start <= aid < grp1_end
    bool inrange = (idx >= grp1idx.first) && (idx < grp1idx.second);
    if (exclgrp1.empty() || !inrange) return inrange;
    else
    {
        for (const int& ai : exclgrp1)
            if (ai == idx)
            {
                inrange = false;
                break;
            }
        return inrange;
    }
}

bool ForceAnalysis::in_grp2(const int idx)
{
    // For atom in group 2, its ID should be in range grp2_start <= aid < grp2_end
    bool inrange = (idx >= grp2idx.first) && (idx < grp2idx.second);
    if (exclgrp2.empty() || !inrange) return inrange;
    else
    {
        for (const int& ai : exclgrp2)
            if (ai == idx)
            {
                inrange = false;
                break;
            }
        return inrange;
    }
}

bool ForceAnalysis::in_grp(const int i, const int j)
{
    return in_grp1(i) && in_grp2(j);
}

void ForceAnalysis::mapping_index(int &i, int &j)
{
    switch (forceunit)
    {
        case ForceAnal::FORCE_UNIT::Atom:
            break;
        case ForceAnal::FORCE_UNIT::Residue:
            i = resmap[i];
            j = resmap[j];
            break;
        case ForceAnal::FORCE_UNIT::Molecule:
            i = molmap[i];
            j = molmap[j];
            break;
    }
}

void ForceAnalysis::add_pairforce(int i, int j, ForceAnal::InteractionType type, rvec f_ij)
{
    if (datamode == ForceAnal::DATA_MODE::None) return;
    switch (forceunit)
    {
        case ForceAnal::FORCE_UNIT::Atom:
            break;
        case ForceAnal::FORCE_UNIT::Residue:
            i = resmap[i];
            j = resmap[j];
            break;
        case ForceAnal::FORCE_UNIT::Molecule:
            i = molmap[i];
            j = molmap[j];
            break;
    }
    if (in_grp(i, j))
    {
        switch (datamode)
        {
            case ForceAnal::DATA_MODE::SummedMode:
                summed_forces.add_detailed_force(i, j, type, f_ij);
                break;
            case ForceAnal::DATA_MODE::DetailedMode:
                detailed_forces.add_detailed_force(i, j, type, f_ij);
                break;
            case ForceAnal::DATA_MODE::ListMode:
                listed_forces.add_detailed_force(i, j, type, f_ij);
                break;
            default:
                // Do Nothing
                break;
        }
    }
    rvec f_ji;
    rvec_opp(f_ij, f_ji);
    if (in_grp(j, i))
    {
        switch (datamode)
        {
            case ForceAnal::DATA_MODE::SummedMode:
                summed_forces.add_detailed_force(j, i, type, f_ji);
                break;
            case ForceAnal::DATA_MODE::DetailedMode:
                detailed_forces.add_detailed_force(j, i, type, f_ji);
                break;
            case ForceAnal::DATA_MODE::ListMode:
                listed_forces.add_detailed_force(j, i, type, f_ji);
                break;
            default:
                // Do Nothing
                break;
        }
    }
}

void ForceAnalysis::add_nonbonded(int i, int j, real pf_coul, real pf_vdw, real dx, real dy, real dz)
{
    add_nonbonded_coulomb(i, j, pf_coul, dx, dy, dz);
    add_nonbonded_vdw(i, j, pf_vdw, dx, dy, dz);
}

void ForceAnalysis::add_nonbonded_coulomb(int i, int j, real pf_coul, real dx, real dy, real dz)
{
    rvec coul_force;
    coul_force[XX] = dx * pf_coul;
    coul_force[YY] = dy * pf_coul;
    coul_force[ZZ] = dz * pf_coul;
    add_pairforce(i, j, ForceAnal::Interact_COULOMB, coul_force);
}

void ForceAnalysis::add_nonbonded_vdw(int i, int j, real pf_vdw, real dx, real dy, real dz)
{
    rvec lj_force;
    lj_force[XX] = dx * pf_vdw;
    lj_force[YY] = dy * pf_vdw;
    lj_force[ZZ] = dz * pf_vdw;
    add_pairforce(i, j, ForceAnal::Interact_VDW, lj_force);
}

int ForceAnalysis::dbres(rvec f_i, rvec r_ij, rvec r_ik, rvec f_ij, rvec f_ik)
{
    int n, m, nrhs, lda, ldb, info, *ipiv;
    double **R, *c;
    const char trans = 'N';
    m = lda = ldb = n = 2;
    nrhs = 1;
    snew(ipiv, n);
    snew(c, n);
    copy_rvec_to_2dvec(f_i, c);
    R = alloc_matrix(n, n);
    copy_rvec_to_2dvec(r_ij, R[XX]);
    copy_rvec_to_2dvec(r_ik, R[YY]);
    F77_FUNC(dgetrf, DGETRF) (&n, &m, R[0], &lda, ipiv, &info);
    if (info != 0) return info;
    F77_FUNC(dgetrs, DGETRS) (&trans, &n, &nrhs, R[0], &lda, ipiv, c, &ldb, &info);
    if (info != 0) return info;
    svmul(c[XX], r_ij, f_ij);
    svmul(c[YY], r_ik, f_ik);
    sfree(ipiv);
    sfree(c);
    return 0;
}

int ForceAnalysis::trires(rvec f_i, rvec r_ij, rvec r_ik, rvec r_il, rvec f_ij, rvec f_ik, rvec f_il)
{
    int n, m, nrhs, lda, ldb, info, *ipiv;
    double **R, *c;
    const char trans = 'N';
    m = lda = ldb = n = 3;
    nrhs = 1;
    snew(ipiv, n);
    snew(c, n);
    copy_rvec_to_3dvec(f_i, c);
    R = alloc_matrix(n, n);
    copy_rvec_to_3dvec(r_ij, R[XX]);
    copy_rvec_to_3dvec(r_ik, R[YY]);
    copy_rvec_to_3dvec(r_il, R[ZZ]);
    F77_FUNC(dgetrf, DGETRF) (&n, &m, R[0], &lda, ipiv, &info);
    if (info != 0) return info;
    F77_FUNC(dgetrs, DGETRS) (&trans, &n, &nrhs, R[0], &lda, ipiv, c, &ldb, &info);
    if (info != 0) return info;
    svmul(c[XX], r_ij, f_ij);
    svmul(c[YY], r_ik, f_ik);
    svmul(c[ZZ], r_il, f_il);
    sfree(ipiv);
    sfree(c);
    return 0;
}

void ForceAnalysis::add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k, rvec r_ij, rvec r_kj, rvec r_ik)
{
    /*
    Fi, Fj and Fk are dependent, Fi + Fj + Fk = 0.
    Therefore, not all 3 of Fi, Fj and Fk are required.
    In this case, only Fi and Fk are required to recover pairwise forces.
    */

    // Clear force component perpendicular to the plane formed by position vector r_ij and r_ik
    real rvnorm2;
    rvec r_ki, f_ij, f_ik, f_ji, f_jk, f_ki, f_kj;
    rvec rv, proj_f_i, proj_f_k, f_i_err, f_j_err, f_k_err;
    cprod(r_ij, r_ik, rv);
    rvnorm2 = norm2(rv);
    svmul(iprod(f_i, rv) / rvnorm2, rv, proj_f_i);
    svmul(iprod(f_k, rv) / rvnorm2, rv, proj_f_k);

    rvec_opp(r_ik, r_ki);
    if (dbres(proj_f_i, r_ij, r_ik, f_ij, f_ik) != 0)
    {
        clear_rvec(f_ij);
        clear_rvec(f_ik);
        gmx_warning("Failed to resolve force at atom %d on angle potential\n"
                    "(%d, %d, %d).", ai, ai, aj, ak);
    }
    if (dbres(proj_f_k, r_ki, r_kj, f_ki, f_kj) != 0)
    {
        clear_rvec(f_ki);
        clear_rvec(f_kj);
        gmx_warning("Failed to resolve force at atom %d on angle potential\n"
                    "(%d, %d, %d).", ak, ai, aj, ak);
    }
    
#ifdef FORCEANAL_DEBUG
    bool ferr = false;
    rvec f_i_s, f_k_s, f_ik_s;
    rvec_add(f_ij, f_ik, f_i_s);
    rvec_dec(f_i_s, f_i);
    if (!rvec_abs_small(f_i_s, force_threshold))
    {
        ferr = true;
        gmx_warning("Force resolution at atom %d on angle potential (%d, %d, %d)\n"
                    "shows a large deviation [%.3f, %.3f, %.3f].",
                    ai, ai, aj, ak, f_i_s[XX], f_i_s[YY], f_i_s[ZZ]);
    }
    rvec_add(f_ki, f_kj, f_k_s);
    rvec_dec(f_k_s, f_k);
    if (!rvec_abs_small(f_k_s, force_threshold))
    {
        ferr = true;
        gmx_warning("Force resolution at atom %d on angle potential (%d, %d, %d)\n"
                    "shows a large deviation [%.3f, %.3f, %.3f].",
                    ak, ai, aj, ak, f_k_s[XX], f_k_s[YY], f_k_s[ZZ]);
    }
    rvec_add(f_ik, f_ki, f_ik_s);
    if (!rvec_abs_small(f_ik_s, force_threshold))
    {
        ferr = true;
        gmx_warning("Pairwise forces resolved at atom %d and %d on angle\n"
                    "potential (%d, %d, %d) are asymmetric with a large\n"
                    "deviation [%.3f, %.3f, %.3f].",
                    ai, ak, ai, aj, ak, f_ik_s[XX], f_ik_s[YY], f_ik_s[ZZ]);
    }
    if (ferr)
    {
        gmx_warning("Angle Potential (%d, %d, %d)\n"
                    "Fi [%.3f, %.3f, %.3f], Fj [%.3f, %.3f, %.3f], Fk [%.3f, %.3f, %.3f]\n"
                    "Rij [%.3f, %.3f, %.3f], Rkj [%.3f, %.3f, %.3f], Rik [%.3f, %.3f, %.3f]",
                    ai, aj, ak,
                    f_i[XX], f_i[YY], f_i[ZZ], f_j[XX], f_j[YY], f_j[ZZ],
                    f_k[XX], f_k[YY], f_k[ZZ],
                    r_ij[XX], r_ij[YY], r_ij[ZZ], r_kj[XX], r_kj[YY], r_kj[ZZ],
                    r_ik[XX], r_ik[YY], r_ik[ZZ]);
    }
#endif

    add_pairforce(ai, aj, ForceAnal::Interact_ANGLE, f_ij);
    add_pairforce(ai, ak, ForceAnal::Interact_ANGLE, f_ik);
    add_pairforce(ak, aj, ForceAnal::Interact_ANGLE, f_kj);

    // Save force decomposition deviations to storage
    rvec_opp(f_ij, f_ji);
    rvec_opp(f_kj, f_jk);
    rvec_sub(f_i, f_ij, f_i_err);
    rvec_dec(f_i_err, f_ik);
    rvec_sub(f_j, f_ji, f_j_err);
    rvec_dec(f_j_err, f_jk);
    rvec_sub(f_k, f_ki, f_k_err);
    rvec_dec(f_k_err, f_kj);
    force_deviation[ai] += f_i_err;
    force_deviation[aj] += f_j_err;
    force_deviation[ak] += f_k_err;
}

void ForceAnalysis::add_dihedral(int ai, int aj, int ak, int al, rvec f_i, rvec f_j, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl)
{
    /*
    Fi, Fj, Fk and Fl are dependent, Fi + Fj + Fk + Fl = 0.
    Therefore, not all 4 of Fi, Fj, Fk and Fl are required.
    In this case, only Fi, Fj and Fk are required to recover pairwise forces.
    Note that f_j and f_k here indicate -Fj and -Fk, respectively.
    */
    rvec_opp(f_j);
    rvec_opp(f_k);

    rvec r_ji, r_jk, r_ki, r_ik, r_il, r_jl;
    rvec f_ij, f_ik, f_il, f_ji, f_jk, f_jl, f_ki, f_kj, f_kl, f_li, f_lj, f_lk;
    rvec f_i_err, f_j_err, f_k_err, f_l_err;
    rvec_opp(r_ij, r_ji);
    rvec_opp(r_kj, r_jk);
    rvec_sub(r_ij, r_kj, r_ik);
    rvec_opp(r_ik, r_ki);
    rvec_add(r_ik, r_kl, r_il);
    rvec_sub(r_kj, r_kl, r_jl);
    if (trires(f_i, r_ij, r_ik, r_il, f_ij, f_ik, f_il) != 0)
    {
        clear_rvec(f_ij);
        clear_rvec(f_ik);
        clear_rvec(f_il);
        gmx_warning("Failed to resolve force at atom %d on dihedral potential\n"
                    "(%d, %d, %d, %d).", ai, ai, aj, ak, al);
    }
    if (trires(f_j, r_ji, r_jk, r_jl, f_ji, f_jk, f_jl) != 0)
    {
        clear_rvec(f_ji);
        clear_rvec(f_jk);
        clear_rvec(f_jl);
        gmx_warning("Failed to resolve force at atom %d on dihedral potential\n"
                    "(%d, %d, %d, %d).", aj, ai, aj, ak, al);
    }
    if (trires(f_k, r_ki, r_kj, r_kl, f_ki, f_kj, f_kl) != 0)
    {
        clear_rvec(f_ki);
        clear_rvec(f_kj);
        clear_rvec(f_kl);
        gmx_warning("Failed to resolve force at atom %d on dihedral potential\n"
                    "(%d, %d, %d, %d).", ak, ai, aj, ak, al);
    }
    
#ifdef FORCEANAL_DEBUG
    bool ferr = false;
    rvec f_i_s, f_j_s, f_k_s, f_ij_s, f_ik_s, f_jk_s;
    rvec_add(f_ij, f_ik, f_i_s);
    rvec_inc(f_i_s, f_il);
    rvec_dec(f_i_s, f_i);
    if (!rvec_abs_small(f_i_s, force_threshold))
    {
        ferr = true;
        gmx_warning("Force resolution at atom %d on dihedral potential (%d, %d, %d, %d)\n"
                    "shows a large deviation [%.3f, %.3f, %.3f].",
                    ai, ai, aj, ak, al, f_i_s[XX], f_i_s[YY], f_i_s[ZZ]);
    }
    rvec_add(f_ji, f_jk, f_j_s);
    rvec_inc(f_j_s, f_jl);
    rvec_dec(f_j_s, f_j);
    if (!rvec_abs_small(f_j_s, force_threshold))
    {
        ferr = true;
        gmx_warning("Force resolution at atom %d on dihedral potential (%d, %d, %d, %d)\n"
                    "shows a large deviation [%.3f, %.3f, %.3f].",
                    aj, ai, aj, ak, al, f_j_s[XX], f_j_s[YY], f_j_s[ZZ]);
    }
    rvec_add(f_ki, f_kj, f_k_s);
    rvec_inc(f_k_s, f_kl);
    rvec_dec(f_k_s, f_k);
    if (!rvec_abs_small(f_k_s, force_threshold))
    {
        ferr = true;
        gmx_warning("Force resolution at atom %d on dihedral potential (%d, %d, %d, %d)\n"
                    "shows a large deviation [%.3f, %.3f, %.3f].",
                    ak, ai, aj, ak, al, f_k_s[XX], f_k_s[YY], f_k_s[ZZ]);
    }
    rvec_add(f_ij, f_ji, f_ij_s);
    if (!rvec_abs_small(f_ij_s, force_threshold))
    {
        ferr = true;
        gmx_warning("Pairwise forces resolved at atom %d and %d on dihedral\n"
                    "potential (%d, %d, %d, %d) are asymmetric with a large\n"
                    "deviation [%.3f, %.3f, %.3f].",
                    ai, aj, ai, aj, ak, al, f_ij_s[XX], f_ij_s[YY], f_ij_s[ZZ]);
    }
    rvec_add(f_ik, f_ki, f_ik_s);
    if (!rvec_abs_small(f_ik_s, force_threshold))
    {
        ferr = true;
        gmx_warning("Pairwise forces resolved at atom %d and %d on dihedral\n"
                    "potential (%d, %d, %d, %d) are asymmetric with a large\n"
                    "deviation [%.3f, %.3f, %.3f].",
                    ai, ak, ai, aj, ak, al, f_ik_s[XX], f_ik_s[YY], f_ik_s[ZZ]);
    }
    rvec_add(f_jk, f_kj, f_jk_s);
    if (!rvec_abs_small(f_jk_s, force_threshold))
    {
        ferr = true;
        gmx_warning("Pairwise forces resolved at atom %d and %d on dihedral\n"
                    "potential (%d, %d, %d, %d) are asymmetric with a large\n"
                    "deviation [%.3f, %.3f, %.3f].",
                    aj, ak, ai, aj, ak, al, f_jk_s[XX], f_jk_s[YY], f_jk_s[ZZ]);
    }
    if (ferr)
    {
        gmx_warning("Dihedral Potential (%d, %d, %d, %d)\n"
                    "Fi [%.3f, %.3f, %.3f], Fj [%.3f, %.3f, %.3f]\n"
                    "Fk [%.3f, %.3f, %.3f], Fl [%.3f, %.3f, %.3f]\n"
                    "Rij [%.3f, %.3f, %.3f], Rkj [%.3f, %.3f, %.3f], Rkl [%.3f, %.3f, %.3f]",
                    ai, aj, ak, al,
                    f_i[XX], f_i[YY], f_i[ZZ], f_j[XX], f_j[YY], f_j[ZZ],
                    f_k[XX], f_k[YY], f_k[ZZ], f_l[XX], f_l[YY], f_l[ZZ],
                    r_ij[XX], r_ij[YY], r_ij[ZZ], r_kj[XX], r_kj[YY], r_kj[ZZ],
                    r_kl[XX], r_kl[YY], r_kl[ZZ]);
    }
#endif

    add_pairforce(ai, aj, ForceAnal::Interact_DIHEDRAL, f_ij);
    add_pairforce(ai, al, ForceAnal::Interact_DIHEDRAL, f_il);
    add_pairforce(aj, ak, ForceAnal::Interact_DIHEDRAL, f_jk);
    add_pairforce(aj, al, ForceAnal::Interact_DIHEDRAL, f_jl);
    add_pairforce(ak, ai, ForceAnal::Interact_DIHEDRAL, f_ki);
    add_pairforce(ak, al, ForceAnal::Interact_DIHEDRAL, f_kl);

    // Save force decomposition deviations to storage
    rvec_opp(f_il, f_li);
    rvec_opp(f_jl, f_lj);
    rvec_opp(f_kl, f_lk);
    rvec_sub(f_i, f_ij, f_i_err);
    rvec_dec(f_i_err, f_ik);
    rvec_dec(f_i_err, f_il);
    rvec_sub(f_j, f_ji, f_j_err);
    rvec_dec(f_j_err, f_jk);
    rvec_dec(f_j_err, f_jl);
    rvec_sub(f_k, f_ki, f_k_err);
    rvec_dec(f_k_err, f_kj);
    rvec_dec(f_k_err, f_kl);
    rvec_sub(f_l, f_li, f_l_err);
    rvec_dec(f_l_err, f_lj);
    rvec_dec(f_l_err, f_lk);
    force_deviation[ai] += f_i_err;
    force_deviation[aj] += f_j_err;
    force_deviation[ak] += f_k_err;
    force_deviation[al] += f_l_err;
}

void ForceAnalysis::write_frame(bool write_last_frame)
{
    if (write_last_frame)
    {
        frame_count = 0;   // Reset frame counter
        return ;
    }

    if (datamode == ForceAnal::DATA_MODE::None) return;

    // Can not average forces in Listed Forces mode
    if (datamode == ForceAnal::DATA_MODE::ListMode)
        Naverage = 1;

    if ((frame_count % Naverage) == 0)
        write_forces();
}

void ForceAnalysis::write_forces()
{
    if (!res_bin_fn.empty())
    {
        std::ofstream binstream(res_bin_fn, std::ios::binary | std::ios::ate | std::ios::in);
        
        if (!binstream.is_open())
            gmx_fatal(FARGS, "GROMACS Force Analysis module can not write force data to file.\n");

        binstream.seekp(0, std::ios::end);
        uint32_t forces_count = 0;
        binstream.write((char*)&frame_count, sizeof(uint32_t));
        binstream.write((char*)&forces_count, sizeof(uint32_t));

        int64_t saddr = binstream.tellp();
        int64_t eaddr = 0;
        switch (datamode)
        {
            case ForceAnal::DATA_MODE::SummedMode:
                summed_forces.write_forces_bin(binstream, forces_count, saddr, eaddr);
                break;
            case ForceAnal::DATA_MODE::DetailedMode:
                detailed_forces.write_forces_bin(binstream, forces_count, saddr, eaddr);
                break;
            case ForceAnal::DATA_MODE::ListMode:
                listed_forces.write_forces_bin(binstream, forces_count, saddr, eaddr);
                break;
            default:
                // Do Nothing
                break;
        }

        // Get forces count and back to write at the start of the frame
        binstream.seekp(saddr - sizeof(uint32_t));
        binstream.write((char*)&forces_count, sizeof(uint32_t));
        binstream.seekp(eaddr);

        binstream.close();
    }

    if (!res_txt_fn.empty())
    {
        std::ofstream txtstream(res_txt_fn, std::ios::app);

        if (!txtstream.is_open())
            gmx_fatal(FARGS, "GROMACS Force Analysis module can not write force data to file.\n");
        
        txtstream << "START FRAME " << frame_count << std::endl;

        switch (datamode)
        {
            case ForceAnal::DATA_MODE::SummedMode:
                summed_forces.write_forces_txt(txtstream);
                break;
            case ForceAnal::DATA_MODE::DetailedMode:
                detailed_forces.write_forces_txt(txtstream);
                break;
            case ForceAnal::DATA_MODE::ListMode:
                listed_forces.write_forces_txt(txtstream);
                break;
            default:
                // Do Nothing
                break;
        }

        txtstream << "END FRAME " << frame_count << std::endl;

        txtstream.close();
    }

    switch (datamode)
    {
        case ForceAnal::DATA_MODE::SummedMode:
            summed_forces.clear();
            break;
        case ForceAnal::DATA_MODE::DetailedMode:
            detailed_forces.clear();
            break;
        case ForceAnal::DATA_MODE::ListMode:
            listed_forces.clear();
            break;
        default:
            // Do Nothing
            break;
    }
}

void ForceAnalysis::write_atom_forces(const char* fnm, const rvec* f)
{
    rvec fi;
    uint32_t forces_count = 0;
    uint64_t fvec_len = 4 * sizeof(real);   // Each force vector block has 4 real

    std::ofstream forcestream(fnm, std::ios::binary | std::ios::app);
    if (forcestream.is_open())
    {
        forcestream.write((char*)&frame_count, sizeof(uint32_t));
        if (forceunit == ForceAnal::FORCE_UNIT::Atom)
        {
            real* buff = new real[4 * atomn];
            uint32_t ai, idx;
            for (ai = 0; ai < atomn; ++ai)
            {
                if (!in_grp1(ai)) continue;
                idx = 4 * forces_count;
                copy_rvec(f[ai], fi);
                buff[idx + XX] = fi[XX];
                buff[idx + YY] = fi[YY];
                buff[idx + ZZ] = fi[ZZ];
                buff[idx + ForceAnal::XYZ] = norm(fi);
                ++forces_count;
            }
            forcestream.write((char*)&forces_count, sizeof(uint32_t));
            forcestream.write((char*)buff, forces_count * fvec_len);
            delete[] buff;
        }
        else
        {
            std::vector<ForceAnal::atomindex>* pmap = (forceunit == ForceAnal::FORCE_UNIT::Residue) ? &resmap : &molmap;
            uint32_t nres = (forceunit == ForceAnal::FORCE_UNIT::Residue) ? resn : moln;
            uint32_t* rescount = new uint32_t[nres];
            uint32_t forcelen = 4 * nres;
            real* resforces = new real[forcelen];
            real* buff = new real[forcelen];
            uint32_t resi, idx;
            for (resi = 0; resi < nres; ++resi) rescount[resi] = 0;
            for (idx = 0; idx < forcelen; ++idx) resforces[idx] = 0.;
            for (uint32_t ai = 0; ai < atomn; ++ai)
            {
                resi = pmap->at(ai);
                if (!in_grp1(resi)) continue;
                idx = 4 * resi;
                copy_rvec(f[ai], fi);
                resforces[idx + XX] += fi[XX];
                resforces[idx + YY] += fi[YY];
                resforces[idx + ZZ] += fi[ZZ];
                ++rescount[resi];
            }
            for (resi = 0; resi < nres; ++resi)
            {
                if (rescount[resi] == 0) continue;
                idx = 4 * resi;
                fi[XX] = resforces[idx + XX];
                fi[YY] = resforces[idx + YY];
                fi[ZZ] = resforces[idx + ZZ];
                idx = 4 * forces_count;
                buff[idx + XX] = fi[XX];
                buff[idx + YY] = fi[YY];
                buff[idx + ZZ] = fi[ZZ];                
                buff[idx + ForceAnal::XYZ] = norm(fi);
                ++forces_count;
            }
            forcestream.write((char*)&forces_count, sizeof(uint32_t));
            forcestream.write((char*)buff, forces_count * fvec_len);
            delete[] rescount;
            delete[] resforces;
            delete[] buff;
            pmap = nullptr;
        }
        forcestream.close();
    }
}

void ForceAnalysis::write_atom_forces(ForceAnal::OUT_FORCE_TYPE out_ftype, rvec *f)
{
    if (out_ftype == ForceAnal::OUT_FORCE_TYPE::AtomForceTotal && !totf_bin_fn.empty())
        write_atom_forces(totf_bin_fn.c_str(), f);
    else if (out_ftype == ForceAnal::OUT_FORCE_TYPE::AtomForceNonbonded && !atomf_nb_bin_fn.empty())
        write_atom_forces(atomf_nb_bin_fn.c_str(), f);
    else if (out_ftype == ForceAnal::OUT_FORCE_TYPE::AtomForceNonbondedBonded && !atomf_nb_b_bin_fn.empty())
        write_atom_forces(atomf_nb_b_bin_fn.c_str(), f);
}

void ForceAnalysis::write_dev_forces(bool clearvec)
{
    if (!fdev_bin_fn.empty())
        write_atom_forces(fdev_bin_fn.c_str(), as_rvec_array(force_deviation.data()));
    if (clearvec)
    {
        force_deviation.clear();
        force_deviation.resizeWithPadding(atomn);
    }
}

void FA_add_nonbonded(class ForceAnalysis *FA, int i, int j, real pf_coul, real pf_vdw, real dx, real dy, real dz)
{
    FA->add_nonbonded(i, j, pf_coul, pf_vdw, dx, dy, dz);
}

void FA_add_nonbonded_coulomb(class ForceAnalysis *FA, int i, int j, real pf_coul, real dx, real dy, real dz)
{
    FA->add_nonbonded_coulomb(i, j, pf_coul, dx, dy, dz);
}

void FA_add_nonbonded_vdw(class ForceAnalysis *FA, int i, int j, real pf_vdw, real dx, real dy, real dz)
{
    FA->add_nonbonded_vdw(i, j, pf_vdw, dx, dy, dz);
}
