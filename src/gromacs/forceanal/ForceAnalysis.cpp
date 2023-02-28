/*
    ForceAnalysis.cpp
    Author: Pengbo Song [pbsong-ccme2019@pku.edu.cn]
    Date created: 2020/09/19
    Description: Core module for Force Analysis.
*/

// C++ STL
#include <algorithm>

// ForceAnal module
#include "ForceAnalysis.h"

ForceAnalysis::ForceAnalysis()
 : frame_count(0)
{
}

ForceAnalysis::ForceAnalysis(int nfile, const t_filenm fnm[], gmx_mtop_t *mtop)
 : ForceAnal::ForceParaSet(nfile, fnm, mtop),
   frame_count(0)
{
    init_outfiles();

    switch (datamode)
    {
        case ForceAnal::DATA_MODE::SummedMode:
            summed_forces = ForceAnal::SummedData(grp1idx, grp2idx, threshold, Naverage);
            break;
        case ForceAnal::DATA_MODE::DetailedMode:
            detailed_forces = ForceAnal::DetailedData(grp1idx, grp2idx, threshold, Naverage);
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
            totfile.close();
        }
    }

#ifdef FORCEANAL_DEBUG
    if (!totf_bin_fn.empty())
    {
        std::ofstream nb_totfile("internal_nb_forces.far", std::ios::out | std::ios::trunc);
        if (nb_totfile.is_open())
        {
            uint8_t filecode = static_cast<uint8_t>(ForceAnal::DATA_MODE::AtomForceMode);
            nb_totfile.write((char*)&filecode, sizeof(uint8_t));
            nb_totfile.close();
        }
    }
    if (!totf_bin_fn.empty())
    {
        std::ofstream nb_b_totfile("internal_nb+b_forces.far", std::ios::out | std::ios::trunc);
        if (nb_b_totfile.is_open())
        {
            uint8_t filecode = static_cast<uint8_t>(ForceAnal::DATA_MODE::AtomForceMode);
            nb_b_totfile.write((char*)&filecode, sizeof(uint8_t));
            nb_b_totfile.close();
        }
    }
#endif
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

bool ForceAnalysis::in_grp(int& i, int& j, rvec f_ij)
{
    bool ij_inrange = in_grp1(i) && in_grp2(j);
    bool ji_inrange = in_grp1(j) && in_grp2(i);
    if (ji_inrange && !ij_inrange)
    {
        std::swap<int>(i, j);
        rvec_opp(f_ij);
    }
    return ij_inrange || ji_inrange;
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
    if (!in_grp(i, j, f_ij)) return;
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

void ForceAnalysis::add_nonbonded(int i, int j, real pf_coul, real pf_vdw, real dx, real dy, real dz)
{
    add_nonbonded_coulomb(i, j, pf_coul, dx, dy, dz);
    add_nonbonded_vdw(i, j, pf_vdw, dx, dy, dz);
}

void ForceAnalysis::add_nonbonded_coulomb(int i, int j, real pf_coul, real dx, real dy, real dz)
{
    rvec coul_force;
    coul_force[0] = dx * pf_coul;
    coul_force[1] = dy * pf_coul;
    coul_force[2] = dz * pf_coul;
    add_pairforce(i, j, ForceAnal::Interact_COULOMB, coul_force);
}

void ForceAnalysis::add_nonbonded_vdw(int i, int j, real pf_vdw, real dx, real dy, real dz)
{
    rvec lj_force;
    lj_force[0] = dx * pf_vdw;
    lj_force[1] = dy * pf_vdw;
    lj_force[2] = dz * pf_vdw;
    add_pairforce(i, j, ForceAnal::Interact_VDW, lj_force);
}

void ForceAnalysis::tri_res(rvec f_i, rvec r_ij, rvec r_ik, rvec f_ij, rvec f_ik)
{
    /*
    FLoating point calculations: iprod for 5, cprod for 9, norm for 15, svmul for 3
    Total: 80
    */
    rvec cpij, cpik, cpjk;
    cprod(f_i, r_ij, cpij);
    cprod(f_i, r_ik, cpik);
    cprod(r_ij, r_ik, cpjk);
    real ncpij = norm(cpij);
    real ncpik = norm(cpik);
    real ncpjk = norm(cpjk);
    ncpjk += 1.0E-6;   // Used to avoid division by zero
    svmul(ncpik / ncpjk, r_ij, f_ij);
    svmul(ncpij / ncpjk, r_ik, f_ik);
    
    /*
    Floating point calculations: rvec_add/rvec_sub/rvec_opp for 3, rvec_abs_small for 6
    Total: 12 (best) 54 (worst)
    */
    rvec f_i_a, f_i_pp;
    rvec_add(f_ij, f_ik, f_i_a);
    rvec_sub(f_i_a, f_i, f_i_pp);
    if (rvec_abs_small(f_i_pp, 1.0E-6))
        return;
    rvec f_i_nn;
    rvec_add(f_i_a, f_i, f_i_nn);
    if (rvec_abs_small(f_i_nn, 1.0E-6))
    {
        rvec_opp(f_ij);
        rvec_opp(f_ik);
        return;
    }
    rvec f_i_b, f_i_pn;
    rvec_sub(f_ij, f_ik, f_i_b);
    rvec_sub(f_i_b, f_i, f_i_pn);
    if (rvec_abs_small(f_i_pn, 1.0E-6))
    {
        rvec_opp(f_ik);
        return;
    }
    rvec f_i_np;
    rvec_add(f_i_b, f_i, f_i_np);
    if (rvec_abs_small(f_i_np, 1.0E-6))
    {
        rvec_opp(f_ij);
        return;
    }
    gmx_warning("Force resolution failed to find a valid solution with a force\n"
                "vector [%.3f, %.3f, %.3f] and displacement vectors\n"
                "[%.4f, %.4f, %.4f] and [%.4f, %.4f, %.4f]",
                f_i[XX], f_i[YY], f_i[ZZ], r_ij[XX], r_ij[YY], r_ij[ZZ],
                r_ik[XX], r_ik[YY], r_ik[ZZ]);
}

void ForceAnalysis::add_angle(int ai, int aj, int ak, rvec f_i, rvec gmx_unused f_j, rvec f_k, rvec r_ij, rvec r_kj, rvec r_ik)
{
    /*
    Fi, Fj and Fk are dependent, Fi + Fj + Fk = 0.
    Therefore, not all 3 of Fi, Fj and Fk are required.
    In this case, only Fi and Fk are required to recover pairwise forces.
    */
    rvec r_ki, f_ij, f_ik, f_ki, f_kj;
    rvec_opp(r_ik, r_ki);
    tri_res(f_i, r_ij, r_ik, f_ij, f_ik);
    tri_res(f_k, r_ki, r_kj, f_ki, f_kj);
    rvec f_ik_s;
    rvec_add(f_ik, f_ki, f_ik_s);
    if (rvec_abs_small(f_ik_s, 1.0E-6))
        gmx_warning("Angle potential %d-%d-%d force resolution failed. Pairwise\n"
                    "forces are asymmetric.", ai, aj, ak);
    add_pairforce(ai, aj, ForceAnal::Interact_ANGLE, f_ij);
    add_pairforce(ai, ak, ForceAnal::Interact_ANGLE, f_ik);
    add_pairforce(ak, aj, ForceAnal::Interact_ANGLE, f_kj);
}

void ForceAnalysis::add_dihedral(int ai, int aj, int ak, int al, rvec f_i, rvec f_j, rvec f_k, rvec gmx_unused f_l, rvec r_ij, rvec r_kj, rvec r_kl)
{
    /*
    Fi, Fj, Fk and Fl are dependent, Fi + Fj + Fk + Fl = 0.
    Therefore, not all 4 of Fi, Fj, Fk and Fl are required.
    In this case, only Fi, Fj and Fk are required to recover pairwise forces.
    */
    rvec r_jk, r_ki, r_ik, r_il, r_jl;
    rvec f_ij, f_il, f_jk, f_jl, f_ki, f_kl;
    rvec_opp(r_kj, r_jk);
    rvec_sub(r_ij, r_kj, r_ik);
    rvec_opp(r_ik, r_ki);
    rvec_add(r_ik, r_kl, r_il);
    rvec_sub(r_kj, r_kl, r_jl);
    tri_res(f_i, r_ij, r_il, f_ij, f_il);
    tri_res(f_j, r_jk, r_jl, f_jk, f_jl);
    tri_res(f_k, r_ki, r_kl, f_ki, f_kl);
    add_pairforce(ai, aj, ForceAnal::Interact_ANGLE, f_ij);
    add_pairforce(ai, al, ForceAnal::Interact_ANGLE, f_il);
    add_pairforce(aj, ak, ForceAnal::Interact_ANGLE, f_jk);
    add_pairforce(aj, al, ForceAnal::Interact_ANGLE, f_jl);
    add_pairforce(ak, ai, ForceAnal::Interact_ANGLE, f_ki);
    add_pairforce(ak, al, ForceAnal::Interact_ANGLE, f_kl);
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
    switch (datamode)
    {
        case ForceAnal::DATA_MODE::SummedMode:
            summed_forces.average_forces();
            break;
        case ForceAnal::DATA_MODE::DetailedMode:
            detailed_forces.average_forces();
            break;
        default:
            // Do Nothing
            break;
    }

    if (!res_bin_fn.empty())
    {
        std::ofstream binstream(res_bin_fn, std::ios::binary | std::ios::ate | std::ios::in);
        
        if (!binstream.is_open())
            gmx_fatal(FARGS, "GROMACS Force Analysis module can not write force data to file.\n");

        binstream.seekp(0, std::ios::end);
        uint32_t forces_count = 0;
        binstream.write((char*)&frame_count, sizeof(uint32_t));
        binstream.write((char*)&forces_count, sizeof(uint32_t));

        uint64_t saddr = binstream.tellp();
        uint64_t eaddr = 0;
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
    real force;

    std::ofstream forcestream(fnm, std::ios::binary | std::ios::app);
    if (forcestream.is_open())
    {
        forcestream.write((char*)&frame_count, sizeof(uint32_t));
        if (forceunit == ForceAnal::FORCE_UNIT::Atom)
        {
            forcestream.write((char*)&atomn, sizeof(uint32_t));
            for (uint32_t i = 0; i < atomn; ++i)
            {
                copy_rvec(f[i], fi);
                force = norm(fi);
                forcestream.write((char*)&fi, sizeof(real) * 3);
                forcestream.write((char*)&force, sizeof(real));
            }
        }
        else
        {
            std::vector<ForceAnal::atomindex>* pmap = (forceunit == ForceAnal::FORCE_UNIT::Residue) ? &resmap : &molmap;
            uint32_t nres = (forceunit == ForceAnal::FORCE_UNIT::Residue) ? resn : moln;
            forcestream.write((char*)&nres, sizeof(uint32_t));
            uint32_t forcelen = 4 * nres;
            real* resforces = new real[forcelen];
            uint32_t idx;
            for (idx = 0; idx < forcelen; ++idx) resforces[idx] = 0.;
            for (uint32_t ai = 0; ai < atomn; ++ai)
            {
                idx = 4 * pmap->at(ai);
                copy_rvec(f[ai], fi);
                resforces[idx + XX] += fi[XX];
                resforces[idx + YY] += fi[YY];
                resforces[idx + ZZ] += fi[ZZ];
            }
            real fx, fy, fz;
            for (idx = 0; idx < forcelen; )
            {
                fx = resforces[idx++];
                fy = resforces[idx++];
                fz = resforces[idx++];
                resforces[idx++] = ForceAnal::real_norm(fx, fy, fz);
            }
            forcestream.write((char*)resforces, sizeof(real) * forcelen);
            delete[] resforces;
            pmap = nullptr;
        }
        forcestream.close();
    }
}

void ForceAnalysis::write_tot_forces(const rvec* f)
{
    if (!totf_bin_fn.empty())
        write_atom_forces(totf_bin_fn.c_str(), f);
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
