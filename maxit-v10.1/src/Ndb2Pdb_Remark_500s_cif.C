/*
FILE:     Ndb2Pdb_Remark_500s_cif.C
*/
/*
VERSION:  11.200
*/
/*
DATE:     4/20/2024
*/
/*
  Comments and Questions to: sw-help@rcsb.rutgers.edu
*/
/*
COPYRIGHT 1999-2024 Rutgers - The State University of New Jersey

This software is provided WITHOUT WARRANTY OF MERCHANTABILITY OR
FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
IMPLIED.  RUTGERS MAKE NO REPRESENTATION OR WARRANTY THAT THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT OR OTHER
PROPRIETARY RIGHT.

The user of this software shall indemnify, hold harmless and defend
Rutgers, its governors, trustees, officers, employees, students,
agents and the authors against any and all claims, suits,
losses, liabilities, damages, costs, fees, and expenses including
reasonable attorneys' fees resulting from or arising out of the
use of this software.  This indemnification shall include, but is
not limited to, any and all claims alleging products liability.
*/
/*
               RCSB PDB SOFTWARE LICENSE AGREEMENT

BY CLICKING THE ACCEPTANCE BUTTON OR INSTALLING OR USING 
THIS "SOFTWARE, THE INDIVIDUAL OR ENTITY LICENSING THE  
SOFTWARE ("LICENSEE") IS CONSENTING TO BE BOUND BY AND IS 
BECOMING A PARTY TO THIS AGREEMENT.  IF LICENSEE DOES NOT 
AGREE TO ALL OF THE TERMS OF THIS AGREEMENT
THE LICENSEE MUST NOT INSTALL OR USE THE SOFTWARE.

1. LICENSE AGREEMENT

This is a license between you ("Licensee") and the Protein Data Bank (PDB) 
at Rutgers, The State University of New Jersey (hereafter referred to 
as "RUTGERS").   The software is owned by RUTGERS and protected by 
copyright laws, and some elements are protected by laws governing 
trademarks, trade dress and trade secrets, and may be protected by 
patent laws. 

2. LICENSE GRANT

RUTGERS grants you, and you hereby accept, non-exclusive, royalty-free 
perpetual license to install, use, modify, prepare derivative works, 
incorporate into other computer software, and distribute in binary 
and source code format, or any derivative work thereof, together with 
any associated media, printed materials, and on-line or electronic 
documentation (if any) provided by RUTGERS (collectively, the "SOFTWARE"), 
subject to the following terms and conditions: (i) any distribution 
of the SOFTWARE shall bind the receiver to the terms and conditions 
of this Agreement; (ii) any distribution of the SOFTWARE in modified 
form shall clearly state that the SOFTWARE has been modified from 
the version originally obtained from RUTGERS.  

2. COPYRIGHT; RETENTION OF RIGHTS.  

The above license grant is conditioned on the following: (i) you must 
reproduce all copyright notices and other proprietary notices on any 
copies of the SOFTWARE and you must not remove such notices; (ii) in 
the event you compile the SOFTWARE, you will include the copyright 
notice with the binary in such a manner as to allow it to be easily 
viewable; (iii) if you incorporate the SOFTWARE into other code, you 
must provide notice that the code contains the SOFTWARE and include 
a copy of the copyright notices and other proprietary notices.  All 
copies of the SOFTWARE shall be subject to the terms of this Agreement.  

3. NO MAINTENANCE OR SUPPORT; TREATMENT OF ENHANCEMENTS 

RUTGERS is under no obligation whatsoever to: (i) provide maintenance 
or support for the SOFTWARE; or (ii) to notify you of bug fixes, patches, 
or upgrades to the features, functionality or performance of the 
SOFTWARE ("Enhancements") (if any), whether developed by RUTGERS 
or third parties.  If, in its sole discretion, RUTGERS makes an 
Enhancement available to you and RUTGERS does not separately enter 
into a written license agreement with you relating to such bug fix, 
patch or upgrade, then it shall be deemed incorporated into the SOFTWARE 
and subject to this Agreement. You are under no obligation whatsoever 
to provide any Enhancements to RUTGERS or the public that you may 
develop over time; however, if you choose to provide your Enhancements 
to RUTGERS, or if you choose to otherwise publish or distribute your 
Enhancements, in source code form without contemporaneously requiring 
end users or RUTGERS to enter into a separate written license agreement 
for such Enhancements, then you hereby grant RUTGERS a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare
derivative works, incorporate into the SOFTWARE or other computer
software, distribute, and sublicense your Enhancements or derivative
works thereof, in binary and source code form.

4. FEES.  There is no license fee for the SOFTWARE.  If Licensee
wishes to receive the SOFTWARE on media, there may be a small charge
for the media and for shipping and handling.  Licensee is
responsible for any and all taxes.

5. TERMINATION.  Without prejudice to any other rights, Licensor
may terminate this Agreement if Licensee breaches any of its terms
and conditions.  Upon termination, Licensee shall destroy all
copies of the SOFTWARE.

6. PROPRIETARY RIGHTS.  Title, ownership rights, and intellectual
property rights in the Product shall remain with RUTGERS.  Licensee 
acknowledges such ownership and intellectual property rights and will 
not take any action to jeopardize, limit or interfere in any manner 
with RUTGERS' ownership of or rights with respect to the SOFTWARE.  
The SOFTWARE is protected by copyright and other intellectual 
property laws and by international treaties.  Title and related 
rights in the content accessed through the SOFTWARE is the property 
of the applicable content owner and is protected by applicable law.  
The license granted under this Agreement gives Licensee no rights to such
content.

7. DISCLAIMER OF WARRANTY.  THE SOFTWARE IS PROVIDED FREE OF 
CHARGE, AND, THEREFORE, ON AN "AS IS" BASIS, WITHOUT WARRANTY OF 
ANY KIND, INCLUDING WITHOUT LIMITATION THE WARRANTIES THAT IT 
IS FREE OF DEFECTS, MERCHANTABLE, FIT FOR A PARTICULAR PURPOSE 
OR NON-INFRINGING.  THE ENTIRE RISK AS TO THE QUALITY AND 
PERFORMANCE OF THE SOFTWARE IS BORNE BY LICENSEE.  SHOULD THE 
SOFTWARE PROVE DEFECTIVE IN ANY RESPECT, THE LICENSEE AND NOT 
LICENSOR ASSUMES THE ENTIRE COST OF ANY SERVICE AND REPAIR.  
THIS DISCLAIMER OF WARRANTY CONSTITUTES AN ESSENTIAL PART OF 
THIS AGREEMENT.  NO USE OF THE PRODUCT IS AUTHORIZED HEREUNDER 
EXCEPT UNDER THIS DISCLAIMER.

8. LIMITATION OF LIABILITY.  TO THE MAXIMUM EXTENT PERMITTED BY
APPLICABLE LAW,  IN NO EVENT WILL LICENSOR BE LIABLE FOR ANY 
INDIRECT, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING 
OUT OF THE USE OF OR INABILITY TO USE THE SOFTWARE, INCLUDING, 
WITHOUT LIMITATION, DAMAGES FOR LOSS OF GOODWILL, WORK 
STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, OR ANY AND ALL 
OTHER COMMERCIAL DAMAGES OR LOSSES, EVEN IF ADVISED OF THE
POSSIBILITY THEREOF. 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Maxit.h"
#include "NdbToken.h"
#include "Remark_define.h"
#include "Remark_extern.h"
#include "utillib.h"

void Maxit::_ndb_to_pdb_get_remark_500()
{
       if (_molecules.empty()) return;
       if (!_CifObj) return;
       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       _ndb_to_pdb_get_remark_500_asym_contact(_cifblock);
       _ndb_to_pdb_get_remark_500_symm_contact(_cifblock);
       _ndb_to_pdb_get_remark_500_bond_deviation(_cifblock);
       _ndb_to_pdb_get_remark_500_angle_deviation(_cifblock);
       _ndb_to_pdb_get_remark_500_Ramachandran_outliers(_cifblock);
       _ndb_to_pdb_get_remark_500_non_cis_trans_torsions(_cifblock);
       _ndb_to_pdb_get_remark_500_side_chain_plane(_cifblock);
       _ndb_to_pdb_get_remark_500_main_chain_planarity(_cifblock);
       // alpha_carbon_chiral();
}

void Maxit::_ndb_to_pdb_get_remark_500_asym_contact(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_validate_close_contact");
       if (!t) return;

       if (_experiment_type & EXPERIMENT_TYPE_MODEL || _experiment_type & EXPERIMENT_TYPE_NMR_SOLID || _experiment_type & EXPERIMENT_TYPE_NMR)
            _ndb_to_pdb_get_general_remark(Num_Remark_500_NMR, Remark_500_NMR, 0, 0, 0);
       else _ndb_to_pdb_get_general_remark(Num_Remark_500_A, Remark_500_A, 0, 0, 0);

       _ndb_to_pdb_get_remark_500_close_contact(t, false); 
}

void Maxit::_ndb_to_pdb_get_remark_500_symm_contact(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_validate_symm_contact");
       if (!t) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_500_S, Remark_500_S, 0, 0, 0);

       _ndb_to_pdb_get_remark_500_close_contact(t, true);
}

void Maxit::_ndb_to_pdb_get_remark_500_close_contact(ISTable *t, const bool& sym_flag)
{
       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();
       data.clear();
       data.push_back("auth_asym_id_1");
       data.push_back("auth_comp_id_1");
       data.push_back("auth_seq_id_1");
       data.push_back("PDB_ins_code_1");
       data.push_back("auth_atom_id_1");
       data.push_back("label_alt_id_1");
       items.push_back(data);
       data.clear();
       data.push_back("auth_asym_id_2");
       data.push_back("auth_comp_id_2");
       data.push_back("auth_seq_id_2");
       data.push_back("PDB_ins_code_2");
       data.push_back("auth_atom_id_2");
       data.push_back("label_alt_id_2");
       items.push_back(data);

       std::string cs, symmetry;
       int rowNo = t->GetNumRows();
       int count = 0, mol_index = 0;
       std::vector<RCSB::Atom*> atoms;
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs, t, i, "PDB_model_num");
            if (atoi(cs.c_str()) != _molecules[0]->Mol_ID()) continue;
            if (count >= 50) {
                 cs = "THIS ENTRY HAS " + FloatToString((double) rowNo, 7, 0, false, true);
                 if (sym_flag)
                      cs += " SYMMETRY CONTACTS";
                 else cs += " CLOSE CONTACTS";
                 _addNewRemark(500, "");
                 _addNewRemark(500, cs);
                 break;
            }

            get_values(t, i, items, names);
            _getAtomsfromNames(0, names, mol_index, atoms, data); 
            if (atoms.size() < 2) continue;
            count++;

            get_value(cs, t, i, "dist");
            double dist = atof(cs.c_str());

            symmetry.clear();
            if (sym_flag) {
                 get_value(symmetry, t, i, "site_symmetry_2");
                 std::string::size_type p = symmetry.find("_");
                 if (p != std::string::npos) symmetry.erase(p, 1);
            }

            _remark_500_close_contact(atoms, dist, symmetry);
       }
       _addNewRemark(500, "");
       _addNewRemark(500, "REMARK: NULL");
}

void Maxit::_ndb_to_pdb_get_remark_500_bond_deviation(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_validate_rmsd_bond");
       if (!t) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();
       data.clear();
       data.push_back("auth_asym_id_1");
       data.push_back("auth_comp_id_1");
       data.push_back("auth_seq_id_1");
       data.push_back("PDB_ins_code_1");
       data.push_back("auth_atom_id_1");
       data.push_back("label_alt_id_1");
       items.push_back(data);
       data.clear();
       data.push_back("auth_asym_id_2");
       data.push_back("auth_comp_id_2");
       data.push_back("auth_seq_id_2");
       data.push_back("PDB_ins_code_2");
       data.push_back("auth_atom_id_2");
       data.push_back("label_alt_id_2");
       items.push_back(data);

       _ndb_to_pdb_get_general_remark(Num_Remark_500_BOND, Remark_500_BOND, 0, 0, 0);

       std::string cs;
       int count = 0, mol_index = 0;
       std::vector<RCSB::Atom*> atoms;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            if (count >= 50) {
                 std::string remark = "THIS ENTRY HAS "
                     + FloatToString((double) rowNo, 7, 0, false, true)
                     + " BOND DEVIATIONS.";
                 _addNewRemark(500, "");
                 _addNewRemark(500, remark);
                 break;
            }

            get_value(cs, t, i, "PDB_model_num");
            int mol_id = atoi(cs.c_str());
            get_values(t, i, items, names);
            _getAtomsfromNames(mol_id, names, mol_index, atoms, data); 
            if (atoms.size() < 2) continue;
            count++;

            get_value(cs, t, i, "bond_deviation");
            double val = atof(cs.c_str());
            if (_molecules.size() < 2) mol_id = 0;
            _remark_500_bond_deviation(atoms, val, mol_id);
       }
       _addNewRemark(500, "");
       _addNewRemark(500, "REMARK: NULL");
}

void Maxit::_ndb_to_pdb_get_remark_500_angle_deviation(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_validate_rmsd_angle");
       if (!t) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_500_ANGLE, Remark_500_ANGLE, 0, 0, 0);
       if (!t) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();
       data.clear();
       data.push_back("auth_asym_id_1");
       data.push_back("auth_comp_id_1");
       data.push_back("auth_seq_id_1");
       data.push_back("PDB_ins_code_1");
       data.push_back("auth_atom_id_1");
       data.push_back("label_alt_id_1");
       items.push_back(data);
       data.clear();
       data.push_back("auth_asym_id_2");
       data.push_back("auth_comp_id_2");
       data.push_back("auth_seq_id_2");
       data.push_back("PDB_ins_code_2");
       data.push_back("auth_atom_id_2");
       data.push_back("label_alt_id_2");
       items.push_back(data);
       data.clear();
       data.push_back("auth_asym_id_3");
       data.push_back("auth_comp_id_3");
       data.push_back("auth_seq_id_3");
       data.push_back("PDB_ins_code_3");
       data.push_back("auth_atom_id_3");
       data.push_back("label_alt_id_3");
       items.push_back(data);

       std::string cs;
       int count = 0, mol_index = 0;
       std::vector<RCSB::Atom*> atoms;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            if (count >= 50) {
                 std::string remark = "THIS ENTRY HAS "
                     + FloatToString((double) rowNo, 7, 0, false, true)
                     + " ANGLE DEVIATIONS.";
                 _addNewRemark(500, "");
                 _addNewRemark(500, remark);
                 break;
            }

            get_value(cs, t, i, "PDB_model_num");
            int mol_id = atoi(cs.c_str());
            get_values(t, i, items, names);
            _getAtomsfromNames(mol_id, names, mol_index, atoms, data); 
            if (atoms.size() < 3) continue;
            count++;

            get_value(cs, t, i, "angle_deviation");
            double val = atof(cs.c_str());
            if (_molecules.size() < 2) mol_id = 0;
            _remark_500_angle_deviation(atoms, val, mol_id);
       }
       _addNewRemark(500, "");
       _addNewRemark(500, "REMARK: NULL");
}

void Maxit::_ndb_to_pdb_get_remark_500_Ramachandran_outliers(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_validate_torsion");
       if (!t) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_500_TORSION, Remark_500_TORSION, 0, 0, 0);

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();
       data.clear();
       data.push_back("auth_asym_id");
       data.push_back("auth_comp_id");
       data.push_back("auth_seq_id");
       data.push_back("PDB_ins_code");
       items.push_back(data);

       std::string cs;
       int count = 0, mol_index = 0;
       std::vector<RCSB::Residue*> residues;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            if (count >= 50) {
                 std::string remark = "THIS ENTRY HAS " + FloatToString((double) rowNo, 7, 0, false, true) + " RAMACHANDRAN OUTLIERS.";
                 _addNewRemark(500, "");
                 _addNewRemark(500, remark);
                 break;
            }

            get_value(cs, t, i, "PDB_model_num");
            int mol_id = atoi(cs.c_str());
            get_values(t, i, items, names);
            _getResiduesfromNames(mol_id, names, mol_index, residues, data); 
            if (residues.empty()) continue;
            count++;

            get_value(cs, t, i, "psi");
            double psi = atof(cs.c_str());
            get_value(cs, t, i, "phi");
            double phi = atof(cs.c_str());
            if (_molecules.size() < 2) mol_id = 0;
            _remark_500_Ramachandran_outliers(residues[0]->GetFirstAtom(), psi, phi, mol_id);
       }
       _addNewRemark(500, "");
       _addNewRemark(500, "REMARK: NULL");
}

void Maxit::_ndb_to_pdb_get_remark_500_non_cis_trans_torsions(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_validate_peptide_omega");
       if (!t) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_500_NON_CIS, Remark_500_NON_CIS, 0, 0, 0);

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();
       data.clear();
       data.push_back("auth_asym_id_1");
       data.push_back("auth_comp_id_1");
       data.push_back("auth_seq_id_1");
       data.push_back("PDB_ins_code_1");
       items.push_back(data);
       data.clear();
       data.push_back("auth_asym_id_2");
       data.push_back("auth_comp_id_2");
       data.push_back("auth_seq_id_2");
       data.push_back("PDB_ins_code_2");
       items.push_back(data);

       std::string cs;
       int count = 0, mol_index = 0;
       std::vector<RCSB::Residue*> residues;
       std::vector<RCSB::Atom*> atoms;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            if (count >= 50) {
                 std::string remark = "THIS ENTRY HAS " + FloatToString((double) rowNo, 7, 0, false, true) + " NON CIS, NON-TRANS OMEGA OUTLIERS.";
                 _addNewRemark(500, "");
                 _addNewRemark(500, remark);
                 break;
            }

            get_value(cs, t, i, "PDB_model_num");
            int mol_id = atoi(cs.c_str());
            get_values(t, i, items, names);
            _getResiduesfromNames(mol_id, names, mol_index, residues, data); 
            if (residues.size() < 2) continue;
            count++;

            atoms.clear();
            atoms.push_back(residues[0]->GetFirstAtom());
            atoms.push_back(residues[1]->GetFirstAtom());

            get_value(cs, t, i, "omega");
            double omega = atof(cs.c_str());
            if (_molecules.size() < 2) mol_id = 0;
            _remark_500_non_cis_trans_torsions(atoms, omega, mol_id);
       }
       _addNewRemark(500, "");
       _addNewRemark(500, "REMARK: NULL");
}

void Maxit::_ndb_to_pdb_get_remark_500_side_chain_plane(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_validate_planes");
       if (!t) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_500_SPLANE, Remark_500_SPLANE, 0, 0, 0);

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();
       data.clear();
       data.push_back("auth_asym_id");
       data.push_back("auth_comp_id");
       data.push_back("auth_seq_id");
       data.push_back("PDB_ins_code");
       items.push_back(data);

       std::string cs;
       int count = 0, mol_index = 0;
       std::vector<RCSB::Residue*> residues;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            if (count >= 50) {
                 std::string remark = "THIS ENTRY HAS " + FloatToString((double) rowNo, 7, 0, false, true) + " PLANE DEVIATIONS.";
                 _addNewRemark(500, "");
                 _addNewRemark(500, remark);
                 break;
            }

            get_value(cs, t, i, "PDB_model_num");
            int mol_id = atoi(cs.c_str());
            get_values(t, i, items, names);
            _getResiduesfromNames(mol_id, names, mol_index, residues, data); 
            if (residues.empty()) continue;
            count++;

            get_value(cs, t, i, "rmsd");
            double rmsd = atof(cs.c_str());
            get_value(cs, t, i, "type");
            if (_molecules.size() < 2) mol_id = 0;
            _remark_500_side_chain_plane(residues[0]->GetFirstAtom(), rmsd, mol_id, cs);
       }
       _addNewRemark(500, "");
       _addNewRemark(500, "REMARK: NULL");
}

void Maxit::_ndb_to_pdb_get_remark_500_main_chain_planarity(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_validate_main_chain_plane");
       if (!t) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_500_MPLANE, Remark_500_MPLANE, 0, 0, 0);

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();
       data.clear();
       data.push_back("auth_asym_id");
       data.push_back("auth_comp_id");
       data.push_back("auth_seq_id");
       data.push_back("PDB_ins_code");
       items.push_back(data);

       std::string cs;
       int count = 0, mol_index = 0;
       std::vector<RCSB::Residue*> residues;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            if (count >= 50) {
                 std::string remark = "THIS ENTRY HAS " + FloatToString((double) rowNo, 7, 0, false, true) + " MAIN CHAIN PLANARITY DEVIATIONS.";
                 _addNewRemark(500, "");
                 _addNewRemark(500, remark);
                 break;
            }

            get_value(cs, t, i, "PDB_model_num");
            int mol_id = atoi(cs.c_str());
            get_values(t, i, items, names);
            _getResiduesfromNames(mol_id, names, mol_index, residues, data); 
            if (residues.empty()) continue;
            count++;

            get_value(cs, t, i, "improper_torsion_angle");
            double pseudo_torsion = atof(cs.c_str());
            if (_molecules.size() < 2) mol_id = 0;
            _remark_500_main_chain_planarity(residues[0]->GetFirstAtom(), pseudo_torsion, mol_id);
       }
       _addNewRemark(500, "");
       _addNewRemark(500, "REMARK: NULL");
}

void Maxit::_ndb_to_pdb_get_remark_525()
{
       if (_molecules.empty()) return;
       if (!_CifObj) return;
       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);
       ISTable *t = getTablePtr(_cifblock, "pdbx_distant_solvent_atoms");
       if (!t) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_525, Remark_525, 0, 0, 0);

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();
       data.clear();
       data.push_back("auth_asym_id");
       data.push_back("auth_comp_id");
       data.push_back("auth_seq_id");
       data.push_back("PDB_ins_code");
       data.push_back("auth_atom_id");
       data.push_back("label_alt_id");
       items.push_back(data);

       std::string cs;
       int rowNo = t->GetNumRows();
       int mol_index = 0;
       std::vector<RCSB::Atom*> atoms;
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs, t, i, "PDB_model_num");
            int mol_id = atoi(cs.c_str());

            get_values(t, i, items, names);
            _getAtomsfromNames(mol_id, names, mol_index, atoms, data); 
            if (atoms.empty()) continue;

            get_value(cs, t, i, "neighbor_macromolecule_distance");
            if (cs.empty()) get_value(cs, t, i, "neighbor_ligand_distance");
            double dist = atof(cs.c_str());

            _remark_525(atoms[0], dist);
       }
}
