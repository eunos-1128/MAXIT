/*
FILE:     ContactUtil.C
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
#include <math.h>

#include "BondUtil.h"
#include "GetPairList.h"
#include "GetPairList.C"
#include "SeqCodeUtil.h"
#include "ValidateObj.h"
#include "utillib.h"
#include "xtal.h"

void ValidateObj::calculate_all_contact()
{
       if (_molecules.empty()) return;

       bool /* exclude_hydrogen = true;
       if (_experiment_type & EXPERIMENT_TYPE_NMR ||
           _experiment_type & EXPERIMENT_TYPE_NMR_SOLID) */
            exclude_hydrogen = false;

       _clear_contact();
       _create_link_mapping();

       _cal_a_contact_flag = true;
       _cal_s_contact_flag = true;
       _cal_special_position_flag = true;

       Value value;
       std::set<std::string> chain_types;
       chain_types.clear();
       std::set<std::string> residue_set;
       residue_set.clear();
       std::vector<CONTACT> contact;
       std::set<std::string> UniIndex, SpcIndex;
       std::multimap<double, Value> _Index;
       std::string cs;

       std::map<unsigned int, std::string> sym_occ_map;
       sym_occ_map.clear();
       sym_occ_map.insert(std::make_pair(2, "0.50"));
       sym_occ_map.insert(std::make_pair(3, "0.33"));
       sym_occ_map.insert(std::make_pair(4, "0.25"));
       sym_occ_map.insert(std::make_pair(6, "0.16"));

       // first key: residue ID - pdbchnid_resname_resnum_inscode
       // second key: atom ID - atmnam_altloc
       // values: symmetry operators
       std::map<std::string, std::map<std::string, std::set<std::string> > > special_residues;

       for (unsigned int i = 0; i < _molecules.size(); ++i) {
            xtal mycrys;
            if (!set_mycrys(mycrys, _cell, _molecules[i], chain_types, residue_set, exclude_hydrogen, false, false)) return;

            get_contact(mycrys, contact, 0.0, 2.2, "asym", "all_asym", "all_asym", 1, 1);
            UniIndex.clear();
            _Index.clear();
            for (unsigned int j = 0; j < contact.size(); ++j) {
                 if (!_cell.is_artifical()) {
                      if (atoi(contact[j].a_atm->atnum().c_str()) < atoi(contact[j].b_atm->atnum().c_str()))
                           cs = contact[j].a_atm->atnum() + "_" + contact[j].b_atm->atnum();
                      else cs = contact[j].b_atm->atnum() + "_" + contact[j].a_atm->atnum();
                      UniIndex.insert(cs);
                 }

                 // filter low occupancy contacts
                 if ((atof(contact[j].a_atm->occ().c_str()) < 1.0 && atof(contact[j].b_atm->occ().c_str()) < 1.0) ||
                     (atof(contact[j].a_atm->occ().c_str()) < 1.0 && (contact[j].a_atm->restype() == "HOH" || contact[j].a_atm->restype() == "DOD")) ||
                     (atof(contact[j].b_atm->occ().c_str()) < 1.0 && (contact[j].b_atm->restype() == "HOH" || contact[j].b_atm->restype() == "DOD")))
                      continue;

                 if (_is_a_link(contact[j].a_atm, contact[j].b_atm, contact[j].dist)) continue;
     
                 value.clear();
                 value.set_Mol_ID(_molecules[i]->Mol_ID());
                 value.set_val(contact[j].dist);
                 if (atoi(contact[j].a_atm->atnum().c_str()) < atoi(contact[j].b_atm->atnum().c_str())) {
                      value.insert_atom(contact[j].a_atm);
                      value.insert_atom(contact[j].b_atm);
                 } else {
                      value.insert_atom(contact[j].b_atm);
                      value.insert_atom(contact[j].a_atm);
                 }
                 _Index.insert(std::make_pair(contact[j].dist, value));
            }
            if (!_Index.empty()) {
                 for (std::multimap<double, Value>::const_iterator pos = _Index.begin(); pos != _Index.end(); ++pos) {
                      _a_contact.push_back(pos->second);
                 }
            }

            if (_cell.is_artifical()) continue;

            // get_contact(mycrys, contact, 0.0, 2.2, "sym", "all_asym", "all", 1, 1);
            get_sym_contact_with_special_position_residues(mycrys, contact, special_residues, 0.0, 2.2, "all_asym", "all");

            SpcIndex.clear();
            _Index.clear();
            for (unsigned j = 0; j < contact.size(); ++j) {
                 if (_is_a_link(contact[j].a_atm, contact[j].b_atm, contact[j].dist)) continue;

                 std::string sym = String::IntToString(contact[j].sym) + "_" + String::IntToString(contact[j].lx + 5)
                                 + String::IntToString(contact[j].ly + 5) + String::IntToString(contact[j].lz + 5);

                 if (contact[j].type < 0) {
                      if (contact[j].dist < SPECIAL_WATER_DIST_CUTOFF && contact[j].a_atm == contact[j].b_atm) {
                           cs = contact[j].a_atm->restype() + "_" + contact[j].a_atm->pdb_chnid() + "_" + contact[j].a_atm->pdb_resnum() + "_"
                              + contact[j].a_atm->ins_code() + "_" + contact[j].a_atm->atmtype();

                           if (SpcIndex.find(cs) == SpcIndex.end()) {
                                SpcIndex.insert(cs);
                                value.clear();
                                value.set_Mol_ID(_molecules[i]->Mol_ID());
                                value.set_val(contact[j].dist);
                                value.set_sval(sym);
                                value.insert_atom(contact[j].a_atm);
                                _special_position_atoms.push_back(value);
                           }
                      }
                 } else {
                      // filter low occupancy contacts
                      if ((atof(contact[j].a_atm->occ().c_str()) < 1.0 && atof(contact[j].b_atm->occ().c_str()) < 1.0) ||
                          (atof(contact[j].a_atm->occ().c_str()) < 1.0 && (contact[j].a_atm->restype() == "HOH" ||
                          contact[j].a_atm->restype() == "DOD")) || (atof(contact[j].b_atm->occ().c_str()) < 1.0 &&
                         (contact[j].b_atm->restype() == "HOH" || contact[j].b_atm->restype() == "DOD"))) continue;

                      if (atoi(contact[j].a_atm->atnum().c_str()) < atoi(contact[j].b_atm->atnum().c_str()))
                           cs = contact[j].a_atm->atnum() + "_" + contact[j].b_atm->atnum();
                      else cs = contact[j].b_atm->atnum() + "_" + contact[j].a_atm->atnum();
                      if (UniIndex.find(cs) == UniIndex.end()) {
                           UniIndex.insert(cs);
                           value.clear();
                           value.set_Mol_ID(_molecules[i]->Mol_ID());
                           value.set_val(contact[j].dist);
                           value.set_sval(sym);
                           value.insert_atom(contact[j].a_atm);
                           value.insert_atom(contact[j].b_atm);
                           _Index.insert(std::make_pair(contact[j].dist, value));
                      }
                 }
            }
            if (!_Index.empty()) {
                 for (std::multimap<double, Value>::const_iterator pos = _Index.begin(); pos != _Index.end(); ++pos) {
                      _s_contact.push_back(pos->second);
                 }
            }

            if (special_residues.empty()) continue;

            for (std::map<std::string, std::map<std::string, std::set<std::string> > >::const_iterator
                 mpos = special_residues.begin(); mpos != special_residues.end(); ++mpos) {
                 RCSB::Residue* res = _molecules[i]->find_pdb_residue(mpos->first);
                 if (!res) continue;

                 bool update_flag = false;
                 if (res->ResName() == "HOH" || res->ResName() == "DOD") update_flag = true;
                 else if (res->NonHydrogenAtomNumbers() == 1) {
                      try {
                           const ConnectFormat& drug = _ccDic->find_drug(res->ResName());
                           if (drug.nheavyatoms() == 1) update_flag = true;
                      } catch (const std::exception& exc) {}
                 }
                 if (!update_flag) continue;

                 int sys_op_number = mpos->second.begin()->second.size() + 1;
                 std::map<unsigned int, std::string>::const_iterator sopos = sym_occ_map.find(sys_op_number);
                 if (sopos != sym_occ_map.end()) res->update_occupancy(sopos->second);
            }
       }
}

void ValidateObj::CalculateCovalentBonds(const std::string& link_radii)
{
       if (_molecules.empty()) return;

       ISTable *t = new ISTable("pdbx_struct_link");
       t->AddColumn("id");
       t->AddColumn("pdbx_model");
       t->AddColumn("ptnr1_label_alt_id");
       t->AddColumn("ptnr1_label_asym_id");
       t->AddColumn("ptnr1_label_atom_id");
       t->AddColumn("ptnr1_label_comp_id");
       t->AddColumn("ptnr1_label_seq_id");
       t->AddColumn("ptnr1_label_ins_code");
       t->AddColumn("ptnr1_x");
       t->AddColumn("ptnr1_y");
       t->AddColumn("ptnr1_z");
       t->AddColumn("ptnr2_label_alt_id");
       t->AddColumn("ptnr2_label_asym_id");
       t->AddColumn("ptnr2_label_atom_id");
       t->AddColumn("ptnr2_label_comp_id");
       t->AddColumn("ptnr2_label_seq_id");
       t->AddColumn("ptnr2_label_ins_code");
       t->AddColumn("ptnr2_x");
       t->AddColumn("ptnr2_y");
       t->AddColumn("ptnr2_z");
       t->AddColumn("pdbx_dist_value");

       double cutoff = 3.1;
       if (!link_radii.empty()) cutoff = 3.1 + atof(link_radii.c_str());

       std::set<std::string> chain_types, residue_set, annotated_linkage_set;
       chain_types.clear();
       residue_set.clear();
       annotated_linkage_set.clear();

       for (std::list<_LINK>::const_iterator lpos = _links.begin(); lpos != _links.end(); ++lpos) {
            if ((!lpos->SymOP_1.empty() && (lpos->SymOP_1 != "1_555")) || (!lpos->SymOP_2.empty() && (lpos->SymOP_2 != "1_555"))) continue; 
            std::string idx = String::IntToString(lpos->mol_index) + "_" + lpos->fstAtom->getAtomAllIndex() + "_" + lpos->sndAtom->getAtomAllIndex();
            annotated_linkage_set.insert(idx);
            idx = String::IntToString(lpos->mol_index) + "_" + lpos->sndAtom->getAtomAllIndex() + "_" + lpos->fstAtom->getAtomAllIndex();
            annotated_linkage_set.insert(idx);
       }

       std::map<std::string, std::vector<std::vector<std::string> > > covalent_bonding;
       covalent_bonding.clear();

       std::vector<CONTACT> contact;
       for (unsigned int i = 0; i < _molecules.size(); ++i) {
            xtal mycrys;
            if (!set_mycrys(mycrys, _cell, _molecules[i], chain_types, residue_set, false, false, false)) continue;
            get_contact(mycrys, contact, 0.80, cutoff, "asym", "all_asym", "all_asym", 1, 0, false);
            if (contact.empty()) continue;

            for (std::vector<CONTACT>::const_iterator cpos = contact.begin(); cpos != contact.end(); ++cpos) {
                 if (cpos->sym != 1 || cpos->lx || cpos->ly || cpos->lz) continue;

                 if (cpos->a_atm->pdb_resnam() == "HOH" || cpos->a_atm->pdb_resnam() == "DOD" ||
                     cpos->b_atm->pdb_resnam() == "HOH" || cpos->b_atm->pdb_resnam() == "DOD") continue;

                 if ((SeqCodeUtil::is_a_standard_residue(cpos->a_atm->pdb_resnam()) || cpos->a_atm->pdb_resnam() == "MSE") &&
                     (SeqCodeUtil::is_a_standard_residue(cpos->b_atm->pdb_resnam()) || cpos->b_atm->pdb_resnam() == "MSE")) continue;

                 if (BondUtil::is_a_link_flexible(cpos->a_atm->atom_type().c_str(), cpos->b_atm->atom_type().c_str(), cpos->dist) != 1) continue;

                 int row = t->GetNumRows();
                 t->AddRow();
                 t->UpdateCell(row, "id", String::IntToString(row + 1));
                 t->UpdateCell(row, "pdbx_model", String::IntToString(_molecules[i]->Mol_ID()));
                 t->UpdateCell(row, "ptnr1_label_alt_id",   cpos->a_atm->alt_loc());
                 t->UpdateCell(row, "ptnr1_label_asym_id",  cpos->a_atm->pdb_chnid());
                 t->UpdateCell(row, "ptnr1_label_atom_id",  cpos->a_atm->pdb_atmnam());
                 t->UpdateCell(row, "ptnr1_label_comp_id",  cpos->a_atm->pdb_resnam());
                 t->UpdateCell(row, "ptnr1_label_seq_id",   cpos->a_atm->pdb_resnum());
                 t->UpdateCell(row, "ptnr1_label_ins_code", cpos->a_atm->ins_code());
                 t->UpdateCell(row, "ptnr1_x", FloatToString(cpos->a_atm->orig().x, 0, 3));
                 t->UpdateCell(row, "ptnr1_y", FloatToString(cpos->a_atm->orig().y, 0, 3));
                 t->UpdateCell(row, "ptnr1_z", FloatToString(cpos->a_atm->orig().z, 0, 3));
                 t->UpdateCell(row, "ptnr2_label_alt_id",   cpos->b_atm->alt_loc());
                 t->UpdateCell(row, "ptnr2_label_asym_id",  cpos->b_atm->pdb_chnid());
                 t->UpdateCell(row, "ptnr2_label_atom_id",  cpos->b_atm->pdb_atmnam());
                 t->UpdateCell(row, "ptnr2_label_comp_id",  cpos->b_atm->pdb_resnam());
                 t->UpdateCell(row, "ptnr2_label_seq_id",   cpos->b_atm->pdb_resnum());
                 t->UpdateCell(row, "ptnr2_label_ins_code", cpos->b_atm->ins_code());
                 t->UpdateCell(row, "ptnr2_x", FloatToString(cpos->b_atm->orig().x, 0, 3));
                 t->UpdateCell(row, "ptnr2_y", FloatToString(cpos->b_atm->orig().y, 0, 3));
                 t->UpdateCell(row, "ptnr2_z", FloatToString(cpos->b_atm->orig().z, 0, 3));
                 t->UpdateCell(row, "pdbx_dist_value", FloatToString(cpos->dist, 0, 3));

                 std::string idx = String::IntToString(i) + "_" + cpos->a_atm->getAtomAllIndex() + "_" + cpos->b_atm->getAtomAllIndex();
                 if ((!annotated_linkage_set.empty() && (annotated_linkage_set.find(String::IntToString(i) + "_" + cpos->a_atm->getAtomAllIndex() + "_" +
                       cpos->b_atm->getAtomAllIndex()) == annotated_linkage_set.end())) || (annotated_linkage_set.empty() &&
                      !BondUtil::is_a_link(cpos->a_atm->atom_type().c_str(), cpos->b_atm->atom_type().c_str(), cpos->dist))) continue;

                 std::string a_resIdx = String::IntToString(_molecules[i]->Mol_ID()) + "_" + cpos->a_atm->pdb_chnid() + "_" + cpos->a_atm->pdb_resnam()
                                      + "_" + cpos->a_atm->pdb_resnum() + "_" + cpos->a_atm->ins_code();
                 std::string a_atomIdx = a_resIdx + "_" + cpos->a_atm->pdb_atmnam();
                 if (!cpos->a_atm->alt_loc().empty()) a_atomIdx += "(" + cpos->a_atm->alt_loc() + ")";
                 std::string b_resIdx = String::IntToString(_molecules[i]->Mol_ID()) + "_" + cpos->b_atm->pdb_chnid() + "_" + cpos->b_atm->pdb_resnam()
                                      + "_" + cpos->b_atm->pdb_resnum() + "_" + cpos->b_atm->ins_code();
                 std::string b_atomIdx = b_resIdx + "_" + cpos->b_atm->pdb_atmnam();
                 if (!cpos->b_atm->alt_loc().empty()) b_atomIdx += "(" + cpos->b_atm->alt_loc() + ")";
                 _insert_covalent_bonding_info(a_resIdx, a_atomIdx, b_atomIdx, FloatToString(cpos->dist, 0, 3), covalent_bonding);
                 _insert_covalent_bonding_info(b_resIdx, b_atomIdx, a_atomIdx, FloatToString(cpos->dist, 0, 3), covalent_bonding);
            }
       }

       if (t->GetNumRows() == 0) {
            delete t;
            return;
       }

       CifFile *fobj = create_fobj("", _StructureId);
       Block &block = fobj->GetBlock(_StructureId);
       block.WriteTable(t);
       if (!covalent_bonding.empty()) {
            t = new ISTable("pdbx_covalent_bonding");
            t->AddColumn("inst_id");
            t->AddColumn("first_atom_id");
            t->AddColumn("second_atom_id");
            t->AddColumn("dist");
            for (std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator
                 mpos = covalent_bonding.begin(); mpos != covalent_bonding.end(); ++mpos) {
                 for (std::vector<std::vector<std::string> >::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                      int row = t->GetNumRows();
                      t->AddRow();
                      t->UpdateCell(row, "inst_id", mpos->first);
                      t->UpdateCell(row, "first_atom_id", (*vpos)[0]);
                      t->UpdateCell(row, "second_atom_id", (*vpos)[1]);
                      t->UpdateCell(row, "dist", (*vpos)[2]);
                 }
            }
            block.WriteTable(t);
       }
       fobj->Write(_output_filename);
       delete fobj;
}

std::list<Value> ValidateObj::_get_special_position_atoms(const bool& pdb_format_flag)
{
       double cut_off = SPECIAL_WATER_DIST_CUTOFF;
       if (pdb_format_flag) cut_off = 0.3;

       std::list<Value> sp_atom_list;
       sp_atom_list.clear();

       if (_molecules.empty()) return sp_atom_list;
       if (!_CifObj) return sp_atom_list;
       Block &block = _CifObj->GetBlock(_firstBlockName);
       ISTable *t = getTablePtr(block, "pdbx_struct_special_symmetry");
       if (!t) return sp_atom_list;

       std::vector<RCSB::Residue*> residues;
       residues.clear();

       std::string chnid, resnam, resnum, ins_code, cs;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(chnid, t, i, "auth_asym_id");
            get_value_clean(resnam, t, i, "auth_comp_id");
            get_value_clean(resnum, t, i, "auth_seq_id");
            get_value_clean(ins_code, t, i, "PDB_ins_code");
            RCSB::Residue* res = _molecules[0]->find_pdb_residue(chnid, resnam, resnum, ins_code);
            if (res) residues.push_back(res); 
       }
       if (residues.empty()) return sp_atom_list;

       xtal mycrys;
       set_mycrys(mycrys, _cell, residues, true, false);

       std::set<std::string> SpcIndex;
       SpcIndex.clear();

       std::vector<CONTACT> contact;
       get_contact(mycrys, contact, 0.0, 2.2, "sym", "all_asym", "all", 1, 1);

       Value value;
       for (std::vector<CONTACT>::const_iterator pos = contact.begin(); pos != contact.end(); ++pos) {
            if ((pos->type >= 0) || (pos->dist >= cut_off) || (pos->a_atm != pos->b_atm)) continue;

            cs = pos->a_atm->restype() + "_" + pos->a_atm->pdb_chnid() + "_" + pos->a_atm->pdb_resnum()
               + "_" + pos->a_atm->ins_code() + "_" + pos->a_atm->atmtype();
            if (SpcIndex.find(cs) != SpcIndex.end()) continue;
            SpcIndex.insert(cs);

            cs = String::IntToString(pos->sym) + "_" + String::IntToString(pos->lx + 5) + String::IntToString(pos->ly + 5) + String::IntToString(pos->lz + 5);
            value.clear();
            value.set_Mol_ID(_molecules[0]->Mol_ID());
            value.set_val(pos->dist);
            value.set_sval(cs);
            value.insert_atom(pos->a_atm);
            sp_atom_list.push_back(value);
       }

       return sp_atom_list;
}

double ValidateObj::_cal_distance_between_residues(RCSB::Residue* res1, const std::string& F_atom, RCSB::Residue* res2, const std::string& S_atom)
{
       double dist = 100000.0;

       std::vector<std::pair<RCSB::Residue*, std::string> > res_atom_pair_list;
       res_atom_pair_list.clear();
       res_atom_pair_list.push_back(std::make_pair(res1, F_atom));
       res_atom_pair_list.push_back(std::make_pair(res2, S_atom));

       std::vector<std::vector<RCSB::Atom*> > a_pair_list;
       _get_atom_pair_list(res_atom_pair_list, a_pair_list);
       if (a_pair_list.empty()) return dist;

       for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = a_pair_list.begin(); apos != a_pair_list.end(); ++apos) {
            double dist1 = cal_distance((*apos)[0], (*apos)[1]);
            if (dist1 < dist) dist = dist1;
       }

       return dist;
}

void ValidateObj::_get_atom_pair_list(const std::vector<std::pair<RCSB::Residue*, std::string> >& residue_atom_name_pair,
                                      std::vector<std::vector<RCSB::Atom*> >& pair_lists)
{
       pair_lists.clear();
       if (residue_atom_name_pair.empty()) return;

       std::vector<RCSB::Atom*> atoms;
       std::vector<std::vector<RCSB::Atom*> > atom_lists;
       atom_lists.clear();
       for (std::vector<std::pair<RCSB::Residue*, std::string> >::const_iterator pos = residue_atom_name_pair.begin();
            pos != residue_atom_name_pair.end(); ++pos) {
            pos->first->find_atom(pos->second, atoms);
            if (atoms.empty()) return;
            atom_lists.push_back(atoms);
       }
       GetPairList(atom_lists, pair_lists);
}

void ValidateObj::_insert_covalent_bonding_info(const std::string& resIdx, const std::string& fstAtomIdx, const std::string& sndAtomIdx, const std::string&
                                                dist, std::map<std::string, std::vector<std::vector<std::string> > >& covalent_bonding)
{
       std::vector<std::string> data;
       data.clear();
       data.push_back(fstAtomIdx);
       data.push_back(sndAtomIdx);
       data.push_back(dist);

       std::map<std::string, std::vector<std::vector<std::string> > >::iterator mpos = covalent_bonding.find(resIdx);
       if (mpos != covalent_bonding.end()) mpos->second.push_back(data);
       else {
            std::vector<std::vector<std::string> > t_vec;
            t_vec.clear();
            t_vec.push_back(data);
            covalent_bonding.insert(std::make_pair(resIdx, t_vec));
       }
}
