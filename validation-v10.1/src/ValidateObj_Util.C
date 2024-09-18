/*
FILE:     ValidateObj_Util.C
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

#include "GetPairList.h"
#include "GetPairList.C"
#include "NdbToken.h"
#include "PdbWrite.h"
#include "SeqCodeUtil.h"
#include "TypeDef.h"
#include "utillib.h"
#include "ValidateObj.h"

static void insert_bonded_atoms(std::map<std::string, std::set<std::string> >& bonded_atoms, const int& mol_index,
                                const int& res_index, const std::string& atom_name);

void ValidateObj::_get_missing_or_zero_occupancy_residues_or_atoms()
{
       _unobs_or_zero_occ_residues.clear();
       _unobs_or_zero_occ_atoms.clear();
       _extra_atoms.clear();

       if (_molecules.empty()) return;

       _finding_missing_and_extra_atoms();

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > atom_list;
       std::vector<RCSB::Residue*> residues;
       for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            RCSB::Chain* chain = (*mpos)->GetFirstChain();
            while (chain) {
                 std::string polymer_flag = "N";
                 if (chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP")
                      polymer_flag = "Y";

                 for (unsigned int j = 0; j < chain->SeqLen(); ++j) {
                      _FIELD* SeqRes = chain->SeqRes(j); 
                      if (SeqRes->ResIndex < 0) {
                           data.clear();
                           data.push_back(String::IntToString((*mpos)->Mol_ID()));
                           data.push_back(polymer_flag);
                           data.push_back("1");
                           data.push_back(chain->PDB_ChainID());
                           data.push_back(SeqRes->Field[0]);
                           data.push_back(SeqRes->Field[4]);
                           data.push_back(SeqRes->InsCode);
                           data.push_back(chain->ChainID());
                           data.push_back(SeqRes->Field[0]);
                           data.push_back(SeqRes->Field[2]);
                           _unobs_or_zero_occ_residues.push_back(data);
                           continue;
                      }
                      if (chain->ca_or_p_atom_only()) continue;

                      chain->GetResidueListByIndex(SeqRes->ResIndex, residues);
                      for (std::vector<RCSB::Residue*>::iterator
                           rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                           const std::vector<std::string>& missing = (*rpos)->missing();
                           if (!missing.empty()) {
                                for (std::vector<std::string>::const_iterator
                                     pos = missing.begin(); pos != missing.end(); ++pos) {
                                     data.clear();
                                     data.push_back(String::IntToString((*mpos)->Mol_ID()));
                                     data.push_back(polymer_flag);
                                     data.push_back("1");
                                     data.push_back(chain->PDB_ChainID());
                                     data.push_back(SeqRes->Field[0]);
                                     data.push_back(SeqRes->Field[4]);
                                     data.push_back(SeqRes->InsCode);
                                     data.push_back(*pos);
                                     data.push_back("");
                                     data.push_back(chain->ChainID());
                                     data.push_back(SeqRes->Field[0]);
                                     data.push_back(SeqRes->Field[2]);
                                     data.push_back(*pos);
                                     _unobs_or_zero_occ_atoms.push_back(data);
                                }
                           }
                           const std::vector<std::string>& extras = (*rpos)->extras();
                           if (!extras.empty()) {
                                data.clear();
                                data.push_back(String::IntToString((*mpos)->Mol_ID()));
                                data.push_back(chain->PDB_ChainID());
                                data.push_back(SeqRes->Field[0]);
                                data.push_back(SeqRes->Field[4]);
                                data.push_back(SeqRes->InsCode);
                                _extra_atoms.push_back(std::make_pair(data, extras));
                           }

                           if (!(_experiment_type & EXPERIMENT_TYPE_XRAY)) continue;

                           bool ret = (*rpos)->is_zero_occupancy_residue(atom_list);
                           if (ret) {
                                data.clear();
                                data.push_back(String::IntToString((*mpos)->Mol_ID()));
                                data.push_back(polymer_flag);
                                data.push_back("0");
                                data.push_back((*rpos)->pdb_chnid());
                                data.push_back((*rpos)->ResName());
                                data.push_back((*rpos)->pdb_res_no());
                                data.push_back((*rpos)->ins_code());
                                data.push_back((*rpos)->chnid());
                                data.push_back((*rpos)->ResName());
                                data.push_back((*rpos)->res_no());
                                _unobs_or_zero_occ_residues.push_back(data);
                           } else if (!atom_list.empty()) {
                                data.clear();
                                data.push_back(String::IntToString((*mpos)->Mol_ID()));
                                data.push_back(polymer_flag);
                                data.push_back("0");
                                data.push_back((*rpos)->pdb_chnid());
                                data.push_back((*rpos)->ResName());
                                data.push_back((*rpos)->pdb_res_no());
                                data.push_back((*rpos)->ins_code());
                                data.push_back("");
                                data.push_back("");
                                data.push_back((*rpos)->chnid());
                                data.push_back((*rpos)->ResName());
                                data.push_back((*rpos)->res_no());
                                data.push_back("");
                                for (std::vector<std::vector<std::string> >::const_iterator
                                     pos = atom_list.begin(); pos != atom_list.end(); ++pos) {
                                     data[7] = (*pos)[0];
                                     data[8] = (*pos)[1];
                                     data[12] = (*pos)[0];
                                     _unobs_or_zero_occ_atoms.push_back(data);
                                }
                           }
                      }
                 }
                 chain = (*mpos)->GetNextChain();
            }
       }
}

void ValidateObj::_finding_missing_and_extra_atoms(const bool& exclude_hydrogen)
{
       if (_molecules.empty()) return;
/*
       // for ASP, DAS, DGL & GLU
       std::set<std::string> non_standard_linkage;
       non_standard_linkage.clear();
*/
       // key: molecule index + "_" + residue index
       // value: bonded atom set
       std::map<std::string, std::set<std::string> > bonded_atoms;
       bonded_atoms.clear();

       // key: residue index
       std::set<int> c_n_bonds;

       for (std::list<_LINK>::const_iterator lpos = _links.begin(); lpos != _links.end(); ++lpos) {
/*
            if (lpos->type != "covale") continue;
            if (((lpos->fstAtom->pdb_resnam() == "ASP" || lpos->fstAtom->pdb_resnam() == "DAS") && lpos->fstAtom->pdb_atmnam() == "CG") ||
                ((lpos->fstAtom->pdb_resnam() == "DGL" || lpos->fstAtom->pdb_resnam() == "GLU") && lpos->fstAtom->pdb_atmnam() == "CD"))
                 non_standard_linkage.insert(lpos->fstAtom->pdb_chnid() + "_" + lpos->fstAtom->pdb_resnam() + "_" +
                                             lpos->fstAtom->pdb_resnum() + "_" + lpos->fstAtom->ins_code());
            if (((lpos->sndAtom->pdb_resnam() == "ASP" || lpos->sndAtom->pdb_resnam() == "DAS") && lpos->sndAtom->pdb_atmnam() == "CG") ||
                ((lpos->sndAtom->pdb_resnam() == "DGL" || lpos->sndAtom->pdb_resnam() == "GLU") && lpos->sndAtom->pdb_atmnam() == "CD"))
                 non_standard_linkage.insert(lpos->sndAtom->pdb_chnid() + "_" + lpos->sndAtom->pdb_resnam() + "_" +
                                             lpos->sndAtom->pdb_resnum() + "_" + lpos->sndAtom->ins_code());
*/
            for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
                 RCSB::Residue* res = (*mpos)->find_pdb_residue(lpos->fstAtom->pdb_chnid(), lpos->fstAtom->pdb_resnam(), lpos->fstAtom->pdb_resnum(),
                                                                lpos->fstAtom->ins_code());
                 if (res) insert_bonded_atoms(bonded_atoms, (*mpos)->index(), res->index(), lpos->fstAtom->pdb_atmnam());
                 res = (*mpos)->find_pdb_residue(lpos->sndAtom->pdb_chnid(), lpos->sndAtom->pdb_resnam(), lpos->sndAtom->pdb_resnum(),
                                                       lpos->sndAtom->ins_code());
                 if (res) insert_bonded_atoms(bonded_atoms, (*mpos)->index(), res->index(), lpos->sndAtom->pdb_atmnam());
            }
       }

       int need_hydrogen = 0;
       if (_experiment_type & EXPERIMENT_TYPE_NMR ||
           _experiment_type & EXPERIMENT_TYPE_NMR_SOLID) need_hydrogen = 1;
       if (exclude_hydrogen) need_hydrogen = 0;

       std::set<std::string> EmptyAtoms, AA_N_Terminal_Allowed_Atoms;
       std::set<std::string> AA_N_Terminal_Leaving_Atoms, AA_N_Terminal_Not_Allowed_Atoms;
       std::set<std::string> AA_Leaving_Atoms, AA_Not_Allowed_Atoms, NA_5_Terminal_Allowed_Atoms;
       std::set<std::string> NA_5_Terminal_Leaving_Atoms, NA_5_Terminal_Not_Allowed_Atoms;
       std::set<std::string> NA_C_Terminal_Leaving_Atoms, NA_Leaving_Atoms, NA_Not_Allowed_Atoms;
       std::set<std::string> CYS_N_Terminal_Leaving_Atoms, CYS_Leaving_Atoms;
       std::set<std::string> HIS_N_Terminal_Leaving_Atoms, HIS_Leaving_Atoms;
       std::set<std::string> CYS_C_Terminal_Leaving_Atoms, HIS_C_Terminal_Leaving_Atoms;

       EmptyAtoms.clear();
       AA_N_Terminal_Allowed_Atoms.clear();
       CYS_C_Terminal_Leaving_Atoms.clear();
       HIS_C_Terminal_Leaving_Atoms.clear();
       for (int i = 0; i < NUM_AA_N_TERMINAL_ATOMS; ++i) {
            AA_N_Terminal_Allowed_Atoms.insert(_AA_N_Terminal_Atoms[i]);
            CYS_C_Terminal_Leaving_Atoms.insert(_AA_N_Terminal_Atoms[i]);
            HIS_C_Terminal_Leaving_Atoms.insert(_AA_N_Terminal_Atoms[i]);
       }
       CYS_C_Terminal_Leaving_Atoms.insert("HG");
       HIS_C_Terminal_Leaving_Atoms.insert("HE2");
       HIS_C_Terminal_Leaving_Atoms.insert("HD1");

       AA_N_Terminal_Leaving_Atoms.clear();
       CYS_N_Terminal_Leaving_Atoms.clear();
       HIS_N_Terminal_Leaving_Atoms.clear();
       for (int i = 0; i < NUM_AA_N_LEAVING_ATOMS; ++i) {
            AA_N_Terminal_Leaving_Atoms.insert(_AA_N_Leaving_Atoms[i]);
            CYS_N_Terminal_Leaving_Atoms.insert(_AA_N_Leaving_Atoms[i]);
            HIS_N_Terminal_Leaving_Atoms.insert(_AA_N_Leaving_Atoms[i]);
       }
       CYS_N_Terminal_Leaving_Atoms.insert("HG");
       HIS_N_Terminal_Leaving_Atoms.insert("HE2");
       HIS_N_Terminal_Leaving_Atoms.insert("HD1");

       AA_N_Terminal_Not_Allowed_Atoms.clear();
       for (int i = 0; i < NUM_AA_C_TERMINAL_ATOMS; ++i) {
            AA_N_Terminal_Not_Allowed_Atoms.insert(_AA_C_Terminal_Atoms[i]);
       }
 
       AA_Leaving_Atoms.clear();
       CYS_Leaving_Atoms.clear();
       HIS_Leaving_Atoms.clear();
       for (int i = 0; i < NUM_AA_LEAVING_ATOMS; ++i) {
            AA_Leaving_Atoms.insert(_AA_Leaving_Atoms[i]);
            CYS_Leaving_Atoms.insert(_AA_Leaving_Atoms[i]);
            HIS_Leaving_Atoms.insert(_AA_Leaving_Atoms[i]);
       }
       CYS_Leaving_Atoms.insert("HG");
       HIS_Leaving_Atoms.insert("HE2");
       HIS_Leaving_Atoms.insert("HD1");

       AA_Not_Allowed_Atoms.clear();
       for (int i = 0; i < NUM_AA_TERMINAL_ATOMS; ++i) {
            AA_Not_Allowed_Atoms.insert(_AA_Terminal_Atoms[i]);
       }
 
       NA_5_Terminal_Allowed_Atoms.clear();
       for (int i = 0; i < NUM_NA_5_TERMINAL_ATOMS; ++i) {
            NA_5_Terminal_Allowed_Atoms.insert(_NA_5_Terminal_Atoms[i]);
       }

       NA_5_Terminal_Leaving_Atoms.clear();
       for (int i = 0; i < NUM_NA_5_LEAVING_ATOMS; ++i) {
            NA_5_Terminal_Leaving_Atoms.insert(_NA_5_Leaving_Atoms[i]);
       }

       NA_5_Terminal_Not_Allowed_Atoms.clear();
       for (int i = 0; i < NUM_NA_3_TERMINAL_ATOMS; ++i) {
            NA_5_Terminal_Not_Allowed_Atoms.insert(_NA_3_Terminal_Atoms[i]);
       }

       NA_C_Terminal_Leaving_Atoms.clear();
       for (int i = 0; i < NUM_NA_3_LEAVING_ATOMS; ++i) {
            NA_C_Terminal_Leaving_Atoms.insert(_NA_3_Leaving_Atoms[i]);
       }
       
       NA_Leaving_Atoms.clear();
       for (int i = 0; i < NUM_NA_LEAVING_ATOMS; ++i) {
            NA_Leaving_Atoms.insert(_NA_Leaving_Atoms[i]);
       }
 
       NA_Not_Allowed_Atoms.clear();
       for (int i = 0; i < NUM_NA_TERMINAL_ATOMS; ++i) {
            NA_Not_Allowed_Atoms.insert(_NA_Terminal_Atoms[i]);
       }
 
       std::vector<RCSB::Residue*> residues, next;
       std::vector<std::vector<RCSB::Residue*> > r_pair_list, residue_lists;
       for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            RCSB::Chain* chain = (*mpos)->GetFirstChain();
            while (chain) {
                 int is_connect = 0;
                 if ((chain->chain_type() == "ATOMP" ||
                      chain->chain_type() == "ATOMN") &&
                      chain->SeqLen() > 1) is_connect = 1;

                 for (unsigned int j = 0; j < chain->SeqLen(); ++j) {
                      chain->GetResidueList(j, residues);
                      if (residues.empty()) continue;

                      c_n_bonds.clear();
                      if ((chain->chain_type() == "ATOMP") && ((j + 1) < chain->SeqLen())) {
                           chain->GetResidueList(j + 1, next);
                           if (!next.empty()) {
                                residue_lists.clear();
                                residue_lists.push_back(residues);
                                residue_lists.push_back(next);
                                GetPairList(residue_lists, r_pair_list);
                                if (!r_pair_list.empty()) {
                                     for (std::vector<std::vector<RCSB::Residue*> >::const_iterator rpos = r_pair_list.begin(); rpos != r_pair_list.end(); ++rpos) {
                                          if (SeqCodeUtil::is_standard_aa_residue_plus_MSE((*rpos)[0]->ResName()) &&
                                              SeqCodeUtil::is_standard_aa_residue_plus_MSE((*rpos)[1]->ResName()))
                                               c_n_bonds.insert((*rpos)[0]->index());
                                          else {
                                               std::string index = String::IntToString((*mpos)->index()) + "_" + String::IntToString((*rpos)[0]->index());
                                               std::map<std::string, std::set<std::string> >::const_iterator bpos = bonded_atoms.find(index);
                                               if (bpos != bonded_atoms.end() && bpos->second.find("C") != bpos->second.end())
                                                    c_n_bonds.insert((*rpos)[0]->index());
                                          }
                                     }
                                }
                           } else {
                                for (std::vector<RCSB::Residue*>::const_iterator pos = residues.begin(); pos != residues.end(); ++pos) {
                                     c_n_bonds.insert((*pos)->index());
                                }
                           }
                      }

                      for (std::vector<RCSB::Residue*>::iterator pos = residues.begin(); pos != residues.end(); ++pos) {
                           if (chain->chain_type() == "ATOMP") {
                                bool has_c_n_bond = false;
                                if (c_n_bonds.find((*pos)->index()) != c_n_bonds.end()) has_c_n_bond = true;
                                // std::string index = (*pos)->pdb_chnid() + "_" + (*pos)->ResName() + "_" + (*pos)->pdb_res_no() + "_" + (*pos)->ins_code();
                                // if (non_standard_linkage.find(index) != non_standard_linkage.end()) is_connect = 2;
                                if (j == 0) {
                                     if (has_c_n_bond) {
                                          if ((*pos)->ResName() == "CYS")
                                               (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, CYS_N_Terminal_Leaving_Atoms,
                                                    AA_N_Terminal_Allowed_Atoms, AA_N_Terminal_Not_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                          else if ((*pos)->ResName() == "HIS")
                                               (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, HIS_N_Terminal_Leaving_Atoms,
                                                    AA_N_Terminal_Allowed_Atoms, AA_N_Terminal_Not_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                          else (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, AA_N_Terminal_Leaving_Atoms,
                                                    AA_N_Terminal_Allowed_Atoms, AA_N_Terminal_Not_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                     } else {
                                          if ((*pos)->ResName() == "CYS")
                                               (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, CYS_N_Terminal_Leaving_Atoms,
                                                    AA_Not_Allowed_Atoms, EmptyAtoms, _exclude_extra_hydrogen_flag);
                                          else if ((*pos)->ResName() == "HIS")
                                               (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, HIS_N_Terminal_Leaving_Atoms,
                                                    AA_Not_Allowed_Atoms, EmptyAtoms, _exclude_extra_hydrogen_flag);
                                          else (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, AA_N_Terminal_Leaving_Atoms,
                                                    AA_Not_Allowed_Atoms, EmptyAtoms, _exclude_extra_hydrogen_flag);
                                     }
                                } else if (j == (chain->SeqLen() - 1)) {
                                     if ((*pos)->ResName() == "CYS") 
                                          (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, CYS_C_Terminal_Leaving_Atoms,
                                               EmptyAtoms, AA_N_Terminal_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                     else if ((*pos)->ResName() == "HIS") 
                                          (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, HIS_C_Terminal_Leaving_Atoms,
                                               EmptyAtoms, AA_N_Terminal_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                     else (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, AA_N_Terminal_Allowed_Atoms,
                                               EmptyAtoms, AA_N_Terminal_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                } else {
                                     if (has_c_n_bond) {
                                          if ((*pos)->ResName() == "CYS")
                                               (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, CYS_Leaving_Atoms,
                                                    EmptyAtoms, AA_Not_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                          else if ((*pos)->ResName() == "HIS")
                                               (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, HIS_Leaving_Atoms,
                                                    EmptyAtoms, AA_Not_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                          else (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, AA_Leaving_Atoms,
                                                    EmptyAtoms, AA_Not_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                     } else {
                                          if ((*pos)->ResName() == "CYS")
                                               (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, CYS_Leaving_Atoms,
                                                    AA_N_Terminal_Not_Allowed_Atoms, AA_N_Terminal_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                          else if ((*pos)->ResName() == "HIS")
                                               (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, HIS_Leaving_Atoms,
                                                    AA_N_Terminal_Not_Allowed_Atoms, AA_N_Terminal_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                          else (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, AA_Leaving_Atoms,
                                                    AA_N_Terminal_Not_Allowed_Atoms, AA_N_Terminal_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                     }
                                }
                           } else if (chain->chain_type() == "ATOMN") {
                                if (j == 0)
                                     (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, NA_5_Terminal_Leaving_Atoms,
                                              NA_5_Terminal_Allowed_Atoms, NA_5_Terminal_Not_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                else if (j == (chain->SeqLen() - 1))
                                     (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, NA_C_Terminal_Leaving_Atoms,
                                              EmptyAtoms, NA_5_Terminal_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                                else (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, NA_Leaving_Atoms,
                                              EmptyAtoms, NA_Not_Allowed_Atoms, _exclude_extra_hydrogen_flag);
                           } else {
                                (*pos)->find_missing_or_extra_atoms(is_connect, need_hydrogen, EmptyAtoms, EmptyAtoms, EmptyAtoms,
                                                                    _exclude_extra_hydrogen_flag);
                           }
                           std::string index = String::IntToString((*mpos)->index()) + "_" + String::IntToString((*pos)->index());
                           std::map<std::string, std::set<std::string> >::const_iterator bpos = bonded_atoms.find(index);
                           if (bpos != bonded_atoms.end()) (*pos)->refine_missing_atoms(bpos->second);
                      }
                 }
                 chain = (*mpos)->GetNextChain();
            }
       }
}

void ValidateObj::_checking_dbref(Block& block)
{
       ISTable *t = getTablePtr(block, "entity_poly");
       if (!t) return;

       std::set<std::string> polymer_chain_set;
       polymer_chain_set.clear();

       std::string cs;
       std::vector<std::string> data;

       bool has_pdbx_target_identifier = false;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs, t, i, "pdbx_target_identifier");
            if (!cs.empty()) has_pdbx_target_identifier = true;
            get_value_clean(cs, t, i, "pdbx_strand_id");
            if (cs.empty()) continue;
            get_wordarray(data, cs, ",");
            for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
                 polymer_chain_set.insert(*pos);
            }
       }
       if (has_pdbx_target_identifier) {
            ISTable *t1 = getTablePtr(block, "pdbx_SG_project");
            if (!t1) _error_messages.insertMessage("sg_target_identifier_print", "warning",
               "Target DB filled in, but not a structural genomics entry. See _entity_poly.pdbx_target_identifer.");
       }

       if (polymer_chain_set.empty()) return;

       t = getTablePtr(block, "struct_ref");
       if (!t) _error_messages.insertMessage("dbref_print", "error", "Missing 'struct_ref' category.");

       std::set<std::string> dbref_chain_set;
       dbref_chain_set.clear();
       t = getTablePtr(block, "struct_ref_seq");
       if (t) {
            rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "pdbx_strand_id");
                 if (!cs.empty()) dbref_chain_set.insert(cs);
            }
       } else {
            _error_messages.insertMessage("dbref_print", "error", "Missing 'struct_ref_seq' category.");
            return;
       }

       for (std::set<std::string>::const_iterator spos = polymer_chain_set.begin(); spos != polymer_chain_set.end(); ++spos) {
            if (dbref_chain_set.find(*spos) != dbref_chain_set.end()) continue;
             _error_messages.insertMessage("dbref_print", "error", "Chain " + *spos + " does not present in DBREF mapping.");
       }
}

bool ValidateObj::_checking_assembly(Block& block)
{
       std::vector<std::string> categories;
       categories.clear();
       categories.push_back("pdbx_struct_assembly");
       categories.push_back("pdbx_struct_assembly_gen");
       categories.push_back("pdbx_struct_oper_list");

       for (std::vector<std::string>::const_iterator pos = categories.begin(); pos != categories.end(); ++pos) {
            ISTable *t = getTablePtr(block, *pos);
            if (t) continue;
            _error_messages.insertMessage("assembly_print", "warning", "Missing " + *pos + " category.");
       }

       bool flag = false;

       std::string cs, cs1;
       std::vector<std::string> data, data1;
       // map.first:  assembly_id
       // map.second: oper_expression list
       std::map<std::string, std::vector<std::string> > assembly_oper_id_mapping;
       assembly_oper_id_mapping.clear();

       ISTable *t1 = getTablePtr(block, "pdbx_struct_assembly_gen");
       if (t1 && !_molecules.empty()) {
            std::set<std::string> biol_assembly_asym_ids, missing_asym_ids;
            biol_assembly_asym_ids.clear();
            missing_asym_ids.clear();

            for (unsigned int i = 0; i < t1->GetNumRows(); ++i) {
                 get_value_clean(cs, t1, i, "assembly_id");
                 if (!cs.empty()) {
                      std::map<std::string, std::vector<std::string> >::iterator mpos = assembly_oper_id_mapping.find(cs);
                      if (mpos == assembly_oper_id_mapping.end()) {
                           data.clear();
                           assembly_oper_id_mapping.insert(std::make_pair(cs, data));
                           mpos = assembly_oper_id_mapping.find(cs);
                      }
                      get_value_clean(cs1, t1, i, "oper_expression");
                      if (!cs1.empty()) mpos->second.push_back(cs1);
                 }
                 get_value_clean(cs, t1, i, "asym_id_list");
                 if (cs.empty()) continue;
                 get_wordarray(data, cs, ",");
                 for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
                      biol_assembly_asym_ids.insert(*pos);
                 }
            }

            RCSB::Chain* chain = _molecules[0]->GetFirstChain();
            while (chain) {
                 if (biol_assembly_asym_ids.find(chain->ChainID()) == biol_assembly_asym_ids.end()) {
                      missing_asym_ids.insert(chain->ChainID());
                 }
                 chain = _molecules[0]->GetNextChain();
            }

            if (!missing_asym_ids.empty()) {
                 cs.clear();
                 for (std::set<std::string>::const_iterator spos = missing_asym_ids.begin(); spos != missing_asym_ids.end(); ++spos) {
                      if (!cs.empty()) cs += ",";
                      cs += *spos;
                 }

                 std::string error = "Asym ";
                 if (missing_asym_ids.size() > 1)
                      error += "IDs '" + cs + "' are";
                 else error += "ID '" + cs + "' is";
                 error += " not defined in biological assembly.";
                 _error_messages.insertMessage("assembly_print", "warning", error);
            }
       }

       std::set<std::string> biol_assembly_type_set, oper_id_set, auth_assembly_id_set, missing_id_set;
       biol_assembly_type_set.clear();
       biol_assembly_type_set.insert("author_and_software_defined_assembly");
       biol_assembly_type_set.insert("author_defined_assembly");
       biol_assembly_type_set.insert("complete icosahedral assembly");
       biol_assembly_type_set.insert("complete point assembly");
       biol_assembly_type_set.insert("representative helical assembly");
       biol_assembly_type_set.insert("software_defined_assembly");
       bool found_biol_assembly_type = false;

       oper_id_set.clear();
       ISTable *operT = getTablePtr(block, "pdbx_struct_oper_list");
       if (operT) {
            for (unsigned int i = 0; i < operT->GetNumRows(); ++i) {
                 get_value_clean(cs, operT, i, "id");
                 if (!cs.empty()) oper_id_set.insert(cs);
            }
       }

       auth_assembly_id_set.clear();
       ISTable *t = getTablePtr(block, "pdbx_struct_assembly_auth_evidence");
       if (t) {
            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean(cs, t, i, "assembly_id");
                 if (cs.empty()) continue;

                 get_value_clean_lower(cs1, t, i, "experimental_support");
                 if (!cs1.empty() && (cs1 != "none")) auth_assembly_id_set.insert(cs);
            }
       }

       t = getTablePtr(block, "pdbx_struct_assembly");
       if (t) {
            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean(cs, t, i, "id");
                 if (cs.empty()) continue;

                 get_value_clean(cs1, t, i, "details");
                 if (biol_assembly_type_set.find(cs1) != biol_assembly_type_set.end()) found_biol_assembly_type = true;
                 if ((cs1 == "software_defined_assembly") && (auth_assembly_id_set.find(cs) != auth_assembly_id_set.end())) {
                      _error_messages.insertMessage("assembly_print", "warning", "For '_pdbx_struct_assembly.id=" + cs + 
                            "' it has experimental evidence provided by author but is marked as 'software_defined_assembly' only.");
                 }

                 std::map<std::string, std::vector<std::string> >::const_iterator mpos = assembly_oper_id_mapping.find(cs);      
                 if (mpos != assembly_oper_id_mapping.end()) {
                      missing_id_set.clear();
                      for (std::vector<std::string>::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                           data.clear();
                           if (parseString(*vpos, data) != 0) continue;

                           for (std::vector<std::string>::const_iterator vpos1 = data.begin(); vpos1 != data.end(); ++vpos1) {
                                get_wordarray(data1, *vpos1, " ");
                                for (std::vector<std::string>::const_iterator vpos2 = data1.begin(); vpos2 != data1.end(); ++vpos2) {
                                     if (oper_id_set.find(*vpos2) == oper_id_set.end()) missing_id_set.insert(*vpos2);
                                }
                           }
                      }
                      if (missing_id_set.empty()) continue;

                      cs = "The following value(s) [ " + join_string(missing_id_set, ",")
                         + " ] from '_pdbx_struct_assembly_gen.oper_expression' item for assembly ID '" + cs
                         + "' is(are) not defined in '_pdbx_struct_oper_list.id' item.";
                      _error_messages.insertMessage("assembly_error_print", "error", cs); 
                 } else {
                      _error_messages.insertMessage("assembly_error_print", "error", "Missing 'pdbx_struct_assembly_gen' information for assembly ID '" + cs + "'.");
                 }
            }
       }
       if (!found_biol_assembly_type) {
            cs = "Missing biological assembly information. The '_pdbx_struct_assembly.details' value needs to be one of the following values:\n";
            cs += "[ 'author_and_software_defined_assembly', 'author_defined_assembly', 'complete icosahedral assembly', 'complete point assembly',\n";
            cs += "'representative helical assembly', 'software_defined_assembly' ].";
            _error_messages.insertMessage("assembly_print", "warning", cs);
       }

       if (operT) {
            std::string identity_id = "";
            for (unsigned int i = 0; i < operT->GetNumRows(); ++i) {
                 get_value_clean_lower(cs, operT, i, "type");
                 if (cs == "identity operation") {
                      get_value_clean(identity_id, operT, i, "id");
                      break;
                 }
            }
            if (identity_id.empty()) {
                 _error_messages.insertMessage("assembly_print", "warning", "Missing unit matrix in biological assembly.");
            } else if (t1) {
                 std::set<std::string> data_set, data_set_1;
                 std::map<std::string, std::set<std::string> > tmp_map;
                 // first key: assembly_id
                 // second key: asym_id
                 // value: _pdbx_struct_oper_list.id set
                 std::map<std::string, std::map<std::string, std::set<std::string> > > mapping;
                 mapping.clear();

                 for (unsigned int i = 0; i < t1->GetNumRows(); ++i) {
                      get_value_clean(cs, t1, i, "assembly_id");
                      if (cs.empty()) continue;
                      get_value_clean(cs1, t, i, "details");
                      // skip the unit matrix checking for non-biological assemblies
                      if (biol_assembly_type_set.find(cs1) == biol_assembly_type_set.end()) continue;
                      get_value_clean(cs1, t1, i, "oper_expression");
                      if (cs1.empty()) continue;
                      get_wordset(data_set, cs1, ",");
                      get_value_clean(cs1, t1, i, "asym_id_list");
                      if (cs1.empty()) continue;
                      get_wordset(data_set_1, cs1, ",");
                      std::map<std::string, std::map<std::string, std::set<std::string> > >::iterator mpos = mapping.find(cs);
                      if (mpos == mapping.end()) {
                           tmp_map.clear();
                           for (std::set<std::string>::const_iterator spos = data_set_1.begin(); spos != data_set_1.end(); ++spos) {
                                tmp_map.insert(std::make_pair(*spos, data_set));
                           }
                           mapping.insert(std::make_pair(cs, tmp_map));
                      } else {
                           for (std::set<std::string>::const_iterator spos = data_set_1.begin(); spos != data_set_1.end(); ++spos) {
                                std::map<std::string, std::set<std::string> >::iterator mpos1 = mpos->second.find(*spos);
                                if (mpos1 == mpos->second.end())
                                     mpos->second.insert(std::make_pair(*spos, data_set));
                                else {
                                     for (std::set<std::string>::const_iterator spos1 = data_set.begin(); spos1 != data_set.end(); ++spos1) {
                                          mpos1->second.insert(*spos1);
                                     }
                                }
                           }
                      }
                 }

                 for (std::map<std::string, std::map<std::string, std::set<std::string> > >::const_iterator
                      mpos = mapping.begin(); mpos != mapping.end(); ++mpos) {
                      data_set.clear();
                      for (std::map<std::string, std::set<std::string> >::const_iterator
                           mpos1 = mpos->second.begin(); mpos1 != mpos->second.end(); ++mpos1) {
                           if (mpos1->second.find(identity_id) != mpos1->second.end()) continue;
                           data_set.insert(mpos1->first);
                      }
                      if (data_set.empty()) continue;
                      if (data_set.size() == mpos->second.size()) {
                           _error_messages.insertMessage("assembly_print", "warning", "Missing unit matrix for assembly '" + mpos->first + "'.");
                           flag = true;
                      }
                 }
            }
       }

       return flag;
}

std::string ValidateObj::_write_atom_list(const std::list<RCSB::Atom*>& atom_list)
{
       if (atom_list.empty()) return "";

       NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       std::set<std::string> affectedFieldSet;
       affectedFieldSet.clear();
       for (std::list<RCSB::Atom*>::const_iterator lpos = atom_list.begin(); lpos != atom_list.end(); ++lpos) {
            if ((*lpos)->pdb_chnid().size() > 1) affectedFieldSet.insert("Strand_ID");
            if ((*lpos)->pdb_resnam().size() > 3) affectedFieldSet.insert("Residue_Name");
       }

       PdbWrite writer;
       writer.setLog(_logIo);
       writer.setCCDic(_ccDic);

       if (!affectedFieldSet.empty()) {
            for (std::set<std::string>::const_iterator spos = affectedFieldSet.begin(); spos != affectedFieldSet.end(); ++spos) {
                 writer.setFieldLabel(*spos);
            }
       }

       std::vector<std::string> records, FieldInfo;
       std::string content = "";
       for (std::list<RCSB::Atom*>::const_iterator lpos = atom_list.begin(); lpos != atom_list.end(); ++lpos) {
            writer.WriteAtom(records, "ATOMP", *lpos, FieldInfo, false);
            content += records[0].substr(11) + "\n";
       }

       NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);

       return content;
}

std::string ValidateObj::_check_cell_info(Block& block, const bool& comment_flag)
{
       std::string cs, error;
       error.clear();

       ISTable *t = getTablePtr(block, "cell");
       if (t) {
            std::vector<std::pair<std::string, string> > check_values;
            check_values.clear();
            check_values.push_back(std::make_pair("length_a", "1.0"));
            check_values.push_back(std::make_pair("length_b", "1.0"));
            check_values.push_back(std::make_pair("length_c", "1.0"));
            check_values.push_back(std::make_pair("angle_alpha", "90.0"));
            check_values.push_back(std::make_pair("angle_beta", "90.0"));
            check_values.push_back(std::make_pair("angle_gamma", "90.0"));
            check_values.push_back(std::make_pair("Z_PDB", "1.0"));

            for (std::vector<std::pair<std::string, string> >::const_iterator pos = check_values.begin(); pos != check_values.end(); ++pos) {
                 get_value_clean(cs, t, 0, pos->first);
                 if (cs.empty()) continue;
                 bool ok = true;
                 if (!String::IsNumber(cs)) ok = false;
                 else {
                      double val = atof(cs.c_str());
                      double std = atof(pos->second.c_str());
                      if (fabs(val - std) > 0.0001) ok = false;
                 }
                 if (ok) continue;
                 if (!error.empty()) error += "\n";
                 error += "The value '" + cs + "' in _cell." + pos->first + " is not " + pos->second  + ".";
            }
       }

       t = getTablePtr(block, "symmetry");
       if (t) {
            get_value_clean(cs, t, 0, "space_group_name_H-M");
            if (!cs.empty() && cs != "P 1") {
                 if (!error.empty()) error += "\n";
                 error += "The space group is '" + cs + "'.";
            }
       }

       if (!error.empty() and comment_flag) {
            error = "Structure models solved by non-diffraction method should have nominal unit cell parameters a=b=c=1, alpha=beta=gamma=90\n";
            error += "and space group='P 1'. But it is not the case in the deposited model. We have changed the parameters to nominal values.";
       }

       return error;
}

std::string ValidateObj::_write_link_list(std::list<std::vector<std::string> >& links)
{
       if (links.empty()) return "";

       std::vector<std::string> data, leaving_atom_flag_value, link_ID;
       std::set<std::string> affectedFieldSet;

       PdbWrite writer;
       writer.setLog(_logIo);
       writer.setCCDic(_ccDic);
       writer.setOutputFormat(NDB_FILE_FORMAT_PDB);

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("LINK");

       affectedFieldSet.clear();
       writer.getAffectedFieldSet(ndbformat, links, affectedFieldSet);
       if (affectedFieldSet.empty()) {
            writer.writeGeneralRecords(data, ndbformat, links);
       } else {
            ndb_token_format updated_ndbformat = ndbformat;
            writer.updateNdbTokenFormat(affectedFieldSet, updated_ndbformat);
            writer.writeGeneralRecords(data, updated_ndbformat, links);
       }

       if (data.empty()) return "";

       unsigned int link_id_size = 0;
       leaving_atom_flag_value.clear();
       link_ID.clear();
       for (std::list<std::vector<std::string> >::const_iterator lpos = links.begin(); lpos != links.end(); ++lpos) {
            if (lpos->size() == 18) {
                 leaving_atom_flag_value.push_back((*lpos)[16]);
                 link_ID.push_back((*lpos)[17]);
                 if ((*lpos)[17].size() > link_id_size) link_id_size = (*lpos)[17].size();
            }
       }

       std::string text = "";
       for (unsigned int i = 0; i < data.size(); ++i) {
            if (!text.empty()) text += "\n";
            if ((data.size() == link_ID.size()) && (link_id_size > 0)) text += FormattedString(link_ID[i], link_id_size + 3, true, false);
            text += data[i];
            if (data.size() == leaving_atom_flag_value.size()) text += "    " + leaving_atom_flag_value[i];
       }

       return text;
}

std::string ValidateObj::_get_zero_occupancy_atom_list_summary(std::list<RCSB::Atom*>& zero_occupancy_atom_list)
{
       std::list<RCSB::Atom*> non_hydrogen_zero_occupancy_atom_list;
       non_hydrogen_zero_occupancy_atom_list.clear();
       std::map<std::string, int> zero_occupancy_atom_count;
       zero_occupancy_atom_count.clear();
       for (std::list<RCSB::Atom*>::const_iterator apos = zero_occupancy_atom_list.begin(); apos != zero_occupancy_atom_list.end(); ++apos) {
            std::map<std::string, int>::iterator mpos = zero_occupancy_atom_count.find((*apos)->atom_type());
            if (mpos != zero_occupancy_atom_count.end()) mpos->second += 1;
            else zero_occupancy_atom_count.insert(std::make_pair((*apos)->atom_type(), 1)); 
            if (((*apos)->atom_type() != "H") && ((*apos)->atom_type() != "D")) non_hydrogen_zero_occupancy_atom_list.push_back(*apos);
       }
       std::string summary = "";
       for (std::map<std::string, int>::const_iterator mpos = zero_occupancy_atom_count.begin(); mpos != zero_occupancy_atom_count.end(); ++mpos) {
            if (!summary.empty()) summary += "\n";
            summary += mpos->first + " atoms with zero occupancy: " + String::IntToString(mpos->second);
       }
       zero_occupancy_atom_list = non_hydrogen_zero_occupancy_atom_list;
       return summary;
}

static void insert_bonded_atoms(std::map<std::string, std::set<std::string> >& bonded_atoms, const int& mol_index,
                                const int& res_index, const std::string& atom_name)
{
      std::string index = String::IntToString(mol_index) + "_" + String::IntToString(res_index);
      std::map<std::string, std::set<std::string> >::iterator mpos = bonded_atoms.find(index);
       if (mpos != bonded_atoms.end())
            mpos->second.insert(atom_name);
       else {
            std::set<std::string> atom_set;
            atom_set.clear();
            atom_set.insert(atom_name);
            bonded_atoms.insert(std::make_pair(index, atom_set));
       }
}
