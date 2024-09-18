/*
FILE:     AnnotationObj_IO_Util.C
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

#include "AnnotationObj.h"
#include "CategoryMapping.h"
#include "SeqCodeUtil.h"
#include "SgCenter.h"
#include "utillib.h"

void AnnotationObj::_read_modres(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_struct_mod_residue");
       if (!t) return;

       _MODRES modres;
       std::string cs;
       std::vector<std::string> data, items;
       std::vector<std::vector<std::string> > names;
       std::vector<RCSB::Residue*> residues;

       items.clear();
       items.push_back("auth_asym_id");
       items.push_back("auth_comp_id");
       items.push_back("auth_seq_id");
       items.push_back("PDB_ins_code");

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            names.clear();
            data.clear();
            for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
                 get_value_clean(cs, t, i, *pos);
                 data.push_back(cs);
            }
            names.push_back(data);

            _getResiduesfromNames(0, names, modres.mol_index, residues, data);
            if (residues.empty()) continue;

            modres.res_index = residues[0]->index();
            get_value_clean(modres.Standard_Name, t, i, "parent_comp_id");
            get_value_clean(modres.details, t, i, "details");
            _modres.push_back(modres);
       }
}

void AnnotationObj::_write_modres(Block& block)
{
       if (_modres.empty() || _molecules.empty()) {
            deleteTable(block, "pdbx_struct_mod_residue");
            return;
       }

       ISTable *t = _newTablePtr("pdbx_struct_mod_residue");

       int row = 0;
       for (std::list<_MODRES>::const_iterator lpos = _modres.begin(); lpos != _modres.end(); ++lpos) {
            bool is_removed = false;
            RCSB::Residue* res = _molecules[0]->find_residue(lpos->res_index, is_removed);
            if (!res) continue;
            t->AddRow();
            t->UpdateCell(row, "id", String::IntToString(row + 1));
            t->UpdateCell(row, "label_asym_id", res->chnid());
            t->UpdateCell(row, "label_comp_id", res->ResName());
            t->UpdateCell(row, "label_seq_id", res->res_no());
            t->UpdateCell(row, "auth_comp_id", res->ResName());
            t->UpdateCell(row, "auth_asym_id", res->pdb_chnid());
            t->UpdateCell(row, "auth_seq_id", res->pdb_res_no());
            t->UpdateCell(row, "PDB_ins_code", res->ins_code());
            t->UpdateCell(row, "parent_comp_id", lpos->Standard_Name);
            t->UpdateCell(row, "details", lpos->details);
            row++;
       }

       if (row == 0) {
            deleteTable(block, "pdbx_struct_mod_residue");
            delete t;
       } else block.WriteTable(t);
}

void AnnotationObj::_read_prd_features(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_molecule_features");
       if (!t) return;

       const std::vector<std::string>& itemNames = t->GetColumnNames();

       std::map<std::string, std::string> tmp_map;
       std::string prd_id, cs;

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(prd_id, t, i, "prd_id");
            if (prd_id.empty()) continue;
            tmp_map.clear();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 if (*pos == "prd_id") continue;
                 get_value(cs, t, i, *pos);
                 if (cs.empty()) continue;
                 tmp_map.insert(std::make_pair(*pos, cs));
            }
            if (tmp_map.empty()) continue;
            _prd_features.insert(std::make_pair(prd_id, tmp_map));
       }
}

void AnnotationObj::_write_prd_features(Block& block)
{
       if (_prd_features.empty()) {
            deleteTable(block, "pdbx_molecule_features");
            return;
       }

       ISTable *t = _newTablePtr("pdbx_molecule_features");
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       int row = 0;
       for (std::map<std::string, std::map<std::string, std::string> >::const_iterator mpos = _prd_features.begin(); mpos != _prd_features.end(); ++mpos) {
            t->AddRow();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 if (*pos == "prd_id") t->UpdateCell(row, "prd_id", mpos->first);
                 else {
                      std::map<std::string, std::string>::const_iterator mmpos = mpos->second.find(*pos);
                      if (mmpos != mpos->second.end()) t->UpdateCell(row, *pos, mmpos->second);
                 }
            }
            row++;
       }

       block.WriteTable(t);
}

void AnnotationObj::_read_prd_instances(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "pdbx_molecule");
       if (!t) return;

       std::set<int> chain_set, residue_set;
       chain_set.clear();
       residue_set.clear();

       std::string prd_id = "";
       std::string pre_instance_id = "";

       std::string instance_id, asym_id;

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(instance_id, t, i, "instance_id");
            if (instance_id != pre_instance_id) {
                 if (!prd_id.empty() && (!chain_set.empty() || !residue_set.empty())) {
                      _prd_instances.push_back(std::make_pair(prd_id, std::make_pair(chain_set, residue_set)));
                 }
                 prd_id.clear();
                 chain_set.clear();
                 residue_set.clear();
                 pre_instance_id = instance_id;
            }

            get_value_clean(prd_id, t, i, "prd_id");
            get_value_clean(asym_id, t, i, "asym_id");
            if (prd_id.empty() || asym_id.empty()) continue;

            RCSB::Chain* chain = _molecules[0]->GetAsymChain(asym_id);
            if (!chain) {
                 chain = _molecules[0]->GetOrigAsymSugarChain(asym_id);
                 if (!chain) {
                      std::string error = "Asym ID " + asym_id + " does not exist in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                      _logIo->message(error.c_str());
                      continue;
                 }
            }

            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN")
                 chain_set.insert(chain->index());
            else {
                 RCSB::Residue* residue = chain->GetFirstResidue();
                 while (residue) {
                      residue_set.insert(residue->index());
                      residue = chain->GetNextResidue();
                 }
            }
       }

       if (!prd_id.empty() && (!chain_set.empty() || !residue_set.empty())) {
            _prd_instances.push_back(std::make_pair(prd_id, std::make_pair(chain_set, residue_set)));
       }
}

void AnnotationObj::_write_prd_instances(Block& block)
{
       if (_prd_instances.empty() || _molecules.empty()) {
            deleteTable(block, "pdbx_molecule");
            return;
       }

       ISTable *t = _newTablePtr("pdbx_molecule");
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       bool has_auth_asym_id = false;
       bool has_auth_seq_id = false;
       for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
            if (*pos == "auth_asym_id") has_auth_asym_id = true;
            if (*pos == "auth_seq_id") has_auth_seq_id = true;
       }

       std::set<std::string> unique_asym_id_set;
       unique_asym_id_set.clear();

       int row = 0;
       int instance_id = 0;
       for (std::vector<std::pair<std::string, std::pair<std::set<int>, std::set<int> > > >::const_iterator
            pos = _prd_instances.begin(); pos != _prd_instances.end(); ++pos) {
            bool found_info = false;
            for (std::set<int>::const_iterator spos = pos->second.first.begin(); spos != pos->second.first.end(); ++spos) {
                 bool is_removed = false;
                 RCSB::Chain* chain = _molecules[0]->GetIndexChain(*spos, is_removed);
                 if (!chain) {
                      if (!is_removed) {
                           std::string error = "Chain index " + String::IntToString(*spos) + " does not exist in Model "
                                             + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                           _logIo->message(error.c_str());
                      }
                      continue;
                 }
                 found_info = true;
            }
            for (std::set<int>::const_iterator spos = pos->second.second.begin(); spos != pos->second.second.end(); ++spos) {
                 bool is_removed = false;
                 RCSB::Residue* residue = _molecules[0]->find_residue(*spos, is_removed);
                 if (!residue) {
                      if (!is_removed) {
                           std::string error = "Residue index " + String::IntToString(*spos) + " does not exist in Model "
                                             + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                           _logIo->message(error.c_str());
                      }
                      continue;
                 }
                 found_info = true;
            }
            if (!found_info) continue;

            instance_id++;
            std::string instance_ID = String::IntToString(instance_id);
            for (std::set<int>::const_iterator spos = pos->second.first.begin(); spos != pos->second.first.end(); ++spos) {
                 bool is_removed = false;
                 RCSB::Chain* chain = _molecules[0]->GetIndexChain(*spos, is_removed);
                 if (!chain) continue;

                 if (unique_asym_id_set.find(chain->ChainID()) != unique_asym_id_set.end()) continue;
                 unique_asym_id_set.insert(chain->ChainID());

                 t->AddRow();
                 t->UpdateCell(row, "instance_id", instance_ID);
                 t->UpdateCell(row, "prd_id", pos->first);
                 t->UpdateCell(row, "asym_id", chain->ChainID());
                 if (has_auth_asym_id) t->UpdateCell(row, "auth_asym_id", chain->PDB_ChainID());
                 row++;
            }
            for (std::set<int>::const_iterator spos = pos->second.second.begin(); spos != pos->second.second.end(); ++spos) {
                 bool is_removed = false;
                 RCSB::Residue* residue = _molecules[0]->find_residue(*spos, is_removed);
                 if (!residue) continue;

                 if (unique_asym_id_set.find(residue->chnid()) != unique_asym_id_set.end()) continue;
                 unique_asym_id_set.insert(residue->chnid());

                 t->AddRow();
                 t->UpdateCell(row, "instance_id", instance_ID);
                 t->UpdateCell(row, "prd_id", pos->first);
                 t->UpdateCell(row, "asym_id", residue->chnid());
                 if (has_auth_asym_id) t->UpdateCell(row, "auth_asym_id", residue->pdb_chnid());
                 if (has_auth_seq_id) t->UpdateCell(row, "auth_seq_id", residue->pdb_res_no());
                 row++;
            }
       }

       block.WriteTable(t);
}

void AnnotationObj::_read_pdbx_entity_instance_feature(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "pdbx_entity_instance_feature");
       if (!t) return;

       std::map<std::string, RCSB::Residue*> first_residue_mapping;
       first_residue_mapping.clear();
       std::vector<RCSB::Residue*> residue_list;
       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            chain->GetFirstResidueList(residue_list);
            while (!residue_list.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                      if ((chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") &&
                           SeqCodeUtil::is_a_standard_residue((*rpos)->ResName())) continue;

                      if (!(*rpos)->OrigResName().empty() && first_residue_mapping.find((*rpos)->OrigResName()) == first_residue_mapping.end()) {
                           first_residue_mapping.insert(std::make_pair((*rpos)->OrigResName(), *rpos));
                      }
                      if (first_residue_mapping.find((*rpos)->ResName()) == first_residue_mapping.end()) {
                           first_residue_mapping.insert(std::make_pair((*rpos)->ResName(), *rpos));
                      }
                 }
                 chain->GetNextResidueList(residue_list);
            }
            chain = _molecules[0]->GetNextChain();
       }

       std::vector<std::string> items, data, comp_ids;
       items.clear();
       items.push_back("comp_id");
       items.push_back("auth_comp_id");
       items.push_back("auth_asym_id");
       items.push_back("auth_seq_num");
       items.push_back("feature_type");
       items.push_back("details");

       std::string cs;
       RCSB::Residue* res = NULL;
       _PDBX_ENTITY_INSTANCE_FEATURE instance;

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            data.clear();
            for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
                 get_value_clean(cs, t, i, *pos);
                 data.push_back(cs);
            }
            comp_ids.clear();
            for (int j = 0; j < 2; ++j) {
                 if (data[j].empty()) continue;
                 if (!comp_ids.empty() && comp_ids[0] == data[j]) continue;
                 comp_ids.push_back(data[j]);
            }
            if (comp_ids.empty()) {
                 std::string error = "No comp_id defined in pdbx_entity_instance_feature category.\n";
                 _logIo->message(error.c_str());
                 continue;
            }
            if (!data[3].empty()) {
                 res = NULL;
                 for (std::vector<std::string>::const_iterator pos = comp_ids.begin(); pos != comp_ids.end(); ++pos) {
                      res = _molecules[0]->find_pdb_residue(data[2], *pos, data[3], "");
                      if (res) break;
                 }
                 if (!res) {
                      std::string error = "Residue " + data[2] + " " + comp_ids[comp_ids.size() - 1] + " " + data[3]
                                        + " does not exist in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                      _logIo->message(error.c_str());
                      continue;
                 }
                 // instance.residue_index = res->index();
                 instance.is_instance = true;
                 instance.feature_type = data[4];
                 instance.details = data[5];
                 instance.atoms = res->atoms();
                 _pdbx_entity_instance_feature.push_back(instance);
            } else {
                 res = NULL;
                 for (std::vector<std::string>::const_iterator pos = comp_ids.begin(); pos != comp_ids.end(); ++pos) {
                      std::map<std::string, RCSB::Residue*>::const_iterator mpos = first_residue_mapping.find(*pos);
                      if (mpos != first_residue_mapping.end()) {
                           res = mpos->second;
                           break;
                      }
                 }
                 if (!res) {
                      std::string error = "Residue " + comp_ids[comp_ids.size() - 1] + " does not exist in Model "
                                        + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                      _logIo->message(error.c_str());
                      continue;
                 }
                 // instance.residue_index = res->index();
                 instance.is_instance = false;
                 instance.feature_type = data[4];
                 instance.details = data[5];
                 instance.atoms = res->atoms();
                 _pdbx_entity_instance_feature.push_back(instance);
            }
       }
}

void AnnotationObj::_write_pdbx_entity_instance_feature(Block& block)
{
       if (_pdbx_entity_instance_feature.empty() || _molecules.empty()) {
            deleteTable(block, "pdbx_entity_instance_feature");
            return;
       }

       ISTable *t = _newTablePtr("pdbx_entity_instance_feature");

       std::set<std::string> unique_id_set;
       std::vector<std::vector<std::string> > instance_list;
       std::vector<std::string> data;

       int row = 0;
       for (std::list<_PDBX_ENTITY_INSTANCE_FEATURE>::const_iterator pos = _pdbx_entity_instance_feature.begin();
            pos != _pdbx_entity_instance_feature.end(); ++pos) {
/*
            bool is_removed = false;
            RCSB::Residue* res = _molecules[0]->find_residue(pos->residue_index, is_removed);
            if (!res) {
                 if (!is_removed) {
                      std::string error = "Residue index " + String::IntToString(pos->residue_index) + " does not exist in Model "
                                        + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                      _logIo->message(error.c_str());
                 }
                 continue;
            }
*/
            unique_id_set.clear();
            instance_list.clear();
            for (std::vector<RCSB::Atom*>::const_iterator apos = pos->atoms.begin(); apos != pos->atoms.end(); ++apos) {
                 std::string key = (*apos)->pdb_resnam();
                 if (pos->is_instance) key = (*apos)->pdb_resnam() + "_" + (*apos)->pdb_chnid() + "_" + (*apos)->pdb_resnum() + "_" + (*apos)->ins_code();
                 if (unique_id_set.find(key) != unique_id_set.end()) continue;

                 unique_id_set.insert(key);
                 data.clear();
                 data.push_back((*apos)->pdb_resnam());
                 if (pos->is_instance) {
                      data.push_back((*apos)->chnid());
                      data.push_back((*apos)->resnum());
                      data.push_back((*apos)->pdb_chnid());
                      data.push_back((*apos)->pdb_resnum());
                 } else {
                      data.push_back("");
                      data.push_back("");
                      data.push_back("");
                      data.push_back("");
                 }
                 instance_list.push_back(data);
            }
            if (instance_list.empty()) continue;

            for (std::vector<std::vector<std::string> >::const_iterator ipos = instance_list.begin(); ipos != instance_list.end(); ++ipos) {
                 t->AddRow();
                 t->UpdateCell(row, "ordinal", String::IntToString(row + 1));
                 t->UpdateCell(row, "comp_id", (*ipos)[0]);
                 t->UpdateCell(row, "auth_comp_id", (*ipos)[0]);
                 t->UpdateCell(row, "feature_type", pos->feature_type);
                 t->UpdateCell(row, "details", pos->details);
                 if (pos->is_instance) {
                      t->UpdateCell(row, "asym_id", (*ipos)[1]);
                      t->UpdateCell(row, "seq_num", (*ipos)[2]);
                      t->UpdateCell(row, "auth_asym_id", (*ipos)[3]);
                      t->UpdateCell(row, "auth_seq_num", (*ipos)[4]);
                 }
                 row++;
            }
       }

       block.WriteTable(t);
}

void AnnotationObj::_read_pdbx_remediation_atom_site_mapping(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "pdbx_remediation_atom_site_mapping");
       if (!t) return;

       std::vector<std::string> data, pre_items;
       pre_items.clear();
       pre_items.push_back("pre_group_PDB");
       pre_items.push_back("pre_auth_atom_id");
       pre_items.push_back("pre_auth_comp_id");
       pre_items.push_back("pre_auth_asym_id");
       pre_items.push_back("pre_auth_seq_id");
       pre_items.push_back("pre_auth_alt_id");
       pre_items.push_back("pre_PDB_ins_code");
       pre_items.push_back("pre_occupancy");

       std::vector<std::vector<std::string> > items, names;
       items.clear();

       data.clear();
       data.push_back("auth_asym_id");
       data.push_back("auth_comp_id");
       data.push_back("auth_seq_id");
       data.push_back("PDB_ins_code");
       data.push_back("auth_atom_id");
       data.push_back("auth_alt_id");
       items.push_back(data);

       int mol_index = 0;
       std::string cs;
       _PDBX_REMEDIATION_ATOM_SITE_MAPPING _atom_site_mapping;

       std::vector<RCSB::Atom*> atoms;
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_values(t, i, items, names);
            _getAtomsfromNames(0, names, mol_index, atoms, data, true);
            if (atoms.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 continue;
            }
            data.clear();
            for (std::vector<std::string>::const_iterator pos = pre_items.begin(); pos != pre_items.end(); ++pos) {
                 get_value_clean(cs, t, i, *pos);
                 if (cs.empty()) cs = ".";
                 data.push_back(cs);
            }
            _atom_site_mapping.atom = atoms[0];
            _atom_site_mapping.pre_data = data;
            _pdbx_remediation_atom_site_mapping.push_back(_atom_site_mapping);
       }
}

void AnnotationObj::_write_pdbx_remediation_atom_site_mapping(Block& block)
{
       if (_pdbx_remediation_atom_site_mapping.empty() || _molecules.empty()) return;

       ISTable *t = _newTablePtr("pdbx_remediation_atom_site_mapping");

       std::map<std::string, std::string> mapping;
       mapping.clear();
       mapping.insert(std::make_pair("group_PDB", "group_PDB"));
       mapping.insert(std::make_pair("label_atom_id", "label_atom_id"));
       mapping.insert(std::make_pair("label_comp_id", "label_comp_id"));
       mapping.insert(std::make_pair("label_asym_id", "label_asym_id"));
       mapping.insert(std::make_pair("label_seq_id", "label_seq_id"));
       mapping.insert(std::make_pair("label_alt_id", "label_alt_id"));
       mapping.insert(std::make_pair("auth_atom_id", "auth_atom_id"));
       mapping.insert(std::make_pair("auth_comp_id", "auth_comp_id"));
       mapping.insert(std::make_pair("auth_asym_id", "auth_asym_id"));
       mapping.insert(std::make_pair("auth_seq_id", "auth_seq_id"));
       mapping.insert(std::make_pair("auth_alt_id", "label_alt_id"));
       mapping.insert(std::make_pair("PDB_ins_code", "pdbx_PDB_ins_code"));
       mapping.insert(std::make_pair("occupancy", "occupancy"));

       std::vector<std::string> pre_items;
       pre_items.clear();
       pre_items.push_back("pre_group_PDB");
       pre_items.push_back("pre_auth_atom_id");
       pre_items.push_back("pre_auth_comp_id");
       pre_items.push_back("pre_auth_asym_id");
       pre_items.push_back("pre_auth_seq_id");
       pre_items.push_back("pre_auth_alt_id");
       pre_items.push_back("pre_PDB_ins_code");
       pre_items.push_back("pre_occupancy");

       int row = 0;
       for (std::list<_PDBX_REMEDIATION_ATOM_SITE_MAPPING>::const_iterator pos = _pdbx_remediation_atom_site_mapping.begin();
            pos != _pdbx_remediation_atom_site_mapping.end(); ++pos) {
            t->AddRow();
            t->UpdateCell(row, "id", String::IntToString(row + 1));
            for (std::map<std::string, std::string>::const_iterator mpos = mapping.begin(); mpos != mapping.end(); ++mpos) {
                 std::string cs = pos->atom->getValue(mpos->second);
                 if (cs.empty() || cs == "?") cs = ".";
                 t->UpdateCell(row, mpos->first, cs);
            }
            for (unsigned int i = 0; i < pre_items.size(); ++i) {
                 t->UpdateCell(row, pre_items[i], pos->pre_data[i]);
            }
            row++;
       }

       block.WriteTable(t);
}

void AnnotationObj::_read_refine_ls_restr_ncs_list(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "refine_ls_restr_ncs");
       if (!t) return;

       std::string cs, pdbx_asym_id;
       std::map<std::string, std::string> meta_data_map;

       const std::vector<std::string>& itemNames = t->GetColumnNames();

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            pdbx_asym_id.clear();
            meta_data_map.clear();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 if (*pos == "pdbx_auth_asym_id") continue;
                 get_value_clean(cs, t, i, *pos);
                 if (cs.empty()) continue;
                 if (*pos == "pdbx_asym_id") pdbx_asym_id = cs;
                 else meta_data_map.insert(std::make_pair(*pos, cs));
            }
            if (pdbx_asym_id.empty()) continue;
            RCSB::Chain* chain = _molecules[0]->GetAsymChain(pdbx_asym_id);
            if (!chain) continue;
            _refine_ls_restr_ncs_list.push_back(std::make_pair(meta_data_map, chain->index()));
       }
       if (_refine_ls_restr_ncs_list.size() != t->GetNumRows()) _refine_ls_restr_ncs_list.clear();
}

void AnnotationObj::_write_refine_ls_restr_ncs_list(Block& block)
{
       if (_refine_ls_restr_ncs_list.empty() || _molecules.empty()) return;

       ISTable *t = _newTablePtr("refine_ls_restr_ncs");
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       bool is_removed;
       for (std::list<std::pair<std::map<std::string, std::string>, int> >::const_iterator lpos = _refine_ls_restr_ncs_list.begin();
            lpos != _refine_ls_restr_ncs_list.end(); ++lpos) {
            RCSB::Chain* chain = _molecules[0]->GetIndexChain(lpos->second, is_removed);
            if (!chain) continue;
            int row = t->GetNumRows();
            t->AddRow();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 if (*pos == "pdbx_ordinal") t->UpdateCell(row, *pos, String::IntToString(row + 1));
                 else if (*pos == "pdbx_asym_id") t->UpdateCell(row, *pos, chain->ChainID());
                 else if (*pos == "pdbx_auth_asym_id") t->UpdateCell(row, *pos, chain->PDB_ChainID());
                 else {
                      std::map<std::string, std::string>::const_iterator mpos = lpos->first.find(*pos);
                      if (mpos != lpos->first.end()) t->UpdateCell(row, *pos, mpos->second);
                 }
            }
       }
       // if (t->GetNumRows() > 0) block.WriteTable(t);
       if (t->GetNumRows() == _refine_ls_restr_ncs_list.size()) block.WriteTable(t);
       else { delete t; } 
}

void AnnotationObj::_read_struct_ncs_dom_lim_list(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "struct_ncs_dom_lim");
       if (!t) return;

       std::vector<std::string> all_beg_items, all_end_items;
       std::vector<std::vector<std::string> > beg_items, end_items, values;

       _get_struct_ncs_dom_lim_item_names(t, all_beg_items, beg_items, all_end_items, end_items);

       std::set<std::string> skip_item_set;
       skip_item_set.clear();
       for (std::vector<std::vector<std::string> >::const_iterator vpos = beg_items.begin(); vpos != beg_items.end(); ++vpos) {
            for (std::vector<std::string>::const_iterator pos = vpos->begin(); pos != vpos->end(); ++pos) {
                 skip_item_set.insert(*pos);
            }
       }
       skip_item_set.insert("beg_label_alt_id");
       for (std::vector<std::vector<std::string> >::const_iterator vpos = end_items.begin(); vpos != end_items.end(); ++vpos) {
            for (std::vector<std::string>::const_iterator pos = vpos->begin(); pos != vpos->end(); ++pos) {
                 skip_item_set.insert(*pos);
            }
       }
       skip_item_set.insert("end_label_alt_id");

       const std::vector<std::string>& itemNames = t->GetColumnNames();

       std::string cs;
       std::map<std::string, std::string> meta_data_map;

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            RCSB::Residue* begRes = _find_struct_ncs_dom_lim_residue(t, i, beg_items, all_beg_items, cs);
/*
            if (!begRes && !cs.empty()) {
                 std::string err = "Residue ( " + cs + " ) can not be found.";
                 _logIo->message("%s\n", err.c_str());
            }
*/
            RCSB::Residue* endRes = _find_struct_ncs_dom_lim_residue(t, i, end_items, all_end_items, cs);
/*
            if (!endRes && !cs.empty()) {
                 std::string err = "Residue ( " + cs + " ) can not be found.";
                 _logIo->message("%s\n", err.c_str());
            }
*/
            if (!begRes || !endRes) continue;

            meta_data_map.clear();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 if (skip_item_set.find(*pos) != skip_item_set.end()) continue;
                 get_value_clean(cs, t, i, *pos);
                 if (cs.empty()) continue;
                 meta_data_map.insert(std::make_pair(*pos, cs));
            }
            _struct_ncs_dom_lim_list.push_back(std::make_pair(meta_data_map, std::make_pair(begRes->index(), endRes->index())));
       }
       if (_struct_ncs_dom_lim_list.size() != t->GetNumRows()) _struct_ncs_dom_lim_list.clear();
}

void AnnotationObj::_write_struct_ncs_dom_lim_list(Block& block)
{
       if (_struct_ncs_dom_lim_list.empty() || _molecules.empty()) return;

       ISTable *t = _newTablePtr("struct_ncs_dom_lim");
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       bool is_removed;
       for (std::list<std::pair<std::map<std::string, std::string>, std::pair<int, int> > >::const_iterator lpos = _struct_ncs_dom_lim_list.begin();
            lpos != _struct_ncs_dom_lim_list.end(); ++lpos) {
            RCSB::Residue* begRes = _molecules[0]->find_residue(lpos->second.first, is_removed);
            RCSB::Residue* endRes = _molecules[0]->find_residue(lpos->second.second, is_removed);
            if (!begRes || !endRes) continue;

            int row = t->GetNumRows();
            t->AddRow();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 std::string beg_label_seq_id = begRes->res_no(); if (beg_label_seq_id.empty()) beg_label_seq_id = ".";
                 std::string beg_label_alt_id = begRes->alt_loc(); if (beg_label_alt_id.empty()) beg_label_alt_id = ".";
                 std::string end_label_seq_id = endRes->res_no(); if (end_label_seq_id.empty()) end_label_seq_id = ".";
                 std::string end_label_alt_id = endRes->alt_loc(); if (end_label_alt_id.empty()) end_label_alt_id = ".";
                 if (*pos == "beg_label_asym_id") t->UpdateCell(row, *pos, begRes->chnid());
                 else if (*pos == "beg_label_comp_id") t->UpdateCell(row, *pos, begRes->ResName());
                 else if (*pos == "beg_label_seq_id") t->UpdateCell(row, *pos, beg_label_seq_id);
                 else if (*pos == "beg_label_alt_id") t->UpdateCell(row, *pos, beg_label_alt_id);
                 else if (*pos == "beg_auth_asym_id") t->UpdateCell(row, *pos, begRes->pdb_chnid());
                 else if (*pos == "beg_auth_comp_id") t->UpdateCell(row, *pos, begRes->ResName());
                 else if (*pos == "beg_auth_seq_id") t->UpdateCell(row, *pos, begRes->pdb_res_no());
                 else if (*pos == "end_label_asym_id") t->UpdateCell(row, *pos, endRes->chnid());
                 else if (*pos == "end_label_comp_id") t->UpdateCell(row, *pos, endRes->ResName());
                 else if (*pos == "end_label_seq_id") t->UpdateCell(row, *pos, end_label_seq_id);
                 else if (*pos == "end_label_alt_id") t->UpdateCell(row, *pos, end_label_alt_id);
                 else if (*pos == "end_auth_asym_id") t->UpdateCell(row, *pos, endRes->pdb_chnid());
                 else if (*pos == "end_auth_comp_id") t->UpdateCell(row, *pos, endRes->ResName());
                 else if (*pos == "end_auth_seq_id") t->UpdateCell(row, *pos, endRes->pdb_res_no());
                 else {
                      std::map<std::string, std::string>::const_iterator mpos = lpos->first.find(*pos);
                      if (mpos != lpos->first.end()) t->UpdateCell(row, *pos, mpos->second);
                 }
            }
       }
       // if (t->GetNumRows() > 0) block.WriteTable(t);
       if (t->GetNumRows() == _struct_ncs_dom_lim_list.size()) block.WriteTable(t);
       else { delete t; } 
}

void AnnotationObj::_write_entity_scheme_info(Block& block)
{
       if (_molecules.empty()) {
            deleteTable(block, "pdbx_poly_seq_scheme");
            deleteTable(block, "pdbx_branch_scheme");
            deleteTable(block, "pdbx_nonpoly_scheme");
            return;
       }

       _branch_residue_set.clear();

       ISTable *t3 = _newTablePtr("pdbx_poly_seq_scheme");
       ISTable *t4 = _newTablePtr("pdbx_branch_scheme");
       ISTable *t5 = _newTablePtr("pdbx_nonpoly_scheme");
       ISTable *t6 = _newTablePtr("pdbx_entity_branch_list");

       std::vector<std::pair<std::string, std::pair<std::string, std::string> > > data;
       std::vector<std::string> tmp_data;
       std::vector<RCSB::Residue*> residue_list;

       std::map<std::string, std::string> key_val_map, key_val_map1 /* , residue_entity_map, residue_chainid_map */ ;
       std::string code_can;
       std::set<std::string> entity_id_set;
       entity_id_set.clear();

       std::map<unsigned int, std::string> entity_order;
       entity_order.clear();
       std::map<std::string, std::vector<std::vector<std::vector<string> > > > found_entities;
       found_entities.clear();

       std::vector<std::vector<std::vector<string> > > entity_data;
       std::vector<std::vector<std::string> > coor_seq_data, residue_data;
       coor_seq_data.clear();

       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            bool new_entity = false;
            if (found_entities.find(chain->entity_id()) == found_entities.end()) {
                 entity_order.insert(std::make_pair(found_entities.size(), chain->entity_id()));
                 new_entity = true;
            }

            std::string sequence = "";
            std::string sequence_p = "";
            std::string sequence_c = "";
            std::string sequence_c_1 = "";
            int len = 0;
            int len_p = 0;
            int len_c = 0;
            int len_c_1 = 0;
            std::string nstd_monomer = "no";
            entity_data.clear();
            for (unsigned int i = 0; i < chain->SeqLen(); ++i) {
                 if (new_entity) {
                      std::string code = chain->SeqRes(i)->Field[0];
                      if (!SeqCodeUtil::is_a_standard_residue(code)) nstd_monomer = "yes";

                      _get_oneletter_code(code, code_can);

                      if (len + (int) code.size() > 80) { 
                           sequence += "\n";
                           len = 0;
                      }
                      sequence += code;
                      len += code.size();

                      if (len_p + (int) code_can.size() > 80) {
                           sequence_p += "\n";
                           len_p = 0;
                      }
                      sequence_p += code_can;
                      len_p += code_can.size();
                 }

                 // std::string asym_id = chain->ChainID();
                 // residue_entity_map.clear();
                 data.clear();
                 data.push_back(std::make_pair(chain->SeqRes(i)->Field[0], std::make_pair("", "")));
                 // residue_chainid_map.clear();
                 if (chain->SeqRes(i)->ResIndex >= 0) {
                      chain->GetResidueListByIndex(chain->SeqRes(i)->ResIndex, residue_list);
                      if (!residue_list.empty()) {
/*
                           asym_id = residue_list[0]->chnid();
                           if (chain->chain_type() == "ATOMS") {
                                residue_chainid_map.insert(std::make_pair(residue_list[0]->ResName(), residue_list[0]->pdb_chnid()));
                           }
*/
                           RCSB::Atom* atom = residue_list[0]->GetFirstAtom();
                           data[0].second.first = atom->pub_chnid();
                           data[0].second.second = atom->pub_resnum();
                      }
/*
                      if (!residue_list.empty() && !residue_list[0]->entity_id().empty())
                           residue_entity_map.insert(std::make_pair(chain->SeqRes(i)->Field[0], residue_list[0]->entity_id()));
*/
                      if (!residue_list.empty() && (chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP")) {
                           std::string code = residue_list[0]->ResName();
                           _get_oneletter_code(code, code_can);
                           if (len_c + (int) code.size() > 80) { 
                                sequence_c += "\n";
                                len_c = 0;
                           }
                           sequence_c += code;
                           len_c += code.size();

                           if (len_c_1 + (int) code_can.size() > 80) {
                                sequence_c_1 += "\n";
                                len_c_1 = 0;
                           }
                           sequence_c_1 += code_can;
                           len_c_1 += code_can.size();
                      }
                      if (residue_list.size() > 1) {
                           for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {  
                                // if (chain->chain_type() == "ATOMS") residue_chainid_map.insert(std::make_pair((*rpos)->ResName(), (*rpos)->pdb_chnid()));
                                bool found = false;
                                for (std::vector<std::pair<std::string, std::pair<std::string, std::string> > >::const_iterator
                                     pos = data.begin(); pos != data.end(); ++pos) {
                                     if (pos->first == (*rpos)->ResName()) {
                                          found = true;
                                          break;
                                     }
                                } 
                                if (!found) {
                                     RCSB::Atom* atom = (*rpos)->GetFirstAtom();
                                     data.push_back(std::make_pair((*rpos)->ResName(), std::make_pair(atom->pub_chnid(), atom->pub_resnum())));
                                     // if (!(*rpos)->entity_id().empty()) residue_entity_map.insert(std::make_pair((*rpos)->ResName(), (*rpos)->entity_id()));
                                }
                           }
                      }
                 }

                 residue_data.clear();
                 for (std::vector<std::pair<std::string, std::pair<std::string, std::string> > >::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
                      key_val_map.clear();
                      key_val_map.insert(std::make_pair("asym_id", chain->ChainID()));
/*
                      key_val_map.insert(std::make_pair("asym_id", asym_id));
                      std::map<std::string, std::string>::const_iterator mpos = residue_entity_map.find(pos->first);
                      if (mpos != residue_entity_map.end())
                           key_val_map.insert(std::make_pair("entity_id", mpos->second));
                      else */ key_val_map.insert(std::make_pair("entity_id", chain->entity_id()));
                      key_val_map.insert(std::make_pair("mon_id", pos->first));
                      key_val_map.insert(std::make_pair("pdb_seq_num", chain->SeqRes(i)->Field[4]));
                      // key_val_map.insert(std::make_pair("auth_seq_num", chain->SeqRes(i)->Field[6]));
                      key_val_map.insert(std::make_pair("auth_seq_num", pos->second.second));
                      if (data.size() > 1) {
                           key_val_map.insert(std::make_pair("pdb_mon_id", pos->first));
                           key_val_map.insert(std::make_pair("auth_mon_id", pos->first));
                      } else {
                           key_val_map.insert(std::make_pair("pdb_mon_id", chain->SeqRes(i)->Field[3]));
                           key_val_map.insert(std::make_pair("auth_mon_id", chain->SeqRes(i)->Field[5]));
                      }
                      if ( /* ( */ chain->chain_type() == "ATOMS" /* ) && (chain->SeqLen() > 1) */ ) {
                           key_val_map.insert(std::make_pair("num", chain->SeqRes(i)->Field[2]));
                           key_val_map.insert(std::make_pair("pdb_asym_id", chain->PDB_ChainID()));
                           key_val_map.insert(std::make_pair("auth_asym_id", pos->second.first));
                           if (data.size() > 1)
                                key_val_map.insert(std::make_pair("hetero", "y"));
                           else key_val_map.insert(std::make_pair("hetero", "n"));
                           _update_scheme_table(t4, key_val_map);

                           if (entity_id_set.find(chain->entity_id()) == entity_id_set.end()) {
                                key_val_map1.clear();
                                key_val_map1.insert(std::make_pair("entity_id", chain->entity_id()));
                                key_val_map1.insert(std::make_pair("comp_id", pos->first));
                                key_val_map1.insert(std::make_pair("num", chain->SeqRes(i)->Field[2]));
                                if (data.size() > 1)
                                     key_val_map1.insert(std::make_pair("hetero", "y"));
                                else key_val_map1.insert(std::make_pair("hetero", "n"));
                                _update_scheme_table(t6, key_val_map1);
                           }

                           _branch_residue_set.insert(pos->first);
                      } else {
                           key_val_map.insert(std::make_pair("ndb_seq_num", chain->SeqRes(i)->Field[2]));
                           if (!chain->SeqRes(i)->InsCode.empty())
                                key_val_map.insert(std::make_pair("pdb_ins_code", chain->SeqRes(i)->InsCode));
                           else key_val_map.insert(std::make_pair("pdb_ins_code", "."));
/*
                           mpos = residue_chainid_map.find(pos->first);
                           if (mpos != residue_chainid_map.end())
                                key_val_map.insert(std::make_pair("pdb_strand_id", mpos->second));
                           else */ key_val_map.insert(std::make_pair("pdb_strand_id", chain->PDB_ChainID()));
                           if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") {
                                key_val_map.insert(std::make_pair("seq_id", chain->SeqRes(i)->Field[2]));
                                if (data.size() > 1)
                                     key_val_map.insert(std::make_pair("hetero", "y"));
                                else key_val_map.insert(std::make_pair("hetero", "n"));
                                _update_scheme_table(t3, key_val_map);
                                tmp_data.clear();
                                tmp_data.push_back(pos->first);
                                tmp_data.push_back(chain->SeqRes(i)->Field[2]);
                                residue_data.push_back(tmp_data);
                           } else _update_scheme_table(t5, key_val_map);
                      }
                 }
                 entity_data.push_back(residue_data);
            }

            if (!sequence_c.empty()) {
                 tmp_data.clear();
                 tmp_data.push_back(chain->entity_id());
                 tmp_data.push_back(chain->PDB_ChainID());
                 tmp_data.push_back(sequence_c_1);
                 tmp_data.push_back(sequence_c);
                 coor_seq_data.push_back(tmp_data);
            }

            if (new_entity) {
                 std::map<int, Entity>::iterator epos = _entities.find(chain->int_entity_id());
                 if (epos != _entities.end())  {
                      if (chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP") {
                           // if (chain->has_sequence() || _DepUI_Flag)
                           epos->second.insertValue("has_sequence", "yes");
                           if (!chain->EvidenceCode().empty())
                                epos->second.insertValue("pdbx_sequence_evidence_code", chain->EvidenceCode());
                           else if (chain->has_sequence())
                                epos->second.insertValue("pdbx_sequence_evidence_code", "depositor provided");
                           else epos->second.insertValue("pdbx_sequence_evidence_code", "derived from coordinates");
                      }
                      epos->second.insertValue("nstd_monomer", nstd_monomer);
                      if (!sequence.empty())
                           epos->second.insertValue("pdbx_seq_one_letter_code", sequence);
                      if (!sequence_p.empty())
                           epos->second.insertValue("pdbx_seq_one_letter_code_can", sequence_p);
                 }
            }

            if ((chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP") && !entity_data.empty()) {
                 std::map<std::string, std::vector<std::vector<std::vector<string> > > >::iterator mpos = found_entities.find(chain->entity_id());
                 if (mpos == found_entities.end()) found_entities.insert(std::make_pair(chain->entity_id(), entity_data));
                 else if (mpos->second.size() == entity_data.size()) {
                      for (unsigned int i = 0; i < mpos->second.size(); ++i) {
                           for (unsigned int j = 0; j < entity_data[i].size(); ++j) {
                                bool found = false;
                                for (unsigned int k = 0; k < mpos->second[i].size(); ++k) {
                                     if (mpos->second[i][k][0] == entity_data[i][j][0]) {
                                          found = true;
                                          break;
                                     }
                                }
                                if (found) continue;
                                mpos->second[i].push_back(entity_data[i][j]);
                           }
                      }
                 }
            }

            entity_id_set.insert(chain->entity_id());

            chain = _molecules[0]->GetNextChain();
       }

       _write_entity_poly(block);

       if (!found_entities.empty()) {
            std::map<std::string, std::string> hetero_map;
            hetero_map.clear();

            ISTable *t2 = _newTablePtr("entity_poly_seq");
            int row2 = 0;
            for (std::map<unsigned int, std::string>::const_iterator opos = entity_order.begin(); opos != entity_order.end(); ++opos) {
                 std::map<std::string, std::vector<std::vector<std::vector<string> > > >::const_iterator mpos = found_entities.find(opos->second);
                 if (mpos == found_entities.end()) continue;
                 for (std::vector<std::vector<std::vector<string> > >::const_iterator pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                      for (std::vector<std::vector<string> >::const_iterator vpos = pos->begin(); vpos != pos->end(); ++vpos) {
                           t2->AddRow();
                           t2->UpdateCell(row2, "entity_id", opos->second);
                           t2->UpdateCell(row2, "num", (*vpos)[1]);
                           t2->UpdateCell(row2, "mon_id", (*vpos)[0]);
                           std::string index = opos->second + "_" + (*vpos)[1] + "_" + (*vpos)[0];
                           if (pos->size() > 1) {
                                t2->UpdateCell(row2, "hetero", "y");
                                hetero_map.insert(std::make_pair(index, "y"));
                           } else {
                                t2->UpdateCell(row2, "hetero", "n");
                                hetero_map.insert(std::make_pair(index, "n"));
                           }
                           row2++;
                      }
                 }
            }
            block.WriteTable(t2);

            std::string entity_id, seq_id, mon_id;
            for (unsigned int i = 0; i < t3->GetNumRows(); ++i) {
                 get_value_clean(entity_id, t3, i, "entity_id");
                 get_value_clean(seq_id, t3, i, "seq_id");
                 get_value_clean(mon_id, t3, i, "mon_id");
                 std::string index = entity_id + "_" + seq_id + "_" + mon_id;
                 std::map<std::string, std::string>::const_iterator mpos = hetero_map.find(index);
                 if (mpos != hetero_map.end()) t3->UpdateCell(i, "hetero", mpos->second); 
            }
       } else deleteTable(block, "entity_poly_seq");

       if (t3->GetNumRows() == 0) {
            delete t3;
            deleteTable(block, "pdbx_poly_seq_scheme");
       } else block.WriteTable(t3);

       if (t4->GetNumRows() == 0) {
            delete t4;
            deleteTable(block, "pdbx_branch_scheme");
       } else block.WriteTable(t4);

       if (t5->GetNumRows() == 0) {
            delete t5;
            deleteTable(block, "pdbx_nonpoly_scheme");
       } else block.WriteTable(t5);

       if (t6->GetNumRows() == 0) {
            delete t6;
            deleteTable(block, "pdbx_entity_branch_list");
       } else block.WriteTable(t6);

       if (!coor_seq_data.empty()) {
            tmp_data.clear();
            tmp_data.push_back("entity_id");
            tmp_data.push_back("auth_asym_id");
            tmp_data.push_back("one_letter_code");
            tmp_data.push_back("one_letter_code_mod");
            ISTable *t = add_new_table("pdbx_seq_map_depositor_info", tmp_data);
            for (unsigned int i = 0; i < coor_seq_data.size(); ++i) {
                 t->AddRow();
                 for (unsigned int j = 0; j < tmp_data.size(); ++j) {
                      t->UpdateCell(i, tmp_data[j], coor_seq_data[i][j]);
                 }
            }
            block.WriteTable(t);
       }
}

void AnnotationObj::_write_atom_type_and_chem_comp(Block &block)
{
       std::set<std::string> residue_names, atom_types, skip_residue_name_set;
       residue_names.clear();
       atom_types.clear();

       // Adjust _intra_molecular_connectivity_flag flag if "chem_comp_atom" is existed in input file or it's DepUI operation
       ISTable *t = getTablePtr(block, "chem_comp_atom");
       if (t && (t->GetNumRows() > 0)) _intra_molecular_connectivity_flag = true;
       if (_DepUI_Flag) _intra_molecular_connectivity_flag = false;

       skip_residue_name_set.clear();
       skip_residue_name_set.insert("N");
       skip_residue_name_set.insert("UNK");
       skip_residue_name_set.insert("UNL");
       skip_residue_name_set.insert("UNX");

       std::vector<RCSB::Residue*> residues;
       for (std::vector<RCSB::Molecule*>::iterator pos = _molecules.begin(); pos != _molecules.end(); ++pos) {
            RCSB::Chain* chain = (*pos)->GetFirstChain();
            while (chain) {
                 chain->GetFirstResidueList(residues);
                 while (!residues.empty()) {
                      for (std::vector<RCSB::Residue*>::iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                           residue_names.insert((*rpos)->ResName());
                           RCSB::Atom* atom = (*rpos)->GetFirstAtom();
                           while (atom) {
                                atom_types.insert(atom->atom_type());
                                atom = (*rpos)->GetNextAtom();
                           }
                      }
                      chain->GetNextResidueList(residues);
                 }
                 chain = (*pos)->GetNextChain();
            }
       }

       std::string cs;

       if (!atom_types.empty()) {
            std::map<std::string, std::map<std::string, std::string> > atom_info;
            atom_info.clear();
            std::vector<std::string> items;
            items.clear();
            std::map<std::string, std::string> tmp_map;

            t = getTablePtr(block, "atom_type");
            if (t) {
                 items = t->GetColumnNames();
                 for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                      tmp_map.clear();
                      for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
                           if (*pos == "symbol")
                                get_value_clean_upper(cs, t, i, *pos);
                           else get_value_clean(cs, t, i, *pos);
                           if (cs.empty()) continue;
                           tmp_map.insert(std::make_pair(*pos, cs));
                      }
                      if (tmp_map.empty()) continue;
                      std::map<std::string, std::string>::const_iterator mpos = tmp_map.find("symbol");
                      if (mpos != tmp_map.end()) atom_info.insert(std::make_pair(mpos->second, tmp_map));
                 }
            }

            for (std::set<std::string>::const_iterator spos = atom_types.begin(); spos != atom_types.end(); ++spos) {
                 cs = *spos;
                 if (cs.empty()) cs = ".";
                 std::map<std::string, std::map<std::string, std::string> >::const_iterator mpos = atom_info.find(cs);
                 if (mpos != atom_info.end()) continue;
                 tmp_map.clear();
                 tmp_map.insert(std::make_pair("symbol", cs));
                 atom_info.insert(std::make_pair(cs, tmp_map));
            }

            bool found = false;
            for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
                 if (*pos == "symbol") {
                      found = true;
                      break;
                 }
            }
            if (!found) items.push_back("symbol");

            t = add_new_table("atom_type", items);
            int i = 0;
            for (std::map<std::string, std::map<std::string, std::string> >::const_iterator mpos = atom_info.begin(); mpos != atom_info.end(); ++mpos) {
                 t->AddRow();
                 for (std::map<std::string, std::string>::const_iterator mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                      t->UpdateCell(i, mmpos->first, mmpos->second);
                 }
                 i++;
            }
            block.WriteTable(t);
       } else deleteTable(block, "atom_type");

       t = getTablePtr(block, "entity_poly_seq");
       if (t) {
            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean_upper(cs, t, i, "mon_id");
                 if (!cs.empty()) residue_names.insert(cs);
            }
       }

       t = getTablePtr(block, "struct_ref_seq_dif");
       if (t) {
            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean_upper(cs, t, i, "db_mon_id");
                 if (!cs.empty()) residue_names.insert(cs);
            }
       }

       std::vector<std::map<std::string, std::string> > branch_identifier_list, branch_synonyms_list;
       branch_identifier_list.clear();
       branch_synonyms_list.clear();

       std::map<std::string, std::set<std::string> > empty_check_mapping, pdbx_chem_comp_identifier_check_mapping;
       empty_check_mapping.clear();
       pdbx_chem_comp_identifier_check_mapping.clear();
       atom_types.clear();
       atom_types.insert("GMML");
       atom_types.insert("PDB-CARE");
       pdbx_chem_comp_identifier_check_mapping.insert(std::make_pair("program", atom_types));

       std::vector<std::vector<std::string> > chem_comp_atom_list, chem_comp_bond_list;
       chem_comp_atom_list.clear();
       chem_comp_bond_list.clear();

       std::vector<std::string> data_list;

       t = _newTablePtr("chem_comp");
       int i = 0;
       for (std::set<std::string>::const_iterator pos = residue_names.begin(); pos != residue_names.end(); ++pos) {
            t->AddRow();
            t->UpdateCell(i, "id", *pos);

            try {
                 ConnectFormat drug = _ccDic->find_drug(*pos);

                 cs = drug.getMetaData("type");
                 String::LowerCase(cs);
                 if (cs == "d-peptide linking")                    cs = "D-peptide linking";
                 else if (cs == "l-peptide linking")               cs = "L-peptide linking";
                 else if (cs == "d-peptide nh3 amino terminus")    cs = "D-peptide NH3 amino terminus";
                 else if (cs == "l-peptide nh3 amino terminus")    cs = "L-peptide NH3 amino terminus";
                 else if (cs == "d-peptide cooh carboxy terminus") cs = "D-peptide COOH carboxy terminus";
                 else if (cs == "l-peptide cooh carboxy terminus") cs = "L-peptide COOH carboxy terminus";
                 else if (cs == "dna linking")                     cs = "DNA linking";
                 else if (cs == "rna linking")                     cs = "RNA linking";
                 else if (cs.substr(0, 12) == "d-saccharide")      cs[0] = 'D';
                 else if (cs.substr(0, 12) == "l-saccharide")      cs[0] = 'L';
                 t->UpdateCell(i, "type", cs);

                 bool sugar_residue = false;
                 if ((cs.substr(0, 12) == "D-saccharide") || (cs.substr(0, 12) == "L-saccharide")) sugar_residue = true;

                 if (SeqCodeUtil::is_a_standard_residue(*pos) && *pos != "N" && *pos != "UNK") cs = "y";
                 else if (!drug.getMetaData("mon_nstd_parent_comp_id").empty()) cs = "n";
                 else cs = ".";
                 t->UpdateCell(i, "mon_nstd_flag", cs);

                 t->UpdateCell(i, "name", drug.chemical_name()); 
                 t->UpdateCell(i, "pdbx_synonyms", drug.synonym()); 

                 if (!drug.formula().empty() && drug.formula() != "X") {
                      cs = drug.formula();
                      std::string formal_charge = drug.getMetaData("pdbx_formal_charge");
                      if (!formal_charge.empty() && formal_charge != "0") cs += " " + formal_charge;
                      t->UpdateCell(i, "formula", cs); 
                 }

                 cs = drug.getMetaData("formula_weight");
                 if (cs.substr(0, 2) != "0.") t->UpdateCell(i, "formula_weight", cs);

                 // if (_branch_residue_set.find(*pos) != _branch_residue_set.end()) {
                 if (sugar_residue) {
                      ISTable* t1 = drug.getCategory("pdbx_chem_comp_identifier");
                      if (t1) _insert_chem_comp_feature_list(*pos, pdbx_chem_comp_identifier_check_mapping, t1, branch_identifier_list);
                      // t1 = drug.getCategory("pdbx_chem_comp_synonyms");
                      // if (t1) _insert_chem_comp_feature_list(*pos, empty_check_mapping, t1, branch_synonyms_list);
                 }

                 if (_intra_molecular_connectivity_flag && (skip_residue_name_set.find(*pos) == skip_residue_name_set.end())) {
                      const std::vector<AtomFormat>& chemCompAtoms = drug.atoms();
                      for (std::vector<AtomFormat>::const_iterator vpos = chemCompAtoms.begin(); vpos != chemCompAtoms.end(); ++vpos) {
                           data_list.clear();
                           data_list.push_back(*pos);
                           data_list.push_back(vpos->atomname());
                           data_list.push_back(vpos->atomtype());
                           data_list.push_back(vpos->get_value(5));
                           data_list.push_back(vpos->stereo_config());
                           data_list.push_back(String::IntToString(chem_comp_atom_list.size() + 1));
                           chem_comp_atom_list.push_back(data_list);
                      }

                      const std::vector<std::vector<std::string> >& chemCompBonds = drug.bonds();
                      for (std::vector<std::vector<std::string> >::const_iterator vpos = chemCompBonds.begin(); vpos != chemCompBonds.end(); ++vpos) {
                           data_list.clear();
                           data_list.push_back(*pos);
                           for (int j = 0; j < NUM_CHEM_COMP_BOND; ++j) {
                                data_list.push_back((*vpos)[j]);
                           }
                           data_list.push_back(String::IntToString(chem_comp_bond_list.size() + 1));
                           chem_comp_bond_list.push_back(data_list);
                      }
                 }
            } catch (const std::exception& exc) {
                 t->UpdateCell(i, "name", "UNKNOWN LIGAND");
                 t->UpdateCell(i, "type", "non-polymer");
                 t->UpdateCell(i, "mon_nstd_flag", "."); 
            }
            i++;
       }
       block.WriteTable(t);
       if (_intra_molecular_connectivity_flag) {
            if (!chem_comp_atom_list.empty()) {
                 data_list.clear();
                 data_list.push_back("comp_id");
                 data_list.push_back("atom_id");
                 data_list.push_back("type_symbol");
                 data_list.push_back("pdbx_aromatic_flag");
                 data_list.push_back("pdbx_stereo_config");
                 data_list.push_back("pdbx_ordinal");

                 t = add_new_table("chem_comp_atom", data_list, chem_comp_atom_list);
                 block.WriteTable(t);
            }
            if (!chem_comp_bond_list.empty()) {
                 data_list.clear();
                 data_list.push_back("comp_id");
                 data_list.push_back("atom_id_1");
                 data_list.push_back("atom_id_2");
                 data_list.push_back("value_order");
                 data_list.push_back("pdbx_aromatic_flag");
                 data_list.push_back("pdbx_stereo_config");
                 data_list.push_back("pdbx_ordinal");

                 t = add_new_table("chem_comp_bond", data_list, chem_comp_bond_list);
                 block.WriteTable(t);
            }
       }

       _write_chem_comp_feature_categories(block, "pdbx_chem_comp_identifier", branch_identifier_list);
       // _write_chem_comp_feature_categories(block, "pdbx_chem_comp_synonyms", branch_synonyms_list);
}

void AnnotationObj::_insert_chem_comp_feature_list(const std::string& comp_id, const std::map<std::string, std::set<std::string> >& checking_mapping,
                                                   ISTable* t, std::vector<std::map<std::string, std::string> >& branch_list)
{
       if (!t) return;

       std::string cs;
       std::map<std::string, std::string> feature_map;

       const std::vector<std::string>& itemNames = t->GetColumnNames();
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(cs, t, i, "comp_id");
            if (cs != comp_id) continue;

            if (!checking_mapping.empty()) {
                 bool found = false;
                 for (std::map<std::string, std::set<std::string> >::const_iterator mpos = checking_mapping.begin(); mpos != checking_mapping.end(); ++mpos) {
                      get_value_clean_upper(cs, t, i, mpos->first);
                      if (mpos->second.find(cs) != mpos->second.end()) found = true;
                 }
                 if (!found) continue;
            }

            feature_map.clear();
            for (std::vector<std::string>::const_iterator vpos = itemNames.begin(); vpos != itemNames.end(); ++vpos) {
                 get_value_clean(cs, t, i, *vpos);
                 if (!cs.empty()) feature_map.insert(std::make_pair(*vpos, cs));
            }
            if (feature_map.size() > 1) branch_list.push_back(feature_map);
       }
}

void AnnotationObj::_write_chem_comp_feature_categories(Block& block, const std::string& categoryName,
                                                        const std::vector<std::map<std::string, std::string> >& branch_list)
{
       deleteTable(block, categoryName);
       if (branch_list.empty()) return;

       ISTable *t = _newTablePtr(categoryName);
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       int row = 0;
       for (std::vector<std::map<std::string, std::string> >::const_iterator vpos = branch_list.begin(); vpos != branch_list.end(); ++vpos) {
            t->AddRow();
/*
            for (std::map<std::string, std::string>::const_iterator mpos = vpos->begin(); mpos != vpos->end(); ++mpos) {
                 t->UpdateCell(row, mpos->first, mpos->second);
            }
*/
            for (std::vector<std::string>::const_iterator ipos = itemNames.begin(); ipos != itemNames.end(); ++ipos) {
                 std::map<std::string, std::string>::const_iterator mpos = vpos->find(*ipos);
                 if (mpos != vpos->end()) t->UpdateCell(row, mpos->first, mpos->second);
            }
            row++;
       }
       block.WriteTable(t);
}

void AnnotationObj::_read_struct_asym(Block& block)
{
       ISTable *t = getTablePtr(block, "struct_asym");
       if (!t) return;

       std::string asym_id, blank_flag, details;
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(asym_id, t, i, "id");
            if (asym_id.empty()) continue;

            get_value_clean(blank_flag, t, i, "pdbx_blank_PDB_chainid_flag");
            get_value_clean(details, t, i, "details");
            if (blank_flag.empty() && details.empty()) continue;

            for (std::vector<RCSB::Molecule*>::iterator pos = _molecules.begin(); pos != _molecules.end(); ++pos) {
                 RCSB::Chain* chain = (*pos)->GetAsymChain(asym_id);
                 if (chain) {
                      if (!blank_flag.empty()) chain->set_PDB_ChainID_Flag(blank_flag);
                      if (!details.empty())    chain->set_details(details);
                 }
            }
       }
}

void AnnotationObj::_write_struct_asym(Block& block)
{
       if (_molecules.empty()) {
            deleteTable(block, "struct_asym");
            return;
       }

       ISTable *t = _newTablePtr("struct_asym");

       std::vector<RCSB::Residue*> residue_list;
       int row = 0;
       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
/*
            if (chain->chain_type() == "ATOMS") {
                 chain->GetFirstResidueList(residue_list);
                 while (!residue_list.empty()) {
                      t->AddRow();
                      t->UpdateCell(row, "id", residue_list[0]->chnid());
                      t->UpdateCell(row, "pdbx_PDB_id", chain->PDB_ChainID());
                      t->UpdateCell(row, "pdbx_alt_id", chain->PUB_ChainID());
                      if (!chain->PDB_ChainID_Flag().empty())
                           t->UpdateCell(row, "pdbx_blank_PDB_chainid_flag", chain->PDB_ChainID_Flag());
                      else t->UpdateCell(row, "pdbx_blank_PDB_chainid_flag", "N");
                      t->UpdateCell(row, "pdbx_type", chain->chain_type());
                      t->UpdateCell(row, "pdbx_order", String::IntToString(row + 1));
                      if (chain->isModified())
                           t->UpdateCell(row, "pdbx_modified", "Y");
                      else t->UpdateCell(row, "pdbx_modified", "N");
                      if (!residue_list[0]->entity_id().empty())
                           t->UpdateCell(row, "entity_id", residue_list[0]->entity_id());
                      else t->UpdateCell(row, "entity_id", chain->entity_id());
                      t->UpdateCell(row, "details", chain->details());
                      row++;

                      chain->GetNextResidueList(residue_list);
                 }
            } else {
*/
                 t->AddRow();
                 t->UpdateCell(row, "id", chain->ChainID());
                 t->UpdateCell(row, "pdbx_PDB_id", chain->PDB_ChainID());
                 t->UpdateCell(row, "pdbx_alt_id", chain->PUB_ChainID());
                 if (!chain->PDB_ChainID_Flag().empty())
                      t->UpdateCell(row, "pdbx_blank_PDB_chainid_flag", chain->PDB_ChainID_Flag());
                 else t->UpdateCell(row, "pdbx_blank_PDB_chainid_flag", "N");
                 t->UpdateCell(row, "pdbx_type", chain->chain_type());
                 t->UpdateCell(row, "pdbx_order", String::IntToString(row + 1));
                 if (chain->isModified())
                      t->UpdateCell(row, "pdbx_modified", "Y");
                 else t->UpdateCell(row, "pdbx_modified", "N");
                 t->UpdateCell(row, "entity_id", chain->entity_id());
                 t->UpdateCell(row, "details", chain->details());
                 row++;
/*
            }
*/

            chain = _molecules[0]->GetNextChain();
       }
       block.WriteTable(t);
}

void AnnotationObj::_write_entity_poly(Block& block)
{
       ISTable *t1 = _newTablePtr("entity_poly");

       const CIF_CATEGORY& category = CategoryMapping::find_category("entity_poly");

       std::string pdbx_strand_id;
       int row = 0;
       for (std::map<int, Entity>::const_iterator epos = _entities.begin(); epos != _entities.end(); ++epos) {
            if (epos->second.getValue("has_sequence") != "yes") continue;

            t1->AddRow();
            for (std::vector<std::string>::const_iterator pos = category.itemNames.begin(); pos != category.itemNames.end(); ++pos) {
                 if (*pos == "entity_id")
                      t1->UpdateCell(row, *pos, String::IntToString(epos->first));
                 else if (*pos == "nstd_linkage")
                      t1->UpdateCell(row, *pos, "no");
                 else if (*pos == "type") {
                      t1->UpdateCell(row, *pos, epos->second.getValue("poly_type"));
                 } else if (*pos == "pdbx_strand_id") {
                      const std::vector<std::string>& chainids = epos->second.PDB_chainID();
                      pdbx_strand_id.clear();
                      for (std::vector<std::string>::const_iterator cpos = chainids.begin(); cpos != chainids.end(); ++cpos) {
                           if (!pdbx_strand_id.empty()) pdbx_strand_id += ",";
                           pdbx_strand_id += *cpos;
                      }
                      t1->UpdateCell(row, *pos, pdbx_strand_id);
                 } else t1->UpdateCell(row, *pos, epos->second.getValue(*pos));
            }
            row++;
       }

       if (t1->GetNumRows() == 0) {
            delete t1;
            deleteTable(block, "entity_poly");
       } else block.WriteTable(t1);
}

void AnnotationObj::_update_scheme_table(ISTable *t, const std::map<std::string, std::string>& key_val_mapping)
{
       int row = t->GetNumRows();
       t->AddRow();
       for (std::map<std::string, std::string>::const_iterator mpos = key_val_mapping.begin(); mpos != key_val_mapping.end(); ++mpos) {
            t->UpdateCell(row, mpos->first, mpos->second);
       }
}

void AnnotationObj::_update_atom_sites(Block& block)
{
       ISTable *t = getTablePtr(block, "atom_sites");
       if (t) return;

       t = _newTablePtr("atom_sites");
       t->AddRow();
 
       t->UpdateCell(0, "entry_id", _StructureId);
       for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                 std::string item = "fract_transf_matrix[" + String::IntToString(i + 1) + "][" + String::IntToString(j + 1) + "]";
                 if (!_cell.is_artifical())
                      t->UpdateCell(0, item, FloatToString(_cell.o_to_f[i][j], 0, 6, false, false));
                 else if (i == j)
                      t->UpdateCell(0, item, FloatToString(1.0, 0, 6, false, false));
                 else t->UpdateCell(0, item, FloatToString(0.0, 0, 6, false, false));
          
            }
            std::string item = "fract_transf_vector[" + String::IntToString(i + 1) + "]";
            if (!_cell.is_artifical())
                 t->UpdateCell(0, item, FloatToString(_cell.o_to_f_v[i], 0, 5, false, false));
            else t->UpdateCell(0, item, FloatToString(0.0, 0, 5, false, false));
       }
       block.WriteTable(t);
}

void AnnotationObj::_update_pdb_format_compatible(Block &block)
{
       ISTable *t = _getTablePtr(block, "pdbx_database_status");
       if (!t) {
            t = _newTablePtr("pdbx_database_status");
            t->AddRow();
       }

       if (is_large_entry())
            t->UpdateCell(0, "pdb_format_compatible", "N");
       else t->UpdateCell(0, "pdb_format_compatible", "Y");
       block.WriteTable(t);
}

void AnnotationObj::_process_sg_entry(Block &block)
{
       ISTable *t = _getTablePtr(block, "pdbx_SG_project");
       if (!t) return;

       std::string cs;
       bool is_SG_entry = false;
       int is_PSI_entry = 0;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; i++) {
            get_value(cs, t, i, "project_name");
            if (!cs.empty()) is_SG_entry = true;

            get_value(cs, t, i, "full_name_of_center");
            if (cs.empty()) continue;

            std::map<std::string, std::string> sg_info = SgCenter::GetSgInfo(cs);
            if (sg_info.empty()) continue;

            is_SG_entry = true;

            std::map<std::string, std::string>::const_iterator mpos = sg_info.find("center_name");
            if (mpos != sg_info.end()) t->UpdateCell(i, "full_name_of_center", mpos->second);

            mpos = sg_info.find("center_initial");
            if (mpos != sg_info.end()) t->UpdateCell(i, "initial_of_center", mpos->second);

            mpos = sg_info.find("psi_flag");
            if (mpos != sg_info.end()) {
                 int psi_flag = atoi(mpos->second.c_str());
                 if (psi_flag > is_PSI_entry) is_PSI_entry = psi_flag;
                 if (psi_flag >= 2) {
                      get_value(cs, t, i, "project_name");
                      if (cs.empty()) {
                           mpos = sg_info.find("project_name");
                           if (mpos != sg_info.end()) t->UpdateCell(i, "project_name", mpos->second);
                      }
                 }
            }
       }
       block.WriteTable(t);

       if (!is_SG_entry) return;

       ISTable *t1 = _getTablePtr(block, "pdbx_database_status");
       if (t1) {
            t1->UpdateCell(0, "SG_entry", "Y");
            block.WriteTable(t1);
       }

       std::string full_name, FULL_NAME, initials, INITIALS, project, cs1;
       for (int k = 0; k < rowNo; k++) {
            get_value_clean(project, t, k, "project_name");
            get_value_clean(full_name, t, k, "full_name_of_center");
            get_value_clean(initials, t, k, "initial_of_center");
            String::UpperCase(full_name, FULL_NAME);
            String::UpperCase(initials, INITIALS);

            if (!full_name.empty()) {
                 std::string author_name = full_name;
                 if (!initials.empty()) author_name += " (" + initials + ")";

                 t1 = _getTablePtr(block, "audit_author");
                 if (!t1) t1 = _newTablePtr("audit_author");

                 bool found = false;
                 int rowNo1 = t1->GetNumRows();
                 for (int i = 0; i < rowNo1; ++i) {
                      get_value_clean_upper(cs, t1, i, "name");
                      if (cs.find(FULL_NAME) != std::string::npos && cs.find(INITIALS) != std::string::npos) {
                           t1->UpdateCell(i, "name", author_name);
                           found = true;
                           break;
                      }
                 }
                 if (!found) {
                      t1->AddRow();
                      t1->UpdateCell(rowNo1, "name", author_name);
                 }
                 block.WriteTable(t1);
            }

            t1 = _getTablePtr(block, "struct_keywords");
            if (!t1) {
                 t1 = _newTablePtr("struct_keywords");
                 t1->AddRow();
            }

            get_value_clean(cs, t1, 0, "text");
            String::UpperCase(cs, cs1);
            if (cs1.find("STRUCTURAL GENOMICS") == std::string::npos) {
                 if (!cs.empty()) cs += ", ";
                 cs += "Structural Genomics";
            }

            if (project == "PSI, Protein Structure Initiative") {
                 if (cs1.find("PSI") == std::string::npos) {
                      if (!cs.empty()) cs += ", ";
                      cs += "PSI-2";
                 }
                 if (cs1.find("PROTEIN STRUCTURE INITIATIVE") == std::string::npos) {
                      if (!cs.empty()) cs += ", ";
                      cs += "Protein Structure Initiative"; 
                 }
            } else if (project == "PSI:Biology") {
                 if (cs1.find("PSI") == std::string::npos) {
                      if (!cs.empty()) cs += ", ";
                      cs += "PSI-Biology";
                 }
            }
            if (!full_name.empty() && cs1.find(FULL_NAME) == std::string::npos) {
                 if (!cs.empty()) cs += ", ";
                 cs += full_name;
            }
            if (!initials.empty() && cs1.find(initials) == std::string::npos) {
                 if (!cs.empty()) cs += ", ";
                 cs += initials;
            }
            t1->UpdateCell(0, "text", cs);
            block.WriteTable(t1);
       }
 
       t = getTablePtr(block, "entity_poly"); 
       if (!t) return;

       std::vector<string> data;
       std::map<std::string, std::string> target_mapping;
       target_mapping.clear();

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value(cs, t, i, "pdbx_target_identifier");
            if (cs.empty()) continue;

            get_wordarray(data, cs, ",; ");
            for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
                 String::UpperCase(*pos, cs);
                 target_mapping.insert(std::make_pair(cs, *pos));
            }
       }
       if (target_mapping.empty()) return;

       t = _getTablePtr(block, "pdbx_database_related");
       if (t) {
            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean_upper(cs, t, i, "db_name");
                 if ((cs == "TARGETDB") || (cs == "TARGETTRACK")) {
                      get_value_clean_upper(cs, t, i, "db_id");
                      if (!cs.empty() && target_mapping.find(cs) != target_mapping.end())
                           target_mapping.erase(cs);
                 }
            }
       }
       if (target_mapping.empty()) return;

       if (!t) t = _newTablePtr("pdbx_database_related");
 
       for (std::map<std::string, std::string>::const_iterator mpos = target_mapping.begin(); mpos != target_mapping.end(); ++mpos) {
            int rowNo = t->GetNumRows();
            t->AddRow();
            t->UpdateCell(rowNo, "db_name", "TargetTrack");
            t->UpdateCell(rowNo, "db_id", mpos->second); 
       }
       block.WriteTable(t);
}

void AnnotationObj::_update_struct_keywords(Block& block)
{
       ISTable *t = _getTablePtr(block, "struct_keywords");
       if (!t) return;

       std::string header, keyword;
       get_value_clean(header, t, 0, "pdbx_keywords");
       get_value_clean(keyword, t, 0, "text");
       if (header.empty() && keyword.empty()) return;

       std::string stop_word_list = _rcsbroot + "/data/ascii/stop_word_list";
       remove_redundant_keywords(stop_word_list, keyword, header);
       t->UpdateCell(0, "text", keyword);
       block.WriteTable(t);
}

bool AnnotationObj::reorder_model()
{
       if ((_molecules.size() < 2) || !_CifObj) return true;

       std::string mol_id = "";
       Block &block = _CifObj->GetBlock(_firstBlockName);
       ISTable *t = _getTablePtr(block, "pdbx_nmr_representative");
       if (t) get_value_clean(mol_id,  t, 0, "conformer_id");
       if (mol_id.empty() || (mol_id == "1")) return true;

       return reorder_model(mol_id);
}

bool AnnotationObj::reorder_model(const std::string& mol_id)
{
       if (_molecules.size() < 2) return true;

       std::vector<RCSB::Molecule*> tmp_mol;
       tmp_mol.clear();

       for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            std::string Mol_ID = String::IntToString((*mpos)->Mol_ID());
            if (Mol_ID == mol_id) {
                 tmp_mol.push_back(*mpos);
                 break;
            }
       }

       if (tmp_mol.empty()) {
            std::string error = "Model ID '" + mol_id + "' does not exist in coordinate section.\n";
            _logIo->message(error.c_str());
            return false;
       }

       for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            std::string Mol_ID = String::IntToString((*mpos)->Mol_ID());
            if (Mol_ID != mol_id) tmp_mol.push_back(*mpos);
       }

       _molecules = tmp_mol;
       int id = 0;
       for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            id++;
            (*mpos)->set_Mol_ID(id);
       }

       if (_CifObj) {
            Block &block = _CifObj->GetBlock(_firstBlockName);
            ISTable *t = _getTablePtr(block, "pdbx_nmr_representative");
            if (t) {
                 t->UpdateCell(0, "conformer_id", "1");
                 block.WriteTable(t);
            }
       }

       return true;
}

bool AnnotationObj::is_large_entry()
{
       if (_molecules.empty()) return false;

       for (std::list<_SHEET>::const_iterator pos = _sheet.begin(); pos != _sheet.end(); ++pos) {
            if (pos->complicateFlag) return true;
       }

       std::vector<RCSB::Residue*> residue_list;
       for (std::vector<RCSB::Molecule*>::const_iterator pos = _molecules.begin(); pos != _molecules.end(); ++pos) {
            int chain_number_count = 0;
            int atom_numer_count = 0;

            RCSB::Chain* chain = (*pos)->GetFirstChain();
            while (chain) {
                 if ((chain->chain_type() == "ATOMP") || (chain->chain_type() == "ATOMN") || (chain->chain_type() == "ATOMS")) {
                      // Two letter chain ID
                      if (chain->PDB_ChainID().size() > 1) return true;
                      chain_number_count++;
                      // More than 62 polymer chains
                      if (chain_number_count > 62) return true;
                 }
                 atom_numer_count += chain->NumAtoms();
                 // More than 99999 atoms
                 if (atom_numer_count > 99999) return true;

                 chain->GetFirstResidueList(residue_list);
                 while (!residue_list.empty()) {
                      for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                           // Five letter residue name or B-factor > 1000
                           if (!(*rpos)->is_pdb_format_compatible()) return true;
                           // if ((*rpos)->ResName().size() > 3) return true;
                      }
                      chain->GetNextResidueList(residue_list);
                 }

                 chain = (*pos)->GetNextChain();
            }

            // Atom plus TER card numbers more than 99999 
            if ((atom_numer_count + chain_number_count) > 99999) return true;
       }

       return false;
}

void AnnotationObj::Output_public_cif_only()
{
       std::string pdb_id = getUpperCasePDBID(_CifObj);
       if (!pdb_id.empty()) _resetFirstBlockName = pdb_id;
  
       if (!_dictUtil.Read()) return;

       std::set<std::string> keeping_categories;
       keeping_categories.clear();
       keeping_categories.insert("pdbx_branch_scheme");
       keeping_categories.insert("pdbx_chem_comp_identifier");
       keeping_categories.insert("pdbx_chem_comp_synonyms");
       keeping_categories.insert("pdbx_entity_branch");
       keeping_categories.insert("pdbx_entity_branch_descriptor");
       keeping_categories.insert("pdbx_entity_branch_link");
       keeping_categories.insert("pdbx_entity_branch_list");

       Block& block = _CifObj->GetBlock(_firstBlockName);

       std::vector<std::string> TableNames;
       block.GetTableNames(TableNames);
       for (std::vector<std::string>::const_iterator tpos = TableNames.begin(); tpos != TableNames.end(); ++tpos) {
            if (!_dictUtil.isPublicCategory(*tpos) && (keeping_categories.find(*tpos) == keeping_categories.end())) {
                 if (_include_date_original_flag && (*tpos == "database_PDB_rev")) continue;
                 deleteTable(block, *tpos);
                 continue;
            }

            ISTable* t = getTablePtr(block, *tpos);
            std::vector<std::string> ColumnNames = t->GetColumnNames();
            for (std::vector<std::string>::const_iterator cpos = ColumnNames.begin(); cpos != ColumnNames.end(); ++cpos) {
                 if (_dictUtil.isPublicItem(*tpos, *cpos)) continue;
                 t->DeleteColumn(*cpos);
            }
            block.WriteTable(t);
       }
}
