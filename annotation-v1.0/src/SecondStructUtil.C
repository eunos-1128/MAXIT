/*
FILE:     SecondStructUtil.C
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
#include <string.h>
#include <stdlib.h>
#include "AnnotationObj.h"
#include "CategoryMapping.h"
#include "CompositeIndex.h"
#include "Promotif.h"
#include "utillib.h"

void AnnotationObj::CalculateSecondStruct(const std::string& supportfile, const bool& force_flag)
{
       if (_molecules.empty()) return;

       if (_skip_task_set.find("helix") != _skip_task_set.end() && _skip_task_set.find("sheet") != _skip_task_set.end()) return;

       if (supportfile.empty() && !force_flag && (!_helix.empty() || !_sheet.empty())) return;

       Promotif promotifUtil;

       std::string message;
       promotifUtil.setLog(_logIo);
       promotifUtil.setMolecule(_molecules[0]);
       promotifUtil.setSupportFile(supportfile);
       promotifUtil.calculateSecondaryStructure(message);
       if (_skip_task_set.find("helix") == _skip_task_set.end()) _helix = promotifUtil.getHelices();
       if (_skip_task_set.find("sheet") == _skip_task_set.end()) _sheet = promotifUtil.getSheets();

       if (supportfile.empty() || _sheet.empty()) return;

       if (_CifObj) {
            Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
            std::set<std::string> allow_set;
            allow_set.clear();
            allow_set.insert("sheet");
            _cif_update_add_value_to_pdbx_data_processing_status(_cifblock, allow_set);
       }
}

void AnnotationObj::Update_SecondStruct()
{
       if (!_CifObj) return;
       if (_molecules.empty()) return;

       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
       _update_struct_conf(_cifblock);
       _update_sheet_categories(_cifblock);
}

void AnnotationObj::_read_struct_conf(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "struct_conf");
       if (!t) return;

       _HELIX helix;
       helix.mol_index = _molecules[0]->index();
       helix.ID.clear();
       helix.initRes = -1;
       helix.endRes = -1;
       helix.helixClass = 0;
       helix.comment.clear();

       std::string cs, pdb_chnid, res_name, res_num, ins_code;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(pdb_chnid, t, i, "beg_auth_asym_id");
            get_value_clean(res_name, t, i, "beg_auth_comp_id");
            get_value_clean(res_num, t, i, "beg_auth_seq_id");
            get_value_clean(ins_code, t, i, "pdbx_beg_PDB_ins_code");
            RCSB::Residue* res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!res) continue;

            helix.initRes = res->index();

            get_value_clean(pdb_chnid, t, i, "end_auth_asym_id");
            get_value_clean(res_name, t, i, "end_auth_comp_id");
            get_value_clean(res_num, t, i, "end_auth_seq_id");
            get_value_clean(ins_code, t, i, "pdbx_end_PDB_ins_code");
            res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!res) continue;

            helix.endRes= res->index();

            get_value_clean(helix.ID, t, i, "pdbx_PDB_helix_id");
            get_value_clean(cs, t, i, "pdbx_PDB_helix_class");
            if (!cs.empty())
                 helix.helixClass = atoi(cs.c_str());
            else helix.helixClass = 0;
            get_value_clean(helix.comment, t, i, "details");
            get_value_clean_upper(cs, t, i, "conf_type_id");
            if (cs == "HELX_P")
                 _helix.push_back(helix);
            else if (cs == "TURN_P")
                 _turn.push_back(helix);
       }
}

void AnnotationObj::_update_struct_conf(Block& block)
{
       if (_helix.empty() && _turn.empty()) {
            deleteTable(block, "struct_conf_type");
            deleteTable(block, "struct_conf");
            return;
       }

       ISTable *t = _newTablePtr("struct_conf_type");
       int row = 0;
       if (!_helix.empty()) {
            t->AddRow();
            t->UpdateCell(row, "id", "HELX_P");
            row++;
       }
       if (!_turn.empty()) {
            t->AddRow();
            t->UpdateCell(row, "id", "TURN_P");
            row++;
       }
       block.WriteTable(t);

       t = _newTablePtr("struct_conf");
       if (!_helix.empty()) _update_struct_conf_table(t, "HELX_P", _helix);
       if (!_turn.empty()) _update_struct_conf_table(t, "TURN_P", _turn);
       block.WriteTable(t);
}

void AnnotationObj::_update_struct_conf_table(ISTable *t, const std::string& type, const std::list<_HELIX>& helices)
{
       int row = t->GetNumRows();
       int serial_no = 1;
       for (std::list<_HELIX>::const_iterator pos = helices.begin(); pos != helices.end(); ++pos) {
            bool is_removed = false;
            RCSB::Residue* initRes = _molecules[0]->find_residue(pos->initRes, is_removed);
            if (!initRes) continue;
            RCSB::Residue* endRes = _molecules[0]->find_residue(pos->endRes, is_removed);
            if (!endRes) continue;

            t->AddRow();
            t->UpdateCell(row, "conf_type_id", type);
            t->UpdateCell(row, "id", type + String::IntToString(serial_no));
            t->UpdateCell(row, "pdbx_PDB_helix_id", pos->ID);
            t->UpdateCell(row, "beg_label_comp_id",     initRes->ResName());
            t->UpdateCell(row, "beg_label_asym_id",     initRes->chnid());
            t->UpdateCell(row, "beg_label_seq_id",      initRes->res_no());
            t->UpdateCell(row, "pdbx_beg_PDB_ins_code", initRes->ins_code());
            t->UpdateCell(row, "end_label_comp_id",     endRes->ResName());
            t->UpdateCell(row, "end_label_asym_id",     endRes->chnid());
            t->UpdateCell(row, "end_label_seq_id",      endRes->res_no());
            t->UpdateCell(row, "pdbx_end_PDB_ins_code", endRes->ins_code());
            t->UpdateCell(row, "beg_auth_comp_id",      initRes->ResName());
            t->UpdateCell(row, "beg_auth_asym_id",      initRes->pdb_chnid());
            t->UpdateCell(row, "beg_auth_seq_id",       initRes->pdb_res_no());
            t->UpdateCell(row, "end_auth_comp_id",      endRes->ResName());
            t->UpdateCell(row, "end_auth_asym_id",      endRes->pdb_chnid());
            t->UpdateCell(row, "end_auth_seq_id",       endRes->pdb_res_no());
            if (pos->helixClass) {
                 t->UpdateCell(row, "pdbx_PDB_helix_class", String::IntToString(pos->helixClass));
            }
            t->UpdateCell(row, "details", pos->comment);
            if (type == "HELX_P") {
                 int length = endRes->position() - initRes->position() + 1;
                 t->UpdateCell(row, "pdbx_PDB_helix_length", String::IntToString(length));
            }
            row++;
            serial_no++;
       }
}

void AnnotationObj::_read_struct_sheet(Block& block)
{
       if (_molecules.empty()) return;

       std::string cs, cs1;
       bool pdb_format_compatible = true;
       ISTable *t = getTablePtr(block, "pdbx_database_status");
       if (t) {
            get_value_clean_upper(cs, t, 0, "pdb_format_compatible");
            if (cs == "N") pdb_format_compatible = false;
       }

       t = getTablePtr(block, "struct_sheet");
       if (!t) return;

       std::vector<std::string> sheet_id_list;
       sheet_id_list.clear();
       sheet_id_list.reserve(t->GetNumRows());

       // key: sheetID
       std::map<std::string, _SHEET> sheet_mapping;
       sheet_mapping.clear();

       // first key: sheetID
       // second key: strand_id
       std::map<std::string, std::map<std::string, _SHEET_STRAND> > sheet_strand_mapping;
       sheet_strand_mapping.clear();

       // key: strand_id
       std::map<std::string, _SHEET_STRAND> tmp_strand_mapping;
       tmp_strand_mapping.clear();

       _SHEET sheet;

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(cs, t, i, "id");
            if (cs.empty()) continue;
            get_value_clean(cs1, t, i, "number_strands");
            if (cs1.empty()) continue;

            sheet.mol_index = _molecules[0]->index();
            sheet.sheetID = cs;
            sheet.numStrands = atoi(cs1.c_str());
            sheet.complicateFlag = false;
            sheet._strands.clear();
            sheet._strand_orders.clear();
            sheet_id_list.push_back(sheet.sheetID);
            sheet_mapping.insert(std::make_pair(sheet.sheetID, sheet));
            sheet_strand_mapping.insert(std::make_pair(sheet.sheetID, tmp_strand_mapping));
       }
       if (sheet_mapping.empty()) return;

       t = getTablePtr(block, "struct_sheet_range");
       if (!t) return;

       _SHEET_STRAND strand;
       strand.begin_res_index = -1;
       strand.end_res_index = -1;
       strand.sense = 0;
       strand.curr_hbond_res_index = -1;
       strand.curr_hbond_atom_name.clear();
       strand.prev_hbond_res_index = -1;
       strand.prev_hbond_atom_name.clear();

       bool has_residue_nomenclature_error = false;
       std::string pdb_chnid, res_name, res_num, ins_code;
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(pdb_chnid, t, i, "beg_auth_asym_id");
            get_value_clean(res_name, t, i, "beg_auth_comp_id");
            get_value_clean(res_num, t, i, "beg_auth_seq_id");
            get_value_clean(ins_code, t, i, "pdbx_beg_PDB_ins_code");
            RCSB::Residue* begin_res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!begin_res) {
                 has_residue_nomenclature_error = true;
                 continue;
            }

            get_value_clean(pdb_chnid, t, i, "end_auth_asym_id");
            get_value_clean(res_name, t, i, "end_auth_comp_id");
            get_value_clean(res_num, t, i, "end_auth_seq_id");
            get_value_clean(ins_code, t, i, "pdbx_end_PDB_ins_code");
            RCSB::Residue* end_res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!end_res) {
                 has_residue_nomenclature_error = true;
                 continue;
            }

            strand.begin_res_index = begin_res->index();
            strand.end_res_index = end_res->index();

            get_value_clean(cs, t, i, "sheet_id");
            get_value_clean(cs1, t, i, "id");

            std::map<std::string, std::map<std::string, _SHEET_STRAND> >::iterator mpos = sheet_strand_mapping.find(cs);
            if (mpos == sheet_strand_mapping.end()) continue; //
            if (mpos->second.find(cs1) != mpos->second.end()) continue; //

            strand.strand_id = mpos->second.size() + 1;
            mpos->second.insert(std::make_pair(cs1, strand));
       }

       std::map<int, _SHEET_STRAND> strand_mapping;
       for (std::map<std::string, _SHEET>::iterator mpos = sheet_mapping.begin(); mpos != sheet_mapping.end(); ++mpos) {
            std::map<std::string, std::map<std::string, _SHEET_STRAND> >::const_iterator ssmpos = sheet_strand_mapping.find(mpos->first);
            if (ssmpos == sheet_strand_mapping.end()) return; //
            if (mpos->second.numStrands != ssmpos->second.size()) {
                 if (!has_residue_nomenclature_error) mpos->second.numStrands = ssmpos->second.size();
                 else continue;
            }

            strand_mapping.clear();
            for (std::map<std::string, _SHEET_STRAND>::const_iterator smpos = ssmpos->second.begin(); smpos != ssmpos->second.end(); ++smpos) {
                 strand_mapping.insert(std::make_pair(smpos->second.strand_id, smpos->second));
            }
            mpos->second._strands.clear();
            mpos->second._strands.reserve(strand_mapping.size());
            for (std::map<int, _SHEET_STRAND>::const_iterator smpos = strand_mapping.begin(); smpos != strand_mapping.end(); ++smpos) {
                 mpos->second._strands.push_back(smpos->second);
            }
       }

       t = getTablePtr(block, "struct_sheet_order");
       if (t) {
            // key: sheet_id
            // value.first.key: "range_id_1"_"range_id_2"
            // value.first.value: index to value.second
            // value.second: strand_orders
            std::map<std::string, std::pair<std::map<std::string, unsigned int>, std::vector<_SHEET_TOPOLOGY> > > sheet_order_map;
            sheet_order_map.clear();

            std::map<std::string, unsigned int> t_map;
            std::vector<_SHEET_TOPOLOGY> t_vec;

            _SHEET_TOPOLOGY s_order;
            s_order.sense_number = 0;
            s_order.sense_string.clear();
            s_order.range_id_1.clear();
            s_order.range_id_1_res_index = -1;
            s_order.range_id_1_atom_name.clear();
            s_order.range_id_2.clear();
            s_order.range_id_2_res_index = -1;
            s_order.range_id_2_atom_name.clear();
 
            std::string sheet_id, sense;
            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean(sheet_id, t, i, "sheet_id");
                 get_value_clean_lower(sense, t, i, "sense");
                 get_value_clean(cs, t, i, "range_id_1");
                 get_value_clean(cs1, t, i, "range_id_2");
                 if (sheet_id.empty() || sense.empty() || cs.empty() || cs1.empty() || (cs == cs1)) continue;

                 s_order.sense_string = sense;
                 if (sense == "parallel")
                      s_order.sense_number = 1;
                 else if (sense == "anti-parallel")
                      s_order.sense_number = -1;
                 s_order.range_id_1 = cs;
                 s_order.range_id_2 = cs1;

                 std::map<std::string, std::pair<std::map<std::string, unsigned int>, std::vector<_SHEET_TOPOLOGY> > >::iterator
                     mpos = sheet_order_map.find(sheet_id);
                 if (mpos != sheet_order_map.end()) {
                      mpos->second.first.insert(std::make_pair(cs + "_" + cs1, mpos->second.second.size()));
                      mpos->second.second.push_back(s_order);
                 } else {
                      t_map.clear();
                      t_vec.clear();
                      t_map.insert(std::make_pair(cs + "_" + cs1, t_vec.size()));
                      t_vec.push_back(s_order);
                      sheet_order_map.insert(std::make_pair(sheet_id, std::make_pair(t_map, t_vec)));
                 }
            }

            t = getTablePtr(block, "pdbx_struct_sheet_hbond");
            if (t) {
                 std::string range_id_1_atom, range_id_2_atom;
                 for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                      get_value_clean(sheet_id, t, i, "sheet_id");
                      get_value_clean(cs, t, i, "range_id_1");
                      get_value_clean(cs1, t, i, "range_id_2");
                      if (sheet_id.empty() || cs.empty() || cs1.empty() || (cs == cs1)) continue;

                      std::map<std::string, std::pair<std::map<std::string, unsigned int>, std::vector<_SHEET_TOPOLOGY> > >::iterator
                          mpos = sheet_order_map.find(sheet_id);
                      if (mpos == sheet_order_map.end()) continue; //

                      std::map<std::string, unsigned int>::const_iterator impos = mpos->second.first.find(cs + "_" + cs1);
                      if (impos == mpos->second.first.end()) continue; //

                      get_value_clean(pdb_chnid, t, i, "range_1_auth_asym_id");
                      get_value_clean(res_name, t, i, "range_1_auth_comp_id");
                      get_value_clean(res_num, t, i, "range_1_auth_seq_id");
                      get_value_clean(ins_code, t, i, "range_1_PDB_ins_code");
                      get_value_clean(range_id_1_atom, t, i, "range_1_auth_atom_id");
                      RCSB::Residue* range_1_res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
                      if (!range_1_res) continue; //

                      get_value_clean(pdb_chnid, t, i, "range_2_auth_asym_id");
                      get_value_clean(res_name, t, i, "range_2_auth_comp_id");
                      get_value_clean(res_num, t, i, "range_2_auth_seq_id");
                      get_value_clean(ins_code, t, i, "range_2_PDB_ins_code");
                      get_value_clean(range_id_2_atom, t, i, "range_2_auth_atom_id");
                      RCSB::Residue* range_2_res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
                      if (!range_2_res) continue; //

                      mpos->second.second[impos->second].range_id_1_res_index = range_1_res->index();
                      mpos->second.second[impos->second].range_id_1_atom_name = range_id_1_atom;
                      mpos->second.second[impos->second].range_id_2_res_index = range_2_res->index();
                      mpos->second.second[impos->second].range_id_2_atom_name = range_id_2_atom;
                 }
            }

            std::set<int> strand_order_range_ids;
            for (std::map<std::string, std::pair<std::map<std::string, unsigned int>, std::vector<_SHEET_TOPOLOGY> > >::const_iterator
                 sompos = sheet_order_map.begin(); sompos != sheet_order_map.end(); ++sompos) {
                 std::map<std::string, _SHEET>::iterator mpos = sheet_mapping.find(sompos->first);
                 if ((mpos == sheet_mapping.end()) || mpos->second._strands.empty()) continue; //
                 std::map<std::string, std::map<std::string, _SHEET_STRAND> >::const_iterator ssmpos = sheet_strand_mapping.find(mpos->first);
                 if (ssmpos == sheet_strand_mapping.end()) return; //

                 std::vector<_SHEET_TOPOLOGY> strand_orders = sompos->second.second;
                 strand_order_range_ids.clear();
                 for (std::vector<_SHEET_TOPOLOGY>::iterator sopos = strand_orders.begin(); sopos != strand_orders.end(); ++sopos) {
                      std::map<std::string, _SHEET_STRAND>::const_iterator smpos = ssmpos->second.find(sopos->range_id_1);
                      if (smpos == ssmpos->second.end()) continue;
                      // check range
                      if ((sopos->range_id_1_res_index >= 0) && ((sopos->range_id_1_res_index < smpos->second.begin_res_index) ||
                          (sopos->range_id_1_res_index > smpos->second.end_res_index))) continue;

                      sopos->range_id_1 = String::IntToString(smpos->second.strand_id);
                      strand_order_range_ids.insert(smpos->second.strand_id);

                      smpos = ssmpos->second.find(sopos->range_id_2);
                      if (smpos == ssmpos->second.end()) continue;
                      // check range
                      if ((sopos->range_id_2_res_index >= 0) && ((sopos->range_id_2_res_index < smpos->second.begin_res_index) ||
                          (sopos->range_id_2_res_index > smpos->second.end_res_index))) continue;

                      sopos->range_id_2 = String::IntToString(smpos->second.strand_id);
                      strand_order_range_ids.insert(smpos->second.strand_id);
                 }
                 bool unusual_strand_order_flag = false;
                 if (!strand_orders.empty() && (strand_order_range_ids.size() != mpos->second._strands.size())) unusual_strand_order_flag = true;

                 bool is_linear_topology = true;
                 if (!strand_orders.empty() && !unusual_strand_order_flag) {
                      std::string first_strand_id = String::IntToString(mpos->second._strands[0].strand_id);
                      std::string last_strand_id = String::IntToString(mpos->second._strands[mpos->second._strands.size()-1].strand_id);
                      if ((strand_orders[0].range_id_1 != first_strand_id) || (strand_orders[strand_orders.size()-1].range_id_2 != last_strand_id))
                          is_linear_topology = false;
                      for (unsigned int i = 1; i < strand_orders.size(); ++i) {
                           if (strand_orders[i - 1].range_id_2 != strand_orders[i].range_id_1) is_linear_topology = false;
                      }
                 }
                 mpos->second._strand_orders = strand_orders;
                 if (!is_linear_topology && !pdb_format_compatible) mpos->second.complicateFlag = true;
                 else {
                      for (std::vector<_SHEET_TOPOLOGY>::const_iterator sopos = strand_orders.begin(); sopos != strand_orders.end(); ++sopos) {
                           int range_id_2 = atoi(sopos->range_id_2.c_str()) - 1;
                           mpos->second._strands[range_id_2].sense = sopos->sense_number;
                           mpos->second._strands[range_id_2].prev_hbond_res_index = sopos->range_id_1_res_index;
                           mpos->second._strands[range_id_2].prev_hbond_atom_name = sopos->range_id_1_atom_name;
                           mpos->second._strands[range_id_2].curr_hbond_res_index = sopos->range_id_2_res_index;
                           mpos->second._strands[range_id_2].curr_hbond_atom_name = sopos->range_id_2_atom_name;
                      }
                 }
            }
       }

       for (std::vector<std::string>::const_iterator vpos = sheet_id_list.begin(); vpos != sheet_id_list.end(); ++vpos) {
            std::map<std::string, _SHEET>::const_iterator mpos = sheet_mapping.find(*vpos);
            if ((mpos == sheet_mapping.end()) || mpos->second._strands.empty()) continue; //
            _sheet.push_back(mpos->second);
       }
}

void AnnotationObj::_update_sheet_categories(Block& block)
{
       if (_sheet.empty()) {
            deleteTable(block, "struct_sheet");
            deleteTable(block, "struct_sheet_order");
            deleteTable(block, "struct_sheet_range");
            deleteTable(block, "pdbx_struct_sheet_hbond");
            return;
       }

       ISTable *t1 = _newTablePtr("struct_sheet");
       ISTable *t2 = _newTablePtr("struct_sheet_order");
       ISTable *t3 = _newTablePtr("struct_sheet_range");
       ISTable *t4 = _newTablePtr("pdbx_struct_sheet_hbond");

       std::vector<std::map<std::string, std::string> > range_data, order_data, hbond_data;

       for (std::list<_SHEET>::const_iterator pos = _sheet.begin(); pos != _sheet.end(); ++pos) {
            if (!_get_range_data(pos->sheetID, pos->_strands, range_data)) continue;
            if (!_get_order_and_hbond_data(pos->sheetID, pos->_strand_orders, order_data, hbond_data)) continue;

            _insert_sheet(t1, pos->sheetID, String::IntToString(pos->numStrands));
            update_table(t2, order_data);
            update_table(t3, range_data);
            update_table(t4, hbond_data);
       }

       if (t1->GetNumRows() > 0) block.WriteTable(t1);
       else { delete t1; deleteTable(block, "struct_sheet"); }
       if (t2->GetNumRows() > 0) block.WriteTable(t2);
       else { delete t2; deleteTable(block, "struct_sheet_order"); }
       if (t3->GetNumRows() > 0) block.WriteTable(t3);
       else { delete t3; deleteTable(block, "struct_sheet_range"); }
       if (t4->GetNumRows() > 0) block.WriteTable(t4);
       else { delete t4; deleteTable(block, "pdbx_struct_sheet_hbond"); }
}

void  AnnotationObj::_insert_sheet(ISTable *t, const std::string& sheetID, const std::string& num_strands)
{
       int row = t->GetNumRows();
       t->AddRow();
       t->UpdateCell(row, "id", sheetID);
       t->UpdateCell(row, "number_strands", num_strands);
}

bool AnnotationObj::_get_range_data(const std::string& sheetID, const std::vector<_SHEET_STRAND>& strands,
                                    std::vector<std::map<std::string, std::string> >& range_data)
{
       range_data.clear();
       if (strands.empty()) return false;

       range_data.reserve(strands.size());

       bool is_removed = false;
       std::map<std::string, std::string> data_map;
       for (std::vector<_SHEET_STRAND>::const_iterator vpos = strands.begin(); vpos != strands.end(); ++vpos) {
            RCSB::Residue* initRes = _molecules[0]->find_residue(vpos->begin_res_index, is_removed);
            if (!initRes) return false;
            RCSB::Residue* endRes = _molecules[0]->find_residue(vpos->end_res_index, is_removed);
            if (!endRes) return false;

            data_map.clear();
            data_map.insert(std::make_pair("sheet_id",              sheetID));
            data_map.insert(std::make_pair("id",                    String::IntToString(vpos->strand_id)));
            data_map.insert(std::make_pair("beg_label_comp_id",     initRes->ResName()));
            data_map.insert(std::make_pair("beg_label_asym_id",     initRes->chnid()));
            data_map.insert(std::make_pair("beg_label_seq_id",      initRes->res_no()));
            data_map.insert(std::make_pair("pdbx_beg_PDB_ins_code", initRes->ins_code()));
            data_map.insert(std::make_pair("beg_auth_comp_id",      initRes->ResName()));
            data_map.insert(std::make_pair("beg_auth_asym_id",      initRes->pdb_chnid()));
            data_map.insert(std::make_pair("beg_auth_seq_id",       initRes->pdb_res_no()));
            data_map.insert(std::make_pair("end_label_comp_id",     endRes->ResName()));
            data_map.insert(std::make_pair("end_label_asym_id",     endRes->chnid()));
            data_map.insert(std::make_pair("end_label_seq_id",      endRes->res_no()));
            data_map.insert(std::make_pair("pdbx_end_PDB_ins_code", endRes->ins_code()));
            data_map.insert(std::make_pair("end_auth_comp_id",      endRes->ResName()));
            data_map.insert(std::make_pair("end_auth_asym_id",      endRes->pdb_chnid()));
            data_map.insert(std::make_pair("end_auth_seq_id",       endRes->pdb_res_no()));
            range_data.push_back(data_map);
       }

       if (range_data.empty()) return false;
       return true;
}

bool AnnotationObj::_get_order_and_hbond_data(const std::string& sheetID, const std::vector<_SHEET_TOPOLOGY>& strand_orders, std::vector<std::map<std::string,
                                              std::string> >& order_data, std::vector<std::map<std::string, std::string> >& hbond_data)
{
       order_data.clear();
       hbond_data.clear();

       if (strand_orders.empty()) return true;

       order_data.reserve(strand_orders.size());
       hbond_data.reserve(strand_orders.size());

       bool is_removed = false;
       std::map<std::string, std::string> data_map;
       for (std::vector<_SHEET_TOPOLOGY>::const_iterator vpos = strand_orders.begin(); vpos != strand_orders.end(); ++vpos) {
            if (!vpos->sense_string.empty()) {
                 data_map.clear();
                 data_map.insert(std::make_pair("sheet_id",   sheetID));
                 data_map.insert(std::make_pair("range_id_1", vpos->range_id_1));
                 data_map.insert(std::make_pair("range_id_2", vpos->range_id_2));
                 data_map.insert(std::make_pair("sense",      vpos->sense_string));   
                 order_data.push_back(data_map);
            }

            if (vpos->range_id_1_res_index < 0) continue;
            RCSB::Residue* range_id_1_res = _molecules[0]->find_residue(vpos->range_id_1_res_index, is_removed);
            if (!range_id_1_res) return false;
            if (vpos->range_id_2_res_index < 0) continue;
            RCSB::Residue* range_id_2_res = _molecules[0]->find_residue(vpos->range_id_2_res_index, is_removed);
            if (!range_id_2_res) return false;

            data_map.clear();
            data_map.insert(std::make_pair("sheet_id",              sheetID));
            data_map.insert(std::make_pair("range_id_1",            vpos->range_id_1));
            data_map.insert(std::make_pair("range_1_label_atom_id", vpos->range_id_1_atom_name));
            data_map.insert(std::make_pair("range_1_label_comp_id", range_id_1_res->ResName()));
            data_map.insert(std::make_pair("range_1_label_asym_id", range_id_1_res->chnid()));
            data_map.insert(std::make_pair("range_1_label_seq_id",  range_id_1_res->res_no()));
            data_map.insert(std::make_pair("range_1_PDB_ins_code",  range_id_1_res->ins_code()));
            data_map.insert(std::make_pair("range_1_auth_atom_id",  vpos->range_id_1_atom_name));
            data_map.insert(std::make_pair("range_1_auth_comp_id",  range_id_1_res->ResName()));
            data_map.insert(std::make_pair("range_1_auth_asym_id",  range_id_1_res->pdb_chnid()));
            data_map.insert(std::make_pair("range_1_auth_seq_id",   range_id_1_res->pdb_res_no()));
            data_map.insert(std::make_pair("range_id_2",            vpos->range_id_2));
            data_map.insert(std::make_pair("range_2_label_atom_id", vpos->range_id_2_atom_name));
            data_map.insert(std::make_pair("range_2_label_comp_id", range_id_2_res->ResName()));
            data_map.insert(std::make_pair("range_2_label_asym_id", range_id_2_res->chnid()));
            data_map.insert(std::make_pair("range_2_label_seq_id",  range_id_2_res->res_no()));
            data_map.insert(std::make_pair("range_2_PDB_ins_code",  range_id_2_res->ins_code()));
            data_map.insert(std::make_pair("range_2_auth_atom_id",  vpos->range_id_2_atom_name));
            data_map.insert(std::make_pair("range_2_auth_comp_id",  range_id_2_res->ResName()));
            data_map.insert(std::make_pair("range_2_auth_asym_id",  range_id_2_res->pdb_chnid()));
            data_map.insert(std::make_pair("range_2_auth_seq_id",   range_id_2_res->pdb_res_no()));
            hbond_data.push_back(data_map);
       }

       return true;
}
