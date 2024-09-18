/*
FILE:     FileObj_Cif_Util.C
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
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdexcept>

#include "algorithm-util.h"
#include "algorithm-util.C"
#include "BondUtil.h"
#include "CategoryMapping.h"
#include "CifCoordRead.h"
#include "CompositeIndex.h"
#include "FileObj.h"
#include "SeqCodeUtil.h"
#include "SpaceGroup.h"
#include "utillib.h"

#define N_SWITCH_CATEGORY     4
#define N_NAME_CATEGORY       3

bool FileObj::parse_mmcif_file()
{
       if (_CifObj) { delete _CifObj; }
       std::string msg;
       _CifObj = get_fobj(msg, _input_filename, _caseSensitive);
       if (!_CifObj) {
            std::string error = "Read file '" + _input_filename + "' failed";
            if (!msg.empty()) {
                 error += ": ";
                 if (msg.substr(0, 8) == "ERROR - ")
                      error += msg.substr(8);
                 else error += msg;
            } else error += ".";
            throw std::out_of_range(error);
            return false;
       }

       _firstBlockName = _CifObj->GetFirstBlockName();
       Block& block = _CifObj->GetBlock(_firstBlockName);

       _StructureId = _firstBlockName;
       _UppercasePdbId = getUpperCasePDBID(block);
       _LowercasePdbId = getLowerCasePDBID(block);

       _read_branch_polymer_info();

       return true;
}

void FileObj::clear_branch_polymer_info()
{
       _Branch_Seq_Scheme_Mapping.clear();

       if (!_CifObj) return;

       Block& block = _CifObj->GetBlock(_firstBlockName);

       std::set<std::string> remove_categories;
       remove_categories.clear();
       remove_categories.insert("pdbx_branch_scheme");
       remove_categories.insert("pdbx_chem_comp_identifier");
       remove_categories.insert("pdbx_chem_comp_synonyms");
       remove_categories.insert("pdbx_entity_branch");
       remove_categories.insert("pdbx_entity_branch_descriptor");
       remove_categories.insert("pdbx_entity_branch_link");
       remove_categories.insert("pdbx_entity_branch_list");

       for (std::set<std::string>::const_iterator spos = remove_categories.begin(); spos != remove_categories.end(); ++spos) {
            deleteTable(block, *spos);
       }
}

bool FileObj::read_mmcif_file(const bool& coord_and_seq_only, const bool& check_format_flag)
{
       if (!read_mmcif_metadata(coord_and_seq_only)) return false;
       read_mmcif_coordinate(check_format_flag);

       return true;
}

void FileObj::read_mmcif_file_without_parsing(const bool& coord_and_seq_only)
{
       _read_mmcif_metadata(coord_and_seq_only);
       read_mmcif_coordinate();
}

bool FileObj::read_mmcif_metadata(const bool& coord_and_seq_only)
{
       if (!parse_mmcif_file()) return false;

       _read_mmcif_metadata(coord_and_seq_only);

       return true;
}

void FileObj::read_mmcif_coordinate(const bool& check_format_flag)
{
       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
       ISTable *t = getTablePtr(_cifblock, "atom_site");
       if (t) _read_atom_site_record_number = t->GetNumRows();

       CifCoordRead reader;
       reader.setLog(_logIo);
       if (check_format_flag) {
            reader.setMessage(&_error_messages);
            reader.setExperimentType(_experiment_type);
            reader.setFormatChecking();
       }
       reader.setDictObj(&_dictUtil);
       reader.setRecord(&_pdb_records);
       reader.setMolecule(&_molecules);
       reader.setCifObj(_CifObj);
       reader.Read();
       _extra_atom_site_item_list = reader.getExtraItems();
}

bool FileObj::is_pdb_format_compatible()
{
       if (_CifObj) {
            Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
            ISTable *t = getTablePtr(_cifblock, "pdbx_database_status");
            if (t) {
                 std::string cs;
                 get_value_clean_upper(cs, t, 0, "pdb_format_compatible");
                 if (cs == "N") return false;
            }
       }
       return true;
}

void FileObj::write_mmcif_file(const bool& internal_flag, const bool& double_quoting_flag, const bool& start_date_flag)
{
       if (!_CifObj || _output_filename.empty()) return;

       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);

       deleteTable(_cifblock, "pdbx_pdb_compnd");
       deleteTable(_cifblock, "pdbx_pdb_source");
       deleteTable(_cifblock, "struct_biol_gen");
/*
       std::vector<std::string> categories;
       categories.clear();
       categories.push_back("entity_src_nat");
       categories.push_back("entity_src_gen");
       categories.push_back("pdbx_entity_src_syn");

       std::set<std::string> items;
       items.clear();
       items.insert("entity_id");

       for (std::vector<std::string>::const_iterator pos =
            categories.begin(); pos != categories.end(); ++pos) {
            ISTable *t = getTablePtr(_cifblock, *pos);
            if (is_empty_table(t, items)) deleteTable(_cifblock, *pos);
       }
*/

       _reformat_one_letter_code_sequence(_cifblock);
       // _switch_high_low_resolution(_cifblock); /* removed by request from DAOTHER-3616
       _remove_extra_space_in_name(_cifblock);
       _reorder_refine_ls_shell_reflns_shell(_cifblock, "refine_ls_shell");
    
       if (_DepUI_Flag) {
            // DAOTHER-4859, 4755, 4729, 4398
            const std::set<std::string>& remove_categories = CategoryMapping::depUI_remove_categories();
            for (std::set<std::string>::const_iterator spos = remove_categories.begin(); spos != remove_categories.end(); ++spos) {
                 deleteTable(_cifblock, *spos);
            }

            const std::map<std::string, std::set<std::string> >& clear_items = CategoryMapping::depUI_clear_items();
            for (std::map<std::string, std::set<std::string> >::const_iterator mpos = clear_items.begin(); mpos != clear_items.end(); ++mpos) {
                 ISTable *t = getTablePtr(_cifblock, mpos->first);
                 if (!t) continue;

                 int rowNo = t->GetNumRows();
                 for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                      if (!t->IsColumnPresent(*spos)) continue;

                      for (int i = 0; i < rowNo; ++i) {
                           t->UpdateCell(i, *spos, "");
                      }
                 }
                 _cifblock.WriteTable(t);
            }
/*
            std::vector<std::string> categories;
            categories.clear();
            categories.push_back("entity_src_nat");
            categories.push_back("entity_src_gen");
            categories.push_back("pdbx_entity_src_syn");
            for (std::vector<std::string>::const_iterator pos =
                 categories.begin(); pos != categories.end(); ++pos) {
                 ISTable *t = getTablePtr(_cifblock, *pos);
                 if (t) continue;
                 t = _newTablePtr(*pos);
                 t->AddRow();
                 _cifblock.WriteTable(t);
            }
*/
       }

       if (_dictUtil.Read()) {
            std::string caseSensitiveItemName;
            std::vector<std::string> TableNames;
            _cifblock.GetTableNames(TableNames);
            for (std::vector<std::string>::const_iterator tpos = TableNames.begin(); tpos != TableNames.end(); ++tpos) {
                 ISTable *t = getTablePtr(_cifblock, *tpos);
                 if (!t) continue;

                 std::vector<std::string> ColumnNames = t->GetColumnNames();
                 for (std::vector<std::string>::const_iterator cpos = ColumnNames.begin(); cpos != ColumnNames.end(); ++cpos) {
                      if (_dictUtil.isDefinedItem(*tpos, *cpos, caseSensitiveItemName)) {
                           if (!caseSensitiveItemName.empty()) rename_item(t, *cpos, caseSensitiveItemName);
                      }
                 }
                 _cifblock.WriteTable(t);
            }
       }

       if (start_date_flag) _add_start_date(_cifblock);

       std::vector<std::string> order_list;
       order_list.clear();

       if (internal_flag) {
            _legacy_vs_v5_conversion(_cifblock, MASTER_TO_INTERNAL);
       } else {
            std::string order_file = _rcsbroot  + "/data/ascii/order_list";
            get_order_list(order_file, _cifblock, order_list);
       }

       ISTable *t = getTablePtr(_cifblock, "atom_site");
       if (t) _write_atom_site_record_number = t->GetNumRows();

       if (_check_atom_site_record_number_flag && (_read_atom_site_record_number > 0)) {
            if (_read_atom_site_record_number != _write_atom_site_record_number) {
                 _logIo->message("Processing error: Reading %d atom_site records <--> Writing %d atom_site records.\n",
                                  _read_atom_site_record_number, _write_atom_site_record_number);
            }
       }

       if (!_resetFirstBlockName.empty()) {
             reset_entry_id(_cifblock, _resetFirstBlockName);
             _CifObj->RenameBlock(_firstBlockName, _resetFirstBlockName);
       }

       if (double_quoting_flag) _CifObj->SetQuoting(CifFile::eDOUBLE);

       if (!order_list.empty() /* && (_CifObj->GetNumBlocks() == 1) */ )
            _CifObj->Write(_output_filename, order_list);
       else _CifObj->Write(_output_filename);
}

void FileObj::_read_mmcif_metadata(const bool& coord_and_seq_only, const bool& update_database_PDB_rev_flag)
{
       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);

       _legacy_vs_v5_conversion(_cifblock, INTERNAL_TO_MASTER);

       _read_extra_dataBlock();

       _read_mmcif_metadata_process(_cifblock, coord_and_seq_only, update_database_PDB_rev_flag);
}

void FileObj::_legacy_vs_v5_conversion(Block& block, const int& status)
{
       std::string mapping_file = _rcsbroot + "/data/binary/master_internal_cif_mapping.odb";
       master_internal_conversion(mapping_file, block, status);
}

void FileObj::_read_extra_dataBlock()
{
       if (!_CifObj) return;

       std::vector<std::string> blockNames;
       _CifObj->GetBlockNames(blockNames);
       if (blockNames.size() < 2) return;

       bool has_comp_list = false;
       std::set<std::string> comp_block_set;
       comp_block_set.clear();
       for (unsigned int i = 1; i < blockNames.size(); ++i) {
            if (blockNames[i] == "comp_list") has_comp_list = true;
            else if (blockNames[i].substr(0, 4) == "comp") comp_block_set.insert(blockNames[i]);
       }

       std::list<std::string> comp_id_list;
       comp_id_list.clear();
       std::map<std::string, std::string> id_name_mapping, meta_data;
       id_name_mapping.clear();
       if (has_comp_list) {
            Block& block = _CifObj->GetBlock("comp_list");
            ISTable *t = getTablePtr(block, "chem_comp");
            if (t) {
                 std::string id, name;
                 int rowNo = t->GetNumRows();
                 for (int i = 0; i < rowNo; ++i) {
                      get_value_clean(id, t, i, "id");
                      if (id.empty()) continue;
                      std::string block_id = "comp_" + id;
                      if (comp_block_set.find(block_id) == comp_block_set.end()) continue;
                      comp_block_set.erase(block_id);

                      comp_id_list.push_back(id);
                      get_value_clean(name, t, i, "name");
                      if (!name.empty()) id_name_mapping.insert(std::make_pair(id, name));
                 }
            }
       }

       for (std::set<std::string>::const_iterator spos = comp_block_set.begin(); spos != comp_block_set.end(); ++spos) {
            comp_id_list.push_back(spos->substr(5));
       }

       ConnectFormat drug;
       for (std::list<std::string>::const_iterator lpos = comp_id_list.begin(); lpos != comp_id_list.end(); ++lpos) {
            Block& block = _CifObj->GetBlock("comp_" + *lpos);
            if (_verify_chem_comp_info(block, *lpos)) {
                 drug.clear();
                 if (drug.Read(block, false, false)) {
                      meta_data.clear();
                      meta_data.insert(std::make_pair("id", *lpos));
                      meta_data.insert(std::make_pair("three_letter_code", *lpos));
                      std::map<std::string, std::string>::const_iterator mpos = id_name_mapping.find(*lpos);
                      if (mpos != id_name_mapping.end()) meta_data.insert(std::make_pair("name", mpos->second));
                      drug.mergeMetaData(meta_data);
                      _ccDic->add_author_defined_drug(*lpos, drug);
/*
                      // Testing
                      CifFile *fobjOut = create_fobj("", *lpos);
                      drug.Write(*lpos, fobjOut->GetBlock(*lpos));
                      fobjOut->SetQuoting(CifFile::eDOUBLE);
                      fobjOut->Write(*lpos + ".cif");
                      delete fobjOut;
*/
                 }
            }
       }
}

bool FileObj::_verify_chem_comp_info(Block& block, const std::string& comp_id)
{
       std::map<std::string, std::string> error_map;
       error_map.clear();

       std::vector<std::string> items;
       std::vector<std::vector<std::string> > values;

       std::map<std::string, std::list<unsigned int> > atom_id_map;
       atom_id_map.clear();

       ISTable *t = getTablePtr(block, "chem_comp_atom");
       if (t) {
            items.clear();
            items.push_back("comp_id");
            items.push_back("atom_id");
            items.push_back("type_symbol");
            get_values(t, items, values);

            std::string error = _check_chem_comp_error(comp_id, items, values);

            for (unsigned int i = 0; i < values.size(); ++i) {
                 if (values[i][1].empty()) continue;
                 std::map<std::string, std::list<unsigned int> >::iterator mpos = atom_id_map.find(values[i][1]);
                 if (mpos != atom_id_map.end()) mpos->second.push_back(i);
                 else {
                      std::list<unsigned int> uList;
                      uList.clear();
                      uList.push_back(i);
                      atom_id_map.insert(std::make_pair(values[i][1], uList));
                 }
            }

            for (std::map<std::string, std::list<unsigned int> >::const_iterator mpos = atom_id_map.begin(); mpos != atom_id_map.end(); ++mpos) {
                 if (mpos->second.size() < 2) continue;
                 std::string rows = "";
                 for (std::list<unsigned int>::const_iterator lpos = mpos->second.begin(); lpos != mpos->second.end(); ++lpos) {
                      if (!rows.empty()) rows += ", ";
                      rows += String::IntToString(*lpos + 1);
                 }
                 error += "Atom name '" + mpos->first + "' repeats in rows ( " + rows + " ).\n";
            }
            if (!error.empty()) error_map.insert(std::make_pair("chem_comp_atom", error));
       } else error_map.insert(std::make_pair("chem_comp_atom", "missing"));

       t = getTablePtr(block, "chem_comp_bond");
       if (t) {
            items.clear();
            items.push_back("comp_id");
            items.push_back("atom_id_1");
            items.push_back("atom_id_2");
            if (t->IsColumnPresent("value_order")) items.push_back("value_order");
            else if (t->IsColumnPresent("type"))   items.push_back("type");

            get_values(t, items, values);

            std::string error = _check_chem_comp_error(comp_id, items, values);

            std::set<std::string> bonded_atom_set;
            bonded_atom_set.clear();

            for (unsigned int i = 0; i < values.size(); ++i) {
                 for (int j = 1; j < 3; ++j) {
                      if (values[i][j].empty() || (atom_id_map.find(values[i][j]) != atom_id_map.end())) continue;
                      error += "In row '" + String::IntToString(i + 1) + "', atom name '" + values[i][j] + "' in '" + items[j]
                                  + "' field is not defined in '_chem_comp_atom.atom_id' field.\n";
                 }
                 bonded_atom_set.insert(values[i][1]);
                 bonded_atom_set.insert(values[i][2]);
            }
            if (!error.empty()) error_map.insert(std::make_pair("chem_comp_bond", error));

            if (atom_id_map.size() > 1) {
                 error.clear();
                 for (std::map<std::string, std::list<unsigned int> >::const_iterator mpos = atom_id_map.begin(); mpos != atom_id_map.end(); ++mpos) {
                      if (bonded_atom_set.find(mpos->first) != bonded_atom_set.end()) continue;
                      error += "Atom name '" + mpos->first + "' is not linked to other atom(s).\n";
                 }
                 if (!error.empty()) {
                      std::map<std::string, std::string>::iterator mpos = error_map.find("chem_comp_atom");
                      if (mpos != error_map.end()) mpos->second += error;
                      else error_map.insert(std::make_pair("chem_comp_atom", error));
                 }
            }
       } else if (values.size() > 1) error_map.insert(std::make_pair("chem_comp_bond", "missing"));

       if (!error_map.empty()) {
            std::string error = "In 'data_comp_" + comp_id + "' block:\n";
            for (std::map<std::string, std::string>::const_iterator mpos = error_map.begin(); mpos != error_map.end(); ++mpos) {
                 if (mpos->second == "missing") error += "\nMissing '" + mpos->first + "' category.\n";
                 else error += "\nIn '" + mpos->first + "' category:\n" + mpos->second;
            }
            _error_messages.insertMessage("error", "model", error, false);
            return false;
       }
       return true;
}

std::string FileObj::_check_chem_comp_error(const std::string& comp_id, const std::vector<std::string>& items,
                                            const std::vector<std::vector<std::string> >& values) 
{
       std::string error = "";
       if (values.empty()) return error;

       std::set<std::string> comp_id_set;
       comp_id_set.clear();
       for (unsigned int i = 0; i < values.size(); ++i) {
            for (unsigned int j = 0; j < values[i].size(); ++j) {
                 if (!values[i][j].empty()) {
                      if (j == 0) comp_id_set.insert(values[i][j]);
                      continue;
                 }
                 error += "The value for '" + items[j] + "' field is empty in row '" + String::IntToString(i + 1) + "'.\n";
            }
       }

       for (std::set<std::string>::const_iterator spos = comp_id_set.begin(); spos != comp_id_set.end(); ++spos) {
            if (*spos == comp_id) continue;
            error += "The 'comp_id' field value '" + *spos + "' is not same as '" + comp_id + "' defined in data block name.\n";
       }
       return error;
}

void FileObj::_read_mmcif_metadata_process(Block& block, const bool& coord_and_seq_only, const bool& update_database_PDB_rev_flag)
{
       /* if (!coord_and_seq_only) */ _getExperimentTypefrom_exptl(block);

       _getSequenceInformation(block);

       if (!coord_and_seq_only) {
            _getCrySymmetryfrom_symmetry_and_cell(block);

            _getAuthorDefinedScaleMatrixfrom_atom_sites(block);

            if (update_database_PDB_rev_flag) _update_database_PDB_rev(block);

            _updateCrySymmetry();
       }

       std::string cs;

       ISTable *t = getTablePtr(block, "pdbx_data_processing_status");
       if (t) {
            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean_lower(cs, t, i, "status");
                 if (cs != "skip") continue;
                 get_value_clean_lower(cs, t, i, "task_name");
                 if (!cs.empty()) _skip_task_set.insert(cs); 
            }
       }

       t = getTablePtr(block, "struct_conn");
       if (t) {
            std::vector<std::string> items1, items2, data;
            items1.clear();
            items1.push_back("ptnr1_auth_asym_id");
            items1.push_back("ptnr1_auth_comp_id");
            items1.push_back("ptnr1_auth_seq_id");
            items1.push_back("pdbx_ptnr1_PDB_ins_code");
            items1.push_back("ptnr2_auth_asym_id");
            items1.push_back("ptnr2_auth_comp_id");
            items1.push_back("ptnr2_auth_seq_id");
            items1.push_back("pdbx_ptnr2_PDB_ins_code");
            items2.clear();
            items2.push_back("ptnr2_auth_asym_id");
            items2.push_back("ptnr2_auth_comp_id");
            items2.push_back("ptnr2_auth_seq_id");
            items2.push_back("pdbx_ptnr2_PDB_ins_code");
            items2.push_back("ptnr1_auth_asym_id");
            items2.push_back("ptnr1_auth_comp_id");
            items2.push_back("ptnr1_auth_seq_id");
            items2.push_back("pdbx_ptnr1_PDB_ins_code");

            std::string type, atom_id_1, atom_id_2;
            std::vector<unsigned int> delete_rows;
            delete_rows.clear();

            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean(cs, t, i, "ptnr1_symmetry");
                 if (!cs.empty() && cs != "1_555") continue;
                 get_value_clean(cs, t, i, "ptnr2_symmetry");
                 if (!cs.empty() && cs != "1_555") continue;

                 get_value_clean_lower(type, t, i, "conn_type_id");

                 data.clear();
                 for (std::vector<std::string>::const_iterator pos = items1.begin(); pos != items1.end(); ++pos) {
                      get_value_clean(cs, t, i, *pos);
                      data.push_back(cs);
                 }

                 if (type == "covale") {
                      try {
                           get_value_clean(atom_id_1, t, i, "ptnr1_label_atom_id");
                           get_value_clean(atom_id_2, t, i, "ptnr2_label_atom_id");
                           const ConnectFormat& drug_1 = _ccDic->find_drug(data[1]);
                           const ConnectFormat& drug_2 = _ccDic->find_drug(data[5]);
                           const AtomFormat& Atom1 = drug_1.find_atom(atom_id_1);
                           const AtomFormat& Atom2 = drug_2.find_atom(atom_id_2);
                           if ((drug_1.getMetaData("pdbx_type") == "ATOMS") && (drug_2.getMetaData("pdbx_type") == "ATOMS")) {
                                if (((Atom1.atomtype() == "C") && (Atom2.atomtype() == "C")) || ((Atom1.atomtype() == "O") && (Atom2.atomtype() == "O"))) {
                                     delete_rows.push_back(i);
                                     continue;
                                } // else if ((Atom1.atomtype() != "C") && (Atom2.atomtype() != "C")) continue;
                                // only for C-polar_atom bond (in case of D_1000102005 2L65. Should be covered by BondUtil::is_glycosidic_linkage() function.
                           }

                           if (!BondUtil::is_glycosidic_linkage(atom_id_1, Atom1.atomtype(), atom_id_2, Atom2.atomtype())) continue;
                           if (_is_glycosylation_site(CompositeIndex::getIndex(data[1], atom_id_1, data[5], atom_id_2))) continue;

                      } catch (const std::exception& exc) {}
                 } else continue;

                 std::string idx = CompositeIndex::getIndex(data);
                 _link_residue_set.insert(idx);

                 data.clear();
                 for (std::vector<std::string>::const_iterator pos = items2.begin(); pos != items2.end(); ++pos) {
                      get_value_clean(cs, t, i, *pos);
                      data.push_back(cs);
                 }
                 idx = CompositeIndex::getIndex(data);
                 _link_residue_set.insert(idx);
            }
            if (!delete_rows.empty()) {
                 t->DeleteRows(delete_rows);
                 block.WriteTable(t);
            }
       }
}

void FileObj::_getExperimentTypefrom_exptl(Block& block)
{
       ISTable *t = getTablePtr(block, "exptl");
       if (t) {
            std::string cs;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_upper(cs, t, i, "method");
                 if (cs.empty()) continue;
                 _getExperimentType(cs); 
            }
       }
       if (!_experiment_type) _experiment_type = EXPERIMENT_TYPE_BASIC;
       
}

void FileObj::_getSequenceInformation(Block& block)
{
       _getSequenceInformationfrom_entity_poly(block);

       _getSequenceInformationfrom_pdbx_poly_seq_scheme(block);

       _getSequenceType();
}

void FileObj::_getSequenceInformationfrom_pdbx_poly_seq_scheme(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_poly_seq_scheme");
       if (!t) return;

       if (!t->IsColumnPresent("mon_id") || !t->IsColumnPresent("seq_id") || !t->IsColumnPresent("pdb_strand_id") || !t->IsColumnPresent("pdb_mon_id") ||
           !t->IsColumnPresent("pdb_seq_num") || !t->IsColumnPresent("pdb_ins_code")) return;

       std::set<std::string> FirstIndex;
       FirstIndex.clear();

       _Seq_Scheme_Mapping.clear();

       std::string cs, cs1;
       std::vector<std::string> data, items;
       std::vector<std::vector<std::string> > t_vectors;

       items.clear();
       items.push_back("pdb_strand_id");
       items.push_back("mon_id");
       items.push_back("pdb_mon_id");
       items.push_back("pdb_seq_num");
       items.push_back("pdb_ins_code");
       items.push_back("entity_id");

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean_lower(cs, t, i, "hetero");
            if (cs == "y") {
                 get_value(cs1, t, i, "seq_id");
                 cs1 += "_" + cs;
                 if (FirstIndex.find(cs1) != FirstIndex.end()) continue;
                 FirstIndex.insert(cs1);
            }

            data.clear();
            for (unsigned int j = 0; j < items.size(); ++j) {
                 get_value_clean(cs, t, i, items[j]);
                 data.push_back(cs);
            }

            std::map<std::string, std::vector<std::vector<std::string> > >::iterator mpos = _Seq_Scheme_Mapping.find(data[0]);
            if (mpos == _Seq_Scheme_Mapping.end()) {
                 t_vectors.clear();
                 t_vectors.push_back(data);
                 _Seq_Scheme_Mapping.insert(std::make_pair(data[0], t_vectors));
                 FirstIndex.clear();
            } else mpos->second.push_back(data);
       }
}

void FileObj::_getSequenceInformationfrom_entity_poly(Block& block)
{
       ISTable *t = getTablePtr(block, "entity_poly");
       if (!t) {
            // deleteTable(block, "pdbx_poly_seq_scheme");
            return;
       }

       // bool found = false;
       std::string cs;
       int rowNo = t->GetNumRows();
/*
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean_lower(cs, t, i, "pdbx_sequence_evidence_code");
            if (cs == "derived from coordinates") {
                 found = true;
                 break;
            }
       }
       if (found) {
            deleteTable(block, "pdbx_poly_seq_scheme");
            deleteTable(block, "entity_poly");
            return;
       }
*/

       if (/* !t->IsColumnPresent("type") || */ !t->IsColumnPresent("pdbx_strand_id")) {
            _logIo->message("Missing '_entity_poly.pdbx_strand_id' item.\n");
           _error_messages.insertMessage("error", "model", "Missing '_entity_poly.pdbx_strand_id' item.", true);
            return;
       }
       if (!t->IsColumnPresent("pdbx_seq_one_letter_code")) {
            _logIo->message("Missing '_entity_poly.pdbx_seq_one_letter_code' item.\n");
            _error_messages.insertMessage("error", "model", "Missing '_entity_poly.pdbx_seq_one_letter_code' item.", true);
       }
       if (!t->IsColumnPresent("pdbx_seq_one_letter_code") &&
           !t->IsColumnPresent("pdbx_seq_three_letter_code")) return;

       std::string entity_id, cs1, cs2, PolyType, EvidenceCode, error;
       std::vector<std::string> data;
       SEQ seq;

       error.clear();
       _Seqs.clear();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(entity_id, t, i, "entity_id");
            get_value_clean(PolyType, t, i, "type");
            get_value_clean(EvidenceCode, t, i, "pdbx_sequence_evidence_code");
            if (EvidenceCode.empty()) EvidenceCode = "depositor provided";

            get_value_clean_upper(cs, t, i, "pdbx_seq_one_letter_code");
            get_value_clean_upper(cs1, t, i, "pdbx_N_terminal_seq_one_letter_code");
            get_value_clean_upper(cs2, t, i, "pdbx_C_terminal_seq_one_letter_code");
            std::string sequence = cs1 + cs + cs2;

            data.clear();
            if (!sequence.empty()) {
                 std::string err = SeqCodeUtil::GetThreeLetterCodeSeq(sequence, PolyType, data);
                 if (!err.empty()) error += "In 'entity_poly' category, one letter sequence " + err + " for entity '" + entity_id + "'.\n";
            } else if (t->IsColumnPresent("pdbx_seq_three_letter_code")) {
                 get_value_upper(sequence, t, i, "pdbx_seq_three_letter_code");
                 if (!sequence.empty()) get_wordarray(data, sequence, " \t\n");
            }
            if (data.empty()) continue;

            _clearSEQ(seq);
            seq.entity_id = entity_id;
            seq.res = data;
            seq.PolyType = PolyType;
            seq.EvidenceCode = EvidenceCode;
            get_value(cs, t, i, "pdbx_strand_id");
            get_wordarray(data, cs, ", \t\n");
            for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
                 seq.ChainId = *pos;
                 seq.PDB_ChainId = *pos;
                 seq.order = (int) _Seqs.size();
                 _Seqs.insert(std::make_pair(*pos, seq));
            }
       }

       if (!error.empty()) throw std::out_of_range(error);
}

void FileObj::_getCrySymmetryfrom_symmetry_and_cell(Block& block)
{
       // re-set values in cell, symmetry & atom_sites categories for EM & NMR entries, required by DAOTHER-932
       // removed re-set procedure

       std::string cs;
       ISTable *t = NULL;

       if (_update_symmetry_flag)
            t = _getTablePtr(block, "cell");
       else t = getTablePtr(block, "cell");
       if (t) {
            get_value(cs, t, 0, "length_a");    if (!cs.empty()) _cell.setA(cs);
            get_value(cs, t, 0, "length_b");    if (!cs.empty()) _cell.setB(cs);
            get_value(cs, t, 0, "length_c");    if (!cs.empty()) _cell.setC(cs);
            get_value(cs, t, 0, "angle_alpha"); if (!cs.empty()) _cell.setAlpha(cs);
            get_value(cs, t, 0, "angle_beta");  if (!cs.empty()) _cell.setBeta(cs);
            get_value(cs, t, 0, "angle_gamma"); if (!cs.empty()) _cell.setGamma(cs);
       }

       if (_update_symmetry_flag)
            t = _getTablePtr(block, "symmetry");
       else t = getTablePtr(block, "symmetry");
       if (t) {
            get_value(cs, t, 0, "space_group_name_H-M");
            if (cs.empty()) get_value(cs, t, 0, "pdbx_full_space_group_name_H-M");
            if (!cs.empty()) {
                 std::string number;
                 SpaceGroup::CheckSpaceGroup(cs, number);
                 if (_update_symmetry_flag) {
                      t->UpdateCell(0, "space_group_name_H-M", cs);
                      if (!number.empty()) t->UpdateCell(0, "Int_Tables_number", number);
                      block.WriteTable(t);
                 }
                 _cell.setSpaceGroup(cs);
            }
       }
}

void FileObj::_getAuthorDefinedScaleMatrixfrom_atom_sites(Block& block)
{
       const char *_matrix_items[3][3] = {
            { "fract_transf_matrix[1][1]", "fract_transf_matrix[1][2]", "fract_transf_matrix[1][3]" },
            { "fract_transf_matrix[2][1]", "fract_transf_matrix[2][2]", "fract_transf_matrix[2][3]" },
            { "fract_transf_matrix[3][1]", "fract_transf_matrix[3][2]", "fract_transf_matrix[3][3]" }
       };
       const char *_vector_items[3] = {
               "fract_transf_vector[1]", "fract_transf_vector[2]", "fract_transf_vector[3]"
       };

       _exist_scale_matrix = false;

       ISTable *t = getTablePtr(block, "atom_sites");
       if (!t) return;

       std::string cs;
       for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                 get_value_clean(cs, t, 0, _matrix_items[i][j]);
                 if (cs.empty()) return;
                 _scale[i][j] = atof(cs.c_str());
            }
            get_value_clean(cs, t, 0, _vector_items[i]);
            if (cs.empty()) return;
            _vec[i] = atof(cs.c_str());
       }
       _exist_scale_matrix = true;
}

void FileObj::_update_database_PDB_rev(Block& block)
{
       ISTable *t = getTablePtr(block, "database_PDB_rev");
       if (!t) {
            t = getTablePtr(block, "pdbx_database_status");
            if (!t) return;

            std::string cs;
            get_value_clean(cs, t, 0, "recvd_initial_deposition_date");
            if (cs.empty()) return;

            t = _newTablePtr("database_PDB_rev");
            t->AddRow();
            t->UpdateCell(0, "date_original", cs);
            block.WriteTable(t);

            return;
       }

       std::vector<std::string> items, data;
       items.clear();
       items.push_back("pdbx_record_revised_1");
       items.push_back("pdbx_record_revised_2");
       items.push_back("pdbx_record_revised_3");
       items.push_back("pdbx_record_revised_4");

       std::map<int, std::vector<std::string> > rev_records;
       rev_records.clear();

       std::set<std::string> index_set;
       std::string cs;
      
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            index_set.clear();
            data.clear();
            for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
                 get_value_clean(cs, t, i, *pos);
                 if (cs.empty()) continue;
                 if (index_set.find(cs) != index_set.end()) continue;
                 index_set.insert(cs);
                 data.push_back(cs); 
            }
            if (data.empty()) continue;

            get_value_clean(cs, t, i, "num");
            if (cs.empty()) continue;

            rev_records.insert(std::make_pair(atoi(cs.c_str()), data));
       }

       for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
            if (t->IsColumnPresent(*pos)) t->DeleteColumn(*pos);
       }
       block.WriteTable(t);

       if (rev_records.empty()) return;

       t = getTablePtr(block, "database_PDB_rev_record");
       if (t) {
            std::string cs1;
            index_set.clear();
            rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "rev_num");
                 get_value_clean(cs1, t, i, "type");
                 if (index_set.find(cs + "_" + cs1) != index_set.end()) continue;

                 index_set.insert(cs + "_" + cs1);
                 int num = atoi(cs.c_str());
                 std::map<int, std::vector<std::string> >::iterator mpos = rev_records.find(num);
                 if (mpos != rev_records.end())
                      mpos->second.push_back(cs1);
                 else {
                      data.clear();
                      data.push_back(cs1);
                      rev_records.insert(std::make_pair(num, data));
                 }
            }
       }

       t = _newTablePtr("database_PDB_rev_record");
       int row = 0;
       for (std::map<int, std::vector<std::string> >::const_iterator mpos = rev_records.begin(); mpos != rev_records.end(); ++mpos) {
            cs = String::IntToString(mpos->first);
            for (std::vector<std::string>::const_iterator pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                 t->AddRow();
                 t->UpdateCell(row, "rev_num", cs);
                 t->UpdateCell(row, "type", *pos);
                 row++;
            }
       }
       block.WriteTable(t);
}

void FileObj::_reformat_one_letter_code_sequence(Block& block)
{
       std::map<std::string, std::set<std::string> > seq_categories;
       seq_categories.clear();

       std::set<std::string> items;

       items.clear();
       items.insert("pdbx_seq_one_letter_code");
       seq_categories.insert(std::make_pair("struct_ref", items));

       items.clear();
       items.insert("db_seq_one_letter_code");
       seq_categories.insert(std::make_pair("pdbx_struct_ref_seq_depositor_info", items));

       items.clear();
       items.insert("one_letter_code");
       items.insert("one_letter_code_mod");
       seq_categories.insert(std::make_pair("pdbx_seq_map_depositor_info", items));

       std::vector<std::string> data;
       std::string cs;

       std::string error = "";
       for (std::map<std::string, std::set<std::string> >::const_iterator mpos = seq_categories.begin(); mpos != seq_categories.end(); ++mpos) {
            ISTable *t = getTablePtr(block, mpos->first);
            if (!t) continue;

            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                      get_value_clean(cs, t, i, *spos);
                      if (cs.empty()) continue;
                      std::string err = SeqCodeUtil::GetOneLetterCodeSeqArray(cs, data);
                      if (!err.empty()) error += "In '_" + mpos->first + "." + *spos + "', one letter sequence '" + cs + "' " + err + ".\n";
                      cs = SeqCodeUtil::GetOneLetterCodeSeq(data, 80);
                      t->UpdateCell(i, *spos, cs);
                 }
            }
            block.WriteTable(t);
       }

       if (!error.empty()) throw std::out_of_range(error);
}

void FileObj::_switch_high_low_resolution(Block& block)
{
       const char *_switch_categories[N_SWITCH_CATEGORY][3] = {
            { "refine",          "ls_d_res_high",     "ls_d_res_low"     },
            { "reflns",          "d_resolution_high", "d_resolution_low" },
            { "refine_ls_shell", "d_res_high",        "d_res_low"        },
            { "reflns_shell",    "d_res_high",        "d_res_low"        }
       };

       std::string high, low;

       for (int i = 0; i < N_SWITCH_CATEGORY; ++i) {
            ISTable *t = getTablePtr(block, _switch_categories[i][0]);
            if (!t) continue;

            int rowNo = t->GetNumRows();
            for (int j = 0; j < rowNo; ++j) {
                 get_value_clean(high, t, j, _switch_categories[i][1]);
                 if (high.empty()) continue;
                 get_value_clean(low, t, j, _switch_categories[i][2]);
                 if (low.empty()) continue;

                 if (atof(high.c_str()) > atof(low.c_str())) {
                      t->UpdateCell(j, _switch_categories[i][1], low);
                      t->UpdateCell(j, _switch_categories[i][2], high);
                 }
            }
            block.WriteTable(t);
       }
}

void FileObj::_remove_extra_space_in_name(Block& block)
{
       const char *_name_categories[N_NAME_CATEGORY][2] = {
            { "audit_author",    "name"  },
            { "citation_author", "name"  },
            { "citation_editor", "name"  }
       };

       std::string cs;

       for (int i = 0; i < N_NAME_CATEGORY; ++i) {
            ISTable *t = getTablePtr(block, _name_categories[i][0]);
            if (!t) continue;

            int rowNo = t->GetNumRows();
            for (int j = 0; j < rowNo; ++j) {
                 get_value(cs, t, j, _name_categories[i][1]);
                 if (cs.empty()) continue;

                 String::StripAndCompressWs(cs);
                 replace_string(cs, ". ", ".");
                 t->UpdateCell(j, _name_categories[i][1], cs);
            }
            block.WriteTable(t);
       }
}

void FileObj::_reorder_refine_ls_shell_reflns_shell(Block& block, const std::string& category)
{
       ISTable *t = getTablePtr(block, category);
       if (!t) return;

       if (!t->IsColumnPresent("d_res_high")) return;
 
       int rowNo = t->GetNumRows();
       if (rowNo < 2) return;

       std::vector<std::string> pdbx_refine_id;
       pdbx_refine_id.clear();

       // first key: pdbx_refine_id
       // second key: d_res_high
       std::map<std::string, std::multimap<double, std::vector<std::string> > > values;
       values.clear();

       std::multimap<double, std::vector<std::string> > tmp_map;
       std::vector<std::string> data;
       std::string cs, id;

       const std::vector<std::string>& items = t->GetColumnNames();
       for (int i = 0; i < rowNo; ++i) {
            data.clear();
            for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
                 get_value_clean(cs, t, i, *pos);
                 if (cs.empty()) cs = ".";
                 data.push_back(cs);
            }
            get_value_clean(cs, t, i, "d_res_high");
            get_value_clean(id, t, i, "pdbx_refine_id");
            if (id.empty()) id = "NULL";

            std::map<std::string, std::multimap<double, std::vector<std::string> > >::iterator mpos = values.find(id);
            if (mpos != values.end())
                 mpos->second.insert(std::make_pair(atof(cs.c_str()), data));
            else {
                 pdbx_refine_id.push_back(id);
                 tmp_map.clear();
                 tmp_map.insert(std::make_pair(atof(cs.c_str()), data));
                 values.insert(std::make_pair(id, tmp_map));
            }
       }

       int i = 0;
       for (std::vector<std::string>::const_iterator pos = pdbx_refine_id.begin(); pos != pdbx_refine_id.end(); ++pos) {
            std::map<std::string, std::multimap<double, std::vector<std::string> > >::const_iterator mpos = values.find(*pos);
            if (mpos == values.end()) continue;

            for (std::multimap<double, std::vector<std::string> >::const_iterator dpos = mpos->second.begin(); dpos != mpos->second.end(); ++dpos) {
                 for (unsigned int j = 0; j < items.size(); ++j) {
                      t->UpdateCell(i, items[j], dpos->second[j]);
                 }
                 i++;
            }
       }
       block.WriteTable(t);
}

void FileObj::_add_start_date(Block& block)
{
       ISTable *t = _getTablePtr(block, "pdbx_database_status");
       if (!t) return;

       std::string cs;
       get_value_clean(cs, t, 0, "date_begin_processing");
       if (!cs.empty()) return;

       get_value_clean_upper(cs, t, 0, "status_code");
       if (cs != "PROC" && cs != "REPL" && cs != "WAIT") return;

       // find existing start date
       cs.clear();
       ISTable *t1 = getTablePtr(block, "pdbx_database_proc");
       if (t1) {
            int rowNo = t1->GetNumRows(); 
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t1, i, "date_begin_cycle");
                 if (!cs.empty()) break;
            }
       }
       if (cs.empty()) get_date(cs);

       t->UpdateCell(0, "date_begin_processing", cs);
       block.WriteTable(t);
}

ISTable* FileObj::_getTablePtr(Block& block, const std::string& catName)
{
       ISTable *t = getTablePtr(block, catName);
       if (!t) return NULL;

       try {
            const CIF_CATEGORY& category = CategoryMapping::find_category(catName);
            check_missing_item(t, category.itemNames);
       } catch (const std::exception& exc) {}

       return t;
}

ISTable* FileObj::_newTablePtr(const CIF_CATEGORY& cifcategory)
{
       return add_new_table(cifcategory.category, cifcategory.itemNames);
}

ISTable* FileObj::_newTablePtr(const std::string& catName)
{
       try {
            const CIF_CATEGORY& category = CategoryMapping::find_category(catName);
            return _newTablePtr(category);
       } catch (const std::exception& exc) {
            if (_dictUtil.Read()) {
                 const std::vector<std::string>& itemNames = _dictUtil.getItems(catName);
                 if (!itemNames.empty()) return add_new_table(catName, itemNames);
            }
            throw std::out_of_range("Category " + catName + " is unrecognized.\n");
       }
       return NULL;
}

void FileObj::_read_branch_polymer_info()
{
       Block& block = _CifObj->GetBlock(_firstBlockName);

       ISTable *t1 = getTablePtr(block, "pdbx_entity_branch_list");
       if (!t1) return;
       ISTable *t2 = getTablePtr(block, "pdbx_entity_branch_link");
       if (!t2) return;
       ISTable *t3 = getTablePtr(block, "pdbx_entity_branch_descriptor");
       ISTable *t4 = getTablePtr(block, "pdbx_branch_scheme");
       if (!t4) return;

       std::map<std::string, std::vector<std::vector<std::string> > > list_mapping, scheme_mapping;

       std::vector<std::string> items;

       items.clear();
       items.push_back("entity_id");
       items.push_back("comp_id");
       items.push_back("num");
       items.push_back("hetero");
       _read_branch_scheme_list_category(t1, items, list_mapping);
       if (list_mapping.empty()) return;

       std::map<std::string, std::vector<std::map<std::string, std::string> > > link_mapping, descriptor_mapping;

       items.clear();
       items.push_back("entity_branch_list_num_1");
       items.push_back("comp_id_1");
       items.push_back("atom_id_1");
       items.push_back("leaving_atom_id_1");
       items.push_back("entity_branch_list_num_2");
       items.push_back("comp_id_2");
       items.push_back("atom_id_2");
       items.push_back("leaving_atom_id_2");
       items.push_back("value_order");
       items.push_back("details");
       _read_branch_link_descriptor_category(t2, items, link_mapping);
       if (link_mapping.empty()) return;

       _check_branch_list_with_link(link_mapping, list_mapping);
       if (list_mapping.empty()) return;

       descriptor_mapping.clear();
       if (t3) {
            items.clear();
            items.push_back("descriptor");
            items.push_back("type");
            items.push_back("program");
            items.push_back("program_version");
            _read_branch_link_descriptor_category(t3, items, descriptor_mapping);
       }

       items.clear();
       items.push_back("pdb_asym_id");
       items.push_back("pdb_mon_id");
       items.push_back("pdb_seq_num");
       items.push_back("entity_id");
       items.push_back("hetero");
       _read_branch_scheme_list_category(t4, items, scheme_mapping);
       if (scheme_mapping.empty()) return;

       BRANCH_INFO branch_info;
       for (std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator
            smpos = scheme_mapping.begin(); smpos != scheme_mapping.end(); ++smpos) {
            std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator lmpos = list_mapping.find(smpos->second[0][3]);
            if ((lmpos == list_mapping.end()) || (smpos->second.size() != lmpos->second.size())) {
                 _Branch_Seq_Scheme_Mapping.clear();
                 return;
            }
            for (unsigned int i = 0; i < smpos->second.size(); ++i) {
                 if ((smpos->second[i][1] != lmpos->second[i][1]) || (smpos->second[i][1] != lmpos->second[i][1])) {
                      _Branch_Seq_Scheme_Mapping.clear();
                      return;
                 }
            }

            std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator link_mpos = link_mapping.find(lmpos->first);
            if (link_mpos == link_mapping.end()) {
                 _Branch_Seq_Scheme_Mapping.clear();
                 return;
            }

            branch_info.PDB_ChainId = smpos->first;
            branch_info.seqs = smpos->second;
            branch_info.linkages = link_mpos->second;
            branch_info.descriptors.clear();

            std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator dmpos = descriptor_mapping.find(lmpos->first);
            if (dmpos != descriptor_mapping.end()) branch_info.descriptors = dmpos->second;

            _Branch_Seq_Scheme_Mapping.insert(std::make_pair(smpos->first, branch_info));
       }
}

void FileObj::_read_branch_scheme_list_category(ISTable* t, const std::vector<std::string>& items, std::map<std::string,
                                                std::vector<std::vector<std::string> > >& mapping)
{
       mapping.clear();
/*
       std::set<std::string> FirstIndex;
       FirstIndex.clear();
*/
       std::string cs /*, cs1 */ ;
       std::vector<std::string> data;
       std::vector<std::vector<std::string> > t_vectors;

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
/*
            get_value_clean_lower(cs, t, i, "hetero");
            if (cs == "y") {
                 get_value(cs1, t, i, "num");
                 cs1 += "_" + cs;
                 if (FirstIndex.find(cs1) != FirstIndex.end()) continue;
                 FirstIndex.insert(cs1);
            }
*/
            data.clear();
            for (std::vector<std::string>::const_iterator vpos = items.begin(); vpos != items.end(); ++vpos) {
                 get_value_clean(cs, t, i, *vpos);
                 data.push_back(cs);
            }

            std::map<std::string, std::vector<std::vector<std::string> > >::iterator mpos = mapping.find(data[0]);
            if (mpos == mapping.end()) {
                 t_vectors.clear();
                 t_vectors.push_back(data);
                 mapping.insert(std::make_pair(data[0], t_vectors));
/*
                 FirstIndex.clear();
*/
            } else mpos->second.push_back(data);
       }
}

void FileObj::_read_branch_link_descriptor_category(ISTable*t, const std::vector<std::string>& items, std::map<std::string,
                                                    std::vector<std::map<std::string, std::string> > >& mapping)
{
       mapping.clear();

       std::string entity_id, cs;
       std::map<std::string, std::string> tmp_map;
       std::vector<std::map<std::string, std::string> > t_vectors;

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(entity_id, t, i, "entity_id");

            tmp_map.clear();
            for (std::vector<std::string>::const_iterator vpos = items.begin(); vpos != items.end(); ++vpos) {
                 get_value_clean(cs, t, i, *vpos);
                 if (!cs.empty()) tmp_map.insert(std::make_pair(*vpos, cs));
            }

            std::map<std::string, std::vector<std::map<std::string, std::string> > >::iterator mpos = mapping.find(entity_id);
            if (mpos == mapping.end()) {
                 t_vectors.clear();
                 t_vectors.push_back(tmp_map);
                 mapping.insert(std::make_pair(entity_id, t_vectors));
            } else mpos->second.push_back(tmp_map);
       }
}

void FileObj::_check_branch_list_with_link(const std::map<std::string, std::vector<std::map<std::string, std::string> > >& link_mapping,
                                           std::map<std::string, std::vector<std::vector<std::string> > >& list_mapping)
{
       std::set<std::string> remove_entity_set;
       remove_entity_set.clear();

       std::set<int> tmp_set;
       std::vector<std::set<int> > groups;
       std::map<std::string, int> idx_map;
       std::vector<std::pair<int, int> > links;
       for (std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator mpos = list_mapping.begin(); mpos != list_mapping.end(); ++mpos) {
            std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator mpos1 = link_mapping.find(mpos->first);
            if (mpos1 == link_mapping.end()) {
                 remove_entity_set.insert(mpos->first);
                 continue;
            }

            int idx = 0;
            groups.clear();
            idx_map.clear();
            for (std::vector<std::vector<std::string> >::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 idx_map.insert(std::make_pair((*vpos)[1] + "_" + (*vpos)[2], idx));
                 tmp_set.clear();
                 tmp_set.insert(idx);
                 groups.push_back(tmp_set);
                 idx++;
            }

            links.clear();
            for (std::vector<std::map<std::string, std::string> >::const_iterator vpos = mpos1->second.begin(); vpos != mpos1->second.end(); ++vpos) {
                 std::map<std::string, std::string>::const_iterator vmpos1 = vpos->find("comp_id_1");
                 if (vmpos1 == vpos->end()) continue;
                 std::map<std::string, std::string>::const_iterator vmpos2 = vpos->find("entity_branch_list_num_1");
                 if (vmpos2 == vpos->end()) continue;
                 std::map<std::string, int>::const_iterator impos1 = idx_map.find(vmpos1->second + "_" + vmpos2->second);
                 if (impos1 == idx_map.end()) continue;
                
                 vmpos1 = vpos->find("comp_id_2");
                 if (vmpos1 == vpos->end()) continue;
                 vmpos2 = vpos->find("entity_branch_list_num_2");
                 if (vmpos2 == vpos->end()) continue;
                 std::map<std::string, int>::const_iterator impos2 = idx_map.find(vmpos1->second + "_" + vmpos2->second);
                 if (impos2 == idx_map.end()) continue;

                 links.push_back(std::make_pair(impos1->second, impos2->second));
            }

            clustering_with_merging(groups, links);
            if (groups.size() > 1) remove_entity_set.insert(mpos->first);
       }

       if (remove_entity_set.empty()) return;

       for (std::set<std::string>::const_iterator spos = remove_entity_set.begin(); spos != remove_entity_set.end(); ++spos) {
            list_mapping.erase(*spos);
       }
}
