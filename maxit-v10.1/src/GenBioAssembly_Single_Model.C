/*
FILE:     GenBioAssembly_Single_Model.C
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
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>

#include "GenBioAssembly.h"
#include "utillib.h"

void GenBioAssembly::_update_non_atom_site_categories(Block& modelBlock, const std::map<int, std::vector<std::string> >& original_asym_chain_id_mapping,
                                                      Block& biolBlock, RCSB::Molecule* mol)
{
       // key: entity key
       // value: new entity ID
       std::map<std::string, int> entity_key_id_mapping;
       entity_key_id_mapping.clear();

       // key: new entity ID
       // pair.first: old entity ID
       // pair.second: pdbx_number_of_molecules 
       std::map<int, std::pair<std::string, int> > new_old_entity_id_mapping;
       new_old_entity_id_mapping.clear();

       // key: new entity ID
       // pair.first: old entity ID
       // pair.second: PDB chain ID list
       std::map<int, std::pair<std::string, std::vector<std::string> > > polymer_entity_id_mapping; 
       polymer_entity_id_mapping.clear();

       // key: new entity ID
       // pair.first: old entity ID
       // pair.second: updating key/value mapping
       std::map<int, std::pair<std::string, std::map<std::string, std::string> > > branch_entity_id_mapping, nonpolymer_entity_id_mapping;
       branch_entity_id_mapping.clear();
       nonpolymer_entity_id_mapping.clear();

       // key: ordinal order, no meaning
       // pair.first: old asym ID
       // pair.second: updating key/value mapping
       std::map<int, std::pair<std::string, std::map<std::string, std::string> > > polymer_scheme_mapping, branch_scheme_mapping, nonpolymer_scheme_mapping;
       polymer_scheme_mapping.clear();
       branch_scheme_mapping.clear();
       nonpolymer_scheme_mapping.clear();

       std::map<std::string, std::string> key_value_mapping;

       std::vector<std::map<std::string, std::string> > chain_mapping_values;
       chain_mapping_values.clear();

       std::vector<std::string> data;

       int entity_id = 0;
       RCSB::Chain* chain = mol->GetFirstChain();
       while (chain) {
            int number = 1;
            if ((chain->entity_key() == "HOH") || (chain->entity_key() == "DOD")) number = chain->ResidueNumbers();

            std::map<std::string, int>::const_iterator epos = entity_key_id_mapping.find(chain->entity_key());
            if (epos != entity_key_id_mapping.end()) {
                 chain->set_entity_id(epos->second);
                 std::map<int, std::pair<std::string, int> >::iterator nopos = new_old_entity_id_mapping.find(epos->second);
                 if (nopos != new_old_entity_id_mapping.end()) nopos->second.second += number;
                 if ((chain->chain_type() == "ATOMN") || (chain->chain_type() == "ATOMP")) {
                      std::map<int, std::pair<std::string, std::vector<std::string> > >::iterator pepos = polymer_entity_id_mapping.find(epos->second);
                      if (pepos != polymer_entity_id_mapping.end()) pepos->second.second.push_back(chain->PDB_ChainID());
                 }
            } else {
                 entity_id++;
                 entity_key_id_mapping.insert(std::make_pair(chain->entity_key(), entity_id));
                 new_old_entity_id_mapping.insert(std::make_pair(entity_id, std::make_pair(chain->entity_id(), number)));
                 if ((chain->chain_type() == "ATOMN") || (chain->chain_type() == "ATOMP")) {
                      data.clear();
                      data.push_back(chain->PDB_ChainID());
                      polymer_entity_id_mapping.insert(std::make_pair(entity_id, std::make_pair(chain->entity_id(), data)));
                 } else if (chain->chain_type() == "ATOMS") {
                      key_value_mapping.clear();
                      key_value_mapping.insert(std::make_pair("entity_id", String::IntToString(entity_id)));
                      branch_entity_id_mapping.insert(std::make_pair(entity_id, std::make_pair(chain->entity_id(), key_value_mapping)));
                 } else {
                      key_value_mapping.clear();
                      key_value_mapping.insert(std::make_pair("entity_id", String::IntToString(entity_id)));
                      nonpolymer_entity_id_mapping.insert(std::make_pair(entity_id, std::make_pair(chain->entity_id(), key_value_mapping)));
                 }
                 chain->set_entity_id(entity_id);
            }

            std::map<int, std::vector<std::string> >::const_iterator mpos = original_asym_chain_id_mapping.find(chain->index());
            if (mpos != original_asym_chain_id_mapping.end()) {
                 key_value_mapping.clear();
                 key_value_mapping.insert(std::make_pair("entity_id", chain->entity_id()));
                 key_value_mapping.insert(std::make_pair("label_asym_id", chain->ChainID()));
                 key_value_mapping.insert(std::make_pair("auth_asym_id", chain->PDB_ChainID()));
                 key_value_mapping.insert(std::make_pair("orig_label_asym_id", mpos->second[0]));
                 key_value_mapping.insert(std::make_pair("orig_auth_asym_id", mpos->second[1]));
                 key_value_mapping.insert(std::make_pair("applied_operations", mpos->second[2]));
                 chain_mapping_values.push_back(key_value_mapping);

                 key_value_mapping.clear();
                 key_value_mapping.insert(std::make_pair("entity_id", chain->entity_id()));
                 key_value_mapping.insert(std::make_pair("asym_id", chain->ChainID()));
                 if (chain->chain_type() == "ATOMS") {
                      key_value_mapping.insert(std::make_pair("pdb_asym_id", chain->PDB_ChainID()));
                      branch_scheme_mapping.insert(std::make_pair((int) branch_scheme_mapping.size(), std::make_pair(mpos->second[0], key_value_mapping)));
                 } else {
                      key_value_mapping.insert(std::make_pair("pdb_strand_id", chain->PDB_ChainID()));
                      if ((chain->chain_type() == "ATOMN") || (chain->chain_type() == "ATOMP")) {
                           int size = (int) polymer_scheme_mapping.size();
                           polymer_scheme_mapping.insert(std::make_pair(size, std::make_pair(mpos->second[0], key_value_mapping)));
                      } else {
                           int size = (int) nonpolymer_scheme_mapping.size();
                           nonpolymer_scheme_mapping.insert(std::make_pair(size, std::make_pair(mpos->second[0], key_value_mapping)));
                      }
                 }
            }
            chain = mol->GetNextChain();
       }

       _update_entity_id_mapping_category(biolBlock, new_old_entity_id_mapping);
       _update_entity_category(modelBlock, biolBlock, new_old_entity_id_mapping);
       _update_polymer_entity_related_categories(modelBlock, biolBlock, polymer_entity_id_mapping);
       _update_single_key_value_categories(modelBlock, biolBlock, "pdbx_entity_branch", "entity_id", branch_entity_id_mapping);
       _update_multiple_key_value_categories(modelBlock, biolBlock, "pdbx_entity_branch_descriptor", "entity_id", "ordinal", branch_entity_id_mapping);
       _update_multiple_key_value_categories(modelBlock, biolBlock, "pdbx_entity_branch_link", "entity_id", "link_id", branch_entity_id_mapping);
       _update_multiple_key_value_categories(modelBlock, biolBlock, "pdbx_entity_branch_list", "entity_id", "", branch_entity_id_mapping);
       _update_single_key_value_categories(modelBlock, biolBlock, "pdbx_entity_nonpoly", "entity_id", nonpolymer_entity_id_mapping);
       _update_multiple_key_value_categories(modelBlock, biolBlock, "pdbx_poly_seq_scheme", "asym_id", "", polymer_scheme_mapping);
       _update_multiple_key_value_categories(modelBlock, biolBlock, "pdbx_branch_scheme", "asym_id", "", branch_scheme_mapping);
       _update_multiple_key_value_categories(modelBlock, biolBlock, "pdbx_nonpoly_scheme", "asym_id", "", nonpolymer_scheme_mapping);
       _update_chain_id_mapping_category(biolBlock, chain_mapping_values);
}

void GenBioAssembly::_update_struct_conf_categories(Block& modelBlock, const std::map<int, std::vector<std::string> >& original_asym_chain_id_mapping,
                                                    const std::set<std::string>& chain_id_set, const std::vector<std::set<int> >& chain_index_list,
                                                    RCSB::Molecule* mol, Block& biolBlock)
{
       std::vector<std::string> item_names;
       std::map<std::string, std::vector<std::map<std::string, std::string> > > original_values;
       _get_original_struct_conf_values(modelBlock, chain_id_set, item_names, original_values);

       std::vector<std::map<std::string, std::string> > all_values;
       all_values.clear();
       if (!original_values.empty() && !chain_index_list.empty()) {
            std::map<int, std::map<std::string, std::map<std::string, std::string> > > order_key_value_map;
            std::map<std::string, std::map<std::string, std::string> > key_value_map;
            std::map<std::string, std::string> tmp_map;
            std::vector<std::string> data;
            for (std::vector<std::set<int> >::const_iterator vpos = chain_index_list.begin(); vpos != chain_index_list.end(); ++vpos) {
                 order_key_value_map.clear();
                 for (std::set<int>::const_iterator spos = vpos->begin(); spos != vpos->end(); ++spos) {
                      bool is_removed = false;
                      RCSB::Chain* chain = mol->GetIndexChain(*spos, is_removed);
                      if (!chain || (chain->chain_type() != "ATOMP")) continue;
                      std::map<int, std::vector<std::string> >::const_iterator mpos = original_asym_chain_id_mapping.find(*spos);
                      if (mpos == original_asym_chain_id_mapping.end()) continue;

                      tmp_map.clear();
                      tmp_map.insert(std::make_pair("beg_label_asym_id", chain->ChainID()));
                      tmp_map.insert(std::make_pair("end_label_asym_id", chain->ChainID()));
                      tmp_map.insert(std::make_pair("beg_auth_asym_id", chain->PDB_ChainID()));
                      tmp_map.insert(std::make_pair("end_auth_asym_id", chain->PDB_ChainID()));
                      key_value_map.clear();
                      key_value_map.insert(std::make_pair(mpos->second[1], tmp_map));
                      order_key_value_map.insert(std::make_pair(chain->order(), key_value_map)); 
                 }
                 if (order_key_value_map.empty()) continue;

                 for (std::map<int, std::map<std::string, std::map<std::string, std::string> > >::const_iterator
                      ompos = order_key_value_map.begin(); ompos != order_key_value_map.end(); ++ompos) {
                      for (std::map<std::string, std::map<std::string, std::string> >::const_iterator
                           kmpos = ompos->second.begin(); kmpos != ompos->second.end(); ++kmpos) {
                           std::map<std::string, std::vector<std::map<std::string, std::string> > >::iterator ovmpos = original_values.find(kmpos->first);
                           if (ovmpos == original_values.end()) continue;
                           for (std::vector<std::map<std::string, std::string> >::iterator ovpos = ovmpos->second.begin();
                                ovpos != ovmpos->second.end(); ++ovpos) {
                                for (std::map<std::string, std::string>::const_iterator svmpos = kmpos->second.begin(); svmpos != kmpos->second.end(); ++svmpos) {
                                     std::map<std::string, std::string>::iterator mpos = ovpos->find(svmpos->first);
                                     if (mpos != ovpos->end()) mpos->second = svmpos->second;
                                     else ovpos->insert(std::make_pair(svmpos->first, svmpos->second));
                                }
                                all_values.push_back(*ovpos);
                           }
                      } 
                 }
            }
       }

       if (!all_values.empty()) {
            std::string character_1, character_2;
            if (all_values.size() > 6084) {
                 character_1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";
                 character_2 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";
            } else {
                 character_1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
                 character_2 = "123456789";
            }
            std::set<std::string> id_set;
            id_set.clear();
            int count = 0;
            for (std::vector<std::map<std::string, std::string> >::iterator vpos = all_values.begin(); vpos != all_values.end(); ++vpos) {
                 count++;
                 std::map<std::string, std::string>::iterator mpos = vpos->find("id");
                 if (mpos != vpos->end()) mpos->second = "HELX_P" + String::IntToString(count);
                 else vpos->insert(std::make_pair("id", "HELX_P" + String::IntToString(count)));

                 std::string pdb_helix_id = get_next_id(character_1, character_1, character_2, id_set);
                 mpos = vpos->find("pdbx_PDB_helix_id");
                 if (mpos != vpos->end()) mpos->second = pdb_helix_id;
                 else vpos->insert(std::make_pair("pdbx_PDB_helix_id", pdb_helix_id)); 
            }
       }

       _update_values(biolBlock, "struct_conf", item_names, all_values);

       if (!all_values.empty()) {
            ISTable *t = getTableCopy(modelBlock, "struct_conf_type");
            if (t) biolBlock.WriteTable(t);
       } else deleteTable(biolBlock, "struct_conf_type");
}

void GenBioAssembly::_update_struct_sheet_categories(Block& modelBlock, const std::map<int, std::vector<std::string> >& original_asym_chain_id_mapping,
                                                     const std::set<std::string>& chain_id_set, const std::vector<std::set<int> >& chain_index_list,
                                                     RCSB::Molecule* mol, Block& biolBlock)
{
       std::vector<std::string> struct_sheet_items, sheet_order_items, sheet_range_items, sheet_hbond_items, check_item_names;
       // pair.first: sheet_id
       // pair.second: value list
       std::vector<std::pair<std::string, std::vector<std::map<std::string, std::string> > > > orig_struct_sheet_values, orig_sheet_order_values;
       std::vector<std::pair<std::string, std::vector<std::map<std::string, std::string> > > > orig_sheet_range_values, orig_sheet_hbond_values;

       check_item_names.clear();
       check_item_names.push_back("beg_auth_asym_id");
       check_item_names.push_back("end_auth_asym_id");
       _get_sheet_values(modelBlock, "struct_sheet_range", "sheet_id", check_item_names, chain_id_set, sheet_range_items, orig_sheet_range_values);

       std::set<std::string> sheet_id_set;
       sheet_id_set.clear();
       for (std::vector<std::pair<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator
            vpos = orig_sheet_range_values.begin(); vpos != orig_sheet_range_values.end(); ++vpos) {
            sheet_id_set.insert(vpos->first);
       }
       if (orig_sheet_range_values.size() != sheet_id_set.size()) sheet_id_set.clear();

       check_item_names.clear();
       check_item_names.push_back("id");
       _get_sheet_values(modelBlock, "struct_sheet", "id", check_item_names, sheet_id_set, struct_sheet_items, orig_struct_sheet_values);
       if (sheet_id_set.size() != orig_struct_sheet_values.size()) orig_struct_sheet_values.clear();

       check_item_names.clear();
       check_item_names.push_back("sheet_id");
       _get_sheet_values(modelBlock, "struct_sheet_order", "sheet_id", check_item_names,  sheet_id_set, sheet_order_items, orig_sheet_order_values);
       if (sheet_id_set.size() != orig_sheet_order_values.size()) orig_sheet_order_values.clear();

       if (orig_struct_sheet_values.empty() || orig_sheet_order_values.empty() || orig_sheet_range_values.empty() || chain_index_list.empty()) {
            deleteTable(biolBlock, "struct_sheet");
            deleteTable(biolBlock, "struct_sheet_order");
            deleteTable(biolBlock, "struct_sheet_range");
            deleteTable(biolBlock, "pdbx_struct_sheet_hbond");
            return;
       }

       _get_sheet_values(modelBlock, "pdbx_struct_sheet_hbond", "sheet_id", check_item_names,  sheet_id_set, sheet_hbond_items, orig_sheet_hbond_values);

       std::string character_1, character_2;
       if ((orig_struct_sheet_values.size() * chain_index_list.size()) > 6084) {
            character_1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";
            character_2 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";
       } else {
            character_1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
            character_2 = "123456789";
       }
       std::set<std::string> id_set;
       id_set.clear();

       std::vector<std::map<std::string, std::string> > all_struct_sheet_values, all_sheet_order_values, all_sheet_range_values, all_sheet_hbond_values;
       all_struct_sheet_values.clear();
       all_sheet_order_values.clear();
       all_sheet_range_values.clear();
       all_sheet_hbond_values.clear();

       std::map<std::string, std::string> sheet_id_map, empty_chain_id_item_map, range_chain_id_item_map, hbond_chain_id_item_map;
       empty_chain_id_item_map.clear();

       range_chain_id_item_map.clear();
       range_chain_id_item_map.insert(std::make_pair("beg_auth_asym_id", "beg_label_asym_id"));
       range_chain_id_item_map.insert(std::make_pair("end_auth_asym_id", "end_label_asym_id"));

       hbond_chain_id_item_map.clear();
       hbond_chain_id_item_map.insert(std::make_pair("range_1_auth_asym_id", "range_1_label_asym_id"));
       hbond_chain_id_item_map.insert(std::make_pair("range_2_auth_asym_id", "range_2_label_asym_id"));

       // key: original PDB chain ID
       // pair.first: current PDB chain ID
       // pair.second: current asym ID
       std::map<std::string, std::pair<std::string, std::string> > empty_chain_id_map, chain_id_map;
       empty_chain_id_map.clear();

       for (std::vector<std::set<int> >::const_iterator vpos = chain_index_list.begin(); vpos != chain_index_list.end(); ++vpos) {
            chain_id_map.clear();
            for (std::set<int>::const_iterator spos = vpos->begin(); spos != vpos->end(); ++spos) {
                 bool is_removed = false;
                 RCSB::Chain* chain = mol->GetIndexChain(*spos, is_removed);
                 if (!chain || (chain->chain_type() != "ATOMP")) continue;
                 std::map<int, std::vector<std::string> >::const_iterator mpos = original_asym_chain_id_mapping.find(*spos);
                 if (mpos == original_asym_chain_id_mapping.end()) continue;

                 chain_id_map.insert(std::make_pair(mpos->second[1], std::make_pair(chain->PDB_ChainID(), chain->ChainID())));
            }
            if (chain_id_map.empty()) continue;

            sheet_id_map.clear();
            for (std::vector<std::pair<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator
                 pvpos = orig_struct_sheet_values.begin(); pvpos != orig_struct_sheet_values.end(); ++pvpos) {
                 std::string sheet_id = get_next_id(character_1, character_1, character_2, id_set);
                 sheet_id_map.insert(std::make_pair(pvpos->first, sheet_id));

                 for (std::vector<std::map<std::string, std::string> >::const_iterator vspos = pvpos->second.begin(); vspos != pvpos->second.end(); ++vspos) {
                      std::map<std::string, std::string> valMap = *vspos;
                      std::map<std::string, std::string>::iterator mpos = valMap.find("id");
                      if (mpos == valMap.end()) continue;

                      mpos->second = sheet_id;
                      all_struct_sheet_values.push_back(valMap);
                 }
            }

            _get_updated_sheet_value(orig_sheet_order_values, sheet_id_map, empty_chain_id_item_map, empty_chain_id_map, all_sheet_order_values);
            _get_updated_sheet_value(orig_sheet_range_values, sheet_id_map, range_chain_id_item_map, chain_id_map, all_sheet_range_values);
            _get_updated_sheet_value(orig_sheet_hbond_values, sheet_id_map, hbond_chain_id_item_map, chain_id_map, all_sheet_hbond_values);
       }

       _update_values(biolBlock, "struct_sheet", struct_sheet_items, all_struct_sheet_values);
       _update_values(biolBlock, "struct_sheet_order", sheet_order_items, all_sheet_order_values);
       _update_values(biolBlock, "struct_sheet_range", sheet_range_items, all_sheet_range_values);
       _update_values(biolBlock, "pdbx_struct_sheet_hbond", sheet_hbond_items, all_sheet_hbond_values);
}

void GenBioAssembly::_update_entity_id_mapping_category(Block& block, const std::map<int, std::pair<std::string, int> >& entity_id_mapping)
{
       if (entity_id_mapping.empty()) {
            deleteTable(block, "pdbx_entity_remapping");
            return;
       }

       std::vector<std::string> item_names;
       item_names.clear();
       item_names.push_back("entity_id");
       item_names.push_back("orig_entity_id");

       std::map<std::string, std::string> key_value_mapping;

       std::vector<std::map<std::string, std::string> > values;
       values.clear();

       for (std::map<int, std::pair<std::string, int> >::const_iterator mpos = entity_id_mapping.begin(); mpos != entity_id_mapping.end(); ++mpos) {
            key_value_mapping.clear();
            key_value_mapping.insert(std::make_pair("entity_id", String::IntToString(mpos->first)));
            key_value_mapping.insert(std::make_pair("orig_entity_id", mpos->second.first));
            values.push_back(key_value_mapping);
       }
       _update_values(block, "pdbx_entity_remapping", item_names, values);
}

void GenBioAssembly::_update_entity_category(Block& modelBlock, Block& biolBlock, const std::map<int, std::pair<std::string, int> >& entity_id_mapping)
{
       std::map<std::string, std::string> tmp_map;
       std::map<int, std::pair<std::string, std::map<std::string, std::string> > > entity_id_update_mapping;
       entity_id_update_mapping.clear();
       for (std::map<int, std::pair<std::string, int> >::const_iterator empos = entity_id_mapping.begin(); empos != entity_id_mapping.end(); ++empos) {
            tmp_map.clear();
            tmp_map.insert(std::make_pair("id", String::IntToString(empos->first)));
            tmp_map.insert(std::make_pair("pdbx_number_of_molecules", String::IntToString(empos->second.second)));
            entity_id_update_mapping.insert(std::make_pair(empos->first, std::make_pair(empos->second.first, tmp_map)));
       }

       _update_single_key_value_categories(modelBlock, biolBlock, "entity", "id", entity_id_update_mapping);
}

void GenBioAssembly::_update_polymer_entity_related_categories(Block& modelBlock, Block& biolBlock, const std::map<int, std::pair<std::string,
                                                                 std::vector<std::string> > >& polymer_entity_id_mapping)
{
       std::map<std::string, std::string> tmp_map;
       std::map<int, std::pair<std::string, std::map<std::string, std::string> > > entity_id_update_mapping, entity_id_update_mapping1;
       entity_id_update_mapping.clear();
       entity_id_update_mapping1.clear();
       for (std::map<int, std::pair<std::string, std::vector<std::string> > >::const_iterator pmpos = polymer_entity_id_mapping.begin();
            pmpos != polymer_entity_id_mapping.end(); ++pmpos) {
            tmp_map.clear();
            tmp_map.insert(std::make_pair("entity_id", String::IntToString(pmpos->first)));
            entity_id_update_mapping1.insert(std::make_pair(pmpos->first, std::make_pair(pmpos->second.first, tmp_map)));
            tmp_map.insert(std::make_pair("pdbx_strand_id", join_string(pmpos->second.second, ",")));
            entity_id_update_mapping.insert(std::make_pair(pmpos->first, std::make_pair(pmpos->second.first, tmp_map)));
       }

       _update_single_key_value_categories(modelBlock, biolBlock, "entity_poly", "entity_id", entity_id_update_mapping);
       _update_multiple_key_value_categories(modelBlock, biolBlock, "entity_poly_seq", "entity_id", "", entity_id_update_mapping1);
}

void GenBioAssembly::_update_chain_id_mapping_category(Block& block, const std::vector<std::map<std::string, std::string> >& values)
{
       if (values.empty()) {
            deleteTable(block, "pdbx_chain_remapping");
            return;
       }

       std::vector<std::string> item_names;
       item_names.clear();
       item_names.push_back("entity_id");
       item_names.push_back("label_asym_id");
       item_names.push_back("auth_asym_id");
       item_names.push_back("orig_label_asym_id");
       item_names.push_back("orig_auth_asym_id");
       item_names.push_back("applied_operations");

       _update_values(block, "pdbx_chain_remapping", item_names, values);
}

void GenBioAssembly::_update_single_key_value_categories(Block& modelBlock, Block& biolBlock, const std::string& category, const std::string& key_item,
                                         const std::map<int, std::pair<std::string, std::map<std::string, std::string> > >& updated_key_value_mapping)
{
       if (updated_key_value_mapping.empty()) {
            deleteTable(biolBlock, category);
            return;
       }

       std::vector<std::string> item_names;
       std::map<std::string, std::map<std::string, std::string> > single_key_value_map;

       _get_single_key_value_map(modelBlock, category, key_item, item_names, single_key_value_map);
       if (single_key_value_map.empty()) {
            deleteTable(biolBlock, category);
            return;
       }

       std::vector<std::map<std::string, std::string> > values;
       values.clear();

       for (std::map<int, std::pair<std::string, std::map<std::string, std::string> > >::const_iterator empos = updated_key_value_mapping.begin();
            empos != updated_key_value_mapping.end(); ++empos) {
            std::map<std::string, std::map<std::string, std::string> >::iterator avmpos = single_key_value_map.find(empos->second.first);
            if (avmpos == single_key_value_map.end()) continue;

            for (std::map<std::string, std::string>::const_iterator svmpos = empos->second.second.begin(); svmpos != empos->second.second.end(); ++svmpos) {
                 std::map<std::string, std::string>::iterator mpos = avmpos->second.find(svmpos->first);
                 if (mpos != avmpos->second.end()) mpos->second = svmpos->second;
                 else avmpos->second.insert(std::make_pair(svmpos->first, svmpos->second));
            }

            values.push_back(avmpos->second);
       }

       _update_values(biolBlock, category, item_names, values);
}

void GenBioAssembly::_update_multiple_key_value_categories(Block& modelBlock, Block& biolBlock, const std::string& category, const std::string& key_item,
                                                           const std::string& ordinal_item, const std::map<int, std::pair<std::string,
                                                           std::map<std::string, std::string> > >& updated_key_value_mapping)
{
       if (updated_key_value_mapping.empty()) {
            deleteTable(biolBlock, category);
            return;
       }

       std::vector<std::string> item_names;
       std::map<std::string, std::vector<std::map<std::string, std::string> > > multiple_key_value_map;

       _get_multiple_key_value_map(modelBlock, category, key_item, item_names, multiple_key_value_map);
       if (multiple_key_value_map.empty()) {
            deleteTable(biolBlock, category);
            return;
       }

       std::vector<std::map<std::string, std::string> > values;
       values.clear();

       for (std::map<int, std::pair<std::string, std::map<std::string, std::string> > >::const_iterator empos = updated_key_value_mapping.begin();
            empos != updated_key_value_mapping.end(); ++empos) {
            std::map<std::string, std::vector<std::map<std::string, std::string> > >::iterator avmpos = multiple_key_value_map.find(empos->second.first);
            if (avmpos == multiple_key_value_map.end()) continue;

            for (std::vector<std::map<std::string, std::string> >::iterator vpos= avmpos->second.begin(); vpos != avmpos->second.end(); ++vpos) {
                 for (std::map<std::string, std::string>::const_iterator svmpos = empos->second.second.begin(); svmpos != empos->second.second.end(); ++svmpos) {
                      std::map<std::string, std::string>::iterator mpos = vpos->find(svmpos->first);
                      if (mpos != vpos->end()) mpos->second = svmpos->second;
                      else vpos->insert(std::make_pair(svmpos->first, svmpos->second));
                 }
                 if (!ordinal_item.empty()) {
                      std::map<std::string, std::string>::iterator mpos = vpos->find(ordinal_item);
                      if (mpos != vpos->end()) mpos->second = String::IntToString(values.size() + 1);
                      else vpos->insert(std::make_pair(ordinal_item, String::IntToString(values.size() + 1)));
                 }
                 values.push_back(*vpos);
            }
       }

       _update_values(biolBlock, category, item_names, values);
}

void GenBioAssembly::_get_single_key_value_map(Block& block, const std::string& category, const std::string& key_item, std::vector<std::string>& item_names,
                                               std::map<std::string, std::map<std::string, std::string> >& single_key_value_map)
{
       item_names.clear();
       single_key_value_map.clear();

       ISTable *t = getTablePtr(block, category);
       if (!t) return;

       std::vector<std::map<std::string, std::string> > all_values;
       get_values(t, all_values);
       if (all_values.empty()) return;

       item_names = t->GetColumnNames();
       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = all_values.begin(); pos != all_values.end(); ++pos) {
            std::map<std::string, std::string>::const_iterator mpos = pos->find(key_item);
            if (mpos != pos->end()) single_key_value_map.insert(std::make_pair(mpos->second, *pos));
       }
}

void GenBioAssembly::_get_multiple_key_value_map(Block& block, const std::string& category, const std::string& key_item, std::vector<std::string>& item_names,
                                                 std::map<std::string, std::vector<std::map<std::string, std::string> > >& multiple_key_value_map)
{
       item_names.clear();
       multiple_key_value_map.clear();

       ISTable *t = getTablePtr(block, category);
       if (!t) return;

       std::vector<std::map<std::string, std::string> > all_values, tmp_vec;
       get_values(t, all_values);
       if (all_values.empty()) return;

       item_names = t->GetColumnNames();
       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = all_values.begin(); pos != all_values.end(); ++pos) {
            std::map<std::string, std::string>::const_iterator mpos = pos->find(key_item);
            if (mpos == pos->end()) continue;

            std::map<std::string, std::vector<std::map<std::string, std::string> > >::iterator kmpos = multiple_key_value_map.find(mpos->second);
            if (kmpos != multiple_key_value_map.end()) kmpos->second.push_back(*pos);
            else {
                 tmp_vec.clear();
                 tmp_vec.push_back(*pos);
                 multiple_key_value_map.insert(std::make_pair(mpos->second, tmp_vec));
            }
       }
}

void GenBioAssembly::_get_original_struct_conf_values(Block& block, const std::set<std::string>& chain_id_set, std::vector<std::string>& item_names,
                                                      std::map<std::string, std::vector<std::map<std::string, std::string> > >& original_values)
{
       item_names.clear();
       original_values.clear();

       if (chain_id_set.empty()) return;

       std::map<std::string, std::vector<std::map<std::string, std::string> > > multiple_key_value_map;
       _get_multiple_key_value_map(block, "struct_conf", "beg_auth_asym_id", item_names, multiple_key_value_map);
       if (multiple_key_value_map.empty()) return;

       for (std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator
            mpos = multiple_key_value_map.begin(); mpos != multiple_key_value_map.end(); ++mpos) {
            if (chain_id_set.find(mpos->first) == chain_id_set.end()) continue;
            original_values.insert(std::make_pair(mpos->first, mpos->second));
       }
}

void GenBioAssembly::_get_sheet_values(Block& block, const std::string& category, const std::string& sheet_id_item, const std::vector<std::string>&
                                       check_item_names, const std::set<std::string>& check_value_set, std::vector<std::string>& category_item_names,
                                       std::vector<std::pair<std::string, std::vector<std::map<std::string, std::string> > > >& original_values)
{
       category_item_names.clear();
       original_values.clear();

       if (check_value_set.empty()) return;

       ISTable *t = getTablePtr(block, category);
       if (!t) return;

       std::vector<std::map<std::string, std::string> > all_values, selected_values;
       get_values(t, all_values);
       if (all_values.empty()) return;

       category_item_names = t->GetColumnNames();

       std::string prev_sheet_id = "123456789";
       selected_values.clear();

       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = all_values.begin(); pos != all_values.end(); ++pos) {
            std::map<std::string, std::string>::const_iterator mpos = pos->find(sheet_id_item);
            if (mpos == pos->end()) continue;

            if (mpos->second != prev_sheet_id) {
                 if (_is_valid_selected_values(check_item_names, check_value_set, selected_values)) {
                      original_values.push_back(std::make_pair(prev_sheet_id, selected_values));
                 }
                 prev_sheet_id = mpos->second;
                 selected_values.clear();
            }
            selected_values.push_back(*pos);
       }
       if (!selected_values.empty()) {
            if (_is_valid_selected_values(check_item_names, check_value_set, selected_values)) {
                 original_values.push_back(std::make_pair(prev_sheet_id, selected_values));
            }
       }
}

bool GenBioAssembly::_is_valid_selected_values(const std::vector<std::string>& check_item_names, const std::set<std::string>& check_value_set,
                                               const std::vector<std::map<std::string, std::string> >& selected_values)
{
       if (selected_values.empty()) return false;

       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = selected_values.begin(); pos != selected_values.end(); ++pos) {
            for (std::vector<std::string>::const_iterator ipos = check_item_names.begin(); ipos != check_item_names.end(); ++ipos) {
                 std::map<std::string, std::string>::const_iterator mpos = pos->find(*ipos);
                 if (mpos == pos->end()) return false;
                 if (check_value_set.find(mpos->second) == check_value_set.end()) return false;
            }
       }

       return true;
}

void GenBioAssembly::_get_updated_sheet_value(const std::vector<std::pair<std::string, std::vector<std::map<std::string, std::string> > > >& input_values,
                                 const std::map<std::string, std::string>& sheet_id_map, const std::map<std::string, std::string>& chain_id_item_map,
                                 const std::map<std::string, std::pair<std::string, std::string> >& chain_id_map,
                                 std::vector<std::map<std::string, std::string> >& output_values)
{
       if (input_values.empty() || sheet_id_map.empty()) return;

       for (std::vector<std::pair<std::string, std::vector<std::map<std::string, std::string> > > >::const_iterator
            vpos = input_values.begin(); vpos != input_values.end(); ++vpos) {
            std::map<std::string, std::string>::const_iterator idpos = sheet_id_map.find(vpos->first);
            if (idpos == sheet_id_map.end()) continue;

            for (std::vector<std::map<std::string, std::string> >::const_iterator vmpos = vpos->second.begin(); vmpos != vpos->second.end(); ++vmpos) {
                 std::map<std::string, std::string> valMap = *vmpos;
                 std::map<std::string, std::string>::iterator mpos = valMap.find("sheet_id");
                 if (mpos == valMap.end()) continue;

                 mpos->second = idpos->second;
                 if (chain_id_map.empty()) {
                      output_values.push_back(valMap);
                      continue;
                 }

                 bool found_match = true;
                 for (std::map<std::string, std::string>::const_iterator cmpos = chain_id_item_map.begin(); cmpos != chain_id_item_map.end(); ++cmpos) {
                      mpos = valMap.find(cmpos->first);
                      if (mpos == valMap.end()) {
                           found_match = false;
                           break;
                      }
                      std::map<std::string, std::pair<std::string, std::string> >::const_iterator cvmpos = chain_id_map.find(mpos->second);
                      if (cvmpos == chain_id_map.end()) {
                           found_match = false;
                           break;
                      }
                      mpos->second = cvmpos->second.first; // beg_auth_asym_id/end_auth_asym_id

                      mpos = valMap.find(cmpos->second);
                      if (mpos == valMap.end()) {
                           found_match = false;
                           break;
                      }
                      mpos->second = cvmpos->second.second; // beg_label_asym_id/end_label_asym_id
                 }
                 if (found_match) output_values.push_back(valMap);
            }
       }
}
