/*
FILE:     AnnotationObj_Entity_Util.C
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
#include "utillib.h"

#define NUM_SYN_TO_GEN    4
#define NUM_NAT_TO_GEN   19

static const char *_syn_to_gen[NUM_SYN_TO_GEN][2] = {
       { "organism_scientific",   "pdbx_gene_src_scientific_name"  },
       { "organism_common_name",  "gene_src_common_name"           },
       { "ncbi_taxonomy_id",      "pdbx_gene_src_ncbi_taxonomy_id" },
       { "details",               "pdbx_description"               }
};

static const char *_nat_to_gen[NUM_NAT_TO_GEN][2] = {
       { "common_name",              "gene_src_common_name"             },
       { "pdbx_organism_scientific", "pdbx_gene_src_scientific_name"    },
       { "pdbx_ncbi_taxonomy_id",    "pdbx_gene_src_ncbi_taxonomy_id"   },
       { "genus",                    "gene_src_genus"                   },
       { "species",                  "gene_src_species"                 },
       { "strain",                   "gene_src_strain",                 },
       { "tissue",                   "gene_src_tissue",                 },
       { "tissue_fraction",          "gene_src_tissue_fraction"         },
       { "pdbx_fragment",            "pdbx_gene_src_fragment"           },
       { "pdbx_variant",             "pdbx_gene_src_variant"            },
       { "pdbx_cell_line",           "pdbx_gene_src_cell_line"          },
       { "pdbx_atcc",                "pdbx_gene_src_atcc"               },
       { "pdbx_culture_collection",  "pdbx_gene_src_culture_collection" },
       { "pdbx_cellular_location",   "pdbx_gene_src_cellular_location"  },
       { "pdbx_organ",               "pdbx_gene_src_organ"              },
       { "pdbx_organelle",           "pdbx_gene_src_organelle"          },
       { "pdbx_cell",                "pdbx_gene_src_cell"               },
       { "pdbx_plasmid_details",     "plasmid_details"                  },
       { "details",                  "pdbx_description"                 }
};

void AnnotationObj::AssignEntityId()
{
       if (_molecules.empty()) return;

       _annotate_branch_polymers();
/*
       for (std::vector<RCSB::Molecule*>::const_iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            RCSB::Chain* chain = (*mpos)->GetFirstChain();
            while (chain) {
                 if (chain->chain_type() == "ATOMS") _get_branch_entity_links(chain);
                 chain = (*mpos)->GetNextChain();
            }
       }
*/
       std::map<std::string, std::string> old_new_entity_mapping;
       _assign_new_entity_ids(_new_chain_entity_id_mapping, old_new_entity_mapping);
       _updateEntityIDs(old_new_entity_mapping);
}

void AnnotationObj::_updateEntityIDs(const std::map<std::string, std::string>& old_new_entity_mapping)
{
       if (old_new_entity_mapping.empty()) return;
       if (!_CifObj) return;

       bool need_update = false;
       for (std::map<std::string, std::string>::const_iterator mpos = old_new_entity_mapping.begin(); mpos != old_new_entity_mapping.end(); ++mpos) {
            if (mpos->first != mpos->second) {
                 need_update = true;
                 break;
            }
       }
       if (!need_update) return;

       Block& block = _CifObj->GetBlock(_firstBlockName);
       
       std::vector<std::string> categories, data;
       categories.clear();
       categories.push_back("pdbx_entity_func_bind_mode");
       categories.push_back("pdbx_entity_poly_domain");
       categories.push_back("pdbx_struct_ref_seq_depositor_info");
       categories.push_back("struct_ref");

       std::string cs;

       for (std::vector<std::string>::const_iterator vpos = categories.begin(); vpos != categories.end(); ++vpos) {
            ISTable *t = getTablePtr(block, *vpos);
            if (!t) continue;
            if (!t->IsColumnPresent("entity_id")) continue;

            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "entity_id");
                 std::map<std::string, std::string>::const_iterator mpos = old_new_entity_mapping.find(cs);
                 if (mpos == old_new_entity_mapping.end()) continue;
                 t->UpdateCell(i, "entity_id", mpos->second);
            }
            block.WriteTable(t);
       }

       ISTable *t = getTablePtr(block, "em_entity_assembly");
       if (t) {
            std::set<std::string> data_set;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "entity_id_list");
                 if (cs.empty()) continue;

                 data.clear();
                 data_set.clear();
                 get_wordarray(categories, cs, ", \t\n");
                 for (std::vector<std::string>::const_iterator vpos = categories.begin(); vpos != categories.end(); ++vpos) {
                      std::map<std::string, std::string>::const_iterator mpos = old_new_entity_mapping.find(*vpos);
                      if (mpos == old_new_entity_mapping.end()) continue;
                      if (data_set.find(mpos->second) == data_set.end()) {
                           data_set.insert(mpos->second);
                           data.push_back(mpos->second);
                      }
                 }
                 if (data.empty()) continue;
                 cs.clear();
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      if (!cs.empty()) cs += ",";
                      cs += *vpos;
                 }
                 t->UpdateCell(i, "entity_id_list", cs);
            }
            block.WriteTable(t);
       }
}

void AnnotationObj::_updateStructRef(const std::map<std::string, std::string>& chain_entity_id_mapping)
{
       if (!_CifObj) return;

       Block& block = _CifObj->GetBlock(_firstBlockName);
       ISTable *ref = getTablePtr(block, "struct_ref");
       ISTable *seq = getTablePtr(block, "struct_ref_seq");
       if (!ref || !seq) return;

       // add missing entity_id item to struct_ref category
       if (!ref->IsColumnPresent("entity_id")) {
            std::vector<std::string> missing_items;
            missing_items.clear();
            missing_items.push_back("entity_id");
            check_missing_item(ref, missing_items);
       }

       // check if struct_ref.id & struct_ref.entity_id exist
       const std::vector<std::string>& ref_items = ref->GetColumnNames();
       int id_idx = -1;
       int entity_idx = -1;
       for (unsigned int i = 0; i < ref_items.size(); ++i) {
            if (ref_items[i] == "id") id_idx = i;
            if (ref_items[i] == "entity_id") entity_idx = i;
       }
       if (id_idx < 0 || entity_idx < 0) return;

       // check if struct_ref_seq.ref_id & struct_ref_seq.pdbx_strand_id exist
       const std::vector<std::string>& seq_items = seq->GetColumnNames();
       int ref_idx = -1;
       int chn_idx = -1;
       for (unsigned int i = 0; i < seq_items.size(); ++i) {
            if (seq_items[i] == "ref_id") ref_idx = i;
            if (seq_items[i] == "pdbx_strand_id") chn_idx = i;
       }
       if (ref_idx < 0 || chn_idx < 0) return;

       std::vector<std::vector<std::string> > ref_values, seq_values;

       // read struct_ref category
       get_values(ref, ref_values);
       if (ref_values.empty()) return;

       std::map<std::string, std::vector<std::string> > ref_values_mapping;
       ref_values_mapping.clear();
       for (std::vector<std::vector<std::string> >::const_iterator pos = ref_values.begin(); pos != ref_values.end(); ++pos) {
            ref_values_mapping.insert(std::make_pair((*pos)[id_idx], *pos));
       }

       // read struct_ref_seq category
       get_values(seq, seq_values);
       if (seq_values.empty()) return;

       std::vector<std::vector<std::string> > new_ref_values, new_seq_values;
       new_ref_values.clear();
       new_seq_values.clear();

       // get PDB chain IDs from struct_ref_seq.pdbx_strand_id
       std::set<std::string> chain_set;
       chain_set.clear();
       for (std::vector<std::vector<std::string> >::const_iterator pos = seq_values.begin(); pos != seq_values.end(); ++pos) {
            if ((*pos)[chn_idx].empty()) continue;
            chain_set.insert((*pos)[chn_idx]);
       }
       if (chain_set.empty()) return;

       // re-generate struct_ref & struct_ref_seq
       for (std::vector<std::vector<std::string> >::iterator pos = seq_values.begin(); pos != seq_values.end(); ++pos) {
            std::map<std::string, std::vector<std::string> >::const_iterator mpos = ref_values_mapping.find((*pos)[ref_idx]);
            if (mpos == ref_values_mapping.end()) continue;

            std::vector<std::string> data = mpos->second;
            if (data[entity_idx].empty()) {
                 std::map<std::string, std::string>::const_iterator epos = chain_entity_id_mapping.find((*pos)[chn_idx]);
                 if (epos != chain_entity_id_mapping.end()) data[entity_idx] = epos->second;
            }

            bool found = false;
            for (std::vector<std::vector<std::string> >::const_iterator npos = new_ref_values.begin(); npos != new_ref_values.end(); ++npos) {
                 bool mismatch = false;
                 for (unsigned int i = 0; i < data.size(); ++i) {
                      if ((int) i == id_idx || data[i] == (*npos)[i]) continue;

                      mismatch = true;
                      break;
                 }
                 if (mismatch) continue;

                 found = true;
                 (*pos)[ref_idx] = (*npos)[id_idx];
                 new_seq_values.push_back(*pos);
                 break;
            }
            if (found) continue;

            data[id_idx] = String::IntToString(new_ref_values.size() + 1);
            new_ref_values.push_back(data);

            (*pos)[ref_idx] = data[id_idx];
            new_seq_values.push_back(*pos);
       }

       if (new_ref_values.empty() || new_seq_values.empty()) return;

       ISTable *t = add_new_table("struct_ref", ref_items);
       for (unsigned int i = 0; i < new_ref_values.size(); ++i) {
            t->AddRow();
            for (unsigned int j = 0; j < ref_items.size(); ++j) {
                 t->UpdateCell(i, ref_items[j], new_ref_values[i][j]);
            }
       }
       block.WriteTable(t);

       t = add_new_table("struct_ref_seq", seq_items);
       for (unsigned int i = 0; i < new_seq_values.size(); ++i) {
            t->AddRow();
            for (unsigned int j = 0; j < seq_items.size(); ++j) {
                 t->UpdateCell(i, seq_items[j], new_seq_values[i][j]);
            }
       }
       block.WriteTable(t);
}

void AnnotationObj::_split_entity_info(const int& entity_id, const std::map<int, std::pair<std::vector<std::string>, std::vector<std::string> > >& entityids)
{
       if (entity_id == 0 || entityids.empty()) return;

       std::map<int, Entity>::const_iterator epos = _entities.find(entity_id);
       if (epos == _entities.end()) return;
       if (epos->second.empty()) return;

       std::vector<std::string> items;
       items.clear();
       items.push_back("pdbx_beg_seq_num");
       items.push_back("pdbx_end_seq_num");

       Entity entity;
       
       for (std::map<int, std::pair<std::vector<std::string>, std::vector<std::string> > >::const_iterator
            mpos = entityids.begin(); mpos != entityids.end(); ++mpos) {
            entity.clear();
            entity.setEntityID(String::IntToString(mpos->first));
            for (unsigned int i = 2; i < mpos->second.second.size(); ++i) entity.insertPDBChainID(mpos->second.second[i]);
            entity.setSeqs(mpos->second.first);
            entity.CopyMetaData(epos->second);
            std::vector<std::pair<std::string, std::map<std::string, std::string> > > source = epos->second.source();
            if (source.size() == 1) {
                 for (unsigned int i = 0; i < items.size(); ++i) {
                      std::map<std::string, std::string>::iterator spos = source[0].second.find(items[i]);
                      if (spos != source[0].second.end())
                           spos->second = mpos->second.second[i];
                      else source[0].second.insert(std::make_pair(items[i], mpos->second.second[i]));
                 }
                 entity.setSource(source);
            }
            _entities.insert(std::make_pair(mpos->first, entity));
       }
}

void AnnotationObj::_merge_entity_info(const int& entity_id, const std::string& pdb_chainid, const std::vector<std::string>& seqs,
                                       const std::vector<std::vector<int> >& entityids)
{
       if (entity_id == 0 || entityids.empty()) return;

       std::vector<std::string> items;
       items.clear();
       items.push_back("pdbx_beg_seq_num");
       items.push_back("pdbx_end_seq_num");

       Entity entity;
       entity.clear();
       for (std::vector<std::vector<int> >::const_iterator vpos = entityids.begin(); vpos != entityids.end(); ++vpos) {
            std::map<int, Entity>::const_iterator epos = _entities.find((*vpos)[0]);
            if (epos == _entities.end()) continue;

            entity.MergeMetaData(epos->second);

            Entity old_entity = epos->second;
            std::vector<std::pair<std::string, std::map<std::string, std::string> > > source = old_entity.source();
            if (source.size() == 1) {
                 for (unsigned int i = 0; i < items.size(); ++i) {
                      std::map<std::string, std::string>::iterator mpos = source[0].second.find(items[i]);
                      if (mpos != source[0].second.end())
                           mpos->second = String::IntToString((*vpos)[i + 1]);
                      else source[0].second.insert(std::make_pair(items[i], String::IntToString((*vpos)[i + 1])));
                 }
                 old_entity.setSource(source);
            }
            entity.MergeSource(old_entity);
       }

       if (entity.empty()) return;

       entity.setEntityID(String::IntToString(entity_id));
       entity.insertPDBChainID(pdb_chainid);
       entity.setSeqs(seqs);
       _entities.insert(std::make_pair(entity_id, entity));
}

void AnnotationObj::_get_entity_category_definition(std::vector<std::pair<std::string, std::pair<std::string, std::map<std::string, std::string> > > >&
                                                    entity_category)
{
       entity_category.clear();

       std::map<std::string, std::string> key_item_mapping;

       key_item_mapping.clear();
       entity_category.push_back(std::make_pair("entity", std::make_pair("id", key_item_mapping)));

       entity_category.push_back(std::make_pair("entity_keywords", std::make_pair("entity_id", key_item_mapping)));

       key_item_mapping.clear();
       key_item_mapping.insert(std::make_pair("name", "com_name"));
       entity_category.push_back(std::make_pair("entity_name_com", std::make_pair("entity_id", key_item_mapping)));

       key_item_mapping.clear();
       key_item_mapping.insert(std::make_pair("name", "sys_name"));
       entity_category.push_back(std::make_pair("entity_name_sys", std::make_pair("entity_id", key_item_mapping)));

       key_item_mapping.clear();
       key_item_mapping.insert(std::make_pair("type", "poly_type"));
       key_item_mapping.insert(std::make_pair("pdbx_target_identifier", "pdbx_target_identifier"));
       entity_category.push_back(std::make_pair("entity_poly", std::make_pair("entity_id", key_item_mapping)));
}

void AnnotationObj::_get_entity_source_category_definition(std::vector<std::pair<std::string, std::pair<std::string, std::map<std::string, std::string> > > >&
                                                           entity_category)
{
       entity_category.clear();

       std::map<std::string, std::string> key_item_mapping;

       key_item_mapping.clear();
       entity_category.push_back(std::make_pair("entity_src_gen", std::make_pair("man", key_item_mapping)));

       for (int i = 0; i < NUM_NAT_TO_GEN; ++i) {
            key_item_mapping.insert(std::make_pair(_nat_to_gen[i][0], _nat_to_gen[i][1]));
       }
       entity_category.push_back(std::make_pair("entity_src_nat", std::make_pair("nat", key_item_mapping)));

       key_item_mapping.clear();
       for (int i = 0; i < NUM_SYN_TO_GEN; ++i) {
            key_item_mapping.insert(std::make_pair(_syn_to_gen[i][0], _syn_to_gen[i][1]));
       }
       entity_category.push_back(std::make_pair("pdbx_entity_src_syn", std::make_pair("syn", key_item_mapping)));
}

void AnnotationObj::_read_entity_info(Block& block)
{
       // pair.first: category
       // pair.second.first: entity ID item name
       // pair.second.second: key/cif_item mapping
       std::vector<std::pair<std::string, std::pair<std::string, std::map<std::string, std::string> > > > entity_category;
       _get_entity_category_definition(entity_category);

       std::map<int, std::vector<std::map<std::string, std::string> > > info_mapping;
       for (std::vector<std::pair<std::string, std::pair<std::string, std::map<std::string, std::string> > > >::const_iterator
            pos = entity_category.begin(); pos != entity_category.end(); ++pos) {
            _read_info_mapping(block, pos->first, pos->second.first, pos->second.second, info_mapping);
            if (info_mapping.empty()) continue;

            for (std::map<int, std::vector<std::map<std::string, std::string> > >::const_iterator
                 ipos = info_mapping.begin(); ipos != info_mapping.end(); ++ipos) {
                 std::map<int, Entity>::iterator epos = _entities.find(ipos->first);
                 if (epos == _entities.end()) continue;
                 for (std::map<std::string, std::string>::const_iterator mpos = ipos->second[0].begin(); mpos != ipos->second[0].end(); ++mpos) {
                      epos->second.insertValue(mpos->first, mpos->second);
                 }
            }
       }

       _get_entity_source_category_definition(entity_category);

       for (std::vector<std::pair<std::string, std::pair<std::string, std::map<std::string, std::string> > > >::const_iterator
            pos = entity_category.begin(); pos != entity_category.end(); ++pos) {
            _read_info_mapping_1(block, pos->first, "entity_id", pos->second.second, info_mapping);
            if (info_mapping.empty()) continue;

            for (std::map<int, std::vector<std::map<std::string, std::string> > >::const_iterator
                 ipos = info_mapping.begin(); ipos != info_mapping.end(); ++ipos) {
                 std::map<int, Entity>::iterator epos = _entities.find(ipos->first);
                 if (epos == _entities.end()) continue;
                 for (unsigned int i = 0; i < ipos->second.size(); ++i) {
                      epos->second.insertSource(pos->second.first, ipos->second[i]);
                 }
            }
       }
/*
       std::map<std::string, std::string> key_item_mapping;
       key_item_mapping.clear();
       _read_info_mapping_1(block, "pdbx_entity_branch_descriptor", "entity_id", key_item_mapping, info_mapping, false);
       if (!info_mapping.empty()) {
            for (std::map<int, std::vector<std::map<std::string, std::string> > >::const_iterator
                 ipos = info_mapping.begin(); ipos != info_mapping.end(); ++ipos) {
                 std::map<int, Entity>::iterator epos = _entities.find(ipos->first);
                 if (epos == _entities.end()) continue;
                 epos->second.setLinearDescriptors(ipos->second);
            }
       }
*/
}

void AnnotationObj::_read_info_mapping(Block& block, const std::string& category, const std::string& item_id, const std::map<std::string, std::string>&
                                       item_mapping, std::map<int, std::vector<std::map<std::string, std::string> > >& info_mapping)
{
       info_mapping.clear();

       ISTable *t = getTablePtr(block, category);
       if (!t) return;

       std::string id, cs;
       std::map<std::string, std::string> tmp_map;
       std::vector<std::map<std::string, std::string> > tmp_vec;

       const std::vector<std::string>& item_names = t->GetColumnNames();

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(id, t, i, item_id);
            if (id.empty()) continue;

            tmp_map.clear();
            if (!item_mapping.empty()) {
                 for (std::map<std::string,std::string>::const_iterator mpos = item_mapping.begin(); mpos != item_mapping.end(); ++mpos) {
                      get_value(cs, t, i, mpos->first);
                      if (!cs.empty()) tmp_map.insert(std::make_pair(mpos->second, cs));
                 }
            } else {
                 for (std::vector<std::string>::const_iterator pos = item_names.begin(); pos != item_names.end(); ++pos) {
                      if (*pos == item_id) continue;
                      if ((category == "entity") && (*pos == "src_method")) {
                           get_value_clean_lower(cs, t, i, *pos);
                           if (cs == "recombinant") {
                                tmp_map.insert(std::make_pair(*pos, "man"));
                                continue;
                           } else if (cs == "synthetic") {
                                tmp_map.insert(std::make_pair(*pos, "syn"));
                                continue;
                           } else if ((cs != "man") && (cs != "nat") && (cs != "syn"))
                                continue; 
                      }
                      get_value(cs, t, i, *pos);
                      if (!cs.empty()) tmp_map.insert(std::make_pair(*pos, cs));
                 }
            }
            if (tmp_map.empty()) continue;

            int int_id = atoi(id.c_str());
            std::map<int, std::vector<std::map<std::string, std::string> > >::iterator mpos = info_mapping.find(int_id);
            if (mpos != info_mapping.end()) mpos->second.push_back(tmp_map);
            else {
                 tmp_vec.clear();
                 tmp_vec.push_back(tmp_map);
                 info_mapping.insert(std::make_pair(int_id, tmp_vec));
            }
       }
}

void AnnotationObj::_read_info_mapping_1(Block& block, const std::string& category, const std::string& item_id, const std::map<std::string, std::string>&
                                   item_mapping, std::map<int, std::vector<std::map<std::string, std::string> > >& info_mapping, const bool& source_flag)
{
       info_mapping.clear();

       ISTable *t = getTablePtr(block, category);
       if (!t) return;

       std::string id, cs;
       std::map<std::string, std::string> tmp_map;
       std::vector<std::map<std::string, std::string> > tmp_vec;

       const std::vector<std::string>& item_names = t->GetColumnNames();

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(id, t, i, item_id);
            if (id.empty()) continue;

            tmp_map.clear();
            for (std::vector<std::string>::const_iterator pos = item_names.begin(); pos != item_names.end(); ++pos) {
                 if ((*pos == item_id) || (*pos == "ordinal")) continue;
                 get_value(cs, t, i, *pos);
                 if (cs.empty()) continue;
                 std::string key = *pos;
                 std::map<std::string, std::string>::const_iterator mpos = item_mapping.find(*pos);
                 if (mpos != item_mapping.end()) key = mpos->second;
                 tmp_map.insert(std::make_pair(key, cs));
            }
            if (tmp_map.empty()) continue;

            if (source_flag && (tmp_map.find("pdbx_alt_source_flag") == tmp_map.end())) tmp_map.insert(std::make_pair("pdbx_alt_source_flag", "sample"));

            int int_id = atoi(id.c_str());
            std::map<int, std::vector<std::map<std::string, std::string> > >::iterator mpos = info_mapping.find(int_id);
            if (mpos != info_mapping.end()) mpos->second.push_back(tmp_map);
            else {
                 tmp_vec.clear();
                 tmp_vec.push_back(tmp_map);
                 info_mapping.insert(std::make_pair(int_id, tmp_vec));
            }
       }
}

void AnnotationObj::_write_entity_info(Block& block)
{
       if (_molecules.empty()) return;

       std::map<std::string, std::string> key_item_mapping;

       key_item_mapping.clear();
       _write_entity_category(block, "entity", "id", key_item_mapping);

       _write_entity_category(block, "entity_keywords", "entity_id", key_item_mapping);

       key_item_mapping.clear();
       key_item_mapping.insert(std::make_pair("name", "com_name"));
       _write_entity_category(block, "entity_name_com", "entity_id", key_item_mapping);

       key_item_mapping.clear();
       key_item_mapping.insert(std::make_pair("name", "sys_name"));
       _write_entity_category(block, "entity_name_sys", "entity_id", key_item_mapping);

       _write_entity_nonpoly_and_branch_categores(block);

       key_item_mapping.clear();
       for (int i = 0; i < NUM_NAT_TO_GEN; ++i) {
            key_item_mapping.insert(std::make_pair(_nat_to_gen[i][0], _nat_to_gen[i][1]));
       }
       _write_entity_source_category(block, "entity_src_nat", "nat", key_item_mapping);
       key_item_mapping.clear();
       _write_entity_source_category(block, "entity_src_gen", "man", key_item_mapping);
       key_item_mapping.clear();
       for (int i = 0; i < NUM_SYN_TO_GEN; ++i) {
            key_item_mapping.insert(std::make_pair(_syn_to_gen[i][0], _syn_to_gen[i][1]));
       }
       _write_entity_source_category(block, "pdbx_entity_src_syn", "syn", key_item_mapping);

       _write_pdbx_entity_branch_descriptor(block);
       _write_pdbx_entity_branch_link(block);
}

void AnnotationObj::_write_entity_category(Block& block, const std::string& category, const std::string& id,
                                           const std::map<std::string, std::string>& item_mapping)
{
       ISTable *t = _newTablePtr(category);
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       int row = 0;
       for (std::map<int, Entity>::const_iterator epos = _entities.begin(); epos != _entities.end(); ++epos) {
            bool found = false;
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 if (*pos == id) continue;
                 else {
                      std::string key = *pos;
                      std::map<std::string, std::string>::const_iterator mpos = item_mapping.find(*pos);
                      if (mpos != item_mapping.end()) key = mpos->second;
                      std::string val = epos->second.getValue(key);
                      if (!val.empty()) {
                           if (!found) t->AddRow();
                           found = true;
                           t->UpdateCell(row, *pos, val);
                      }
                 }
            }
            if (found) {
                 t->UpdateCell(row, id, String::IntToString(epos->first));
                 row++;
            }
       }

       if (row) block.WriteTable(t);
       else {
            delete t;
            deleteTable(block, category);
       }
}

void AnnotationObj::_write_entity_nonpoly_and_branch_categores(Block& block)
{
       ISTable *t1 = _newTablePtr("pdbx_entity_branch");
       ISTable *t2 = _newTablePtr("pdbx_entity_nonpoly");

       int row1 = 0;
       int row2 = 0;
       for (std::map<int, Entity>::const_iterator epos = _entities.begin(); epos != _entities.end(); ++epos) {
            std::string val = epos->second.getValue("type");
            if (val == "polymer") continue;
            if (val == "branched") {
                 t1->AddRow();
                 t1->UpdateCell(row1, "entity_id", String::IntToString(epos->first));
                 t1->UpdateCell(row1, "type", "oligosaccharide");
                 row1++;
            } else {
                 t2->AddRow();
                 t2->UpdateCell(row2, "entity_id", String::IntToString(epos->first));
                 t2->UpdateCell(row2, "name", epos->second.getValue("pdbx_description"));
                 t2->UpdateCell(row2, "comp_id", epos->second.getFirstResName());
                 row2++;
            }
       }

       if (row1) block.WriteTable(t1);
       else {
            delete t1;
            deleteTable(block, "pdbx_entity_branch");
       }
       if (row2) block.WriteTable(t2);
       else {
            delete t2;
            deleteTable(block, "pdbx_entity_nonpoly");
       }
}

void AnnotationObj::_write_entity_source_category(Block& block, const std::string& category, const std::string& type,
                                                  const std::map<std::string, std::string>& item_mapping)
{
       ISTable *t = _newTablePtr(category);
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       int row = 0;
       for (std::map<int, Entity>::const_iterator epos = _entities.begin(); epos != _entities.end(); ++epos) {
            const std::vector<std::pair<std::string, std::map<std::string, std::string> > >& source_info = epos->second.source();
            int count = 0;
            for (std::vector<std::pair<std::string, std::map<std::string, std::string> > >::const_iterator
                 pos = source_info.begin(); pos != source_info.end(); ++pos) {
                 if (pos->first != type) continue;

                 count++;

                 t->AddRow();
                 for (std::vector<std::string>::const_iterator vpos = itemNames.begin(); vpos != itemNames.end(); ++vpos) {
                      if (*vpos == "entity_id")
                           t->UpdateCell(row, "entity_id", String::IntToString(epos->first));
                      else {
                           std::string key = *vpos;
                           std::map<std::string, std::string>::const_iterator mpos = item_mapping.find(*vpos);
                           if (mpos != item_mapping.end()) key = mpos->second;

                           mpos = pos->second.find(key);
                           if (mpos != pos->second.end())
                                t->UpdateCell(row, *vpos, mpos->second);
                           else if (*vpos == "pdbx_src_id")
                                t->UpdateCell(row, *vpos, String::IntToString(count));
                      }
                 }
                 row++;
            }
       }

       if (row) block.WriteTable(t);
       else {
            delete t;
            deleteTable(block, category);
       }
}

void AnnotationObj::_write_pdbx_entity_branch_descriptor(Block& block)
{
       ISTable *t = _newTablePtr("pdbx_entity_branch_descriptor");
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       int row = 0;
       for (std::map<int, Entity>::const_iterator epos = _entities.begin(); epos != _entities.end(); ++epos) {
            const std::vector<std::map<std::string, std::string> >& linear_descriptors = epos->second.linear_descriptors();
            if (linear_descriptors.empty()) continue;

            for (std::vector<std::map<std::string, std::string> >::const_iterator pos = linear_descriptors.begin(); pos != linear_descriptors.end(); ++pos) {
                 t->AddRow();
                 for (std::vector<std::string>::const_iterator vpos = itemNames.begin(); vpos != itemNames.end(); ++vpos) {
                      if (*vpos == "ordinal") t->UpdateCell(row, "ordinal", String::IntToString(row + 1));
                      if (*vpos == "entity_id") t->UpdateCell(row, "entity_id", String::IntToString(epos->first));
                      else {
                           std::map<std::string, std::string>::const_iterator mpos = pos->find(*vpos);
                           if (mpos != pos->end()) t->UpdateCell(row, *vpos, mpos->second);
                      }
                 }
                 row++;
            }
       }

       if (row) block.WriteTable(t);
       else {
            delete t;
            deleteTable(block, "pdbx_entity_branch_descriptor");
       }
}

void AnnotationObj::_write_pdbx_entity_branch_link(Block& block)
{
       ISTable *t = _newTablePtr("pdbx_entity_branch_link");
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       int row = 0;
       for (std::map<int, Entity>::const_iterator epos = _entities.begin(); epos != _entities.end(); ++epos) {
            const std::vector<std::map<std::string, std::string> >& entity_links = epos->second.entity_links();
            if (entity_links.empty()) continue;

            for (std::vector<std::map<std::string, std::string> >::const_iterator pos = entity_links.begin(); pos != entity_links.end(); ++pos) {
                 t->AddRow();
                 for (std::vector<std::string>::const_iterator vpos = itemNames.begin(); vpos != itemNames.end(); ++vpos) {
                      if (*vpos == "link_id") t->UpdateCell(row, "link_id", String::IntToString(row + 1));
                      if (*vpos == "entity_id") t->UpdateCell(row, "entity_id", String::IntToString(epos->first));
                      else {
                           std::map<std::string, std::string>::const_iterator mpos = pos->find(*vpos);
                           if (mpos != pos->end()) t->UpdateCell(row, *vpos, mpos->second);
                      }
                 }
                 row++;
            }
       }

       if (row) block.WriteTable(t);
       else {
            delete t;
            deleteTable(block, "pdbx_entity_branch_link");
       }
}

void AnnotationObj::_assign_new_entity_ids(std::map<std::string, std::string>& new_chain_entity_id_mapping,
                                           std::map<std::string, std::string>& old_new_entity_mapping)
{
       if (_molecules.empty()) return;

       _entity_ids.clear();
       new_chain_entity_id_mapping.clear();
       old_new_entity_mapping.clear();

       std::map<int, Entity> new_entities;
       new_entities.clear();

       Entity entity;
       std::vector<std::string> seqs;

       std::string descriptor, cs;

       std::vector<RCSB::Residue*> residue_list;

       int entity_id = 0;
       for (int i = 0; i < 2; ++i) {
            bool firstModel = true;
            for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
                 if (i == 0) (*mpos)->InternalOrder(_numbering_based_order_flag);
                 RCSB::Chain* chain = (*mpos)->GetFirstChain();
                 while (chain) {
                      if (i == 0) {
                           if (chain->entity_key() == "HOH" || chain->entity_key() == "DOD") {
                                chain = (*mpos)->GetNextChain();
                                continue;
                           }
                      } else {
                           if (chain->entity_key() != "HOH" && chain->entity_key() != "DOD") {
                                chain = (*mpos)->GetNextChain();
                                continue;
                           }
                      }
                      std::map<std::string, int>::const_iterator epos = _entity_ids.find(chain->entity_key());
                      if (epos != _entity_ids.end()) {
                           std::map<int, Entity>::const_iterator opos = _entities.find(chain->int_prev_entity_id());
                           chain->set_entity_id(epos->second);
                           std::map<int, Entity>::iterator ipos = new_entities.find(epos->second);
                           if (ipos != new_entities.end() && firstModel) {
                                if (opos != _entities.end()) {
                                     ipos->second.MergeMetaData(opos->second, false);
                                     ipos->second.MergeSource_1(opos->second);
                                }
                                ipos->second.insertPDBChainID(chain->PDB_ChainID());
                                if (chain->entity_key() == "HOH" || chain->entity_key() == "DOD")
                                     ipos->second.addNumber(chain->ResidueNumbers());
                                else ipos->second.addNumber(1);
                           }
                           old_new_entity_mapping.insert(std::make_pair(chain->prev_entity_id(), String::IntToString(epos->second)));
                           if (chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMS")
                                new_chain_entity_id_mapping.insert(std::make_pair(chain->PDB_ChainID(), String::IntToString(epos->second)));
/*
                      } else if (chain->chain_type() == "ATOMS") {
                           chain->GetFirstResidueList(residue_list);
                           while (!residue_list.empty()) {
                                for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                                     epos = _entity_ids.find((*rpos)->ResName());
                                     if (epos != _entity_ids.end()) {
                                          (*rpos)->set_entity_id(String::IntToString(epos->second));
                                          std::map<int, Entity>::iterator ipos = new_entities.find(epos->second);
                                          if (ipos != new_entities.end()) {
                                               ipos->second.insertPDBChainID(chain->PDB_ChainID());
                                               ipos->second.addNumber(1);
                                          }
                                     } else {
                                          entity_id++;
                                          (*rpos)->set_entity_id(String::IntToString(entity_id));
                                          _entity_ids.insert(std::make_pair((*rpos)->ResName(), entity_id));

                                          entity.clear();
                                          entity.setEntityID(String::IntToString(entity_id));
                                          entity.setChainType(chain->chain_type());
                                          entity.insertPDBChainID(chain->PDB_ChainID());
                                          seqs.clear();
                                          seqs.push_back((*rpos)->ResName());
                                          entity.setSeqs(seqs);
                                          entity.addNumber(1);

                                          float weight = _ccDic->get_molecule_weight((*rpos)->ResName());
                                          if (weight > 0) entity.insertValue("formula_weight", FloatToString(weight, 0, 3));

                                          try {
                                                const ConnectFormat& drug = _ccDic->find_drug((*rpos)->ResName());
                                                entity.insertValue("pdbx_description", drug.chemical_name());
                                          } catch (const std::exception& exc) {}

                                          new_entities.insert(std::make_pair(entity_id, entity));
                                          _entity_ids.insert(std::make_pair((*rpos)->ResName(), entity_id));
                                     }
                                }
                                chain->GetNextResidueList(residue_list);
                           }
*/
                      } else {
                           entity_id++;
                           chain->set_entity_id(entity_id);
                           _entity_ids.insert(std::make_pair(chain->entity_key(), entity_id));

                           chain->get_seq(seqs);

                           entity.clear();
                           entity.setEntityID(String::IntToString(entity_id));
                           entity.setChainType(chain->chain_type());
                           entity.insertPDBChainID(chain->PDB_ChainID());
                           entity.setSeqs(seqs);
                           if (seqs.size() == 1 && (seqs[0] == "HOH" || seqs[0] == "DOD"))
                                entity.addNumber(chain->ResidueNumbers());
                           else entity.addNumber(1);

                           std::map<int, Entity>::const_iterator ipos = _entities.find(chain->int_prev_entity_id());
                           if (ipos != _entities.end()) entity.CopyAllMetaData(ipos->second);

                           std::string poly_type = entity.getValue("poly_type");
                           if (poly_type.empty()) poly_type = chain->get_poly_type();
                           if (!poly_type.empty()) entity.insertValue("poly_type", poly_type);

                           if (chain->weight() > 0) entity.insertValue("formula_weight", FloatToString(chain->weight(), 0, 3));

                           chain->get_descriptor(descriptor);
                           if (!descriptor.empty()) {
                                cs.clear();
                                if (seqs.size() == 1 && (seqs[0] == "HOH" || seqs[0] == "DOD")) descriptor = WATER_TEXT;
                                /* if (chain->chain_type() == "ATOMS")
                                     cs = "SUGAR (" + descriptor + ")";
                                else */ if (chain->chain_type() == "ATOMN") {
                                     if (chain->na_type() == ATOMN_TYPE_RNA_ONLY)
                                          cs = RNA_TEXT;
                                     else if (chain->na_type() == ATOMN_TYPE_TRNA)
                                          cs = T_RNA_TEXT;
                                     else if (chain->na_type() == ATOMN_TYPE_DNA_RNA)
                                          cs = HYBRID_TEXT;
                                     else cs = DNA_TEXT;
                                     cs += " (" + descriptor + ")";
                                }
                                if (!cs.empty()) descriptor = cs;
                                if ((seqs.size() == 1) || (chain->chain_type() == "ATOMS"))
                                     entity.insertValue("pdbx_description", descriptor, true);
                                else entity.insertValue("pdbx_description", descriptor);
                           }

                           const std::vector<std::map<std::string, std::string> >& linear_descriptors = chain->get_linear_descriptors();
                           if (!linear_descriptors.empty()) entity.setLinearDescriptors(linear_descriptors);

                           const std::vector<std::map<std::string, std::string> >& branch_links = chain->get_branch_links();
                           if (!branch_links.empty()) entity.setEntityLinks(branch_links);

                           new_entities.insert(std::make_pair(entity_id, entity));
                           _entity_ids.insert(std::make_pair(chain->entity_key(), entity_id));

                           old_new_entity_mapping.insert(std::make_pair(chain->prev_entity_id(), String::IntToString(entity_id)));
                           if (chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP")
                                new_chain_entity_id_mapping.insert(std::make_pair(chain->PDB_ChainID(), String::IntToString(entity_id)));
                      }
                      chain = (*mpos)->GetNextChain();
                 }
                 firstModel = false;
            }
       }

       _entities = new_entities;
       for (std::map<int, Entity>::iterator pos = _entities.begin(); pos != _entities.end(); ++pos) {
            pos->second.update();
       }
}

void AnnotationObj::_get_branch_entity_links(RCSB::Chain* chain, const std::list<_LINK>& input_links,
                                    const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& bad_links)
{
       if (input_links.empty()) return;

       std::set<std::string> bad_link_set;
       bad_link_set.clear();
       for (std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >::const_iterator lpos = bad_links.begin(); lpos != bad_links.end(); ++lpos) {
            bad_link_set.insert(lpos->first->getAtomAllIndex(false) + "-" + lpos->second->getAtomAllIndex(false));
            bad_link_set.insert(lpos->second->getAtomAllIndex(false) + "-" + lpos->first->getAtomAllIndex(false));
       }

       std::multimap<int, std::multimap<int, std::map<std::string, std::string> > > link_maps;
       link_maps.clear();
       std::set<std::string> link_indices;
       link_indices.clear();

       for (std::list<_LINK>::const_iterator pos = input_links.begin(); pos != input_links.end(); ++pos) {
            if ((!pos->SymOP_1.empty() && pos->SymOP_1 != "1_555") || (!pos->SymOP_2.empty() && pos->SymOP_2 != "1_555")) continue;
            // if ((pos->fstAtom->pdb_chnid() != chain->PDB_ChainID()) || (pos->sndAtom->pdb_chnid() != chain->PDB_ChainID())) continue;

            if (bad_link_set.find(pos->fstAtom->getAtomAllIndex(false) + "-" + pos->sndAtom->getAtomAllIndex(false)) != bad_link_set.end()) continue;

            std::pair<int, RCSB::Residue*> index_pair_1 = chain->find_pdb_residue_with_sort_index(pos->fstAtom->pdb_chnid(), pos->fstAtom->pdb_resnam(),
                                                                  pos->fstAtom->pdb_resnum(), pos->fstAtom->ins_code());
            if (index_pair_1.first < 0) continue;
            RCSB::Atom* atom_1 = index_pair_1.second->find_atom(pos->fstAtom->pdb_atmnam(), pos->fstAtom->alt_loc());
            if (atom_1 == NULL) continue;

            std::pair<int, RCSB::Residue*> index_pair_2 = chain->find_pdb_residue_with_sort_index(pos->sndAtom->pdb_chnid(), pos->sndAtom->pdb_resnam(),
                                                                  pos->sndAtom->pdb_resnum(), pos->sndAtom->ins_code());
            if (index_pair_2.first < 0) continue;
            RCSB::Atom* atom_2 = index_pair_2.second->find_atom(pos->sndAtom->pdb_atmnam(), pos->sndAtom->alt_loc());
            if (atom_2 == NULL) continue;

            if (index_pair_1.first == index_pair_2.first) continue;
            // else if (index_pair_1.first > index_pair_2.first)

            if (atom_1->atom_type() == atom_2->atom_type()) continue;
            if (atom_2->atom_type() == "C")
                 _insert_link_maps(pos->bondtype, index_pair_2.first, atom_2, index_pair_1.first, atom_1, link_maps, link_indices);
            else _insert_link_maps(pos->bondtype, index_pair_1.first, atom_1, index_pair_2.first, atom_2, link_maps, link_indices);
       }

       if (link_maps.empty()) return;

       std::vector<std::map<std::string, std::string> > branch_links;
       branch_links.clear();
       for (std::multimap<int, std::multimap<int, std::map<std::string, std::string> > >::const_iterator
            mpos = link_maps.begin(); mpos != link_maps.end(); ++mpos) {
            for (std::multimap<int, std::map<std::string, std::string> >::const_iterator mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                 branch_links.push_back(mmpos->second);
            }
       }

       chain->set_branch_links(branch_links);
}

void AnnotationObj::_insert_link_maps(const std::string& bondtype, const int& first_index, RCSB::Atom* first_atom, const int& second_index, RCSB::Atom* second_atom,
                           std::multimap<int, std::multimap<int, std::map<std::string, std::string> > >& link_maps, std::set<std::string>& link_indices)
{
       std::string cs = first_atom->getAtomAllIndex(false) + "_" + second_atom->getAtomAllIndex(false);
       if (link_indices.find(cs) != link_indices.end()) return;

       link_indices.insert(cs);

       std::string value_order = "sing";
       if (!bondtype.empty()) {
            value_order = bondtype;
            String::LowerCase(value_order);
       }

       std::map<std::string, std::string> key_val_map;
       key_val_map.clear();
       key_val_map.insert(std::make_pair("value_order", value_order));
       key_val_map.insert(std::make_pair("entity_branch_list_num_1", String::IntToString(first_index + 1)));
       key_val_map.insert(std::make_pair("comp_id_1", first_atom->pdb_resnam()));
       key_val_map.insert(std::make_pair("atom_id_1", first_atom->pdb_atmnam()));
       key_val_map.insert(std::make_pair("leaving_atom_id_1", _get_leaving_atom(first_atom->pdb_resnam(), first_atom->pdb_atmnam())));
       // key_val_map.insert(std::make_pair("atom_stereo_config_1", _get_stereo_config(first_atom->pdb_resnam(), first_atom->pdb_atmnam())));
       key_val_map.insert(std::make_pair("entity_branch_list_num_2", String::IntToString(second_index + 1)));
       key_val_map.insert(std::make_pair("comp_id_2", second_atom->pdb_resnam()));
       key_val_map.insert(std::make_pair("atom_id_2", second_atom->pdb_atmnam()));
       key_val_map.insert(std::make_pair("leaving_atom_id_2", _get_leaving_atom(second_atom->pdb_resnam(), second_atom->pdb_atmnam())));
       // key_val_map.insert(std::make_pair("atom_stereo_config_2", _get_stereo_config(second_atom->pdb_resnam(), second_atom->pdb_atmnam())));

       std::multimap<int, std::multimap<int, std::map<std::string, std::string> > >::iterator mpos = link_maps.find(first_index);
       if (mpos != link_maps.end())
            mpos->second.insert(std::make_pair(second_index, key_val_map));
       else {
            std::multimap<int, std::map<std::string, std::string> > tmp_map;
            tmp_map.clear();
            tmp_map.insert(std::make_pair(second_index, key_val_map));
            link_maps.insert(std::make_pair(first_index, tmp_map));
       }
}

std::string AnnotationObj::_get_leaving_atom(const std::string& compId, const std::string& atomId)
{
       try {
            const ConnectFormat& drug = _ccDic->find_drug(compId);
            return drug.getLeavingdAtom(atomId);
       } catch (const std::exception& exc) {}

       std::string error = "Can not find leaving atom for atom '" + atomId + "' in '" + compId + "' definition.\n";
       _logIo->message(error.c_str()); 

       return "";
}

std::string AnnotationObj::_get_stereo_config(const std::string& compId, const std::string& atomId)
{
       try {
            const ConnectFormat& drug = _ccDic->find_drug(compId);
            const AtomFormat& atom = drug.find_atom(atomId);
            return atom.stereo_config();
       } catch (const std::exception& exc) {}

       return "";
}
