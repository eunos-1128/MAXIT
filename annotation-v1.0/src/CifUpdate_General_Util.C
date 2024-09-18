/*
FILE:     CifUpdate_General_Util.C
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
#include <sys/stat.h>
#include <math.h>

#include "AnnotationObj.h"
#include "CategoryMapping.h"
#include "V4V5Mapping.h"
#include "utillib.h"

#define NUM_EM_ROWS      4
#define NUM_EM_COLUMNS  10

static const char *_em_speciman_related_categories[NUM_EM_ROWS][NUM_EM_COLUMNS] = {
       { "embedding_applied",     "em_embedding",     "details",             "material",     "",         "",          "",           "",       "",     ""                    },
       { "shadowing_applied",     "em_shadowing",     "angle",               "details",      "material", "thickness", "",           "",       "",     ""                    },
       { "staining_applied",      "em_staining",      "details",             "material",     "type",     "",          "",           "",       "",     ""                    },
       { "vitrification_applied", "em_vitrification", "chamber_temperature", "cryogen_name", "details",  "humidity",  "instrument", "method", "temp", "time_resolved_state" }
};

void AnnotationObj::V4_V5_Cleanup()
{
       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       V4V5Mapping cleanupMaping;
       cleanupMaping.ReadMapping(_rcsbroot + "/data/binary/v4_v5_clean_up.odb");
       std::string error = cleanupMaping.Update(_cifblock, OTHER_MAPPING_TYPE);
       if (!error.empty()) _logIo->message(error.c_str());
}

void AnnotationObj::cif_update_EM_V4_V5()
{
       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       _cif_update_em_V4_V5_conversion(_cifblock);
       _cif_update_em_categories(_cifblock);
}

bool AnnotationObj::old_to_new_audit_revision_conversion(const std::string& depid, const std::string& version)
{
       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       std::string status_code = "";
       ISTable *t = getTablePtr(_cifblock, "pdbx_database_status");
       if (t) get_value_clean_upper(status_code, t, 0, "status_code");
       if (status_code != "REL" && status_code != "OBS") { 
            std::string error = "Status Code=" + status_code + "\n";
            // _logIo->message(error.c_str());
            return false;
       }

       std::string revision_date, revision_type;
       std::set<std::string> t_set;
       
       std::string initial_release_date = "";
       // key: release date
       // values: _pdbx_audit_revision_group.group
       std::map<std::string, std::set<std::string> > revision_type_data;
       revision_type_data.clear();

       t = getTablePtr(_cifblock, "pdbx_version");
       if (t) {
            std::map<std::string, std::string> revision_type_mapping = _read_old_new_version_mapping();

            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(revision_date, t, i, "revision_date");
                 if (revision_date.empty()) continue;
                 get_value_clean_upper(revision_type, t, i, "revision_type");
                 if (revision_type.empty()) continue;
                 if (revision_type == "INITIAL RELEASE") initial_release_date = revision_date;

                 std::map<std::string, std::string>::const_iterator mpos = revision_type_mapping.find(revision_type);
                 if (mpos == revision_type_mapping.end()) {
                      get_value_clean(revision_type, t, i, "revision_type");
                      std::string error = "Can't find group value for '" + revision_type + "'\n";
                      _logIo->message(error.c_str());
                      continue;
                 }

                 std::map<std::string, std::set<std::string> >::iterator rpos = revision_type_data.find(revision_date);
                 if (rpos != revision_type_data.end()) {
                      rpos->second.insert(mpos->second);
                 } else {
                      t_set.clear();
                      t_set.insert(mpos->second);
                      revision_type_data.insert(std::make_pair(revision_date, t_set));
                 }
            }
       }

       std::string obsolete_date = "";
       if (initial_release_date.empty() || status_code == "OBS") {
            std::string cs;
            t = getTablePtr(_cifblock, "database_PDB_rev");
            if (t) {
                 if (initial_release_date.empty()) {
                      std::string num, mod_type;
                      get_value_clean(num, t, 0, "num");
                      get_value_clean(mod_type, t, 0, "mod_type");
                      if (num == "1" && mod_type == "0") get_value_clean(initial_release_date, t, 0, "date");
                      if (revision_type_data.empty() && t->GetNumRows() > 1 && status_code == "REL") {
                           // _logIo->message("No pdbx_version category.\n");
                           return false;
                      }
                 }
                 if (status_code == "OBS") {
                      std::string rev_num = "";
                      ISTable *t1 = getTablePtr(_cifblock, "database_PDB_rev_record");
                      if (t1) {
                           int rowNo = t1->GetNumRows();
                           for (int i = 0; i < rowNo; ++i) {
                                get_value_clean_upper(cs, t1, i, "type");
                                if (cs == "OBSLTE") {
                                     get_value_clean(rev_num, t1, i, "rev_num");
                                     break;
                                }
                           }
                      }
                      int rowNo = t->GetNumRows();
                      for (int i = 0; i < rowNo; ++i) {
                           get_value_clean(cs, t, i, "num");
                           if (cs == rev_num) {
                                get_value_clean(obsolete_date, t, i, "date");
                                break;
                           }
                      }
                 }
            }
            if (obsolete_date.empty() && status_code == "OBS") {
                 t = getTablePtr(_cifblock, "pdbx_database_PDB_obs_spr");
                 int rowNo = t->GetNumRows();
                 for (int i = 0; i < rowNo; ++i) {
                      get_value_clean_upper(cs, t, i, "id");
                      if (cs == "OBSLTE") {
                           get_value_clean(obsolete_date, t, i, "date");
                           break;
                      }
                 }
            }
       }

       if (initial_release_date.empty() && status_code == "REL") {
            _logIo->message("Can't find initial release date.\n");
            return false;
       }
       if (obsolete_date.empty() && status_code == "OBS") {
            _logIo->message("Can't find obsolete date for obsolete entry.\n");
            return false;
       }

       if (!obsolete_date.empty()) {
            std::map<std::string, std::set<std::string> >::const_iterator mpos = revision_type_data.find(obsolete_date);
            if (mpos == revision_type_data.end()) {
                 t_set.clear();
                 revision_type_data.insert(std::make_pair(obsolete_date, t_set));
            }
       }

       std::map<int, _Audit_Revision> audit_revision_records;
       audit_revision_records.clear();

       _Audit_Revision audit_revision;
       audit_revision.major_revision = 1;
       audit_revision.minor_revision = (int) audit_revision_records.size();
       audit_revision.data_content_type = "Structure model";
       audit_revision.revision_date = initial_release_date;
       audit_revision.internal_deposition_id.clear();
       audit_revision.internal_version.clear();
       audit_revision.revision_group.clear();
       audit_revision.revision_details.clear();
       audit_revision.revision_category.clear();
       audit_revision.revision_item.clear();

       std::vector<std::string> data;
       data.clear();
       data.push_back("repository");
       data.push_back("Initial release");
       data.push_back("");
       data.push_back("");
       audit_revision.revision_details.push_back(data);

       int index = (int) audit_revision_records.size() + 1;
       audit_revision_records.insert(std::make_pair(index, audit_revision));

       for (std::map<std::string, std::set<std::string> >::const_iterator mpos = revision_type_data.begin(); mpos != revision_type_data.end(); ++mpos) {
            if (mpos->first <= initial_release_date) continue;

            audit_revision.minor_revision = (int) audit_revision_records.size();
            audit_revision.revision_date = mpos->first;
            audit_revision.revision_group.clear();
            audit_revision.revision_details.clear();
            if (mpos->first == obsolete_date) {
                 data.clear();
                 data.push_back("repository");
                 data.push_back("Obsolete");
                 data.push_back("");
                 data.push_back("");
                 audit_revision.revision_details.push_back(data);
            }
            for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                 audit_revision.revision_group.insert(*spos);
            }
            index = (int) audit_revision_records.size() + 1;
            audit_revision_records.insert(std::make_pair(index, audit_revision));
       }

       if (!depid.empty() && !version.empty()) {
            std::map<int, _Audit_Revision>::iterator mpos = audit_revision_records.find(index);
            mpos->second.internal_deposition_id = depid;
            mpos->second.internal_version = version;
       }

       _write_audit_revision_records(_cifblock, audit_revision_records);

       return true;
}

std::map<std::string, std::string> AnnotationObj::_read_old_new_version_mapping()
{
       std::map<std::string, std::string> revision_type_mapping;
       revision_type_mapping.clear();

       struct stat stattext;
       std::string mapping_support_file = _rcsbroot + "/data/binary/old_new_version_mapping.odb";
       if (stat(mapping_support_file.c_str(), &stattext) == 0) {
            CifFile* fobj = new CifFile(READ_MODE, mapping_support_file);
            if (fobj) {
                 Block &block = fobj->GetBlock(fobj->GetFirstBlockName());
                 ISTable* t = getTablePtr(block, "version_mapping");
                 if (t) {
                      std::string cs1, cs2;
                      int rowNo = t->GetNumRows();
                      for (int i = 0; i < rowNo; ++i) {
                           get_value_clean_upper(cs1, t, i, "old_version_type");
                           if (cs1.empty()) continue;
                           get_value_clean(cs2, t, i, "new_group");
                           if (cs2.empty()) continue;
                           revision_type_mapping.insert(std::make_pair(cs1, cs2));
                      }
                 }
                 delete fobj;
            }
       }
       return revision_type_mapping;
}

void AnnotationObj::_write_audit_revision_records(Block& block, const std::map<int, _Audit_Revision>& audit_revision_records)
{
       if (audit_revision_records.empty()) return;

       ISTable *history = _newTablePtr("pdbx_audit_revision_history");
       ISTable *details = NULL;
       ISTable *group = NULL;
       ISTable *category = NULL;
       ISTable *item = NULL;
       for (std::map<int, _Audit_Revision>::const_iterator mpos = audit_revision_records.begin(); mpos != audit_revision_records.end(); ++mpos) {
            int history_row = history->GetNumRows();
            history->AddRow();
            history->UpdateCell(history_row, "ordinal", String::IntToString(mpos->first));
            history->UpdateCell(history_row, "data_content_type", mpos->second.data_content_type);
            history->UpdateCell(history_row, "major_revision", String::IntToString(mpos->second.major_revision));
            history->UpdateCell(history_row, "minor_revision", String::IntToString(mpos->second.minor_revision));
            history->UpdateCell(history_row, "revision_date", mpos->second.revision_date);
            history->UpdateCell(history_row, "internal_deposition_id", mpos->second.internal_deposition_id);
            history->UpdateCell(history_row, "internal_version", mpos->second.internal_version);
            for (std::vector<std::vector<std::string> >::const_iterator vpos = mpos->second.revision_details.begin();
                 vpos != mpos->second.revision_details.end(); ++vpos) {
                 if (details == NULL) details = _newTablePtr("pdbx_audit_revision_details");
                 int details_row = details->GetNumRows();
                 details->AddRow();
                 details->UpdateCell(details_row, "ordinal", String::IntToString(details_row + 1));
                 details->UpdateCell(details_row, "revision_ordinal", String::IntToString(mpos->first));
                 details->UpdateCell(details_row, "data_content_type", mpos->second.data_content_type);
                 details->UpdateCell(details_row, "provider", (*vpos)[0]);
                 details->UpdateCell(details_row, "type", (*vpos)[1]);
                 details->UpdateCell(details_row, "description", (*vpos)[2]);
                 details->UpdateCell(details_row, "details", (*vpos)[3]);
            }
            for (std::set<std::string>::const_iterator spos = mpos->second.revision_group.begin(); spos != mpos->second.revision_group.end(); ++spos) {
                 if (group == NULL) group = _newTablePtr("pdbx_audit_revision_group");
                 int group_row = group->GetNumRows();
                 group->AddRow();
                 group->UpdateCell(group_row, "ordinal", String::IntToString(group_row + 1));
                 group->UpdateCell(group_row, "revision_ordinal", String::IntToString(mpos->first));
                 group->UpdateCell(group_row, "data_content_type", mpos->second.data_content_type);
                 group->UpdateCell(group_row, "group", *spos);
            }
            for (std::set<std::string>::const_iterator spos = mpos->second.revision_category.begin(); spos != mpos->second.revision_category.end(); ++spos) {
                 if (category == NULL) category = _newTablePtr("pdbx_audit_revision_category");
                 int category_row = category->GetNumRows();
                 category->AddRow();
                 category->UpdateCell(category_row, "ordinal", String::IntToString(category_row + 1));
                 category->UpdateCell(category_row, "revision_ordinal", String::IntToString(mpos->first));
                 category->UpdateCell(category_row, "data_content_type", mpos->second.data_content_type);
                 category->UpdateCell(category_row, "category", *spos);
            }
            for (std::set<std::string>::const_iterator spos = mpos->second.revision_item.begin(); spos != mpos->second.revision_item.end(); ++spos) {
                 if (item == NULL) item = _newTablePtr("pdbx_audit_revision_item");
                 int item_row = item->GetNumRows();
                 item->AddRow();
                 item->UpdateCell(item_row, "ordinal", String::IntToString(item_row + 1));
                 item->UpdateCell(item_row, "revision_ordinal", String::IntToString(mpos->first));
                 item->UpdateCell(item_row, "data_content_type", mpos->second.data_content_type);
                 item->UpdateCell(item_row, "item", *spos);
            }
       }
       block.WriteTable(history);
       if (details) block.WriteTable(details);
       if (group) block.WriteTable(group);
       if (category) block.WriteTable(category);
       if (item) block.WriteTable(item);
}

void AnnotationObj::_cif_update_xray_V4_V5_conversion(Block& block)
{
       V4V5Mapping xrayMaping;
       xrayMaping.ReadMapping(_rcsbroot + "/data/binary/xray_v4_v5_mapping.odb");
       std::string error = xrayMaping.Update(block, XRAY_MAPPING_TYPE);
       if (!error.empty()) _logIo->message(error.c_str());
}

void AnnotationObj::_cif_update_em_V4_V5_conversion(Block &block)
{
       V4V5Mapping emMaping;
       emMaping.ReadMapping(_rcsbroot + "/data/binary/em_v4_v5_mapping.odb");
       emMaping.Update(block, EM_MAPPING_TYPE);

       std::list<std::string> cat_list;
       CategoryMapping::GetCategoryNameListByType("EM", cat_list);

       for (std::list<std::string>::const_iterator lpos = cat_list.begin(); lpos != cat_list.end(); ++lpos) {
            ISTable *t = _getTablePtr(block, *lpos);
            if (!t) deleteTable(block, *lpos);
            else block.WriteTable(t);
       }
}

void AnnotationObj::_cif_update_nmr_V4_V5_conversion(Block &block)
{
       V4V5Mapping nmrMaping;
       nmrMaping.ReadMapping(_rcsbroot + "/data/binary/nmr_v4_v5_mapping.odb");
       std::string error = nmrMaping.Update(block, NMR_MAPPING_TYPE);
       if (!error.empty()) _logIo->message(error.c_str());
}

void AnnotationObj::_cif_update_nmr_ensemble(Block &block)
{
       ISTable *t = _getTablePtr(block, "pdbx_nmr_ensemble");
       if (!t) {
            t = _newTablePtr("pdbx_nmr_ensemble");
            t->AddRow();
            t->UpdateCell(0, "entry_id", _StructureId);
       }
       t->UpdateCell(0, "conformers_submitted_total_number", String::IntToString(_molecules.size()));
       block.WriteTable(t);
}

void AnnotationObj::_cif_update_em_categories(Block &block)
{
       std::string cs;

       std::vector<std::string> items, data;
       std::map<std::string, std::vector<std::string> > em_speciman_related_category_mapping;
       em_speciman_related_category_mapping.clear();
       for (int i = 0; i < NUM_EM_ROWS; ++i) {
            data.clear();
            for (int j = 2; j < NUM_EM_COLUMNS; ++j) {
                 if (strcmp(_em_speciman_related_categories[i][j], "")) { data.push_back(_em_speciman_related_categories[i][j]);
                   
                 }
            }
            em_speciman_related_category_mapping.insert(std::make_pair(_em_speciman_related_categories[i][1], data));
       }

       for (std::map<std::string, std::vector<std::string> >::const_iterator mpos = em_speciman_related_category_mapping.begin();
            mpos != em_speciman_related_category_mapping.end(); ++mpos) {
            ISTable *t = getTablePtr(block, mpos->first);
            if (!t) {
                 deleteTable(block, mpos->first);
                 continue;
            }

            bool found_value = false;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 for (std::vector<std::string>::const_iterator pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                      get_value_clean(cs, t, i, *pos);
                      if (!cs.empty()) {
                           found_value = true;
                           break;
                      }
                 }
                 if (found_value) break;
            }
            if (!found_value) deleteTable(block, mpos->first);
       }

       ISTable *t = getTablePtr(block, "em_3d_reconstruction");
       if (t) {
            std::string resolution;
            get_value_clean(resolution, t, 0, "resolution");
            if (!resolution.empty() && String::IsNumber(resolution)) {
                 t = getTablePtr(block, "refine");
                 if (t) {
                      get_value_clean(cs, t, 0, "ls_d_res_high");
                      bool changed = false;
                      if (!cs.empty()) {
                           if (!String::IsNumber(cs)) {
                                t->UpdateCell(0, "ls_d_res_high", resolution);
                                changed = true;
                           } else if (fabs(atof(cs.c_str()) - atof(resolution.c_str())) > 0.001) {
                                t->UpdateCell(0, "ls_d_res_high", resolution);
                                changed = true;
                           }
                      }
                      if (changed) {
                           get_value_clean(cs, t, 0, "ls_d_res_low");
                           if (!cs.empty()) t->UpdateCell(0, "ls_d_res_low", resolution);
                           block.WriteTable(t);
                      }
                 }
            }
       }

       t = _getTablePtr(block, "em_3d_fitting_list");
       if (!t) return;

       items.clear();
       items.push_back("id");
       items.push_back("type");
       items.push_back("source_name");
       items.push_back("accession_code");

       ISTable *t1 = getTablePtr(block, "pdbx_initial_refinement_model");
       if (!t1) {
            std::vector<std::vector<std::string> > starting_model_list;
            starting_model_list.clear();
            std::set<std::string> unique_accession_code_set, unique_data_set;
            unique_accession_code_set.clear();
            unique_data_set.clear();

            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 data.clear();
                 data.push_back("");
                 for (unsigned int j = 1; j < items.size(); ++j) {
                      get_value_clean(cs, t, i, items[j]);
                      data.push_back(cs);
                 }
                 if (data[1].empty() || data[2].empty()) continue;

                 String::UpperCase(data[3], cs);
                 std::string unique_key = data[1] + "_" + data[2] + "_" + data[3];
                 if (unique_data_set.find(unique_key) != unique_data_set.end()) continue;

                 unique_data_set.insert(unique_key);
                 if (!cs.empty()) unique_accession_code_set.insert(cs);
                 data[0] = String::IntToString(starting_model_list.size() + 1);
                 starting_model_list.push_back(data);
            }           

            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {     
                 get_value_clean_upper(cs, t, i, "pdb_entry_id");
                 if (!cs.empty() && (unique_accession_code_set.find(cs) == unique_accession_code_set.end())) {
                      unique_accession_code_set.insert(cs);
                      data.clear();
                      data.push_back(String::IntToString(starting_model_list.size() + 1));
                      data.push_back("experimental model");
                      data.push_back("PDB");
                      data.push_back(cs);
                      starting_model_list.push_back(data);
                 }
            }
            if (starting_model_list.empty()) return;

            t1 = add_new_table("pdbx_initial_refinement_model", items);
            for (unsigned int i = 0; i < starting_model_list.size(); ++i) {
                 t1->AddRow();
                 for (unsigned int j = 0; j < items.size(); ++j) {
                      t1->UpdateCell(i, items[j], starting_model_list[i][j]);
                 }
            }
            block.WriteTable(t1);
       }

       std::map<std::string, std::string> initial_model_id_mapping;
       initial_model_id_mapping.clear();

       for (unsigned int i = 0; i < t1->GetNumRows(); ++i) {
            data.clear();
            for (std::vector<std::string>::const_iterator vpos = items.begin(); vpos != items.end(); ++vpos) {
                 get_value_clean_upper(cs, t1, i, *vpos);
                 if (!cs.empty()) data.push_back(cs);
            }
            if (data.size() < items.size()) continue;

            initial_model_id_mapping.insert(std::make_pair(data[1] + "_" + data[2] + "_" + data[3], data[0]));
       }
       if (initial_model_id_mapping.empty()) return;

       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            data.clear();
            for (unsigned int j = 1; j < items.size(); ++j) {
                 get_value_clean_upper(cs, t, i, items[j]);
                 if (!cs.empty()) data.push_back(cs);
            }
            if (data.size() < 3) {
                 get_value_clean_upper(cs, t, i, "pdb_entry_id");
                 if (!cs.empty()) {
                      t->UpdateCell(i, "type", "experimental model");
                      t->UpdateCell(i, "source_name", "PDB");
                      t->UpdateCell(i, "accession_code", cs);
                      data.clear();
                      data.push_back("EXPERIMENTAL MODEL");
                      data.push_back("PDB");
                      data.push_back(cs);
                 } else continue;
            }

            std::map<std::string, std::string>::const_iterator mpos = initial_model_id_mapping.find(data[0] + "_" + data[1] + "_" + data[2]);
            if (mpos != initial_model_id_mapping.end()) t->UpdateCell(i, "initial_refinement_model_id", mpos->second);
       }
}

void AnnotationObj::_cif_update_category(ISTable *t)
{
       std::string cs, cs1;

       int rowNo = t->GetNumRows();

       if (t->GetName() == "diffrn") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "crystal_id");
                 if (cs.empty()) t->UpdateCell(i, "crystal_id", "1");
                 _remove_k(t, "ambient_temp", i);
            }
       } else if (t->GetName() == "diffrn_radiation_wavelength") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "wt");
                 if (cs.empty()) t->UpdateCell(i, "wt", "1.0");
            }
       } else if (t->GetName() == "em_imaging") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "date");
                 if ((cs.size() == 11 || cs.size() == 9) && cs[2] == '-' && cs[6] == '-') {
                      convert_time_format(10, cs);
                      t->UpdateCell(i, "date", cs);
                 }
            }
       } else if (t->GetName() == "exptl_crystal_grow") {
            for (int i = 0; i < rowNo; ++i) {
                 _remove_k(t, "temp", i);

                 get_value_clean(cs, t, i, "crystal_id");
                 if (cs.empty()) t->UpdateCell(i, "crystal_id", "1");
     
                 get_value_clean(cs, t, i, "pH");
                 get_value_clean(cs1, t, i, "pdbx_pH_range");
                 if (!cs.empty() && !String::IsNumber(cs) && cs1.empty()) {
                      t->UpdateCell(i, "pH", ".");
                      t->UpdateCell(i, "pdbx_pH_range", cs);
                 }
            }
       } else if (t->GetName() == "exptl_crystal_grow_comp") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "crystal_id");
                 if (cs.empty()) t->UpdateCell(i, "crystal_id", "1");
            }
       } else if (t->GetName() == "pdbx_nmr_exptl") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "conditions_id");
                 if (cs.empty()) t->UpdateCell(i, "conditions_id", "1");
                 get_value_clean(cs, t, i, "solution_id");
                 if (cs.empty()) t->UpdateCell(i, "solution_id", "1");
            }
       } else if (t->GetName() == "pdbx_nmr_exptl_sample_conditions") {
            for (int i = 0; i < rowNo; ++i) {
                 _remove_k(t, "temperature", i);

                 get_value_clean(cs, t, i, "pressure");
                 if (cs == "1 ATM") {
                      t->UpdateCell(i, "pressure", "1");
                      if (t->IsColumnPresent("pressure_units")) t->UpdateCell(i, "pressure_units", "ATM");
                 } else if (!cs.empty() && t->IsColumnPresent("pressure_units") && isxdigit(cs[0]) && isxdigit(cs[cs.size() - 1])) {
                      get_value_clean(cs, t, i, "pressure_units");
                      if (cs.empty()) t->UpdateCell(i, "pressure_units", "ATM");
                 }
            }
       } else if (t->GetName() == "pdbx_nmr_spectrometer") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs, t, i, "field_strength");
                 if (!cs.empty()) {
                      String::StripAndCompressWs(cs);
                      t->UpdateCell(i, "field_strength", cs);
                 }
            }
       } else if (t->GetName() == "phasing_MAD_set") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "clust_id");
                 if (cs.empty()) t->UpdateCell(i, "clust_id", "1");

                 get_value_clean(cs, t, i, "expt_id");
                 if (cs.empty()) t->UpdateCell(i, "expt_id", "1");
            }
       } else if (t->GetName() == "phasing_MIR_der") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "der_set_id");
                 if (cs.empty()) t->UpdateCell(i, "der_set_id", "1");

                 get_value_clean(cs, t, i, "native_set_id");
                 if (cs.empty()) t->UpdateCell(i, "native_set_id", ".");
            }
       } else if (t->GetName() == "refine" || t->GetName() == "reflns" /* || t->GetName() == "reflns_shell" */ ) {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "pdbx_diffrn_id");
                 if (!cs.empty()) continue;

                 get_value_clean_upper(cs, t, i, "pdbx_refine_id");
                 if (!cs.empty()) {
                      std::map<std::string, std::string>::const_iterator mpos = _scattering_type_mapping.find(cs);
                      if (mpos != _scattering_type_mapping.end()) {
                           t->UpdateCell(i, "pdbx_diffrn_id", mpos->second);
                           continue;
                      }
                 }

                 if (rowNo == 1 &&  _scattering_type_mapping.size() <= 1 && _experimental_method_count == 1 && !_default_diffrn_id.empty()) {
                      t->UpdateCell(i, "pdbx_diffrn_id", _default_diffrn_id);
                 }
            }
       } else if ((t->GetName() == "refine_analyze") || (t->GetName() == "refine_ls_shell")) {
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs, t, i, "pdbx_refine_id");
                 if (cs.empty() || cs == "1") t->UpdateCell(i, "pdbx_refine_id", _default_exp_type);
            }
       } else if (t->GetName() == "refine_B_iso") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_upper(cs, t, i, "treatment");
                 if (cs == "I")
                      t->UpdateCell(i, "treatment", "isotropic");
                 else if (cs == "A")
                      t->UpdateCell(i, "treatment", "anisotropic");
            }
       } else if (t->GetName() == "refine_ls_shell") {
            std::string text = "";
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs, t, i, "pdbx_total_number_of_bins_used");
                 if (!cs.empty()) {
                      text = String::IntToString(atoi(cs.c_str()));
                      t->UpdateCell(i, "pdbx_total_number_of_bins_used", text);
                 } else if (!text.empty()) t->UpdateCell(i, "pdbx_total_number_of_bins_used", text);

                 get_value(cs, t, i, "number_reflns_R_free");
                 if (!cs.empty()) t->UpdateCell(i, "number_reflns_R_free", String::IntToString(atoi(cs.c_str())));
            }
       } else if (t->GetName() == "refine_occupancy") {
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_upper(cs, t, i, "treatment");
                 if (cs == "OF")
                      t->UpdateCell(i, "treatment", "fix");
                 else if (cs == "OR")
                      t->UpdateCell(i, "treatment", "ref");
            }
       }
}

void AnnotationObj::_update_table_index(ISTable *source_t, const std::string& source_item, ISTable *target_t, const std::string& target_item)
{
       std::set<std::string> data_set;
       data_set.clear();

       std::string cs;

       int rowNo = source_t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs, source_t, i, source_item);
            if (!cs.empty()) data_set.insert(cs);
       }
       if (data_set.empty()) return;

       _update_table_index(data_set, target_t, target_item);
}

void AnnotationObj::_update_table_index(const std::set<std::string>& data_set, ISTable *t, const std::string& item)
{
       std::vector<std::string> data;
       data.clear();
       std::set<int> int_set;
       int_set.clear();

       bool is_number = true;
       for (std::set<std::string>::const_iterator
            spos = data_set.begin(); spos != data_set.end(); ++spos) {
            data.push_back(*spos);
            if (String::IsNumber(*spos) && spos->find(".") == std::string::npos)
                 int_set.insert(atoi(spos->c_str()));
            else is_number = false;
       }
       if (is_number) {
            data.clear();
            for (std::set<int>::const_iterator spos = int_set.begin();
                 spos != int_set.end(); ++spos) {
                 data.push_back(String::IntToString(*spos));
            }
       }

       std::string cs;

       int rowNo = t->GetNumRows();
       for (std::vector<std::string>::const_iterator
            vpos = data.begin(); vpos != data.end(); ++vpos) {
            bool found = false;
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, item);
                 if (cs == *vpos) {
                      found = true;
                      break;
                 } else if (cs.empty()) {
                      t->UpdateCell(i, item, *vpos);
                      found = true;
                      break;
                 }
            }
            if (found) continue;

            t->AddRow();
            t->UpdateCell(rowNo, item, *vpos);
            rowNo++;
       }
}

void AnnotationObj::_remove_k(ISTable *t, const std::string &item, const int& row)
{
       std::string cs;
       get_value_clean(cs, t, row, item);
       if (cs.empty()) return;
       std::string::size_type p = cs.find("K");
       if (p != std::string::npos) cs.erase(p);
       p = cs.find("k");
       if (p != std::string::npos) cs.erase(p);
       t->UpdateCell(row, item, cs);
}
