/*
FILE:     Cif2Ndb.C
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

#include "CategoryMapping.h"
#include "Maxit.h"
#include "NdbToken.h"
#include "RmsRecord.h"
#include "SgCenter.h"
#include "utillib.h"

// static void ndb_name_conversion_first_last(std::string &str);

void Maxit::cif_to_ndb()
{
       if (!_CifObj) return;

       bool is_refmac5 = false;
       bool is_phenix = false;
       bool is_cns_xplor = false;
       bool is_buster = false;
       _get_refinement_program(is_refmac5, is_phenix, is_cns_xplor, is_buster);

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       _cif_update_deposit_process_sites(_cifblock);

       _cif_to_ndb_update_exptl_crystal_grow(_cifblock);

       _cif_to_ndb_update_struct_ref_seq_and_diff(_cifblock);

       _cif_to_ndb_mapping(_cifblock, is_refmac5 || is_phenix);

       _cif_to_ndb_post_processing(is_cns_xplor);
}

void Maxit::_cif_to_ndb_update_exptl_crystal_grow(Block& block)
{
       ISTable *t = _getTablePtr(block, "exptl_crystal_grow");
       if (!t) return;

       std::string details, U_details, pH, temp, method, cs;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(details, t, i, "pdbx_details");
            get_value_clean_upper(U_details, t, i, "pdbx_details");
            get_value_clean(pH, t, i, "pH");
            get_value_clean(temp, t, i, "temp");
            get_value_clean(method, t, i, "method");
      
            if (!pH.empty()) {
                 std::string cs1 = "PH " + pH;
                 std::string cs2 = "PH" + pH;
                 if (U_details.find(cs1) == std::string::npos && U_details.find(cs2) == std::string::npos) {
                      if (!details.empty()) details += ", ";
                      details += "pH " + pH;
                 }
            }

            if (!method.empty()) {
                 String::UpperCase(method, cs);
                 if (U_details.find(cs) == std::string::npos) {
                      if (!details.empty()) details += ", ";
                      details += method;
                 }
            }

            if (!temp.empty() && details.find(temp) == std::string::npos) {
                 if (!details.empty()) details += ", ";
                 details += "temperature " + temp + "K";
            }

            t->UpdateCell(i, "pdbx_details", details);
       }
       block.WriteTable(t);
}

void Maxit::_cif_to_ndb_update_struct_ref_seq_and_diff(Block& block)
{
       std::string cs;

       ISTable *t = _getTablePtr(block, "struct_ref_seq");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "pdbx_auth_seq_align_beg");
                 if (!cs.empty()) t->UpdateCell(i, "seq_align_beg", cs);
                 get_value_clean(cs, t, i, "pdbx_auth_seq_align_end");
                 if (!cs.empty()) t->UpdateCell(i, "seq_align_end", cs);
            }
            block.WriteTable(t);
       }

       t = _getTablePtr(block, "struct_ref_seq_dif");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_upper(cs, t, i, "details");
                 if (cs == "SEE _PDBX_ENTRY_DETAILS.SEQUENCE_DETAILS") t->UpdateCell(i, "details", "SEE REMARK 999");

                 if (cs == "DELETION") continue;
                 get_value_clean(cs, t, i, "pdbx_auth_seq_num");
                 if (!cs.empty()) t->UpdateCell(i, "seq_num", cs);
            }
            block.WriteTable(t);
       }
}

void Maxit::_cif_to_ndb_mapping(Block& block, const bool& refmac_phenix_flag)
{
       ISTable *Struct_Ref = NULL;

       static const std::map<std::string, CIF_CATEGORY>& CifCategory =
                           CategoryMapping::CifCategory();
       for (std::map<std::string, CIF_CATEGORY>::const_iterator
            cpos = CifCategory.begin(); cpos != CifCategory.end(); ++cpos) {
            if (!cpos->second.transfer) continue;

            ISTable *t = getTablePtr(block, cpos->second.category);
            if (!t) continue;

            if (cpos->second.category == "database_PDB_remark")
                 _cif_to_ndb_process_remark(t, "PREMRK", 2);
            else if (cpos->second.category == "database_PDB_rev")
                 _cif_to_ndb_process_database_pdb_rev(block, t);
            else if (cpos->second.category == "pdbx_database_related")
                 _cif_to_ndb_process_rcsb_database_related(t);
            else if (cpos->second.category == "pdbx_database_remark")
                 _cif_to_ndb_process_remark(t, "PDBRMK", 3);
            else if (cpos->second.category == "pdbx_pdb_compnd" ||
                     cpos->second.category == "pdbx_pdb_source")
                 _cif_to_ndb_process_general_category(t, cpos->second, false);
            else if (cpos->second.category == "refine_ls_restr") 
                 _cif_to_ndb_process_refine_ls_restr(t);
            else if (cpos->second.category == "refine_ls_restr_ncs") {
                 if (refmac_phenix_flag)
                      _cif_to_ndb_process_refmac_refine_ls_restr_ncs(t);
                 else _cif_to_ndb_process_category(t, cpos->second, false);
            } else if (cpos->second.category == "struct_ref") {
                 Struct_Ref = t;
                 continue;
/*
            } else if (cpos->second.category == "struct_site_gen")  {
                 if (!__need_file_version_update) cifndb_process_struct_site_gen(Table);
*/
            } else if (cpos->second.general)
                 _cif_to_ndb_process_general_category(t, cpos->second, true);
            else if (cpos->second.xray_or_nmr & 4)
                 _cif_to_ndb_process_matrix(t, cpos->second);
            else _cif_to_ndb_process_category(t, cpos->second, false);
       }

       if (Struct_Ref) {
            const CIF_CATEGORY& category = CategoryMapping::find_category("struct_ref");
            _cif_to_ndb_process_category(Struct_Ref, category, true);
       }
}

void Maxit::_cif_to_ndb_post_processing(const bool& is_cns_xplor)
{
       // reorder citation records
       _orderRecordWithIntegerKey("AUTH", 1);
       _orderRecordWithIntegerKey("TITL", 1);
       _orderRecordWithIntegerKey("EDIT", 1);
       _orderRecordWithIntegerKey("REF", 1);
       _orderRecordWithIntegerKey("PUBL", 1);
       _orderRecordWithIntegerKey("REFN", 1);

       // Strip and compress white spaces in TITLE record
       _cleanString("TITLE", 2);

       // copy initial deposition date to STATUS record
       _copyValue("HEADER", 2, "STATUS", 6, false);
/*
       // update process site
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("PROCES");
       if (ppos == _pdb_records.end()) _addNewRecord("PROCES");
       _updateRecordFront("PROCES", 5, "RCSB", false);
*/
       // reorder RADIAT records
       _orderRecordWithIntegerKey("RADIAT", 1);

       // Copy R_Value_Obs to R_Value_Work and vs versa for CNS & X-PLOR program
       if (is_cns_xplor) _updateRValueForCNS_XPLOR();

       // remove empty records
       _removeEmptyRecords();
}

void Maxit::_cif_to_ndb_process_remark(ISTable *t, const std::string& tokenid,
                                        const int& field_no)
{
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat(tokenid);
       int width = ndbformat.FieldList[field_no].FieldWidth;

       std::string id, text;
       std::vector<std::string> data;

       std::string old_id = "-1";
       bool first = true;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(id, t, i, "id");
            get_value(text, t, i, "text");
            if (text.empty()) continue; 
 
            if (id != old_id) {
                 old_id = id;
                 first = true;
            } else first = false;

            get_text_array_from_block(data, text, width);

            for (std::vector<std::string>::const_iterator
                 pos = data.begin(); pos != data.end(); ++pos) {
                 if (first && pos->empty()) continue;
                 _addNewRecord(tokenid);
                 _updateRecordBack(tokenid, 1, id);
                 _updateRecordBack(tokenid, field_no, *pos); 
            }
       }
}

void Maxit::_cif_to_ndb_process_database_pdb_rev(Block& block, ISTable *t)
{
       std::map<int, std::vector<std::string> > pdb_rev;
       pdb_rev.clear();

       std::vector<std::string> vector_data;
       std::string cs;

       int rowNo = t->GetNumRows(); 
       for (int i = 0; i < rowNo; ++i) {
            vector_data.clear();
            get_value_clean(cs, t, i, "date");
            vector_data.push_back(cs);
            get_value_clean(cs, t, i, "replaces");
            vector_data.push_back(cs);
            get_value_clean(cs, t, i, "mod_type");
            vector_data.push_back(cs);

            get_value_clean(cs, t, i, "num");
            int num = atoi(cs.c_str());
            pdb_rev.insert(std::make_pair(num, vector_data));

            get_value_clean(cs, t, i, "date_original");
            if (!cs.empty()) _updateRecordFront("HEADER", 2, cs);
       }

       ISTable *t1 = getTablePtr(block, "database_PDB_rev_record");
       if (t1) {
            std::set<std::string> tmp_set;

            std::map<int, std::set<std::string> > pdb_record;
            pdb_record.clear();

            rowNo = t1->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t1, i, "rev_num");
                 if (cs.empty()) continue;;

                 int num = atoi(cs.c_str());

                 get_value_clean_upper(cs, t1, i, "type");
                 if (cs.empty()) continue;;

                 std::map<int, std::vector<std::string> >::iterator rev_pos
                           = pdb_rev.find(num);
                 if (rev_pos == pdb_rev.end()) continue;

                 std::map<int, std::set<std::string> >::iterator record_pos
                           = pdb_record.find(num);
                 if (record_pos == pdb_record.end()) {
                      tmp_set.clear();
                      tmp_set.insert(cs);
                      pdb_record.insert(std::make_pair(num, tmp_set));
                 } else {
                      if (record_pos->second.find(cs) != record_pos->second.end())
                           continue;
                      record_pos->second.insert(cs);
                 }
                 rev_pos->second.push_back(cs);
            }
       }

       if (pdb_rev.empty()) return;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REVDAT");

       for (std::map<int, std::vector<std::string> >::reverse_iterator
            rev_pos = pdb_rev.rbegin(); rev_pos != pdb_rev.rend(); ++rev_pos) {
            if (rev_pos->first == 0) continue;
            _addNewRecord("REVDAT");
            _updateRecordBack("REVDAT", 1, String::IntToString(rev_pos->first));
            _updateRecordBack("REVDAT", 3, rev_pos->second[0]);
            _updateRecordBack("REVDAT", 4, rev_pos->second[1]);
            _updateRecordBack("REVDAT", 5, rev_pos->second[2]);
            int serial_no = 2;
            int pos = 6;
            for (unsigned int i = 3; i < rev_pos->second.size(); ++i) {
                 if (pos >= ndbformat.NumField) {
                      _addNewRecord("REVDAT");
                      _updateRecordBack("REVDAT", 1, String::IntToString(rev_pos->first));
                      _updateRecordBack("REVDAT", 2, String::IntToString(serial_no));
                      _updateRecordBack("REVDAT", 5, rev_pos->second[2]);
                      serial_no++;
                      pos = 6;
                 }
                 _updateRecordBack("REVDAT", pos, rev_pos->second[i]);
                 pos++;
            }
       }
}

void Maxit::_cif_to_ndb_process_rcsb_database_related(ISTable *t)
{
       std::set<std::string> split_ids;
       split_ids.clear();

       std::string cs, cs1, cs2;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs1, t, i, "db_id");
            if (cs1.empty()) continue;

            get_value_clean_lower(cs, t, i, "content_type");
            if (cs == "split") {
                 split_ids.insert(cs1);
                 continue;
            } else if (cs == "re-refinement") {
                 _original_entry_ids.push_back(cs1);
            }
            get_value_clean(cs, t, i, "db_name");
            get_value(cs2, t, i, "details");
            if (cs.empty() && cs1.empty() && cs2.empty()) continue;

            _addNewRecord("RENTRY");
            _updateRecordBack("RENTRY", 2, cs);
            _updateRecordBack("RENTRY", 3, cs1);
            _updateRecordBack("RENTRY", 4, cs2);
       }

       if (split_ids.empty()) return;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("SPLIT");

       _addNewRecord("SPLIT");
       int pos = 2;
       int serial_no = 2;
       for (std::set<std::string>::const_iterator
            spos = split_ids.begin(); spos != split_ids.end(); ++spos) { 
            if (pos >= ndbformat.NumField) {
                 _addNewRecord("SPLIT");
                 _updateRecordBack("SPLIT", 1, String::IntToString(serial_no));
                 serial_no++;
                 pos = 2;
            }
            _updateRecordBack("SPLIT", pos, *spos);
            pos++;
       }
}

void Maxit::_cif_to_ndb_process_general_category(ISTable* t, const CIF_CATEGORY &Category, 
                                                 const bool& clean_string_flag)
{
       std::string tokenid = "";
       std::map<int, std::string> field_item_mapping;
       field_item_mapping.clear();

       for (std::vector<CIF_KEYWORD>::const_iterator
            kpos = Category.keywords.begin(); kpos != Category.keywords.end(); ++kpos) {
            if (kpos->NdbInfo[0].NdbTokenName.empty() ||
                kpos->NdbInfo[0].NdbTokenName == "HEADER") continue;
            int ndb_fieldno = kpos->NdbInfo[0].NdbFieldNo - 1;
            if (ndb_fieldno < 0) continue;
            field_item_mapping.insert(std::make_pair(ndb_fieldno, kpos->name));
            tokenid = kpos->NdbInfo[0].NdbTokenName;
       }
       if (tokenid.empty() || field_item_mapping.empty()) return;

       _cif_to_ndb_process_single_token(t, tokenid, field_item_mapping, clean_string_flag);
}

void Maxit::_cif_to_ndb_process_single_token(ISTable* t, const std::string& tokenid, const
                                             std::map<int, std::string>& field_item_mapping,
                                             const bool& clean_string_flag)
{
       int rowNo = t->GetNumRows();

       std::string cs;
       for (int i = 0; i < rowNo; i++) {
            _addNewRecord(tokenid);
            for (std::map<int, std::string>::const_iterator mpos =
                 field_item_mapping.begin(); mpos != field_item_mapping.end(); ++mpos) {
                 get_value(cs, t, i, mpos->second);
                 if (cs.empty()) continue;
                 if (clean_string_flag) String::StripAndCompressWs(cs);
                 if (cs.empty()) continue;
                 _updateRecordBack(tokenid, mpos->first, cs);
            }
       }
}

void Maxit::_cif_to_ndb_process_refine_ls_restr(ISTable * t)
{
       if (!t->IsColumnPresent("type")) return;

       RmsRecord rmsrecord;
       if (!rmsrecord.Read(*_logIo, _rcsbroot)) return;

       std::string type, refine_code, cs;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(type, t, i, "type");
            if (type.empty()) continue;

            std::map<std::string, std::pair<std::string, int> > item_mapping
                = rmsrecord.findItemTokenMapping(type);
            if (item_mapping.empty()) continue;

            get_value_clean(refine_code, t, i, "pdbx_refine_id");

            for (std::map<std::string, std::pair<std::string, int> >::const_iterator
                 mpos = item_mapping.begin(); mpos != item_mapping.end(); ++mpos) {
                 get_value_clean(cs, t, i, mpos->first);
                 if (cs.empty()) continue;

                 const ndb_token_format& ndbformat =
                                NdbToken::getTokenFormat(mpos->second.first);

                 bool update_flag = false;
                 std::map<std::string, std::list<std::vector<std::string> > >::iterator
                     ppos = _pdb_records.find(mpos->second.first);
                 if (ppos != _pdb_records.end() && ndbformat.RefineField > 0) {
                      for (std::list<std::vector<std::string> >::iterator lppos =
                           ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                           if (refine_code == (*lppos)[ndbformat.RefineField - 1]) {
                                if ((*lppos)[mpos->second.second].empty()) {
                                     (*lppos)[mpos->second.second] = cs;
                                     update_flag = true;
                                }
                                break;
                           }
                      }
                 }
                 if (update_flag) continue;

                 _addNewRecord(mpos->second.first);
                 _updateRecordBack(mpos->second.first, mpos->second.second, cs);
                 if (ndbformat.RefineField > 0)
                      _updateRecordBack(mpos->second.first, ndbformat.RefineField - 1,
                                        refine_code);
            }
       }
}

void Maxit::_cif_to_ndb_process_refmac_refine_ls_restr_ncs(ISTable* t)
{
       int n_keywords = 8;
       const char *keywords[8] = { "", "dom_id", "pdbx_auth_asym_id", "pdbx_number",
             "rms_dev_position", "weight_position", "pdbx_type", "pdbx_ens_id" };

       const std::vector<std::string>& ColumnNames = t->GetColumnNames();

       std::map<int, std::string> field_item_mapping;
       field_item_mapping.clear();

       for (std::vector<std::string>::const_iterator
            pos = ColumnNames.begin(); pos != ColumnNames.end(); ++pos) {
            for (int i = 1; i < n_keywords; ++i) {
                 if (String::IsEqual(*pos, keywords[i], Char::eCASE_INSENSITIVE)) {
                      field_item_mapping.insert(std::make_pair(i, *pos));
                      break;
                 }
            }
       }
       if (field_item_mapping.empty()) return;

       _cif_to_ndb_process_single_token(t, "TMLPTL", field_item_mapping, true);
}

void Maxit::_cif_to_ndb_process_category(ISTable *t, const CIF_CATEGORY& Category,
                                         const bool& is_special)
{
       std::vector<std::string> SerialNo, idValue;

       std::string cs;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            for (std::vector<CIF_KEYWORD>::const_iterator kpos = Category.keywords.begin(); kpos != Category.keywords.end(); ++kpos) {
                 get_value(cs, t, i, kpos->name);
                 if (cs.empty()) continue;

                 for (std::vector<ndb_token_match>::const_iterator ipos = kpos->NdbInfo.begin(); ipos != kpos->NdbInfo.end(); ++ipos) {
                      std::string tokenid = ipos->NdbTokenName;
                      if (tokenid.empty()) continue;
                      std::string jrnlid  = ipos->JrnlTokenName;
                      int fieldno = ipos->NdbFieldNo - 1;
                      if (fieldno <= 0) continue;
                  
                      _cif_to_ndb_find_serial_no(t, Category, i, tokenid, idValue, SerialNo);

                      if (!Category.CifKeyField.empty()) {
                           bool found = true;
                           for (unsigned int k = 0; k < Category.CifKeyField.size(); k++) {
                                if (!idValue[k].empty() && !String::IsEqual(idValue[k],
                                     ipos->JrnlTokenName, Char::eCASE_INSENSITIVE)) {
                                     found = false;
                                     break;
                                }
                           }
                           if (!found) continue;
                      }

                      if (tokenid == "JRNL") {
                           if (String::IsEqual(jrnlid, "AUTH", Char::eCASE_INSENSITIVE) &&
                               String::IsEqual(Category.category, "citation", Char::eCASE_INSENSITIVE)) continue;

                           String::StripAndCompressWs(cs);
                           _cif_to_ndb_process_general(jrnlid, fieldno, ipos->NdbTokenKeyFieldNo, SerialNo, cs);
                      } else if (is_special)
                           _cif_to_ndb_process_special(tokenid, fieldno, ipos->NdbTokenKeyFieldNo, SerialNo, cs);
                      else _cif_to_ndb_process_general(tokenid, fieldno, ipos->NdbTokenKeyFieldNo, SerialNo, cs);
                 }
            }
       }
}

void Maxit::_cif_to_ndb_find_serial_no(ISTable* t, const CIF_CATEGORY& Category, const int& irow,
                                       const std::string& tokenid, std::vector<std::string> &idvalue,
                                       std::vector<std::string> &SerialNo)
{
       idvalue.clear();
       SerialNo.clear();

       std::string cs;

       if (!Category.CifKeyField.empty()) {
            for (std::vector<int>::const_iterator pos =
                 Category.CifKeyField.begin(); pos != Category.CifKeyField.end(); ++pos) {
                 get_value(cs, t, irow, Category.keywords[*pos].name);
                 if (cs.empty()) {
                      idvalue.push_back("");
                      if (tokenid == "JRNL")
                           SerialNo.push_back(String::IntToString(irow));
                      else SerialNo.push_back(String::IntToString(irow + 1));
                 } else {
                      std::string value = cs;
                      std::string ivalue = "";
                      if (tokenid == "JRNL")
                           idvalue.push_back("");
                      else {
                           for (std::vector<ndb_token_match>::const_iterator
                                ipos = Category.keywords[*pos].NdbInfo.begin();
                                ipos != Category.keywords[*pos].NdbInfo.end(); ++ipos) {
                                if (Category.category == "database_2") {
                                     if (cs == ipos->JrnlTokenName) {
                                          value.clear();
                                          ivalue = ipos->JrnlTokenName;
                                          break;
                                     }
                                } else {
                                     unsigned int str_len = ipos->JrnlTokenName.size();
                                     if (str_len != 0 && cs.size() >= str_len &&
                                         cs.substr(0, str_len) == ipos->JrnlTokenName) {
                                          value = cs.substr(str_len);
                                          ivalue = ipos->JrnlTokenName;
                                          break;
                                     }
                                }
                           }
                      }
                      idvalue.push_back(ivalue);
                      SerialNo.push_back(value);
                 }
            }
       } SerialNo.push_back(String::IntToString(irow + 1));
}

void Maxit::_cif_to_ndb_process_general(const std::string& tokenid, const int& fieldno,
                                        const std::vector<int>& SerialNoField,
                                        const std::vector<std::string>& SerialNo,
                                        const std::string& value)
{
       if (tokenid == "HEADER" && String::IsEqual(value, "TMP_ID", Char::eCASE_INSENSITIVE)) return;

       if (tokenid == "NDBFIL" && String::IsEqual(value, "DEPOSIT-0", Char::eCASE_INSENSITIVE)) {
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
                ppos = _pdb_records.find(tokenid);
            if (ppos == _pdb_records.end()) _addNewRecord(tokenid);

            _updateRecordFront(tokenid, fieldno, value, false);
            return;
       }

       if (tokenid == "AUTHOR") {
            std::string cs = value;
            ndb_name_conversion_first_last(cs);
            _updateRecordBack("AUTHOR", fieldno, cs, ", ", true);
            return;
       }

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat(tokenid);

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(tokenid);
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator
                 lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                 bool found_match_record = false;

                 if (!SerialNoField.empty()) {
                      bool found = false;
                      for (unsigned int i = 0; i < SerialNoField.size(); ++i) {
                           if (SerialNoField[i] > 0 && !String::IsEqual(SerialNo[i],
                               (*lppos)[SerialNoField[i]], Char::eCASE_INSENSITIVE)) {
                                found = true; break;
                           }
                      }
                      if (!found) found_match_record = true;
                 } else if (ndbformat.SeqField != 0 && String::IsEqual(SerialNo[0],
                            (*lppos)[ndbformat.SeqField - 1], Char::eCASE_INSENSITIVE))
                      found_match_record = true;
                 else if ((tokenid == "HEADER" || tokenid == "PDBFIL" || tokenid == "NDBFIL") &&
                          String::IsEqual(value, (*lppos)[fieldno], Char::eCASE_INSENSITIVE))
                      found_match_record = true;
                 else if ((*lppos)[fieldno].empty()) found_match_record = true;

                 if (!found_match_record) continue;

                 _cif_to_ndb_process_value(tokenid, fieldno, SerialNoField, value,
                                           ndbformat.SeqField - 1, *lppos);
                 return;
            }
       }

       _addNewRecord(tokenid);
       ppos = _pdb_records.find(tokenid);

       std::vector<std::string>& Field = ppos->second.back();
       if (!SerialNoField.empty()) {
            for (unsigned int i = 0; i < SerialNoField.size(); ++i) {
                 if (SerialNoField[i] > 0) Field[SerialNoField[i]] = SerialNo[i];
            }
       } else if (!SerialNo[0].empty() && ndbformat.SeqField != 0)
            Field[ndbformat.SeqField - 1] = SerialNo[0];

       _cif_to_ndb_process_value(tokenid, fieldno, SerialNoField, value,
                                 ndbformat.SeqField - 1, Field);
}

void Maxit::_cif_to_ndb_process_value(const std::string& tokenid, const int& fieldno, const
                                      std::vector<int>& SerialNoField, const std::string& value,
                                      const int& SeqField, std::vector<std::string>& Field)
{
       if (fieldno == SeqField) return;

       bool logical = false;
       for (std::vector<int>::const_iterator
            pos = SerialNoField.begin(); pos != SerialNoField.end(); ++pos) {
            if (fieldno == *pos) logical = true;
       }

       if (NdbToken::IsNdbToken(tokenid) && logical) return;

       if (tokenid == "AUTH" || tokenid == "EDIT") {
            std::string cs = value;
            ndb_name_conversion_first_last(cs);
            if (Field[fieldno].empty())
                 Field[fieldno] = cs;
            else Field[fieldno] += ", " + cs;
       } else if (tokenid == "REF") {
            int vField = 5;
            if (fieldno == vField + 1) Field[vField] = "V.";
            Field[fieldno] = value;
       } else Field[fieldno] = value;
}

void Maxit::_cif_to_ndb_process_special(const std::string& tokenid, const int& fieldno,
                                        const std::vector<int>& SerialNoField,
                                        const std::vector<std::string>& SerialNo,
                                        const std::string& value)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(tokenid);
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            bool found_match_record = false;
            if (!SerialNoField.empty()) {
                 bool found_mismatch = false;
                 for (unsigned int i = 0; i < SerialNoField.size(); ++i) {
                      if (SerialNoField[i] > 0 &&
                          SerialNo[i] != (*lppos)[SerialNoField[i]]) {
                           found_mismatch = true;
                           break;
                      }
                 }
                 if (!found_mismatch) found_match_record = true;
            }
            if (!found_match_record) continue;

            (*lppos)[fieldno] = value;
       }
}

void Maxit::_cif_to_ndb_process_matrix(ISTable *t, const CIF_CATEGORY& Category)
{
       std::string tokenid = "";
       std::map<std::string, int> tmp_map1, tmp_map2, tmp_map3;
       tmp_map1.clear();
       tmp_map2.clear();
       tmp_map3.clear();

       for (std::vector<CIF_KEYWORD>::const_iterator
            kpos = Category.keywords.begin(); kpos != Category.keywords.end(); ++kpos) {
            if (kpos->NdbInfo[0].NdbTokenName.empty() ||
                kpos->NdbInfo[0].NdbTokenName == "HEADER") continue;
            int ndb_fieldno = kpos->NdbInfo[0].NdbFieldNo - 1;
            if (ndb_fieldno < 0) continue;

            tokenid = kpos->NdbInfo[0].NdbTokenName;

            std::string::size_type idx = kpos->name.find_first_of('[');
            if (idx != std::string::npos) {
                 if (kpos->name.substr(idx + 1, 1) == "1")
                      tmp_map1.insert(std::make_pair(kpos->name, ndb_fieldno));
                 else if (kpos->name.substr(idx + 1, 1) == "2")
                      tmp_map2.insert(std::make_pair(kpos->name, ndb_fieldno));
                 else if (kpos->name.substr(idx + 1, 1) == "3")
                      tmp_map3.insert(std::make_pair(kpos->name, ndb_fieldno));
            } else {
                 tmp_map1.insert(std::make_pair(kpos->name, ndb_fieldno));
                 tmp_map2.insert(std::make_pair(kpos->name, ndb_fieldno));
                 tmp_map3.insert(std::make_pair(kpos->name, ndb_fieldno));
            }
       }
       if (tokenid.empty() || tmp_map1.empty() || tmp_map2.empty() ||
           tmp_map3.empty()) return;

       std::vector<std::map<std::string, int> > field_item_mapping;
       field_item_mapping.clear();
       field_item_mapping.push_back(tmp_map1);
       field_item_mapping.push_back(tmp_map2);
       field_item_mapping.push_back(tmp_map3);

       std::string cs, cs1;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; i++) {
            for (std::vector<std::map<std::string, int> >::const_iterator pos =
                 field_item_mapping.begin(); pos != field_item_mapping.end(); ++pos) {
                 _addNewRecord(tokenid);
                 for (std::map<std::string, int>::const_iterator
                      mpos = pos->begin(); mpos != pos->end(); ++mpos) {
                      get_value_clean(cs, t, i, mpos->first);
                      if (cs.empty()) continue;

                      String::LowerCase(cs, cs1);
                      if (cs1 == "given")
                           cs = "1";
                      else if (cs1 == "generate")
                           cs.clear();
                      _updateRecordBack(tokenid, mpos->second, cs);
                 }
            }
       }
}
/*
static void ndb_name_conversion_first_last(std::string &str)
{
       if (str.empty()) return;

       if (SgCenter::IsSgCenter(str)) return;

       std::string::size_type pos = str.find(",");
       if (pos != std::string::npos) {
            std::string last = str.substr(0, pos);
            std::string first = str.substr(pos+1);
            String::StripAndCompressWs(last);
            String::StripAndCompressWs(first);
            if (first.find(".") != std::string::npos)
                 str = first + last;
            else str = first + " " + last;
       }
}
*/
