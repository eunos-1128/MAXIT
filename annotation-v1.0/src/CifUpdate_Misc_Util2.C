/*
FILE:     CifUpdate_Misc_Util2.C
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
#include "NdbRefn.h"
#include "utillib.h"

#define NUM_CAT_ITEM_TO_LOWER     5

static const char *_cat_items_to_lower[NUM_CAT_ITEM_TO_LOWER][2] = {
       { "pdbx_database_status",    "author_approval_type" },
       { "pdbx_nmr_representative", "selection_criteria"   },
       { "pdbx_nmr_software",       "classification"       },
       { "software",                "classification"       },
       { "symmetry",                "cell_setting"         }
};

#define NUM_SC   14

static const char *skip_categories[NUM_SC] = {
       "atom_site",
       "atom_site_anisotrop",
       "database_PDB_remark",
       "entity_poly",
       "entity_src_gen",
       "entity_src_nat",
       "ndb_original_ndb_coordinates",
       "ndb_original_pdb_coordinates",
       "pdbx_database_remark",
       "pdbx_entity_src_syn",
       "pdbx_original_pdb_coordinates",
       "pdbx_pdb_compnd",
       "pdbx_pdb_source",
       "struct_ref"
};

void AnnotationObj::_cif_update_coordinate_model_type(Block &block)
{
       if (_molecules.empty()) return;

       ISTable *t = _newTablePtr("pdbx_coordinate_model");

       // bool has_ATOMN = false;
       int irow = 0;
       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") {
                 // if (chain->chain_type() == "ATOMN") has_ATOMN = true;
                 if (chain->ca_or_p_atom_only()) {
                      t->AddRow();
                      t->UpdateCell(irow, "asym_id", chain->ChainID());
                      if (chain->chain_type() == "ATOMP")
                           t->UpdateCell(irow, "type", "CA ATOMS ONLY");
                      else t->UpdateCell(irow, "type", "P ATOMS ONLY");
                      irow++;
                 }
            }
            chain = _molecules[0]->GetNextChain(); 
       }
       if (t->GetNumRows() > 0)
            block.WriteTable(t);
       else {
            delete t;
            deleteTable(block, "pdbx_coordinate_model");
       }

       // if (!has_ATOMN) return;

       t = getTablePtr(block, "database_2");
       if (!t) return;

       std::string pdbid = getUpperCasePDBID(block);
       // std::string depid = getUpperCaseDEPID(block);
       std::string rcsbid = getUpperCaseRCSBID(block);
       std::string cs1, cs2;

       std::vector<unsigned int> deleterows;
       deleterows.clear();

       // bool has_NDBID = false;
       unsigned int rowNo = t->GetNumRows();
       for (unsigned int i = 0; i < rowNo; ++i) {
            get_value_clean_upper(cs1, t, i, "database_id");
            get_value_clean_upper(cs2, t, i, "database_code");
            if (cs1 == pdbid) {
                 deleterows.push_back(i);
                 continue;
            } else if ((cs1 == "NDB") && (cs2 == pdbid || cs2 == rcsbid)) {
                 deleterows.push_back(i);
                 continue;
            } else if ((cs1 == "RCSB") && (cs2 == pdbid)) {
                 deleterows.push_back(i);
                 continue;
            }
/*
            if (cs1 == "NDB") {
                 has_NDBID = true;
                 if (cs2.empty() || (cs2 == depid && (cs2.substr(0, 4) == "RCSB" || cs2.substr(0, 2) == "D_")))
                      t->UpdateCell(i, "database_code", pdbid);
            }
*/
       }
       if (!deleterows.empty()) t->DeleteRows(deleterows);

/*
       if (!has_NDBID) {
            rowNo = t->GetNumRows();
            t->AddRow();
            t->UpdateCell(rowNo, "database_id", "NDB");
            t->UpdateCell(rowNo, "database_code", pdbid);
       }
*/

       block.WriteTable(t);

}

void AnnotationObj::_cif_update_deposit_process_sites(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_database_status");
       if (!t) return; 

       std::vector<std::string> sites;
       sites.clear();
       sites.push_back("deposit_site");
       sites.push_back("process_site");

       std::string cs;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            for (std::vector<std::string>::const_iterator pos = sites.begin(); pos != sites.end(); ++pos) {
                 get_value_clean_upper(cs, t, i, *pos);
                 if (cs == "OSAKA")
                      t->UpdateCell(i, *pos, "PDBJ");
                 else if (cs == "EBI")
                      t->UpdateCell(i, *pos, "PDBE");
            }
       }
       block.WriteTable(t);
}

void AnnotationObj::_cif_update_prerelease_seq(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_database_status");
       if (!t) return; 

       std::string status;
       get_value_clean_upper(status, t, 0, "status_code");
       if (status == "REL" || status == "OBS") {
            deleteTable(block, "pdbx_prerelease_seq");
            return;
       }

       get_value_clean_upper(status, t, 0, "dep_release_code_sequence");
       if (status != "RELEASE NOW") return;

       t = getTablePtr(block, "entity_poly");
       if (!t) return;

       std::string entity_id, seq;
       int row = 0;
       ISTable *t1 = _newTablePtr("pdbx_prerelease_seq");

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(entity_id, t, i, "entity_id");
            get_value(seq, t, i, "pdbx_seq_one_letter_code");
            if (entity_id.empty() || seq.empty()) continue;

            t1->AddRow();
            t1->UpdateCell(row, "entity_id", entity_id);
            t1->UpdateCell(row, "seq_one_letter_code", seq);
            row++;
       }
       if (row > 0)
            block.WriteTable(t1);
       else {
            delete t1;
            deleteTable(block, "pdbx_prerelease_seq");
       }
}

void AnnotationObj::_cif_update_miscellaneous_categories(Block& block)
{
       // delete pdbx_audit category
       // deleteTable(block, "pdbx_audit");

       // delete empty pdbx_SG_project category
       ISTable *t = getTablePtr(block, "pdbx_SG_project");
       if (is_empty_table(t, "id")) deleteTable(block, "pdbx_SG_project");

       // add entry category
       t = _getTablePtr(block, "entry");
       if (!t) {
            t = new ISTable("entry");
            t->AddColumn("id");
            t->AddRow();
       }       
       t->UpdateCell(0, "id", _StructureId);
       block.WriteTable(t);

       std::string cs;

       // update cell category
/*
       t = _getTablePtr(block, "cell");
       if (t) {
            get_value_clean(cs, t, 0, "Z_PDB");
            if (cs.empty() || cs == "0") {
                 t->UpdateCell(0, "Z_PDB", "1");
                 block.WriteTable(t);
            }
       }
*/

       // update database_PDB_rev category
       t = _getTablePtr(block, "database_PDB_rev");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_lower(cs, t, i, "status");
                 if (cs == "y") t->UpdateCell(i, "status", "full release");
            }
            block.WriteTable(t);
       }

       // update to lowercase
       for (int k = 0; k < NUM_CAT_ITEM_TO_LOWER; ++k) {
            t = _getTablePtr(block, _cat_items_to_lower[k][0]);
            if (!t) continue;
            if (!t->IsColumnPresent(_cat_items_to_lower[k][1])) continue;

            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_lower(cs, t, i, _cat_items_to_lower[k][1]);
                 if (cs.empty()) continue;
                 t->UpdateCell(i, _cat_items_to_lower[k][1], cs);
            }
            block.WriteTable(t);
       }

       // update database_PDB_caveat category
       t = _getTablePtr(block, "database_PDB_caveat");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 t->UpdateCell(i, "id", String::IntToString(i + 1));
            }
            block.WriteTable(t);
       }

       // update pdbx_database_related category
       t = _getTablePtr(block, "pdbx_database_related");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_lower(cs, t, i, "content_type");
                 if (cs.empty() || cs == "other") t->UpdateCell(i, "content_type", "unspecified");
                 get_value_clean(cs, t, i, "details");
                 if (cs.empty()) t->UpdateCell(i, "details", ".");
            }
            block.WriteTable(t);
       }

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > related_ids;
       related_ids.clear();

       cs = getUpperCaseBMRBID(block);
       if (!cs.empty()) {
            data.clear();
            data.push_back("BMRB");
            data.push_back(cs);
            data.push_back("unspecified");
            related_ids.push_back(data);
       }

       cs = getUpperCaseEMDBID(block);
       if (!cs.empty()) {
            data.clear();
            data.push_back("EMDB");
            data.push_back(cs);
            data.push_back("associated EM volume");
            related_ids.push_back(data);
       }

       if (!related_ids.empty()) {
            bool found = false;
            std::string db_name, db_id;
            if (!t) t = _newTablePtr("pdbx_database_related");
            int rowNo = t->GetNumRows();
            for (unsigned i = 0; i < related_ids.size(); ++i) {
                 for (int j = 0; j < rowNo; ++j) {
                      get_value_clean_upper(db_name, t, j, "db_name");
                      get_value_clean_upper(db_id, t, j, "db_id");
                      get_value_clean(cs, t, j, "content_type");
                      if (db_name == related_ids[i][0] && db_id == related_ids[i][1] && cs == related_ids[i][2]) {
                           found = true;
                           break;
                      }
                 }
                 if (!found) {
                      t->AddRow();
                      t->UpdateCell(rowNo, "db_name",      related_ids[i][0]);
                      t->UpdateCell(rowNo, "db_id",        related_ids[i][1]);
                      t->UpdateCell(rowNo, "content_type", related_ids[i][2]);
                      t->UpdateCell(rowNo, "details",      ".");
                      rowNo++;
                 }
            }
            block.WriteTable(t);
       }

       // update reflns category
       t = _getTablePtr(block, "reflns");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 _remove_infinite(t, "d_resolution_low", i);
                 _remove_comma(t, "number_obs", i);
                 _remove_comma(t, "number_all", i);
            }
            block.WriteTable(t);
       }

       // update cifation category
       t = _getTablePtr(block, "citation");
       if (t) {
            NdbRefn::initialize();
            NdbRefn::Read(*_logIo, _rcsbroot);

            std::map<std::string, std::string> jrnlInfo;

            std::string cs1, cs2;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 // get_value_clean(cs, t, i, "id");
                 // if (cs.empty()) continue;
                 get_value_clean(cs, t, i, "journal_abbrev");
                 if (cs.empty()) {
                      t->UpdateCell(i, "journal_abbrev", "To Be Published");
                      cs = "To Be Published";
                 }
                 bool is_epub = true;
                 get_value_clean(cs1, t, i, "page_first");
                 if (!cs1.empty()) is_epub = false;
                 get_value_clean(cs1, t, i, "journal_id_ISSN");
                 get_value_clean(cs2, t, i, "year");
                 jrnlInfo.clear();
                 if (!cs1.empty()) jrnlInfo = NdbRefn::GetJournalInfo(cs, cs1, cs2);
                 else jrnlInfo = NdbRefn::GetJournalInfo(cs, is_epub, cs2);
                 if (jrnlInfo.empty()) continue;

                 if (cs1.empty()) t->UpdateCell(i, "journal_id_ISSN", jrnlInfo["issn_code"]);
                 get_value_clean(cs2, t, i, "journal_id_ASTM");
                 if (cs2.empty()) t->UpdateCell(i, "journal_id_ASTM", jrnlInfo["astm_list"]);
                 get_value_clean(cs2, t, i, "country");
                 if (cs2.empty()) t->UpdateCell(i, "country", jrnlInfo["country"]);
                 get_value_clean(cs2, t, i, "journal_id_CSD");
                 if (cs2.empty()) t->UpdateCell(i, "journal_id_CSD", jrnlInfo["coden"]);
            }
            block.WriteTable(t);
       }
}

void AnnotationObj::_cif_update_remove_empty_row_and_tables(Block& block)
{
       std::set<std::string> skip_category_set, alloweditems;
       skip_category_set.clear();
       for (int i = 0; i < NUM_SC; ++i) skip_category_set.insert(skip_categories[i]);

       alloweditems.clear();
       alloweditems.insert("entry_id");
       alloweditems.insert("entity_id");

       std::vector<unsigned int> deleterows;

       std::string cs, item;

       std::vector<std::string> TableNames;

       block.GetTableNames(TableNames);
       for (std::vector<std::string>::const_iterator tpos = TableNames.begin(); tpos != TableNames.end(); ++tpos) {
            if (skip_category_set.find(*tpos) != skip_category_set.end())
                 continue;

            ISTable *t = _getTablePtr(block, *tpos);
            if (!t) {
                 deleteTable(block, *tpos); // remove empty table
                 continue;
            }

            unsigned int rowNo = t->GetNumRows();
            const std::vector<std::string>& ColumnNames = t->GetColumnNames();
            deleterows.clear();
            for (unsigned int j = 0; j < rowNo; ++j) {
                 if (is_empty_row(t, ColumnNames, j, alloweditems)) {
                      deleterows.push_back(j);
                      continue;
                 }
                 for (std::vector<std::string>::const_iterator cpos = ColumnNames.begin(); cpos != ColumnNames.end(); ++cpos) {
                      cs = (*t)(j, *cpos);
                      if (*cpos == "entry_id")
                           cs = _StructureId;
                      else {
                           String::StripLeadingWs(cs);
                           String::StripTrailingWs(cs);
                      }
                      t->UpdateCell(j, *cpos, cs);
                 }
            }

            // remove empty row(s)
            if (!deleterows.empty()) t->DeleteRows(deleterows);

            rowNo = t->GetNumRows();
            if (rowNo == 0)
                 deleteTable(block, *tpos); // remove empty table
            else {
                 item.clear();
                 if (t->IsColumnPresent("pdbx_ordinal"))
                      item = "pdbx_ordinal";
                 else if (t->IsColumnPresent("ordinal"))
                      item = "ordinal";
                 else if ((*tpos == "pdbx_database_proc") && t->IsColumnPresent("cycle_id"))
                      item = "cycle_id";
                 // update ordinal/pdbx_ordinal
                 if (!item.empty() && *tpos != "pdbx_nmr_software") {
                      bool update_flag = true;
                      if (*tpos == "pdbx_nmr_software") {
                           for (unsigned int j = 0; j < rowNo; ++j) {
                                get_value_clean(cs, t, j, item);
                                if (!cs.empty()) {
                                     update_flag = false;
                                     break;
                                }
                           }
                      }
                      if (update_flag) {
                           for (unsigned int j = 0; j < rowNo; ++j) {
                                t->UpdateCell(j, item, String::IntToString(j + 1));
                           }
                      }
                 }
                 if (*tpos == "audit_author" || *tpos == "citation_author") {
                      for (unsigned int j = 0; j < rowNo; ++j) {
                           get_value_clean(cs, t, j, "name");
                           if (cs.empty() || cs.find(", ") != std::string::npos) continue;
                           replace_string(cs, ",", ", ");
                           t->UpdateCell(j, "name", cs);
                      }
                 }
                 block.WriteTable(t);
            }
       }

       deleteTable(block, "database");
       deleteTable(block, "atom_sites_alt");
}
