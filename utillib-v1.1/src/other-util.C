/*
FILE:     other-util.C
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

#include "utillib.h"

#define NUM_STOP_WORD  65

static const char *_stop_words[NUM_STOP_WORD] = {
       "3D CRYSTAL STRUCTURE",
       "3D-STRUCTURE",
       "ACT_SITE",
       "BINDING",
       "CA_BIND",
       "CARBOHYD",
       "CHAIN",
       "COMPLETE PROTEOME",
       "CONFLICT",
       "CRYO-ELECTRON MICROSCOPY",
       "CRYSTALLOGRAPHY",
       "CRYSTAL STRUCTURE",
       "DIRECT PROTEIN SEQUENCING",
       "DISULFID",
       "DNA_BIND",
       "DOMAIN",
       "ELECTRON CRYSTALLOGRAPHY",
       "ELECTRON DIFFRACTION",
       "ELECTRON MICROSCOPY",
       "ELECTRON TOMOGRAPHY",
       "FIBER DIFFRACTION",
       "FIBRE DIFFRACTION",
       "HELIX TURN",
       "INIT_MET",
       "MECHANISM",
       "METAL",
       "MOD_RES",
       "MUTAGEN",
       "NATIVE CRYSTAL",
       "NEUTRON DIFFRACTION",
       "NMR",
       "NMR STRUCTURE",
       "NON_CONS",
       "NON_TER",
       "NP_BIND",
       "NUCLEAR MAGNETIC RESONANCE",
       "PEPTIDE",
       "PROPEP",
       "REPEAT",
       "REPEAT",
       "SE_CYS",
       "SEQUENCING",
       "SIGNAL",
       "SIGNAL",
       "SITE",
       "SOLID-STATE NMR",
       "SOLUTION NMR",
       "SOLUTION SCATTERING",
       "STRAND",
       "STRUCTURE",
       "THEORETICAL MODEL",
       "THIOETH",
       "THIOLEST",
       "TRANSIT",
       "TRANSMEM",
       "ULTRA-HIGH RESOLUTION",
       "UNSURE",
       "VARIANT",
       "VARSPLIC",
       "X-RAY",
       "X-RAY CRYSTALLOGRAPHY",
       "X-RAY DIFFRACTION",
       "X-RAY POWDER DIFFRACTION",
       "X-RAY STRUCTURE",
       "ZN_FING"
};

static void get_cif_mapping(const std::string& filename, std::map<std::string, std::string>& category_mapping,
                            std::map<std::string, std::map<std::string, std::string> >& item_mapping,
                            std::list<std::pair<std::string, std::string > >& rename_categories,
                            std::map<std::string, std::list<std::pair<std::string, std::string> > >& rename_items,
                            std::map<std::string, std::set<std::string> >& remove_items, const int& status);
static void process_one_cif_mapping(const std::string& item_v5, const std::string& item_v4, std::map<std::string, std::string>& category_mapping,
                                    std::map<std::string, std::map<std::string, std::string> >& item_mapping, const int& status);
static void _multiple_sort(const std::vector<std::vector<int> >& matrix, const std::vector<unsigned int>& start_index,
                           const unsigned int field_index, std::vector<int>& score);
static void separate_string_by_parenthesis(std::string &cifstring, std::vector<std::string> &list);
static void separate_string_by_comma(std::string &cifstring, std::vector<std::string> &list);

bool reset_entry_id(Block& block, const std::string& entry_id, const bool& include_entry)
{
       bool change_flag = false;
       std::vector<std::string> tableNames;
       block.GetTableNames(tableNames);
       for (std::vector<std::string>::const_iterator pos = tableNames.begin(); pos != tableNames.end(); ++pos) {
            ISTable *t = getTablePtr(block, *pos);
            if (!t) {
                 deleteTable(block,*pos);
                 continue;
            }
            if (!t->IsColumnPresent("entry_id")) continue;

            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 if (update_value(t, "entry_id", i, entry_id, true)) change_flag = true;
            }
            block.WriteTable(t);
       }

       ISTable *t = getTablePtr(block, "entry");
       if (!t && include_entry) {
            t = new ISTable("entry");
            t->AddColumn("id");
            t->AddRow();
       }
       if (t) {
            if (update_value(t, "id", 0, entry_id, true)) change_flag = true;
            block.WriteTable(t);
       }

       return change_flag;
}

void master_internal_conversion(const std::string& mapping_file, Block& block, const int& status)
{
       std::map<std::string, std::string> category_mapping;
       std::map<std::string, std::map<std::string, std::string> > item_mapping;
       std::list<std::pair<std::string, std::string > > rename_categories;
       std::map<std::string, std::list<std::pair<std::string, std::string> > > rename_items;
       std::map<std::string, std::set<std::string> > remove_items;

       get_cif_mapping(mapping_file, category_mapping, item_mapping, rename_categories, rename_items, remove_items, status);
       if (category_mapping.empty()) return; 

       std::vector<std::string> tableNames;
       block.GetTableNames(tableNames);
       for (std::vector<std::string>::const_iterator vpos = tableNames.begin(); vpos != tableNames.end(); ++vpos) {
            ISTable *t = getTablePtr(block, *vpos);
            if (t == NULL) {
                 deleteTable(block, *vpos);
                 continue;
            }

            std::map<std::string, std::string>::const_iterator pos1 = category_mapping.find(*vpos);
            std::map<std::string, std::map<std::string, std::string> >::const_iterator pos2 = item_mapping.find(*vpos);
            if (pos1 == category_mapping.end() && pos2 == item_mapping.end()) continue;

            if (pos2 != item_mapping.end()) {
                 for (std::map<std::string, std::string>::const_iterator mpos = pos2->second.begin(); mpos != pos2->second.end(); ++mpos) {
                      rename_item(t, mpos->first, mpos->second);
                 }
            }
            block.WriteTable(t);

            if (pos1 != category_mapping.end()) {
                 t = getTablePtr(block, pos1->second);
                 if (t == NULL) {
                      deleteTable(block, pos1->second);
                      block.RenameTable(pos1->first, pos1->second);
                 }
            }
       }

       if (status == INTERNAL_TO_MASTER) {
            for (std::list<std::pair<std::string, std::string > >::const_iterator pos = rename_categories.begin(); pos != rename_categories.end(); ++pos) {
                 if (block.IsTablePresent(pos->first) && !block.IsTablePresent(pos->second)) {
                      block.RenameTable(pos->first, pos->second);
                 }
            }

            for (std::map<std::string, std::list<std::pair<std::string, std::string> > >::const_iterator
                 mpos = rename_items.begin(); mpos != rename_items.end(); ++mpos) {
                 ISTable *t = getTablePtr(block, mpos->first);
                 if (!t) continue;

                 for (std::list<std::pair<std::string, std::string> >::const_iterator lpos = mpos->second.begin(); lpos != mpos->second.end(); ++lpos) {
                      rename_item(t, lpos->first, lpos->second);
                 }
                 block.WriteTable(t);
            }

            for (std::map<std::string, std::set<std::string> >::const_iterator mpos = remove_items.begin(); mpos != remove_items.end(); ++mpos) {
                 ISTable *t = getTablePtr(block, mpos->first);
                 if (!t) continue;

                 for (std::set<std::string>::const_iterator
                      spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                      if (t->IsColumnPresent(*spos)) t->DeleteColumn(*spos);
                 }
                 block.WriteTable(t);
            }

            ISTable *t1 = getTablePtr(block, "entity_keywords");
            ISTable *t2 = getTablePtr(block, "entity");
            if (!t1 || !t2) return;

            std::vector<std::string> items;
            items.clear();
            items.push_back("pdbx_ec");
            items.push_back("pdbx_mutation");
            items.push_back("pdbx_fragment");
            check_missing_item(t2, items);

            std::string cs1, cs2;
            std::map<std::string, std::string> mapping;

            int rowNo1 = t1->GetNumRows();
            int rowNo2 = t2->GetNumRows();
            for (std::vector<std::string>::const_iterator
                 pos = items.begin(); pos != items.end(); ++pos) {
                 mapping.clear();
                 for (int i = 0; i < rowNo1; ++i) {
                      get_value_clean(cs1, t1, i, "entity_id");
                      if (cs1.empty()) continue;
                      get_value_clean(cs2, t1, i, *pos);
                      if (cs2.empty()) continue;
                      mapping.insert(std::make_pair(cs1, cs2));
                 }
                 if (t1->IsColumnPresent(*pos)) t1->DeleteColumn(*pos);
                 if (mapping.empty()) continue;

                 for (int i = 0; i < rowNo2; ++i) {
                      get_value_clean(cs1, t2, i, "id");
                      std::map<std::string, std::string>::const_iterator
                          mpos = mapping.find(cs1);
                      if (mpos != mapping.end()) t2->UpdateCell(i, *pos, mpos->second);
                 }
            }
            block.WriteTable(t2);
       } else {
            ISTable *t = getTablePtr(block, "entity");
            if (!t) return;

            std::map<std::string, std::string> mapping, item_mapping;

            item_mapping.clear();
            item_mapping.insert(std::make_pair("pdbx_ec", "ndb_ec"));
            item_mapping.insert(std::make_pair("pdbx_mutation", "ndb_mutation"));
            item_mapping.insert(std::make_pair("pdbx_fragment", "ndb_fragment"));

            std::map<int, std::map<std::string, std::string> > value_mapping;
            value_mapping.clear();

            std::string cs;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 mapping.clear();
                 for (std::map<std::string, std::string>::const_iterator
                      mpos = item_mapping.begin(); mpos != item_mapping.end(); ++mpos) {
                      get_value_clean(cs, t, i, mpos->first);
                      if (!cs.empty()) mapping.insert(std::make_pair(mpos->second, cs));
                 }
                 if (mapping.empty()) continue;

                 get_value_clean(cs, t, i, "id");
                 value_mapping.insert(std::make_pair(atoi(cs.c_str()), mapping));
            }

            if (t->IsColumnPresent("pdbx_ec")) t->DeleteColumn("pdbx_ec");
            if (t->IsColumnPresent("pdbx_mutation")) t->DeleteColumn("pdbx_mutation");
            if (t->IsColumnPresent("pdbx_fragment")) t->DeleteColumn("pdbx_fragment");

            block.WriteTable(t);

            if (value_mapping.empty()) return;

            mapping.clear();
            mapping.insert(std::make_pair("ndb_ec", "?"));
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, "id");
                 int entity_id = atoi(cs.c_str());
                 if (value_mapping.find(entity_id) != value_mapping.end()) continue;
                 value_mapping.insert(std::make_pair(entity_id, mapping));
            }

            std::vector<std::string> items;
            items.clear();

            t = getTablePtr(block, "entity_keywords");
            if (t) {
                 items = t->GetColumnNames();
                 rowNo = t->GetNumRows();
                 for (int i = 0; i < rowNo; ++i) {
                      mapping.clear();
                      for (std::vector<std::string>::const_iterator
                           pos = items.begin(); pos != items.end(); ++pos) {
                           if (*pos == "entity_id") continue;
                           get_value_clean(cs, t, i, *pos);
                           if (!cs.empty()) mapping.insert(std::make_pair(*pos, cs));
                      }
                      if (mapping.empty()) continue;

                      get_value_clean(cs, t, i, "entity_id");
                      int entity_id = atoi(cs.c_str());
                      std::map<int, std::map<std::string, std::string> >::iterator
                          ipos = value_mapping.find(entity_id);
                      if (ipos == value_mapping.end())
                           value_mapping.insert(std::make_pair(entity_id, mapping));
                      else {
                           for (std::map<std::string, std::string>::const_iterator
                                mpos = mapping.begin(); mpos != mapping.end(); ++mpos) {
                                ipos->second.insert(std::make_pair(mpos->first, mpos->second));
                           }
                      }
                 }
            } else items.push_back("entity_id");

            items.push_back("ndb_ec");
            items.push_back("ndb_mutation");
            items.push_back("ndb_fragment");

            t = add_new_table("entity_keywords", items);

            rowNo = 0;
            for (std::map<int, std::map<std::string, std::string> >::const_iterator
                 ipos = value_mapping.begin(); ipos != value_mapping.end(); ++ipos) {
                 t->AddRow();
                 t->UpdateCell(rowNo, "entity_id", String::IntToString(ipos->first));
                 for (std::map<std::string, std::string>::const_iterator
                      mpos = ipos->second.begin(); mpos != ipos->second.end(); ++mpos) {
                      t->UpdateCell(rowNo, mpos->first, mpos->second);
                 }
                 rowNo++;
            }
            block.WriteTable(t);
       }
}

void get_order_list(const std::string& filename, Block& block, std::vector<std::string>& order_list)
{
       order_list.clear();

       struct stat stattext;
       if (stat(filename.c_str(), &stattext) != 0) return;

       int count = 0;
       std::map<std::string, int> order_index_mapping;
       order_index_mapping.clear();

       std::string cs;

       FILE *fp = fopen(filename.c_str(), "r");
       while (!feof(fp)) {
            get_one_line(fp, cs);
            if (feof(fp)) break;

            String::StripAndCompressWs(cs);
            if (cs.empty() || cs == "atom_site" || cs == "atom_site_anisotrop") continue;

            order_index_mapping.insert(std::make_pair(cs, count));
            count++;
       }
       fclose (fp);

       std::map<int, std::string> order_list_mapping;
       order_list_mapping.clear();
       bool has_atom_site = false;
       bool has_atom_site_anisotrop = false;
       std::set<std::string> extra_list;
       extra_list.clear();

       std::vector<std::string> tableNames;
       block.GetTableNames(tableNames);
       for (std::vector<std::string>::const_iterator
            vpos = tableNames.begin(); vpos != tableNames.end(); ++vpos) {
            if (*vpos == "atom_site")
                 has_atom_site = true;
            else if (*vpos == "atom_site_anisotrop")
                 has_atom_site_anisotrop = true;
            else {
                 std::map<std::string, int>::const_iterator
                     mpos = order_index_mapping.find(*vpos);
                 if (mpos != order_index_mapping.end())
                      order_list_mapping.insert(std::make_pair(mpos->second, mpos->first));
                 else extra_list.insert(*vpos);
            }
       }

       if (!order_list_mapping.empty()) {
            for (std::map<int, std::string>::const_iterator
                 mpos = order_list_mapping.begin(); mpos != order_list_mapping.end(); ++mpos) {
                 order_list.push_back(mpos->second);
            }
       }
       if (!extra_list.empty()) {
            for (std::set<std::string>::const_iterator
                 spos = extra_list.begin(); spos != extra_list.end(); ++spos) {
                 order_list.push_back(*spos);
            }
       }
       if (has_atom_site) order_list.push_back("atom_site");
       if (has_atom_site_anisotrop) order_list.push_back("atom_site_anisotrop");
}

void get_validation_letter_template(const std::string& filename, const bool is_xray_entry,
                std::vector<std::vector<std::string> > &head, std::map<std::string, std::string>& mapping,
                std::vector<std::vector<std::string> >& validation_tokens)
{
       head.clear();
       mapping.clear();
       validation_tokens.clear();

       struct stat stattext;
       if (stat(filename.c_str(), &stattext) != 0) return;

       std::string /* odbfile, */ cs, cs1;
       // get_temp_filename(odbfile);
       CifFile *fobj = read_cif_file("", filename, cs);
       Block &block = fobj->GetBlock(fobj->GetFirstBlockName());

       std::map<std::string, int> head_index;
       head_index.clear();
       std::vector<std::string> data;
       ISTable *t = getTablePtr(block, "letter_order");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs, t, i, "token");
                 head_index.insert(std::make_pair(cs, head.size()));
                 data.clear();
                 data.push_back(cs);
                 head.push_back(data);
            }
       }

       t = getTablePtr(block, "letter_sub_order");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs, t, i, "token");
                 get_value(cs1, t, i, "sub_token");
                 if (is_xray_entry && cs1 == "close_contact_other") continue;
                 else if (!is_xray_entry && (cs1 == "close_contact_asym" ||
                          cs1 == "close_contact_symm")) continue;
                 std::map<std::string, int>::iterator pos = head_index.find(cs);
                 if (pos != head_index.end()) head[pos->second].push_back(cs1); 
            }
       }

       t = getTablePtr(block, "letter_template");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs, t, i, "token");
                 get_value(cs1, t, i, "text");
                 mapping.insert(std::make_pair(cs, cs1));
            }
       }

       t = getTablePtr(block, "validation");
       if (t) {
            const std::vector<std::string>& ColumnNames = t->GetColumnNames();
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 data.clear();
                 for (unsigned int j = 0; j < ColumnNames.size(); ++j) {
                      get_value(cs, t, i, ColumnNames[j]);
                      data.push_back(cs);
                 }
                 validation_tokens.push_back(data);
            }
       }

       if (fobj) delete fobj;

       // remove(odbfile.c_str());
}

void get_no_loop_category_set(const std::string& filePath, std::set<std::string>& no_loop_categories)
{
       no_loop_categories.clear();

       struct stat statbuf;
       if (stat(filePath.c_str(), &statbuf) != 0) return;

       std::string cs;
       FILE *fp = fopen(filePath.c_str(), "r");
       while (!feof(fp)) {
            get_one_line(fp, cs);
            if (feof(fp)) break;
            String::StripAndCompressWs(cs);
            no_loop_categories.insert(cs);
       }
       fclose (fp);
}

std::set<std::string> get_stop_words(const std::string& filename)
{
       std::set<std::string> Stop_words;
       Stop_words.clear();

       std::string cs;
       struct stat statbuf;
       if (stat(filename.c_str(), &statbuf) == 0) {
            FILE *fp = fopen(filename.c_str(), "r");
            while (!feof(fp)) {
                 get_one_line(fp, cs);
                 if (feof(fp)) break;
                 String::StripAndCompressWs(cs);
                 String::UpperCase(cs);
                 Stop_words.insert(cs);
            }
            fclose (fp);
       } else {
            for (int i = 0; i < NUM_STOP_WORD; ++i) {
                 Stop_words.insert(_stop_words[i]);
            }
       }
       return Stop_words;
}

void remove_redundant_keywords(const std::string& filename, std::string& keyword, std::string& header)
{
       if (!header.empty()) {
            bool is_complex = false;
            for (unsigned int i = 0; i < header.size(); ++i) {
                 if (header[i] == '/') {
                      header[i] = '-';
                      is_complex = true;
                 }
            }
            if (is_complex) header += " complex";
       }

       std::set<std::string> Stop_words = get_stop_words(filename);

       std::string cs;
       std::vector<std::string> data, header_and_keywords;
       header_and_keywords.clear();
       if (!keyword.empty()) header_and_keywords.push_back(keyword);
       if (!header.empty()) header_and_keywords.push_back(header);
       if (header_and_keywords.empty()) return;

       keyword.clear();
       for (unsigned int l = 0; l < header_and_keywords.size(); ++l) {
            get_wordarray(data, header_and_keywords[l], ",");
            for (unsigned int i = 0; i < data.size(); ++i) {
                 String::UpperCase(data[i], cs);
                 if (cs.find("/") != std::string::npos && cs.find("COMPLEX") != std::string::npos) {
                      for (unsigned int j = 0; j < cs.size(); ++j) {
                           if (cs[j] == '/') cs[j] = '-';
                      }
                      for (unsigned int j = 0; j < data[i].size(); ++j) {
                           if (data[i][j] == '/') data[i][j] = '-';
                      }
                 }
                 if (Stop_words.find(cs) != Stop_words.end()) continue;
                 Stop_words.insert(cs);
                 if (keyword != "") keyword += ", ";
                 keyword += data[i];
            }
       }
}

void _MultipleSort(const std::vector<std::vector<int> >& score_matrix, std::vector<int>& index)
{
       index.clear();
       if (score_matrix.size() == 1) {
            index.push_back(0);
            return;
       }

       std::vector<int> score;
       score.clear();
       std::vector<unsigned int> idx;
       idx.clear();
       for (unsigned int i = 0; i < score_matrix.size(); ++i) {
            idx.push_back(i);
            score.push_back(0);
       }

       _multiple_sort(score_matrix, idx, 0, score);

       std::multimap<int, int> mapping; 
       mapping.clear();
       for (unsigned int i = 0; i < score.size(); ++i) {
            mapping.insert(std::make_pair(score[i], i));
       }

       for (std::multimap<int, int>::iterator pos = mapping.begin();
                                         pos != mapping.end(); ++pos) {
            index.push_back(pos->second);
       }
}

void add_type(std::map<std::string, int>& all_types, const std::string& type)
{
       std::map<std::string, int>::iterator
           pos = all_types.find(type);
       if (pos != all_types.end()) {
            pos->second++;
            return;
       }

       all_types.insert(std::make_pair(type, 1));
}

std::string get_type(const std::map<std::string, int>& all_types, const std::string& d_value)
{
       std::multimap<int, std::string> type_index;
       type_index.clear();
       for (std::map<std::string, int>::const_iterator
            mpos = all_types.begin(); mpos != all_types.end(); ++mpos) {
            type_index.insert(std::make_pair(mpos->second, mpos->first));
       }

       std::string type = d_value;
       std::multimap<int, std::string>::reverse_iterator
           rpos = type_index.rbegin();
       if (rpos != type_index.rend()) {
            // if ((rpos->second == "ATOMS") || ((rpos->first > 1) && (rpos->second != "unknown")))
            if (rpos->second != "unknown") type = rpos->second;
       }
       return type;
}

int parseString(const std::string &orginal_string, std::vector<std::string> &list)
{
       std::string op_string;
       String::RemoveWhiteSpace(orginal_string, op_string);
       int has_special_character = 0;
       for (unsigned int i = 0; i < op_string.size(); ++i) {
            if (op_string[i] == '(' || op_string[i] == ')' || op_string[i] == ',' || op_string[i] == '-') {
                 has_special_character = 1;
                 break;
            }
       }
       if (!has_special_character) {
            list.push_back(op_string);
            return 0;
       }

       int numlp = 0;
       int numrp = 0;
       int numd = 0;
       int firstl = -1;
       int firstr = -1;
       int firstc = -1;
      
       for (unsigned int i = 0; i < op_string.size(); ++i) {
            if (op_string[i] == '(') {
                 numlp++;
                 if (firstl < 0) firstl = i;
            } else if (op_string[i] == ')') {
                 numrp++;
                 if (firstr < 0) firstr = i;
            } else if (op_string[i] == ',') {
                 if (firstc < 0) firstc = i;
            } else if (op_string[i] == '-') {
                 numd++;
            }
       } 

       if (numlp != numrp) return 2;
       if (numlp && numrp && firstr < firstl) return 2;

       int choice = 0;
       if (firstl >= 0) {
            if (firstc >= 0 && firstc < firstl)
                 choice = 2;
            else choice = 1;
       } else if (firstc >= 0) {
            if (firstl >= 0 && firstl < firstc)
                 choice = 1;
            else choice = 2;
       } else if (numd == 1) choice = 3;

       std::vector<std::string> list1, list2;
       if (choice == 1) {
            std::vector<std::vector<std::string> > lists;
            lists.clear();
            separate_string_by_parenthesis(op_string, list1); 
            if (list1.empty()) return 2;
            for (unsigned int i = 0; i < list1.size(); ++i) {
                 list2.clear();
                 int error = parseString(list1[i], list2);
                 if (error) return error;
                 lists.push_back(list2);
            }
            if (lists.empty()) return 2;

            list2 = lists[lists.size() - 1];
            for (int i = (int) lists.size() - 2; i >= 0; i--) {
                 list1.clear();
                 for (unsigned int j = 0; j < lists[i].size(); ++j) {
                      for (unsigned int k = 0; k < list2.size(); ++k) {
                           std::string cs = lists[i][j] + " " + list2[k];
                           list1.push_back(cs);
                      }
                 }
                 list2 = list1;
            }
            for (unsigned int i = 0; i < list2.size(); ++i) {
                 list.push_back(list2[i]);
            }
            return 0;
       } else if (choice == 2) {
            separate_string_by_comma(op_string, list1);
            if (list1.empty()) return 3;
            for (unsigned int i = 0; i < list1.size(); ++i) {
                 list2.clear();
                 int error = parseString(list1[i], list2);
                 if (error) return error;
                 for (unsigned int j = 0; j < list2.size(); ++j) {
                      list.push_back(list2[j]);
                 }
            }
            return 0;
       } else if (choice == 3) {
            std::string::size_type pos = op_string.find("-");
            if (pos == std::string::npos) return 3;
            int start = atoi(op_string.substr(0, pos).c_str());
            int end = atoi(op_string.substr(pos + 1).c_str());
            if (end <= start) return 3;
            if (end == 0) return 3;
            for (int i = start; i <= end; ++i) {
                 list.push_back(String::IntToString(i));
            }
            return 0;
       } else return 4;
}

static void get_cif_mapping(const std::string& filename, std::map<std::string, std::string>& category_mapping,
                            std::map<std::string, std::map<std::string, std::string> >& item_mapping,
                            std::list<std::pair<std::string, std::string > >& rename_categories,
                            std::map<std::string, std::list<std::pair<std::string, std::string> > >& rename_items,
                            std::map<std::string, std::set<std::string> >& remove_items, const int& status)
{
       category_mapping.clear();
       item_mapping.clear();
       rename_categories.clear();
       rename_items.clear();
       remove_items.clear();

       struct stat stattext;
       if (stat(filename.c_str(), &stattext) != 0) return;

       CifFile* fobjR = new CifFile(READ_MODE, filename);
       if (!fobjR) return;

       Block &block = fobjR->GetBlock(fobjR->GetFirstBlockName());

       ISTable *t = getTablePtr(block, "master_internal_cif_mapping");
       if (t) {
            std::string item_v4, item_v5;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(item_v5, t, i, "item_v5");
                 get_value_clean(item_v4, t, i, "item_v4");
                 if (item_v5 == item_v4) continue;
                 process_one_cif_mapping(item_v5, item_v4, category_mapping, item_mapping, status);
            }
       }

       t = getTablePtr(block, "rename_category");
       if (t) {
            std::string category_v4, category_v5;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(category_v5, t, i, "category_v5");
                 get_value_clean(category_v4, t, i, "category_v4");
                 rename_categories.push_back(std::make_pair(category_v4, category_v5));
            }
       }

       t = getTablePtr(block, "rename_item");
       if (t) {
            std::string category, old_item, new_item;
            std::list<std::pair<std::string, std::string> > t_list;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(category, t, i, "category");
                 get_value_clean(old_item, t, i, "old_item");
                 get_value_clean(new_item, t, i, "new_item");
                 std::map<std::string, std::list<std::pair<std::string, std::string> > >::iterator mpos = rename_items.find(category);
                 if (mpos != rename_items.end())
                      mpos->second.push_back(std::make_pair(old_item, new_item));
                 else {
                      t_list.clear();
                      t_list.push_back(std::make_pair(old_item, new_item));
                      rename_items.insert(std::make_pair(category, t_list));
                 }
            }
       }

       t = getTablePtr(block, "remove_item");
       if (t) {
            std::string category, item;
            std::set<std::string> t_set;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(category, t, i, "category");
                 get_value_clean(item, t, i, "item");
                 std::map<std::string, std::set<std::string> >::iterator mpos = remove_items.find(category);
                 if (mpos != remove_items.end())
                      mpos->second.insert(item);
                 else {
                      t_set.clear();
                      t_set.insert(item);
                      remove_items.insert(std::make_pair(category, t_set));
                 }
            }
       }

       delete fobjR;
}

static void process_one_cif_mapping(const std::string& item_v5, const std::string& item_v4, std::map<std::string, std::string>& category_mapping,
                                    std::map<std::string, std::map<std::string, std::string> >& item_mapping, const int& status)
{
       std::vector<std::string> data1, data2;
       std::map<std::string, std::string> items;

       get_wordarray(data1, item_v5, ".");
       get_wordarray(data2, item_v4, ".");
       std::string old_category = "";
       std::string new_category = "";
       std::string category = "";
       std::string old_item = "";
       std::string new_item = "";

       if (data2[0] != data1[0]) {
            if (status == INTERNAL_TO_MASTER) {
                 old_category = data2[0].substr(1);
                 new_category = data1[0].substr(1);
            } else if (status == MASTER_TO_INTERNAL) {
                 old_category = data1[0].substr(1);
                 new_category = data2[0].substr(1);
            }
       }

       if (data2[1] != data1[1]) {
            if (status == INTERNAL_TO_MASTER) {
                 category = data2[0].substr(1);
                 old_item = data2[1];
                 new_item = data1[1];
            } else if (status == MASTER_TO_INTERNAL) {
                 category = data1[0].substr(1);
                 old_item = data1[1];
                 new_item = data2[1];
            }
       }

       if (!old_category.empty() && !new_category.empty()) {
            category_mapping.insert(std::make_pair(old_category, new_category));
       }

       if (!category.empty() && !old_item.empty() && !new_item.empty()) {
            std::map<std::string, std::map<std::string, std::string> >::iterator
                pos = item_mapping.find(category);
            if (pos == item_mapping.end()) {
                 items.clear();
                 items.insert(std::make_pair(old_item, new_item));
                 item_mapping.insert(std::make_pair(category, items));
            } else pos->second.insert(std::make_pair(old_item, new_item));
       }
}

static void _multiple_sort(const std::vector<std::vector<int> >& matrix, const std::vector<unsigned int>& start_index,
                           const unsigned int field_index, std::vector<int>& score)
{
       if (field_index == matrix[0].size()) return;

       int d = matrix[0].size() - field_index;

       std::multimap<int, unsigned int> mapping;
       mapping.clear();

       std::set<int> unique_values;
       unique_values.clear();

       for (unsigned int i = 0; i < start_index.size(); ++i) {
            mapping.insert(std::make_pair(matrix[start_index[i]][field_index],
                                     start_index[i]));
            unique_values.insert(matrix[start_index[i]][field_index]);
       }

       int count = 0;
       std::vector<unsigned int> index;

       std::set<int>::iterator spt;
       std::multimap<int, unsigned int>::iterator mpt;
       for (spt = unique_values.begin(); spt != unique_values.end(); ++spt) {
            std::pair<std::multimap<int, unsigned int>::iterator,
                 std::multimap<int, unsigned int>::iterator> range
                 = mapping.equal_range(*spt);

            index.clear();
            for (mpt = range.first; mpt != range.second; ++mpt) {
                 score[mpt->second] += count * start_index.size() * ((int) pow(10, (double) d));
                 index.push_back(mpt->second);
            }
            _multiple_sort(matrix, index, field_index + 1, score);
            count++;
       }
}

void get_chain_start_and_end(const std::set<unsigned int>& chain_break_set, unsigned int size,
                             std::vector<std::pair<unsigned int, unsigned int> >& pair_array)
{
       pair_array.clear();
       if (chain_break_set.empty()) {
            pair_array.push_back(std::make_pair(0, size - 1));
            return;
       }

       std::vector<unsigned int> chain_break_array;
       chain_break_array.clear();
       chain_break_array.reserve(chain_break_set.size());

       for (std::set<unsigned int>::const_iterator
            pos = chain_break_set.begin(); pos != chain_break_set.end(); ++pos) {
            chain_break_array.push_back(*pos);
       }

       pair_array.reserve(chain_break_array.size());
       pair_array.push_back(std::make_pair(0, chain_break_array[0]));
       for (unsigned int i = 0; i < chain_break_array.size() - 1; ++i) {
            pair_array.push_back(std::make_pair(chain_break_array[i]+1,
                                                chain_break_array[i+1]));
       }
       unsigned int i = chain_break_array.size() - 1;
       if (chain_break_array[i] < (size - 1)) {
            pair_array.push_back(std::make_pair(chain_break_array[i] + 1, size - 1));
       }
}

std::string get_next_id(const std::string& firstcolumn, const std::string& secondcolumn,
                        const std::string& thirdcolumn, std::set<std::string>& id_set)
{
       for (unsigned int i = 0; i < firstcolumn.size(); i++) {
            for (unsigned int j = 0; j < secondcolumn.size(); j++) {
                 for (unsigned int k = 0; k < thirdcolumn.size(); k++) {
                      std::string id = "";
                      id += firstcolumn[i];
                      id += secondcolumn[j];
                      id += thirdcolumn[k];
                      if (id_set.find(id) == id_set.end()) {
                           id_set.insert(id);
                           return id;
                      }
                 }
            }
       }

       id_set.insert("AA0");
       return "AA0";
}

static void separate_string_by_parenthesis(std::string &cifstring, std::vector<std::string> &list)
{
       list.clear();
       std::string cs;
       cs.clear();
       int count = 0;
       for (unsigned int i = 0; i < cifstring.size(); i++) {
            if (cifstring[i] == ' ' || cifstring[i] == '\t' ||
                cifstring[i] == '\n') continue;
            if (cifstring[i] == '(') {
                 if (count != 0) cs += cifstring[i];
                 count++;
            } else if (cifstring[i] == ')') {
                 count--;
                 if (count == 0) {
                      if (cs != "") list.push_back(cs);
                      cs.clear();
                 } else cs += cifstring[i];
            } else cs += cifstring[i];
       }
}

static void separate_string_by_comma(std::string &cifstring, std::vector<std::string> &list)
{
       list.clear();
       std::string cs;
       cs.clear();
       int count = 0;
       for (unsigned int i = 0; i < cifstring.size(); ++i) {
            if (cifstring[i] == ' ' || cifstring[i] == '\t' ||
                cifstring[i] == '\n') continue;
            if (cifstring[i] == ',' && count == 0) {
                 if (cs != "") list.push_back(cs);
                 cs.clear();
                 i++;
            }
            cs += cifstring[i];
            if (cifstring[i] == '(') count++;
            if (cifstring[i] == ')') count--;
       }
       if (cs != "") list.push_back(cs);
}
