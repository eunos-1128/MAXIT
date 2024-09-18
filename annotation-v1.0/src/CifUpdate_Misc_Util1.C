/*
FILE:     CifUpdate_Misc_Util1.C
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
#include <sys/stat.h>

#include "AnnotationObj.h"
#include "CategoryMapping.h"
#include "CompositeIndex.h"
#include "SeqCodeUtil.h"
#include "utillib.h"

static bool is_defualt_setting_diffrn_radiation(ISTable *t);
static void get_program_name_and_version(const std::string& val, std::string& program, std::string& version);
static bool is_same(const std::string& name1, const std::string& name2);
#if 0
static void read_cleanup_support_file(const std::string& filename, std::map<std::string, std::list<std::string> >& renumbering_items,
                                      std::list<std::string>& remove_categories, std::map<std::string, std::list<std::string> >& remove_items,
                                      std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >& replace_values);

void AnnotationObj::V4_V5_Cleanup()
{
       std::string cleanup_support_file = _rcsbroot + "/data/binary/v4_v5_clean_up.odb";
       struct stat stattext;
       if (stat(cleanup_support_file.c_str(), &stattext) != 0) return;

       std::list<std::string> remove_categories;
       // key: category name
       // values: item list
       std::map<std::string, std::list<std::string> > renumbering_items, remove_items;
       // primary key: category name
       // second key: item name
       // value: old vs. new value mapping
       std::map<std::string, std::map<std::string, std::map<std::string, std::string> > > replace_values;

       read_cleanup_support_file(cleanup_support_file, renumbering_items, remove_categories, remove_items, replace_values);

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       if (!renumbering_items.empty()) {
            for (std::map<std::string, std::list<std::string> >::const_iterator mpos = renumbering_items.begin(); mpos != renumbering_items.end(); ++mpos) {
                 ISTable *t = getTablePtr(_cifblock, mpos->first);
                 if (!t) continue;

                 int rowNo = t->GetNumRows();
                 for (int i = 0; i < rowNo; ++i) {
                      for (std::list<std::string>::const_iterator lpos = mpos->second.begin(); lpos != mpos->second.end(); ++lpos) {
                           if (t->IsColumnPresent(*lpos)) t->UpdateCell(i, *lpos, String::IntToString(i + 1));
                      }
                 }
                 _cifblock.WriteTable(t);
            }
       }

       if (!remove_categories.empty()) {
            for (std::list<std::string>::const_iterator lpos = remove_categories.begin(); lpos != remove_categories.end(); ++lpos) {
                 deleteTable(_cifblock, *lpos);
            }
       }

       if (!remove_items.empty()) {
            for (std::map<std::string, std::list<std::string> >::const_iterator mpos = remove_items.begin(); mpos != remove_items.end(); ++mpos) {
                 ISTable *t = getTablePtr(_cifblock, mpos->first);
                 if (!t) continue;

                 for (std::list<std::string>::const_iterator lpos = mpos->second.begin(); lpos != mpos->second.end(); ++lpos) {
                      if (t->IsColumnPresent(*lpos)) t->DeleteColumn(*lpos);
                 }
                 _cifblock.WriteTable(t);
            }
       }

       if (!replace_values.empty()) {
            std::string cs;
            for (std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator mpos = replace_values.begin();
                 mpos != replace_values.end(); ++mpos) {
                 ISTable *t = getTablePtr(_cifblock, mpos->first);
                 if (!t) continue;

                 int rowNo = t->GetNumRows();
                 for (int i = 0; i < rowNo; ++i) {
                      for (std::map<std::string, std::map<std::string, std::string> >::const_iterator mmpos = mpos->second.begin();
                           mmpos != mpos->second.end(); ++mmpos) {
                           get_value_clean(cs, t, i, mmpos->first);
                           std::map<std::string, std::string>::const_iterator mmmpos = mmpos->second.find(cs);
                           if (mmmpos != mmpos->second.end())  t->UpdateCell(i, mmpos->first, mmmpos->second);
                      }
                 }
                 _cifblock.WriteTable(t);
            }
       }
}
#endif

void AnnotationObj::_cif_update_remove_exp_spcific_categories(Block &block)
{
       if (_experiment_type == 0 || _experiment_type == EXPERIMENT_TYPE_BASIC) return;

       std::set<std::string> remove_exp_type_set;
       remove_exp_type_set.clear();
       remove_exp_type_set.insert("XRAY");
       remove_exp_type_set.insert("NMR");
       remove_exp_type_set.insert("EM");
       unsigned int size = remove_exp_type_set.size();
   

       if ((_experiment_type & EXPERIMENT_TYPE_NMR) || (_experiment_type & EXPERIMENT_TYPE_NMR_SOLID)) remove_exp_type_set.erase("NMR");
       if ((_experiment_type & EXPERIMENT_TYPE_EM) || (_experiment_type & EXPERIMENT_TYPE_EC) || (_experiment_type & EXPERIMENT_TYPE_ET))
            remove_exp_type_set.erase("EM");
       if ((_experiment_type & EXPERIMENT_TYPE_XRAY) || (_experiment_type & EXPERIMENT_TYPE_NEUTRON) || (_experiment_type & EXPERIMENT_TYPE_FIBER) ||
           ( _experiment_type & EXPERIMENT_TYPE_ELECTRON)) remove_exp_type_set.erase("XRAY");

       if (remove_exp_type_set.empty() || (remove_exp_type_set.size() == size)) return;

       std::list<std::string> cat_list;
       std::set<std::string> ignore_items;

       for (std::set<std::string>::const_iterator spos = remove_exp_type_set.begin(); spos != remove_exp_type_set.end(); ++spos) {
            CategoryMapping::GetCategoryNameListByType(*spos, cat_list);
            if (cat_list.empty()) continue;
            for (std::list<std::string>::const_iterator lpos = cat_list.begin(); lpos != cat_list.end(); ++lpos) {
                 ISTable *t = getTablePtr(block, *lpos);
                 if (!t) deleteTable(block, *lpos);
                 else if (*lpos == "diffrn_radiation") {
                      if (is_defualt_setting_diffrn_radiation(t)) deleteTable(block, *lpos);
                 } else if (*lpos == "refine_ls_restr") {
                      ignore_items.clear();
                      ignore_items.insert("type");
                      if (is_empty_table(t, ignore_items)) deleteTable(block, *lpos);
                 }
            }
       }
}

void AnnotationObj::_cif_update_software(Block &block)
{
       ISTable *t = getTablePtr(block, "computing");
       ISTable *t1 = _getTablePtr(block, "software");

       std::string cs;
       if (t1) {
            int rowNo = t1->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t1, i, "version");
                 if (cs.empty()) t1->UpdateCell(i, "version", ".");
                 t1->UpdateCell(i, "citation_id", "?");
            }
            block.WriteTable(t1);
            return;
       }

       if (!t) return;

       if (!t1) t1 = _newTablePtr("software");

       std::vector<std::string> classification, item;
       item.clear();
       item.push_back("data_collection");
       // item.push_back("data_reduction");
       item.push_back("pdbx_data_reduction_ii");
       item.push_back("pdbx_data_reduction_ds");
       item.push_back("structure_solution");
       item.push_back("structure_refinement");

       classification.clear();
       classification.push_back("data collection");
       classification.push_back("data reduction");
       classification.push_back("data scaling");
       // classification.push_back("model building");
       classification.push_back("phasing");
       classification.push_back("refinement");

       std::string cs1, cs2, cs3, program, version;
       std::vector<std::string> data;
       for (unsigned int i = 0; i < item.size(); ++i) {
            get_value(cs, t, 0, item[i]);
            if (cs.empty()) continue;
            std::string::size_type p = cs.find(" AND ");
            if (p != std::string::npos) cs.replace(p, 5, ",");
            p = cs.find(" and ");
            if (p != std::string::npos) cs.replace(p, 5, ",");
            p = cs.find(" And ");
            if (p != std::string::npos) cs.replace(p, 5, ",");
            get_wordarray(data, cs, ",;");
            for (std::vector<std::string>::iterator pos = data.begin(); pos != data.end(); ++pos) {
                 get_program_name_and_version(*pos, program, version);
                 bool found = false;
                 int rowNo = t1->GetNumRows();
                 for (int j = 0; j < rowNo; ++j) {
                      get_value_clean_lower(cs, t1, j, "classification");
                      get_value_clean_lower(cs1, t1, j, "name");
                      get_value_clean_lower(cs2, t1, j, "version");
                      String::LowerCase(program, cs3);
                      if (cs == classification[i] && is_same(cs1, cs3)) {
                           found = true;
                           if (cs2.empty()) t1->UpdateCell(j, "version", version);
                           break;
                      }
                 }
                 if (found) continue;

                 t1->AddRow();
                 t1->UpdateCell(rowNo, "name", program);
                 t1->UpdateCell(rowNo, "classification", classification[i]);
                 t1->UpdateCell(rowNo, "version", version);
            }
       }
       block.WriteTable(t1);
}

void AnnotationObj::_cif_update_repeat_row_categories(Block& block)
{
       std::list<std::string> categories;
       categories.clear();
       categories.push_back("database_PDB_rev_record");
       categories.push_back("pdbx_entity_name");

       std::string cs;
       std::vector<unsigned int> delete_rows;
       std::set<std::string> value_index;
       std::vector<std::string> data;

       for (std::list<std::string>::const_iterator lpos = categories.begin(); lpos != categories.end(); ++lpos) {
            ISTable *t = getTablePtr(block, *lpos);
            if (!t) {
                 deleteTable(block, *lpos);
                 continue;
            }

            const std::vector<std::string>& ColumnNames = t->GetColumnNames();
            delete_rows.clear();
            value_index.clear();
            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 data.clear();
                 for (std::vector<std::string>::const_iterator pos = ColumnNames.begin(); pos != ColumnNames.end(); ++pos) {
                      get_value_clean(cs, t, i, *pos);
                      data.push_back(cs); 
                 }
                 cs = CompositeIndex::getIndex(data);
                 if (value_index.find(cs) != value_index.end())
                      delete_rows.push_back(i);
                 else value_index.insert(cs);
            }
            if (!delete_rows.empty()) {
                 t->DeleteRows(delete_rows);
                 block.WriteTable(t);
            }
       }
}

void AnnotationObj::_cif_update_matthew_and_solvent(Block& block)
{
       if (!(_experiment_type & EXPERIMENT_TYPE_XRAY)) return;

       std::string matthew, solvent;
       matthew.clear();
       solvent.clear();

       ISTable *t = _getTablePtr(block, "exptl_crystal");
       if (t) {
            get_value_clean(matthew, t, 0, "density_Matthews");
            get_value_clean(solvent, t, 0, "density_percent_sol");
            if (!matthew.empty() && !solvent.empty()) return;
       }

       std::vector<double> return_values;
       if (!_cal_matthew_and_solvent(return_values)) return;

       double d_matthew = return_values[0];
       double d_solvent = return_values[1];
       d_solvent *= 100.0; 

       if (!t) {
            t = _newTablePtr("exptl_crystal");
            t->AddRow();
       }

       if (matthew.empty() && d_matthew >= return_values[2] && d_matthew <= return_values[3]) {
            t->UpdateCell(0, "density_Matthews", FloatToString(d_matthew, 0, 2));
       }

       if (solvent.empty() && d_solvent >= return_values[4] && d_solvent <= return_values[5]) {
            t->UpdateCell(0, "density_percent_sol", FloatToString(d_solvent, 0, 2));
       }

       block.WriteTable(t);
}

static bool is_defualt_setting_diffrn_radiation(ISTable *t)
{
       bool is_defualt_setting = true;

       if (!t) return is_defualt_setting;

       int rowNo = t->GetNumRows();

       if (rowNo == 0) return is_defualt_setting;

       std::string cs;
       const std::vector<std::string>& ColumnNames = t->GetColumnNames();
       
       for (int i = 0; i < rowNo; ++i) {
            for (std::vector<std::string>::const_iterator pos = ColumnNames.begin(); pos != ColumnNames.end(); ++pos) {
                 get_value_clean_upper(cs, t, i, *pos);
                 if (cs.empty()) continue;
                 else if (*pos == "diffrn_id") continue;
                 else if (*pos == "pdbx_monochromatic_or_laue_m_l" && cs == "M") continue;
                 else if (*pos == "pdbx_diffrn_protocol" && cs == "SINGLE WAVELENGTH") continue;
                 is_defualt_setting = false;
                 break;
            }
            if (!is_defualt_setting) break;
       }
       return is_defualt_setting;
}

static void get_program_name_and_version(const std::string& val, std::string& program, std::string& version)
{
       program.clear();
       version.clear();

       std::vector<std::string> data; 
       get_wordarray(data, val, " ");
       program = data[0];
       if (data.size() == 1) {
            version = ".";
            return;
       }

       for (unsigned int i = 1; i < data.size(); ++i) {
            if (!version.empty()) version += " ";
            version += data[i];
       }

       std::string cs;
       cs.clear();
       String::LowerCase(version, cs);
       if (cs.find("phenix") != std::string::npos) {
            replace_string(cs, "(", "");
            replace_string(cs, ")", "");
            get_wordarray(data, cs, " ");
            version.clear();
            for (unsigned int i = 0; i < data.size(); ++i) {
                 if (data[i].find("phenix") != std::string::npos) continue;
                 if (!version.empty()) version += " ";
                 version += data[i];
            }
       }
}

static bool is_same(const std::string& name1, const std::string& name2)
{
       if (String::IsEqual(name1, name2, Char::eCASE_INSENSITIVE)) return true;
       if ((String::IsEqual(name1, "D*trek", Char::eCASE_INSENSITIVE) && String::IsEqual(name2, "Dtrek", Char::eCASE_INSENSITIVE)) ||
           (String::IsEqual(name1, "Dtrek", Char::eCASE_INSENSITIVE) && String::IsEqual(name2, "D*trek", Char::eCASE_INSENSITIVE)))
            return true;

       return false;
}
#if 0
static void read_cleanup_support_file(const std::string& filename, std::map<std::string, std::list<std::string> >& renumbering_items,
                                      std::list<std::string>& remove_categories, std::map<std::string, std::list<std::string> >& remove_items,
                                      std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >& replace_values)
{
       renumbering_items.clear();
       remove_categories.clear();
       remove_items.clear();
       replace_values.clear();

       CifFile* fobj = new CifFile(READ_MODE, filename);
       if (!fobj) return;

       Block &block = fobj->GetBlock(fobj->GetFirstBlockName());

       std::string category, item;
       std::list<std::string> t_list;

       ISTable* t = getTablePtr(block, "renumbering");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(category, t, i, "category");
                 get_value_clean(item, t, i, "item");
                 std::map<std::string, std::list<std::string> >::iterator mpos = renumbering_items.find(category);
                 if (mpos != renumbering_items.end())
                      mpos->second.push_back(item);
                 else {
                      t_list.clear();
                      t_list.push_back(item);
                      renumbering_items.insert(std::make_pair(category, t_list));
                 }
            }
       }

       t = getTablePtr(block, "remove_category");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(category, t, i, "category");
                 remove_categories.push_back(category);
            }
       }

       t = getTablePtr(block, "remove_item");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(category, t, i, "category");
                 get_value_clean(item, t, i, "item");
                 std::map<std::string, std::list<std::string> >::iterator mpos = remove_items.find(category);
                 if (mpos != remove_items.end())
                      mpos->second.push_back(item);
                 else {
                      t_list.clear();
                      t_list.push_back(item);
                      remove_items.insert(std::make_pair(category, t_list));
                 }
            }
       }

       t = getTablePtr(block, "replace_value");
       if (t) {
            std::string old_value, new_value;
            std::map<std::string, std::map<std::string, std::string> > t_map;
            std::map<std::string, std::string> tt_map;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(category, t, i, "category");
                 get_value_clean(item, t, i, "item");
                 get_value_clean(old_value, t, i, "old_value");
                 get_value_clean(new_value, t, i, "new_value");
                 std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::iterator mpos = replace_values.find(category);
                 if (mpos != replace_values.end()) {
                      std::map<std::string, std::map<std::string, std::string> >::iterator mmpos = mpos->second.find(item);
                      if (mmpos != mpos->second.end()) {
                           mmpos->second.insert(std::make_pair(old_value, new_value));
                      } else {
                           tt_map.clear();
                           tt_map.insert(std::make_pair(old_value, new_value));
                           mpos->second.insert(std::make_pair(item, tt_map));
                      }
                 } else {
                      tt_map.clear();
                      tt_map.insert(std::make_pair(old_value, new_value));
                      t_map.clear();
                      t_map.insert(std::make_pair(item, tt_map));
                      replace_values.insert(std::make_pair(category, t_map));
                 }
            }
       }

       delete fobj;
}
#endif
