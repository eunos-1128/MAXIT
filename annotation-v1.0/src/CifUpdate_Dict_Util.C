/*
FILE:     CifUpdate_Dict_Util.C
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
// #include "regex.h"
#include "utillib.h"

#define NS           10
#define NUM_SC       5

static const char *skip_categories[NUM_SC] = {
       "adit_ligand_info",
       "atom_site",
       "atom_site_anisotrop",
       "citation",
       "pdbx_struct_link"
};

void AnnotationObj::_cif_update_dictionary_compliance(Block& block)
{
       if (!_dictUtil.Read()) return;

       std::set<std::string> skip_category_set;
       skip_category_set.clear();
       for (int i = 0; i < NUM_SC; ++i) skip_category_set.insert(skip_categories[i]);

       const std::map<int, std::set<std::string> >& childParantHierarchyOrder = _dictUtil.getHierarchyOrder();
       for (std::map<int, std::set<std::string> >::const_iterator mpos = childParantHierarchyOrder.begin(); mpos != childParantHierarchyOrder.end(); ++mpos) {
            for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                 ISTable *t = getTablePtr(block, *spos);
                 if (!t) {
                      deleteTable(block, *spos);
                      continue;
                 } else if (_is_empty_table_with_dictionary_check(block, *spos, t)) continue;

                 t = _getTablePtr(block, *spos);
                 _cif_update_category(t);
                 _cif_update_mandatory_item(t);
                 block.WriteTable(t);

                 _cif_update_parent_categories(block, *spos, t);
            }
       }

       std::string cs;
       std::vector<std::string> TableNames;

       block.GetTableNames(TableNames);
       for (std::vector<std::string>::const_iterator tpos = TableNames.begin(); tpos != TableNames.end(); ++tpos) {
            if (skip_category_set.find(*tpos) != skip_category_set.end()) continue;

            ISTable *t = _getTablePtr(block, *tpos);
            if (!t) {
                 deleteTable(block, *tpos);
                 continue;
            }

            _cif_update_mandatory_item(t);

            const std::vector<std::string>& ColumnNames = t->GetColumnNames();
            for (std::vector<std::string>::const_iterator cpos = ColumnNames.begin(); cpos != ColumnNames.end(); ++cpos) {
                 const std::map<std::string, std::string>& enums = _dictUtil.getEnumeration(*tpos, *cpos);
                 _cif_update_enumeration_value(t, *cpos, enums);
                 std::string regular_expression = _dictUtil.getRegularExpression(*tpos, *cpos);
                 _cif_update_check_value(t, *cpos, regular_expression);
            }
            block.WriteTable(t);
       }

       std::vector<unsigned int> deleterows;
       const std::vector<std::vector<std::string> >& categories = _dictUtil.getSingleValuableItemCategores();
       for (std::vector<std::vector<std::string> >::const_iterator pos = categories.begin(); pos != categories.end(); ++pos) {
            ISTable *t = _getTablePtr(block, (*pos)[0]);
            if (!t) {
                 deleteTable(block, (*pos)[0]);
                 continue;
            }

            deleterows.clear();
            unsigned int rowNo = t->GetNumRows();
            for (unsigned int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, t, i, (*pos)[1]);
                 if (cs.empty()) deleterows.push_back(i);
            }
            if (!deleterows.empty()) t->DeleteRows(deleterows);
            if (t->GetNumRows() == 0)
                 deleteTable(block, (*pos)[0]);
            else block.WriteTable(t);
       }
}

bool AnnotationObj::_is_empty_table_with_dictionary_check(Block& block, const std::string& catName, ISTable* t)
{
       const std::set<std::string>& child_category_set = _dictUtil.getChildCategories(catName);
       for (std::set<std::string>::const_iterator spos = child_category_set.begin(); spos != child_category_set.end(); ++spos) {
            ISTable *child = getTablePtr(block, *spos);
            if (child) return false;
       }

       std::set<std::string> skip_item_set;
       skip_item_set.clear();
       skip_item_set.insert("entry_id");

       const std::set<std::string>& keyItems = _dictUtil.getKeyItems(catName);

       for (std::set<std::string>::const_iterator spos = keyItems.begin(); spos != keyItems.end(); ++spos) {
            if (_dictUtil.isIntegerKeyItem(*spos) || *spos == "pdbx_refine_id") skip_item_set.insert(*spos);
       }
       if (is_empty_table(t, skip_item_set)) {
            deleteTable(block, catName);
            return true;
       }
       return false;
}

void AnnotationObj::_cif_update_parent_categories(Block& block, const std::string& catName, ISTable* child)
{
       const std::map<std::string, std::map<std::string, std::string> >& pItemMapping = _dictUtil.getParentItemMapping(catName);
       if (pItemMapping.empty()) return;

       std::vector<std::string> missing_items;

       for (std::map<std::string, std::map<std::string, std::string> >::const_iterator mpos = pItemMapping.begin(); mpos != pItemMapping.end(); ++mpos) {
            for (std::map<std::string, std::string>::const_iterator mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                 bool new_table = false;
                 ISTable *t = _getTablePtr(block, mmpos->first);
                 if (!t) {
                      t = _newTablePtr(mmpos->first);
                      new_table = true;
                 }
                 missing_items.clear();
                 missing_items.push_back(mmpos->second);
                 check_missing_item(t, missing_items);
                 _update_table_index(child, mpos->first, t, mmpos->second);
                 if (new_table && t->GetNumRows() == 0) {
                      delete t;
                      deleteTable(block, mmpos->first);
                      continue;
                 }
                 block.WriteTable(t);
            }
       }
}

void AnnotationObj::_cif_update_mandatory_item(ISTable *t)
{
       if (_mandatory_updated_categories.find(t->GetName()) != _mandatory_updated_categories.end()) return;
       _mandatory_updated_categories.insert(t->GetName());
    
       const std::map<std::string, bool>& mandItemsNames = _dictUtil.getMandatoryItemsWithKeyInfo(t->GetName());
       if (mandItemsNames.empty()) return;

       int numRows = t->GetNumRows();

       std::vector<std::string> data;
       std::string cs;

       for (std::map<std::string, bool>::const_iterator mpos = mandItemsNames.begin(); mpos != mandItemsNames.end(); ++mpos) { 
            if (!t->IsColumnPresent(mpos->first)) {
                 data.clear();
                 for (int l = 0; l < numRows; ++l) data.push_back(_cif_update_get_default_mandatory_item_value(mpos->first, mpos->second, l));
                 t->AddColumn(mpos->first, data);
            } else {
                 for (int l = 0; l < numRows; ++l) {
                      get_value_clean(cs, t, l, mpos->first);
                      if (!cs.empty()) continue;
                      t->UpdateCell(l, mpos->first, _cif_update_get_default_mandatory_item_value(mpos->first, mpos->second, l));
                 }
            }
       }
}

void AnnotationObj::_cif_update_enumeration_value(ISTable *t, const std::string& itemName, const std::map<std::string, std::string>& enums)
{
       if (enums.empty()) return;

       int numRows = t->GetNumRows();

       std::map<std::string, std::string> enumValues_with_Hyphen;
       std::string cs;

       bool sign_flag = false;
       enumValues_with_Hyphen.clear();
       for (std::map<std::string, std::string>::const_iterator mpos = enums.begin(); mpos != enums.end(); ++mpos) {
            cs = mpos->first;
            RemoveNonAlnum(cs, sign_flag);
            if (enumValues_with_Hyphen.find(cs) != enumValues_with_Hyphen.end()) {
                 enumValues_with_Hyphen.clear();
                 break;
            }
            enumValues_with_Hyphen.insert(std::make_pair(cs, mpos->second));
       }
       if (enumValues_with_Hyphen.empty()) {
            sign_flag = true;
            for (std::map<std::string, std::string>::const_iterator mpos = enums.begin(); mpos != enums.end(); ++mpos) {
                 cs = mpos->first;
                 RemoveNonAlnum(cs, sign_flag);
                 if (enumValues_with_Hyphen.find(cs) != enumValues_with_Hyphen.end()) {
                      enumValues_with_Hyphen.clear();
                      break;
                 }
                 enumValues_with_Hyphen.insert(std::make_pair(cs, mpos->second));
            }
       }

       for (int l = 0; l < numRows; ++l) {
            get_value_clean_upper(cs, t, l, itemName);
            std::map<std::string, std::string>::const_iterator mpos = enums.find(cs);
            if (mpos != enums.end()) t->UpdateCell(l, itemName, mpos->second);
            else if (!enumValues_with_Hyphen.empty()) {
                 RemoveNonAlnum(cs, sign_flag);
                 mpos = enumValues_with_Hyphen.find(cs);
                 if (mpos != enumValues_with_Hyphen.end()) t->UpdateCell(l, itemName, mpos->second);
            }
       }
}

std::string AnnotationObj::_cif_update_get_default_mandatory_item_value(const std::string& itemName, const bool& keyFlag, const int& row)
{
       std::string value = ".";
       if (itemName == "entry_id")
            value = _StructureId;
       else if (itemName == "pdbx_refine_id")
            value = _default_exp_type;
       else if (keyFlag && _dictUtil.isIntegerKeyItem(itemName))
            value = String::IntToString(row + 1);
       return value;
}

void AnnotationObj::_cif_update_check_value(ISTable *t, const std::string& itemName, const std::string& regular_expression)
{
       // regex_t preg;
       // regmatch_t pmatch[NS];
       std::string cs;

       int numRows = t->GetNumRows();
       for (int l = 0; l < numRows; ++l) {
            get_value_clean(cs, t, l, itemName);
            if (cs.empty()) continue;
            if (IsNotValidatedValue(cs)) t->UpdateCell(l, itemName, ".");
/*
            else if (!regular_expression.empty()) {
                 rcsb_regcomp(&preg, regular_expression.c_str(), REG_EXTENDED);
                 int ret = rcsb_regexec(&preg, cs.c_str(), NS, pmatch, 0);
                 rcsb_regfree(&preg);
                 if (ret != 0) t->UpdateCell(l, itemName, ".");
            }
*/
       }
}
