/*
FILE:     DictSdbUtil.C
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
#include <math.h>

#include "CifFile.h"
#include "DicFile.h"
#include "DictSdbUtil.h"
#include "utillib.h"

#define dictbinfile  "mmcif_pdbx.sdb"
#define dictname     "mmcif_pdbx.dic"
#define dictconfigfile "dict_config.odb"

DictUtil::DictUtil()
{
       clear();
}

DictUtil::~DictUtil()
{
       clear();
}

void DictUtil::clear()
{
       _rcsbroot.clear();
       _dict_sdb_path.clear();
       _debug_flag = false;
       _read_flag = false; 
       _read_suceessful = false;
       _init();
}

void DictUtil::setRCSBROOT(const std::string& path)
{
       _rcsbroot = path;
}

void DictUtil::setDictSdbPath(const std::string& path)
{
       _dict_sdb_path = path;
}

void DictUtil::setDebugFlag()
{
       _debug_flag = true;
}

bool DictUtil::Read()
{
       if (_read_flag) return _read_suceessful;

       _read_flag = true;

       struct stat stattext;
       std::string binaryfilepath = _rcsbroot + "/data/binary/";

       std::string binaryconfigfile = binaryfilepath + dictconfigfile;
       if (stat(binaryconfigfile.c_str(), &stattext) == 0) _read_dict_config(binaryconfigfile);

       if (!_dict_sdb_path.empty() && (stat(_dict_sdb_path.c_str(), &stattext) == 0)) _read_dict(_dict_sdb_path);
       else {
            std::string binarydictfile = binaryfilepath + dictbinfile;
            if (stat(binarydictfile.c_str(), &stattext) == 0) _read_dict(binarydictfile);
       }

       if (!_read_suceessful) _init();

       if (_debug_flag) _print_dict();

       return _read_suceessful;
}

bool DictUtil::isPublicCategory(const std::string& catName)
{
       std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            // if (cpos->second.context_type.find("LOCAL") != std::string::npos) return false;
            if ((cpos->second.context_type == "WWPDB_LOCAL") || (cpos->second.context_type == "WWPDB_DEPRECATED")) return false;
            return true;
       }
       return false;
}

bool DictUtil::isDefinedCategory(const std::string& catName)
{
       std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) return true;
       return false;
}

bool DictUtil::isPublicItem(const std::string& catName, const std::string& itemName)
{
       std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.find(itemName);
            if (ipos != cpos->second.items_def_mapping.end()) {
                 if (ipos->second.context_type.find("LOCAL") != std::string::npos) return false;
                 return true;
            }
            return false;
       }
       return false;
}

bool DictUtil::isDefinedItem(const std::string& catName, const std::string& itemName, std::string& caseSensitiveItemName)
{
       caseSensitiveItemName.clear();
       std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.find(itemName);
            if (ipos != cpos->second.items_def_mapping.end()) return true;

            std::string uppercase_value;
            String::UpperCase(itemName, uppercase_value);

            _updateCategoryItems(cpos->second);
            std::map<std::string, std::string>::const_iterator mpos = cpos->second.uppercase_items_map.find(uppercase_value);
            if (mpos != cpos->second.uppercase_items_map.end()) {
                 caseSensitiveItemName = mpos->second;
                 return true;
            }
       }
       return false;
}

const std::vector<std::string>& DictUtil::getItems(const std::string& catName)
{
       std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            _updateCategoryItems(cpos->second);
            return cpos->second.items;
       }
       return _default_vector;
}

const std::set<std::string>& DictUtil::getKeyItems(const std::string& catName)
{
       std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            _updateCategoryKeyItems(cpos->second);
            return cpos->second.key_items;
       }
       return _default_set;
}

const std::set<std::string>& DictUtil::getMandatoryItems(const std::string& catName)
{
       std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            _updateCategoryMandatoryItems(cpos->second);
            return cpos->second.mandatory_items;
       }
       return _default_set;
}

const std::map<std::string, bool>& DictUtil::getMandatoryItemsWithKeyInfo(const std::string& catName)
{
       std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            _updateCategoryMandatoryItemsWithKeyInfo(cpos->second);
            return cpos->second.mandatory_items_with_info;
       }
       return _default_mandatory_key_map;
}

bool DictUtil::isIntegerKeyItem(const std::string& itemName)
{
       return (_integer_key_items.find(itemName) != _integer_key_items.end());
}

const std::map<std::string, std::string>& DictUtil::getEnumeration(const std::string& catName, const std::string& itemName) const
{
       std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.find(itemName);
            if (ipos != cpos->second.items_def_mapping.end()) return ipos->second.enumeration;
       }
       return _default_map;
}

const std::string& DictUtil::getTypeCode(const std::string& catName, const std::string& itemName) const
{
       std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.find(itemName);
            if (ipos != cpos->second.items_def_mapping.end()) return ipos->second.code;
       }
       return _default_string;
}

const std::string& DictUtil::getRegularExpression(const std::string& catName, const std::string& itemName) const
{
       std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.find(itemName);
            if (ipos != cpos->second.items_def_mapping.end()) {
                 std::map<std::string, std::string>::const_iterator mpos = _code_expression_mapping.find(ipos->second.code);
                 if (mpos != _code_expression_mapping.end()) return mpos->second;
            }
       }
       return _default_string;
}

bool DictUtil::hasItemRange(const std::string& catName, const std::string& itemName)
{
       return _hasRangeDefinition(catName, itemName, "item_range");
}

const std::string& DictUtil::getItemRange(const std::string& catName, const std::string& itemName, const std::string& valueType) const
{
       return _getRangeBoundaryValue(catName, itemName, "item_range", valueType);
}

const std::string& DictUtil::checkItemRange(const std::string& catName, const std::string& itemName, const std::string& dataValue)
{
       return _checkValueIsInRangeBoundaries(catName, itemName, "item_range", dataValue);
}

bool DictUtil::hasPdbxItemRange(const std::string& catName, const std::string& itemName)
{
       return _hasRangeDefinition(catName, itemName, "pdbx_item_range");
}

const std::string& DictUtil::getPdbxItemRange(const std::string& catName, const std::string& itemName, const std::string& valueType) const
{
       return _getRangeBoundaryValue(catName, itemName, "pdbx_item_range", valueType);
}

const std::string& DictUtil::checkPdbxItemRange(const std::string& catName, const std::string& itemName, const std::string& dataValue)
{
       return _checkValueIsInRangeBoundaries(catName, itemName, "pdbx_item_range", dataValue);
}

const std::map<std::string, std::map<std::string, std::string> >& DictUtil::getParentItemMapping(const std::string& catName) const
{
       std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator mpos = _child_parant_mapping.find(catName);
       if (mpos != _child_parant_mapping.end()) return mpos->second;

       return _empty_parent_item_map;
}

const std::map<std::string, std::string>& DictUtil::getChildItemMapping(const std::string& catName, const std::string& itemName) const
{
       std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator mpos = _parant_child_mapping.find(catName);
       if (mpos != _child_parant_mapping.end()) {
            std::map<std::string, std::map<std::string, std::string> >::const_iterator mmpos = mpos->second.find(itemName);
            if (mmpos != mpos->second.end()) return mmpos->second;
       }

       return _default_map;
}

const std::set<std::string>& DictUtil::getChildCategories(const std::string& catName) const
{
       std::map<std::string, std::set<std::string> >::const_iterator mpos = _parant_child_category_mapping.find(catName);
       if (mpos != _parant_child_category_mapping.end()) return mpos->second;

       return _default_set;
}

const std::map<int, std::set<std::string> >& DictUtil::getHierarchyOrder() const
{
       return _parant_child_hierarchy_mapping;
}

const std::vector<std::vector<std::string> >& DictUtil::getSingleValuableItemCategores() const
{
       return _single_valuable_item_categores;
}

const std::set<std::string>& DictUtil::getIgnoreKeyItemCategores() const
{
       return _ignore_key_item_categories;
}

void DictUtil::getChildCifItems(const std::string& parentCifItem, std::list<std::string>& childCifItems)
{
       childCifItems.clear();

       std::vector<std::vector<std::string> > linked_group_list_values;
       _get_pdbx_item_linked_group_list(linked_group_list_values);

       for (std::vector<std::vector<std::string> >::const_iterator vpos = linked_group_list_values.begin(); vpos != linked_group_list_values.end(); ++vpos) {
            if ((*vpos)[3] == parentCifItem) childCifItems.push_back((*vpos)[2]);
       }
}

void DictUtil::getChildCifItems(const std::set<std::string>& parentCifItems, std::map<std::string, std::list<std::string> >& childCifItemMap)
{
       childCifItemMap.clear();

       std::list<std::string> childCifItems;

       std::vector<std::vector<std::string> > linked_group_list_values;
       _get_pdbx_item_linked_group_list(linked_group_list_values);

       for (std::vector<std::vector<std::string> >::const_iterator vpos = linked_group_list_values.begin(); vpos != linked_group_list_values.end(); ++vpos) {
            if (parentCifItems.find((*vpos)[3]) == parentCifItems.end()) continue;

            std::map<std::string, std::list<std::string> >::iterator mpos = childCifItemMap.find((*vpos)[3]);
            if (mpos != childCifItemMap.end()) mpos->second.push_back((*vpos)[2]);
            else {
                 childCifItems.clear();
                 childCifItems.push_back((*vpos)[2]);
                 childCifItemMap.insert(std::make_pair((*vpos)[3], childCifItems));
            }
       }
}

const std::map<std::string, std::string>& DictUtil::getCrossCheckingEnumerationMapping(const std::string& catName, const std::string& itemName) const
{
       std::string cifItem = "_" + catName + "." + itemName;
       std::map<std::string, std::map<std::string, std::string> >::const_iterator mpos = _cross_checking_enumeration_mapping.find(cifItem);
       if (mpos != _cross_checking_enumeration_mapping.end()) return mpos->second;

       return _default_map;
}

void DictUtil::checkEnumerationAndRange()
{
       for (std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.begin(); cpos != _dict_category_mapping.end(); ++cpos) {
            for (std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.begin();
                 ipos != cpos->second.items_def_mapping.end(); ++ipos) {
                 if (ipos->second.item_range_mapping.empty() /* && ipos->second.enumeration.empty() */ ) continue;

                 fprintf(stdout, "_%s.%s (code=%s):\n", cpos->first.c_str(), ipos->first.c_str(), ipos->second.code.c_str());
                 for (std::map<std::string, _ITEM_RANGE>::const_iterator rpos = ipos->second.item_range_mapping.begin();
                      rpos != ipos->second.item_range_mapping.end(); ++rpos) {
                      if (ipos->second.code == "float")
                           fprintf(stdout, "\t%s: minimum=%s (%.3f %d) maximum=%s (%.3f %d)\n", rpos->first.c_str(), rpos->second.minimum.c_str(),
                                   rpos->second.d_minimum, rpos->second.include_equal_minimum, rpos->second.maximum.c_str(),
                                   rpos->second.d_maximum, rpos->second.include_equal_maximum);
                      else fprintf(stdout, "\t%s: minimum=%s (%d %d) maximum=%s (%d %d)\n", rpos->first.c_str(), rpos->second.minimum.c_str(),
                                   rpos->second.i_minimum, rpos->second.include_equal_minimum, rpos->second.maximum.c_str(),
                                   rpos->second.i_maximum, rpos->second.include_equal_maximum);
                 }
/*
                 if (!ipos->second.enumeration.empty()) {
                      fprintf(stdout, "\tenumeration:");
                      for (std::map<std::string, std::string>::const_iterator epos = ipos->second.enumeration.begin();
                           epos != ipos->second.enumeration.end(); ++epos) {
                           fprintf(stdout, " %s,", epos->second.c_str());
                      }
                      fprintf(stdout, "\n");
                 }
*/
                 fprintf(stdout, "\n");
            }
       }
}

void DictUtil::_init()
{
       _default_string.clear();
       _default_set.clear();
       _default_vector.clear();
       _default_map.clear();
       _default_mandatory_key_map.clear();
       _empty_parent_item_map.clear();
       _selected_groups.clear();
       _special_selected_categories.clear();
       _additional_v4_key_items.clear();
       _additional_v4_mandatory_items.clear();
       _integer_key_items.clear();
       _additional_v4_child_parent_mapping.clear();
       _child_parant_mapping.clear();
       _parant_child_mapping.clear();
       _parant_child_category_mapping.clear();
       _parant_child_hierarchy_mapping.clear();
       _code_expression_mapping.clear();
       _dict_group_mapping.clear();
       _dict_category_mapping.clear();
       _single_valuable_item_categores.clear();
       _ignore_key_item_categories.clear();
       _cross_checking_enumeration_mapping.clear();
}

void DictUtil::_read_dict_config(const std::string& binaryfile)
{
       CifFile* fobj = new CifFile(READ_MODE, binaryfile);
       if (!fobj) return;

       Block &block = fobj->GetBlock(fobj->GetFirstBlockName());
       _read_single_item_as_set(block, "selected_groups", "id", _selected_groups);
       _read_single_item_as_set(block, "special_selected_categories", "id", _special_selected_categories);
       _read_single_item_as_set(block, "additional_key_items", "cifitem", _additional_v4_key_items);
       _read_single_item_as_set(block, "additional_mandatory_items", "cifitem", _additional_v4_mandatory_items);
       _read_single_item_as_set(block, "integer_key_items", "name", _integer_key_items);
       _read_two_items_as_map(block, "additional_child_parent_link", "child_item", "parent_item", _additional_v4_child_parent_mapping);

       std::vector<std::string> items;
       items.clear();
       items.push_back("category");
       items.push_back("item");
       _read_table(block, "single_valuable_item_category", items, true, _single_valuable_item_categores);

       _read_single_item_as_set(block, "ignore_key_item_categories", "name", _ignore_key_item_categories);

       delete fobj;
}

void DictUtil::_read_single_item_as_set(Block& block, const std::string& category, const std::string& item, std::set<std::string>& val_set)
{
       val_set.clear();

       std::vector<std::vector<std::string> > values;
       std::vector<std::string> items;
       items.clear();
       items.push_back(item);
       _read_table(block, category, items, true, values);

       for (std::vector<std::vector<std::string> >::const_iterator pos = values.begin(); pos != values.end(); ++pos) {
            val_set.insert((*pos)[0]);
       }
}

void DictUtil::_read_two_items_as_map(Block& block, const std::string& category, const std::string& key_item, const std::string& value_item,
                                         std::map<std::string, std::string>& val_map)
{
       val_map.clear();

       std::vector<std::vector<std::string> > values;
       std::vector<std::string> items;
       items.clear();
       items.push_back(key_item);
       items.push_back(value_item);
       _read_table(block, category, items, true, values);

       for (std::vector<std::vector<std::string> >::const_iterator pos = values.begin(); pos != values.end(); ++pos) {
            val_map.insert(std::make_pair((*pos)[0], (*pos)[1]));
       }
}

void DictUtil::_read_table(Block& block, const std::string& category, const std::vector<std::string>& items, const bool& skip_flag,
                              std::vector<std::vector<std::string> >& values)
{
       values.clear();

       ISTable *t = getTablePtr(block, category);
       if (!t) return;

       std::string cs;
       std::vector<std::string> data;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            data.clear();
            for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
                 get_value_clean(cs, t, i, *pos);
                 if (cs.empty()) {
                      if (_debug_flag) printf("empty value in _%s.%s at row %d.\n", category.c_str(), pos->c_str(), i + 1);
                      if (skip_flag) continue;
                 }
                 data.push_back(cs);
            }
            if (data.size() == items.size()) values.push_back(data);
       }
}

void DictUtil::_read_dict(const std::string& binarydictfile)
{
       DicFile* dictObj = new DicFile(READ_MODE, binarydictfile, false);
       if (!dictObj) return;

       _read_suceessful = true;
       Block& block = dictObj->GetBlock(dictObj->GetFirstBlockName());

       // key: _item_type_list.code
       // value: _item_type_list.construct
       _read_two_items_as_map(block, "item_type_list", "code", "construct", _code_expression_mapping);

       // _category_key.name
       std::set<std::string> all_key_items;
       _read_single_item_as_set(block, "category_key", "name", all_key_items);

       // key: _pdbx_item.name
       // value: _pdbx_item.mandatory_code
       std::map<std::string, std::string> pdbx_items;
       _read_two_items_as_map(block, "pdbx_item", "name", "mandatory_code", pdbx_items);

       // key: _item_type.name
       // value: _item_type.code
       std::map<std::string, std::string> item_codes, pdbx_item_codes;
       _read_two_items_as_map(block, "item_type", "name", "code", item_codes);
       _read_two_items_as_map(block, "pdbx_item_type", "name", "code", pdbx_item_codes);
       for (std::map<std::string, std::string>::const_iterator mpos = pdbx_item_codes.begin(); mpos != pdbx_item_codes.end(); ++mpos) {
            std::map<std::string, std::string>::iterator mpos1 = item_codes.find(mpos->first);
            if (mpos1 != item_codes.end())
                 mpos1->second = mpos->second;
            else item_codes.insert(std::make_pair(mpos->first, mpos->second));
       }

       // key: _pdbx_item_context.item_name
       // value: _pdbx_item_context.type
       std::map<std::string, std::string> item_contexts;
       _read_two_items_as_map(block, "pdbx_item_context", "item_name", "type", item_contexts);

       // key: _category.id
       // value: _category.mandatory_code
       std::map<std::string, std::string> category_mandatory_code;
       _read_two_items_as_map(block, "category", "id", "mandatory_code", category_mandatory_code);

       // key: _pdbx_category_context.category_id
       // value: _pdbx_category_context.type
       std::map<std::string, std::string> category_contexts;
       _read_two_items_as_map(block, "pdbx_category_context", "category_id", "type", category_contexts);

       // key: item_name
       // second key: Uppercase enumeration value
       // value: real enumeration value
       std::map<std::string, std::map<std::string, std::string> > enumeration;
       enumeration.clear();

       _read_enumeration(block, "item_enumeration", enumeration);
       _read_enumeration(block, "pdbx_item_enumeration", enumeration);

       // key: item_name
       // value[0]: minimum
       // value[1]: equal_minimum
       // value[2]: maximum
       // value[3]: equal_maximum
       std::map<std::string, std::vector<std::string> > item_range, pdbx_item_range;

       _read_item_range(block, "item_range", item_range);
       _read_item_range(block, "pdbx_item_range", pdbx_item_range);

       _read_items(block, all_key_items, pdbx_items, item_codes, item_contexts, enumeration,
                   category_mandatory_code, category_contexts, item_range, pdbx_item_range);

       for (std::map<std::string, std::string>::const_iterator mpos = _additional_v4_child_parent_mapping.begin();
            mpos != _additional_v4_child_parent_mapping.end(); ++mpos) {
            _update_selected_parent_child_relationship(mpos->first, mpos->second);
       }

       _read_category_group(block);

       std::vector<std::vector<std::string> > values;
       std::vector<std::string> items;
       items.clear();
       items.push_back("child_category_id");
       items.push_back("parent_category_id");
       items.push_back("child_name");
       items.push_back("parent_name");
       _read_table(block, "pdbx_item_linked_group_list", items, true, values);

       if (_debug_flag) {
            for (std::set<std::string>::const_iterator spos = all_key_items.begin(); spos != all_key_items.end(); ++spos) {
                 printf("Key item %s is not used.\n", spos->c_str());
            }
            for (std::map<std::string, std::string>::const_iterator mpos = pdbx_items.begin(); mpos != pdbx_items.end(); ++mpos) {
                 printf("pdbx_item %s mandatory_code %s is not used.\n", mpos->first.c_str(), mpos->second.c_str());
            }
            for (std::map<std::string, std::string>::const_iterator mpos = item_codes.begin(); mpos != item_codes.end(); ++mpos) {
                 printf("item_type %s code %s is not used.\n", mpos->first.c_str(), mpos->second.c_str());
            }
            for (std::map<std::string, std::string>::const_iterator mpos = item_contexts.begin(); mpos != item_contexts.end(); ++mpos) {
                 printf("pdbx_item_context %s type %s is not used.\n", mpos->first.c_str(), mpos->second.c_str());
            }
            for (std::map<std::string, std::map<std::string, std::string> >::const_iterator epos = enumeration.begin(); epos != enumeration.end(); ++epos) {
                 printf("item %s has following enumeration value(s) that are not used.\n", epos->first.c_str());
                 for (std::map<std::string, std::string>::const_iterator mpos = epos->second.begin(); mpos != epos->second.end(); ++mpos) {
                      printf("\t%s\n", mpos->second.c_str());
                 }
            }
            for (std::map<std::string, std::vector<std::string> >::const_iterator ipos = item_range.begin(); ipos != item_range.end(); ++ipos) {
                 printf("item_range %s ( '%s' '%s' '%s' '%s' ) is not used.\n", ipos->first.c_str(), ipos->second[0].c_str(),
                        ipos->second[1].c_str(), ipos->second[2].c_str(), ipos->second[3].c_str());
            }
            for (std::map<std::string, std::vector<std::string> >::const_iterator ipos = pdbx_item_range.begin(); ipos != pdbx_item_range.end(); ++ipos) {
                 printf("pdbx_item_range %s ( '%s' '%s' '%s' '%s' ) is not used.\n", ipos->first.c_str(), ipos->second[0].c_str(),
                        ipos->second[1].c_str(), ipos->second[2].c_str(), ipos->second[3].c_str());
            }
       }

       delete dictObj;

       std::set<std::string> selected_categories;
       selected_categories.clear();
       for (std::set<std::string>::const_iterator gpos = _selected_groups.begin(); gpos != _selected_groups.end(); ++gpos) {
            std::map<std::string, std::set<std::string> >::const_iterator mpos = _dict_group_mapping.find(*gpos);
            if (mpos == _dict_group_mapping.end()) continue;
            for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                 selected_categories.insert(*spos);
            }
       }

       for (std::set<std::string>::const_iterator spos = _special_selected_categories.begin(); spos != _special_selected_categories.end(); ++spos) {
            selected_categories.insert(*spos);
       }

       for (std::vector<std::vector<std::string> >::const_iterator vpos = values.begin(); vpos != values.end(); ++vpos) {
            if (selected_categories.find((*vpos)[0]) == selected_categories.end() || selected_categories.find((*vpos)[1]) == selected_categories.end()) continue;
            _update_selected_parent_child_relationship((*vpos)[2], (*vpos)[3]);
       }

       _update_parant_child_hierarchy(selected_categories);
}

void DictUtil::_read_enumeration(Block& block, const std::string& category, std::map<std::string, std::map<std::string, std::string> >& enumeration)
{
       std::vector<std::vector<std::string> > values;
       std::vector<std::string> items;
       items.clear();
       items.push_back("name");
       items.push_back("value");
       items.push_back("detail");
       _read_table(block, category, items, false, values);

       std::map<std::string, std::string> t_map;
       std::string uppercase_value;
       for (std::vector<std::vector<std::string> >::const_iterator pos = values.begin(); pos != values.end(); ++pos) {
            if ((*pos)[2] == "deprecated") continue;

            if ((category == "pdbx_item_enumeration") && !(*pos)[2].empty() && ((*pos)[2] != ".") && ((*pos)[2] != "?")) {
                 std::map<std::string, std::map<std::string, std::string> >::iterator mpos = _cross_checking_enumeration_mapping.find((*pos)[0]);
                 if (mpos != _cross_checking_enumeration_mapping.end()) mpos->second.insert(std::make_pair((*pos)[1], (*pos)[2]));
                 else {
                      t_map.clear();
                      t_map.insert(std::make_pair((*pos)[1], (*pos)[2]));
                      _cross_checking_enumeration_mapping.insert(std::make_pair((*pos)[0], t_map));
                 }
            }
            String::UpperCase((*pos)[1], uppercase_value);
            std::map<std::string, std::map<std::string, std::string> >::iterator mpos = enumeration.find((*pos)[0]);
            if (mpos != enumeration.end()) mpos->second.insert(std::make_pair(uppercase_value, (*pos)[1]));
            else {
                 t_map.clear();
                 t_map.insert(std::make_pair(uppercase_value, (*pos)[1]));
                 enumeration.insert(std::make_pair((*pos)[0], t_map));
            }
       }
}

void DictUtil::_read_item_range(Block& block, const std::string& category, std::map<std::string, std::vector<std::string> >& item_range)
{
       item_range.clear();

       std::vector<std::vector<std::string> > values;
       std::vector<std::string> items;
       items.clear();
       items.push_back("name");
       items.push_back("minimum");
       items.push_back("maximum");
       _read_table(block, category, items, false, values);

       std::map<std::string, std::vector<std::pair<std::string, std::string> > > range_mapping;
       range_mapping.clear();
       std::map<std::string, std::set<std::string> > equal_mapping;
       equal_mapping.clear();
       std::vector<std::pair<std::string, std::string> > tmp_vec;
       std::set<std::string> tmp_set;

       for (std::vector<std::vector<std::string> >::const_iterator pos = values.begin(); pos != values.end(); ++pos) {
            if ((*pos)[1] == (*pos)[2]) {
                 if ((*pos)[1].empty()) continue;
                 std::map<std::string, std::set<std::string> >::iterator smpos = equal_mapping.find((*pos)[0]);
                 if (smpos != equal_mapping.end()) smpos->second.insert((*pos)[1]);
                 else {
                      tmp_set.clear();
                      tmp_set.insert((*pos)[1]);
                      equal_mapping.insert(std::make_pair((*pos)[0], tmp_set));
                 }
            } else {
                 std::map<std::string, std::vector<std::pair<std::string, std::string> > >::iterator rmpos = range_mapping.find((*pos)[0]);
                 if (rmpos != range_mapping.end()) rmpos->second.push_back(std::make_pair((*pos)[1], (*pos)[2]));
                 else {
                      tmp_vec.clear();
                      tmp_vec.push_back(std::make_pair((*pos)[1], (*pos)[2]));
                      range_mapping.insert(std::make_pair((*pos)[0], tmp_vec));
                 }
            }
       }

       if (range_mapping.empty()) return;

       for (std::map<std::string, std::vector<std::pair<std::string, std::string> > >::const_iterator
            rpos = range_mapping.begin(); rpos != range_mapping.end(); ++rpos) {
            std::string minimum = rpos->second[0].first;
            std::string maximum = rpos->second[0].second;

            for (unsigned int i = 1; i < rpos->second.size(); ++i) {
                 bool has_same_minimum = false;
                 if (minimum == rpos->second[i].first) has_same_minimum = true;
                 int use_prev_minimum = 0;
                 if (minimum.empty()) {
                      use_prev_minimum = 1;
                      if (rpos->second[i].first.empty()) use_prev_minimum = 0;
                 } else {
                      if (rpos->second[i].first.empty()) use_prev_minimum = -1;
                      else {
                           double d_previous = String::StringToDouble(minimum);
                           double d_current = String::StringToDouble(rpos->second[i].first);
                           if (fabs(d_previous - d_current) < 0.0000001) {
                                has_same_minimum = true;
                                use_prev_minimum = 0;
                           } else if (d_previous > d_current) use_prev_minimum = -1;
                           else use_prev_minimum = 1;
                      }
                 }
                 bool has_same_maximum = false;
                 if (maximum == rpos->second[i].second) has_same_maximum = true;
                 int use_prev_maximum = 0;
                 if (maximum.empty()) {
                      use_prev_maximum = 1;
                      if (rpos->second[i].second.empty()) use_prev_maximum = 0;
                 } else {
                      if (rpos->second[i].second.empty()) use_prev_maximum = -1;
                      else {
                           double d_previous = String::StringToDouble(maximum);
                           double d_current = String::StringToDouble(rpos->second[i].second);
                           if (fabs(d_previous - d_current) < 0.0000001) {
                                has_same_maximum = true;
                                use_prev_maximum = 0;
                           } else if (d_previous < d_current) use_prev_maximum = -1;
                           else use_prev_maximum = 1;
                      }
                 }

                 if (use_prev_minimum == 0) {
                      if (use_prev_maximum == 0) {
                           if (has_same_minimum && has_same_maximum) continue;
                      } else if (use_prev_maximum < 0) {
                           continue;
                      } else {
                           minimum = rpos->second[i].first;
                           maximum = rpos->second[i].second;
                           continue;
                      }
                 } else if (use_prev_minimum < 0) {
                      if (use_prev_maximum <= 0) continue;
                 } else {
                      if (use_prev_maximum >= 0) {
                           minimum = rpos->second[i].first;
                           maximum = rpos->second[i].second;
                           continue;
                      }
                 }
                 fprintf(stdout, "Incorret range definitions in %s (%s): [ '%s', '%s'] vs [ '%s', '%s']\n", rpos->first.c_str(), category.c_str(),
                                 minimum.c_str(), maximum.c_str(), rpos->second[i].first.c_str(), rpos->second[i].second.c_str());
            }

            std::string equal_minimum = "";
            std::string equal_maximum = "";
            std::map<std::string, std::set<std::string> >::const_iterator epos = equal_mapping.find(rpos->first);
            if (epos != equal_mapping.end()) {
                 if (epos->second.find(minimum) != epos->second.end()) equal_minimum = minimum;
                 if (epos->second.find(maximum) != epos->second.end()) equal_maximum = maximum;
            }
            items.clear();
            items.push_back(minimum);
            items.push_back(equal_minimum);
            items.push_back(maximum);
            items.push_back(equal_maximum);
            item_range.insert(std::make_pair(rpos->first, items));
       }
}

void DictUtil::_read_items(Block& block, std::set<std::string>& all_key_items, std::map<std::string, std::string>& pdbx_items,
                              std::map<std::string, std::string>& item_codes, std::map<std::string, std::string>& item_contexts,
                              std::map<std::string, std::map<std::string, std::string> >& enumeration, std::map<std::string, std::string>& 
                              category_mandatory_code, std::map<std::string, std::string>& category_contexts,
                              std::map<std::string, std::vector<std::string> >& item_range,
                              std::map<std::string, std::vector<std::string> >& pdbx_item_range)
{
       std::vector<std::vector<std::string> > values;
       std::vector<std::string> items;
       items.clear();
       items.push_back("name");
       items.push_back("mandatory_code");
       _read_table(block, "item", items, false, values);

       _ITEM_DEF item_def;
       _CATEGORY_DEF category_def;
       std::string CatName, ItemName;
       for (std::vector<std::vector<std::string> >::const_iterator pos = values.begin(); pos != values.end(); ++pos) {
            item_def.is_key = false;
            if (all_key_items.find((*pos)[0]) != all_key_items.end()) {
                 item_def.is_key = true;
                 all_key_items.erase((*pos)[0]);
            } else if (_additional_v4_key_items.find((*pos)[0]) != _additional_v4_key_items.end()) item_def.is_key = true;

            item_def.is_mandatory = false;
            if (((*pos)[1] == "yes") || (_additional_v4_mandatory_items.find((*pos)[0]) != _additional_v4_mandatory_items.end())) item_def.is_mandatory = true;

            item_def.is_pdbx_item = false;
            item_def.is_pdbx_mandatory = false;
            std::map<std::string, std::string>::const_iterator mpos = pdbx_items.find((*pos)[0]);
            if (mpos != pdbx_items.end()) {
                 item_def.is_pdbx_item = true;
                 if (mpos->second == "yes") item_def.is_pdbx_mandatory = true;
                 pdbx_items.erase((*pos)[0]);
            }

            item_def.code.clear();
            mpos = item_codes.find((*pos)[0]);
            if (mpos != item_codes.end()) {
                 item_def.code = mpos->second;
                 item_codes.erase((*pos)[0]);
            } else if (_debug_flag) printf("No code defined for item '%s'.\n", (*pos)[0].c_str());

            item_def.context_type.clear();
            mpos = item_contexts.find((*pos)[0]);
            if (mpos != item_contexts.end()) {
                 item_def.context_type = mpos->second;
                 item_contexts.erase((*pos)[0]);
            }

            _ITEM_RANGE _item_range;
            item_def.item_range_mapping.clear();

            std::map<std::string, std::vector<std::string> >::const_iterator ipos = item_range.find((*pos)[0]);
            if (ipos != item_range.end()) {
                 _get_item_range(ipos->second, _item_range);
                 item_def.item_range_mapping.insert(std::make_pair("item_range", _item_range));

                 item_range.erase((*pos)[0]);
            }

            ipos = pdbx_item_range.find((*pos)[0]);
            if (ipos != pdbx_item_range.end()) {
                 _get_item_range(ipos->second, _item_range);
                 item_def.item_range_mapping.insert(std::make_pair("pdbx_item_range", _item_range));

                 pdbx_item_range.erase((*pos)[0]);
            }

            item_def.enumeration.clear();
            std::map<std::string, std::map<std::string, std::string> >::const_iterator epos = enumeration.find((*pos)[0]);
            if (epos != enumeration.end()) {
                 item_def.enumeration = epos->second;
                 enumeration.erase((*pos)[0]);
            }

            CifString::GetCategoryFromCifItem(CatName, (*pos)[0]);
            CifString::GetItemFromCifItem(ItemName, (*pos)[0]);

            std::map<std::string, _CATEGORY_DEF>::iterator dpos = _dict_category_mapping.find(CatName);
            if (dpos != _dict_category_mapping.end()) {
                 dpos->second.items_def_mapping.insert(std::make_pair(ItemName, item_def));
            } else {
                 category_def.is_mandatory = false;
                 mpos = category_mandatory_code.find(CatName);
                 if (mpos != category_mandatory_code.end()) {
                      if (mpos->second == "yes") category_def.is_mandatory = true;
                      category_mandatory_code.erase(CatName);
                 } else if (_debug_flag) printf("No mandatory_code defined for category '%s'.\n", CatName.c_str());

                 category_def.context_type.clear();
                 mpos = category_contexts.find(CatName);
                 if (mpos != category_contexts.end()) {
                      category_def.context_type = mpos->second;
                      category_contexts.erase(CatName);
                 }

                 // category_def.parent_groups.clear();
                 category_def.items.clear();
                 category_def.key_items.clear();
                 category_def.mandatory_items.clear();
                 category_def.mandatory_items_with_info.clear();
                 category_def.items_def_mapping.clear();
                 category_def.items_def_mapping.insert(std::make_pair(ItemName, item_def));
                 category_def.uppercase_items_map.clear();
                 _dict_category_mapping.insert(std::make_pair(CatName, category_def));
            }
       }

       if (!_debug_flag) return;

       for (std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.begin(); cpos != _dict_category_mapping.end(); ++cpos) {
            for (std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.begin();
                 ipos != cpos->second.items_def_mapping.end(); ++ipos) {
                 if (ipos->second.is_key) cpos->second.key_items.insert(ipos->first);
                 if (ipos->second.is_key || ipos->second.is_mandatory) cpos->second.mandatory_items.insert(ipos->first);
            }
            if (cpos->second.key_items.empty() && _debug_flag) printf("No key item defined for %s\n", cpos->first.c_str());
       }
}

void DictUtil::_get_item_range(const std::vector<std::string>& values, _ITEM_RANGE& _item_range)
{
       _item_range.minimum = values[0];
       _item_range.maximum = values[2];
       _item_range.rangeTxt = "[" + values[0] + ", " + values[2] + "]";

       _item_range.include_equal_minimum = false;
       if (!_item_range.minimum.empty() && (_item_range.minimum == values[1])) _item_range.include_equal_minimum = true;
       _item_range.include_equal_maximum = false;
       if (!_item_range.maximum.empty() && (_item_range.maximum == values[3])) _item_range.include_equal_maximum = true;

       _item_range.i_minimum = 0;
       _item_range.d_minimum = 0;
       if (!_item_range.minimum.empty()) {
            try {
                 _item_range.i_minimum = String::StringToInt(_item_range.minimum);
            } catch (const std::exception& exc) {
                 _item_range.i_minimum = 0;
            }
            try {
                 _item_range.d_minimum = String::StringToDouble(_item_range.minimum);
            } catch (const std::exception& exc) {
                 _item_range.d_minimum = 0;
            }
       }

       _item_range.i_maximum = 0;
       _item_range.d_maximum = 0;
       if (!_item_range.maximum.empty()) {
            try {
                 _item_range.i_maximum = String::StringToInt(_item_range.maximum);
            } catch (const std::exception& exc) {
                 _item_range.i_maximum = 0;
            }
            try {
                 _item_range.d_maximum = String::StringToDouble(_item_range.maximum);
            } catch (const std::exception& exc) {
                 _item_range.d_maximum = 0;
            }
       }
}

void DictUtil::_update_selected_parent_child_relationship(const std::string& child_cifitem, const std::string& parent_cifitem)
{
       // break em_imaging -> em_image_recording -> em_image_scans -> em_imaging loop
       if (child_cifitem == "_em_image_scans.image_recording_id") return;

       std::string childCat, childItem, parentCat, parentItem;
       CifString::GetCategoryFromCifItem(childCat, child_cifitem);
       CifString::GetItemFromCifItem(childItem, child_cifitem);

       std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.find(childCat);
       if (cpos == _dict_category_mapping.end()) return;

       std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.find(childItem);
       if (ipos == cpos->second.items_def_mapping.end()) return;
   
       CifString::GetCategoryFromCifItem(parentCat, parent_cifitem);
       CifString::GetItemFromCifItem(parentItem, parent_cifitem);

       cpos = _dict_category_mapping.find(parentCat);
       if (cpos == _dict_category_mapping.end()) return;

       ipos = cpos->second.items_def_mapping.find(parentItem);
       if (ipos == cpos->second.items_def_mapping.end()) return;

       _insert_parent_child_item_mapping(childCat, childItem, parentCat, parentItem, _child_parant_mapping);
       _insert_parent_child_item_mapping(parentCat, parentItem, childCat, childItem, _parant_child_mapping);
       _insert_parent_child_category_mapping(parentCat, childCat);
}

void DictUtil::_read_category_group(Block& block)
{
       std::vector<std::vector<std::string> > values;
       std::vector<std::string> items;
       items.clear();
       items.push_back("id");
       items.push_back("category_id");
       _read_table(block, "category_group", items, true, values);

       std::set<std::string> t_set;
       for (std::vector<std::vector<std::string> >::const_iterator pos = values.begin(); pos != values.end(); ++pos) {
            std::map<std::string, _CATEGORY_DEF>::iterator cpos = _dict_category_mapping.find((*pos)[1]);
            if (cpos == _dict_category_mapping.end()) {
                 if (_debug_flag) printf("_category_group.category_id %s can not be found in _dict_category_mapping.\n", (*pos)[1].c_str());
                 continue;
            }

            // cpos->second.parent_groups.insert((*pos)[0]);

            std::map<std::string, std::set<std::string> >::iterator mpos = _dict_group_mapping.find((*pos)[0]);
            if (mpos != _dict_group_mapping.end()) mpos->second.insert((*pos)[1]);
            else {
                 t_set.clear();
                 t_set.insert((*pos)[1]);
                 _dict_group_mapping.insert(std::make_pair((*pos)[0], t_set));
            }
       }
/*
       if (!_debug_flag) return;

       for (std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.begin(); cpos != _dict_category_mapping.end(); ++cpos) {
            if (cpos->second.parent_groups.empty()) printf("No group defined for category %s.\n", cpos->first.c_str());
       }
*/
}

void DictUtil::_update_parant_child_hierarchy(const std::set<std::string>& selected_categories)
{
       std::map<std::string, int> order_mapping;
       order_mapping.clear();
       for (std::set<std::string>::const_iterator spos = selected_categories.begin(); spos != selected_categories.end(); ++spos) {
            std::map<std::string, std::set<std::string> >::const_iterator mpos = _parant_child_category_mapping.find(*spos);
            if (mpos == _parant_child_category_mapping.end()) order_mapping.insert(std::make_pair(*spos, 0));
       }

       while (order_mapping.size() < selected_categories.size()) {
            for (std::set<std::string>::const_iterator spos = selected_categories.begin(); spos != selected_categories.end(); ++spos) {
                 std::map<std::string, int>::const_iterator opos = order_mapping.find(*spos);
                 if (opos != order_mapping.end()) continue;
                 std::map<std::string, std::set<std::string> >::const_iterator ppos = _parant_child_category_mapping.find(*spos);
                 if (ppos == _parant_child_category_mapping.end()) continue;
                 int order = 0;
                 bool found_missing_value = false;
                 for (std::set<std::string>::const_iterator cpos = ppos->second.begin(); cpos != ppos->second.end(); ++cpos) {
                      opos = order_mapping.find(*cpos);
                      if (opos == order_mapping.end()) {
                           found_missing_value = true;
                           break;
                      }
                      if (opos->second > order) order = opos->second;
                 }
                 if (found_missing_value) continue;
                 order_mapping.insert(std::make_pair(*spos, order + 1)); 
            }
       }

       std::set<std::string> t_set;
       for (std::map<std::string, int>::const_iterator mpos = order_mapping.begin(); mpos != order_mapping.end(); ++mpos) {
            std::map<int, std::set<std::string> >::iterator mmpos = _parant_child_hierarchy_mapping.find(mpos->second);
            if (mmpos != _parant_child_hierarchy_mapping.end())
                 mmpos->second.insert(mpos->first);
            else {
                 t_set.clear(); t_set.insert(mpos->first);
                 _parant_child_hierarchy_mapping.insert(std::make_pair(mpos->second, t_set));
            }
       }
}

void DictUtil::_insert_parent_child_item_mapping(const std::string& childCat, const std::string& childItem, const std::string& parentCat, const std::string&
                                                 parentItem, std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >& mapping)
{
       std::map<std::string, std::string> t_map;

       std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::iterator mpos = mapping.find(childCat);
       if (mpos != mapping.end()) {
            std::map<std::string, std::map<std::string, std::string> >::iterator mmpos = mpos->second.find(childItem);
            if (mmpos != mpos->second.end()) mmpos->second.insert(std::make_pair(parentCat, parentItem));
            else {
                 t_map.clear();
                 t_map.insert(std::make_pair(parentCat, parentItem));
                 mpos->second.insert(std::make_pair(childItem, t_map));
            }
       } else {
            t_map.clear();
            t_map.insert(std::make_pair(parentCat, parentItem));
            std::map<std::string, std::map<std::string, std::string> > tmp_map;
            tmp_map.clear();
            tmp_map.insert(std::make_pair(childItem, t_map));
            mapping.insert(std::make_pair(childCat, tmp_map));
       }
}

void DictUtil::_insert_parent_child_category_mapping(const std::string& parentCat, const std::string& childCat)
{
       std::set<std::string> t_set;

       std::map<std::string, std::set<std::string> >::iterator mpos = _parant_child_category_mapping.find(parentCat);
       if (mpos != _parant_child_category_mapping.end())
            mpos->second.insert(childCat);
       else {
            t_set.clear(); t_set.insert(childCat);
            _parant_child_category_mapping.insert(std::make_pair(parentCat, t_set));
       }
}

void DictUtil::_updateCategoryItems(_CATEGORY_DEF& category)
{
       if (!category.items.empty()) return;

       std::string uppercase_value;
       for (std::map<std::string, _ITEM_DEF>::const_iterator mpos = category.items_def_mapping.begin(); mpos != category.items_def_mapping.end(); ++mpos) {
            category.items.push_back(mpos->first);
            String::UpperCase(mpos->first, uppercase_value);
            category.uppercase_items_map.insert(std::make_pair(uppercase_value, mpos->first));
       }
}

void DictUtil::_updateCategoryKeyItems(_CATEGORY_DEF& category)
{
       if (!category.key_items.empty()) return;

       for (std::map<std::string, _ITEM_DEF>::const_iterator mpos = category.items_def_mapping.begin(); mpos != category.items_def_mapping.end(); ++mpos) {
            if (mpos->second.is_key) category.key_items.insert(mpos->first);
       }
}

void DictUtil::_updateCategoryMandatoryItems(_CATEGORY_DEF& category)
{
       if (!category.mandatory_items.empty()) return;

       for (std::map<std::string, _ITEM_DEF>::const_iterator mpos = category.items_def_mapping.begin(); mpos != category.items_def_mapping.end(); ++mpos) {
            if (mpos->second.is_key || mpos->second.is_mandatory) category.mandatory_items.insert(mpos->first);
       }
}

void DictUtil::_updateCategoryMandatoryItemsWithKeyInfo(_CATEGORY_DEF& category)
{
       if (!category.mandatory_items_with_info.empty()) return;

       for (std::map<std::string, _ITEM_DEF>::const_iterator mpos = category.items_def_mapping.begin(); mpos != category.items_def_mapping.end(); ++mpos) {
            if (mpos->second.is_key) category.mandatory_items_with_info.insert(std::make_pair(mpos->first, true));
            else if (mpos->second.is_mandatory) category.mandatory_items_with_info.insert(std::make_pair(mpos->first, false));
       }
}

bool DictUtil::_hasRangeDefinition(const std::string& catName, const std::string& itemName, const std::string& rangeType)
{
       std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.find(itemName);
            if (ipos != cpos->second.items_def_mapping.end()) {
                 std::map<std::string, _ITEM_RANGE>::const_iterator rpos = ipos->second.item_range_mapping.find(rangeType);
                 if (rpos != ipos->second.item_range_mapping.end()) return true;
            }
       }
       return false;
}

const std::string& DictUtil::_getRangeBoundaryValue(const std::string& catName, const std::string& itemName, const std::string& rangeType,
                                                    const std::string& valueType) const
{
       std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.find(itemName);
            if (ipos != cpos->second.items_def_mapping.end()) {
                 std::map<std::string, _ITEM_RANGE>::const_iterator rpos = ipos->second.item_range_mapping.find(rangeType);
                 if (rpos != ipos->second.item_range_mapping.end()) {
                      if (valueType == "minimum") return rpos->second.minimum;
                      else if (valueType == "maximum") return rpos->second.maximum;
                 }
            }
       }
       return _default_string;
}

const std::string& DictUtil::_checkValueIsInRangeBoundaries(const std::string& catName, const std::string& itemName, const std::string& rangeType,
                                                            const std::string& dataValue)
{
       if (dataValue.empty()) return _default_string;

       std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.find(catName);
       if (cpos != _dict_category_mapping.end()) {
            std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.find(itemName);
            if (ipos != cpos->second.items_def_mapping.end()) {
                 std::map<std::string, _ITEM_RANGE>::const_iterator rpos = ipos->second.item_range_mapping.find(rangeType);
                 if (rpos != ipos->second.item_range_mapping.end()) {
                      bool out_of_range_flag = false;
                      if (ipos->second.code == "float") {
                           try {
                                double value = String::StringToDouble(dataValue);
                                out_of_range_flag = _checkFloatValueIsInRangeBoundaries(rpos->second, value);
                           } catch (const std::exception& exc) {}
                      } else if ((ipos->second.code == "int") || (ipos->second.code == "positive_int")) {
                           try {
                                int value = String::StringToInt(dataValue);
                                out_of_range_flag = _checkIntValueIsInRangeBoundaries(rpos->second, value);
                           } catch (const std::exception& exc) {}
                      }
                      if (out_of_range_flag) return rpos->second.rangeTxt;
                 }
            }
       }

       return _default_string;
}

void DictUtil::_get_pdbx_item_linked_group_list(std::vector<std::vector<std::string> >& values)
{
       values.clear();

       struct stat stattext;
       std::string binaryfilepath = _rcsbroot + "/data/binary/";

       std::string binarydictfile = binaryfilepath + dictbinfile;
       if (stat(binarydictfile.c_str(), &stattext) != 0) return;

       DicFile* dictObj = new DicFile(READ_MODE, binarydictfile, false);
       if (!dictObj) return;

       Block& block = dictObj->GetBlock(dictObj->GetFirstBlockName());

       std::vector<std::string> items;
       items.clear();
       items.push_back("child_category_id");
       items.push_back("parent_category_id");
       items.push_back("child_name");
       items.push_back("parent_name");
       _read_table(block, "pdbx_item_linked_group_list", items, true, values);

       delete dictObj;
}

bool DictUtil::_checkFloatValueIsInRangeBoundaries(const _ITEM_RANGE& item_range, const double& value)
{
       bool out_of_range_flag = false;
       if (!item_range.minimum.empty()) {
            if (item_range.include_equal_minimum) {
                 if (value < item_range.d_minimum) out_of_range_flag = true;
            } else {
                 if (value <= item_range.d_minimum) out_of_range_flag = true;
            }
       }
       if (!item_range.maximum.empty()) {
            if (item_range.include_equal_maximum) {
                 if (value > item_range.d_maximum) out_of_range_flag = true;
            } else {
                 if (value >= item_range.d_maximum) out_of_range_flag = true;
            }
       }
       return out_of_range_flag;
}

bool DictUtil::_checkIntValueIsInRangeBoundaries(const _ITEM_RANGE& item_range, const int& value)
{
       bool out_of_range_flag = false;
       if (!item_range.minimum.empty()) {
            if (item_range.include_equal_minimum) {
                 if (value < item_range.i_minimum) out_of_range_flag = true;
            } else {
                 if (value <= item_range.i_minimum) out_of_range_flag = true;
            }
       }
       if (!item_range.maximum.empty()) {
            if (item_range.include_equal_maximum) {
                 if (value > item_range.i_maximum) out_of_range_flag = true;
            } else {
                 if (value >= item_range.i_maximum) out_of_range_flag = true;
            }
       }
       return out_of_range_flag;
}

void DictUtil::_print_dict()
{
       printf("Code/Regular expression map:\n");
       for (std::map<std::string, std::string>::const_iterator mpos = _code_expression_mapping.begin(); mpos != _code_expression_mapping.end(); ++mpos) {
            printf("\t%s <--> %s\n", mpos->first.c_str(), mpos->second.c_str());
       }

       printf("Category groups:\n");
       for (std::map<std::string, std::set<std::string> >::const_iterator mpos = _dict_group_mapping.begin(); mpos != _dict_group_mapping.end(); ++mpos) {
            printf("\tGroup %s:\n", mpos->first.c_str());
            for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos !=mpos->second.end(); ++spos) {
                 printf("\t\t%s\n", spos->c_str());
            }
            printf("\n");
       }

       for (std::map<std::string, _CATEGORY_DEF>::const_iterator cpos = _dict_category_mapping.begin(); cpos != _dict_category_mapping.end(); ++cpos) {
            printf("Category %s: %d %s\n", cpos->first.c_str(), cpos->second.is_mandatory, cpos->second.context_type.c_str());
            if (!cpos->second.key_items.empty()) {
                 printf("\tKey item(s):\n");
                 for (std::set<std::string>::const_iterator spos = cpos->second.key_items.begin(); spos != cpos->second.key_items.end(); ++spos) {
                      printf("\t\t%s\n", spos->c_str());
                 }
                 printf("\n");
            }
            if (!cpos->second.mandatory_items.empty()) {
                 printf("\tMandatory item(s):\n");
                 for (std::set<std::string>::const_iterator spos = cpos->second.mandatory_items.begin(); spos != cpos->second.mandatory_items.end(); ++spos) {
                      printf("\t\t%s\n", spos->c_str());
                 }
                 printf("\n");
            }
/*
            if (!cpos->second.parent_groups.empty()) {
                 printf("\tGroup(s):\n");
                 for (std::set<std::string>::const_iterator spos = cpos->second.parent_groups.begin(); spos != cpos->second.parent_groups.end(); ++spos) {
                      printf("\t\t%s\n", spos->c_str());
                 }
                 printf("\n");
            }
*/
            for (std::map<std::string, _ITEM_DEF>::const_iterator ipos = cpos->second.items_def_mapping.begin();
                 ipos != cpos->second.items_def_mapping.end(); ++ipos) {
                 printf("\tItem '%s': %d %d %d %d -%s- -%s-\n", ipos->first.c_str(), ipos->second.is_key, ipos->second.is_mandatory, ipos->second.is_pdbx_item,
                        ipos->second.is_pdbx_mandatory, ipos->second.code.c_str(), ipos->second.context_type.c_str());
                 if (ipos->second.enumeration.empty()) continue;
                 printf("\t\tEnumeration:\n");
                 for (std::map<std::string, std::string>::const_iterator mpos = ipos->second.enumeration.begin(); mpos != ipos->second.enumeration.end(); ++mpos) {
                      printf("\t\t%s\n", mpos->second.c_str());
                 }
                 printf("\n");
            }
            printf("\n\n");
       }

       printf("\n_child_parant_mapping:\n");
       for (std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator mpos = _child_parant_mapping.begin();
            mpos != _child_parant_mapping.end(); ++mpos) {
            printf("%s:\n", mpos->first.c_str());
            for (std::map<std::string, std::map<std::string, std::string> >::const_iterator mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                 printf("\t%s:\n", mmpos->first.c_str());
                 for (std::map<std::string, std::string>::const_iterator vpos = mmpos->second.begin(); vpos != mmpos->second.end(); ++vpos) {
                      printf("\t\t_%s.%s\n", vpos->first.c_str(), vpos->second.c_str());
                 } 
            }
       }

       printf("\n_parant_child_category_mapping:\n");
       for (std::map<std::string, std::set<std::string> >::const_iterator mpos = _parant_child_category_mapping.begin();
            mpos != _parant_child_category_mapping.end(); ++mpos) {
            printf("%s:\n", mpos->first.c_str());
            for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                 printf("\t%s\n", spos->c_str());
            }
       }

       printf("\n_parant_child_hierarchy_mapping:\n");
       for (std::map<int, std::set<std::string> >::const_iterator mpos = _parant_child_hierarchy_mapping.begin();
            mpos != _parant_child_hierarchy_mapping.end(); ++mpos) {
            printf("%d:\n", mpos->first);
            for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                 printf("\t%s\n", spos->c_str());
            }
       }
}
