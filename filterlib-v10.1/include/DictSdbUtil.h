/*
FILE:     DictSdbUtil.h
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
#ifndef _H_DICT_UTIL_USING_SDB_FILE_H_
#define _H_DICT_UTIL_USING_SDB_FILE_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "TableFile.h"

#define DIFFRN_GROUP_TYPE      1
#define EM_GROUP_TYPE          2
#define EXPTL_GROUP_TYPE       3
#define NMR_GROUP_TYPE         4
#define PHASING_GROUP_TYPE     5
#define REFINE_GROUP_TYPE      6
#define REFLN_GROUP_TYPE       7

class DictUtil
{
   private:
       typedef struct {
            std::string minimum;                   // allowed lower boundary: _item_range.minimum, _pdbx_item_range.minimum
            std::string maximum;                   // allowed upper boundary: _item_range.maximum, _pdbx_item_range.maximum
            std::string rangeTxt;
            bool include_equal_minimum;
            bool include_equal_maximum;
            int i_minimum;
            int i_maximum;
            double d_minimum;
            double d_maximum;
       } _ITEM_RANGE;
            
       typedef struct {
            bool is_key;                           // _category_key.name
            bool is_mandatory;                     // _item.mandatory_code
            bool is_pdbx_item;                     // _pdbx_item.name
            bool is_pdbx_mandatory;                // _pdbx_item.mandatory_code
            std::string code;                      // _item_type.code ( int, float etc. )
            std::string context_type;              // _pdbx_item_context.type ( WWPDB_LOCAL etc. )
            // key: item_range, pdbx_item_range
            std::map<std::string, _ITEM_RANGE> item_range_mapping;
            // concatenation of _item_enumeration.value & _pdbx_item_enumeration.value
            // key: uppercase enumeration value
            // value real enumeration value
            std::map<std::string, std::string> enumeration;
       } _ITEM_DEF;

       typedef struct {
            bool is_mandatory;                     // _category.mandatory_code
            std::string context_type;              // _pdbx_category_context.type
            // std::set<std::string> parent_groups;   // _category_group.id
            std::vector<std::string> items;
            std::set<std::string> key_items;
            std::set<std::string> mandatory_items;
            std::map<std::string, bool> mandatory_items_with_info;
            std::map<std::string, _ITEM_DEF> items_def_mapping; // key: item name
            std::map<std::string, std::string> uppercase_items_map; // key: upppercase item name
       } _CATEGORY_DEF;

       std::string _rcsbroot, _dict_sdb_path;
       bool _debug_flag;
       bool _read_flag;
       bool _read_suceessful;

       std::string _default_string;
       std::set<std::string> _default_set;
       std::vector<std::string> _default_vector;
       std::map<std::string, bool> _default_mandatory_key_map;
       std::map<std::string, std::string> _default_map;
       std::map<std::string, std::map<std::string, std::string> > _empty_parent_item_map;

       // _selected_groups.id
       std::set<std::string> _selected_groups;

       // _special_selected_categories.id
       std::set<std::string> _special_selected_categories;

       // _additional_key_items.cifitem
       std::set<std::string> _additional_v4_key_items;

       // _additional_mandatory_items.cifitem
       std::set<std::string> _additional_v4_mandatory_items;

       // _integer_key_items.name
       std::set<std::string> _integer_key_items;

       // key: _additional_child_parent_link.child_item
       // value: _additional_child_parent_link.parent_item
       std::map<std::string, std::string> _additional_v4_child_parent_mapping;

       // first: child category
       // second.first: child item
       // second.second.first: parent category
       // second.second.second: parent item
       std::map<std::string, std::map<std::string, std::map<std::string, std::string> > > _child_parant_mapping;

       // first: parent category
       // second.first: parent item
       // second.second.first: child category
       // second.second.second: child item
       std::map<std::string, std::map<std::string, std::map<std::string, std::string> > > _parant_child_mapping;

       // first: parent category
       // second: child categories
       std::map<std::string, std::set<std::string> > _parant_child_category_mapping;

       // first: order 0 - has no childs
       //              1 - has childs
       //              2 - has childs and grandchilds
       //              ....
       // second: category set
       std::map<int, std::set<std::string> > _parant_child_hierarchy_mapping;

       // key: code
       // value: regular expression
       std::map<std::string, std::string> _code_expression_mapping;

       // key: group_id
       // value: category_id
       std::map<std::string, std::set<std::string> > _dict_group_mapping;

       // key: category name
       std::map<std::string, _CATEGORY_DEF> _dict_category_mapping;

       std::vector<std::vector<std::string> > _single_valuable_item_categores;

       // _ignore_key_item_categories.name
       std::set<std::string> _ignore_key_item_categories;

       // map.first: cif     values from _pdbx_item_enumeration.name
       // map.second.first:  values from _pdbx_item_enumeration.value
       // map.second.second: values from _pdbx_item_enumeration.detail
       std::map<std::string, std::map<std::string, std::string> > _cross_checking_enumeration_mapping;

       void _init();
       void _read_dict_config(const std::string& binaryfile);
       void _read_single_item_as_set(Block& block, const std::string& category, const std::string& item, std::set<std::string>& val_set);
       void _read_two_items_as_map(Block& block, const std::string& category, const std::string& key_item, const std::string& value_item,
                                   std::map<std::string, std::string>& val_map);
       void _read_table(Block& block, const std::string& category, const std::vector<std::string>& items, const bool& skip_flag,
                        std::vector<std::vector<std::string> >& values);
       void _read_dict(const std::string& binarydictfile);
       void _read_enumeration(Block& block, const std::string& category, std::map<std::string, std::map<std::string, std::string> >& enumeration);
       void _read_item_range(Block& block, const std::string& category, std::map<std::string, std::vector<std::string> >& item_range);
       void _read_items(Block& block, std::set<std::string>& all_key_items, std::map<std::string, std::string>& pdbx_items,
                        std::map<std::string, std::string>& item_codes, std::map<std::string, std::string>& item_contexts,
                        std::map<std::string, std::map<std::string, std::string> >& enumeration, std::map<std::string, std::string>&
                        category_mandatory_code, std::map<std::string, std::string>& category_contexts,
                        std::map<std::string, std::vector<std::string> >& item_range, std::map<std::string, std::vector<std::string> >& pdbx_item_range);
       void _get_item_range(const std::vector<std::string>& values, _ITEM_RANGE& _item_range);
       void _update_selected_parent_child_relationship(const std::string& child_cifitem, const std::string& parent_cifitem);
       void _read_category_group(Block& block);
       void _update_parant_child_hierarchy(const std::set<std::string>& selected_categories);
       void _insert_parent_child_item_mapping(const std::string& childCat, const std::string& childItem, const std::string& parentCat, const std::string&
                                              parentItem, std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >& mapping);
       void _insert_parent_child_category_mapping(const std::string& parentCat, const std::string& childCat);
       void _updateCategoryItems(_CATEGORY_DEF& category);
       void _updateCategoryKeyItems(_CATEGORY_DEF& category);
       void _updateCategoryMandatoryItems(_CATEGORY_DEF& category);
       void _updateCategoryMandatoryItemsWithKeyInfo(_CATEGORY_DEF& category);
       bool _hasRangeDefinition(const std::string& catName, const std::string& itemName, const std::string& rangeType);
       const std::string& _getRangeBoundaryValue(const std::string& catName, const std::string& itemName, const std::string& rangeType,
                                                 const std::string& valueType) const;
       const std::string& _checkValueIsInRangeBoundaries(const std::string& catName, const std::string& itemName, const std::string& rangeType,
                                                         const std::string& dataValue);
       void _get_pdbx_item_linked_group_list(std::vector<std::vector<std::string> >& values);
       bool _checkFloatValueIsInRangeBoundaries(const _ITEM_RANGE& item_range, const double& value);
       bool _checkIntValueIsInRangeBoundaries(const _ITEM_RANGE& item_range, const int& value);
       void _print_dict();
   public:
       DictUtil();
       ~DictUtil();
       void clear();
       void setRCSBROOT(const std::string& path);
       void setDictSdbPath(const std::string& path); // Full path dictionary sdb file name: /wwpdb_da/da_top/reference/dict-v4.0/mmcif_pdbx_v5_next.sdb
       void setDebugFlag();
       bool Read();
       bool isPublicCategory(const std::string& catName);
       bool isDefinedCategory(const std::string& catName);
       bool isPublicItem(const std::string& catName, const std::string& itemName);
       bool isDefinedItem(const std::string& catName, const std::string& itemName, std::string& caseSensitiveItemName);
       const std::vector<std::string>& getItems(const std::string& catName);
       const std::set<std::string>& getKeyItems(const std::string& catName);
       const std::set<std::string>& getMandatoryItems(const std::string& catName);
       const std::map<std::string, bool>& getMandatoryItemsWithKeyInfo(const std::string& catName);
       bool isIntegerKeyItem(const std::string& itemName);
       const std::map<std::string, std::string>& getEnumeration(const std::string& catName, const std::string& itemName) const;
       const std::string& getTypeCode(const std::string& catName, const std::string& itemName) const;
       const std::string& getRegularExpression(const std::string& catName, const std::string& itemName) const;
       bool hasItemRange(const std::string& catName, const std::string& itemName);
       const std::string& getItemRange(const std::string& catName, const std::string& itemName, const std::string& valueType) const;
       const std::string& checkItemRange(const std::string& catName, const std::string& itemName, const std::string& dataValue);
       bool hasPdbxItemRange(const std::string& catName, const std::string& itemName);
       const std::string& getPdbxItemRange(const std::string& catName, const std::string& itemName, const std::string& valueType) const;
       const std::string& checkPdbxItemRange(const std::string& catName, const std::string& itemName, const std::string& dataValue);
       const std::map<std::string, std::map<std::string, std::string> >& getParentItemMapping(const std::string& catName) const;
       const std::map<std::string, std::string>& getChildItemMapping(const std::string& catName, const std::string& itemName) const;
       const std::set<std::string>& getChildCategories(const std::string& catName) const;
       const std::map<int, std::set<std::string> >& getHierarchyOrder() const;
       const std::vector<std::vector<std::string> >& getSingleValuableItemCategores() const;
       const std::set<std::string>& getIgnoreKeyItemCategores() const;
       void getChildCifItems(const std::string& parentCifItem, std::list<std::string>& childCifItems);
       void getChildCifItems(const std::set<std::string>& parentCifItems, std::map<std::string, std::list<std::string> >& childCifItemMap);
       const std::map<std::string, std::string>& getCrossCheckingEnumerationMapping(const std::string& catName, const std::string& itemName) const;
       void checkEnumerationAndRange();
};

#endif
