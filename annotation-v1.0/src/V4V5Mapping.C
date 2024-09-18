/*
FILE:     V4V5Mapping.C
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
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdexcept>

#include "V4V5Mapping.h"
#include "utillib.h"

V4V5Mapping::V4V5Mapping()
{
       _clear();
}

V4V5Mapping::~V4V5Mapping()
{
       _clear();
}

void V4V5Mapping::_clear()
{
       _message.clear();
       _mappings.clear();
       _rename_categories.clear();
       _delete_categories.clear();
       _delete_items.clear();
       _renumbering_items.clear();
       _replace_values.clear();
}

void V4V5Mapping::ReadMapping(const std::string& binaryfile)
{
       struct stat statbuf;
       if (stat(binaryfile.c_str(), &statbuf) != 0) return;

       CifFile* fobjR = new CifFile(READ_MODE, binaryfile);
       if (!fobjR) return;

       Block &block = fobjR->GetBlock(fobjR->GetFirstBlockName());

       std::string cs;
       std::vector<std::string> data;

       ISTable *t = getTablePtr(block, "mapping");
       if (t) {
            std::vector<std::pair<std::string, std::string> > t_vec;
            std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > tt_vec;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs, t, i, "item_v4");
                 get_wordarray(data, cs, ".");
                 t_vec.clear();
                 t_vec.push_back(std::make_pair(data[0].substr(1), data[1]));
                 get_value(cs, t, i, "item_v5");
                 get_wordarray(data, cs, ".");
                 std::string category = data[0].substr(1);
                 std::string item = data[1];
                 std::map<std::string, std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > >::iterator
                     mpos = _mappings.find(category);
                 if (mpos != _mappings.end()) {
                      unsigned int size = mpos->second.size() - 1;
                      if (mpos->second[size].first == item)
                           mpos->second[size].second.push_back(std::make_pair(t_vec[0].first, t_vec[0].second));
                      else mpos->second.push_back(std::make_pair(item, t_vec));
                 } else {
                      tt_vec.clear();
                      tt_vec.push_back(std::make_pair(item, t_vec));
                      _mappings.insert(std::make_pair(category, tt_vec));
                 }
            }
       }

       t = getTablePtr(block, "rename_category");
       if (t) {
            std::string cs1, cs2;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs1, t, i, "old_category");
                 get_value(cs2, t, i, "new_category");
                 _rename_categories.insert(std::make_pair(cs1, cs2));
            }
       }

       t = getTablePtr(block, "delete_category");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs, t, i, "name");
                 _delete_categories.push_back(cs);
            }
       }

       std::list<std::string> t_list;

       t = getTablePtr(block, "delete_item");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value(cs, t, i, "name");
                 get_wordarray(data, cs, ".");
                 std::string category = data[0].substr(1);
                 std::map<std::string, std::list<std::string> >::iterator mpos = _delete_items.find(category);
                 if (mpos != _delete_items.end())
                      mpos->second.push_back(data[1]);
                 else {
                      t_list.clear();
                      t_list.push_back(data[1]);
                      _delete_items.insert(std::make_pair(category, t_list));
                 }
            }
       }

       std::string category, item;

       t = getTablePtr(block, "renumbering");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(category, t, i, "category");
                 get_value_clean(item, t, i, "item");
                 std::map<std::string, std::list<std::string> >::iterator mpos = _renumbering_items.find(category);
                 if (mpos != _renumbering_items.end())
                      mpos->second.push_back(item);
                 else {
                      t_list.clear();
                      t_list.push_back(item);
                      _renumbering_items.insert(std::make_pair(category, t_list));
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
                 std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::iterator mpos = _replace_values.find(category);
                 if (mpos != _replace_values.end()) {
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
                      _replace_values.insert(std::make_pair(category, t_map));
                 }
            }
       }

       delete fobjR;
}

std::string V4V5Mapping::Update(Block& block, const int& mapping_type)
{
       for (std::map<std::string, std::string>::const_iterator mpos = _rename_categories.begin(); mpos != _rename_categories.end(); ++mpos) {
            ISTable *old_t = getTablePtr(block, mpos->first);
            if (!old_t) {
                 deleteTable(block, mpos->first);
                 continue;
            }
            ISTable *new_t = getTablePtr(block, mpos->second);
            if (new_t) {
                 // _message += "Both '" + mpos->first + "' and '" + mpos->second + "' categories existed.\n";
                 deleteTable(block, mpos->first);
            } else {
                 deleteTable(block, mpos->second);
                 block.RenameTable(mpos->first, mpos->second);
            }
       }

       for (std::map<std::string, std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > >::const_iterator
            mpos = _mappings.begin(); mpos != _mappings.end(); ++mpos) {
            if (mapping_type == EM_MAPPING_TYPE) _update_em_categories(block, mpos->first, mpos->second);
            else if (mapping_type == NMR_MAPPING_TYPE || mapping_type == XRAY_MAPPING_TYPE || mapping_type == OTHER_MAPPING_TYPE)
                 _update_categories(block, mpos->first, mpos->second);
       }

       for (std::list<std::string>::const_iterator pos = _delete_categories.begin(); pos != _delete_categories.end(); ++pos) {
            deleteTable(block, *pos);
       }

       for (std::map<std::string, std::list<std::string> >::const_iterator mpos = _delete_items.begin(); mpos != _delete_items.end(); ++mpos) {
            ISTable *t = getTablePtr(block, mpos->first);
            if (!t) {
                 deleteTable(block, mpos->first);
                 continue;
            }
            for (std::list<std::string>::const_iterator pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                 if (t->IsColumnPresent(*pos)) t->DeleteColumn(*pos);
            }
            block.WriteTable(t);
       }

       for (std::map<std::string, std::list<std::string> >::const_iterator mpos = _renumbering_items.begin(); mpos != _renumbering_items.end(); ++mpos) {
            ISTable *t = getTablePtr(block, mpos->first);
            if (!t) continue;

            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 for (std::list<std::string>::const_iterator lpos = mpos->second.begin(); lpos != mpos->second.end(); ++lpos) {
                      if (t->IsColumnPresent(*lpos)) t->UpdateCell(i, *lpos, String::IntToString(i + 1));
                 }
            }
            block.WriteTable(t);
       }

       std::string cs;
       for (std::map<std::string, std::map<std::string, std::map<std::string, std::string> > >::const_iterator mpos = _replace_values.begin();
            mpos != _replace_values.end(); ++mpos) {
            ISTable *t = getTablePtr(block, mpos->first);
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
            block.WriteTable(t);
       }

       return _message;
}

void V4V5Mapping::_update_em_categories(Block& block, const std::string& category, const std::vector<std::pair<std::string,
                                        std::vector<std::pair<std::string, std::string> > > >& items)
{
       // key: old category name
       // pair.first: ISTable
       // pair.second.first: old item name
       // pair.second.second: ISTable rowNo
       std::map<std::string, std::pair<ISTable*, std::pair<std::string, int> > > old_categories_mapping;
       old_categories_mapping.clear();
       for (unsigned int i = 0; i < items[0].second.size(); ++i) {
            ISTable *t = getTablePtr(block, items[0].second[i].first);
            if (!t) continue;
            old_categories_mapping.insert(std::make_pair(items[0].second[i].first, std::make_pair(t,
                                          std::make_pair(items[0].second[i].first, t->GetNumRows()))));
       }
       if (old_categories_mapping.empty()) return;

       std::string key, cs;
       std::vector<std::string> data;
       std::map<std::string, std::vector<std::string> > values; 
       values.clear();
       for (unsigned int i = 1; i < items.size(); ++i) {
            for (std::vector<std::pair<std::string, std::string> >::const_iterator pos = items[i].second.begin(); pos != items[i].second.end(); ++pos) {
                 std::map<std::string, std::pair<ISTable*, std::pair<std::string, int> > >::const_iterator mpos = old_categories_mapping.find(pos->first);
                 if (mpos == old_categories_mapping.end()) continue;
                 for (int j = 0; j < mpos->second.second.second; j++) {
                      get_value_clean(cs, mpos->second.first, j, pos->second);
                      if (cs.empty()) continue;
                      get_value_clean(key, mpos->second.first, j, mpos->second.second.first);
                      std::map<std::string, std::vector<std::string> >::iterator mmpos = values.find(key);
                      if (mmpos != values.end()) mmpos->second[i] = cs;
                      else {
                           data.clear();
                           for (unsigned int k = 0; k < items.size(); ++k) data.push_back("");
                           data[0] = key;
                           data[i] = cs;
                           values.insert(std::make_pair(key, data));
                      }
                 }
            }
       }
       if (values.empty()) return;

       ISTable* t = new ISTable(category);
       for (std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >::const_iterator
            pos = items.begin(); pos != items.end(); ++pos) {
            t->AddColumn(pos->first);
       }

       int row = 0;
       for (std::map<std::string, std::vector<std::string> >::const_iterator mpos = values.begin(); mpos != values.end(); ++mpos) {
            t->AddRow();
            for (unsigned int i = 0; i < items.size(); ++i) {
                 t->UpdateCell(row, items[i].first, mpos->second[i]);
            }
            row++;
       }
       block.WriteTable(t);
}

void V4V5Mapping::_update_categories(Block& block, const std::string& category, const std::vector<std::pair<std::string,
                                     std::vector<std::pair<std::string, std::string> > > >& items)
{
       ISTable *t = getTablePtr(block, items[0].second[0].first);
       if (!t) {
            deleteTable(block, items[0].second[0].first);
            return;
       }

       std::string cs1, cs2;
       int rowNo = t->GetNumRows();
       for (std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >::const_iterator
            ipos = items.begin(); ipos != items.end(); ++ipos) {
            for (std::vector<std::pair<std::string, std::string> >::const_iterator ppos = ipos->second.begin(); ppos != ipos->second.end(); ++ppos) {
                 if (ipos->first == ppos->second || !t->IsColumnPresent(ppos->second)) continue;

                 if (!t->IsColumnPresent(ipos->first)) {
                      t->RenameColumn(ppos->second, ipos->first);
                 } else {
                      for (int i = 0; i < rowNo; ++i) {
                           get_value_clean(cs1, t, i, ppos->second);
                           get_value_clean(cs2, t, i, ipos->first);
                           if (!cs2.empty()) continue;
                           if (!cs1.empty()) t->UpdateCell(i, ipos->first, cs1);
                      }
                      t->DeleteColumn(ppos->second);
                 }
            }
       }
       block.WriteTable(t);

       if (category != items[0].second[0].first) {
            ISTable *new_t = getTablePtr(block, category);
            if (new_t) {
                 // _message += "Both '" + items[0].second[0].first + "' and '" + category + "' categories existed.\n";
                 deleteTable(block, items[0].second[0].first);
            } else {
                 deleteTable(block, category);
                 block.RenameTable(items[0].second[0].first, category);
            }
       }
}
