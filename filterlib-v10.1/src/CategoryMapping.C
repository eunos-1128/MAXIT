/*
FILE:     CategoryMapping.C
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

#include "CategoryMapping.h"
#include "CategoryMapping_global.h"
#include "NdbToken.h"
#include "utillib.h"

void CategoryMapping::initialize()
{
       _all_categories.clear();
       _type_categories.clear();
       _depUI_remove_categories.clear();
       _skip_non_ascii_checking_categories.clear();
       _depUI_clear_items.clear();
}

bool CategoryMapping::Read(LogUtil& logIo, const std::string& path)
{
       std::string binaryfile = path + "/data/binary/";
       binaryfile += cifbinfile;
       struct stat statbuf;
       if (stat(binaryfile.c_str(), &statbuf) != 0) {
            std::string error = "CategoryMapping::Read: Can't find " + binaryfile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       CifFile* fobjR = new CifFile(READ_MODE, binaryfile);
       if (!fobjR) {
            std::string error = "CategoryMapping::Read: Can't open " + binaryfile + "\n";
            logIo.message(error.c_str());
            return false;
       }

       Block &block = fobjR->GetBlock(fobjR->GetFirstBlockName());

       bool found_error = false;
       if (!block.IsTablePresent("ndb_cif_category")) {
            std::string error = "CategoryMapping::Read: Can't find ndb_cif_category category in " + binaryfile + "\n";
            logIo.message(error.c_str());
            found_error = true;
       }
       if (!block.IsTablePresent("ndb_cif_keyword")) {
            std::string error = "CategoryMapping::Read: Can't find ndb_cif_keyword category in " + binaryfile + "\n";
            logIo.message(error.c_str());
            found_error = true;
       }
       if (found_error) {
            delete fobjR; return false;
       }

       read_ndb_cif_category(block);

       read_ndb_cif_keyword(block);

       delete fobjR;

       update_key_index();

       return true;
}	

bool CategoryMapping::IsValidCategory(const std::string& category)
{
       if (category.empty()) return false;

       std::map<std::string, CIF_CATEGORY>::const_iterator mpos = _all_categories.find(category);
       if (mpos != _all_categories.end()) return true;
       return false;
}

const CIF_CATEGORY& CategoryMapping::find_category(const std::string& category)
{
       std::map<std::string, CIF_CATEGORY>::const_iterator mpos = _all_categories.find(category);
       if (mpos == _all_categories.end()) {
            throw std::out_of_range("Category " + category + " is unrecognized.\n");
       }
       return mpos->second;
}

const std::map<std::string, CIF_CATEGORY>& CategoryMapping::CifCategory()
{
       return _all_categories;
}

void CategoryMapping::GetCategoryNameListByType(const std::string& type, std::list<std::string>& categoryNameList)
{
       categoryNameList.clear();

       std::map<std::string, std::list<std::string> >::const_iterator mpos = _type_categories.find(type);
       if (mpos != _type_categories.end()) categoryNameList = mpos->second;
}

const std::set<std::string>& CategoryMapping::depUI_remove_categories()
{
       return _depUI_remove_categories;
}

const std::map<std::string, std::set<std::string> >& CategoryMapping::depUI_clear_items()
{
       return _depUI_clear_items;
}

bool CategoryMapping::IsSkipNonAsciiCheckingCategory(const std::string& category)
{
       if (category.empty()) return false;
       if (_skip_non_ascii_checking_categories.find(category) != _skip_non_ascii_checking_categories.end()) return true;
       return false;
}

void CategoryMapping::read_ndb_cif_category(Block& block)
{
       ISTable *t = getTablePtr(block, "ndb_cif_category");

       std::list<std::string> t_list;

       CIF_CATEGORY cifcategory;
       cifcategory.CifKeyField.clear();
       cifcategory.keywords.clear();
       cifcategory.itemNames.clear();

       std::string cs;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value(cifcategory.category, t, i, "category");

            get_value(cs, t, i, "xyz");
            if (!cs.empty())
                 cifcategory.xyz = atoi(cs.c_str());
            else cifcategory.xyz = 0;

            get_value(cs, t, i, "sf");
            if (!cs.empty())
                 cifcategory.sf = atoi(cs.c_str());
            else cifcategory.sf = 0;

            get_value(cs, t, i, "correspondence");
            if (!cs.empty())
                 cifcategory.coordence = atoi(cs.c_str());
            else cifcategory.coordence = 0;

            get_value(cs, t, i, "deposition");
            if (!cs.empty())
                 cifcategory.deposition = atoi(cs.c_str());
            else cifcategory.deposition = 0;

            get_value(cs, t, i, "xray_or_nmr");
            if (!cs.empty())
                 cifcategory.xray_or_nmr = atoi(cs.c_str());
            else cifcategory.xray_or_nmr = 0;

            get_value(cs, t, i, "general");
            if (!cs.empty())
                 cifcategory.general = atoi(cs.c_str());
            else cifcategory.general = 0;

            get_value(cs, t, i, "transfer");
            if (!cs.empty())
                 cifcategory.transfer = atoi(cs.c_str());
            else cifcategory.transfer = 0;

            get_value(cs, t, i, "type");
            if (!cs.empty())
                 cifcategory.type = cs;
            else cifcategory.type.clear();

            get_value_clean_upper(cs, t, i, "remove_by_DepUI_Flag");
            if (cs == "Y") _depUI_remove_categories.insert(cifcategory.category);

            get_value_clean_upper(cs, t, i, "skip_non_ascii_checking");
            if (cs == "Y") _skip_non_ascii_checking_categories.insert(cifcategory.category);

            _all_categories.insert(std::make_pair(cifcategory.category, cifcategory));

            if (!cifcategory.type.empty()) {
                 std::map<std::string, std::list<std::string> >::iterator mpos = _type_categories.find(cifcategory.type);
                 if (mpos != _type_categories.end())
                      mpos->second.push_back(cifcategory.category);
                 else {
                      t_list.clear();
                      t_list.push_back(cifcategory.category);
                      _type_categories.insert(std::make_pair(cifcategory.type, t_list));
                 }
            }
       }
}

void CategoryMapping::read_ndb_cif_keyword(Block& block)
{
       ISTable *t = getTablePtr(block, "ndb_cif_keyword");

       CIF_KEYWORD cif_keyword;

       ndb_token_match token_match;
       token_match.NdbTokenKeyFieldNo.clear();
       std::string cs;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value(cs, t, i, "category");
            std::map<std::string, CIF_CATEGORY>::iterator mpos = _all_categories.find(cs);
            if (mpos == _all_categories.end()) continue;

            get_value(token_match.NdbTokenName, t, i, "ndb_token");
            get_value(token_match.JrnlTokenName, t, i, "special_type");
            get_value(cs, t, i, "ndb_field");
            if (!cs.empty())
                 token_match.NdbFieldNo = atoi(cs.c_str());
            else token_match.NdbFieldNo = 0;
     
            get_value(cs, t, i, "key");
            if (atoi(cs.c_str()) == 1)
                 token_match.IsKey = true;
            else token_match.IsKey = false;

            get_value(cs, t, i, "deposition");
            if (!cs.empty())
                 token_match.deposition = atoi(cs.c_str());
            else token_match.deposition = 0;

            bool clear_flag = false;
            get_value_clean_upper(cs, t, i, "clear_by_DepUI_Flag");
            if (cs == "Y") clear_flag = true;

            get_value(cs, t, i, "keyword");
            bool found = false;
            for (unsigned int j = 0; j < mpos->second.keywords.size(); ++j) {
                 if (mpos->second.keywords[j].name == cs) {
                      found = true;
                      mpos->second.keywords[j].NdbInfo.push_back(token_match);
                      break;
                 }
            }
            if (!found) {
                 cif_keyword.name = cs;
                 cif_keyword.NdbInfo.clear();
                 cif_keyword.NdbInfo.push_back(token_match);
                 mpos->second.keywords.push_back(cif_keyword);
            }

            if (!clear_flag) continue;

            std::map<std::string, std::set<std::string> >::iterator ipos = _depUI_clear_items.find(mpos->first);
            if (ipos != _depUI_clear_items.end()) ipos->second.insert(cs);
            else {
                 std::set<std::string> t_set;
                 t_set.clear();
                 t_set.insert(cs);
                 _depUI_clear_items.insert(std::make_pair(mpos->first, t_set));
            }
       }
}

void CategoryMapping::update_key_index()
{
       for (std::map<std::string, CIF_CATEGORY>::iterator mpos = _all_categories.begin(); mpos != _all_categories.end(); ++mpos) {
            mpos->second.CifKeyField.clear();
            mpos->second.itemNames.clear();
            mpos->second.itemNames.reserve(mpos->second.keywords.size());
            for (unsigned int i = 0; i < mpos->second.keywords.size(); ++i) {
                 mpos->second.itemNames.push_back(mpos->second.keywords[i].name);
                 bool is_key = false;
                 if (mpos->second.keywords[i].name != "entry_id") {
                      for (unsigned int j = 0; j < mpos->second.keywords[i].NdbInfo.size(); ++j) {
                           if (mpos->second.keywords[i].NdbInfo[j].IsKey) {
                                is_key = true;
                                break;
                           }
                      }
                 }
                 if (is_key) mpos->second.CifKeyField.push_back(i);
            }
       
            if (mpos->second.CifKeyField.empty()) continue;

            for (unsigned int i = 0; i < mpos->second.keywords.size(); ++i) {
                 for (unsigned int j = 0; j < mpos->second.keywords[i].NdbInfo.size(); ++j) {
                      for (unsigned int l = 0; l < mpos->second.CifKeyField.size(); ++l) {
                           int keyField = mpos->second.CifKeyField[l];

                           std::string keyTokenId = "";
                           int KeyFieldNo = -1;
                           for (unsigned int k = 0; k < mpos->second.keywords[keyField].NdbInfo.size(); ++k) {
                                if (mpos->second.keywords[keyField].NdbInfo[k].JrnlTokenName == mpos->second.keywords[i].NdbInfo[j].JrnlTokenName) {
                                     keyTokenId = mpos->second.keywords[keyField].NdbInfo[k].NdbTokenName;
                                     KeyFieldNo = mpos->second.keywords[keyField].NdbInfo[k].NdbFieldNo - 1;
                                     break;
                                }
                           }

                           std::string TokenId = mpos->second.keywords[i].NdbInfo[j].NdbTokenName;
                           if (TokenId == "JRNL") {
                                TokenId = mpos->second.keywords[i].NdbInfo[0].JrnlTokenName;
                                keyTokenId = mpos->second.keywords[keyField].NdbInfo[0].JrnlTokenName;
                                KeyFieldNo = mpos->second.keywords[keyField].NdbInfo[0].NdbFieldNo - 1;
                           }

                           if (TokenId.empty() || keyTokenId.empty() || KeyFieldNo < 0) {
                                mpos->second.keywords[i].NdbInfo[j].NdbTokenKeyFieldNo.push_back(0);
                                continue;
                           }

                           try {
                                const ndb_token_format& ndbformat = NdbToken::getTokenFormat(TokenId);
                                const ndb_token_format& ndbformat1 = NdbToken::getTokenFormat(keyTokenId);
                                for (int k = 0; k < ndbformat.NumField; ++k) {
                                     if (ndbformat.FieldList[k].FieldLabel == ndbformat1.FieldList[KeyFieldNo].FieldLabel) {
                                          mpos->second.keywords[i].NdbInfo[j].NdbTokenKeyFieldNo.push_back(k);
                                          break;
                                     }
                                }
                           } catch (const std::exception& exc) {}
                      }
                 }
            }
       }
}
