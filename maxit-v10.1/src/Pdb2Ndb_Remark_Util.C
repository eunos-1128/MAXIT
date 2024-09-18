/*
FILE:     Pdb2Ndb_Remark_Util.C
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
#include <ctype.h>

#include "Maxit.h"
#include "NdbToken.h"
#include "Remark_extern.h"
#include "utillib.h"

#define NUM_LOCAL_FIELD 10
#define NUM_SHELL_FIELD  9
#define NUM_STR_FORMAT   3

static const int _long_format[NUM_SHELL_FIELD][2] = {
       {  0, 5 },
       {  5, 8 },
       { 14, 1 },
       { 15, 8 },
       { 25, 6 },
       { 33, 7 },
       { 41, 5 },
       { 48, 6 },
       { 55, 6 }
};

static const int _short_format[NUM_SHELL_FIELD][2] = {
       {  0, 5 },
       {  5, 8 },
       { 14, 1 },
       { 15, 8 },
       { 24, 6 },
       { 32, 7 },
       { 40, 5 },
       { 47, 5 },
       { 53, 5 }
};

static const int _string_format[NUM_STR_FORMAT][2] = {
       { 4, 2 },
       { 7, 4 },
       { 8, 4 }
};

bool Maxit::_pdb_to_ndb_process_remark(const int& remark_no, const std::vector<std::string>&
                                       remarkRecord, const int& remark_num, REMARKS *remark,
                                       const int& repeat)
{
       std::vector<int> index;
       int count = _pdb_to_ndb_get_match_index(remark_no, remarkRecord, remark_num,
                                               remark, repeat, index);

       if (count && remark == Remark_3_Phenix_BIN) {
            _pdb_to_ndb_process_Phenix_Shell_BIN(remarkRecord, index); 
            return true;
       }

       if (count == remark_num && remark == Remark_3_REFMAC_NCS_LOCAL) {
            std::vector<std::string> wordlist;
            for (unsigned int i = 0; i < remarkRecord.size(); ++i) {
                 if (index[i]) continue;
                 get_wordarray(wordlist, remarkRecord[i], " ");
                 if (wordlist.size() != NUM_LOCAL_FIELD) continue;

                 _addNewRecord("NCSLOC");
                 for (int j = 0; j < NUM_LOCAL_FIELD; ++j) {
                      _updateRecordBack("NCSLOC", j + 1, wordlist[j]);
                 }
            }
            return true;
       }

       int count1 = 0;
       for (int i = 0; i < remark_num; ++i) {
            if (remark[i].NumField > 1) count1++;
            else if (strcmp(remark[i].FieldList[0].TokenName, "") ||
                     strcmp(remark[i].FieldList[0].Text, "")) count1++;
       }
       if (remark == Remark_3_BValue_PHENIX && count >= 6) {}
       else if (!count || (((count < (int) remarkRecord.size() && remarkRecord.size() < 8) ||
           (count < 4 * (int) remarkRecord.size() / 5 && remarkRecord.size() >= 8)) &&
           count < 2 * count1 / 3)) return false;

       std::set<std::string> token_Index, order_Index;
       token_Index.clear();
       order_Index.clear();

       std::string block_remark;

       // first key: token ID
       // second key: field_no
       // value: remark context
       std::map<std::string, std::map<int, std::string> > value_pairs;

       for (unsigned int k = 0; k < remarkRecord.size(); ++k) {
            if (!index[k]) continue;

            block_remark.clear();
            for (unsigned int l = k + 1; l < remarkRecord.size(); ++l) {
                 if (index[l]) break;
                 if (remarkRecord[l].empty() && block_remark.empty()) continue;

                 std::string::size_type p = remarkRecord[l].find(":");
                 if (p != std::string::npos) {
                      if ((p > 0) && isdigit(remarkRecord[l][p - 1])) {
                           if (!block_remark.empty()) block_remark += "\n";
                           std::string cs = remarkRecord[l];
                           String::StripLeadingWs(cs);
                           block_remark += cs;
                      } else {
                           std::string cs = remarkRecord[l];
                           String::StripLeadingWs(cs);
                           if (cs[0] == ':') {
                                cs.erase(0, 1);
                                String::StripLeadingWs(cs);
                                if (!block_remark.empty()) block_remark += "\n";
                                block_remark += cs;
                           } else if (cs.substr(0, 3) == "ID:" || cs.find(" ID:") != std::string::npos) {
                                if (!block_remark.empty()) block_remark += "\n";
                                block_remark += cs;
                           }
                      }
                 } else {
                      if (!block_remark.empty()) block_remark += "\n";
                      block_remark += remarkRecord[l];
                 }
            }
         
            bool is_REFMAC_NCS_GROUP = false;
            if (remark == Remark_3_REFMAC_NCS_GROUP && index[k] == 5)
                 is_REFMAC_NCS_GROUP = true;

            _pdb_to_ndb_get_value_pairs(remarkRecord[k], remark[index[k] - 1],
                        is_REFMAC_NCS_GROUP, repeat, block_remark, token_Index,
                        value_pairs);

            if (value_pairs.empty()) continue;

            for (std::map<std::string, std::map<int, std::string> >::const_iterator
                 mpos = value_pairs.begin(); mpos != value_pairs.end(); ++mpos) {
                 const ndb_token_format& ndbformat =
                            NdbToken::getTokenFormat(mpos->first);
                 if (remark[index[k] - 1].Multiple_Treatment < 0) {
                      _addNewRecord(mpos->first);
                      _pdb_to_ndb_insert_back(mpos->first, mpos->second);
                 } else if (repeat > 0) {
                      if ((mpos->first == "TLSGRO" || mpos->first == "ASSEMB" ||
                           mpos->first == "VIRPAR")) {
                           if (mpos->second.find(1) != mpos->second.end())
                                _addNewRecord(mpos->first);
                      } else if (order_Index.find(mpos->first) == order_Index.end()) {
                           order_Index.insert(mpos->first);
                           _addNewRecord(mpos->first);
                      }
                      _pdb_to_ndb_insert_back(mpos->first, mpos->second);
                 } else if (remark[index[k] - 1].Multiple_Treatment > 0)
                      _pdb_to_ndb_insert_with_Multiple_Treatment(remark[index[k]
                                 - 1].Multiple_Treatment, ndbformat.RefineField,
                                 mpos->first, mpos->second);
                 else if (ndbformat.RefineField)
                      _pdb_to_ndb_insert_with_RefineField(ndbformat.RefineField - 1,
                               mpos->first, mpos->second);
                 else _pdb_to_ndb_insert_front(mpos->first, mpos->second);
            }
       }

       return true;
}

int Maxit::_pdb_to_ndb_get_match_index(const int& remark_no, const std::vector<std::string>&
                                       remarkRecord, const int& remark_num, REMARKS *remark,
                                       const int& repeat, std::vector<int>& index)
{
       index.reserve(remarkRecord.size());
       index.clear();
       for (unsigned int i = 0; i < remarkRecord.size(); ++i) index.push_back(0);

       int j = 0;
       int count = 0;
       for (unsigned int i = 0; i < remarkRecord.size(); ++i) {
            if (remark_no == 280 && remarkRecord[i].empty()) continue;

            while (j < remark_num) {
                 if (_pdb_to_ndb_find_a_match(remarkRecord[i], remark[j])) {
                      index[i] = j + 1;
                      count++;
                      break;
                 }
                 j++;
            }
            if (index[i]) {
                 if (repeat) j = 0;
                 else if (j < remark_num && !remark[j].Multiple_Treatment) j++;
            } else j = 0;
       }

       if (remark != Remark_3_REFMAC_NCS_GROUP) return count;

       int range_start = -1;
       int range_end = -1;
       for (unsigned int i = 0; i < remarkRecord.size(); ++i) {
            if (index[i] == 4) range_start = i;
            if (index[i] == 6) range_end = i;
       }
       bool not_found_range = true;
       for (int i = range_start + 1; i < range_end; ++i) {
            if (index[i]) not_found_range = false;
       }
       if (not_found_range) {
            std::vector<std::string> data;
            for (int i = range_start + 1; i < range_end; ++i) {
                 get_wordarray(data, remarkRecord[i], " ");
                 if ((int) data.size() * 2 == remark[4].NumField) {
                      index[i] = 5;
                      count++;
                 }
            }
       }
       return count;
}

bool Maxit::_pdb_to_ndb_find_a_match(const std::string& context, const REMARKS& remark)
{
       bool strip_blanks = true;
       bool label = false;
       std::string cs1 = context;

       for (int i = 0; i < remark.NumField; ++i) {
            if (!strcmp(remark.FieldList[i].TokenName, "") &&
                 strcmp(remark.FieldList[i].Text, "")) {
                 std::string cs = remark.FieldList[i].Text;
                 String::UpperCase(cs);               
                 String::StripLeadingWs(cs);
                 if (cs.empty()) strip_blanks = false;
                 break;
            }
       }
       if (strip_blanks) String::StripLeadingWs(cs1);

       std::vector<std::string> text_tokens;
       text_tokens.clear();
       bool first_match = false;
       for (int i = 0; i < remark.NumField; ++i) {
            if (!strcmp(remark.FieldList[i].TokenName, "") &&
                 strcmp(remark.FieldList[i].Text, "")) {
                 label = true;
                 std::string cs = remark.FieldList[i].Text;
                 String::UpperCase(cs);               
                 String::StripLeadingWs(cs);
                 if (cs == "CRYSTAL" && cs != cs1) {
                      label = false;
                      break;
                 }
                 if (cs1.substr(0, 5) == "SITE_" && cs == "SITE") {
                      label = false;
                      break;
                 }
                 bool not_empty = true;
                 if (cs.empty()) {
                      cs = remark.FieldList[i].Text;
                      not_empty = false;
                 } else text_tokens.push_back(cs);
                 int pos = remark.FieldList[i].FieldStCol;
                 if (strip_blanks)
                      pos = remark.FieldList[i].FieldStCol
                          - remark.FieldList[0].FieldStCol;
                 if (pos < (int) cs1.size()) {
                      if (cs1.substr(pos, cs.size()) != cs) {
                           if (!not_empty) {
                                label = false;
                                break;
                           }
                           String::StripAndCompressWs(cs);
                           std::string cs2 = cs1;
                           String::StripAndCompressWs(cs2);
                           if (cs2.substr(0, cs.size()) != cs) {
                                label = false;
                                break;
                           }
                      } else if (text_tokens.size() == 1) first_match = true;
                 } else {
                      label = false;
                      break;
                 }
            }
       }

       if (!label && first_match) {
            text_tokens.clear();
            for (int i = 0; i < remark.NumField; ++i) {
                 if (!strcmp(remark.FieldList[i].TokenName, "") &&
                      strcmp(remark.FieldList[i].Text, "")) {
                      std::string cs = remark.FieldList[i].Text;
                      String::UpperCase(cs);               
                      String::StripLeadingWs(cs);
                      if (!cs.empty()) text_tokens.push_back(cs);
                 }
            }

            if (text_tokens.size() > 1) {
                 std::string::size_type pos = cs1.find(text_tokens[0]);
                 if (pos != std::string::npos) {
                      label = true;
                      for (unsigned int i = 1; i < text_tokens.size(); ++i) {
                           pos = cs1.find(text_tokens[i], pos + text_tokens[i - 1].size());
                           if (pos == std::string::npos) {
                                label = false;
                                break;
                           }
                      }
                 }
            }
       }

       return label;
}

void Maxit::_pdb_to_ndb_process_Phenix_Shell_BIN(const std::vector<std::string>& remarkRecord,
                                                 const std::vector<int>& index)
{
       std::vector<std::string> wordlist;
       for (unsigned int i = 0; i < remarkRecord.size(); ++i) {
            if (index[i]) continue;
            get_wordarray(wordlist, remarkRecord[i], " ");
            if (wordlist.size() == NUM_SHELL_FIELD) {}
            else if (remarkRecord[i].size() >= 61) {
                 wordlist.clear();
                 for (int j = 0; j < NUM_SHELL_FIELD; ++j) {
                      std::string cs = remarkRecord[i].substr(_long_format[j][0],
                                                              _long_format[j][1]);
                      String::StripAndCompressWs(cs);
                      if (!cs.empty()) wordlist.push_back(cs);
                 }
            } else if (remarkRecord[i].size() >= 58) {
                 wordlist.clear();
                 for (int j = 0; j < NUM_SHELL_FIELD; ++j) {
                      std::string cs = remarkRecord[i].substr(_short_format[j][0],
                                                              _short_format[j][1]);
                      String::StripAndCompressWs(cs);
                      if (!cs.empty()) wordlist.push_back(cs);
                 }
            }
            if (wordlist.size() == NUM_SHELL_FIELD) {
                 for (int j = 0; j < NUM_STR_FORMAT; ++j) {
                      double f = atof(wordlist[_string_format[j][0]].c_str());
                      if (j == 0) {
                           if (f > 1) continue;
                           wordlist[_string_format[j][0]] = FloatToString(f * 100.0, 0,
                                                            _string_format[j][1]);
                      } else if (f > 1) {
                           wordlist[_string_format[j][0]] = FloatToString(f / 100.0, 0,
                                                            _string_format[j][1]);
                      }
                 }

                 _addNewRecord("SHELL");
                 _updateRecordBack("SHELL", 1, wordlist[3]);
                 _updateRecordBack("SHELL", 2, wordlist[1]);
                 _updateRecordBack("SHELL", 3, wordlist[4]);
                 _updateRecordBack("SHELL", 4, wordlist[5]);
                 _updateRecordBack("SHELL", 5, wordlist[7]);
                 _updateRecordBack("SHELL", 6, wordlist[8]);
                 _updateRecordBack("SHELL", 8, wordlist[6]);
            }
       }
       _orderRecordWithFloatKey("SHELL", 1);
}

void Maxit::_pdb_to_ndb_get_value_pairs(const std::string& remarkRecord, const REMARKS&
                        remark, const bool& is_REFMAC_NCS_GROUP, const int& repeat,
                        const std::string& block_remark, std::set<std::string>& token_Index,
                        std::map<std::string, std::map<int, std::string> >& value_pairs)
{
       value_pairs.clear();

       std::vector<std::string> wordlist, data;
       wordlist.clear();

       if (is_REFMAC_NCS_GROUP) {
            get_wordarray(data, remarkRecord, " ");
            if ((int) data.size() * 2 == remark.NumField) {
                 for (unsigned int t = 0; t < data.size(); ++t) {
                      wordlist.push_back(data[t]);
                      wordlist.push_back(data[t]);
                 }
            }
       } else {
            int count = 0;
            bool found = true;
            std::vector<std::string> text_tokens;
            text_tokens.clear();
            for (int i = 0; i < remark.NumField; ++i) {
                 if (remark.FieldList[i].FieldId) count++;
                 if (!strcmp(remark.FieldList[i].TokenName, "") &&
                      strcmp(remark.FieldList[i].Text, "")) {
                      std::string cs = remark.FieldList[i].Text;
                      String::UpperCase(cs);               
                      String::StripLeadingWs(cs);
                      if (cs.empty()) continue;
                      text_tokens.push_back(cs);
                      if (remarkRecord.substr(remark.FieldList[i].FieldStCol - 1, cs.size()) != cs) found = false;
                 }
            }

            if (text_tokens.size() > 1 && !found) {
                 std::string cs = remarkRecord;
                 std::string::size_type pos = cs.find(text_tokens[0]);
                 if (pos != std::string::npos) {
                      found = true;
                      std::string blank_space = "";
                      for (unsigned int j = 0; j < text_tokens[0].size(); ++j) blank_space += " ";
                      cs.replace(pos, text_tokens[0].size(), blank_space);
                      for (unsigned int i = 1; i < text_tokens.size(); ++i) {
                           pos = cs.find(text_tokens[i], pos + text_tokens[i - 1].size());
                           if (pos == std::string::npos) {
                                found = false;
                                break;
                           }
                           blank_space = "";
                           for (unsigned int j = 0; j < text_tokens[i].size(); ++j) blank_space += " ";
                           cs.replace(pos, text_tokens[i].size(), blank_space);
                      }
                      if (found) {
                           get_wordarray(data, cs, " ");
                           if ((int) data.size() == count) {
                                count = 0;
                                for (int i = 0; i < remark.NumField; ++i) {
                                     if (!remark.FieldList[i].FieldId) continue;
                                     std::string token = remark.FieldList[i].TokenName;
                                     if (token.empty()) continue;
                                     int field_no = remark.FieldList[i].FieldId - 1;
                                     _pdb_to_ndb_insert_value_pairs(value_pairs, token, field_no, data[count]);
                                     count++;
                                }
                                return;
                           }
                      }
                 }
            }
       }

       std::string last_token = "";
       int last_fieldno = 0;
       int last_type = 0;

       std::string value;
       for (int i = 0; i < remark.NumField; ++i) {
            if (!remark.FieldList[i].FieldId) continue;

            std::string token = remark.FieldList[i].TokenName;
            if (token.empty()) continue;

            int field_no = remark.FieldList[i].FieldId - 1;
            last_token = token;
            last_fieldno = field_no;
            last_type = remark.FieldList[i].FieldType;

            if (!strcmp(remark.FieldList[i].Text, "RESIDUE RANGE :") &&
                remarkRecord.find("RESIDUE RANGE :") != std::string::npos) {
                 get_wordarray(data, remarkRecord, " ");
                 for (unsigned int t = 2; t < data.size(); ++t) {
                      wordlist.push_back(data[t]);
                 }
            }
            std::string cs = "";
            if (!token.empty() && remark.FieldList[i].FieldId) {
                 bool moving_flag = true;
                 if (token == "CCPNCS" || token == "TLSRNG" ||
                    (token == "RFACTR" && field_no == 6) ||
                    (token == "TLSGRO" && (field_no == 3 || field_no == 4 ||
                     field_no == 5))) moving_flag = false;
                 if ((int) remarkRecord.size() >
                       remark.FieldList[i].FieldStCol -
                       remark.FieldList[0].FieldStCol)
                      _pdb_to_ndb_get_remark_answer(cs, remarkRecord, remark,
                                                    i, moving_flag);
            }

            if ((int) wordlist.size() == remark.NumField) cs = wordlist[i];

            _pdb_to_ndb_process_record_value(cs, false);
            if (cs.empty()) continue;

            if (value_pairs.find(token) == value_pairs.end()) {
                 if (remark.Multiple_Treatment < 0) {
                      if (token == "CCPNCS") {
                           _getRecordBack("NCSGRO", 4, value);
                           if (!value.empty())
                                _pdb_to_ndb_insert_value_pairs(value_pairs, token, 1, value);
                      } else if (token == "TLSRNG") {
                           _getRecordBack("TLSGRO", 1, value);
                           if (!value.empty())
                                _pdb_to_ndb_insert_value_pairs(value_pairs, token, 1, value);
                      }
                 } else if (repeat > 0) {
                      if (token_Index.find(token) == token_Index.end() &&
                          token != "TLSGRO" && token != "ASSEMB" && token != "VIRPAR") {
                           token_Index.insert(token);
                           if (token == "TTENSR" || token == "LTENSR" ||
                               token == "STENSR") {
                                _getRecordBack("TLSGRO", 1, value);
                                if (!value.empty())
                                     _pdb_to_ndb_insert_value_pairs(value_pairs, token, 1, value);
                           } else if (token == "PHENCS") {
                                _getRecordBack("NCSGRO", 4, value);
                                if (!value.empty())
                                     _pdb_to_ndb_insert_value_pairs(value_pairs, token, 1, value);
                           }
                      }
                 }
            }

            _pdb_to_ndb_insert_value_pairs(value_pairs, token, field_no, cs);
       }

       if (block_remark.empty() || last_token.empty() || last_fieldno <= 0) return; 

       std::map<std::string, std::map<int, std::string> >::iterator
           mpos = value_pairs.find(last_token);
       if (mpos == value_pairs.end())
            _pdb_to_ndb_insert_value_pairs(value_pairs, last_token, last_fieldno, block_remark);
       else {
            std::map<int, std::string>::iterator
                mmpos = mpos->second.find(last_fieldno);
            if (mmpos != mpos->second.end()) {
                 if (last_type == 3) mmpos->second += "\n" + block_remark;
            } else mpos->second.insert(std::make_pair(last_fieldno, block_remark));
       }
}

void Maxit::_pdb_to_ndb_get_remark_answer(std::string& answer, const std::string&
                                          remarkCard, const REMARKS& remark, const
                                          int& field_no, const bool& moving_flag)
{
       answer.clear();
       if (remarkCard.empty()) return;

       bool strip_blanks = true;
       std::string cs = remark.FieldList[0].Text;
       String::UpperCase(cs);               
       String::StripLeadingWs(cs);
       if (cs.empty()) strip_blanks = false;

       cs = remarkCard;
       if (strip_blanks) String::StripLeadingWs(cs);

       if (field_no == 0) {
            if (cs != "NULL" && cs != "N/A") 
                 answer = cs.substr(0, remark.FieldList[field_no].FieldWidth);
            else answer = "NULL";
       } else {
            int tmp = remark.FieldList[field_no].FieldStCol;
            int end = remark.FieldList[field_no].FieldStCol
                    + remark.FieldList[field_no].FieldWidth;
            if (strip_blanks) {
                 tmp = remark.FieldList[field_no].FieldStCol
                     - remark.FieldList[0].FieldStCol;
                 end = remark.FieldList[field_no].FieldStCol 
                     - remark.FieldList[0].FieldStCol
                     + remark.FieldList[field_no].FieldWidth;
            }
            if (field_no == remark.NumField - 1) {
                 const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
                 end = ndbformat.FieldList[2].FieldWidth;
            }
            if (end > (int) cs.size()) end = cs.size();
            int i = 0;
            for (i = tmp; i < end; ++i)
                 if (cs[i] != ' ' && cs[i] != ':' &&
                     cs[i] != ';' && cs[i] != '(' &&
                     cs[i] != ')') break;
            if (moving_flag) {
                 while (i && cs[i-1] != ' ' && cs[i-1] != ':' &&
                             cs[i-1] != ';') i--;
            }
            int width = end - i;
            if (i >= 0 && width > 0) answer = cs.substr(i, width);
       }
}

void Maxit::_pdb_to_ndb_insert_value_pairs(std::map<std::string, std::map<int, std::string> >&
                                           value_pairs, const std::string& token, const
                                           int& field_no, const std::string& value)
{
       std::map<std::string, std::map<int, std::string> >::iterator
           mpos = value_pairs.find(token);
       if (mpos != value_pairs.end()) {
            std::map<int, std::string>::iterator
                mmpos = mpos->second.find(field_no);
            if (mmpos != mpos->second.end())
                 mmpos->second = value;
            else mpos->second.insert(std::make_pair(field_no, value));
       } else {
            std::map<int, std::string> tmp_map;
            tmp_map.clear();
            tmp_map.insert(std::make_pair(field_no, value));
            value_pairs.insert(std::make_pair(token, tmp_map));
       }
}

void Maxit::_pdb_to_ndb_insert_front(const std::string& token, const std::map<int,
                                      std::string>& value_pairs)
{
       for (std::map<int, std::string>::const_iterator
            mpos = value_pairs.begin(); mpos != value_pairs.end(); ++mpos) {
            _updateRecordFront(token, mpos->first, mpos->second);
       }
}

void Maxit::_pdb_to_ndb_insert_back(const std::string& token, const std::map<int,
                                    std::string>& value_pairs)
{
       for (std::map<int, std::string>::const_iterator
            mpos = value_pairs.begin(); mpos != value_pairs.end(); ++mpos) {
            _updateRecordBack(token, mpos->first, mpos->second);
       }
}

void Maxit::_pdb_to_ndb_insert_with_Multiple_Treatment(const int& Multiple_Treatment,
                                  const int& RefineFieldNo, const std::string& token,
                                  const std::map<int, std::string>& value_pairs)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) {
            _addNewRecord(token);
            _pdb_to_ndb_insert_back(token, value_pairs);
            return;
       }

       int field_no = -1;
       std::string value = "";
       // bool XFILES_flag = false;
       for (std::map<int, std::string>::const_iterator
            // if (mpos->first == 1 && mpos->second == "1" && token == "XFILES")
            //     XFILES_flag = true;
            mpos = value_pairs.begin(); mpos != value_pairs.end(); ++mpos) {
            if (Multiple_Treatment == (mpos->first + 1)) {
                 field_no = mpos->first;
                 value = mpos->second;
            }
       }
/*
       if (token == "XFILES") {
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[RefineFieldNo - 1].empty()) {
                      for (std::map<int, std::string>::const_iterator
                           mpos = value_pairs.begin(); mpos != value_pairs.end(); ++mpos) {
                           (*lpos)[mpos->first] = mpos->second;
                      }
                      return;
                 }
            }
       } else */ if (field_no >= 0) {
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if (RefineFieldNo && !(*lpos)[RefineFieldNo - 1].empty()) continue;
                 if ((*lpos)[field_no] == value || (token == "XFILES" && (*lpos)[field_no].empty())) {
                      for (std::map<int, std::string>::const_iterator
                           mpos = value_pairs.begin(); mpos != value_pairs.end(); ++mpos) {
                           (*lpos)[mpos->first] = mpos->second;
                      }
                      return;
                 }
            }
       }

       _addNewRecord(token);
       _pdb_to_ndb_insert_back(token, value_pairs);
}

void Maxit::_pdb_to_ndb_insert_with_RefineField(const int& RefineFieldNo, const std::string&
                                       token, const std::map<int, std::string>& value_pairs)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) {
            _addNewRecord(token);
            _pdb_to_ndb_insert_back(token, value_pairs);
            return;
       }

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[RefineFieldNo].empty()) {
                 for (std::map<int, std::string>::const_iterator
                      mpos = value_pairs.begin(); mpos != value_pairs.end(); ++mpos) {
                      (*lpos)[mpos->first] = mpos->second;
                 }
                 return;
            }
       }

       _addNewRecord(token);
       _pdb_to_ndb_insert_back(token, value_pairs);
}
