/*
FILE:     Ndb2Pdb_Remark_Util.C
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

#include "Maxit.h"
#include "NdbToken.h"
#include "Remark_define.h"
#include "Remark_extern.h"
#include "utillib.h"

void Maxit::_ndb_to_pdb_get_general_remark(const int& num_remarks, REMARKS *Remarks, const int& repeat, const int& index, const int& empty_flag,
                              const std::string& sctype)
{
       _ndb_to_pdb_update_template_FieldWidth(num_remarks, Remarks);

       if (Remarks == Remark_3_REFMAC_TWIN) {
            _ndb_to_pdb_get_REFMAC_TWIN_remark(num_remarks, Remarks, sctype);
            return;
       }

       std:: string card_id = "";
       for (int i = 0; i < num_remarks; ++i) {
            card_id.clear();
            for (int j = 0; j < Remarks[i].NumField; ++j) {
                 if (strcmp(Remarks[i].FieldList[j].TokenName, "") && card_id.empty()) card_id = Remarks[i].FieldList[j].TokenName;
            }
       }

       if (repeat && !card_id.empty())
            _ndb_to_pdb_get_repeat_remark(card_id, num_remarks, Remarks, empty_flag);
       else _ndb_to_pdb_get_non_repeat_remark(num_remarks, Remarks, index, sctype);

       if (Remarks == Remark_3_REFMAC_NCS_LOCAL) _ndb_to_pdb_get_REFMAC_NCS_LOCAL(Remarks->Remark_No);

       if (Remarks == Remark_3_Phenix_BIN) _ndb_to_pdb_get_Phenix_BIN_remark(Remarks->Remark_No);

       if (Remarks == Remark_3_PHENIX_NCS) _ndb_to_pdb_get_PHENIX_NCS_remark();
}

void Maxit::_ndb_to_pdb_update_template_FieldWidth(const int& num_remarks, REMARKS *Remarks)
{
       for (int i = 0; i < num_remarks; ++i) {
            for (int j = 0; j < Remarks[i].NumField; ++j) {
                 if (!strcmp(Remarks[i].FieldList[j].TokenName, "") && Remarks[i].FieldList[j].FieldId == 0 && Remarks[i].FieldList[j].FieldWidth == 0) {
                      Remarks[i].FieldList[j].FieldWidth = strlen(Remarks[i].FieldList[j].Text);
                 }
            }
       }
}

void Maxit::_ndb_to_pdb_get_REFMAC_TWIN_remark(const int& num_remarks, REMARKS *Remarks, const std::string& sctype)
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("TWIN");
       if (ppos == _pdb_records.end()) return;

       _addNewRemark(Remarks->Remark_No, "");

       _ndb_to_pdb_update_template_FieldWidth(Num_Remark_3_REFMAC_TWIN_OP, Remark_3_REFMAC_TWIN_OP);

       for (int i = 0; i < num_remarks; ++i) {
            _ndb_to_pdb_get_general_remark(Remarks[i], sctype);
       }

       for (std::list<std::vector<std::string> >::const_iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            for (int i = 0; i < Num_Remark_3_REFMAC_TWIN_OP; ++i) {
                 _ndb_to_pdb_get_general_remark(Remark_3_REFMAC_TWIN_OP[i], "", "", *lppos);
            }
       }
}

void Maxit::_ndb_to_pdb_get_repeat_remark(const std::string& card_id, const int& num_remarks, REMARKS *Remarks, const int& empty_flag)
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find(card_id);
       if (ppos == _pdb_records.end()) return;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat(card_id);

       bool found = false;
       for (std::list<std::vector<std::string> >::const_iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            if (!_current_method.empty() && ndbformat.RefineField &&
                 _current_method != (*lppos)[ndbformat.RefineField - 1]) continue;
            found = true;
            break;
       }

       if (Remarks != Remark_3_PHENIX_TLS_GROUP && Remarks != &Remark_3_XPlor4 && found)
            _addNewRemark(Remarks->Remark_No, "");

       std::vector<int> single;
       single.clear();
       for (int i = 0; i < num_remarks; ++i) single.push_back(1);

       for (std::list<std::vector<std::string> >::const_iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            if (!_current_method.empty() && ndbformat.RefineField && _current_method != (*lppos)[ndbformat.RefineField - 1]) continue;

            if (empty_flag && lppos != ppos->second.begin()) _addNewRemark(Remarks->Remark_No, "");

            for (int i = 0; i < num_remarks; ++i) {
                 if (Remarks[i].NumField == 1) {
                      if (single[i] || Remarks == Remark_265_2) {
                           _ndb_to_pdb_get_general_remark(Remarks[i], "", "", *lppos);
                           single[i] = 0;
                      } else if (Remarks == Remark_800) {
                           _addNewRemark(Remarks->Remark_No, "");
                      }
                 } else _ndb_to_pdb_get_general_remark(Remarks[i], "", "", *lppos);
            }
       }
}

void Maxit::_ndb_to_pdb_get_non_repeat_remark(const int& num_remarks, REMARKS *Remarks, const int& index, const std::string& sctype)
{
       if (Remarks != Remark_3_PHENIX_TLS_GROUP) _addNewRemark(Remarks->Remark_No, "");

       for (int i = 0; i < num_remarks; ++i) {
            std::string card_id = "";
            int field_no = 1;
            for (int j = 0; j < Remarks[i].NumField; ++j) {
                 if (strcmp(Remarks[i].FieldList[j].TokenName, "")) {
                      card_id = Remarks[i].FieldList[j].TokenName;
                      if (card_id == "NCSGRO")
                           field_no = 4;
                      else field_no = 1;
                 }
                 if (!card_id.empty()) break;
            }

            if ((Remarks[i].Multiple_Treatment || index) && !card_id.empty()) {
                 std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find(card_id);
                 if (ppos == _pdb_records.end()) continue;

                 const ndb_token_format& ndbformat = NdbToken::getTokenFormat(card_id);

                 for (std::list<std::vector<std::string> >::const_iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                      if (!_current_method.empty() && ndbformat.RefineField && _current_method != (*lppos)[ndbformat.RefineField - 1]) continue;

                      if (!index || (index && atoi((*lppos)[field_no].c_str()) == index)) _ndb_to_pdb_get_general_remark(Remarks[i], "", "", *lppos);
                 }
            } else _ndb_to_pdb_get_general_remark(Remarks[i], sctype);
       }
}

void Maxit::_ndb_to_pdb_get_REFMAC_NCS_LOCAL(const int& Remark_No)
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("NCSLOC");
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::const_iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            std::string buffer = FormattedString((*lppos)[1], 4, false, true) + " " + FormattedString((*lppos)[2], 5, false, true) + " "
                               + FormattedString((*lppos)[3], 5, false, true) + " " + FormattedString((*lppos)[4], 6, false, true) + " "
                               + FormattedString((*lppos)[5], 7, false, true) + " " + FormattedString((*lppos)[6], 5, false, true) + " "
                               + FormattedString((*lppos)[7], 6, false, true) + " " + FormattedString((*lppos)[8], 7, false, true) + " "
                               + FormattedString((*lppos)[9], 5, false, true) + " " + FormattedString((*lppos)[10], 5, false, true);
            _addNewRemark(Remark_No, buffer);
       }
}

void Maxit::_ndb_to_pdb_get_Phenix_BIN_remark(const int& Remark_No)
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("SHELL");
       if (ppos == _pdb_records.end()) return;

       _orderRecordWithFloatKey("SHELL", 1, true);

       ppos = _pdb_records.find("SHELL");
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("SHELL");

       int count = 0;
       for (std::list<std::vector<std::string> >::const_iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            if (!_current_method.empty() && ndbformat.RefineField && _current_method != (*lppos)[ndbformat.RefineField - 1]) continue;

            count++;
            std::string buffer = FloatToString(count, 5, 0, false, true) + FloatToString(atof((*lppos)[2].c_str()), 8, 4, false, true) + " -"
                               + FloatToString(atof((*lppos)[1].c_str()), 8, 4, false, true) + "    "
                               + FloatToString(atof((*lppos)[3].c_str()) / 100.0, 4, 2, false, true) + "   "
                               + FloatToString(atoi((*lppos)[4].c_str()), 6, 0, false, true) + " "
                               + FloatToString(atoi((*lppos)[8].c_str()), 5, 0, false, true) + "  "
                               + FloatToString(atof((*lppos)[5].c_str()), 6, 4, false, true) + " "
                               + FloatToString(atof((*lppos)[6].c_str()), 6, 4, false, true);
            _addNewRemark(Remark_No, buffer);
       }
}

void Maxit::_ndb_to_pdb_get_PHENIX_NCS_remark()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("NCSGRO");
       if (ppos == _pdb_records.end()) return;

       _ndb_to_pdb_update_template_FieldWidth(1, &Remark_3_PHENIX_NCS_GROUP);
       _ndb_to_pdb_update_template_FieldWidth(Num_Remark_3_PHENIX_NCS_GROUP_OP, Remark_3_PHENIX_NCS_GROUP_OP);

       for (std::list<std::vector<std::string> >::const_iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            _ndb_to_pdb_get_general_remark(Remark_3_PHENIX_NCS_GROUP, "", "", *lppos);
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator qpos = _pdb_records.find("PHENCS");
            if (qpos == _pdb_records.end()) continue;

            for (std::list<std::vector<std::string> >::const_iterator lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                 if ((*lqpos)[1] != (*lppos)[4]) continue;

                 for (int i = 0; i < Num_Remark_3_PHENIX_NCS_GROUP_OP; ++i) {
                      _ndb_to_pdb_get_general_remark(Remark_3_PHENIX_NCS_GROUP_OP[i], "", "", *lqpos);
                 }
            }
       }
}

void Maxit::_ndb_to_pdb_get_general_remark(const REMARKS& Remark, const std::string& sctype)
{
       std::string diffrn_id = "";
       if (!sctype.empty()) {
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("SCTYPE");
            if (ppos != _pdb_records.end()) {
                 for (std::list<std::vector<std::string> >::const_iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                      if ((*lppos)[2] != sctype) continue;
                      diffrn_id = (*lppos)[1];
                      break;
                 }
            }
       }

       bool is_skip_remark = false;
       bool is_null_value = true;
       for (int i = 0; i < Remark.NumField; i++) {
            if (!strcmp(Remark.FieldList[i].Text, "B VALUE TYPE :")) {
                 is_skip_remark = true;
            }
            if (strcmp(Remark.FieldList[i].TokenName, "") &&
                Remark.FieldList[i].FieldId) {
                 std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find(Remark.FieldList[i].TokenName);
                 if (ppos == _pdb_records.end()) continue;

                 const std::vector<std::string>& Field = ppos->second.front();
                 if (!Field[Remark.FieldList[i].FieldId-1].empty()) is_null_value = false;
            }
       }
       if (is_skip_remark && is_null_value) return;

       std::vector<std::string> pCardField;
       pCardField.clear();
       _ndb_to_pdb_get_general_remark(Remark, sctype, diffrn_id, pCardField);
}

void Maxit::_ndb_to_pdb_get_general_remark(const REMARKS& Remark, const std::string& sctype, const std::string& diffrn_id,
                                           const std::vector<std::string>& pCardField)
{
       std::vector<std::string> remark_array;
       remark_array.clear();

       std::string remark_string = "";
       std::string value = "";

       for (int i = 0; i < Remark.NumField; ++i) {
            if (!strcmp(Remark.FieldList[i].TokenName, "") && !strcmp(Remark.FieldList[i].Text, "")) continue;

            remark_string += get_space_between_remark_field(Remark, i);

            bool last_field = false;
            if (i == Remark.NumField - 1) last_field = true;

            bool left_adjust = false;
            if (Remark.FieldList[i].FieldJustification == PDB_REMARK_LEFT_JUSTIFIED) left_adjust = true;

            if (strcmp(Remark.FieldList[i].TokenName, "") && Remark.FieldList[i].FieldId) {
                 if (!pCardField.empty())
                      value = pCardField[Remark.FieldList[i].FieldId - 1];
                 else value = _ndb_to_pdb_get_field_value(Remark.FieldList[i], Remark.Remark_No, sctype, diffrn_id, last_field);

                 bool clean_flag = true;
                 if ((!strcmp(Remark.FieldList[i].TokenName, "CRMET1") || !strcmp(Remark.FieldList[i].TokenName, "NMREMK") ||
                      !strcmp(Remark.FieldList[i].TokenName, "REFREM")) && Remark.FieldList[i].FieldId == 3 &&
                      (int) value.size() > Remark.FieldList[i].FieldWidth) clean_flag = false;
                 if (clean_flag) String::StripAndCompressWs(value);
               
                 bool force_reformat = false;
                 if ((!strcmp(Remark.FieldList[i].TokenName, "TLSGRO") || !strcmp(Remark.FieldList[i].TokenName, "TTENSR") ||
                      !strcmp(Remark.FieldList[i].TokenName, "LTENSR") || !strcmp(Remark.FieldList[i].TokenName, "STENSR")) && 
                     (int) value.size() > Remark.FieldList[i].FieldWidth && Remark.FieldList[i].FieldType == 2)
                      force_reformat = true;

                 if (value.empty())
                      remark_string += FormattedFieldValue("NULL", 3, Remark.FieldList[i].FieldWidth, 0, left_adjust);
                 else if (force_reformat)
                      remark_string += FormattedFieldValue(value, Remark.FieldList[i].FieldType, Remark.FieldList[i].FieldWidth,
                                                           Remark.FieldList[i].FieldPrec, left_adjust);
                 else if (!last_field || (int) value.size() <= Remark.FieldList[i].FieldWidth)
                      remark_string += FormattedFieldValue(value, 3, Remark.FieldList[i].FieldWidth, 0, left_adjust, true);
                 else _ndb_to_pdb_format_field_value(Remark.Remark_No, Remark.FieldList[0], Remark.FieldList[i], value, remark_string, remark_array);
            } else remark_string += FormattedFieldValue(Remark.FieldList[i].Text, Remark.FieldList[i].FieldType,
                                         Remark.FieldList[i].FieldWidth, Remark.FieldList[i].FieldPrec, left_adjust);
       }

       if (remark_array.empty()) {
            if (Remark.Remark_No == 100) {
                 String::StripTrailingWs(remark_string);
                 remark_string += ".";
            }
            _addNewRemark(Remark.Remark_No, remark_string);
       } else {
            for (unsigned int i = 0; i < remark_array.size(); ++i) {
                 if (Remark.Remark_No == 100 && i == (remark_array.size() - 1)) {
                      String::StripTrailingWs(remark_array[i]);
                      remark_array[i] += ".";
                 }
                 _addNewRemark(Remark.Remark_No, remark_array[i]);
            }
       }
}

std::string Maxit::_ndb_to_pdb_get_field_value(const REMARK_FIELD& FieldList, const int& Remark_No, const std::string& sctype,
                                               const std::string& diffrn_id, const bool& last_field)
{
       std::string value = "";

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find(FieldList.TokenName);
       if (ppos == _pdb_records.end()) return value;

       std::string card_id = FieldList.TokenName;
       int field_no = FieldList.FieldId - 1;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat(card_id);

       int serial_no_field = 0;
       if (card_id == "DIFPTL" || card_id == "DTMEAS" || card_id == "DTMETH" || card_id == "DTTEMP" || card_id == "DTWAVE" || card_id == "PERCOM" ||
           card_id == "PHVAL"  || card_id == "RADIAT" || card_id == "RESTOT" || card_id == "REFLEC" || card_id == "SHELC"  || card_id == "WAVLEN")
            serial_no_field = ndbformat.SeqField - 1;
       else if (card_id == "DTWAVE")
            serial_no_field = 2;

       int FieldWidth = get_remark_field_width(FieldList, last_field);

       for (std::list<std::vector<std::string> >::const_iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            if ((*lppos)[field_no].empty()) continue;

            std::string cs = (*lppos)[field_no];
            String::StripTrailingWs(cs);
            if (!ndbformat.FieldList[field_no].MixedCase) String::UpperCase(cs);

            if (Remark_No == 245) {
                 if (FieldList.FieldType < 3) _ndb_to_pdb_reformat_numeric_value(cs, FieldWidth, FieldList);
                 String::StripAndCompressWs(cs);
                 if (!value.empty()) value += "; ";
                 value += cs;
            } else {
                 if (!_current_method.empty() && ndbformat.RefineField && (*lppos)[ndbformat.RefineField - 1] != _current_method) continue;
             
                 if (!sctype.empty()) {
                      if (card_id == "EXPDTA" && (*lppos)[2] != sctype) continue;
                      else if (card_id == "CRMETH" && (*lppos)[5] != sctype) continue;
                      else if (!diffrn_id.empty() && serial_no_field && !(*lppos)[serial_no_field].empty() &&
                               !IsSameStringCombination(diffrn_id, (*lppos)[serial_no_field], " ,")) continue;
                 }

                 if (FieldList.FieldType < 3) _ndb_to_pdb_reformat_numeric_value(cs, FieldWidth, FieldList);
                 value = cs;
                 break;
            }
       }

       return value;
}

void Maxit::_ndb_to_pdb_reformat_numeric_value(std::string& value, const int& width, const REMARK_FIELD& FieldList)
{
       if (!String::IsNumber(value)) return;

       if (_remark_decimal_precision == UNFORCE_DECIMAL_PRECISION && (int) value.size() <= width) return;

       int left_adjust = false;
       if (FieldList.FieldJustification == PDB_REMARK_LEFT_JUSTIFIED) left_adjust = true;

       value = FormattedFieldValue(value, FieldList.FieldType, FieldList.FieldWidth, FieldList.FieldPrec, left_adjust, true);
}

void Maxit::_ndb_to_pdb_format_field_value(const int& Remark_No, const REMARK_FIELD& FieldList0, const REMARK_FIELD& FieldList, std::string& value,
                                           std::string& remark_string, std::vector<std::string>& remark_array)
{
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
       int FieldLen = ndbformat.FieldList[2].FieldWidth;

       remark_array.clear();

       std::vector<std::string> data;

       if ((!strcmp(FieldList.TokenName, "CRMET1") || !strcmp(FieldList.TokenName, "NMREMK") || !strcmp(FieldList.TokenName, "REFREM")) &&
             FieldList.FieldId == 3 && value.find("\n") != std::string::npos) {
            get_wordarray_with_space(data, value, "\n");
            int count = 0;
            for (unsigned int i = 0; i < data.size(); ++i) {
                 String::StripTrailingWs(data[i]);
                 if (data[i].empty() && i == 0) continue;
                 count++;
            }
            if (count > 1) {
                 std::vector<std::string> data1;
                 remark_array.push_back(remark_string);
                 bool first = true;
                 for (unsigned int i = 0; i < data.size(); ++i) {
                      if (data[i].empty() && first) continue;
                      if ((int) data[i].size() < FieldLen) {
                           if (first && ((int) data[i].size() + (int) remark_string.size() + 1) < FieldLen)
                                remark_array[0] += data[i];
                           else remark_array.push_back(" " + data[i]);
                      } else {
                           get_max_length_words(data1, data[i], FieldLen - 1);
                           for (unsigned int j = 0; j < data1.size(); ++j) {
                                remark_array.push_back(" " + data1[j]);
                           }
                      }
                      first = false;
                 }
            }
       }
       if (!remark_array.empty()) return;

       String::StripAndCompressWs(value);

       if (((int) remark_string.size() + (int) value.size() + 1) < FieldLen) {
            remark_string += value;
            return;
       }

       int max_length = FieldLen - (remark_string.size() + 2);

       std::string::size_type r = std::string::npos;
       std::string::size_type p = value.find_first_of(" -");
       while (p != std::string::npos &&
             (int) value.substr(0, p).size() <= max_length) {
            r = p;
            p = value.find_first_of(" -", p + 1);
       }
       if (r != std::string::npos) {
            remark_string += value.substr(0, r);
            value = value.substr(r);
            String::StripAndCompressWs(value);
       }

       max_length = FieldLen - 3;
       std::string space = " ";

       int align = check_align(Remark_No, FieldList);
       if (align) {
            int length = strlen(FieldList0.Text) + FieldList0.FieldStCol;
            max_length = FieldLen - length - 1;
            space = PrintSpace(length);
            get_wordarray(data, value, " ");
            for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
                 if ((int) pos->size() > (max_length - 1)) {
                      max_length = FieldLen - 3;
                      space = " ";
                      break;
                 }
            }
       }

       get_max_length_words(data, value, max_length);
       remark_array.push_back(remark_string);
       for (std::vector<std::string>::iterator pos = data.begin(); pos != data.end(); ++pos) {
            String::StripLeadingWs(*pos);
            remark_array.push_back(space + *pos);
       }
}

void Maxit::_ndb_to_pdb_add_remark(const int& Remark_No, const std::vector<std::string>& remark_array)
{
       if (remark_array.empty()) return;

       for (std::vector<std::string>::const_iterator pos = remark_array.begin(); pos != remark_array.end(); ++pos) {
            _addNewRemark(Remark_No, *pos);
       }
}
