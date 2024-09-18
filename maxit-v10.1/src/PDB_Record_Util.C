/*
FILE:     PDB_Record_Util.C
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
#include <ctype.h>

#include "Maxit.h"
#include "NdbToken.h"
#include "ReservedWord.h"
#include "TypeDef.h"
#include "utillib.h"

void Maxit::_addNewRecord(const std::string& token)
{
       std::vector<std::string> FieldInfo;
       FieldInfo.clear();

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat(token);
       for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
       if (NdbToken::IsJrnlToken(token)) {
            FieldInfo[0] = "JRNL";
            FieldInfo[2] = token;
       } else FieldInfo[0] = token;

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) {
            std::list<std::vector<std::string> > rlist;
            rlist.clear();
            rlist.push_back(FieldInfo);
            _pdb_records.insert(std::make_pair(token, rlist));
       } else ppos->second.push_back(FieldInfo);
}

void Maxit::_updateRecordFront(const std::string& token, const int& field_no,
                               const std::string& value, const bool& force_flag)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) {
            if (force_flag) {
                 _addNewRecord(token);
                 ppos = _pdb_records.find(token);
            } else return;
       }

       std::vector<std::string>& Field = ppos->second.front();
       if (force_flag || Field[field_no].empty()) Field[field_no] = value;
}

void Maxit::_updateRecordFront(const std::string& token, const int& field_no,
                               const std::string& value, const std::string& delimiter,
                               const bool& force_flag)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) {
            if (force_flag) {
                 _addNewRecord(token);
                 ppos = _pdb_records.find(token);
            } else return;
       }

       std::vector<std::string>& Field = ppos->second.front();
       if (Field[field_no].empty())
            Field[field_no] = value;
       else Field[field_no] += delimiter + value;
}

void Maxit::_updateRecordBack(const std::string& token, const int& field_no,
                              const std::string& value, const bool& force_flag)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) {
            if (force_flag) {
                 _addNewRecord(token);
                 ppos = _pdb_records.find(token);
            } else return;
       }

       std::vector<std::string>& Field = ppos->second.back();
       if (force_flag || Field[field_no].empty()) Field[field_no] = value;
}

void Maxit::_updateRecordBack(const std::string& token, const int& field_no,
                              const std::string& value, const std::string& delimiter,
                              const bool& force_flag)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) {
            if (force_flag) {
                 _addNewRecord(token);
                 ppos = _pdb_records.find(token);
            } else return;
       }

       std::vector<std::string>& Field = ppos->second.back();
       if (Field[field_no].empty())
            Field[field_no] = value;
       else Field[field_no] += delimiter + value;
}

void Maxit::_getRecordFront(const std::string& token, const int& field_no,
                            std::string& value)
{
       value.clear();

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       const std::vector<std::string>& Field = ppos->second.front();
       if (field_no >= 1 && field_no < (int) Field.size()) value = Field[field_no];
}

void Maxit::_getRecordBack(const std::string& token, const int& field_no,
                           std::string& value)
{
       value.clear();

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       const std::vector<std::string>& Field = ppos->second.back();
       if (field_no >= 1 && field_no < (int) Field.size()) value = Field[field_no];
}

void Maxit::_copyValue(const std::string& source_token, const int& source_field,
                       const std::string& target_token, const int& target_field,
                       const bool& force_flag)
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(source_token);
       if (ppos == _pdb_records.end()) return;

       const std::vector<std::string>& source = ppos->second.front();
       if (source[source_field].empty()) return;

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           qpos = _pdb_records.find(target_token);
       if (qpos == _pdb_records.end()) _addNewRecord(target_token);

       _updateRecordFront(target_token, target_field, source[source_field], force_flag);
}

void Maxit::_copyValue(const std::string& source_token, const std::string& target_token,
                       const std::map<int, int>& field_mapping, const bool& force_flag)
{
       if (field_mapping.empty()) return;

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(source_token);
       if (ppos == _pdb_records.end()) return;

       bool has_value = false;
       const std::vector<std::string>& source = ppos->second.front();
       for (std::map<int, int>::const_iterator
            mpos = field_mapping.begin(); mpos != field_mapping.end(); ++mpos) {
            if (!source[mpos->first].empty()) {
                 has_value = true;
                 break;
            }
       }
       if (!has_value) return;

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           qpos = _pdb_records.find(target_token);
       if (qpos == _pdb_records.end()) _addNewRecord(target_token);

       for (std::map<int, int>::const_iterator
            mpos = field_mapping.begin(); mpos != field_mapping.end(); ++mpos) {
            if (source[mpos->first].empty()) continue;
            _updateRecordFront(target_token, mpos->second, source[mpos->first], force_flag);
       }
}

void Maxit::_updateValue(const std::string& token, const int& key_field,  const std::string&
                         key_value, const std::map<int, std::string>& field_value_pair)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[key_field] == key_value) {
                      for (std::map<int, std::string>::const_iterator mpos =
                           field_value_pair.begin(); mpos != field_value_pair.end(); ++mpos) {
                           (*lpos)[mpos->first] = mpos->second;
                      }
                      return;
                 }
            }
       }

       _addNewRecord(token);
       ppos = _pdb_records.find(token);
       std::vector<std::string>& Field = ppos->second.back();
       Field[key_field] = key_value;
       for (std::map<int, std::string>::const_iterator
            mpos = field_value_pair.begin(); mpos != field_value_pair.end(); ++mpos) {
            Field[mpos->first] = mpos->second;
       }
}

void Maxit::_getValue(const std::string& token, const int& key_field,  const std::string&
                      key_value, const int& fieldno, std::string& value)
{
       value.clear();
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[key_field] == key_value) {
                 value = (*lpos)[fieldno];
                 return;
            }
       }
}

void Maxit::_orderRecordWithIntegerKey(const std::string& token, const int& fieldno)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       std::multimap<int, std::vector<std::string> > order_mapping;
       order_mapping.clear();

       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            order_mapping.insert(std::make_pair(atoi((*lpos)[fieldno].c_str()), *lpos));
       }

       ppos->second.clear();
       for (std::multimap<int, std::vector<std::string> >::const_iterator
            mpos = order_mapping.begin(); mpos != order_mapping.end(); ++mpos) {
            ppos->second.push_back(mpos->second);
       }
}

void Maxit::_orderRecordWithFloatKey(const std::string& token, const int& fieldno,
                                     const bool& reverse_order)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       std::multimap<double, std::vector<std::string> > order_mapping;
       order_mapping.clear();

       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            order_mapping.insert(std::make_pair(atof((*lpos)[fieldno].c_str()), *lpos));
       }

       ppos->second.clear();
       if (reverse_order) {
            for (std::multimap<double, std::vector<std::string> >::const_reverse_iterator
                 mpos = order_mapping.rbegin(); mpos != order_mapping.rend(); ++mpos) {
                 ppos->second.push_back(mpos->second);
            }
       } else {
            for (std::multimap<double, std::vector<std::string> >::const_iterator
                 mpos = order_mapping.begin(); mpos != order_mapping.end(); ++mpos) {
                 ppos->second.push_back(mpos->second);
            }
       }
}

void Maxit::_orderRecordWithOrderKey(const std::string& token, const int& fieldno,
                                     const std::vector<std::string>& order)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       std::map<std::string, std::vector<std::vector<std::string> > > order_mapping;
       order_mapping.clear();

       std::vector<std::vector<std::string> > tmp_vector;
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            std::map<std::string, std::vector<std::vector<std::string> > >::iterator
                mpos = order_mapping.find((*lpos)[fieldno]);
            if (mpos != order_mapping.end())
                 mpos->second.push_back(*lpos);
            else {
                 tmp_vector.clear();
                 tmp_vector.push_back(*lpos);
                 order_mapping.insert(std::make_pair((*lpos)[fieldno], tmp_vector));
            }
       }

       ppos->second.clear();
       for (std::vector<std::string>::const_iterator
            vpos1 = order.begin(); vpos1 != order.end(); ++vpos1) {
            std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator
                mpos = order_mapping.find(*vpos1);
            if (mpos == order_mapping.end()) continue;
            for (std::vector<std::vector<std::string> >::const_iterator
                 vpos2 = mpos->second.begin(); vpos2 != mpos->second.end(); ++vpos2) {
                 ppos->second.push_back(*vpos2);
            }
            order_mapping.erase(*vpos1);
       }

       for (std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator
            mpos = order_mapping.begin(); mpos != order_mapping.end(); ++mpos) {
            for (std::vector<std::vector<std::string> >::const_iterator
                 vpos2 = mpos->second.begin(); vpos2 != mpos->second.end(); ++vpos2) {
                 ppos->second.push_back(*vpos2);
            }
       }
}

void Maxit::_cleanString(const std::string& token, const int& fieldno)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            String::StripAndCompressWs((*lpos)[fieldno]);
       }
}

void Maxit::_removeWhiteSpace(const std::string& token, const int& fieldno)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       std::string cs;
       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            String::RemoveWhiteSpace((*lpos)[fieldno], cs);
            (*lpos)[fieldno] = cs;
       }
}

void Maxit::_removeWhiteSpaceBetweenName(const std::string& token, const int& fieldno)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            delete_space_between_names((*lpos)[fieldno]);
       }
}

void Maxit::_mixCase(const std::string& token, const int& fieldno)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            ReservedWord::SetStringMixedCase((*lpos)[fieldno]);
       }
}

void Maxit::_updateRValueForCNS_XPLOR()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("RFACTR");
       if (ppos == _pdb_records.end()) return;

       std::vector<std::string>& Field = ppos->second.front();
       if (!Field[7].empty() && Field[10].empty())
            Field[10] = Field[7];
       else if (Field[7].empty() && !Field[10].empty())
            Field[7] = Field[10];
}

void Maxit::_updateRefineMethod(const std::string& method)
{
       if (method.empty()) return;

       const std::vector<std::string>& ndbTokens = NdbToken::getNdbTokens();
       for (std::vector<std::string>::const_iterator
            tpos = ndbTokens.begin(); tpos != ndbTokens.end(); ++tpos) {
            const ndb_token_format& ndbformat = NdbToken::getTokenFormat(*tpos);
            if (!ndbformat.RefineField) continue;

            std::map<std::string, std::list<std::vector<std::string> > >::iterator
                ppos = _pdb_records.find(*tpos);
            if (ppos == _pdb_records.end()) continue;

            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[ndbformat.RefineField - 1].empty())
                     (*lpos)[ndbformat.RefineField - 1] = method;
            }
       }
}

#if 0
void Maxit::_removeEmptyRecord(const std::string& token, const std::set<unsigned int>& skip_field)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       std::vector<std::vector<std::string> > tmp_vector;
       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            bool found = false;
            for (unsigned int i = 1; i < lpos->size(); ++i) {
                 if (skip_field.find(i) != skip_field.end()) continue;
                 if (!(*lpos)[i].empty()) {
                      found = true;
                      break;
                 }
            }
            if (found) continue;
            lpos = ppos->second.erase(lpos);
            --lpos;
       }
}
#endif

void Maxit::_removeEmptyRecords()
{
       std::set<std::string> empty_tokens;
       empty_tokens.clear();

       for (std::map<std::string, std::list<std::vector<std::string> > >::iterator
            ppos = _pdb_records.begin(); ppos != _pdb_records.end(); ++ppos) {
            if (ppos->first == "NCSGRO" || ppos->first == "PDBRMK" ||
                ppos->first == "PREMRK" || ppos->first == "REMARK") continue;

            int start = 1;
            if (NdbToken::IsJrnlToken(ppos->first)) start = 4;

            const ndb_token_format& ndbformat = NdbToken::getTokenFormat(ppos->first);

            for (std::list<std::vector<std::string> >::iterator
                 lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                 bool found = false;
                 for (int i = start; i < ndbformat.NumField; ++i) {
                      if (i != (ndbformat.SeqField - 1) && i  != (ndbformat.ConField - 1)
                          && !(*lppos)[i].empty()) {
                           found = true;
                           break;
                      }
                 }
                 if (!found) {
                      lppos = ppos->second.erase(lppos);
                      --lppos;
                 }
            }

            if (ppos->second.empty()) empty_tokens.insert(ppos->first);
       }

       if (empty_tokens.empty()) return;

       for (std::set<std::string>::const_iterator
            spos = empty_tokens.begin(); spos != empty_tokens.end(); ++spos) {
            _pdb_records.erase(*spos);
       }
}

void Maxit::_updateRecords(const int& format)
{
       int length = 10;
       if (format == NDB_FILE_FORMAT_PDB) length = 9;

       for (std::map<std::string, std::list<std::vector<std::string> > >::iterator
            ppos = _pdb_records.begin(); ppos != _pdb_records.end(); ++ppos) {
            if (ppos->first == "NCSGRO" || ppos->first == "PDBRMK" ||
                ppos->first == "PREMRK" || ppos->first == "REMARK") continue;

            const ndb_token_format& ndbformat = NdbToken::getTokenFormat(ppos->first);

            int MatrixNo = 1;
            int SerialNo = 0;
            int field_no = ndbformat.SeqField - 1;
            for (std::list<std::vector<std::string> >::iterator
                 lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                 if (!ndbformat.MatrixType && ndbformat.SeqField) {
                      SerialNo++;
                      if ((*lppos)[field_no].empty())
                           (*lppos)[field_no] = String::IntToString(SerialNo);
                 } else if (ndbformat.MatrixType) {
                      int len = (*lppos)[0].size();
                      if (!isdigit((*lppos)[0][len - 1]))
                           (*lppos)[0] = (*lppos)[0] + String::IntToString(MatrixNo);
                      if (MatrixNo % 3 == 0) MatrixNo = 1;
                      else MatrixNo++;
                 }
                 int start = 1;
                 if (NdbToken::IsJrnlToken(ppos->first)) start = 4;
                 for (int i = start; i < ndbformat.NumField; ++i) {
                      if ((*lppos)[i].empty()) continue;
                      if (ndbformat.FieldList[i].FieldType == NDB_DATE_FIELD) {
                           convert_time_format(length, (*lppos)[i]);
                      }
                      if (format == NDB_FILE_FORMAT_PDB && ((ppos->first != "COMPND" &&
                         !ndbformat.FieldList[i].MixedCase && ppos->first != "SOURCE") ||
                         NdbToken::IsJrnlToken(ppos->first))) {
                           String::UpperCase((*lppos)[i]);
                      }
                 }
            }
       }
}

void Maxit::_collapseRecord(const std::string& token, const std::set<int>& field_no_set,
                            const std::string& delimiter)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end() || ppos->second.size() < 2) return;

       std::map<int, std::string> field_value_pair;
       field_value_pair.clear();
       for (std::list<std::vector<std::string> >::iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            for (std::set<int>::const_iterator
                 spos = field_no_set.begin(); spos != field_no_set.end(); ++spos) {
                 if ((*lppos)[*spos].empty()) (*lppos)[*spos] = "NULL";
                 std::map<int, std::string>::iterator
                     mpos = field_value_pair.find(*spos);
                 if (mpos == field_value_pair.end())
                      field_value_pair.insert(std::make_pair(*spos, (*lppos)[*spos]));
                 else mpos->second += delimiter + (*lppos)[*spos];
            }
       }

       std::vector<std::string>& Field = ppos->second.front();
       for (std::map<int, std::string>::const_iterator
            mpos = field_value_pair.begin(); mpos != field_value_pair.end(); ++mpos) {
            Field[mpos->first] = mpos->second;
       }
}

void Maxit::_checkRecord(const std::string& token)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) {
            printf("No record for %s\n", token.c_str());
            return;
       }

       for (std::list<std::vector<std::string> >::const_iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            for (std::vector<std::string>::const_iterator
                 vpos = lppos->begin(); vpos != lppos->end(); ++vpos) {
                 printf("-%s- ", vpos->c_str());
            }
            printf("\n");
       }
}

void Maxit::_get_remark_map(const std::string& token, const int& field_no, std::map<int,
                            std::vector<std::string> >& Remarks)
{
       Remarks.clear();

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       std::string id = "------1";

       std::vector<std::string> data;
       data.clear();

       for (std::list<std::vector<std::string> >::const_iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            if ((*lppos)[1] != id) {
                 if (!data.empty() && String::IsNumber(id))
                      Remarks.insert(std::make_pair(atoi(id.c_str()), data));
                 id = (*lppos)[1];
                 data.clear();
                 data.push_back((*lppos)[field_no]);
            } else data.push_back((*lppos)[field_no]);
       }
       if (!data.empty() && String::IsNumber(id))
            Remarks.insert(std::make_pair(atoi(id.c_str()), data));
}

void Maxit::_addNewRemark(const int& remark_no, const std::string& value)
{
       std::vector<std::string> FieldInfo;
       FieldInfo.clear();

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
       for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
       FieldInfo[0] = "REMARK";
       FieldInfo[1] = String::IntToString(remark_no);
       FieldInfo[2] = value;

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("REMARK");
       if (ppos == _pdb_records.end()) {
            std::list<std::vector<std::string> > rlist;
            rlist.clear();
            rlist.push_back(FieldInfo);
            _pdb_records.insert(std::make_pair("REMARK", rlist));
       } else ppos->second.push_back(FieldInfo);
}
