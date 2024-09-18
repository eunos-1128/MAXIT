/*
FILE:     Ndb2Pdb_Remark_NCSTLS.C
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

#include "Maxit.h"
#include "NdbToken.h"
#include "utillib.h"

#define NUM_TML 6

static const char *_tml_tokens[NUM_TML] = {
       "RFTPOS", "RFMPOS", "RFLPOS", "RFTTHR", "RFMTHR", "RFLTHR"
};

static const char *_tml_types[NUM_TML] = {
       "tight positional", "medium positional", "loose positional",
       "tight thermal",    "medium thermal",    "loose thermal"
};

int Maxit::_ndb_to_pdb_update_TLSGRO()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("TLSGRO");
       if (ppos == _pdb_records.end()) return 0;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            int count = 0;
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
                qpos = _pdb_records.find("TLSRNG");
            if (qpos != _pdb_records.end()) {
                 for (std::list<std::vector<std::string> >::const_iterator
                      lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                      if ((*lpos)[1] == (*lqpos)[1]) count++;
                 }
            }
            (*lpos)[2] = String::IntToString(count);
       }

       return ppos->second.size();
}

int Maxit::_ndb_to_pdb_processing_refmac5()
{
       std::vector<std::vector<std::string> > ncs_data;
       ncs_data.clear();

       std::vector<std::string> data;
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("NCSGRO");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 bool found = false;
                 for (std::vector<std::vector<std::string> >::iterator
                      vpos = ncs_data.begin(); vpos != ncs_data.end(); ++vpos) {
                      if ((*vpos)[0] == (*lpos)[4]) {
                           found = true;
                           if ((*vpos)[1].find((*lpos)[3]) == std::string::npos) {
                                (*vpos)[1] += " " + (*lpos)[3];
                           }
                      }
                 }
                 if (found) continue;

                 data.clear();
                 data.push_back((*lpos)[4]);
                 data.push_back((*lpos)[3]);
                 ncs_data.push_back(data); 
            }
       }

       if (!ncs_data.empty()) {
            _pdb_records.erase("NCSGRO");
            for (std::vector<std::vector<std::string> >::iterator
                 vpos = ncs_data.begin(); vpos != ncs_data.end(); ++vpos) {
                 _addNewRecord("NCSGRO");
                 _updateRecordBack("NCSGRO", 3, (*vpos)[1]);
                 _updateRecordBack("NCSGRO", 4, (*vpos)[0]);
            }
       }

       int num_cns = 0;
       ppos = _pdb_records.find("NCSGRO");
       if (ppos != _pdb_records.end()) {
            num_cns = ppos->second.size();

            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 int count = 0;
                 std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
                     qpos = _pdb_records.find("CCPNCS");
                 if (qpos != _pdb_records.end()) {
                      for (std::list<std::vector<std::string> >::const_iterator
                           lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                           if ((*lpos)[4] == (*lqpos)[1] && atoi((*lqpos)[6].c_str()) > count)
                                count = atoi((*lqpos)[6].c_str());
                      }
                 }
                 (*lpos)[2] = String::IntToString(count);
            }
       }

       ppos = _pdb_records.find("TMLPTL");
       if (ppos != _pdb_records.end()) {
            _pdb_records.erase("NCSLOC");
            for (int i = 0; i < NUM_TML; ++i) _pdb_records.erase(_tml_tokens[i]);

            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 (*lpos)[1] = (*lpos)[7];
                 String::LowerCase((*lpos)[6]);
                 String::StripAndCompressWs((*lpos)[6]);
                 for (int i = 0; i < NUM_TML; ++i) {
                      if ((*lpos)[6] != _tml_types[i]) continue;

                      const ndb_token_format& ndbformat =
                                           NdbToken::getTokenFormat(_tml_tokens[i]);
                      _addNewRecord(_tml_tokens[i]);
                      for (int j = 1; j < ndbformat.NumField; ++j) {
                           if (j == 2 && (*lpos)[j].empty())
                                _updateRecordBack(_tml_tokens[i], j, " ");
                           else _updateRecordBack(_tml_tokens[i], j, (*lpos)[j]);
                      }
                      break;
                 }

                 if ((*lpos)[6] != "interatomic distance") continue;

                 bool found = false;
                 std::map<std::string, std::list<std::vector<std::string> > >::iterator
                     qpos = _pdb_records.find("NCSLOC");
                 if (qpos != _pdb_records.end()) {
                      for (std::list<std::vector<std::string> >::iterator
                           lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                           if (((*lqpos)[1] == (*lpos)[1]) && (*lqpos)[5].empty()) {
                                found = true;
                                (*lqpos)[5] = (*lpos)[2];
                           }
                      }
                 }
                 if (found) continue;

                 _addNewRecord("NCSLOC");
                 _updateRecordBack("NCSLOC", 1, (*lpos)[1]);
                 _updateRecordBack("NCSLOC", 2, (*lpos)[2]);
                 _updateRecordBack("NCSLOC", 8, (*lpos)[3]);
                 _updateRecordBack("NCSLOC", 9, (*lpos)[4]);
                 _updateRecordBack("NCSLOC", 10, (*lpos)[5]);
            }
       }

       ppos = _pdb_records.find("TLSRNG");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[2].empty()) (*lpos)[2] = " ";
                 if ((*lpos)[4].empty()) (*lpos)[4] = " ";
            }
       }

       ppos = _pdb_records.find("CCPNCS");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 std::map<std::string, std::list<std::vector<std::string> > >::iterator
                     qpos = _pdb_records.find("NCSLOC");
                 if (qpos != _pdb_records.end()) {
                      for (std::list<std::vector<std::string> >::iterator
                           lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                           if ((*lqpos)[1] != (*lpos)[1]) continue;
                           if ((*lqpos)[2] == (*lpos)[2]) {
                                (*lqpos)[3] = (*lpos)[3];
                                (*lqpos)[4] = (*lpos)[5];
                           } else if ((*lqpos)[5] == (*lpos)[2]) {
                                (*lqpos)[6] = (*lpos)[3];
                                (*lqpos)[7] = (*lpos)[5];
                           }
                      }
                 }
                 if ((*lpos)[2].empty()) (*lpos)[2] = " ";
                 if ((*lpos)[4].empty()) (*lpos)[4] = " ";
            }
       }

       ppos = _pdb_records.find("NCSLOC");
       if (ppos != _pdb_records.end()) num_cns = ppos->second.size();

       return num_cns;
}

int Maxit::_ndb_to_pdb_processing_phenix(const bool& is_buster)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("TLSRNG");
       if (ppos != _pdb_records.end()) {
            std::map<int, std::string> field_value_pair;
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if (!(*lpos)[2].empty() && !(*lpos)[3].empty() && !(*lpos)[4].empty() &&
                     !(*lpos)[5].empty() && (*lpos)[11].empty() && is_buster) { 
                      (*lpos)[11] = " " + FormattedString((*lpos)[2], 2, false, true)
                                  + " " + FormattedString((*lpos)[3], 4, false, true)
                                  + "   " + FormattedString((*lpos)[4], 2, false, true)
                                  + " " + FormattedString((*lpos)[5], 4, false, true);
                 }
                 if ((*lpos)[11].empty()) continue;

                 field_value_pair.clear();
                 field_value_pair.insert(std::make_pair(7, (*lpos)[11]));
                 _updateValue("TLSGRO", 1, (*lpos)[1], field_value_pair);
            }
       }

       ppos = _pdb_records.find("CCPNCS");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[15].empty()) continue;
                 bool found = false;
                 std::map<std::string, std::list<std::vector<std::string> > >::iterator
                     qpos = _pdb_records.find("NCSGRO");
                 if (qpos != _pdb_records.end()) {
                      for (std::list<std::vector<std::string> >::iterator
                           lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                           if ((*lpos)[1] == (*lqpos)[4] && (*lpos)[14] == (*lqpos)[1]) {
                                (*lqpos)[3] = (*lpos)[15];
                                found = true;
                                break;
                           }
                      }
                 }
                 if (found) continue;

                 _addNewRecord("NCSGRO");
                 _updateRecordBack("NCSGRO", 1, (*lpos)[14]);
                 _updateRecordBack("NCSGRO", 3, (*lpos)[15]);
                 _updateRecordBack("NCSGRO", 4, (*lpos)[1]);
            }
       }

       std::vector<std::string> Ens_IDs;
       Ens_IDs.clear();

       std::string cs = "";

       ppos = _pdb_records.find("NCSGRO");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 bool found = false;
                 for (std::vector<std::string>::const_iterator
                      vpos = Ens_IDs.begin(); vpos != Ens_IDs.end(); ++vpos) {
                      if (*vpos == (*lpos)[4]) {
                           found = true;
                           break;
                      }
                 }
                 if (!found) Ens_IDs.push_back((*lpos)[4]);

                 if ((*lpos)[1] == "1") {
                      cs = (*lpos)[3]; 
                      std::string::size_type p = (*lpos)[3].find("REFERENCE SELECTION:");
                      if (p != std::string::npos) cs = (*lpos)[3].substr(p + 20);
                      String::StripAndCompressWs(cs);
                 } else {
                      std::string cs2 = (*lpos)[3];
                      std::string::size_type p = (*lpos)[3].find("SELECTION :");
                      if (p != std::string::npos) cs2 = (*lpos)[3].substr(p + 11);
                      String::StripAndCompressWs(cs2);

                      _addNewRecord("PHENCS");
                      _updateRecordBack("PHENCS", 1, (*lpos)[4]);
                      _updateRecordBack("PHENCS", 3, cs);
                      _updateRecordBack("PHENCS", 4, cs2);
                      
                      std::map<std::string, std::list<std::vector<std::string> > >::iterator
                          qpos = _pdb_records.find("TMLPTL");
                      if (qpos != _pdb_records.end()) {
                           for (std::list<std::vector<std::string> >::iterator
                                lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                                if ((*lpos)[1] == (*lqpos)[1] && (*lpos)[4] == (*lqpos)[7]) {
                                     _updateRecordBack("PHENCS", 5, (*lqpos)[3]);
                                     _updateRecordBack("PHENCS", 6, (*lqpos)[4]);
                                     break;
                                }
                           }
                      }
                 }
            }
       }

       _pdb_records.erase("NCSGRO");
       for (std::vector<std::string>::const_iterator
            vpos = Ens_IDs.begin(); vpos != Ens_IDs.end(); ++vpos) {
            _addNewRecord("NCSGRO");
            _updateRecordBack("NCSGRO", 4, *vpos);
       }

       ppos = _pdb_records.find("PHENCS");
       if (ppos != _pdb_records.end()) {
            cs = "xxx";
            int i = 1;
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[1] != cs) {
                      i = 1;
                      cs = (*lpos)[1];
                 }
                 (*lpos)[2] = String::IntToString(i);
                 i++;
            } 
       }

       return (Ens_IDs.size());
}
