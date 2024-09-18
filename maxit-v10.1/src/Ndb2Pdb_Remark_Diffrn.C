/*
FILE:     Ndb2Pdb_Remark_Diffrn.C
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
#include <math.h>

#include "Maxit.h"
#include "NdbToken.h"
#include "Remark_define.h"
#include "Remark_extern.h"
#include "utillib.h"

void Maxit::_ndb_to_pdb_update_SCTYPE()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("SCTYPE");
       if (ppos == _pdb_records.end()) return;

       std::string cs;
       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            String::UpperCase((*lpos)[2], cs);
            if (cs == "X-RAY")
                 (*lpos)[2] = "X-RAY DIFFRACTION";
            else if (cs == "NEUTRON")
                 (*lpos)[2] = "NEUTRON DIFFRACTION";
            else if (cs == "ELECTRON")
                 (*lpos)[2] = "ELECTRON CRYSTALLOGRAPHY";
       }
}

void Maxit::_ndb_to_pdb_processing_shelx_and_tnt()
{
       std::string value;
       _getRecordFront("REFMET", 3, value);

       if (value.empty()) return;

       String::UpperCase(value);

       if (value.substr(0, 5) == "SHELX") {
            std::map<int, int> field_mapping;
            field_mapping.clear();
            field_mapping.insert(std::make_pair(7, 2));
            field_mapping.insert(std::make_pair(8, 1));
            _copyValue("RFACTR", "RNOCUT", field_mapping, false);
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
                ppos = _pdb_records.find("RFACTR");
            if (ppos != _pdb_records.end()) {
                 field_mapping.clear();
                 field_mapping.insert(std::make_pair(1, 3));
                 field_mapping.insert(std::make_pair(3, 4));
                 field_mapping.insert(std::make_pair(4, 5));
                 _copyValue("FREERF", "RNOCUT", field_mapping, false);
            }
       } else if (value.substr(0, 3) == "TNT" || value.substr(0, 6) == "NUCLSQ" ||
                  value.substr(0, 6) == "PROLSQ") {
            _getRecordFront("RFACTR", 4, value);
            if (value.empty() || fabs(atof(value.c_str())) < 0.1) {
                 _copyValue("RFACTR", 8, "RNOCUT", 1, false);
            }
       }
}

void Maxit::_ndb_to_pdb_processing_synchrotron()
{
       _ndb_to_pdb_update_DTMETH_DTMEAS_RADIAT();

       _ndb_to_pdb_update_WAVLEN();

       std::set<int> diffrn_ids;
       _ndb_to_pdb_get_diffrn_ids(diffrn_ids);
       if (diffrn_ids.empty()) return;

       std::string field_4, field_5;
       _getRecordFront("DTMEAS", 4, field_4);
       _getRecordFront("DTMEAS", 5, field_5);

       std::map<int, std::string> diffrn_sctype;
       _ndb_to_pdb_get_diffrn_value("SCTYPE", 1, 2, diffrn_sctype);

       std::map<std::string, std::set<int> > sctype_diffrn_ids;
       _ndb_to_pdb_get_sctype_diffrn_ids(diffrn_ids, diffrn_sctype, sctype_diffrn_ids);

       std::map<int, int> diffrn_crystal_ids;
       _ndb_to_pdb_get_diffrn_crystal_ids(diffrn_crystal_ids);

       std::map<int, std::string> crystal_id_ph;
       _ndb_to_pdb_get_diffrn_value("PHVAL", 2, 1, crystal_id_ph);

       // Update SCTYPE & PHVAL tokens
       std::set<int> ids;
       std::string cs, cs1;
       for (std::map<std::string, std::set<int> >::const_iterator
            mvpos = sctype_diffrn_ids.begin(); mvpos != sctype_diffrn_ids.end(); ++mvpos) {
            cs.clear();
            ids.clear();
            for (std::set<int>::const_iterator
                 spos = mvpos->second.begin(); spos != mvpos->second.end(); ++spos) {
                 if (!cs.empty()) cs += ",";
                 cs += String::IntToString(*spos);
                 std::map<int, int>::const_iterator
                     mpos = diffrn_crystal_ids.find(*spos);
                 if (mpos == diffrn_crystal_ids.end()) continue;
                 ids.insert(mpos->second);
            }
            cs1.clear();
            for (std::set<int>::const_iterator
                 spos = ids.begin(); spos != ids.end(); ++spos) {
                 std::map<int, std::string>::const_iterator
                     mpos = crystal_id_ph.find(*spos);
                 if (mpos == crystal_id_ph.end()) continue;
                 if (!cs1.empty()) cs1 += "; ";
                 cs1 += mpos->second;
            }

            _addNewRecord("SCTYPE");
            _updateRecordBack("SCTYPE", 1, cs);
            _updateRecordBack("SCTYPE", 2, mvpos->first);

            _addNewRecord("PHVAL");
            _updateRecordBack("PHVAL", 1, cs1);
            _updateRecordBack("PHVAL", 2, cs);
       }

       _ndb_to_pdb_update_200_align(sctype_diffrn_ids);

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
            ppos = _pdb_records.find("DTMEAS");
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            (*lpos)[4] = field_4;
            (*lpos)[5] = field_5;
       }
}

void Maxit::_ndb_to_pdb_update_DTMETH_DTMEAS_RADIAT()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
            ppos = _pdb_records.find("DTMETH");
       if (ppos == _pdb_records.end()) return;

       std::string site, beamline;
       std::map<int, std::string> field_value_pair;
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("DTMETH");
       const ndb_token_format& ndbformat1 = NdbToken::getTokenFormat("DTMEAS");
       const ndb_token_format& ndbformat2 = NdbToken::getTokenFormat("RADIAT");

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if (!String::IsEqual((*lpos)[4], "SYNCHROTRON", Char::eCASE_INSENSITIVE)) {
                 (*lpos)[1] = "N";
                 continue;
            }

            (*lpos)[1] = "Y";
            _ndb_to_pdb_get_beamline((*lpos)[ndbformat.SeqField-1], site, beamline);

            if (!site.empty() || !beamline.empty()) {
                 field_value_pair.clear();
                 if (!site.empty()) field_value_pair.insert(std::make_pair(2, site));
                 if (!beamline.empty()) field_value_pair.insert(std::make_pair(3, beamline));
                 _updateValue("DTMEAS", ndbformat1.SeqField-1, (*lpos)[ndbformat.SeqField-1],
                              field_value_pair);
            }
            if (_pdb_records.find("RADIAT") != _pdb_records.end()) {
                 field_value_pair.clear();
                 field_value_pair.insert(std::make_pair(5, ""));
                 _updateValue("RADIAT", ndbformat2.SeqField-1, (*lpos)[ndbformat.SeqField-1],
                              field_value_pair);
            }
            _getValue("DTMEAS", ndbformat1.SeqField-1, (*lpos)[ndbformat.SeqField-1], 2, site);
            (*lpos)[4] = site;
       }
}

void Maxit::_ndb_to_pdb_update_WAVLEN()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
            ppos = _pdb_records.find("DTMETH");
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
            qpos = _pdb_records.find("WAVLEN");
       if (ppos != _pdb_records.end() && qpos == _pdb_records.end()) {
            std::string value = "";
            for (std::list<std::vector<std::string> >::const_iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[3].empty()) continue;
                 if (!value.empty()) value += ", ";
                 value += (*lpos)[3];
            }
            if (!value.empty()) {
                 _addNewRecord("WAVLEN");
                 _updateRecordBack("WAVLEN", 3, value);
            }
       }

       qpos = _pdb_records.find("WAVLEN");
       if (qpos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = qpos->second.begin(); lpos != qpos->second.end(); ++lpos) {
            if ((*lpos)[2].empty() && !(*lpos)[3].empty()) (*lpos)[2] = (*lpos)[3];
       }
}

void Maxit::_ndb_to_pdb_get_beamline(const std::string& key_value, std::string& site,
                                      std::string& beamline)
{
       site.clear(); beamline.clear();

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
            ppos = _pdb_records.find("RADIAT");
       if (ppos == _pdb_records.end()) return;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("RADIAT");

       std::vector<std::string> data;
       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[ndbformat.SeqField-1] == key_value) {
                 String::UpperCase((*lpos)[5]);
                 get_wordarray_delimit_by_string(data, (*lpos)[5], "BEAMLINE");
                 if (!data.empty()) {
                      site = data[0];
                      String::StripAndCompressWs(site);
                      if (data.size() > 1) {
                           beamline = data[1];
                           String::StripAndCompressWs(beamline);
                      }
                 }
                 (*lpos)[5].clear();
                 return;
            }
       }
}

void Maxit::_ndb_to_pdb_get_diffrn_ids(std::set<int>& diffrn_ids)
{
       diffrn_ids.clear();

       std::set<int> ids;
       ids.clear();
       for (int i = 0; i < NUM_200_ALIGN; i++) {
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
                ppos = _pdb_records.find(_remark_200_align[i].TokenName);
            if (ppos == _pdb_records.end()) continue;

            for (std::list<std::vector<std::string> >::const_iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[_remark_200_align[i].Field_ID].empty()) continue;
                 int id = atoi((*lpos)[_remark_200_align[i].Field_ID].c_str());
                 diffrn_ids.insert(id);
                 if (strcmp(_remark_200_align[i].TokenName, "DTTEMP")) ids.insert(id);
            }
       }
       if (!ids.empty()) diffrn_ids = ids;
}

void Maxit::_ndb_to_pdb_get_diffrn_value(const std::string& token, const int& key_field,
                       const int& value_field, std::map<int, std::string>& diffrn_value)
{
       diffrn_value.clear();

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(token);
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[key_field].empty() || (*lpos)[value_field].empty()) continue;
            diffrn_value.insert(std::make_pair(atoi((*lpos)[key_field].c_str()),
                                                     (*lpos)[value_field]));
       }
       _pdb_records.erase(token);
}

void Maxit::_ndb_to_pdb_get_sctype_diffrn_ids(const std::set<int>& diffrn_ids, const
                        std::map<int, std::string>& diffrn_sctype, std::map<std::string,
                        std::set<int> >& sctype_diffrn_ids)
{
       sctype_diffrn_ids.clear();

       if (!diffrn_sctype.empty()) {
            std::set<int> ids;
            for (std::set<int>::const_iterator
                 spos = diffrn_ids.begin(); spos != diffrn_ids.end(); ++spos) {
                 std::map<int, std::string>::const_iterator
                     mpos = diffrn_sctype.find(*spos);
                 if (mpos == diffrn_sctype.end()) mpos = diffrn_sctype.begin();
                 std::map<std::string, std::set<int> >::iterator
                     mvpos = sctype_diffrn_ids.find(mpos->second);
                 if (mvpos != sctype_diffrn_ids.end()) {
                      mvpos->second.insert(*spos);
                 } else {
                      ids.clear();
                      ids.insert(*spos);
                      sctype_diffrn_ids.insert(std::make_pair(mpos->second, ids));
                 }
            }
       } else {
            std::string cs = "X-RAY DIFFRACTION";
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
                ppos = _pdb_records.find("SCTYPE");
            if (ppos != _pdb_records.end()) {
                 for (std::list<std::vector<std::string> >::const_iterator
                      lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                      if ((*lpos)[2] == "X-RAY DIFFRACTION" ||
                          (*lpos)[2] == "FIBER DIFFRACTION" ||
                          (*lpos)[2] == "NEUTRON DIFFRACTION" ||
                          (*lpos)[2] == "ELECTRON CRYSTALLOGRAPHY") {
                           cs = (*lpos)[2];
                           break;
                      }
                 }
            }
            sctype_diffrn_ids.insert(std::make_pair(cs, diffrn_ids));
       }
}

void Maxit::_ndb_to_pdb_get_diffrn_crystal_ids(std::map<int, int>& diffrn_crystal_ids)
{
       diffrn_crystal_ids.clear();

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("DTTEMP");
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[4].empty() || (*lpos)[5].empty()) continue;
            diffrn_crystal_ids.insert(std::make_pair(atoi((*lpos)[4].c_str()),
                                                     atoi((*lpos)[5].c_str())));
       }
}

void Maxit::_ndb_to_pdb_update_200_align(const std::map<std::string, std::set<int> >&
                                         sctype_diffrn_ids)
{
       std::map<int, std::vector<std::string> > records;
       std::vector<std::string> data;
       std::string cs;

       for (int i = 0; i < NUM_200_ALIGN; i++) {
            records.clear();
            std::map<std::string, std::list<std::vector<std::string> > >::iterator
                ppos = _pdb_records.find(_remark_200_align[i].TokenName);
            if (ppos != _pdb_records.end()) {
                for (std::list<std::vector<std::string> >::iterator
                     lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                     int id = atoi((*lpos)[_remark_200_align[i].Field_ID].c_str());
                     records.insert(std::make_pair(id, *lpos));
                }
                _pdb_records.erase(_remark_200_align[i].TokenName);
            }

            const ndb_token_format& ndbformat =
                      NdbToken::getTokenFormat(_remark_200_align[i].TokenName);
            for (std::map<std::string, std::set<int> >::const_iterator mvpos =
                 sctype_diffrn_ids.begin(); mvpos != sctype_diffrn_ids.end(); ++mvpos) {
                 cs.clear();
                 data.clear();
                 for (int j = 0; j < ndbformat.NumField; ++j) data.push_back("");
                 for (std::set<int>::const_iterator
                      spos = mvpos->second.begin(); spos != mvpos->second.end(); ++spos) {
                      if (!cs.empty()) cs += ",";
                      cs += String::IntToString(*spos);

                      std::map<int, std::vector<std::string> >::const_iterator
                          mvpos1 = records.find(*spos);
                      for (unsigned int k = 0; k < data.size(); ++k) {
                           if (!data[k].empty()) data[k] += "; ";
                           if (mvpos1 != records.end() && !mvpos1->second[k].empty())
                                data[k] += mvpos1->second[k];
                           else data[k] += "NULL";
                      }
                 }

                 _addNewRecord(_remark_200_align[i].TokenName);
                 _updateRecordBack(_remark_200_align[i].TokenName,
                                   _remark_200_align[i].Field_ID, cs);
                 for (int j = 0; j < _remark_200_align[i].NumField; j++) {
                      _updateRecordBack(_remark_200_align[i].TokenName,
                                        _remark_200_align[i].Field_No[j],
                                        data[_remark_200_align[i].Field_No[j]]);
                 }
            }
       }
}
