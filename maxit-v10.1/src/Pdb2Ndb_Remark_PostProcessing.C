/*
FILE:     Pdb2Ndb_Remark_PostProcessing.C
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
#include "Remark_define.h"
#include "Remark_extern.h"
#include "TypeDef.h"
#include "utillib.h"

void Maxit::_pdb_to_ndb_postprocessing()
{
       _pdb_to_ndb_update_EXPDTA();

       bool is_refmac5 = false;
       bool is_phenix = false;
       bool is_cns_xplor = false;
       bool is_buster = false;

       _get_refinement_program(is_refmac5, is_phenix, is_cns_xplor, is_buster);

       if (is_refmac5) {
            if (!_pdb_to_ndb_update_local_NCS()) _pdb_to_ndb_update_NCSGRO();
            _pdb_to_ndb_update_CCPNCS();
       }

       if (is_phenix) _pdb_to_ndb_update_PHENCS();

       if (is_phenix || is_buster) _pdb_to_ndb_update_TLSGRO_and_TLSRNG();

       _mixCase("REF", 4);
       _mixCase("AUTH", 4);
       _mixCase("TITL", 4);
       _mixCase("PUBL", 4);
       _mixCase("EDIT", 4);
       _removeWhiteSpace("DOI", 4);
       _mixCase("AUTHOR", 2);

       _pdb_to_ndb_update_RFACTR();

       // _pdb_to_ndb_update_200_align();

       _pdb_to_ndb_update_WAVLEN();

       _pdb_to_ndb_update_DTMETH_DTMEAS_RADIAT();

       _pdb_to_ndb_update_EMDATA();

       _pdb_to_ndb_update_TLSGRO();

       _pdb_to_ndb_update_TWIN();

       _pdb_to_ndb_update_SOLMOD();

       _removeEmptyRecords();
       _updateRecords(NDB_FILE_FORMAT_NDB);
}

void Maxit::_pdb_to_ndb_update_EXPDTA()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("EXPDTA");
       if (ppos == _pdb_records.end()) {
            for (std::map<int, std::string>::const_iterator mpos = _TypeExpMapping.begin(); mpos != _TypeExpMapping.end(); ++mpos) {
                 if (!(_experiment_type & mpos->first)) continue;
                 _addNewRecord("EXPDTA");
                 _updateRecordBack("EXPDTA", 2, mpos->second);
            }
       }

       std::string value;
       _getRecordFront("EXPDTA", 2, value);
       if (!value.empty()) _updateRefineMethod(value);
}

bool Maxit::_pdb_to_ndb_update_local_NCS()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("NCSLOC");
       if (ppos == _pdb_records.end()) return false;

       int count = 0;
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            count++;

            _addNewRecord("TMLPTL");
            _updateRecordBack("TMLPTL", 1, "1");
            _updateRecordBack("TMLPTL", 2, (*lpos)[2]);
            _updateRecordBack("TMLPTL", 3, (*lpos)[8]);
            _updateRecordBack("TMLPTL", 4, (*lpos)[9]);
            _updateRecordBack("TMLPTL", 5, (*lpos)[10]);
            _updateRecordBack("TMLPTL", 6, "interatomic distance");
            _updateRecordBack("TMLPTL", 7, String::IntToString(count));

            _addNewRecord("TMLPTL");
            _updateRecordBack("TMLPTL", 1, "2");
            _updateRecordBack("TMLPTL", 2, (*lpos)[5]);
            _updateRecordBack("TMLPTL", 3, (*lpos)[8]);
            _updateRecordBack("TMLPTL", 4, (*lpos)[9]);
            _updateRecordBack("TMLPTL", 5, (*lpos)[10]);
            _updateRecordBack("TMLPTL", 6, "interatomic distance");
            _updateRecordBack("TMLPTL", 7, String::IntToString(count));

            _addNewRecord("CCPNCS");
            _updateRecordBack("CCPNCS", 1, String::IntToString(count));
            _updateRecordBack("CCPNCS", 2, (*lpos)[2]);
            _updateRecordBack("CCPNCS", 3, (*lpos)[3]);
            _updateRecordBack("CCPNCS", 4, (*lpos)[2]);
            _updateRecordBack("CCPNCS", 5, (*lpos)[4]);
            _updateRecordBack("CCPNCS", 6, "0");
            _updateRecordBack("CCPNCS", 7, "0");
            _updateRecordBack("CCPNCS", 14, "1");

            _addNewRecord("CCPNCS");
            _updateRecordBack("CCPNCS", 1, String::IntToString(count));
            _updateRecordBack("CCPNCS", 2, (*lpos)[5]);
            _updateRecordBack("CCPNCS", 3, (*lpos)[6]);
            _updateRecordBack("CCPNCS", 4, (*lpos)[5]);
            _updateRecordBack("CCPNCS", 5, (*lpos)[7]);
            _updateRecordBack("CCPNCS", 6, "0");
            _updateRecordBack("CCPNCS", 7, "0");
            _updateRecordBack("CCPNCS", 14, "2");

            _addNewRecord("NCSGRO");
            _updateRecordBack("NCSGRO", 1, "1");
            _updateRecordBack("NCSGRO", 3, (*lpos)[2]);
            _updateRecordBack("NCSGRO", 4, String::IntToString(count));

            _addNewRecord("NCSGRO");
            _updateRecordBack("NCSGRO", 1, "2");
            _updateRecordBack("NCSGRO", 3, (*lpos)[5]);
            _updateRecordBack("NCSGRO", 4, String::IntToString(count));
       }

       _updateRecordFront("NCSTLS", 1, String::IntToString(count));

       return true;
}

void Maxit::_pdb_to_ndb_update_NCSGRO()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("NCSGRO");
       if (ppos == _pdb_records.end()) return;

       std::vector<std::string> data, data1;

       std::vector<std::vector<std::string> > ncs_ens;
       ncs_ens.clear();
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            get_wordarray(data, (*lpos)[3], " ");
            for (unsigned int i = 0; i < data.size(); ++i) {
                 data1.clear();
                 data1.push_back(String::IntToString(i+1));
                 data1.push_back(data[i]);
                 data1.push_back((*lpos)[4]);
                 ncs_ens.push_back(data1);
            }
       }
       if (!ncs_ens.empty()) {
            _pdb_records.erase("NCSGRO");
            for (std::vector<std::vector<std::string> >::const_iterator
                 pos = ncs_ens.begin(); pos != ncs_ens.end(); ++pos) {
                 _addNewRecord("NCSGRO");
                 _updateRecordBack("NCSGRO", 1, (*pos)[0]);
                 _updateRecordBack("NCSGRO", 3, (*pos)[1]);
                 _updateRecordBack("NCSGRO", 4, (*pos)[2]);
            }
       }
}

void Maxit::_pdb_to_ndb_update_CCPNCS()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("NCSGRO");
       if (ppos == _pdb_records.end()) return;

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           qpos = _pdb_records.find("CCPNCS");
       if (qpos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
            for (std::list<std::vector<std::string> >::const_iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lqpos)[1] == (*lpos)[4] && (*lqpos)[2] == (*lpos)[3]) {
                      (*lqpos)[14] = (*lpos)[1];
                 }
            }
       }
}

void Maxit::_pdb_to_ndb_update_PHENCS()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("PHENCS");
       if (ppos == _pdb_records.end()) return;

       std::string cs = "xxx";
       int count = 1;
       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[1] != cs) {
                 count = 1;
                 cs = (*lpos)[1];
            }
            (*lpos)[2] = String::IntToString(count);
            count++;
       }
}

void Maxit::_pdb_to_ndb_update_TLSGRO_and_TLSRNG()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("TLSGRO");
       if (ppos == _pdb_records.end()) return;

       int i = 0;
       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            i++;
            (*lpos)[1] =  String::IntToString(i);
            if ((*lpos)[7].empty()) continue;

            bool found = false;
            std::map<std::string, std::list<std::vector<std::string> > >::iterator
                qpos = _pdb_records.find("TLSRNG");
            if (qpos != _pdb_records.end()) {
                 for (std::list<std::vector<std::string> >::iterator
                      lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                      if ((*lpos)[1] == (*lqpos)[1]) {
                           (*lqpos)[11] = (*lpos)[7];
                           found = true;
                           break;
                      }
                 }
            }

            if (found) continue;

            _addNewRecord("TLSRNG");
            _updateRecordBack("TLSRNG", 1, (*lpos)[1]);
            _updateRecordBack("TLSRNG", 6, (*lpos)[1]);
            _updateRecordBack("TLSRNG", 11, (*lpos)[7]);
       }
}

void Maxit::_pdb_to_ndb_update_RFACTR()
{
       std::string program;
       _getRecordFront("REFMET", 3, program);
       if (program.empty()) return;

       String::UpperCase(program);

       if (program.substr(0, 3) == "CNS" || program.substr(0, 6) == "X-PLOR") {
            std::string value;
            _getRecordFront("RFACTR", 10, value);
            if (!value.empty()) _updateRecordFront("RFACTR", 7, value);
       } else if (program.substr(0, 5) == "SHELX")
            _pdb_to_ndb_update_RFACTR_and_FREERF_for_SHELX();
       else if (program.substr(0, 3) == "TNT" || program.substr(0, 6) == "NUCLSQ" ||
                program.substr(0, 6) == "PROLSQ")
            _pdb_to_ndb_update_RFACTR_for_TNT_NUCLSQ_PROLSQ();
}

void Maxit::_pdb_to_ndb_update_RFACTR_and_FREERF_for_SHELX()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("RNOCUT");
       if (ppos == _pdb_records.end()) return;

       const std::vector<std::string>& Field = ppos->second.front();

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           qpos = _pdb_records.find("RFACTR");
       if (qpos != _pdb_records.end()) {
            std::vector<std::string>& Field1 = qpos->second.front();
            Field1[7] = Field[2];
            Field1[8] = Field[1];
            Field1[13] = Field[8];
       } else {
            _addNewRecord("RFACTR");
            _updateRecordBack("RFACTR", 7, Field[2]);
            _updateRecordBack("RFACTR", 8, Field[1]);
            _updateRecordBack("RFACTR", 13, Field[8]);
       }

       qpos = _pdb_records.find("FREERF");
       if (qpos != _pdb_records.end()) {
            std::vector<std::string>& Field1 = qpos->second.front();
            Field1[1] = Field[3];
            Field1[3] = Field[4];
            Field1[4] = Field[5];
            Field1[7] = Field[8];
       } else {
            _addNewRecord("FREERF");
            _updateRecordBack("FREERF", 1, Field[3]);
            _updateRecordBack("FREERF", 3, Field[4]);
            _updateRecordBack("FREERF", 4, Field[5]);
            _updateRecordBack("FREERF", 7, Field[8]);
       }
}

void Maxit::_pdb_to_ndb_update_RFACTR_for_TNT_NUCLSQ_PROLSQ()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("RNOCUT");
       if (ppos == _pdb_records.end()) return;

       const std::vector<std::string>& Field = ppos->second.front();


       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           qpos = _pdb_records.find("RFACTR");
       if (qpos != _pdb_records.end()) {
            std::vector<std::string>& Field1 = qpos->second.front();
            if (Field1[8].empty() && (Field1[4].empty() || (!Field1[4].empty() &&
                fabs(atof(Field1[4].c_str())) < 0.1))) {
                 Field1[8] = Field[1];
            }
       } else {
            _addNewRecord("RFACTR");
            _updateRecordBack("RFACTR", 4, "0.000");
            _updateRecordBack("RFACTR", 8, Field[1]);
       }
}

int Maxit::_pdb_to_ndb_update_200_align(const int& start_id)
{
       std::vector<std::string> data;

       int finish_id = start_id;
       for (int i = 0; i < NUM_200_ALIGN; ++i) {
            std::map<std::string, std::list<std::vector<std::string> > >::iterator ppos = _pdb_records.find(_remark_200_align[i].TokenName);
            if (ppos == _pdb_records.end()) continue;

            std::vector<std::string>& Field = ppos->second.front();

            for (int j = 0; j < _remark_200_align[i].NumField; ++j) {
                 if (Field[_remark_200_align[i].Field_No[j]].empty()) continue;
                 get_wordlist_from_string_separated_by_delimits(data, Field[_remark_200_align[i].Field_No[j]], ";", "", true);
                 for (unsigned int k = 0; k < data.size(); ++k) {
                      std::string cs = data[k];
                      if (!strcmp(_remark_200_align[i].TokenName, "DTMEAS") && _remark_200_align[i].Field_No[j] == 1)
                           String::RemoveWhiteSpace(data[k], cs);
                      bool found = false;
                      for (std::list<std::vector<std::string> >::iterator lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                           if ((*lpos)[_remark_200_align[i].Field_ID].empty() ||
                               atoi((*lpos)[_remark_200_align[i].Field_ID].c_str()) == (int) (k + start_id)) {
                                (*lpos)[_remark_200_align[i].Field_ID] = String::IntToString(k + start_id);
                                (*lpos)[_remark_200_align[i].Field_No[j]] = cs;
                                if ((int) (k + start_id) > finish_id) finish_id = k + start_id;
                                found = true;
                                break;
                           }
                      }
                      if (found) continue;

                      _addNewRecord(_remark_200_align[i].TokenName);
                      _updateRecordBack(_remark_200_align[i].TokenName, _remark_200_align[i].Field_ID, String::IntToString(k + start_id));
                      _updateRecordBack(_remark_200_align[i].TokenName, _remark_200_align[i].Field_No[j], cs);
                      if ((int) (k + start_id) > finish_id) finish_id = k + start_id;
                 }
            }

            ppos = _pdb_records.find(_remark_200_align[i].TokenName);
            if (ppos->second.size() != 1) continue;
            Field = ppos->second.front();
            if (Field[_remark_200_align[i].Field_ID].empty()) Field[_remark_200_align[i].Field_ID] = String::IntToString(start_id);
       }
       return finish_id;
}

void Maxit::_pdb_to_ndb_update_WAVLEN()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("WAVLEN");
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[2].find(",") != std::string::npos) {
                 (*lpos)[3] = (*lpos)[2];
                 (*lpos)[2].clear();
            }
       }
}

void Maxit::_pdb_to_ndb_update_DTMETH_DTMEAS_RADIAT()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator ppos = _pdb_records.find("DTMETH");
       if (ppos == _pdb_records.end()) return;

       std::map<int, std::string> field_value_pair;
       for (std::list<std::vector<std::string> >::iterator lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[1] != "Y") continue;

            std::string beamline = "";
            std::string Synchrotron_Site = (*lpos)[4];
            (*lpos)[4] = "SYNCHROTRON";
            if (Synchrotron_Site.find(" BEAMLINE") != std::string::npos) {
                 std::string::size_type p = Synchrotron_Site.find(" BEAMLINE");
                 if (p != std::string::npos) Synchrotron_Site.erase(p);
            }

            bool found = false;
            std::map<std::string, std::list<std::vector<std::string> > >::iterator qpos = _pdb_records.find("DTMEAS");
            if (qpos != _pdb_records.end()) {
                 for (std::list<std::vector<std::string> >::iterator kpos = qpos->second.begin(); kpos != qpos->second.end(); ++kpos) {
                      if ((*kpos)[6] == (*lpos)[6]) {
                           found = true;
                           (*kpos)[2] = Synchrotron_Site;
                           if (!(*kpos)[2].empty() && !(*kpos)[3].empty())
                                beamline = (*kpos)[2] + " BEAMLINE " + (*kpos)[3];
                           else if (!(*kpos)[2].empty())
                                beamline = (*kpos)[2];
                           else if (!(*kpos)[3].empty())
                                beamline = "BEAMLINE " + (*kpos)[3];
                           break;
                      }
                 }
            }
            if (!found) {
                 _addNewRecord("DTMEAS");
                 _updateRecordBack("DTMEAS", 2, Synchrotron_Site);
                 beamline = Synchrotron_Site;
            }

            String::StripAndCompressWs(beamline);
            if (beamline.empty()) continue;

            field_value_pair.clear();
            field_value_pair.insert(std::make_pair(5, beamline));
            _updateValue("RADIAT", 1, (*lpos)[6], field_value_pair);
       }
}

void Maxit::_pdb_to_ndb_update_EMDATA()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("EMDATA");
       if (ppos == _pdb_records.end()) return;

       const std::vector<std::string>& Field = ppos->second.front();

       std::vector<std::vector<std::string> > all_data;
       all_data.clear();

       std::vector<std::string> data;
       for (int i = 1; i < 19; ++i) { 
            get_wordlist_from_string_separated_by_delimits(data, Field[i], ";", "", false);
            all_data.push_back(data);
       }

       int max_number = 0;
       for (std::vector<std::vector<std::string> >::const_iterator
            pos = all_data.begin(); pos != all_data.end(); ++pos) {
            if ((int) pos->size() > max_number)
                 max_number = pos->size();
       }

       _pdb_records.erase("EMDATA");
       
       for (int i = 0; i < max_number; ++i) {
            _addNewRecord("EMDATA");
            for (unsigned int j = 0; j < all_data.size(); ++j) {
                 if (i >= (int) all_data[j].size()) continue;
                 if (j == 11) {
                      // re-scale _em_image_recording.avg_electron_dose_per_image: divide 100 for PDB -> CIF
                      double f = atof(all_data[j][i].c_str());
                      all_data[j][i] = FloatToString(f / 100.0, 0, 2);
                 }
                 _updateRecordBack("EMDATA", j  + 1, all_data[j][i]);
            }
            _updateRecordBack("EMDATA", 19, String::IntToString(i + 1));
       }
}

void Maxit::_pdb_to_ndb_update_TLSGRO()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("TLSGRO");
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            (*lpos)[6] = "refined";
       }
}

void Maxit::_pdb_to_ndb_update_TWIN()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("TWIN");
       if (ppos == _pdb_records.end()) return;

       int count = 1;
       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[4].empty()) (*lpos)[4] = String::IntToString(count);
            (*lpos)[5] = "1";
            (*lpos)[6] = "1";
            count++;
       }
}

void Maxit::_pdb_to_ndb_update_SOLMOD()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("SOLMOD");
       if (ppos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[3].empty()) continue;
            std::string::size_type p = (*lpos)[3].find("(A**2)");
            if (p == std::string::npos) continue;
            (*lpos)[3].erase(p);
            String::StripAndCompressWs((*lpos)[3]);
       }
}
