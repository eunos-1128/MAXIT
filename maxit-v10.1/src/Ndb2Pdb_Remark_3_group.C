/*
FILE:     Ndb2Pdb_Remark_3_group.C
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
#include "RefineProgram.h"
#include "Remark_define.h"
#include "Remark_extern.h"
#include "utillib.h"

void Maxit::_ndb_to_pdb_get_group_remark_3(const std::string& refinement_program)
{
       _ndb_to_pdb_get_remark_3_top();

       int program_index = _pdb_to_ndb_get_refinement_index(refinement_program);

       if (program_index == 0 || program_index == 1 || program_index == 9 || program_index == 11)
            _ndb_to_pdb_get_remark_3_update_Xplor();

       _ndb_to_pdb_update_joint_methods();

       GROUP_REMARKS *group_remark = group_remarks[program_index];
       if (program_index == 5 && _pdb_records.find("NCSLOC") != _pdb_records.end())
            group_remark = &Remark_Refmac_5_0_local;

       if (_joint_methods.size() > 1) {
            for (std::vector<std::string>::const_iterator
                 pos = _joint_methods.begin(); pos != _joint_methods.end(); ++pos) {
                  _current_method = *pos;
                 for (int i = 0; i < NUM_METHOD; ++i) {
                      if (_current_method == _method_name[i]) {
                           std::string cs = " "; cs += _method_remark[i];
                           _addNewRemark(3, cs);
                           break;
                      }
                 }
                 _ndb_to_pdb_get_remark_3(group_remark);
            }
       } else {
            _current_method.clear();
            _ndb_to_pdb_get_remark_3(group_remark);
       }
}

void Maxit::_ndb_to_pdb_get_remark_3_top()
{
       _addNewRemark(3, "REFINEMENT.");

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("REFMET");
       if (ppos == _pdb_records.end()) {
            _addNewRemark(3, "  PROGRAM     : NULL");
            _addNewRemark(3, "  AUTHORS     : NULL");
            return;
       }

       int ContNo = 1;
       bool one_flag = true;
       if (ppos->second.size() > 1) one_flag = false;

       RefineProgram refprogram;
       refprogram.Read(*_logIo, _rcsbroot);

       std::vector<std::string> data;

       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            if ((*lpos)[3].empty() || (*lpos)[3].size() > 62) continue;
            std::string remark = "  PROGRAM     : " + (*lpos)[3];
            if (!one_flag)
                 remark = "  PROGRAM " + FloatToString(ContNo, 3, 0, true)
                        + " : " + (*lpos)[3];
            String::UpperCase(remark);
            _addNewRemark(3, remark);
                        
            std::string program = (*lpos)[3];
            String::UpperCase(program);
            if (program.find("BUSTER-TNT 2") != std::string::npos)
                 program = "BUSTER-TNT 2";
            else if (program.find("BUSTER 2") != std::string::npos)
                 program = "BUSTER 2";
            std::string author = refprogram.findProgramAuthor(program);
            String::UpperCase(author);
            delete_space_between_names(author);

            if (author.empty()) author = "NULL";
            get_max_length_words(data, author, 42);
            
            remark = "  AUTHORS     : " + data[0];
            if (!one_flag)
                 remark = "  AUTHORS " + FloatToString(ContNo, 3, 0, true)
                        + " : " + data[0];
            _addNewRemark(3, remark);
            for (unsigned int i = 1; i < data.size(); ++i) {
                 _addNewRemark(3, "              : " + data[i]);
            }
            ContNo++;
       }
}

void Maxit::_ndb_to_pdb_get_remark_3_update_Xplor()
{
       std::string cs, cs1;
       _getRecordFront("RFACTR", 7, cs);
       _getRecordFront("RFACTR", 10, cs1);
       if (!cs.empty() && cs1.empty()) _updateRecordFront("RFACTR", 10, cs);

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("NCSMOD");
       if (ppos == _pdb_records.end()) _updateRecordFront("NCSMOD", 1, "1");

       _getRecordBack("XFILES", 1, cs);
       int idx = atoi(cs.c_str()) + 1;
       cs = String::IntToString(idx);
       _addNewRecord("XFILES");
       _updateRecordBack("XFILES", 1, cs);

       ppos = _pdb_records.find("XFILES");
       if (ppos != _pdb_records.end()) {
            cs.clear();
            for (std::list<std::vector<std::string> >::const_iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[4].empty()) continue;
                 if (!cs.empty()) cs += ", ";
                 cs += (*lpos)[4];
            }
            if (!cs.empty()) {
                 std::vector<std::string>& Field = ppos->second.front();
                 Field[4] = cs;
            }
       }     
}

void Maxit::_ndb_to_pdb_update_joint_methods()
{
       _joint_methods.clear();

       std::vector<std::string> data;
       std::set<std::string> data_set;
       for (std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
            ppos = _pdb_records.begin(); ppos != _pdb_records.end(); ++ppos) {
            const ndb_token_format& ndbformat = NdbToken::getTokenFormat(ppos->first);
            if (!ndbformat.RefineField) continue;

            data.clear();
            data_set.clear();
            for (std::list<std::vector<std::string> >::const_iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[ndbformat.RefineField - 1].empty()) continue;
                 if (data_set.find((*lpos)[ndbformat.RefineField - 1]) != data_set.end())
                      continue;
                 data_set.insert((*lpos)[ndbformat.RefineField - 1]);
                 data.push_back((*lpos)[ndbformat.RefineField - 1]);
            }
            if (data.size() > _joint_methods.size()) _joint_methods = data;
       }
}

void Maxit::_ndb_to_pdb_get_remark_3(GROUP_REMARKS *group_remark)
{
       for (int i = 0; i < group_remark->Remarks_No; ++i) {
            if (strcmp(group_remark->TokenName[i], "")) {
                 std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
                     ppos = _pdb_records.find(group_remark->TokenName[i]);
                 if (ppos == _pdb_records.end()) continue;

                 for (std::list<std::vector<std::string> >::const_iterator
                      lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                      _ndb_to_pdb_get_general_remark(group_remark->num_remarks[i],
                             group_remark->remarks[i], 0, 
                             atoi((*lpos)[group_remark->repeat[i]-1].c_str()), 0);
                 }
            } else if (group_remark->repeat[i] < 0)
                 _ndb_to_pdb_get_general_remark(group_remark->num_remarks[i],
                            group_remark->remarks[i], 1, 0, 0);
            else _ndb_to_pdb_get_general_remark(group_remark->num_remarks[i],
                            group_remark->remarks[i], 0, 0, 0);
       }
}
