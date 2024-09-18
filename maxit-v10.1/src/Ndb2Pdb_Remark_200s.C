/*
FILE:     Ndb2Pdb_Remark_200s.C
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

/*
void Maxit::_ndb_to_pdb_get_remark_100(const bool& is_unprocessed_model, const std::string& EBI_ID)
{
       std::string Process_Site;
       _getRecordFront("PROCES", 5, Process_Site);

       if (is_unprocessed_model)
            _ndb_to_pdb_get_general_remark(Num_Remark_100_model, Remark_100_model, 0, 0, 0);
       else if (!EBI_ID.empty()) {
            _updateRecordFront("NDBFIL", 1, EBI_ID);
            _updateRecordFront("PROCES", 5, "PDBE");
            _ndb_to_pdb_get_general_remark(Num_Remark_100_EBI, Remark_100_EBI, 0, 0, 0); 
       } else if (!Process_Site.empty()) {
            std::string wwpdb_id, rcsb_id, pdb_id, bmrb_id;
            _getRecordFront("WPDBID", 1, wwpdb_id);
            _getRecordFront("RCSBID", 1, rcsb_id);
            if (rcsb_id.find("RCSB") != std::string::npos)
                 _updateRecordFront("NDBFIL", 1, rcsb_id); 

            _getRecordFront("NDBFIL", 1, rcsb_id);
            _getRecordFront("PDBFIL", 1, pdb_id);
            _getRecordFront("BMRBID", 1, bmrb_id);
            // if (rcsb_id.find("RCSB") != std::string::npos)
            if (!wwpdb_id.empty())
                 _ndb_to_pdb_get_general_remark(Num_Remark_100, Remark_100, 0, 0, 0);
            // else if (Process_Site == "RCSB")
            else if (rcsb_id.substr(0, 4) == "RCSB")
                 _ndb_to_pdb_get_general_remark(Num_Remark_100_RCSB, Remark_100_RCSB, 0, 0, 0);
            else if (pdb_id.empty() && !bmrb_id.empty())
                 _ndb_to_pdb_get_general_remark(Num_Remark_100_BMRB, Remark_100_BMRB, 0, 0, 0);
            else _ndb_to_pdb_get_general_remark(1, &Remark_100_BNL, 0, 0, 0);
       } else _ndb_to_pdb_get_general_remark(1, &Remark_100_BNL, 0, 0, 0); 
}
*/

void Maxit::_ndb_to_pdb_get_remark_100(const bool& is_unprocessed_model, const std::string& EBI_ID)
{
       std::string Process_Site, Process_Date, Dep_ID;
       _getRecordFront("PROCES", 5, Process_Site);
       _getRecordFront("NDBPRO", 4, Process_Date);
       _getRecordFront("WPDBID", 1, Dep_ID);

       std::string rcsb_id, pdb_id, bmrb_id;
       _getRecordFront("RCSBID", 1, rcsb_id);
       if (rcsb_id.find("RCSB") != std::string::npos) _updateRecordFront("NDBFIL", 1, rcsb_id); 

       _getRecordFront("NDBFIL", 1, rcsb_id);
       _getRecordFront("PDBFIL", 1, pdb_id);
       _getRecordFront("BMRBID", 1, bmrb_id);
       if (Dep_ID.empty() && rcsb_id.empty() && pdb_id.empty() && bmrb_id.empty()) return;

       if (!Dep_ID.empty()) {
            if (!EBI_ID.empty()) _updateRecordFront("PROCES", 5, "PDBE");
            else if (Process_Site.empty()) _updateRecordFront("PROCES", 5, "PDB");
            if (!Process_Date.empty())
                 _ndb_to_pdb_get_general_remark(Num_Remark_100, Remark_100, 0, 0, 0);
            else _ndb_to_pdb_get_general_remark(Num_Remark_100_NoDate, Remark_100_NoDate, 0, 0, 0);
       } else if (is_unprocessed_model)
            _ndb_to_pdb_get_general_remark(Num_Remark_100_model, Remark_100_model, 0, 0, 0);
       else if (!EBI_ID.empty()) {
            _updateRecordFront("NDBFIL", 1, EBI_ID);
            _updateRecordFront("PROCES", 5, "PDBE");
            _ndb_to_pdb_get_general_remark(Num_Remark_100_EBI, Remark_100_EBI, 0, 0, 0); 
       } else if (!Process_Site.empty()) {
            if (rcsb_id.substr(0, 4) == "RCSB")
                 _ndb_to_pdb_get_general_remark(Num_Remark_100_RCSB, Remark_100_RCSB, 0, 0, 0);
            else if (pdb_id.empty() && !bmrb_id.empty())
                 _ndb_to_pdb_get_general_remark(Num_Remark_100_BMRB, Remark_100_BMRB, 0, 0, 0);
            else _ndb_to_pdb_get_general_remark(1, &Remark_100_BNL, 0, 0, 0);
       } else _ndb_to_pdb_get_general_remark(1, &Remark_100_BNL, 0, 0, 0); 
}

void Maxit::_ndb_to_pdb_get_remark_265()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("SOLEXP");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if (String::IsEqual((*lpos)[2], "modelling", Char::eCASE_INSENSITIVE) ||
                     String::IsEqual((*lpos)[2], "THEORETICAL MODELLING",
                                      Char::eCASE_INSENSITIVE)) {
                      lpos = ppos->second.erase(lpos);
                      --lpos;
                 } else if (String::IsEqual((*lpos)[2], "x-ray", Char::eCASE_INSENSITIVE))
                      (*lpos)[2] = "SMALL ANGLE X-RAY SCATTERING";
                 else if (String::IsEqual((*lpos)[2], "neutron", Char::eCASE_INSENSITIVE))
                      (*lpos)[2] = "SMALL ANGLE NEUTRON SCATTERING";
            }
            if (ppos->second.empty()) _pdb_records.erase("SOLEXP");
       }

       for (int i = 0; i < Remark_265.Remarks_No; ++i) {
            if (strcmp(Remark_265.TokenName[i], "")) {
                 ppos = _pdb_records.find(Remark_265.TokenName[i]);
                 if (ppos == _pdb_records.end()) continue;
            }
            int j = 0;
            if (Remark_265.repeat[i]) j = 1;
            _ndb_to_pdb_get_general_remark(Remark_265.num_remarks[i], Remark_265.remarks[i],
                                           j, 0, 1);
       }
}

void Maxit::_ndb_to_pdb_get_remark_280() 
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("CRDTLS");
       if (ppos == _pdb_records.end()) {
            ppos = _pdb_records.find("CRSOLU");
            if (ppos != _pdb_records.end()) {
                 std::string cs = "";
                 for (std::list<std::vector<std::string> >::const_iterator
                      lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                      if ((*lpos)[3] == "1" && !String::IsEqual((*lpos)[4], "WATER",
                           Char::eCASE_INSENSITIVE)) {
                           if (!cs.empty()) cs += ", ";
                           cs += (*lpos)[4];
                      }
                 }
                 if (!cs.empty()) _updateRecordFront("CRDTLS", 2, cs);
            }
       }

       ppos = _pdb_records.find("CRDTLS");
       if (ppos != _pdb_records.end()) {
            std::string cs = "";
            for (std::list<std::vector<std::string> >::const_iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[2].empty()) continue;
                 if (!cs.empty()) {
                      if (cs[cs.size() - 1] == '.')
                           cs += " ";
                      else cs += ". ";
                 }
                 cs += (*lpos)[2];
            }
            if (!cs.empty()) {
                 std::vector<std::string>& Field = ppos->second.front();
                 Field[2] = cs;
            }
       }
       _ndb_to_pdb_get_general_remark(Num_Remark_280, Remark_280, 0, 0, 0);
}

void Maxit::_ndb_to_pdb_get_remark_290()
{
       if (_cell.is_artifical()) return;

       _ndb_to_pdb_update_template_FieldWidth(Num_Remark_290, Remark_290);

       double result[3][3], tr[3];
       std::string remark_string, cs;
       for (int j = 0; j < Num_Remark_290; ++j) {
            int lineNo = 1;
            if (Remark_290[j].Special_Treatment == PDB_REMARK_NNNMMM_SYMOP)
                 lineNo = _cell.symops().size();
            else if (Remark_290[j].Special_Treatment == PDB_REMARK_XYZT) 
                 lineNo = _cell.symops().size() * 3;
      
            for (int k = 0; k < lineNo; ++k) {
                 remark_string.clear();

                 if (Remark_290[j].Special_Treatment == PDB_REMARK_XYZT && k % 3 == 0)
                      _cell.ndb_get_symmetry_matrix(k/3, 0, 0, 0, result, tr);

                 for (int i = 0; i < Remark_290[j].NumField; ++i) {
                      bool left_adjust = false;
                      if (Remark_290[j].FieldList[i].FieldJustification ==
                          PDB_REMARK_LEFT_JUSTIFIED) left_adjust = true;

                      remark_string += get_space_between_remark_field(Remark_290[j], i);
                      if (strcmp(Remark_290[j].FieldList[i].TokenName,"") ||
                                 Remark_290[j].FieldList[i].FieldId != 0) {
                           std::string card_id = Remark_290[j].FieldList[i].TokenName;
                           int field_no = Remark_290[j].FieldList[i].FieldId-1;

                           if (!strcmp(Remark_290[j].FieldList[i].Text,"NNNMMM")) {
                                cs = String::IntToString(k + 1) + "555";
                                remark_string += FormattedFieldValue(cs,
                                                Remark_290[j].FieldList[i].FieldType,
                                                Remark_290[j].FieldList[i].FieldWidth,
                                                Remark_290[j].FieldList[i].FieldPrec,
                                                left_adjust);
                           } else if (!strcmp(Remark_290[j].FieldList[i].Text, "SYMOP")) {
                                _cell.ndb_get_symmetry_operation_name(k, 0, 0, 0, cs);
                                String::UpperCase(cs);
                                remark_string += FormattedFieldValue(cs,
                                                Remark_290[j].FieldList[i].FieldType,
                                                Remark_290[j].FieldList[i].FieldWidth,
                                                Remark_290[j].FieldList[i].FieldPrec,
                                                left_adjust);
                           } else if (!strcmp(Remark_290[j].FieldList[i].Text, "SMTRY")) {
                                cs = Remark_290[j].FieldList[i].Text
                                   + String::IntToString( k % 3 + 1);;
                                remark_string += FormattedFieldValue(cs,
                                                Remark_290[j].FieldList[i].FieldType,
                                                Remark_290[j].FieldList[i].FieldWidth,
                                                Remark_290[j].FieldList[i].FieldPrec,
                                                left_adjust);
                           } else if (!strcmp(Remark_290[j].FieldList[i].Text, "Serial_No")) {
                                remark_string += FloatToString((double) (k / 3 + 1),
                                                  Remark_290[j].FieldList[i].FieldWidth,
                                                  Remark_290[j].FieldList[i].FieldPrec,
                                                  left_adjust);
                           } else if (!strcmp(Remark_290[j].FieldList[i].Text, "X")) { 
                                remark_string += FloatToString(result[k % 3][0],
                                                  Remark_290[j].FieldList[i].FieldWidth,
                                                  Remark_290[j].FieldList[i].FieldPrec,
                                                  left_adjust);
                           } else if (!strcmp(Remark_290[j].FieldList[i].Text, "Y")) {
                                remark_string += FloatToString(result[k % 3][1],
                                                  Remark_290[j].FieldList[i].FieldWidth,
                                                  Remark_290[j].FieldList[i].FieldPrec,
                                                  left_adjust);
                           } else if (!strcmp(Remark_290[j].FieldList[i].Text, "Z")) {
                                remark_string += FloatToString(result[k % 3][2],
                                                  Remark_290[j].FieldList[i].FieldWidth,
                                                  Remark_290[j].FieldList[i].FieldPrec,
                                                  left_adjust);
                           } else if (!strcmp(Remark_290[j].FieldList[i].Text, "T")) {
                                remark_string += FloatToString(tr[k % 3],
                                                  Remark_290[j].FieldList[i].FieldWidth,
                                                  Remark_290[j].FieldList[i].FieldPrec,
                                                  left_adjust);
                           } else if (!card_id.empty() && field_no > 0) {
                                cs.clear();
                                if (card_id == "CRYST1") {
                                     _getRecordFront("CRYST1", field_no, cs);
                                     if (cs.empty()) _getRecordFront("CRYST1", 7, cs);
                                }
                                if (!cs.empty())
                                     remark_string += FormattedFieldValue(cs,
                                                     Remark_290[j].FieldList[i].FieldType,
                                                     Remark_290[j].FieldList[i].FieldWidth,
                                                     Remark_290[j].FieldList[i].FieldPrec,
                                                     left_adjust);
                                else remark_string += FormattedFieldValue("NULL", 3,
                                   Remark_290[j].FieldList[i].FieldWidth, 0, left_adjust);
                           } else remark_string += FormattedFieldValue("NULL", 3,
                                   Remark_290[j].FieldList[i].FieldWidth, 0, left_adjust);
                      } else {
                           remark_string += FormattedFieldValue(
                                                Remark_290[j].FieldList[i].Text,
                                                Remark_290[j].FieldList[i].FieldType,
                                                Remark_290[j].FieldList[i].FieldWidth,
                                                Remark_290[j].FieldList[i].FieldPrec,
                                                left_adjust);
                      }
                 }
                 _addNewRemark(Remark_290[j].Remark_No, remark_string);
            }
       }
}
