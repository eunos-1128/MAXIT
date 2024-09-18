/*
FILE:     Ndb2Cif.C
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

#include "CategoryMapping.h"
#include "CompositeIndex.h"
#include "Maxit.h"
#include "NdbToken.h"
#include "RmsRecord.h"
#include "SgCenter.h"
#include "utillib.h"

// static void ndb_name_conversion_last_first(std::string &str);

void Maxit::ndb_to_cif()
{
       bool is_refmac5 = false;
       bool is_phenix = false;
       bool is_cns_xplor = false;
       bool is_buster = false;
       _get_refinement_program(is_refmac5, is_phenix, is_cns_xplor, is_buster);

       // Copy R_Value_Obs to R_Value_Work and vs versa for CNS & X-PLOR program
       if (is_cns_xplor) _updateRValueForCNS_XPLOR();

       // Copy RCSBID
       _copyValue("NDBFIL", 1, "RCSBID", 1);

       _ndb_to_cif_update_Ref_ID_and_Align_ID();
       if (is_refmac5) _ndb_to_cif_add_TMLPTL_token_for_REFMAC_5();
       if (is_phenix) _ndb_to_cif_update_TMLPTL_and_NCSGRO_tokens_for_PHENIX();

       // _ndb_to_cif_capture_pdbx_database_remark();
 
       if (!_CifObj) {
            _CifObj = new CifFile(false, Char::eCASE_SENSITIVE, 250, CifString::UnknownValue);
            if (_StructureId.empty()) {
                 _StructureId = _input_filename;
                 std::string::size_type p = _input_filename.find_last_of("/");
                 if (p != std::string::npos) _StructureId = _input_filename.substr(p + 1);
                 p = _StructureId.find_first_of(".");
                 if (p != std::string::npos) _StructureId.erase(p);
                 String::StripAndCompressWs(_StructureId);
                 p = _StructureId.find_first_of(" ");
                 if (p != std::string::npos) _StructureId.erase(p);
                 String::UpperCase(_StructureId);
            }
            _CifObj->AddBlock(_StructureId);
            _firstBlockName = _StructureId;
       }

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       _ndb_to_cif_update_BVALUE_Refine_ID();

       _ndb_to_cif_mapping(_cifblock, (is_refmac5 || is_phenix));

       // cif_update();

       // _get_missing_or_zero_occupancy_residues_or_atoms();
}

void Maxit::_ndb_to_cif_update_Ref_ID_and_Align_ID()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("DBREF");
       if (ppos == _pdb_records.end()) return;

       std::map<std::string, int> dbref_id_mapping;
       dbref_id_mapping.clear();
       int dbref_id = 0;

       for (std::list<std::vector<std::string> >::iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            std::string cs = CompositeIndex::getIndex((*lppos)[7], (*lppos)[9]);
            std::map<std::string, int>::iterator mpos = dbref_id_mapping.find(cs);
            if (mpos != dbref_id_mapping.end())
                 (*lppos)[15] = String::IntToString(mpos->second);
            else {
                 dbref_id++;
                 (*lppos)[15] = String::IntToString(dbref_id);
                 dbref_id_mapping.insert(std::make_pair(cs, dbref_id));
            }
       }

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           qpos = _pdb_records.find("SEQADV");
       if (qpos == _pdb_records.end()) return;

       for (std::list<std::vector<std::string> >::iterator
            lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
            for (std::list<std::vector<std::string> >::const_iterator
                 lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                 if ((*lppos)[2] == (*lqpos)[3]) {
                      (*lqpos)[11] = (*lppos)[14];
                      break;
                 }
            }
       }
}

void Maxit::_ndb_to_cif_add_TMLPTL_token_for_REFMAC_5()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("TMLPTL");
       if (ppos != _pdb_records.end()) return;
       
       std::vector<std::pair<std::string, std::string> > tml_pair;
       tml_pair.clear();

       tml_pair.push_back(std::make_pair("RFTPOS", "tight positional"));
       tml_pair.push_back(std::make_pair("RFMPOS", "medium positional"));
       tml_pair.push_back(std::make_pair("RFLPOS", "loose positional"));
       tml_pair.push_back(std::make_pair("RFTTHR", "tight thermal"));
       tml_pair.push_back(std::make_pair("RFMTHR", "medium thermal"));
       tml_pair.push_back(std::make_pair("RFLTHR", "loose thermal"));
 
       for (std::vector<std::pair<std::string, std::string> >::const_iterator
            pos = tml_pair.begin(); pos != tml_pair.end(); ++pos) {
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
                qpos = _pdb_records.find(pos->first);
            if (qpos == _pdb_records.end()) continue;

            for (std::list<std::vector<std::string> >::const_iterator
                 lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                 _addNewRecord("TMLPTL");
                 ppos = _pdb_records.find("TMLPTL");
                 std::vector<std::string>& field = ppos->second.back();
                 for (unsigned int i = 1; i < lqpos->size(); ++i) field[i] = (*lqpos)[i];
                 field[lqpos->size()] = pos->second;
            }
       }

       ppos = _pdb_records.find("TMLPTL");
       if (ppos == _pdb_records.end()) return;

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           qpos = _pdb_records.find("NCSGRO");

       for (std::list<std::vector<std::string> >::iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            (*lppos)[7] = (*lppos)[1];
            if (qpos == _pdb_records.end()) continue;

            for (std::list<std::vector<std::string> >::const_iterator
                 lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                 if ((*lppos)[7] == (*lqpos)[4] && (*lppos)[2] == (*lqpos)[3]) {
                      (*lppos)[1] = (*lqpos)[1];
                      break;
                 }
            }
       }
}

void Maxit::_ndb_to_cif_update_TMLPTL_and_NCSGRO_tokens_for_PHENIX()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("PHENCS");
       if (ppos == _pdb_records.end()) return;

       _pdb_records.erase("NCSGRO");
       _pdb_records.erase("TMLPTL");
       _pdb_records.erase("CCPNCS");

       std::vector<std::string> data;
       std::string cs;

       for (std::list<std::vector<std::string> >::const_iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            _addNewRecord("NCSGRO");
            _updateRecordBack("NCSGRO", 4, (*lppos)[1]);

            _addNewRecord("TMLPTL");
            _updateRecordBack("TMLPTL", 7, (*lppos)[1]);
            _updateRecordBack("TMLPTL", 6, "POSITIONAL");
            _updateRecordBack("TMLPTL", 3, (*lppos)[5]);

            if ((*lppos)[2] == "1") {
                 _updateRecordBack("NCSGRO", 3, (*lppos)[3]);

                 _addNewRecord("NCSGRO");
                 _updateRecordBack("NCSGRO", 4, (*lppos)[1]);
                 _updateRecordBack("NCSGRO", 3, (*lppos)[4]);

                 String::UpperCase((*lppos)[3], cs);
                 std::string::size_type p = cs.find("CHAIN ");
                 if (p != std::string::npos) {
                      get_wordarray(data, cs.substr(p + 6), " ");
                      if (!data.empty()) _updateRecordBack("TMLPTL", 2, data[0]);
                 }

                 _addNewRecord("TMLPTL");
                 _updateRecordBack("TMLPTL", 7, (*lppos)[1]);
                 _updateRecordBack("TMLPTL", 6, "POSITIONAL");
                 _updateRecordBack("TMLPTL", 3, (*lppos)[5]);
                 _updateRecordBack("TMLPTL", 4, (*lppos)[6]);

                 String::UpperCase((*lppos)[4], cs);
                 p = cs.find("CHAIN ");
                 if (p != std::string::npos) {
                      get_wordarray(data, cs.substr(p + 6), " ");
                      if (!data.empty()) _updateRecordBack("TMLPTL", 2, data[0]);
                 }
            } else {
                 _updateRecordBack("NCSGRO", 3, (*lppos)[4]);

                 _updateRecordBack("TMLPTL", 4, (*lppos)[6]);

                 String::UpperCase((*lppos)[4], cs);
                 std::string::size_type p = cs.find("CHAIN ");
                 if (p != std::string::npos) {
                      get_wordarray(data, cs.substr(p + 6), " ");
                      if (!data.empty()) _updateRecordBack("TMLPTL", 2, data[0]);
                 }
            }
       }

       ppos = _pdb_records.find("NCSGRO");
       if (ppos != _pdb_records.end()) {
            cs = "xxx";
            int i = 1;
            for (std::list<std::vector<std::string> >::iterator
                 lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                 if ((*lppos)[4] != cs) {
                      cs = (*lppos)[4];
                      i = 1;
                 }

                 (*lppos)[1] = String::IntToString(i);
                 i++;

                 if ((*lppos)[3].empty()) continue;

                 bool found = false;
                 std::map<std::string, std::list<std::vector<std::string> > >::iterator
                     qpos = _pdb_records.find("CCPNCS");
                 if (qpos != _pdb_records.end()) {
                      for (std::list<std::vector<std::string> >::iterator
                           lqpos = qpos->second.begin(); lqpos != qpos->second.end(); ++lqpos) {
                           if ((*lppos)[4] == (*lqpos)[1] && (*lppos)[1] == (*lqpos)[14]) {
                                (*lqpos)[15] = (*lppos)[3];
                                (*lqpos)[6] = "1";
                                found = true;
                                break;
                           }
                      }
                 }
                 if (!found) {
                      _addNewRecord("CCPNCS");
                      _updateRecordBack("CCPNCS", 1, (*lppos)[4]);
                      _updateRecordBack("CCPNCS", 6, "1");
                      _updateRecordBack("CCPNCS", 14, (*lppos)[1]);
                      _updateRecordBack("CCPNCS", 15, (*lppos)[3]);
                 }
                 (*lppos)[3].clear();
            }
       }

       ppos = _pdb_records.find("TMLPTL");
       if (ppos != _pdb_records.end()) {
            cs = "xxx";
            int i = 1;
            for (std::list<std::vector<std::string> >::iterator
                 lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                 if ((*lppos)[7] != cs) {
                      cs = (*lppos)[7];
                      i = 1;
                 }
                 (*lppos)[1] = String::IntToString(i);
                 i++;
            }
       }
}

void Maxit::_ndb_to_cif_capture_pdbx_database_remark()
{
       std::vector<int> remark_no;
       remark_no.clear();
       remark_no.push_back(0);
       remark_no.push_back(650);
       remark_no.push_back(700);

       std::string cs;
       for (std::vector<int>::const_iterator
            ipos = remark_no.begin(); ipos != remark_no.end(); ++ipos) {
            std::map<int, std::vector<std::string> >::const_iterator
                mpos = _remarks.find(*ipos);
            if (mpos == _remarks.end()) continue;
            if (mpos->second.empty()) continue;

            int last = mpos->second.size() - 1;
            if (*ipos == 0) {
                 for (unsigned int i = 0; i < mpos->second.size(); ++i) {
                      if (mpos->second[i].find("ORIGINAL DATA REFERENCE") !=
                          std::string::npos) {
                           last = i - 1;
                           break;
                      }
                 }
            }
            for (int i = last; i >= 0; --i) {
                 if (!mpos->second[i].empty()) {
                      last = i;
                      break;
                 } else last--;
            }

            cs.clear();
            for (int i = 0; i <= last; ++i) {
                 if (!cs.empty()) cs += "\n";
                 cs += mpos->second[i];
            }
            if (cs.empty()) continue;

            _addNewRecord("PDBRMK");
            _updateRecordBack("PDBRMK", 1, String::IntToString(*ipos));
            _updateRecordBack("PDBRMK", 3, cs);
       }
}

void Maxit::_ndb_to_cif_update_BVALUE_Refine_ID()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator ppos = _pdb_records.find("BVALUE");
       if (ppos == _pdb_records.end()) return;

       std::string exp_type = _get_diffraction_experiment_method();
       if (exp_type.empty()) return;

       for (std::list<std::vector<std::string> >::iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            if ((*lppos)[10].empty()) (*lppos)[10] = exp_type;
       }
}

void Maxit::_ndb_to_cif_mapping(Block& block, const bool& refmac_phenix_flag)
{
       int xray_or_nmr = 1;
       if (_experiment_type & EXPERIMENT_TYPE_NMR ||
           _experiment_type & EXPERIMENT_TYPE_NMR_SOLID) xray_or_nmr = 2;

       static const std::map<std::string, CIF_CATEGORY>& CifCategory =
                           CategoryMapping::CifCategory();
       for (std::map<std::string, CIF_CATEGORY>::const_iterator cpos = CifCategory.begin(); cpos != CifCategory.end(); ++cpos) {
            if (!cpos->second.transfer) continue;
            if (!(cpos->second.xray_or_nmr & xray_or_nmr)) continue;

            std::string tokenid = "";
            std::string jrnlid = "";
            if (cpos->second.keywords.size() > 1) {
                 for (unsigned int i = 1; i < cpos->second.keywords.size(); ++i) {
                      if (!cpos->second.keywords[i].NdbInfo[0].NdbTokenName.empty()) {
                           tokenid = cpos->second.keywords[i].NdbInfo[0].NdbTokenName;
                           if (tokenid == "JRNL") {
                                jrnlid = cpos->second.keywords[i].NdbInfo[0].JrnlTokenName;
                                if (!jrnlid.empty()) break;
                           } else if (!tokenid.empty()) break;
                      }
                 }
            } else if (!cpos->second.keywords.empty())
                 tokenid = cpos->second.keywords[0].NdbInfo[0].NdbTokenName;

            if (cpos->second.category == "refine_ls_restr") 
                 _ndb_to_cif_process_refine_ls_restr(block);      
/*
            else if (cpos->second.category == "struct_site_gen")
                 ndbcif_process_site_gen(block, cpos->second);
*/
            else if (cpos->second.category == "database_PDB_rev")
                 _ndb_to_cif_process_database_pdb_rev(block);
            else if (cpos->second.category == "database_PDB_rev_record")
                 continue;
            else if (cpos->second.category == "pdbx_database_related")
                 _ndb_to_cif_process_pdbx_database_related(block);
            else if (!tokenid.empty() && NdbToken::IsMatrixToken(tokenid))
                 _ndb_to_cif_process_matrix(block, cpos->second, tokenid);
            else if (cpos->second.category == "audit_author" ||
                     cpos->second.category == "citation_author" ||
                     cpos->second.category == "citation_editor")
                 _ndb_to_cif_process_name(block, cpos->second);
            else if (cpos->second.category == "refine_ls_restr_ncs") {
                 if (refmac_phenix_flag)
                      _ndb_to_cif_process_refmac_refine_ls_restr_ncs(block, cpos->second);
                 else _ndb_to_cif_process_general(block, cpos->second);
	    } else if (tokenid == "PREMRK")
                 _ndb_to_cif_process_remark(block, cpos->second /*, tokenid */);
            else if (cpos->second.general)
                 _ndb_to_cif_process_general_category(block, cpos->second, tokenid);
            else _ndb_to_cif_process_general(block, cpos->second);
       }
}

void Maxit::_ndb_to_cif_process_refine_ls_restr(Block& block)
{
       RmsRecord rmsrecord;
       if (!rmsrecord.Read(*_logIo, _rcsbroot)) return;

       std::string cs;

       std::string abbreviation = "o";
       _getRecordFront("REFMET", 3, cs);
       if (!cs.empty()) abbreviation = rmsrecord.findAbbreviation(cs);

       std::vector<std::string> types = rmsrecord.findTypes(abbreviation);
       if (types.empty()) return;

       std::vector<std::string> data;

       _joint_methods.clear();
       for (std::vector<std::string>::const_iterator
            pos = types.begin(); pos != types.end(); ++pos) {
            std::map<std::string, std::pair<std::string, int> > item_mapping
                = rmsrecord.findItemTokenMapping(*pos);
            if (item_mapping.empty()) continue;

            for (std::map<std::string, std::pair<std::string, int> >::const_iterator
                 mpos = item_mapping.begin(); mpos != item_mapping.end(); ++mpos) {
                 std::map<std::string, std::list<std::vector<std::string> > >
                    ::const_iterator ppos = _pdb_records.find(mpos->second.first);
                 if (ppos == _pdb_records.end()) continue;

                 const ndb_token_format& ndbformat =
                                NdbToken::getTokenFormat(mpos->second.first);
                 int RefineField = ndbformat.RefineField - 1;
                 if (RefineField < 0) continue;

                 data.clear();
                 for (std::list<std::vector<std::string> >::const_iterator lppos =
                      ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                      if (!(*lppos)[RefineField].empty()) {
                           bool found = false;
                           for (std::vector<std::string>::const_iterator
                                dpos = data.begin(); dpos != data.end(); ++dpos) {
                                if (*dpos == (*lppos)[RefineField]) {
                                     found = true;
                                     break;
                                }
                           }
                           if (!found) data.push_back((*lppos)[RefineField]);
                      }
                 }
                 if (data.size() > _joint_methods.size()) _joint_methods = data;
            }
       }
       if (_joint_methods.empty()) _joint_methods.push_back("");

       ISTable *t = _newTablePtr("refine_ls_restr");

       bool has_value = false;

       int irow = 0;
       for (std::vector<std::string>::const_iterator
            jpos = _joint_methods.begin(); jpos != _joint_methods.end(); ++jpos) {
            for (std::vector<std::string>::const_iterator
                 tpos = types.begin(); tpos != types.end(); ++tpos) {
                 std::map<std::string, std::pair<std::string, int> > item_mapping
                     = rmsrecord.findItemTokenMapping(*tpos);
                 if (item_mapping.empty()) continue;

                 t->AddRow();
                 t->UpdateCell(irow, "type", *tpos);

                 for (std::map<std::string, std::pair<std::string, int> >::const_iterator
                      mpos = item_mapping.begin(); mpos != item_mapping.end(); ++mpos) {
                      std::map<std::string, std::list<std::vector<std::string> > >
                         ::const_iterator ppos = _pdb_records.find(mpos->second.first);
                      if (ppos == _pdb_records.end()) continue;

                      const ndb_token_format& ndbformat =
                                     NdbToken::getTokenFormat(mpos->second.first);
                      int RefineField = ndbformat.RefineField - 1;

                      cs.clear();
                      for (std::list<std::vector<std::string> >::const_iterator lppos =
                           ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                           if (RefineField >= 0 && (*lppos)[RefineField] == *jpos) {
                                cs = (*lppos)[mpos->second.second];
                                break;
                           }
                      }
                      t->UpdateCell(irow, mpos->first, cs);
                      if (!cs.empty()) has_value = true;
                 }
                 t->UpdateCell(irow, "pdbx_refine_id", *jpos);
                 irow++;
            }
       }

       if (has_value)
            block.WriteTable(t);
       else {
            delete t;
            deleteTable(block, "refine_ls_restr");
       }
}

void Maxit::_ndb_to_cif_process_database_pdb_rev(Block& block)
{
       std::string date_original;
       _getRecordFront("HEADER", 2, date_original);

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("REVDAT");
       if (ppos == _pdb_records.end()) {
            deleteTable(block, "database_PDB_rev");
            deleteTable(block, "database_PDB_rev_record");

            if (_include_date_original_flag && !date_original.empty()) {
                 ISTable *t = _newTablePtr("database_PDB_rev");
                 t->AddRow();
                 t->UpdateCell(0, "date_original", date_original);
                 block.WriteTable(t);
            }
            return;
       }

       std::vector<std::string> data;

       std::map<int, std::vector<std::string> > pdb_rev;
       pdb_rev.clear();

       for (std::list<std::vector<std::string> >::const_iterator lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            int num = atoi((*lppos)[1].c_str());
            std::map<int, std::vector<std::string> >::iterator pos = pdb_rev.find(num);
            if (pos == pdb_rev.end()) {
                 data.clear();
                 data.push_back((*lppos)[3]);
                 data.push_back((*lppos)[4]);
                 data.push_back((*lppos)[5]);
                 pdb_rev.insert(std::make_pair(num, data));
                 pos = pdb_rev.find(num);
            }
            for (unsigned int i = 6; i < lppos->size(); ++i) {
                 if (!(*lppos)[i].empty()) pos->second.push_back((*lppos)[i]);
            }
       }

       ISTable *t = _newTablePtr("database_PDB_rev");
       ISTable *t1 = _newTablePtr("database_PDB_rev_record");

       int irow = 0;
       int irow1 = 0;
       for (std::map<int, std::vector<std::string> >::const_iterator pos = pdb_rev.begin(); pos != pdb_rev.end(); ++pos) {
            t->AddRow();
            t->UpdateCell(irow, "num", String::IntToString(pos->first));
            t->UpdateCell(irow, "date", pos->second[0]);
            t->UpdateCell(irow, "replaces", pos->second[1]);
            t->UpdateCell(irow, "mod_type", pos->second[2]);
            if (pos->second[2] == "0") t->UpdateCell(irow, "date_original", date_original);
            irow++;

            for (unsigned int i = 3; i < pos->second.size(); ++i) {
                 t1->AddRow();
                 t1->UpdateCell(irow1, "rev_num", String::IntToString(pos->first));
                 t1->UpdateCell(irow1, "type", pos->second[i]);
                 irow1++;
            }
       }

       block.WriteTable(t);

       if (irow1 > 0) block.WriteTable(t1);
       else {
            delete t1;
            deleteTable(block, "database_PDB_rev_record");
       }
}

void Maxit::_ndb_to_cif_process_pdbx_database_related(Block& block)
{
       ISTable *t = _newTablePtr("pdbx_database_related");
       int row = 0;

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("RENTRY");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::const_iterator
                 lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                 t->AddRow();
                 t->UpdateCell(row, "db_name", (*lppos)[2]);
                 t->UpdateCell(row, "db_id", (*lppos)[3]);
                 t->UpdateCell(row, "content_type", "unspecified");
                 t->UpdateCell(row, "details", (*lppos)[4]);
                 row++;
            }
       }

       ppos = _pdb_records.find("SPLIT");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::const_iterator
                 lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                 for (unsigned int i = 2; i < lppos->size(); ++i) {
                      if ((*lppos)[i].empty()) continue;
                      t->AddRow();
                      t->UpdateCell(row, "db_name", "PDB");
                      t->UpdateCell(row, "db_id", (*lppos)[i]);
                      t->UpdateCell(row, "content_type", "split");
                      row++;
                 }
            }
       }

       if (!_original_entry_ids.empty()) {
            std::string db_name, db_id, content_type;
            for (std::vector<std::string>::iterator pos =
                 _original_entry_ids.begin(); pos != _original_entry_ids.end(); ++pos) {
                 bool found = false;
                 for (int i = 0; i < row; ++i) {
                      get_value_clean(db_name, t, i, "db_name");
                      get_value_clean(db_id, t, i, "db_id");
                      get_value_clean(content_type, t, i, "content_type");
                      if (db_name == "PDB" && db_id == *pos) {
                           if (content_type.empty())
                                t->UpdateCell(i, "content_type", "re-refinement");
                           found = true;
                           break;
                      }
                 }
                 if (found) continue;
                
                 t->AddRow();
                 t->UpdateCell(row, "db_name", "PDB");
                 t->UpdateCell(row, "db_id", *pos);
                 t->UpdateCell(row, "content_type", "re-refinement");
                 row++;
            }
       }

       if (row > 0)
            block.WriteTable(t);
       else {
            delete t;
            deleteTable(block, "pdbx_database_related");
       }
}

void Maxit::_ndb_to_cif_process_matrix(Block& block, const CIF_CATEGORY &Category,
                                       const std::string& tokenid)
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(tokenid);
       if (ppos == _pdb_records.end()) {
            deleteTable(block, Category.category);
            return;
       }

       ISTable *t = _newTablePtr(Category);
       int irow = -1;

       for (std::list<std::vector<std::string> >::const_iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            std::string TokenName = (*lppos)[0];
            std::string TokenIndex = TokenName.substr(TokenName.size() - 1);
            TokenName.erase(TokenName.size() - 1);
            if (TokenIndex == "1") {
                 t->AddRow();
                 irow++;
            }

            for (std::vector<CIF_KEYWORD>::const_iterator kpos =
                 Category.keywords.begin(); kpos != Category.keywords.end(); ++kpos) {
                 if (kpos->name == "entry_id")
                      t->UpdateCell(irow, kpos->name, _StructureId);
                 else if (kpos->NdbInfo[0].NdbTokenName.empty() ||
                         !kpos->NdbInfo[0].NdbFieldNo) continue;
                 else if (kpos->NdbInfo[0].NdbTokenName == tokenid) {
                      int fieldno = kpos->NdbInfo[0].NdbFieldNo - 1;
                      std::string ItemIndex = "";
                      std::string::size_type idx = kpos->name.find_first_of('[');
                      if (idx != std::string::npos)
                           ItemIndex = kpos->name.substr(idx + 1, 1);
                      if (!ItemIndex.empty()) {
                           if (ItemIndex == TokenIndex)
                                t->UpdateCell(irow, kpos->name, (*lppos)[fieldno]);
                      } else {
                           std::string cs = (*lppos)[fieldno];
                           if (kpos->name == "code") {
                                if (cs.empty())     cs = "generate";
                                else if (cs == "1") cs = "given";
                           }
                           t->UpdateCell(irow, kpos->name, cs);
                      }
                 }
            }
       }
       block.WriteTable(t);
}

void Maxit::_ndb_to_cif_process_name(Block& block, const CIF_CATEGORY &Category)
{
       std::string tokenid = Category.keywords[0].NdbInfo[0].NdbTokenName;
       if (tokenid == "JRNL") tokenid = Category.keywords[0].NdbInfo[0].JrnlTokenName;

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(tokenid);
       if (ppos == _pdb_records.end()) {
            deleteTable(block, Category.category);
            return;
       }

       std::vector<std::string> data;

       ISTable *t = _newTablePtr(Category);


       std::string ordinal = "";
       int id_fieldno = -1;
       for (std::vector<CIF_KEYWORD>::const_iterator
            kpos = Category.keywords.begin(); kpos != Category.keywords.end(); ++kpos) {
            if (kpos->name.find("ordinal") != std::string::npos) {
                 ordinal = kpos->name;
            }
            if (kpos->name == "citation_id") {
                 for (std::vector<ndb_token_match>::const_iterator
                      ipos = kpos->NdbInfo.begin(); ipos != kpos->NdbInfo.end(); ++ipos) { 
                      id_fieldno = ipos->NdbFieldNo - 1;
                 }
            }
       }
      
       int irow = 0;
       for (std::list<std::vector<std::string> >::const_iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            std::string citation_id = "";
            if (id_fieldno >= 0) {
                 citation_id = (*lppos)[id_fieldno];
                 if (citation_id.empty() || citation_id == "0") citation_id = "primary";
            }
            for (std::vector<CIF_KEYWORD>::const_iterator kpos =
                 Category.keywords.begin(); kpos != Category.keywords.end(); ++kpos) {
                 if (kpos->name == "citation_id") continue;
                 for (std::vector<ndb_token_match>::const_iterator
                      ipos = kpos->NdbInfo.begin(); ipos != kpos->NdbInfo.end(); ++ipos) { 
                      if (!ipos->NdbFieldNo) continue;
                      int fieldno = ipos->NdbFieldNo - 1;
                      if ((*lppos)[fieldno].empty()) continue;

                      get_wordarray(data, (*lppos)[fieldno], ",");
                      for (std::vector<std::string>::const_iterator
                           vpos = data.begin(); vpos != data.end(); ++vpos) {
                           std::string cs = *vpos;
                           String::StripAndCompressWs(cs);
                           ndb_name_conversion_last_first(cs);
                           t->AddRow();
                           t->UpdateCell(irow, kpos->name, cs);
                           if (!citation_id.empty())
                                t->UpdateCell(irow, "citation_id", citation_id);
                           if (!ordinal.empty())
                                t->UpdateCell(irow, ordinal, String::IntToString(irow+1));
                           irow++;
                      }
                 }
            }
       }

       if (irow > 0)
            block.WriteTable(t);
       else {
            delete t;
            deleteTable(block, Category.category);
       }
}

void Maxit::_ndb_to_cif_process_refmac_refine_ls_restr_ncs(Block& block,
                                            const CIF_CATEGORY &Category)
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("TMLPTL");
       if (ppos == _pdb_records.end()) {
            deleteTable(block, Category.category);
            return;
       }

       ISTable *t = _newTablePtr(Category);

       int n_keywords = 7;
       const char *keywords[7] = { "dom_id", "pdbx_auth_asym_id", "pdbx_number",
                       "rms_dev_position", "weight_position", "pdbx_type", "pdbx_ens_id" };

       int irow = 0;
       for (std::list<std::vector<std::string> >::const_iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            t->AddRow();
            for (int i = 0; i < n_keywords; ++i) {
                 t->UpdateCell(irow, keywords[i], (*lppos)[i + 1]);
            }
            irow++;
       }
       block.WriteTable(t);
}

void Maxit::_ndb_to_cif_process_general(Block& block, const CIF_CATEGORY &Category)
{
       std::vector<std::string> Value;
       std::string cs, cs1;

       ISTable *t = _newTablePtr(Category);
       int row = 0;

       for (std::vector<CIF_KEYWORD>::const_iterator
            kpos = Category.keywords.begin(); kpos != Category.keywords.end(); ++kpos) {
            for (std::vector<ndb_token_match>::const_iterator
                 ipos = kpos->NdbInfo.begin(); ipos != kpos->NdbInfo.end(); ++ipos) { 
                 std::string tokenid = ipos->NdbTokenName;
                 int fieldno = ipos->NdbFieldNo - 1;
                 bool IsJrnl = false;

                 if (tokenid == "JRNL") { 
                      tokenid = ipos->JrnlTokenName;
                      IsJrnl = true;
                 }
                 if (tokenid.empty() || fieldno < 0) continue;

                 std::map<std::string, std::list<std::vector<std::string> > >
                    ::const_iterator ppos = _pdb_records.find(tokenid);
                 if (ppos == _pdb_records.end()) continue;

                 int s_fieldno = -1;
                 const ndb_token_format& ndbformat = NdbToken::getTokenFormat(tokenid);
                 if (ndbformat.SeqField > 0) s_fieldno = ndbformat.SeqField - 1;

                 int SerialNo = 0;
                 for (std::list<std::vector<std::string> >::const_iterator
                      lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
                      if (Category.category == "struct_ncs_ens_gen" &&
                         (*lppos)[9].empty()) continue;

                      cs.clear();

                      bool logical = false;
                      for (std::vector<int>::const_iterator
                           pos = ipos->NdbTokenKeyFieldNo.begin();
                           pos != ipos->NdbTokenKeyFieldNo.end(); ++pos) {
                           if (fieldno == *pos) logical = true;
                      }
                      if (IsJrnl) {
                           if (fieldno == 1 &&
                              ((*lppos)[fieldno].empty() || (*lppos)[fieldno] == "0"))
                                cs = "primary";
                           else cs = (*lppos)[fieldno];
                      } else if (fieldno == 0) {
                           if (!ipos->JrnlTokenName.empty())
                                cs = ipos->JrnlTokenName;
                           else cs = tokenid;
                      } else if (logical) {
                           if (!ipos->JrnlTokenName.empty()) {
                                if (!(*lppos)[fieldno].empty())
                                     cs = ipos->JrnlTokenName + (*lppos)[fieldno];
                                else cs = ipos->JrnlTokenName + "1";
                           } else {
                                if (!(*lppos)[fieldno].empty())
                                     cs = (*lppos)[fieldno];
                                else cs = "1";
                           }
                      } else {
                           if (tokenid == "CRSOLU" && fieldno == 1) {
                                int comp_id = 0;
                                int crystal_id = atoi((*lppos)[2].c_str());
                                for (std::list<std::vector<std::string> >::const_iterator
                                     lqpos = ppos->second.begin();
                                     lqpos != ppos->second.end(); ++lqpos) {
                                     if (atoi((*lqpos)[2].c_str()) == crystal_id)
                                          comp_id++;
                                     if ((*lppos)[1] == (*lqpos)[1]) break;
                                }
                                cs = String::IntToString(comp_id);
                           } else cs = (*lppos)[fieldno];
                      }

                      if (!Category.CifKeyField.empty()) {
                           Value.clear();
                           for (std::vector<int>::const_iterator
                                pos = ipos->NdbTokenKeyFieldNo.begin();
                                pos != ipos->NdbTokenKeyFieldNo.end(); ++pos) {
                                if (*pos > 0) {
                                     if (IsJrnl) {
                                          if (!(*lppos)[*pos].empty() &&
                                               (*lppos)[*pos] != "0")
                                               Value.push_back((*lppos)[*pos]);
                                          else Value.push_back("primary");
                                     } else if (!ipos->JrnlTokenName.empty()) {
                                          if (!(*lppos)[*pos].empty())
                                               cs1 = ipos->JrnlTokenName + (*lppos)[*pos];
                                          else cs1 = ipos->JrnlTokenName + "1";
                                          Value.push_back(cs1);                              
                                     } else if (!(*lppos)[*pos].empty())
                                          Value.push_back((*lppos)[*pos]);
                                     else Value.push_back("1");
                                } else if (*pos == 0)
                                     Value.push_back(ipos->JrnlTokenName);
                           }

                           bool found = false;
                           unsigned int size = Category.CifKeyField.size();
                           if (Value.size() < size) size = Value.size();
                           for (int serialNo1 = 0; serialNo1 < row; ++serialNo1) {
                                bool found1 = true;
                                for (unsigned int k = 0; k < size; ++k) {
                                     int keyField = Category.CifKeyField[k];
                                     get_value_clean(cs1, t, serialNo1,
                                               Category.keywords[keyField].name);
                                     if (cs1 != Value[k] && !Value[k].empty()) {
                                          found1 = false;
                                          break;
                                     }
                                }
                                if (found1) {
                                     found = true;
                                     t->UpdateCell(serialNo1, kpos->name, cs);
                                }
                           }
                           if (!found) {
                                t->AddRow();
                                t->UpdateCell(row, kpos->name, cs);
                                for (unsigned int k = 0; k < size; ++k) {
                                     int keyField = Category.CifKeyField[k];
                                     t->UpdateCell(row,
                                            Category.keywords[keyField].name, Value[k]);
                                }
                                row++;
                           }
                      } else if (cs != "") {
                           if (s_fieldno >= 0) {
                                if (IsJrnl) 
                                     SerialNo = atoi((*lppos)[s_fieldno - 1].c_str());
                                else if (Category.category == "audit_contact_author")
                                     SerialNo = atoi((*lppos)[s_fieldno - 1].c_str()) - 1;
                           }
                           if (SerialNo >= row) {
                                row++;
                                t->AddRow();
                           }
                           t->UpdateCell(SerialNo, kpos->name, cs);
                      }
                      SerialNo++;
                 }
            }
       }

       if (row > 0)
            block.WriteTable(t);
       else {
            delete t;
            deleteTable(block, Category.category);
       }
}

void Maxit::_ndb_to_cif_process_remark(Block& block, const CIF_CATEGORY &Category)
{
/*
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(tokenid);
       if (ppos == _pdb_records.end()) {
            deleteTable(block, Category.category);
            return;
       }

       std::string cs, id;
       cs.clear();
       id = "------1";

       std::vector<std::pair<std::string, std::string> > remarks;
       remarks.clear();

       for (std::list<std::vector<std::string> >::const_iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            if ((*lppos)[1] != id) {
                 if (!cs.empty()) remarks.push_back(std::make_pair(id, cs));
                 id = (*lppos)[1];
                 cs = (*lppos)[2];
            } else cs += "\n" + (*lppos)[2];
       }
       if (!cs.empty()) remarks.push_back(std::make_pair(id, cs));
*/
       if (_remarks.empty()) {
            deleteTable(block, Category.category);
            return;
       }

       std::string cs;
       ISTable *t = _newTablePtr(Category);
       int row = 0;
       // for (std::vector<std::pair<std::string, std::string> >::const_iterator
       for (std::map<int, std::vector<std::string> >::const_iterator
            mpos = _remarks.begin(); mpos != _remarks.end(); ++mpos) {
            cs.clear();
            for (std::vector<std::string>::const_iterator
                 pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                 if (!cs.empty()) cs += "\n";
                 cs += *pos;
            }
            t->AddRow();
            t->UpdateCell(row, "id", String::IntToString(mpos->first));
            t->UpdateCell(row, "text", cs);
            row++;
       }
       block.WriteTable(t);
}

void Maxit::_ndb_to_cif_process_general_category(Block& block, const CIF_CATEGORY&
                                                 Category, const std::string& tokenid)
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find(tokenid);
       if (ppos == _pdb_records.end()) {
            deleteTable(block, Category.category);
            return;
       }

       ISTable *t = _newTablePtr(Category);

       int row = 0;
       for (std::list<std::vector<std::string> >::const_iterator
            lppos = ppos->second.begin(); lppos != ppos->second.end(); ++lppos) {
            t->AddRow();
            for (std::vector<CIF_KEYWORD>::const_iterator kpos =
                 Category.keywords.begin(); kpos != Category.keywords.end(); ++kpos) {
                 if (kpos->name == "citation_id") continue;
                 for (std::vector<ndb_token_match>::const_iterator
                      ipos = kpos->NdbInfo.begin(); ipos != kpos->NdbInfo.end(); ++ipos) { 
                      int fieldno = ipos->NdbFieldNo - 1;
                      if (fieldno < 0) continue;
                      t->UpdateCell(row, kpos->name, (*lppos)[fieldno]);
                 }
            }
            row++;
       }
       block.WriteTable(t);
}
/*
static void ndb_name_conversion_last_first(std::string &str)
{
       if (str.empty()) return;

       if (SgCenter::IsSgCenter(str)) return;

       std::string::size_type pos = str.find_last_of(".");
       if (pos == str.size() - 1)
            pos = str.substr(0, str.size() - 1).find_last_of(".");
       if (pos > 2 && str[pos-2] == 'S' &&
          (str[pos-1] == 'T' || str[pos-1] == 't'))
            pos = str.substr(0, pos).find_last_of(".");
       if (pos != std::string::npos) {
            std::string last = str.substr(pos+1);
            std::string first = str.substr(0, pos+1);
            str = last + ", " + first;
            String::StripAndCompressWs(str);
            return;
       }

       String::StripAndCompressWs(str);
       pos = str.find_last_of(' ');
       if (pos != std::string::npos) {
            std::string last = str.substr(pos+1);
            std::string first = str.substr(0, pos+1);
            str = last + ", " + first;
            String::StripAndCompressWs(str);
       }
}
*/
