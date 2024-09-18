/*
FILE:     Ndb2Pdb_Remark.C
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

void Maxit::_ndb_to_pdb_get_remarks()
{
       _get_missing_or_zero_occupancy_residues_or_atoms_for_remarks();

       std::vector<int> atom_counts;
       _check_atom_number_count(atom_counts);
       if (!atom_counts.empty()) {
            for (unsigned int i = 0; i < atom_counts.size(); ++i) {
                 _updateRecordFront("ATNUMS", i + 2, String::IntToString(atom_counts[i]));
            }
       }

       _pdb_records.erase("REMARK");

       bool is_refmac5 = false;
       bool is_phenix = false;
       bool is_cns_xplor = false;
       bool is_buster = false;
       std::string refinement_program = _get_refinement_program(is_refmac5, is_phenix, is_cns_xplor, is_buster); 

       int num_twin = 0;
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("TWIN");
       if (ppos != _pdb_records.end()) num_twin = ppos->second.size();

       int num_tls = _ndb_to_pdb_update_TLSGRO();

       int num_cns = 0;
       if (is_refmac5) num_cns = _ndb_to_pdb_processing_refmac5();
       if (is_phenix || is_buster) num_cns = _ndb_to_pdb_processing_phenix(is_buster);

       if (num_cns)  _updateRecordFront("NCSTLS", 1, String::IntToString(num_cns), true);
       if (num_tls)  _updateRecordFront("NCSTLS", 2, String::IntToString(num_tls), true);
       if (num_twin) _updateRecordFront("NCSTLS", 3, String::IntToString(num_twin), true);

       _ndb_to_pdb_update_SCTYPE();
       _ndb_to_pdb_processing_shelx_and_tnt();
       _ndb_to_pdb_processing_synchrotron();

       ppos = _pdb_records.find("PROCES");
       if (ppos != _pdb_records.end()) {
            std::vector<std::string>& Field = ppos->second.front();
            if (String::IsEqual(Field[5], "OSAKA", Char::eCASE_INSENSITIVE))
                 Field[5] = "PDBJ";
            if (String::IsEqual(Field[5], "PRAGUE", Char::eCASE_INSENSITIVE))
                 Field[5] = "RCSB";
       }

       // re-scale _em_image_recording.avg_electron_dose_per_image: multiply by 100 for CIF -> PDB
       ppos = _pdb_records.find("EMDATA");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[12].empty()) continue;
                 double f = atof((*lpos)[12].c_str());
                 (*lpos)[12] = FloatToString(f * 100.0, 0, 2);
            }
       }

       std::set<int> SkipRemarkNo;
       bool is_unprocessed_model = _ndb_to_pdb_get_skip_remark_number(SkipRemarkNo);
       if (is_unprocessed_model) {}

       std::string cs;

       bool is_electron_microscopy = false;
       bool is_powder_diffration = false;
       _getRecordFront("EXPDTA", 2, cs);
       if (!cs.empty()) {
            String::UpperCase(cs);
            if (cs.find("ELECTRON MICROSCOPY") != std::string::npos || cs.find("CRYO-EM RECONSTRUCTION") != std::string::npos)
                 is_electron_microscopy = true;
            if (cs.find("POWDER DIFFRACTION") != std::string::npos)
                 is_powder_diffration = true;
       }

       if (_experiment_type & EXPERIMENT_TYPE_NMR || _experiment_type & EXPERIMENT_TYPE_NMR_SOLID || _experiment_type & EXPERIMENT_TYPE_MODEL) {
            _ndb_to_pdb_processing_nmr_software();
            _ndb_to_pdb_processing_nmr_experiment();
       }

       _ndb_to_pdb_processing_em_software();
       _ndb_to_pdb_processing_em_pixel();

       std::string EBI_ID;
       bool is_version_v315 = false;
       _ndb_to_pdb_get_auxiliary_data(EBI_ID, is_version_v315);

       std::map<int, std::vector<std::string> > user_remarks;
       _get_remark_map("PDBRMK", 3, user_remarks);

       std::map<int, std::vector<std::string> >::const_iterator user_pos = user_remarks.find(0);
       if (user_pos != user_remarks.end()) _get_user_defined_remark(0, user_pos->second);
          
       for (int remark_no = 2; remark_no <= 999; ++remark_no) {
            if (SkipRemarkNo.find(remark_no) != SkipRemarkNo.end()) continue;

            if (_pdb_header_only_flag && (remark_no > 350) && (remark_no != 400)) continue;

            user_pos = user_remarks.find(remark_no);
            if (user_pos != user_remarks.end() && remark_no != 350) {
                 _get_user_defined_remark(remark_no, user_pos->second);
                 continue;
            }

            if (remark_no == 2) {
                 if ((_experiment_type & EXPERIMENT_TYPE_XRAY) || (_experiment_type & EXPERIMENT_TYPE_NEUTRON) || (_experiment_type & EXPERIMENT_TYPE_EM) ||
                     (_experiment_type & EXPERIMENT_TYPE_ET) || (_experiment_type & EXPERIMENT_TYPE_FIBER) || (_experiment_type & EXPERIMENT_TYPE_ELECTRON) ||
                     (_experiment_type & EXPERIMENT_TYPE_EC)) {
                      _ndb_to_pdb_processing_RFACTR();
                      _ndb_to_pdb_get_general_remark(1, &Remark_2, 0, 0, 0);
                 } else _ndb_to_pdb_get_general_remark(1, &Remark_2_NMR, 0, 0, 0);
            } else if (remark_no == 3) {
                 if ((_experiment_type & EXPERIMENT_TYPE_XRAY) || (_experiment_type & EXPERIMENT_TYPE_NEUTRON) || (_experiment_type & EXPERIMENT_TYPE_FIBER) ||
                     (_experiment_type & EXPERIMENT_TYPE_ELECTRON) || (_experiment_type & EXPERIMENT_TYPE_EC))
                      _ndb_to_pdb_get_group_remark_3(refinement_program);
                 else if ((_experiment_type & EXPERIMENT_TYPE_NMR) || (_experiment_type & EXPERIMENT_TYPE_NMR_SOLID) ||
                          (_experiment_type & EXPERIMENT_TYPE_MODEL))
                      _ndb_to_pdb_get_general_remark(Num_Remark_3_NMR, Remark_3_NMR, 0, 0, 0);
                 else if ((_experiment_type & EXPERIMENT_TYPE_EM) || (_experiment_type & EXPERIMENT_TYPE_ET))
                      _ndb_to_pdb_get_general_remark(Num_Remark_3_Cryo_Em, Remark_3_Cryo_Em, 0, 0, 0);
                 else if (_experiment_type & EXPERIMENT_TYPE_SOLN_SCT) {
                      for (int i = 0; i < Remark_3_Solution.Remarks_No; ++i) {
                           _ndb_to_pdb_get_general_remark(Remark_3_Solution.num_remarks[i],
                                 Remark_3_Solution.remarks[i], 0, 0, 0);
                      }
                 }
            } else if (remark_no == 4) {
                 if (is_version_v315)
                      _ndb_to_pdb_get_general_remark(1, &Remark_4_V315, 0, 0, 0);
                 else _ndb_to_pdb_get_general_remark(1, &Remark_4, 0, 0, 0);
            } else if (remark_no == 5) {
                 _ndb_to_pdb_get_remark_5();
            } else if (remark_no == 100) {
                 _ndb_to_pdb_get_remark_100(is_unprocessed_model, EBI_ID);
            } else if (remark_no == 200) {
                 if (!is_powder_diffration) {
                      if (_experiment_type & EXPERIMENT_TYPE_XRAY)
                           _ndb_to_pdb_get_general_remark(Num_Remark_200, Remark_200, 0, 0, 0, "X-RAY DIFFRACTION");
                      else if (_experiment_type & EXPERIMENT_TYPE_FIBER)
                           _ndb_to_pdb_get_general_remark(Num_Remark_200, Remark_200, 0, 0, 0, "FIBER DIFFRACTION");
                 }
            } else if (remark_no == 205) {
                 if (_experiment_type & EXPERIMENT_TYPE_FIBER)
                      _ndb_to_pdb_get_general_remark(Num_Remark_205, Remark_205, 0, 0, 0);
            } else if (remark_no == 210) {
                 if (_experiment_type & EXPERIMENT_TYPE_NMR ||
                     _experiment_type & EXPERIMENT_TYPE_NMR_SOLID ||
                     _experiment_type & EXPERIMENT_TYPE_NMR_SOLID)
                      _ndb_to_pdb_get_general_remark(Num_Remark_210, Remark_210, 0, 0, 0);
            } else if (remark_no == 215) {
                 if  (_experiment_type & EXPERIMENT_TYPE_NMR && _cell.is_artifical())
                      _ndb_to_pdb_get_general_remark(Num_Remark_215, Remark_215, 0, 0, 0);
            } else if (remark_no == 217) {
                 if  (_experiment_type & EXPERIMENT_TYPE_NMR_SOLID && _cell.is_artifical())
                      _ndb_to_pdb_get_general_remark(Num_Remark_217, Remark_217, 0, 0, 0);
            } else if (remark_no == 220) {
                 if (_experiment_type & EXPERIMENT_TYPE_MODEL) {
                      std::string context;
                      _getRecordFront("MODTLS", 2, context);
                      if (!context.empty()) {
                           if (context[context.size() - 1] != '.')
                                context += ".";
                           context += " ";
                      }
                      context += "This theoretical model entry was not annotated and ";
                      context += "not validated by the wwPDB staff and therefore may ";
                      context += "not conform to the PDB format.";
                      _updateRecordFront("MODTLS", 2, context);

                      _ndb_to_pdb_get_general_remark(Num_Remark_220, Remark_220, 0, 0, 0);
                 }
            } else if (remark_no == 225) {
                 if (_experiment_type & EXPERIMENT_TYPE_MODEL)
                      _ndb_to_pdb_get_general_remark(Num_Remark_225, Remark_225, 0, 0, 0);
            } else if (remark_no == 230) {
                 if (_experiment_type & EXPERIMENT_TYPE_NEUTRON)
                      _ndb_to_pdb_get_general_remark(Num_Remark_230, Remark_230, 0, 0, 0, "NEUTRON DIFFRACTION");
            } else if (remark_no == 240) {
                 if (_experiment_type & EXPERIMENT_TYPE_EC)
                      _ndb_to_pdb_get_general_remark(Num_Remark_240, Remark_240, 0, 0, 0, "ELECTRON CRYSTALLOGRAPHY");
            } else if (remark_no == 245) {
                 if ((_experiment_type & EXPERIMENT_TYPE_EM) || (_experiment_type & EXPERIMENT_TYPE_ET))
                      _ndb_to_pdb_get_general_remark(Num_Remark_245, Remark_245, 0, 0, 0);
            } else if (remark_no == 247) {
                 if ((_experiment_type & EXPERIMENT_TYPE_EM) || (_experiment_type & EXPERIMENT_TYPE_ET))
                      _ndb_to_pdb_get_general_remark(Num_Remark_247, Remark_247, 0, 0, 0);
            } else if (remark_no == 250) {
                 if (_experiment_type & EXPERIMENT_TYPE_BASIC)
                      _ndb_to_pdb_get_general_remark(Num_Remark_250, Remark_250, 0, 0, 0);
            } else if (remark_no == 265) {
                 if (_experiment_type & EXPERIMENT_TYPE_SOLN_SCT)
                      _ndb_to_pdb_get_remark_265();
            } else if (remark_no == 280) {
                 if ((_experiment_type & EXPERIMENT_TYPE_XRAY) || (_experiment_type & EXPERIMENT_TYPE_NEUTRON))
                      _ndb_to_pdb_get_remark_280();
            } else if (remark_no == 290) {
                 if (is_electron_microscopy || (_experiment_type & EXPERIMENT_TYPE_XRAY) || (_experiment_type & EXPERIMENT_TYPE_NEUTRON) ||
                    (_experiment_type & EXPERIMENT_TYPE_ELECTRON) || (_experiment_type & EXPERIMENT_TYPE_EC)) _ndb_to_pdb_get_remark_290();
            } else if (remark_no == 300) _ndb_to_pdb_get_remark_300();
            else if (remark_no == 350) _ndb_to_pdb_get_remark_350();
            else if (remark_no == 375) {
                 if (_experiment_type & EXPERIMENT_TYPE_XRAY ||
                     _experiment_type & EXPERIMENT_TYPE_NEUTRON)
                      _ndb_to_pdb_get_remark_375();
            } else if (remark_no == 400) {
                 _ndb_to_pdb_get_remark(remark_no, 1, "COMPOUND");
                 _ndb_to_pdb_get_remark_400_group();
            } else if (remark_no == 450) _ndb_to_pdb_get_remark(remark_no, 2, "SOURCE");
            else if (remark_no == 465) _ndb_to_pdb_get_remark_465();
            else if (remark_no == 470) _ndb_to_pdb_get_remark_470();
            else if (remark_no == 475) _ndb_to_pdb_get_remark_475();
            else if (remark_no == 480) _ndb_to_pdb_get_remark_480();
            else if (remark_no == 500) _ndb_to_pdb_get_remark_500();
            else if (remark_no == 525) _ndb_to_pdb_get_remark_525();
            else if (remark_no == 600) _ndb_to_pdb_get_remark(remark_no, 3, "HETEROGEN");
            else if (remark_no == 610) _ndb_to_pdb_get_remark_610();
            else if (remark_no == 615) _ndb_to_pdb_get_remark_615();
            else if (remark_no == 620) _ndb_to_pdb_get_remark_620();
            else if (remark_no == 630) _ndb_to_pdb_get_remark_630();
            else if (remark_no == 800)
                 _ndb_to_pdb_get_general_remark(Num_Remark_800, Remark_800, 1, 0, 0);
            else if (remark_no == 900) _ndb_to_pdb_get_remark_900();
            else if (remark_no == 999) _ndb_to_pdb_get_remark(remark_no, 4, "SEQUENCE");
       }
}

bool Maxit::_ndb_to_pdb_get_skip_remark_number(std::set<int>& RemarkNo)
{
       RemarkNo.clear();

       std::string cs;
       _getRecordFront("PROCES", 7, cs);
       String::UpperCase(cs);
       if (cs[0] == 'Y') {
            RemarkNo.insert(500);
            RemarkNo.insert(525);
       }

       bool is_unprocessed_model = false;
       _getRecordFront("RMKSKP", 1, cs);
       if (!cs.empty()) {
            String::LowerCase(cs);
            std::vector<std::string> data;
            get_wordarray(data, cs, ",; ");
            for (std::vector<std::string>::const_iterator
                 pos = data.begin(); pos != data.end(); ++pos) {
                 if (*pos == "unprocessed") {
                      is_unprocessed_model = true;
                      RemarkNo.insert(500);
                 } else if (*pos == "y") {
                      RemarkNo.insert(500);
                      RemarkNo.insert(525);
                 } else if (String::IsNumber(*pos) && *pos != "350")
                      RemarkNo.insert(atoi(pos->c_str()));
            }
       }
       return is_unprocessed_model;
}

void Maxit::_ndb_to_pdb_get_auxiliary_data(std::string& ebi_id, bool& version_flag)
{
       ebi_id.clear();
       version_flag = false;

       if (!_CifObj) return;

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       std::string cs, dep_id;
       dep_id.clear();

       ISTable *t = getTablePtr(_cifblock, "database_2");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_upper(cs, t, i, "database_id");
                 if (cs == "EBI" || cs == "PDBE") {
                      get_value(ebi_id, t, i, "database_code");
                 } else if (cs == "WWPDB") {
                      get_value(dep_id, t, i, "database_code");
                 }
            }
       }

       t = getTablePtr(_cifblock, "pdbx_version");
       if (t) {
            int rowNo = t->GetNumRows();
            for (int i = rowNo - 1; i >= 0; --i) {
                 get_value(cs, t, i, "details");
                 if (cs == "compliance with PDB format V.3.15") {
                      version_flag = true;
                      if (_molecules.size() <= 1) version_flag = false;
                      break;
                 }
            }
       }

       // update process date
       if (!dep_id.empty()) _updateRecordFront("NDBPRO", 4, "", true);
       t = getTablePtr(_cifblock, "pdbx_database_status");
       if (t) {
            get_value_clean(cs, t, 0, "date_begin_processing");
            if (!cs.empty()) {
                 convert_time_format(9, cs);
                 _updateRecordFront("NDBPRO", 4, cs, true);
            }
       }

       // update computing related software
       t = getTablePtr(_cifblock, "software");
       if (t) {
            int num_mapping = 5;
            const char *_software_computing_mapping[5][4] = {
                 { "data reduction",  "_computing.pdbx_data_reduction_ii",  "DTMEAS", "5" },
                 { "data scaling",    "_computing.pdbx_data_reduction_ds",  "DTMEAS", "6" },
                 { "data collection", "_computing.data_collection",         "DTMEA1", "3" },
                 // { "model building",  "_computing.structure_solution",      "DTMEA1", "2" }, // old mapping
                 { "phasing",         "_computing.structure_solution",      "DTMEA1", "2" }, // new mapping
                 { "refinement",      "_computing.structure_refinement",    "REFMET", "4" }
            };

            std::set<std::string> classification;
            classification.clear();
            for (int i = 0; i < num_mapping; ++i) classification.insert(_software_computing_mapping[i][0]);

            std::map<std::string, std::vector<std::string> > software_value_mapping;
            software_value_mapping.clear();

            std::vector<std::string> data;
            std::string program, version;

            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_lower(cs, t, i, "classification");
                 if (cs.empty()) continue;
                 if (classification.find(cs) == classification.end()) continue;

                 get_value_clean(program, t, i, "name");
                 if (program.empty()) continue;
                 get_value_clean(version, t, i, "version");
                 if (!version.empty()) program += " " + version;

                 std::map<std::string, std::vector<std::string> >::iterator
                     mpos = software_value_mapping.find(cs);
                 if (mpos != software_value_mapping.end())
                      mpos->second.push_back(program);
                 else {
                      data.clear();
                      data.push_back(program);
                      software_value_mapping.insert(std::make_pair(cs, data));
                 }
            }

            if (!software_value_mapping.empty()) {
                 for (int i = 0; i < num_mapping; ++i) {
                      int field_no = atoi(_software_computing_mapping[i][3]) - 1;
                      // _getRecordFront(_software_computing_mapping[i][2], field_no, cs);
                      // if (!cs.empty()) continue;

                      std::map<std::string, std::vector<std::string> >::const_iterator
                          mpos = software_value_mapping.find(_software_computing_mapping[i][0]);
                      if (mpos == software_value_mapping.end()) continue;

                      cs.clear();
                      if (mpos->first == "refinement") {
                           cs = _select_correct_refinement_program(mpos->second);
                      } else {
                           for (std::vector<std::string>::const_iterator
                                pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                                if (!cs.empty()) cs += ", ";
                                cs += *pos;
                                // if (mpos->first == "refinement") break;
                           }
                      }
                      // always use value(s) from software category first.
                      _updateRecordFront(_software_computing_mapping[i][2], field_no, cs);
                 }
            }
       }
}

void Maxit::_get_user_defined_remark(const int& remark_no, const std::vector<std::string>&
                                     remarks)
{
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
       int width = ndbformat.FieldList[2].FieldWidth;

       std::vector<std::string> data;
       for (std::vector<std::string>::const_iterator
            pos = remarks.begin(); pos != remarks.end(); ++pos) {
            get_max_length_words(data, *pos, width);
            for (unsigned int i = 0; i < data.size(); ++i) {
                 if (remark_no != 350 && remark_no != 400 && remark_no != 465 &&
                     remark_no != 470 && remark_no != 500)
                      String::UpperCase(data[i]);
                 _addNewRemark(remark_no, data[i]);
            }
       }
}

void Maxit::_ndb_to_pdb_processing_RFACTR()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("RFACTR");
       if (ppos == _pdb_records.end()) _addNewRecord("RFACTR");

       std::string cs;
       _getRecordFront("RFACTR", 6, cs);
       if (cs.empty()) _getRecordFront("CMRNST", 3, cs);

       if (cs.empty()) return;

       std::string::size_type p = cs.find(".");
       if (p != std::string::npos) {
            std::string::size_type q = std::string::npos;
            for (std::string::size_type r = cs.size()-1; r > p+2; --r) {
                 if (cs[r] == '0') q = r;
                 else break;
            }
            if (q != std::string::npos) cs.erase(q);
       }

       _updateRecordFront("RFACTR", 6, cs);
}

void Maxit::_ndb_to_pdb_get_remark_5()
{
       if (!_CifObj) return;

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       ISTable *t = getTablePtr(_cifblock, "pdbx_database_PDB_obs_spr");
       if (!t) return;

       std::string cs, block_remark;
       block_remark.clear();

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean_upper(cs, t, i, "id");
            if (cs != "OBSLTE") continue;
            get_value_clean_upper(cs, t, i, "details");
            if (cs.empty()) continue;
            if (!block_remark.empty()) block_remark += " ";
            block_remark += cs;
       }
       if (block_remark.empty()) return;
 
       String::StripAndCompressWs(block_remark);
       if (block_remark.empty()) return;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
       int width = ndbformat.FieldList[2].FieldWidth;

       std::vector<std::string> data;
       get_max_length_words(data, block_remark, width);
       for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            _addNewRemark(5, *pos);
       }
       
}
