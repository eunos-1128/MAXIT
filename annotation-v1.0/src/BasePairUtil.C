/*
FILE:     BasePairUtil.C
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

#include <string.h>
#include <stdlib.h>
#include "AnnotationObj.h"
#include "BasePairInfo.h"
#include "CategoryMapping.h"
#include "utillib.h"

void AnnotationObj::CalculateBasePairInfo(const bool& force_flag)
{
       if (_molecules.empty()) return;

       if (!force_flag && !_basebase_params.empty() && !_interbase_params.empty()) return;

       BasePairInfo bpinfo;

       bpinfo.setCell(&_cell);
       bpinfo.setCCDic(_ccDic);
       bpinfo.setMolecule(_molecules[0]);
       bpinfo.calculateBasePairInfo();
       bpinfo.getBasePairs(_bspairs);
       bpinfo.getClassification(_classifications);
       _basebase_params = bpinfo.getBaseBaseParams();
       _interbase_params = bpinfo.getInterBaseParams();
}

void AnnotationObj::Update_Base_Pair_and_Step()
{
       if (!_CifObj) return;

       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
       // _update_ndb_struct_na_base_pair(_cifblock);
       // _update_ndb_struct_na_base_pair_step(_cifblock);
       _update_ndb_struct_conf_na(_cifblock);
}

void AnnotationObj::_read_ndb_struct_na_base_pair(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "ndb_struct_na_base_pair");
       if (!t) return;

       std::string cs, pdb_chnid, res_name, res_num, ins_code;
       _BASEBASE_PARAMS params;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(pdb_chnid, t, i, "i_auth_asym_id");
            get_value_clean(res_name, t, i, "i_label_comp_id");
            get_value_clean(res_num, t, i, "i_auth_seq_id");
            get_value_clean(ins_code, t, i, "i_PDB_ins_code");
            RCSB::Residue* res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name,
                                  res_num, ins_code);
            if (!res) continue;
            params.I_index = res->index();

            get_value_clean(pdb_chnid, t, i, "j_auth_asym_id");
            get_value_clean(res_name, t, i, "j_label_comp_id");
            get_value_clean(res_num, t, i, "j_auth_seq_id");
            get_value_clean(ins_code, t, i, "j_PDB_ins_code");
            res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!res) continue;
            params.J_index = res->index();

            get_value_clean(params.Sym_I, t, i, "i_symmetry");
            get_value_clean(params.Sym_J, t, i, "j_symmetry");
            get_value_clean(params._saenger_classification, t, i, "hbond_type_28");
            get_value_clean(params._leontis_westhof_classification,  t, i, "hbond_type_12");
            
            get_value_clean(cs, t, i, "shear");   params.Shear   = atof(cs.c_str());
            get_value_clean(cs, t, i, "stretch"); params.Stretch = atof(cs.c_str());
            get_value_clean(cs, t, i, "stagger"); params.Stagger = atof(cs.c_str());
            get_value_clean(cs, t, i, "buckle");  params.Buckle  = atof(cs.c_str());
            get_value_clean(cs, t, i, "propeller"); params.Propel  = atof(cs.c_str());
            get_value_clean(cs, t, i, "opening"); params.Opening = atof(cs.c_str());

            _basebase_params.push_back(params);
       }
}

void AnnotationObj::_update_ndb_struct_na_base_pair(Block& block)
{
       if (_basebase_params.empty()) {
            deleteTable(block, "ndb_struct_na_base_pair");
            return;
       }
       if (_molecules.empty()) return;

       ISTable *t = _newTablePtr("ndb_struct_na_base_pair");

       int i = 0;
       for (std::list<_BASEBASE_PARAMS>::const_iterator pos = _basebase_params.begin(); pos != _basebase_params.end(); ++pos) {
            bool is_removed = false;
            RCSB::Residue* I = _molecules[0]->find_residue(pos->I_index, is_removed);
            if (!I) continue;

            RCSB::Residue* J = _molecules[0]->find_residue(pos->J_index, is_removed);
            if (!J) continue;

            std::string buff = I->pdb_chnid() + "_" + I->ResName() + I->pdb_res_no() + I->ins_code() + ":"
                             + J->ResName() + J->pdb_res_no() + J->ins_code() + "_" + J->pdb_chnid();
            t->AddRow();
            t->UpdateCell(i, "model_number",    "1");
            t->UpdateCell(i, "i_label_asym_id", I->chnid());
            t->UpdateCell(i, "i_label_comp_id", I->ResName());
            t->UpdateCell(i, "i_label_seq_id",  I->res_no());
            t->UpdateCell(i, "i_symmetry",      pos->Sym_I);
            t->UpdateCell(i, "j_label_asym_id", J->chnid());
            t->UpdateCell(i, "j_label_comp_id", J->ResName());
            t->UpdateCell(i, "j_label_seq_id",  J->res_no());
            t->UpdateCell(i, "j_symmetry",      pos->Sym_J);
            t->UpdateCell(i, "shear",           FloatToString(pos->Shear, 0, 3));
            t->UpdateCell(i, "stretch",         FloatToString(pos->Stretch, 0, 3));
            t->UpdateCell(i, "stagger",         FloatToString(pos->Stagger, 0, 3));
            t->UpdateCell(i, "buckle",          FloatToString(pos->Buckle, 0, 3));
            t->UpdateCell(i, "propeller",       FloatToString(pos->Propel, 0, 3));
            t->UpdateCell(i, "opening",         FloatToString(pos->Opening, 0, 3));
            t->UpdateCell(i, "pair_number",     String::IntToString(i + 1));
            t->UpdateCell(i, "pair_name",       buff);
            t->UpdateCell(i, "i_auth_asym_id",  I->pdb_chnid());
            t->UpdateCell(i, "i_auth_seq_id",   I->pdb_res_no());
            t->UpdateCell(i, "i_PDB_ins_code",  I->ins_code());
            t->UpdateCell(i, "j_auth_asym_id",  J->pdb_chnid());
            t->UpdateCell(i, "j_auth_seq_id",   J->pdb_res_no());
            t->UpdateCell(i, "j_PDB_ins_code",  J->ins_code());
            t->UpdateCell(i, "hbond_type_28",   pos->_saenger_classification);
            t->UpdateCell(i, "hbond_type_12",   pos->_leontis_westhof_classification);
            i++;
       }
       block.WriteTable(t);
}

void AnnotationObj::_read_ndb_struct_na_base_pair_step(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "ndb_struct_na_base_pair_step");
       if (!t) return;

       std::string cs, pdb_chnid, res_name, res_num, ins_code;
       _INTERBASE_PARAMS params;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(pdb_chnid, t, i, "i_auth_asym_id_1");
            get_value_clean(res_name, t, i, "i_label_comp_id_1");
            get_value_clean(res_num, t, i, "i_auth_seq_id_1");
            get_value_clean(ins_code, t, i, "i_PDB_ins_code_1");
            RCSB::Residue* res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name,
                                  res_num, ins_code);
            if (!res) continue;
            params.I1_index = res->index();

            get_value_clean(pdb_chnid, t, i, "j_auth_asym_id_1");
            get_value_clean(res_name, t, i, "j_label_comp_id_1");
            get_value_clean(res_num, t, i, "j_auth_seq_id_1");
            get_value_clean(ins_code, t, i, "j_PDB_ins_code_1");
            res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!res) continue;
            params.J1_index = res->index();

            get_value_clean(pdb_chnid, t, i, "i_auth_asym_id_2");
            get_value_clean(res_name, t, i, "i_label_comp_id_2");
            get_value_clean(res_num, t, i, "i_auth_seq_id_2");
            get_value_clean(ins_code, t, i, "i_PDB_ins_code_2");
            res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!res) continue;
            params.I2_index = res->index();

            get_value_clean(pdb_chnid, t, i, "j_auth_asym_id_2");
            get_value_clean(res_name, t, i, "j_label_comp_id_2");
            get_value_clean(res_num, t, i, "j_auth_seq_id_2");
            get_value_clean(ins_code, t, i, "j_PDB_ins_code_2");
            res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!res) continue;
            params.J2_index = res->index();

            get_value_clean(params.Sym_I1, t, i, "i_symmetry_1");
            get_value_clean(params.Sym_J1, t, i, "j_symmetry_1");
            get_value_clean(params.Sym_I2, t, i, "i_symmetry_2");
            get_value_clean(params.Sym_J2, t, i, "j_symmetry_2");

            get_value_clean(cs, t, i, "shift");          params.Shift   = atof(cs.c_str());
            get_value_clean(cs, t, i, "slide");          params.Slide   = atof(cs.c_str());
            get_value_clean(cs, t, i, "rise");           params.Rise    = atof(cs.c_str());
            get_value_clean(cs, t, i, "tilt");           params.Tilt    = atof(cs.c_str());
            get_value_clean(cs, t, i, "roll");           params.Roll    = atof(cs.c_str());
            get_value_clean(cs, t, i, "twist");          params.Twist   = atof(cs.c_str());
            get_value_clean(cs, t, i, "x_displacement"); params.X_disp  = atof(cs.c_str());
            get_value_clean(cs, t, i, "y_displacement"); params.Y_disp  = atof(cs.c_str());
            get_value_clean(cs, t, i, "helical_rise");   params.H_rise  = atof(cs.c_str());
            get_value_clean(cs, t, i, "inclination");    params.Incl    = atof(cs.c_str());
            get_value_clean(cs, t, i, "tip");            params.Tip     = atof(cs.c_str());
            get_value_clean(cs, t, i, "helical_twist");  params.H_twist = atof(cs.c_str());

            _interbase_params.push_back(params);
       }
}

void AnnotationObj::_update_ndb_struct_na_base_pair_step(Block& block)
{
       if (_interbase_params.empty()) {
            deleteTable(block, "ndb_struct_na_base_pair_step");
            return;
       }
       if (_molecules.empty()) return;

       ISTable *t = _newTablePtr("ndb_struct_na_base_pair_step");

       int i = 0;
       for (std::list<_INTERBASE_PARAMS>::const_iterator pos = _interbase_params.begin(); pos != _interbase_params.end(); ++pos) {
            bool is_removed = false;
            RCSB::Residue* I1 = _molecules[0]->find_residue(pos->I1_index, is_removed);
            if (!I1) continue;

            RCSB::Residue* J1 = _molecules[0]->find_residue(pos->J1_index, is_removed);
            if (!J1) continue;

            RCSB::Residue* I2 = _molecules[0]->find_residue(pos->I2_index, is_removed);
            if (!I2) continue;

            RCSB::Residue* J2 = _molecules[0]->find_residue(pos->J2_index, is_removed);
            if (!J2) continue;

            std::string buff = I1->pdb_chnid() + I2->pdb_chnid() + "_" + I1->ResName() + I1->pdb_res_no() + I1->ins_code() + I2->ResName()
                             + I2->pdb_res_no() + I2->ins_code() + ":" + J2->ResName() + J2->pdb_res_no() + J2->ins_code() + J1->ResName()
                             + J1->pdb_res_no() + J1->ins_code() + "_" + J2->pdb_chnid() + J1->pdb_chnid();
            t->AddRow();
            t->UpdateCell(i, "model_number", "1");
            t->UpdateCell(i, "i_label_asym_id_1", I1->chnid());
            t->UpdateCell(i, "i_label_comp_id_1", I1->ResName());
            t->UpdateCell(i, "i_label_seq_id_1",  I1->res_no());
            t->UpdateCell(i, "i_symmetry_1",      pos->Sym_I1);
            t->UpdateCell(i, "j_label_asym_id_1", J1->chnid());
            t->UpdateCell(i, "j_label_comp_id_1", J1->ResName());
            t->UpdateCell(i, "j_label_seq_id_1",  J1->res_no());
            t->UpdateCell(i, "j_symmetry_1",      pos->Sym_J1);
            t->UpdateCell(i, "i_label_asym_id_2", I2->chnid());
            t->UpdateCell(i, "i_label_comp_id_2", I2->ResName());
            t->UpdateCell(i, "i_label_seq_id_2",  I2->res_no());
            t->UpdateCell(i, "i_symmetry_2",      pos->Sym_I2);
            t->UpdateCell(i, "j_label_asym_id_2", J2->chnid());
            t->UpdateCell(i, "j_label_comp_id_2", J2->ResName());
            t->UpdateCell(i, "j_label_seq_id_2",  J2->res_no());
            t->UpdateCell(i, "j_symmetry_2",      pos->Sym_J2);
            t->UpdateCell(i, "shift",             FloatToString(pos->Shift, 0, 3));
            t->UpdateCell(i, "slide",             FloatToString(pos->Slide, 0, 3));
            t->UpdateCell(i, "rise",              FloatToString(pos->Rise, 0, 3));
            t->UpdateCell(i, "tilt",              FloatToString(pos->Tilt, 0, 3));
            t->UpdateCell(i, "roll",              FloatToString(pos->Roll, 0, 3));
            t->UpdateCell(i, "twist",             FloatToString(pos->Twist, 0, 3));
            t->UpdateCell(i, "x_displacement",    FloatToString(pos->X_disp, 0, 3));
            t->UpdateCell(i, "y_displacement",    FloatToString(pos->Y_disp, 0, 3));
            t->UpdateCell(i, "helical_rise",      FloatToString(pos->H_rise, 0, 3));
            t->UpdateCell(i, "inclination",       FloatToString(pos->Incl, 0, 3));
            t->UpdateCell(i, "tip",               FloatToString(pos->Tip, 0, 3));
            t->UpdateCell(i, "helical_twist",     FloatToString(pos->H_twist, 0, 3));
            t->UpdateCell(i, "step_number",       String::IntToString(i + 1));
            t->UpdateCell(i, "step_name",         buff);
            t->UpdateCell(i, "i_auth_asym_id_1",  I1->pdb_chnid());
            t->UpdateCell(i, "i_auth_seq_id_1",   I1->pdb_res_no());
            t->UpdateCell(i, "i_PDB_ins_code_1",  I1->ins_code());
            t->UpdateCell(i, "j_auth_asym_id_1",  J1->pdb_chnid());
            t->UpdateCell(i, "j_auth_seq_id_1",   J1->pdb_res_no());
            t->UpdateCell(i, "j_PDB_ins_code_1",  J1->ins_code());
            t->UpdateCell(i, "i_auth_asym_id_2",  I2->pdb_chnid());
            t->UpdateCell(i, "i_auth_seq_id_2",   I2->pdb_res_no());
            t->UpdateCell(i, "i_PDB_ins_code_2",  I2->ins_code());
            t->UpdateCell(i, "j_auth_asym_id_2",  J2->pdb_chnid());
            t->UpdateCell(i, "j_auth_seq_id_2",   J2->pdb_res_no());
            t->UpdateCell(i, "j_PDB_ins_code_2",  J2->ins_code());
            i++;
       }
       block.WriteTable(t);
}

void AnnotationObj::_update_ndb_struct_conf_na(Block& block)
{
       if (_classifications.empty()) {
            deleteTable(block, "ndb_struct_conf_na");
            return;
       }

       ISTable *t = _newTablePtr("ndb_struct_conf_na");

       int i = 0;
       for (std::vector<std::string>::const_iterator
            pos = _classifications.begin(); pos != _classifications.end(); ++pos) {
            t->AddRow();
            t->UpdateCell(i, "entry_id", _firstBlockName);
            t->UpdateCell(i, "feature", *pos);
            i++;
       }
       block.WriteTable(t);
}
