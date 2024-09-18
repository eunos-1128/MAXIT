/*
FILE:     FileObj.C
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
#include <unistd.h>
#include <sys/stat.h>

#include "CompositeIndex.h"
#include "FileObj.h"
#include "PdbWrite.h"
#include "SpaceGroup.h"
#include "SplitUtil.h"
#include "TypeDef.h"
#include "utillib.h"

#define NUM_GLYCO_SITE 28

const static char *__glyco_site_list[NUM_GLYCO_SITE][5] = {
       { "ASN", "ND2", "BGC", "C1", "N-Glycosylation" },
       { "ASN", "ND2", "GAL", "C1", "N-Glycosylation" },
       { "ASN", "ND2", "NAG", "C1", "N-Glycosylation" },
       { "ASN", "ND2", "NGA", "C1", "N-Glycosylation" },
       { "ARG", "NH1", "NAG", "C1", "N-Glycosylation" },
       { "ARG", "NH2", "NAG", "C1", "N-Glycosylation" },
       { "CYS", "SG",  "A2G", "C1", "S-Glycosylation" },
       { "CYS", "SG",  "BGC", "C1", "S-Glycosylation" },
       { "CYS", "SG",  "NAG", "C1", "S-Glycosylation" },
       { "SER", "OG",  "A2G", "C1", "O-Glycosylation" },
       { "SER", "OG",  "B6D", "C1", "O-Glycosylation" },
       { "SER", "OG",  "BGC", "C1", "O-Glycosylation" },
       { "SER", "OG",  "FUC", "C1", "O-Glycosylation" },
       { "SER", "OG",  "GLA", "C1", "O-Glycosylation" },
       { "SER", "OG",  "GLC", "C1", "O-Glycosylation" },
       { "SER", "OG",  "MAN", "C1", "O-Glycosylation" },
       { "SER", "OG",  "NAG", "C1", "O-Glycosylation" },
       { "THR", "OG1", "A2G", "C1", "O-Glycosylation" },
       { "THR", "OG1", "FUC", "C1", "O-Glycosylation" },
       { "THR", "OG1", "GLC", "C1", "O-Glycosylation" },
       { "THR", "OG1", "MAN", "C1", "O-Glycosylation" },
       { "THR", "OG1", "NAG", "C1", "O-Glycosylation" },
       { "THR", "OG1", "NDG", "C1", "O-Glycosylation" },
       { "THR", "OG1", "NGA", "C1", "O-Glycosylation" },
       { "TYR", "OH",  "BGC", "C1", "O-Glycosylation" },
       { "TYR", "OH",  "GLC", "C1", "O-Glycosylation" },
       { "TYR", "OH",  "NDG", "C1", "O-Glycosylation" },
       { "TRP", "CD1", "MAN", "C1", "C-Mannosylation" }
};

FileObj::FileObj()
{
       _initialize_reference_variables();
       _initialize_other_variables();
       _initExpTypeMapping();

       _fill_in_glyco_site_info();
}

FileObj::~FileObj()
{
       _clear_allocated_memories();
       _initialize_reference_variables();
       _initialize_other_variables();
}

void FileObj::_initialize_reference_variables()
{
       _glyco_site_map.clear();

       _logIo = NULL;
       _ccDic = NULL;
       _error_messages.clear();
       _rcsbroot.clear();
       _dictUtil.clear();
       _ExpTypeMapping.clear();
       _TypeExpMapping.clear();
}

void FileObj::_initialize_other_variables()
{
       _read_atom_site_record_number = 0;
       _write_atom_site_record_number = 0;
       _check_atom_site_record_number_flag = true;

       _input_filename.clear();
       _output_filename.clear();
       _StructureId.clear();
       _UppercasePdbId.clear();
       _LowercasePdbId.clear();
       _firstBlockName.clear();
       _resetFirstBlockName.clear();
       _input_format = 0;
       _output_format = 0;
       _experiment_type = 0;
       _DepUI_Flag = false;
       _pdb_records.clear();
       _molecules.clear();
       _caseSensitive = true;
       _rename_residue_flag = false;
       _find_sugar_chain_flag = true;
       _skip_split_flag = false;
       _CifObj = NULL;
       _Seqs.clear();
       _Seq_Scheme_Mapping.clear();
       _Branch_Seq_Scheme_Mapping.clear();
       _cell.clear();
       _exist_scale_matrix = false;
       _SymmMatrices.clear();
       _SymmMatrix_Mapping.clear();
       _assemblies.clear();
       _update_symmetry_flag = true;
       _max_entity_id = 0;
       _entities.clear();
       _entity_ids.clear();
       _skip_task_set.clear();
       _glyco_site_res_pair_set.clear();
       _merge_residue_list.clear();
       _link_residue_set.clear();
       _extra_atom_site_item_list.clear();
}

void FileObj::_initExpTypeMapping()
{
       _ExpTypeMapping.insert(std::make_pair("X-RAY DIFFRACTION", EXPERIMENT_TYPE_XRAY));
       _ExpTypeMapping.insert(std::make_pair("X-RAY POWDER DIFFRACTION", EXPERIMENT_TYPE_XRAY));
       _ExpTypeMapping.insert(std::make_pair("NEUTRON DIFFRACTION", EXPERIMENT_TYPE_NEUTRON));
       _ExpTypeMapping.insert(std::make_pair("FIBER DIFFRACTION", EXPERIMENT_TYPE_FIBER));
       _ExpTypeMapping.insert(std::make_pair("FIBRE DIFFRACTION", EXPERIMENT_TYPE_FIBER));
       _ExpTypeMapping.insert(std::make_pair("ELECTRON DIFFRACTION", EXPERIMENT_TYPE_ELECTRON));
       _ExpTypeMapping.insert(std::make_pair("ELECTRON CRYSTALLOGRAPHY", EXPERIMENT_TYPE_EC));
       _ExpTypeMapping.insert(std::make_pair("CRYO-ELECTRON MICROSCOPY", EXPERIMENT_TYPE_EM));
       _ExpTypeMapping.insert(std::make_pair("ELECTRON TOMOGRAPHY", EXPERIMENT_TYPE_ET));
       _ExpTypeMapping.insert(std::make_pair("ELECTRON MICROSCOPY", EXPERIMENT_TYPE_EM));
       _ExpTypeMapping.insert(std::make_pair("NMR", EXPERIMENT_TYPE_NMR));
       _ExpTypeMapping.insert(std::make_pair("SOLUTION NMR", EXPERIMENT_TYPE_NMR));
       _ExpTypeMapping.insert(std::make_pair("SOLID-STATE NMR", EXPERIMENT_TYPE_NMR_SOLID));
       _ExpTypeMapping.insert(std::make_pair("SOLUTION SCATTERING", EXPERIMENT_TYPE_SOLN_SCT));
       _ExpTypeMapping.insert(std::make_pair("THEORETICAL MODEL", EXPERIMENT_TYPE_MODEL));

       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_XRAY, "X-RAY DIFFRACTION"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_NEUTRON, "NEUTRON DIFFRACTION"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_FIBER, "FIBER DIFFRACTION"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_ELECTRON, "ELECTRON DIFFRACTION"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_EC, "ELECTRON CRYSTALLOGRAPHY"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_EM, "ELECTRON MICROSCOPY"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_ET, "ELECTRON TOMOGRAPHY"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_NMR, "SOLUTION NMR"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_NMR_SOLID, "SOLID-STATE NMR"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_SOLN_SCT, "SOLUTION SCATTERING"));
       _TypeExpMapping.insert(std::make_pair(EXPERIMENT_TYPE_MODEL, "THEORETICAL MODEL"));
}

std::string FileObj::_get_diffraction_experiment_method()
{
       std::string exp_type = "";
       if (_experiment_type & EXPERIMENT_TYPE_XRAY)
            exp_type = "X-RAY DIFFRACTION";
       else if (_experiment_type & EXPERIMENT_TYPE_NEUTRON)
            exp_type = "NEUTRON DIFFRACTION";
       else if (_experiment_type & EXPERIMENT_TYPE_FIBER)
            exp_type = "FIBER DIFFRACTION";
       else if (_experiment_type & EXPERIMENT_TYPE_ELECTRON)
            exp_type = "ELECTRON DIFFRACTION";
       else if (_experiment_type & EXPERIMENT_TYPE_EC)
            exp_type = "ELECTRON CRYSTALLOGRAPHY";
       else if (_experiment_type & EXPERIMENT_TYPE_EM)
            exp_type = "ELECTRON MICROSCOPY";
       else if (_experiment_type & EXPERIMENT_TYPE_ET)
            exp_type = "ELECTRON TOMOGRAPHY";
       else if (_experiment_type & EXPERIMENT_TYPE_SOLN_SCT)
            exp_type = "SOLUTION SCATTERING";
       return exp_type;
}

void FileObj::_clear_allocated_memories()
{
       if (!_molecules.empty()) {
            for (unsigned int i = 0; i < _molecules.size(); ++i) {
                 if (_molecules[i]) delete _molecules[i];
            }
       }
       if(_CifObj) delete _CifObj;
}

void FileObj::_clearSEQ(SEQ& seq)
{
       seq.ChainId.clear();
       seq.PDB_ChainId.clear();
       seq.PolyType.clear();
       seq.EvidenceCode.clear();
       seq.chain_type.clear();
       seq.entity_id.clear();
       seq.never_used = true;
       seq.order = 0;
       seq.res.clear();
}

void FileObj::_fill_in_glyco_site_info()
{
       for (int i = 0; i < NUM_GLYCO_SITE; ++i) {
            _glyco_site_map.insert(std::make_pair(CompositeIndex::getIndex(__glyco_site_list[i][0], __glyco_site_list[i][1],
                                   __glyco_site_list[i][2], __glyco_site_list[i][3]), __glyco_site_list[i][4]));
            _glyco_site_map.insert(std::make_pair(CompositeIndex::getIndex(__glyco_site_list[i][2], __glyco_site_list[i][3],
                                   __glyco_site_list[i][0], __glyco_site_list[i][1]), __glyco_site_list[i][4]));
            _glyco_site_res_pair_set.insert(CompositeIndex::getIndex(__glyco_site_list[i][0], __glyco_site_list[i][2]));
            _glyco_site_res_pair_set.insert(CompositeIndex::getIndex(__glyco_site_list[i][2], __glyco_site_list[i][0]));
       }
}

bool FileObj::_is_glycosylation_site(const std::string& idx)
{
       std::string type = _find_glycosylation_type(idx);
       if (!type.empty()) return true;
       return false;
}

bool FileObj::_is_glycosylation_site(RCSB::Atom* atom1, RCSB::Atom* atom2)
{
       std::string type = _find_glycosylation_type(atom1, atom2);
       if (!type.empty()) return true;
       return false;
}

std::string FileObj::_find_glycosylation_type(const std::string& idx)
{
       std::map<std::string, std::string>::const_iterator mpos = _glyco_site_map.find(idx);
       if (mpos != _glyco_site_map.end()) return mpos->second;
       return "";
}

std::string FileObj::_find_glycosylation_type(RCSB::Atom* atom1, RCSB::Atom* atom2)
{
       return _find_glycosylation_type(CompositeIndex::getIndex(atom1->pdb_resnam(), atom1->atmtype(), atom2->pdb_resnam(), atom2->atmtype()));
}

void FileObj::setLog(LogUtil *logPt)                     { _logIo = logPt; }
void FileObj::setCCDic(ConnectDic *ccdic)                { _ccDic = ccdic; }
void FileObj::setRCSBROOT(const std::string& cs)         { _rcsbroot = cs; _dictUtil.setRCSBROOT(_rcsbroot); }
void FileObj::setDictSdbPath(const std::string& cs)      { _dictUtil.setDictSdbPath(cs); }
void FileObj::set_input_filename(const std::string& cs)  { _input_filename = cs; }
void FileObj::set_output_filename(const std::string& cs) { _output_filename = cs; }
void FileObj::set_input_format(const int format)         { _input_format = format; }
void FileObj::set_output_format(const int format)        { _output_format = format; }
void FileObj::set_depUI_flag(const bool& flag)           { _DepUI_Flag = flag; }
void FileObj::set_CaseSensitive()                        { _caseSensitive = true; }
void FileObj::setRenameResidueFlag(const bool& flag)     { _rename_residue_flag = flag; }
void FileObj::setFindSugarChainFlag(const bool& flag)    { _find_sugar_chain_flag = flag; }
void FileObj::setSkipSplitFlag(const bool& flag)         { _skip_split_flag = flag; }

void FileObj::set_merge_residue_info(const std::vector<std::string>& groups)
{
       if (groups.empty()) return;

       std::vector<std::string> data, data1;
       std::vector<std::vector<std::string> > residue_list;
       for (std::vector<std::string>::const_iterator pos = groups.begin(); pos != groups.end(); ++pos) {
            get_wordarray(data, *pos, ",");
            if ((data.size() < 3) || (data[0].size() > 3)) continue;

            residue_list.clear();
            for (unsigned int i = 1; i < data.size(); ++i) {
                 get_wordarray(data1, data[i], "_");
                 if (data1.size() == 3) data1.push_back("");
                 if (data1.size() != 4) continue;
                 residue_list.push_back(data1);
            }
            if (residue_list.size() < 2) continue;
            _merge_residue_list.push_back(std::make_pair(data[0], residue_list));
       }
}

void FileObj::_getExperimentType(const std::string& exp_type)
{
       std::map<std::string, int>::const_iterator mpos = _ExpTypeMapping.find(exp_type);
       if (mpos != _ExpTypeMapping.end()) _experiment_type |= mpos->second;
}

void FileObj::_updateCrySymmetry()
{
       _cell.check_artifical_flag((_experiment_type & EXPERIMENT_TYPE_XRAY) || (_experiment_type & EXPERIMENT_TYPE_NEUTRON) ||
                                  (_experiment_type & EXPERIMENT_TYPE_FIBER));

       if (((_experiment_type & EXPERIMENT_TYPE_XRAY) || (_experiment_type & EXPERIMENT_TYPE_NEUTRON) ||
            (_experiment_type & EXPERIMENT_TYPE_FIBER)) && _cell.space_group().empty()) {
            _error_messages.insertMessage("space_group_print", "warning", "Missing space group name.");
       }

       std::vector<std::vector<int> > sym_data;
       sym_data.clear();
       if (!_cell.is_artifical() && !_cell.space_group().empty()) {
            std::string IT_Number, cs;
            std::string name  = _cell.space_group();
            SpaceGroup::CheckSpaceGroup(name, IT_Number);
            if (name != _cell.space_group()) _cell.setSpaceGroup(name);
            SpaceGroup::GetSpaceInfo(name, sym_data, cs);
            if (!cs.empty()) {
                 _logIo->message(cs.c_str());
                 _error_messages.insertMessage("space_group_print", "error", cs);
            }
            if (sym_data.empty()) _cell.setArtificalFlag(2);

            if (_cell.a().empty() || _cell.a().substr(0, 2) == "1." || _cell.b().empty() || _cell.b().substr(0, 2) == "1." ||
                _cell.c().empty() || _cell.c().substr(0, 2) == "1." || _cell.alpha().empty() || _cell.beta().empty() || _cell.gamma().empty()) {
                 _cell.setArtificalFlag(2);
                 cs = "Error: Space group=" + _cell.space_group() + " alpha=" + _cell.alpha() + " beta=" + _cell.beta() + " gamma="
                    + _cell.gamma() + " a=" + _cell.a() + " b=" + _cell.b() + " c=" + _cell.c() + ".";
                 _error_messages.insertMessage("space_group", "label", "true");
                 _error_messages.insertMessage("space_group_print", "error", cs);
            }
       }

       _cell.SetCrySymmetry(sym_data);
       if (!_cell.is_artifical() && _exist_scale_matrix) {
            if (!_cell.check_scale_matrix(_scale, _vec)) {
                 if (_CifObj) {
                      Block& block = _CifObj->GetBlock(_firstBlockName);
                      deleteTable(block, "atom_sites");
                 }
                 if (_pdb_records.find("SCALE") != _pdb_records.end()) _pdb_records.erase("SCALE");
            }
       } 
}

void FileObj::_getSequenceType()
{
       std::string cs;

       for (std::map<std::string, SEQ>::iterator mpos = _Seqs.begin(); mpos != _Seqs.end(); ++mpos) {
            String::LowerCase(mpos->second.PolyType, cs);
            if (cs.find("peptide") != std::string::npos)
                 mpos->second.chain_type = "ATOMP";
            else if (cs.find("ribonucleotide") != std::string::npos)
                 mpos->second.chain_type = "ATOMN";
            else if (cs.find("saccharide") != std::string::npos)
                 mpos->second.chain_type = "ATOMS";
            else mpos->second.chain_type = _getSequenceType(mpos->second.res);
       }
}

std::string FileObj::_getSequenceType(const std::vector<std::string>& res)
{
       std::map<std::string, int> all_types, selected_types;
       all_types.clear();
       selected_types.clear();

       bool has_ATOMP = false;
       for (std::vector<std::string>::const_iterator pos = res.begin(); pos != res.end(); ++pos) {
            std::string cs = _ccDic->find_residue_type(*pos);
            if (cs == "ATOMP") has_ATOMP = true;
            add_type(all_types, cs);
            if (cs == "ATOMP" || cs == "ATOMN" || cs == "ATOMS") add_type(selected_types, cs);
       }
       if (selected_types.size() == 1) {
            std::map<std::string, int>::const_iterator mpos = selected_types.begin();
            return mpos->first;
       }

       std::string chain_type = get_type(all_types);
       if ((chain_type != "ATOMP") && has_ATOMP && (res.size() > 1) && (res.size() < 11)) chain_type = "ATOMP";
       return chain_type;
}

void FileObj::_getResiduesfromNames(const int& mol_id, const std::vector<std::vector<std::string> >& data, int& mol_index,
                                    std::vector<RCSB::Residue*>& residues, std::vector<std::string>& errors)
{
       mol_index = 0;
       residues.clear();
       errors.clear();

       if (_molecules.empty()) return;

       RCSB::Molecule* mol = NULL;
       if (mol_id == 0)
            mol = _molecules[0];
       else {
            for (unsigned int i = 0; i < _molecules.size(); ++i) {
                 if (_molecules[i]->Mol_ID() == mol_id) {
                      mol_index = i;
                      mol = _molecules[i];
                      break;
                 }
            }
       }
       if (!mol) return;

       for (std::vector<std::vector<std::string> >::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            RCSB::Residue*res = mol->find_pdb_residue((*pos)[0], (*pos)[1], (*pos)[2], (*pos)[3]);
            if (!res) {
                 std::string cs = "Residue ('" + (*pos)[0] + "' '" + (*pos)[1] + "' '" + (*pos)[2] + (*pos)[3] + "') can not be found.";
                 if (_skip_split_flag) {
                      _SPLIT_RESIDUE_LIST *exist = SplitUtil::split_residue((*pos)[1]);
                      if (exist) continue;
                 }
                 errors.push_back(cs);
                 continue;
            }
            residues.push_back(res);
       }
       if (residues.size() != data.size()) residues.clear();
}

void FileObj::_getAtomsfromNames(const int& mol_id, const std::vector<std::vector<std::string> >& data, int& mol_index,
                                 std::vector<RCSB::Atom*>& atoms, std::vector<std::string>& errors, const bool& origin_flag)
{
       mol_index = 0;
       atoms.clear();
       errors.clear();

       if (_molecules.empty()) return;

       RCSB::Molecule* mol = NULL;
       if (mol_id == 0)
            mol = _molecules[0];
       else {
            for (unsigned int i = 0; i < _molecules.size(); ++i) {
                 if (_molecules[i]->Mol_ID() == mol_id) {
                      mol_index = i;
                      mol = _molecules[i];
                      break;
                 }
            }
       }
       if (!mol) return;

       for (std::vector<std::vector<std::string> >::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            RCSB::Residue*res = mol->find_pdb_residue((*pos)[0], (*pos)[1], (*pos)[2], (*pos)[3]);
            if (!res) {
                 std::string cs = "Residue ('" + (*pos)[0] + "' '" + (*pos)[1] + "' '" + (*pos)[2] + (*pos)[3] + "') can not be found.";
                 if (_skip_split_flag) {
                      _SPLIT_RESIDUE_LIST *exist = SplitUtil::split_residue((*pos)[1]);
                      if (exist) continue;
                 }
                 errors.push_back(cs);
                 continue;
            }
            RCSB::Atom* atom = res->find_atom((*pos)[4], (*pos)[5], origin_flag);
            if (atom == NULL && (*pos)[5].empty()) atom = res->find_atom((*pos)[4]);
            if (!atom) {
                 std::string cs = "Atom ('" + (*pos)[4] + (*pos)[5] + "' '" + (*pos)[0] + "' '" + (*pos)[1] + "' '" + (*pos)[2]
                                + (*pos)[3] + "') can not be found.";
                 errors.push_back(cs);
                 continue;
            }
            atoms.push_back(atom);
       }
       if (atoms.size() != data.size()) atoms.clear();
}

bool FileObj::build_molecule(const bool& skip_update_flag, const bool& split_flag, const std::string& assignfile)
{
       if (_molecules.empty()) return false;

       // key: Instance ID
       // pair.first: Reference Comp_ID
       // pair.second: Instance vs. Reference atom mapping;
       std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > > atom_mapping;
       atom_mapping.clear();

       std::set<std::string> instanceid_set;
       instanceid_set.clear();

       bool status = true;
       if (!assignfile.empty()) status = _readInstanceAssignFile(assignfile, atom_mapping, true);
       if (!status) return false;

       bool no_entity_flag = false;
       bool keep_solvent_position_flag = false;
       if (_skip_task_set.find("solvent position") != _skip_task_set.end()) keep_solvent_position_flag = true;

       bool find_sugar_link_flag = true;
       if (_skip_task_set.find("link") != _skip_task_set.end()) find_sugar_link_flag = false;

       Entity entity;
       std::vector<std::string> seqs;

       for (std::vector<RCSB::Molecule*>::iterator pos = _molecules.begin(); pos != _molecules.end(); ++pos) {
            (*pos)->setCCDic(_ccDic);
            (*pos)->setMessage(&_error_messages);
            (*pos)->setLog(_logIo);
            (*pos)->set_Carbohydrate_Annotation_Info(_Branch_Seq_Scheme_Mapping);
            (*pos)->setRenameResidueFlag(_rename_residue_flag);
            (*pos)->setFindSugarChainFlag(_find_sugar_chain_flag);
            if ((_experiment_type == EXPERIMENT_TYPE_NMR) || (_experiment_type == EXPERIMENT_TYPE_NMR_SOLID)) (*pos)->changed_zero_occupancy();
            if (!(*pos)->find_residues(split_flag, atom_mapping, instanceid_set, !_DepUI_Flag)) status = false;
            if (!_merge_residue_list.empty()) (*pos)->MergeResidue(_merge_residue_list);
            (*pos)->find_chains(_Seq_Scheme_Mapping, _Seqs, _link_residue_set, keep_solvent_position_flag, find_sugar_link_flag);
            if (!skip_update_flag) (*pos)->update_cif_nomenclature();

            // get existing entity IDs
            RCSB::Chain* chain = (*pos)->GetFirstChain();
            while (chain) {
                 if (!chain->prev_entity_id().empty()) {
                      int entity_id = atoi(chain->prev_entity_id().c_str());
                      if (entity_id > _max_entity_id) _max_entity_id = entity_id;
                      _entity_ids.insert(std::make_pair(chain->entity_key(), entity_id));
                      std::map<int, Entity>::iterator ipos = _entities.find(entity_id);
                      if (ipos == _entities.end()) {
                           chain->get_seq(seqs);
                           entity.clear();
                           entity.setEntityID(chain->prev_entity_id());
                           entity.insertPDBChainID(chain->PDB_ChainID());
                           entity.setSeqs(seqs);
                           _entities.insert(std::make_pair(entity_id, entity));
                      } else ipos->second.insertPDBChainID(chain->PDB_ChainID());
                 } else no_entity_flag = true;
                 chain = (*pos)->GetNextChain();
            }
       }

       if (no_entity_flag) assign_prev_entity_id();

       return status;
}

void FileObj::build_molecule_only(const bool& using_coordinate_only)
{
       if (_molecules.empty()) return;

       bool keep_solvent_position_flag = false;
       if (_skip_task_set.find("solvent position") != _skip_task_set.end()) keep_solvent_position_flag = true;

       bool find_sugar_link_flag = true;
       if (_skip_task_set.find("link") != _skip_task_set.end()) find_sugar_link_flag = false;

       for (std::vector<RCSB::Molecule*>::iterator pos = _molecules.begin(); pos != _molecules.end(); ++pos) {
            (*pos)->setCCDic(_ccDic);
            (*pos)->setMessage(&_error_messages);
            (*pos)->setLog(_logIo);
            (*pos)->set_Carbohydrate_Annotation_Info(_Branch_Seq_Scheme_Mapping);
            (*pos)->find_residues(!_DepUI_Flag);
            if (using_coordinate_only)
                 (*pos)->find_chains(_link_residue_set, find_sugar_link_flag);
            else (*pos)->find_chains(_Seq_Scheme_Mapping, _Seqs, _link_residue_set, keep_solvent_position_flag, find_sugar_link_flag);
       }
}

void FileObj::assign_prev_entity_id()
{
       // Assign previous entity IDs

       Entity entity;
       std::vector<std::string> seqs;

       for (std::vector<RCSB::Molecule*>::iterator pos = _molecules.begin(); pos != _molecules.end(); ++pos) {
            RCSB::Chain* chain = (*pos)->GetFirstChain();
            while (chain) {
                 if (chain->prev_entity_id().empty()) {
                      std::map<std::string, int>::const_iterator
                          kpos = _entity_ids.find(chain->entity_key());
                      if (kpos == _entity_ids.end()) {
                           _max_entity_id++;
                           std::string entity_id = String::IntToString(_max_entity_id);
                           chain->set_prev_entity_id(entity_id);
                           chain->set_entity_id(entity_id);
 
                           _entity_ids.insert(std::make_pair(chain->entity_key(),
                                              _max_entity_id));

                           chain->get_seq(seqs);

                           entity.clear();
                           entity.setEntityID(entity_id);
                           entity.insertPDBChainID(chain->PDB_ChainID());
                           entity.setSeqs(seqs);
                           _entities.insert(std::make_pair(_max_entity_id, entity));
                      } else {
                           chain->set_prev_entity_id(String::IntToString(kpos->second));
                           chain->set_entity_id(String::IntToString(kpos->second));
                           std::map<int, Entity>::iterator
                               ipos = _entities.find(kpos->second);
                           if (ipos != _entities.end())
                                ipos->second.insertPDBChainID(chain->PDB_ChainID());
                      }
                 }
                 chain = (*pos)->GetNextChain();
            }
       }
}

void FileObj::print_data()
{
       fprintf(stdout, "_experiment_type=%d\n", _experiment_type);

       fprintf(stdout, "SEQRES:\n");
       for (std::map<std::string, SEQ>::const_iterator 
            mpos = _Seqs.begin(); mpos != _Seqs.end(); ++mpos) {
            fprintf(stdout, "Chain %s: %d\n", mpos->first.c_str(),
                            (int) mpos->second.res.size());
            int count = 0;
            for (std::vector<std::string>::const_iterator vpos =
                 mpos->second.res.begin(); vpos != mpos->second.res.end(); ++vpos) {
                 if (count == 20) {
                      fprintf(stdout, "\n");
                      count = 0;
                 }
                 fprintf(stdout, "%s ", vpos->c_str());
                 count++;
            }
            fprintf(stdout, "\n");
       }

       fprintf(stdout, "SCALE matrix:\n");
       for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) fprintf(stdout, " %.8f", _scale[i][j]);
            fprintf(stdout, " %.8f\n", _vec[i]);
       }

       fprintf(stdout, "CRYST1: %d %s %s %s %s %s %s %s\n", _cell.is_artifical(),
                       _cell.a().c_str(), _cell.b().c_str(), _cell.c().c_str(),
                       _cell.alpha().c_str(), _cell.beta().c_str(),
                       _cell.gamma().c_str(), _cell.space_group().c_str());

       PdbWrite writer;
       writer.setLog(_logIo);
       writer.setCCDic(_ccDic);
       writer.setRecord(&_pdb_records);
       writer.setMolecule(&_molecules);
       writer.setFileName(_output_filename);
       writer.setOutputFormat(NDB_FILE_FORMAT_PDB);
       writer.WriteRecord();
}
