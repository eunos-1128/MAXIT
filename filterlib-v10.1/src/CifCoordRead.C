/*
FILE:     CifCoordRead.C
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

#include "Atom.h"
#include "CifCoordRead.h"
#include "PdbRead.h"
#include "TypeDef.h"
#include "utillib.h"

#define EXPERIMENT_TYPE_XRAY        2
#define EXPERIMENT_TYPE_NEUTRON     4
#define EXPERIMENT_TYPE_FIBER       8
#define EXPERIMENT_TYPE_ELECTRON   16

#define NUM_MISSING_ATOM_SITE_VALUES    8
#define NUM_INCONSISTENT_ATOM_SITE_VALUES 2

static const char *_missing_atom_site_values[NUM_MISSING_ATOM_SITE_VALUES] = {
       "atom name,2,18",
       "residue name,4,16",
       // "chain ID,5,17",
       "residue number,6,15",
       "x-axis coordinate,8",
       "y-axis coordinate,9",
       "z-axis coordinate,10",
       "occupancy,11",
       "B factor,12"
};

#define NUM_MISSING_VALUE_MINIMUM       6
#define NUM_MISSING_VALUE_LABELS        9

static const char *_missing_value_labels[NUM_MISSING_VALUE_LABELS] = {
       "atom name", "residue name", "residue number", "x-axis coordinate", "y-axis coordinate", "z-axis coordinate", "occupancy",
       "B factor", "anisotropic temperature factors"
};

static const char *_atom_records[FORMAT_SIGATM] = {
       "Cartn_x", "Cartn_y", "Cartn_z", "occupancy", "B_iso_or_equiv"
}; 

CifCoordRead::CifCoordRead()
{
       _init();
}

CifCoordRead::~CifCoordRead()
{
       _init();
}

void CifCoordRead::_init()
{
       _logIo = NULL;
       _messageIo = NULL;
       _dictUtil = NULL;
       _pdb_records = NULL;
       _molecules = NULL;
       _cifobj = NULL;
       _experiment_type = 0;
       _format_checking = false;
       _atom_record_set.clear();
       _missing_atom_records.clear();
       _inconsistent_atom_records.clear();
}

void CifCoordRead::setLog(LogUtil *logPt)
{
       _logIo = logPt;
}

void CifCoordRead::setMessage(MessageUtil *message)
{
       _messageIo = message;
}

void CifCoordRead::setDictObj(DictUtil* dictobj)
{
       _dictUtil = dictobj;
}

void CifCoordRead::setRecord(std::map<std::string, std::list<std::vector<std::string> > >* records)
{
       _pdb_records = records;
}

void CifCoordRead::setMolecule(std::vector<RCSB::Molecule*>* mols)
{
       _molecules = mols;
}

void CifCoordRead::setCifObj(CifFile *obj)
{
       _cifobj = obj;
}

void CifCoordRead::setExperimentType(const int& type)
{
       _experiment_type = type;
}

void CifCoordRead::setFormatChecking()
{
       _format_checking = true;
       for (int i = 0; i < FORMAT_SIGATM; ++i) {
            _atom_record_set.insert(_atom_records[i]);
       }

       std::vector<std::string> str_data;
       std::vector<int> int_data;
       for (int i = 0; i < NUM_MISSING_ATOM_SITE_VALUES; ++i) {
            get_wordarray(str_data, _missing_atom_site_values[i], ",");
            int_data.clear();
            for (unsigned int j = 1; j < str_data.size(); ++j) {
                 int_data.push_back(atoi(str_data[j].c_str()));
            }
            _missing_atom_records.push_back(std::make_pair(str_data[0], int_data));
            if (i < NUM_INCONSISTENT_ATOM_SITE_VALUES) _inconsistent_atom_records.push_back(std::make_pair(str_data[0], int_data));
       }
}

void CifCoordRead::Read()
{
       _molecules->clear();

       std::string blockId = _cifobj->GetFirstBlockName();
       Block& block = _cifobj->GetBlock(blockId);

       if (getTablePtr(block, "atom_site"))
            _getCoordinatefrom_atom_site(block);
       else _getCoordinatefrom_encap_coordinates(block);

       if (_molecules->size() == 1) (*_molecules)[0]->set_Mol_ID(1);
}

const std::vector<std::string>& CifCoordRead::getExtraItems() const
{
       return _extra_items;
}

void CifCoordRead::_getCoordinatefrom_atom_site(Block& block)
{
       ISTable *t = getTablePtr(block, "atom_site");
       if (!t) return;

       std::set<std::string> allowed_item_set;

       allowed_item_set.clear();
       allowed_item_set.insert("pdbx_PDB_model_num");
       for (int i = 0; i <= FORMAT_LONG; ++i) {
            allowed_item_set.insert(_format_atom_site_pdbx[i]);
       }

       std::map<std::string, std::string> re_name_mapping;
       re_name_mapping.clear();
       re_name_mapping.insert(std::make_pair("pdbx_PDB_atom_name", "pdbx_auth_atom_name"));
       re_name_mapping.insert(std::make_pair("pdbx_PDB_residue_name", "pdbx_auth_comp_id"));
       re_name_mapping.insert(std::make_pair("pdbx_PDB_residue_no", "pdbx_auth_seq_id"));
       re_name_mapping.insert(std::make_pair("pdbx_PDB_strand_id", "pdbx_auth_asym_id"));

       std::string cs, atom_label;

       std::vector<std::string> itemNames = t->GetColumnNames();
       for (std::vector<std::string>::const_iterator vpos = itemNames.begin(); vpos != itemNames.end(); ++vpos) {
            if (allowed_item_set.find(*vpos) != allowed_item_set.end()) continue;
            std::map<std::string, std::string>::const_iterator mpos = re_name_mapping.find(*vpos);
            if (mpos != re_name_mapping.end()) {
                 if (!t->IsColumnPresent(mpos->second)) t->RenameColumn(mpos->first, mpos->second);
                 continue;
            }
            if ((_dictUtil && _dictUtil->Read() && _dictUtil->isDefinedItem("atom_site", *vpos, cs)) || ((*vpos) == "Q-score")) _extra_items.push_back(*vpos);
       }

       std::map<std::string, std::vector<std::string> > not_a_number_error;
       not_a_number_error.clear();

       std::set<std::string> not_allowed_chain_ids;
       not_allowed_chain_ids.clear();

       std::vector<std::string> data;
       std::map<std::string, int> _AtomNumbers;
       _AtomNumbers.clear();
       ISTable *anisou = getTablePtr(block, "atom_site_anisotrop");
       if (anisou) {
            int rowNo = anisou->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs, anisou, i, "id");
                 if (cs.empty()) continue;
                 _AtomNumbers.insert(std::make_pair(cs, i));

                 if (!_messageIo || !_format_checking) continue;

                 for (int j = 0; j < FORMAT_SIGUIJ; ++j) {
                      get_value_clean(cs, anisou, i, _format_anisotrop_pdbx[8 + j]);
                      if (!cs.empty() && !String::IsNumber(cs)) {
                           std::map<std::string, std::vector<std::string> >::iterator mpos = not_a_number_error.find("atom_site_anisotrop");
                           if (mpos != not_a_number_error.end()) mpos->second.push_back(cs);
                           else {
                                data.clear();
                                data.push_back(cs);
                                not_a_number_error.insert(std::make_pair("atom_site_anisotrop", data));
                           }
                      }
                      get_value_clean(cs, anisou, i, _format_siguij_pdbx[j]);
                      if (!cs.empty() && !String::IsNumber(cs)) {
                           std::map<std::string, std::vector<std::string> >::iterator mpos = not_a_number_error.find("atom_site_anisotrop");
                           if (mpos != not_a_number_error.end()) mpos->second.push_back(cs);
                           else {
                                data.clear();
                                data.push_back(cs);
                                not_a_number_error.insert(std::make_pair("atom_site_anisotrop", data));
                           }
                      }
                 }
            }
       }

       std::list<std::pair<std::pair<std::string, std::string>, std::set<std::string> > > missing_record_error;
       missing_record_error.clear();
       std::set<std::string> t_set;

       std::list<std::pair<std::pair<std::string, std::string>, std::vector<std::vector<std::string> > > > inconsistent_record_error;
       inconsistent_record_error.clear();
       std::vector<std::vector<std::string> > t_vec;
       

       std::vector<std::vector<int> > atom_label_index;
       atom_label_index.clear();
       if (_messageIo && _format_checking) {
            std::vector<int> v_int;
            v_int.clear();
            v_int.push_back(2);
            v_int.push_back(18);
            atom_label_index.push_back(v_int);
            v_int.clear();
            v_int.push_back(4);
            v_int.push_back(16);
            atom_label_index.push_back(v_int);
            v_int.clear();
            v_int.push_back(5);
            v_int.push_back(17);
            atom_label_index.push_back(v_int);
            v_int.clear();
            v_int.push_back(6);
            v_int.push_back(15);
            atom_label_index.push_back(v_int);
       }

       std::vector<std::string> record;
       record.clear();
       record.reserve(FORMAT_LONG + 1);

       bool has_auth_asym_id = false;
       if (t->IsColumnPresent("auth_asym_id")) has_auth_asym_id = true;
       bool has_pdbx_auth_asym_id = false;
       if (t->IsColumnPresent("pdbx_auth_asym_id")) has_pdbx_auth_asym_id = true;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            record.clear();
            for (int j = 0; j <= FORMAT_LONG; ++j) {
                 get_value_clean(cs, t, i, _format_atom_site_pdbx[j]);
                 if (_messageIo && _format_checking && ((j == 5) || (j == 17)) && !IsAlnum(cs)) {
                      not_allowed_chain_ids.insert(cs);
                 }
                 if (_messageIo && _format_checking && _atom_record_set.find(_format_atom_site_pdbx[j]) != _atom_record_set.end() &&
                    !cs.empty() && !String::IsNumber(cs)) {
                      std::map<std::string, std::vector<std::string> >::iterator mpos = not_a_number_error.find("atom_site");
                      if (mpos != not_a_number_error.end()) mpos->second.push_back(cs);
                      else {
                           data.clear();
                           data.push_back(cs);
                           not_a_number_error.insert(std::make_pair("atom_site", data));
                      }
                 }
                 record.push_back(cs);
            }
            if (record[17].empty() && has_auth_asym_id) record[17] = "?";
            if (record[21].empty() && has_pdbx_auth_asym_id) record[21] = "?";
            RCSB::Atom* atom = new RCSB::Atom;
            atom->setValue(record, i + 1);

            for (std::vector<std::string>::const_iterator vpos = _extra_items.begin(); vpos != _extra_items.end(); ++vpos) {
                 get_value_clean(cs, t, i, *vpos);
                 if (!cs.empty()) atom->setExtraValue(*vpos, cs);
            }

            std::string id = record[1];
            if (_messageIo && _format_checking) {
                 atom_label.clear();
                 for (std::vector<std::vector<int> >::const_iterator iipos = atom_label_index.begin(); iipos != atom_label_index.end(); ++iipos) {
                      for (std::vector<int>::const_iterator ipos = iipos->begin(); ipos != iipos->end(); ++ipos) {
                           if (record[*ipos].empty() || record[*ipos] == "?") continue;
                           if (!atom_label.empty()) atom_label += " ";
                           atom_label += record[*ipos];
                           break;
                      }
                 }
            }

            t_set.clear();
            t_vec.clear();
            if (_messageIo && _format_checking) {
                 for (std::vector<std::pair<std::string, std::vector<int> > >::const_iterator vpos = _missing_atom_records.begin();
                      vpos != _missing_atom_records.end(); ++vpos) {
                      bool found_value = false;
                      for (std::vector<int>::const_iterator ipos = vpos->second.begin(); ipos != vpos->second.end(); ++ipos) {
                           if (record[*ipos].empty() || record[*ipos] == "?") continue;
                           found_value = true;
                           break;
                      }
                      if (!found_value) t_set.insert(vpos->first);
                 }
                 for (std::vector<std::pair<std::string, std::vector<int> > >::const_iterator vpos = _inconsistent_atom_records.begin();
                      vpos != _inconsistent_atom_records.end(); ++vpos) {
                      if (record[vpos->second[0]].empty() || record[vpos->second[0]] == "?" || record[vpos->second[1]].empty() ||
                          record[vpos->second[1]] == "?" || record[vpos->second[0]] == record[vpos->second[1]]) continue;
                      data.clear();
                      data.push_back(vpos->first);
                      data.push_back(record[vpos->second[0]]);
                      data.push_back(_format_atom_site_pdbx[vpos->second[0]]);
                      data.push_back(record[vpos->second[1]]);
                      data.push_back(_format_atom_site_pdbx[vpos->second[1]]);
                      t_vec.push_back(data);
                 }
            }

            bool has_value = false;
            record.clear();
            record.push_back("SIGATM");
            for (int j = 0; j < FORMAT_SIGATM; ++j) {
                 get_value_clean(cs, t, i, _format_sigatm_pdbx[j]);
                 if (!cs.empty()) {
                      has_value = true;
                      if (_messageIo && _format_checking && _atom_record_set.find(_format_atom_site_pdbx[j]) != _atom_record_set.end()) {
                           std::map<std::string, std::vector<std::string> >::iterator mpos = not_a_number_error.find("atom_site");
                           if (mpos != not_a_number_error.end()) mpos->second.push_back(cs);
                           else {
                                data.clear();
                                data.push_back(cs);
                                not_a_number_error.insert(std::make_pair("atom_site", data));
                           }
                      }
                 }
                 record.push_back(cs);
            }
            if (has_value) atom->setAuxiliaryValue(record);

            std::map<std::string, int>::const_iterator mpos = _AtomNumbers.find(id);
            if (mpos != _AtomNumbers.end()) {
                 record.clear();
                 record.push_back("ANISOU");
                 has_value = false;
                 for (int j = 0; j < FORMAT_SIGUIJ; ++j) {
                      get_value_clean(cs, anisou, mpos->second, _format_anisotrop_pdbx[8 + j]);
                      if (!cs.empty()) has_value = true;
                      record.push_back(cs);
                 }
                 if (has_value) {
                      if (_messageIo && _format_checking) {
                           for (int j = 0; j < FORMAT_SIGUIJ; ++j) {
                                if (record[j + 1].empty()) t_set.insert("anisotropic temperature factors");
                           }
                      }
                      atom->setAuxiliaryValue(record);
                 }

                 record.clear();
                 record.push_back("SIGUIJ");
                 has_value = false;
                 for (int j = 0; j < FORMAT_SIGUIJ; ++j) {
                      get_value_clean(cs, anisou, mpos->second, _format_siguij_pdbx[j]);
                      if (!cs.empty()) has_value = true;
                      record.push_back(cs);
                 }
                 if (has_value) atom->setAuxiliaryValue(record);
            }
            if (!t_set.empty()) missing_record_error.push_back(std::make_pair(std::make_pair(id, atom_label), t_set));
            if (!t_vec.empty()) inconsistent_record_error.push_back(std::make_pair(std::make_pair(id, atom_label), t_vec));

            get_value_clean(cs, t, i, "pdbx_PDB_model_num");
            if (cs.empty()) cs = "1";
            int model_id = atoi(cs.c_str());
            RCSB::Molecule* mol = NULL;
            for (std::vector<RCSB::Molecule*>::iterator pos = _molecules->begin(); pos != _molecules->end(); ++pos) {
                 if ((*pos)->Mol_ID() == model_id) {
                      mol = *pos;
                      break;
                 }
            }
            if (!mol) {
                 mol = new RCSB::Molecule;
                 mol->set_Mol_ID(model_id);
                 mol->set_index(_molecules->size());
                 _molecules->push_back(mol);
            }
            mol->insert_a_atom(atom);
       }

       if (!not_a_number_error.empty()) {
            for (std::map<std::string, std::vector<std::string> >::const_iterator mpos = not_a_number_error.begin();
                 mpos != not_a_number_error.end(); ++mpos) {
                 std::string error = "In '" + mpos->first + "' category, '" + join_string(mpos->second, " ") + "'";
                 if (mpos->second.size() > 1)
                      error += " are not numbers.";
                 else error += " is not a number.";
                 _messageIo->insertMessage("error", "model", error, true, "", false);
            }
       }

       if (!missing_record_error.empty()) {
            int missing_value_number = NUM_MISSING_VALUE_MINIMUM;
            if ((_experiment_type & EXPERIMENT_TYPE_XRAY) || (_experiment_type & EXPERIMENT_TYPE_NEUTRON) || (_experiment_type & EXPERIMENT_TYPE_FIBER) ||
                (_experiment_type & EXPERIMENT_TYPE_ELECTRON)) missing_value_number = NUM_MISSING_VALUE_LABELS;
            for (std::list<std::pair<std::pair<std::string, std::string>, std::set<std::string> > >::const_iterator
                 lpos = missing_record_error.begin(); lpos != missing_record_error.end(); ++lpos) {
                 cs.clear();
                 int count = 0;
                 for (int i = 0; i < missing_value_number; ++i) {
                      if (lpos->second.find(_missing_value_labels[i]) == lpos->second.end()) continue;
                      if (!cs.empty()) cs += ", ";
                      cs += _missing_value_labels[i];
                      count++;
                 }
                 if (cs.empty()) continue;
                 std::string error = "Missing ( " + cs + " ) ";
                 if (count > 1) error += "values ";
                 else error += "value ";
                 if (!lpos->first.second.empty()) error += "for atom record ( " + lpos->first.second + " ) ";
                 error += " at id row " + lpos->first.first + ".";
                 _messageIo->insertMessage("error", "model", error, true, "", false);
            }
       }

       if (!inconsistent_record_error.empty()) {
            for (std::list<std::pair<std::pair<std::string, std::string>, std::vector<std::vector<std::string> > > >::const_iterator
                 lpos = inconsistent_record_error.begin(); lpos != inconsistent_record_error.end(); ++lpos) {
                 for (std::vector<std::vector<std::string> >::const_iterator vpos = lpos->second.begin(); vpos != lpos->second.end(); ++vpos) {
                      cs = "For atom record ( " + lpos->first.second + " ) at id row " + lpos->first.first + ", " + (*vpos)[0] + " '" + (*vpos)[1]
                         + "' in '_atom_site." + (*vpos)[2] + "' field is inconsistent with '" + (*vpos)[3] + "' in '_atom_site." + (*vpos)[4]
                         + "' field.";
                      _messageIo->insertMessage("error", "model", cs, true, "", false);
                 }
            }
       }

       if (!not_allowed_chain_ids.empty()) {
            std::string value = "";
            for (std::set<std::string>::const_iterator spos = not_allowed_chain_ids.begin(); spos != not_allowed_chain_ids.end(); ++spos) {
                 if (!value.empty()) value += ", ";
                 value += "'" + *spos + "'";
            }
            std::string error = "For '_atom_site.label_asym_id' and '_atom_site.auth_asym_id' fields, only alphanumeric values (A-Z, 0-9, a-z) are allowed.";
            if (not_allowed_chain_ids.size() > 1)
                 error += " The values ( " + value + " ) are not allowed.";
            else error += " The value " + value + " is not allowed.";
            _messageIo->insertMessage("error", "model", error, true, "", false);
       }
}

void CifCoordRead::_getCoordinatefrom_encap_coordinates(Block& block)
{
       std::string cs = "";
       int input_format = NDB_FILE_FORMAT_PDB;
       if (block.IsTablePresent("ndb_original_pdb_coordinates")) {
            input_format = NDB_FILE_FORMAT_PDB;
            cs = "ndb_original_pdb_coordinates";
       } else if (block.IsTablePresent("ndb_original_ndb_coordinates")) {
            input_format = NDB_FILE_FORMAT_MIX;
            cs = "ndb_original_ndb_coordinates";
       } else if (block.IsTablePresent("pdbx_original_pdb_coordinates")) {
            input_format = NDB_FILE_FORMAT_PDB;
            cs = "pdbx_original_pdb_coordinates";
       } else if (block.IsTablePresent("pdb_original_pdb_coordinates")) {
            input_format = NDB_FILE_FORMAT_PDB;
            cs = "pdb_original_pdb_coordinates";
       }

       if (cs.empty()) return;

       ISTable *t = getTablePtr(block, cs);
       if (!t) return;

       const std::vector<std::string>& columnname = t->GetColumnNames();
       get_value(cs, t, 0, columnname[0]);
       if (cs.empty()) return;

       PdbRead reader;
       reader.setLog(_logIo);
       reader.setMessage(_messageIo);
       if (_format_checking) reader.setFormatChecking();
       reader.setInputFormat(input_format);
       reader.setRecord(_pdb_records);
       reader.setMolecule(_molecules);
       reader.ReadRecord(cs);
}
