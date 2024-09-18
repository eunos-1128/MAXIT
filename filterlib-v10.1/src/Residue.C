/*
FILE:     Residue.C
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

#include <math.h>

#include "BondUtil.h"
#include "CompositeIndex.h"
#include "Element.h"
#include "Exceptions.h"
#include "GraphMatch.h"
#include "MatrixUtil.h"
#include "Residue.h"
#include "Residue_global.h"
#include "SeqCodeUtil.h"
#include "SideChainPattern.h"
#include "TypeDef.h"
#include "utillib.h"

using namespace RCSB;

IntegerStringOrder::IntegerStringOrder()
{
       _int_value = 0;
       _str_value.clear();
}

const int IntegerStringOrder::int_value() const { return _int_value; }
const std::string &IntegerStringOrder::str_value() const { return _str_value; }
void IntegerStringOrder::set_int_value(const int& a)   { _int_value = a; }
void IntegerStringOrder::set_str_value(const std::string &a) { _str_value = a; }

IntegerStringOrder& IntegerStringOrder::operator=(const IntegerStringOrder &index)
{
       _int_value = index._int_value;
       _str_value = index._str_value;
       return (*this);
}

bool IntegerStringOrder::operator()(const IntegerStringOrder &index1,
                                    const IntegerStringOrder &index2) const
{
       return (index1 < index2);
}

bool IntegerStringOrder::operator<(const IntegerStringOrder &index) const
{
       if (_int_value < index._int_value)
            return true;
       else if (_int_value == index._int_value) {
            if (_str_value < index._str_value)
                 return true;
            else return false;
       } else return false;
}

bool IntegerStringOrder::operator==(const IntegerStringOrder &index) const
{
       if (_int_value == index._int_value &&
           _str_value == index._str_value)
            return true;
       else return false;
}

Residue::Residue()
{
       _ccDic = NULL;
       _logIo = NULL;
       _messageIo = NULL;
       clear(); 
}

void Residue::clear()
{
       _token.clear();
       _ResName.clear();
       _OrigResName.clear();
       _chnid.clear();
       _res_no.clear();
       _pdb_chnid.clear();
       _pdb_res_no.clear();
       _ins_code.clear();
       _alt_loc.clear();
       _alt_loc_list.clear();
       _ter_flag = 0;
       _atoms.clear();
       _type = MONOMER_TYPE_UNK;
       _chain_type.clear();
       _entity_id.clear();
       _chain_break = 0;
       _wrong_connectivity = false;
       _disorder_flag = false;
       _full_disorder_flag = false;
       _index = 0;
       _tls_group = 0;
       _chain_index = 0;
       _position = 0;
       _uniqueAtomNames.clear();
       _atomNames.clear();
       _atomNameAlts.clear();
       _atomNameAlts_Orininal.clear();
       _atomOrder.clear();
       _missing.clear();
       _extras.clear();
       _atom_pos = _atomOrder.begin();
}

void Residue::_reset()
{
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (*pos) delete *pos;
       }
       clear();
}

bool Residue::_UpdateIndices()
{
       _uniqueAtomNames.clear();
       _atomNames.clear();
       _atomNameAlts.clear();
       int i = 0;
       bool has_hydrogen = false;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if ((*pos)->atom_type() == "H" || (*pos)->atom_type() == "D") has_hydrogen = true;
            _uniqueAtomNames.insert((*pos)->getAtomIndex());
            _atomNames.insert(std::make_pair((*pos)->getAtomIndex(), i));
            _atomNameAlts.insert(std::make_pair((*pos)->getAtomAltIndex(), i));
            i++;
       }
       return has_hydrogen;
}

void Residue::_writeAtomRecord(RCSB::Atom* atom, std::string& record)
{
       record.clear();

       std::string chn_id = atom->pdb_chnid();
       if (chn_id.empty()) chn_id = " ";
       std::string inscode = atom->ins_code();
       if (inscode.empty()) inscode = " ";
       std::string alt_loc = atom->alt_loc();
       if (alt_loc.empty()) alt_loc = " "; 

       char buffer[100];

       sprintf(buffer, FORMAT, atom->atnum().c_str(), atom->pdb_atmnam().c_str(),
            alt_loc[0], atom->pdb_resnam().c_str(), chn_id.c_str(),
            atom->pdb_resnum().c_str(), inscode[0], atom->x(), atom->y(),
            atom->z(), atof(atom->occ().c_str()), atof(atom->t_fct().c_str()));
       record += buffer;
}

void Residue::_writeErrorMessage(RCSB::Atom* atom1, RCSB::Atom* atom2, const std::string& type, const std::string& message, const int& Mol_ID)
{
       std::string cs, buffer;
       cs.clear();
       if (Mol_ID < 0)
            cs += "\n";
       else {
            cs += "\nIn model " + String::IntToString(Mol_ID) + ", ";
       }
       cs += message;

       _writeAtomRecord(atom1, buffer); cs += buffer;
       _writeAtomRecord(atom2, buffer); cs += buffer;

       _messageIo->insertMessage(type, "label", "true");
       _messageIo->insertMessage(type + "_print", "error", cs);
}

void Residue::_writeOccupancyMessage(const std::string &name, const float& occupancy, const int& Mol_ID)
{
       _messageIo->insertMessage("occupancy", "label", "true");

       std::string error = "";
       if (Mol_ID >= 0) error += "In model " + String::IntToString(Mol_ID) + ",";

       error += "Atom '" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + " " + name;
       error += "' has occupancy " + FloatToString(occupancy, 0, 2) + ".";
       _messageIo->insertMessage("occupancy_print", "warning", error);
}

void Residue::_writeOccupancyError(const std::string &name, const float& occupancy, const int& Mol_ID)
{
       std::string error = "";
       if (Mol_ID >= 0) error += "In model " + String::IntToString(Mol_ID) + ",";

       error += "Atom '" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + " " + name;
       error += "' has occupancy " + FloatToString(occupancy, 0, 2) + ".\n";
       _logIo->message(error.c_str());
}

const std::string& Residue::token() const { return _token; }
const std::string& Residue::ResName() const { return _ResName; }
const std::string& Residue::OrigResName() const { return _OrigResName; }
const std::string& Residue::chnid() const { return _chnid; }
const std::string& Residue::res_no() const { return _res_no; }
const std::string& Residue::pdb_chnid() const { return _pdb_chnid; }
const std::string& Residue::pdb_res_no() const { return _pdb_res_no; }
const std::string& Residue::ins_code() const { return _ins_code; }
const std::string& Residue::alt_loc() const { return _alt_loc; }
const std::vector<RCSB::Atom*>& Residue::atoms() const { return _atoms; }
// const bool Residue::has_alt_loc() const { return _alt_loc_list.size() > 0; }
const std::vector<std::string>& Residue::alt_loc_list() const { return _alt_loc_list; }
const bool& Residue::has_alt_loc() const { return _full_disorder_flag; }
const int& Residue::ter_flag() const { return _ter_flag; }
const int& Residue::type() const { return _type; }
const std::string& Residue::chain_type() const { return _chain_type; }
const std::string& Residue::entity_id() const { return _entity_id; }
const int& Residue::chain_break() const { return _chain_break; }
const bool& Residue::wrong_connectivity() const { return _wrong_connectivity; }
const bool& Residue::disorder_flag() const { return _disorder_flag; }
const int& Residue::index() const { return _index; }
const int& Residue::tls_group() const { return _tls_group; }
const int& Residue::chain_index() const { return _chain_index; }
const int& Residue::position() const { return _position; }
const std::vector<std::string>& Residue::missing() const { return _missing; }
const std::vector<std::string>& Residue::extras() const { return _extras; }
const bool Residue::is_single_atom() const { return (_uniqueAtomNames.size() == 1); }
const unsigned int Residue::AtomNumbers() const { return _atomNames.size(); }
const unsigned int Residue::NumAtoms() const { return _atoms.size(); }

const int Residue::NonHydrogenAtomNumbers()
{
       std::set<std::string> index;
       index.clear();
       for (std::vector<RCSB::Atom*>::const_iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if ((*pos)->atom_type() == "H" || (*pos)->atom_type() == "D") continue;
            index.insert((*pos)->atmtype());
       }
       return (index.size());
}

const int Residue::ca_or_p_atom_only() const
{
       std::set<std::string> name_set;
       name_set.clear();
       for (std::vector<RCSB::Atom*>::const_iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            name_set.insert((*pos)->pdb_atmnam());
            if (name_set.size() > 1) return 0;
       }

       if (name_set.size() == 1) {
            std::set<std::string>::const_iterator spos = name_set.begin();
            if (*spos == "CA") return CAAtom_ONLY;
            if (*spos == "P") return PAtom_ONLY;
       }

       return 0;
       
}

void Residue::setCCDic(ConnectDic *ccdic) { _ccDic = ccdic; }
void Residue::setLog(LogUtil *logPt) { _logIo = logPt; }
void Residue::setMessage(MessageUtil *message) { _messageIo = message; }
void Residue::UpdateIndices() { _UpdateIndices(); }
void Residue::set_token(const std::string &a) { _token = a; }
void Residue::set_ResName(const std::string &a) { _ResName = a; }
void Residue::set_chnid(const std::string &a) { _chnid = a; }
void Residue::set_res_no(const std::string& a) { _res_no = a; }
void Residue::set_pdb_chnid(const std::string& a) { _pdb_chnid = a; }
void Residue::set_pdb_res_no(const std::string& a) { _pdb_res_no = a; }
void Residue::set_ins_code(const std::string& a) { _ins_code = a; }
void Residue::set_type(const int& a) { _type = a; }
void Residue::set_chain_type(const std::string& a) { _chain_type = a; }
void Residue::set_entity_id(const std::string& a) { _entity_id = a; }
void Residue::set_chain_break(const int& a) { _chain_break = a; }
void Residue::set_index(const int& a) { _index = a; }
void Residue::set_tls_group(const int& a) { _tls_group = a; }
void Residue::set_chain_index(const int& a) { _chain_index = a; }
void Residue::set_position(const int& a) { _position = a; }

std::string Residue::get_side_chain_pattern()
{
       if (_ResName == "ALA" || _ResName == "GLY") return "ala-gly-like";
       if (!SeqCodeUtil::is_standard_aa_residue_plus_MSE(_ResName) && (_ResName != "UNK")) return "";

       std::set<std::string> atom_name_set;
       atom_name_set.clear();

       for (std::vector<RCSB::Atom*>::const_iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if ((*pos)->is_hydrogen()) continue;
            atom_name_set.insert((*pos)->pdb_atmnam());
       }

       return SideChainPattern::get_pattern(_ResName, atom_name_set); 
}

void Residue::set_alt_loc()
{
       _alt_loc.clear();
       _alt_loc_list.clear();
       _disorder_flag = false;
       _full_disorder_flag = false;

       std::set<std::string> t_set;
       t_set.clear();
       bool has_no_alt_loc = false;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if ((*pos)->alt_loc().empty()) {
                 has_no_alt_loc = true;
                 continue;
            }
            _disorder_flag = true;
            if (t_set.find((*pos)->alt_loc()) != t_set.end()) continue;
            t_set.insert((*pos)->alt_loc());
            _alt_loc_list.push_back((*pos)->alt_loc());
       }
       // if (has_no_alt_loc) _alt_loc_list.clear();
       if (_disorder_flag && !has_no_alt_loc) _full_disorder_flag = true;
       if ((_alt_loc_list.size() == 1) && _full_disorder_flag) _alt_loc = _alt_loc_list[0];
}

RCSB::Atom* Residue::atoms(const int& pos)
{
       if (pos >= 0 && pos < (int) _atoms.size())
            return _atoms[pos];
       else return NULL;
}

RCSB::Atom* Residue::GetFirstAtom()
{
       _atom_pos = _atomOrder.begin();
       if (_atom_pos != _atomOrder.end()) {
            return _atoms[_atom_pos->second];
       } else return NULL;
}

RCSB::Atom* Residue::GetNextAtom()
{
       if (_atom_pos == _atomOrder.end()) return NULL;

       ++_atom_pos;
       if (_atom_pos != _atomOrder.end())
            return _atoms[_atom_pos->second];
       else return NULL;
}

void Residue::GetAtomNameList(std::list<std::string>& name_list)
{
       name_list.clear();

       std::set<std::string> name_set;
       name_set.clear();
       for (std::multimap<IntegerStringOrder, int>::const_iterator pos = _atomOrder.begin(); pos != _atomOrder.end(); ++pos) {
            if (name_set.find(_atoms[pos->second]->pdb_atmnam()) != name_set.end()) continue;
            name_set.insert(_atoms[pos->second]->pdb_atmnam());
            name_list.push_back(_atoms[pos->second]->pdb_atmnam());
       }
}

void Residue::GetAllBonds(std::vector<std::vector<RCSB::Atom*> >& bond_list)
{
       bond_list.clear();

       std::vector<RCSB::Atom*> bond_pair;
       for (unsigned int i = 0; i < _atoms.size() - 1; i++) {
            for (unsigned int j = i + 1; j < _atoms.size(); j++) {
                 if ((_atoms[i]->alt_loc().empty() || _atoms[j]->alt_loc().empty() || _atoms[i]->alt_loc() == _atoms[j]->alt_loc()) &&
                      BondUtil::is_a_bond(_atoms[i], _atoms[j])) {
                      bond_pair.clear();
                      bond_pair.push_back(_atoms[i]);
                      bond_pair.push_back(_atoms[j]);
                      bond_list.push_back(bond_pair);
                 }
            }
       }
}

void Residue::insert_a_atom(RCSB::Atom* atom, const int& Mol_ID)
{
       std::multimap<std::string, int>::iterator pos = _atomNameAlts.find(atom->getAtomAltIndex());
       if (pos != _atomNameAlts.end()) {
            _writeErrorMessage(atom, _atoms[pos->second], "alt_code", "atoms with incorrect alternate conformation ID:\n", Mol_ID);
       }
       _atomNameAlts.insert(std::make_pair(atom->getAtomAltIndex(), _atoms.size()));

       pos = _atomNames.find(atom->getAtomIndex());
       if (pos != _atomNames.end()) {
            if ((!atom->alt_loc().empty() && _atoms[pos->second]->alt_loc().empty()) || (atom->alt_loc().empty() && !_atoms[pos->second]->alt_loc().empty()))
                 _writeErrorMessage(atom, _atoms[pos->second], "alt_code", "atoms with incorrect alternate conformation ID:\n", Mol_ID);
       }
       _uniqueAtomNames.insert(atom->getAtomAltIndex());
       _atomNames.insert(std::make_pair(atom->getAtomIndex(), _atoms.size()));

       IntegerStringOrder idx;
       idx.set_int_value(_atoms.size());
       idx.set_str_value(atom->alt_loc());
       _atomOrder.insert(std::make_pair(idx, _atoms.size()));

       if (atom->ter_flag()) _ter_flag = atom->ter_flag();
       _atoms.push_back(atom);
}

void Residue::delete_a_atom(const std::string& name, std::set<long>& atom_set)
{
       std::string cs = CompositeIndex::getIndex(_ResName, name);
       std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNames.equal_range(cs);
       if (range.first == range.second) return;

       std::set<int> index_set;
       index_set.clear();
       for (std::multimap<std::string, int>::iterator pos = range.first; pos != range.second; ++pos) {
            index_set.insert(pos->second);
       }

       _delete_atoms(index_set, atom_set);

       _atomNameAlts_Orininal.clear();
       _UpdateIndices();

       IntegerStringOrder idx;
       _atomOrder.clear();
       int i = 0;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            idx.set_int_value(i);
            idx.set_str_value((*pos)->alt_loc());
            _atomOrder.insert(std::make_pair(idx, i));
            i++;
       }

       _atom_pos = _atomOrder.begin();
}

void Residue::remove_a_atom(const std::string& name, std::vector<RCSB::Atom*>& atom_list)
{
       atom_list.clear();
       std::string cs = CompositeIndex::getIndex(_ResName, name);
       std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNames.equal_range(cs);
       if (range.first == range.second) return;

       bool is_OXT_atom = false;
       bool has_O_atom = false;

       if (name == "OXT") {
            try {
                 const ConnectFormat& drug = _ccDic->find_drug(_ResName);
                 std::vector<std::string> atom_names;
                 atom_names.clear();
                 atom_names.push_back("OXT");
                 atom_names.push_back("O");
                 for (std::vector<std::string>::const_iterator vpos = atom_names.begin(); vpos != atom_names.end(); ++vpos) {
                      try {
                           const AtomFormat& atom = drug.find_atom(*vpos);
                           if ((*vpos == "OXT") && (atom.leaving_atom_flag() != "Y")) return;
                      } catch (const std::exception& exc) { return; }
                 }
            } catch (const std::exception& exc) { return; }

            is_OXT_atom = true;

            cs = CompositeIndex::getIndex(_ResName, "O");
            std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range1 = _atomNames.equal_range(cs);
            if (range1.first != range1.second) has_O_atom = true;
       }

       std::set<int> index_set;
       index_set.clear();
       for (std::multimap<std::string, int>::iterator pos = range.first; pos != range.second; ++pos) {
            index_set.insert(pos->second);
       }

       if (is_OXT_atom && !has_O_atom) {
            int i = 0;
            for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
                 if (index_set.find(i) != index_set.end()) {
                      (*pos)->set_atmtype("O");
                      (*pos)->set_pdb_atmnam("O");
                 }
                 i++;
            }
       } else {
            std::vector<RCSB::Atom*> atoms;
            atoms.clear();
            int i = 0;
            for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
                 if (index_set.find(i) != index_set.end()) atom_list.push_back(*pos);
                 else atoms.push_back(*pos);
                 i++;
            }
            _atoms = atoms;
       }

       _atomNameAlts_Orininal.clear();
       _UpdateIndices();

       IntegerStringOrder idx;
       _atomOrder.clear();
       int i = 0;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            idx.set_int_value(i);
            idx.set_str_value((*pos)->alt_loc());
            _atomOrder.insert(std::make_pair(idx, i));
            i++;
       }

       _atom_pos = _atomOrder.begin();
}

void Residue::find_atom(const std::string &name, std::vector<RCSB::Atom*>& atom_list)
{
       atom_list.clear();
       std::string cs = CompositeIndex::getIndex(_ResName, name);
       std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNames.equal_range(cs);
       if (range.first == range.second) return;
       for (std::multimap<std::string, int>::iterator mpos = range.first; mpos != range.second; ++mpos) {
            atom_list.push_back(_atoms[mpos->second]);
       }
}

void Residue::find_atom_by_type(const std::string &type, std::vector<RCSB::Atom*>& atom_list)
{
       atom_list.clear();
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if ((*pos)->atom_type() != type) continue;
            atom_list.push_back((*pos));
       }
}

RCSB::Atom* Residue::find_atom(const std::string &name)
{
       std::string cs = CompositeIndex::getIndex(_ResName, name);
       std::multimap<std::string, int>::iterator pos = _atomNames.find(cs);
       if (pos != _atomNames.end()) return _atoms[pos->second];
       else return NULL;
}

RCSB::Atom* Residue::find_atom(const std::string &name, const std::string &alt_loc, const bool& origin_flag)
{
       if (!_atomNameAlts_Orininal.empty() && origin_flag) {
             std::string resName = _ResName;
             if (!_OrigResName.empty()) resName = _OrigResName;
             std::string cs = CompositeIndex::getIndex(resName, name, alt_loc);
             std::multimap<std::string, int>::iterator pos = _atomNameAlts_Orininal.find(cs);
             if (pos != _atomNameAlts_Orininal.end()) return _atoms[pos->second];
       }

       std::string cs = CompositeIndex::getIndex(_ResName, name, "");
       std::multimap<std::string, int>::iterator pos = _atomNameAlts.find(cs);
       if (pos != _atomNameAlts.end()) return _atoms[pos->second];

       if (!alt_loc.empty()) {
            cs = CompositeIndex::getIndex(_ResName, name, alt_loc);
            pos = _atomNameAlts.find(cs);
            if (pos != _atomNameAlts.end()) return _atoms[pos->second];
       }
       return NULL;
}

void Residue::get_atom_type()
{
       if (_atoms.empty()) return;

       bool need_update = false;
       for (std::vector<RCSB::Atom*>::const_iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if ((*pos)->atom_type().empty()) {
                 need_update = true;
                 break;
            }
       }
       if (!need_update) return;

       try {
            const ConnectFormat& drug = _ccDic->find_drug(_ResName);
            for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
                 try {
                      const AtomFormat& atom = drug.find_atom((*pos)->pdb_atmnam());
                      (*pos)->set_atom_type(atom.atomtype());
                 } catch (const std::exception& exc) {}
            }
       } catch (const std::exception& exc) {}

       std::map<std::string, std::string> name_type;
       name_type.clear();
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (!(*pos)->atom_type().empty() &&
                 (*pos)->atom_type() != ".")
                 continue;
            else if ((*pos)->pdb_atmnam() == "1H" || (*pos)->pdb_atmnam() == "2H" || (*pos)->pdb_atmnam() == "3H")
                 (*pos)->set_atom_type("H");
            else if ((*pos)->pdb_atmnam() == "1D" || (*pos)->pdb_atmnam() == "2D" || (*pos)->pdb_atmnam() == "3D")
                 (*pos)->set_atom_type("D");
            else name_type.insert(std::make_pair((*pos)->pdb_atmnam(), ""));
       }

       if (name_type.empty()) return;

       Element::getAtomSymbols(name_type);
       bool is_aa = false;
       if (_ccDic->find_residue_type(_ResName) == "ATOMP") is_aa = true;

       std::vector<unsigned int> check_list;
       check_list.clear();
       for (unsigned int i = 0; i < _atoms.size(); ++i) {
            if (!_atoms[i]->atom_type().empty() && _atoms[i]->atom_type() != ".") continue;
            std::map<std::string, std::string>::iterator mpos = name_type.find(_atoms[i]->pdb_atmnam());
            if (mpos == name_type.end()) continue;
            if ((is_aa && (mpos->second == "HE" || mpos->second == "HF" || mpos->second == "HG")) ||
                 mpos->second == "HO") mpos->second = "H";
            _atoms[i]->set_atom_type(mpos->second);
            if (mpos->second == "HE" || mpos->second == "HF" || mpos->second == "HG") check_list.push_back(i);
       }

       if (check_list.empty()) return;

       for (std::vector<unsigned int>::const_iterator pos = check_list.begin(); pos != check_list.end(); ++pos) {
            double dist = 10000.0;
            for (unsigned int i = 0; i < _atoms.size(); ++i) {
                 if (_atoms[i]->atom_type() == "H" || _atoms[i]->atom_type() == "D" || _atoms[i]->atom_type() == "HE" ||
                     _atoms[i]->atom_type() == "HF" || _atoms[i]->atom_type() == "HG" || *pos == i) continue;

                 if (!_atoms[*pos]->alt_loc().empty() && !_atoms[i]->alt_loc().empty() &&
                     _atoms[*pos]->alt_loc() != _atoms[i]->alt_loc()) continue;
                 double distance = cal_distance(_atoms[*pos], _atoms[i]);
                 if (distance < dist) dist = distance;
            }
            if (dist < 1.6501) _atoms[*pos]->set_atom_type("H");
       }
}

void Residue::reorder_atoms()
{
       if (_atoms.empty()) return;
       if (_ResName == "UNL") return;

       try {
            const ConnectFormat& drug = _ccDic->find_drug(_ResName);
            std::string residueType = _ccDic->find_residue_type(_ResName);

            int heavy_atom = drug.nheavyatoms();
            int extra_atom = drug.natoms() + 4;

            _atomOrder.clear();
            IntegerStringOrder idx;
            for (unsigned i = 0; i < _atoms.size(); ++i) {
                 std::string atom_name = _atoms[i]->pdb_atmnam();
                 int pos = drug.find_atom_pos(atom_name);
                 if (pos < 0 && _atoms[i]->is_hydrogen()) {
                      for (unsigned int k = 0; k < atom_name.size(); ++k) {
                           if (atom_name[k] == 'D') {
                                atom_name[k] = 'H';
                                break;
                           }
                      }
                      pos = drug.find_atom_pos(atom_name);
                 }
                 if (pos >= 0) {
                      if (_atoms[i]->atom_type() == "H" || _atoms[i]->atom_type() == "D") {
                           idx.set_int_value(pos+3);
                      } else idx.set_int_value(pos);
                      if (residueType == "ATOMP" && (_atoms[i]->pdb_atmnam() == "H2" || _atoms[i]->pdb_atmnam() == "D2")) {
                           idx.set_int_value(heavy_atom+2);
                      }
                 } else if (_atoms[i]->pdb_atmnam() == "H1" || _atoms[i]->pdb_atmnam() == "D1" || _atoms[i]->pdb_atmnam() == "H2" ||
                            _atoms[i]->pdb_atmnam() == "D2" || _atoms[i]->pdb_atmnam() == "H3" || _atoms[i]->pdb_atmnam() == "D3") {
                      if (_atoms[i]->pdb_atmnam() == "H1" || _atoms[i]->pdb_atmnam() == "D1")
                           idx.set_int_value(heavy_atom+1);
                      if (_atoms[i]->pdb_atmnam() == "H2" || _atoms[i]->pdb_atmnam() == "D2")
                           idx.set_int_value(heavy_atom+2);
                      if (_atoms[i]->pdb_atmnam() == "H3" || _atoms[i]->pdb_atmnam() == "D3")
                           idx.set_int_value(heavy_atom+3);
                 } else {
                      idx.set_int_value(extra_atom++);
                 }
                 idx.set_str_value(_atoms[i]->alt_loc());
                 _atomOrder.insert(std::make_pair(idx, i));
            }
            _atom_pos = _atomOrder.begin();
       } catch (const std::exception& exc) {}
}

void Residue::find_extra_H_atom()
{
       std::string atomname = _ResName + "_H"; 
       std::multimap<std::string, int>::iterator ret = _atomNames.find(atomname);

       atomname = _ResName + "_H1";
       std::multimap<std::string, int>::iterator ret1 = _atomNames.find(atomname);

       atomname = _ResName + "_H2";
       std::multimap<std::string, int>::iterator ret2 = _atomNames.find(atomname);

       atomname = _ResName + "_H3";
       std::multimap<std::string, int>::iterator ret3 = _atomNames.find(atomname);

       if (ret != _atomNames.end()) {
            if (ret1 != _atomNames.end() || ret2 != _atomNames.end() || ret3 != _atomNames.end())
                 _extras.push_back("H is extra");
            else _extras.push_back("H");
            if (ret1 != _atomNames.end() && ret2 != _atomNames.end() && ret3 != _atomNames.end())
                 _extras.push_back("(H1, H2, H3 already exist)");
            else if (ret1 != _atomNames.end() && ret2 != _atomNames.end())
                 _extras.push_back("(H1, H2 already exist)");

            else if (ret1 != _atomNames.end() && ret3 != _atomNames.end())
                 _extras.push_back("(H1, H3 already exist)");
            else if (ret2 != _atomNames.end() && ret3 != _atomNames.end())
                 _extras.push_back("(H2, H3 already exist)");
            else if (ret1 != _atomNames.end())
                 _extras.push_back("(H1 already exist)");
            else if (ret2 != _atomNames.end())
                 _extras.push_back("(H1 already exist)");
            else if (ret3 != _atomNames.end())
                 _extras.push_back("(H3 already exist)");
       }
}

void Residue::find_missing_or_extra_atoms(const int& is_connect, const int& need_hydrogen, std::set<std::string>& LeavingAtoms, std::set<std::string>&
                                          AllowedTerminalAtoms, std::set<std::string>& NotAllowedTerminalAtoms, const bool& exclude_hydrogen)
{
       if (_ResName == "UNK" || _ResName == "UNL" || _ResName == "N") return;

       _missing.clear();
       _extras.clear();

       try {
            const ConnectFormat& drug = _ccDic->find_drug(_ResName);

            _find_missing_or_extra_atoms(is_connect, need_hydrogen, LeavingAtoms, AllowedTerminalAtoms, NotAllowedTerminalAtoms, drug, exclude_hydrogen);

            if (_extras.empty()) return;

            std::vector<ConnectFormat> variants;
            _ccDic->find_variants(_ResName, variants);
            for (unsigned int i = 0; i < variants.size(); ++i) {
                 if (variants[i].natoms() < drug.natoms()) continue;
                 _missing.clear();
                 _extras.clear();
                 _find_missing_or_extra_atoms(is_connect, need_hydrogen, LeavingAtoms, AllowedTerminalAtoms, NotAllowedTerminalAtoms,
                                              variants[i], exclude_hydrogen);
                 if (_extras.empty()) return;
            }

            _missing.clear();
            _extras.clear();
            _find_missing_or_extra_atoms(is_connect, need_hydrogen, LeavingAtoms, AllowedTerminalAtoms, NotAllowedTerminalAtoms, drug, exclude_hydrogen);
       } catch (const std::exception& exc) {}
}

void Residue::find_missing_or_extra_atoms(const ConnectFormat& drug)
{
       _missing.clear();
       _extras.clear();

       std::set<std::string> empty_set;
       empty_set.clear();

       _find_missing_or_extra_atoms(0, 1, empty_set, empty_set, empty_set, drug, false);
}

void Residue::refine_missing_atoms(const std::set<std::string>& bonded_atoms)
{
       if (_missing.empty()) return;
       try {
            std::set<std::string> connected_atom_set, missing_atom_set, additional_atom_set;
            missing_atom_set.clear();
            for (std::vector<std::string>::const_iterator vpos = _missing.begin(); vpos != _missing.end(); ++vpos) {
                 missing_atom_set.insert(*vpos);
            }

            const ConnectFormat& drug = _ccDic->find_drug(_ResName);
            std::map<std::string, std::vector<std::string> > LinkedAtoms = drug.getLinkedAtoms();

            connected_atom_set.clear();
            std::set<std::string> checking_bonded_atoms = bonded_atoms;
            while (!checking_bonded_atoms.empty()) {
                 additional_atom_set.clear();
                 for (std::set<std::string>::const_iterator spos = checking_bonded_atoms.begin(); spos != checking_bonded_atoms.end(); ++spos) {
                      std::map<std::string, std::vector<std::string> >::const_iterator mpos = LinkedAtoms.find(*spos);
                      if (mpos == LinkedAtoms.end()) continue;
                      for (std::vector<std::string>::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                           if (missing_atom_set.find(*vpos) == missing_atom_set.end()) continue;
                           if (connected_atom_set.find(*vpos) != connected_atom_set.end()) continue;
                           connected_atom_set.insert(*vpos);
                           additional_atom_set.insert(*vpos);
                      }
                 }
                 checking_bonded_atoms = additional_atom_set;
            }

            std::vector<std::string> true_missing;
            true_missing.clear();
            for (std::vector<std::string>::const_iterator vpos = _missing.begin(); vpos != _missing.end(); ++vpos) {
                 if (connected_atom_set.find(*vpos) != connected_atom_set.end()) continue;
                 true_missing.push_back(*vpos);
            }
            _missing = true_missing;
       } catch (const std::exception& exc) {}
}

void Residue::_find_missing_or_extra_atoms(const int& is_connect, const int& need_hydrogen, std::set<std::string>& LeavingAtoms,
                                           std::set<std::string>& AllowedTerminalAtoms, std::set<std::string>& NotAllowedTerminalAtoms,
                                           const ConnectFormat& drug, const bool& exclude_hydrogen)
{
       std::set<std::string> match_index;
       match_index.clear();

       const std::vector<AtomFormat>& atoms = drug.atoms();
       for (unsigned int i = 0; i < atoms.size(); ++i) {
            if (atoms[i].leaving_atom_flag() == "Y" && NotAllowedTerminalAtoms.find(atoms[i].atomname()) != NotAllowedTerminalAtoms.end()) {
                 bool skip = true;
                 if (!(is_connect & 2) && (((_ResName == "ASP" || _ResName == "DAS") && atoms[i].atomname() == "OD2") ||
                    ((_ResName == "GLU" || _ResName == "DGL") && atoms[i].atomname() == "OE2"))) skip = false;
                 else if (!(is_connect & 1) && atoms[i].atomname() == "OXT" && (_ResName == "ASP" || _ResName == "DAS" ||
                     _ResName == "DGL" || _ResName == "GLU")) skip = false;
                 else if (atoms[i].atomname() == "OXT" && !SeqCodeUtil::is_standard_aa_residue(_ResName)) skip = false;
                 if (skip) continue;
            }

            std::string cs = CompositeIndex::getIndex(_ResName, atoms[i].atomname());
            std::multimap<std::string, int>::const_iterator pos = _atomNames.find(cs);
            if ((pos == _atomNames.end()) && (atoms[i].atomtype() == "H")) {
                 std::string atomname = atoms[i].atomname();
                 for (unsigned int j = 0; j < atomname.size(); ++j) {
                      if (atomname[j] == 'H') {
                           atomname[j] = 'D';
                           break;
                      }
                 }
                 cs = CompositeIndex::getIndex(_ResName, atomname);
                 pos = _atomNames.find(cs);
            }

            if (pos != _atomNames.end()) match_index.insert(pos->first);
            else {
                 if (atoms[i].ishydrogen() && !need_hydrogen)
                      continue;
                 else if (LeavingAtoms.find(atoms[i].atomname()) != LeavingAtoms.end())
                      continue;
                 else if (is_connect && atoms[i].leaving_atom_flag() == "Y") {
                      bool skip = true;
                      if (!(is_connect & 2) && (((_ResName == "ASP" || _ResName == "DAS") && atoms[i].atomname() == "OD2") ||
                         ((_ResName == "GLU" || _ResName == "DGL") && atoms[i].atomname() == "OE2"))) skip = false;
                      if (!(is_connect & 1) && atoms[i].atomname() == "OXT" && (_ResName == "ASP" || _ResName == "DAS" ||
                          _ResName == "DGL" || _ResName == "GLU")) skip = false;
                      if (skip) continue;
                 } else if ((_ResName == "HOH" || _ResName == "DOD") &&
                             atoms[i].ishydrogen()) continue;

                 _missing.push_back(atoms[i].atomname());
            }
       }

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (exclude_hydrogen && (*pos)->is_hydrogen()) continue;
            if (match_index.find((*pos)->getAtomIndex()) != match_index.end()) continue;

            if (NotAllowedTerminalAtoms.find((*pos)->pdb_atmnam()) != NotAllowedTerminalAtoms.end()) {
                 _extras.push_back((*pos)->pdb_atmnam());
                 continue;
            }

            int p = drug.find_atom_pos((*pos)->pdb_atmnam());
            if (p >= 0) continue;

            std::string atomname = (*pos)->pdb_atmnam();
            if ((*pos)->is_hydrogen() && (*pos)->atom_type() == "D") {
                 for (unsigned int j = 0; j < atomname.size(); ++j) {
                      if (atomname[j] == 'D') {
                           atomname[j] = 'H';
                           break;
                      }
                 }
                 p = drug.find_atom_pos(atomname);
                 if (p >= 0) continue;
            }

            if (AllowedTerminalAtoms.find(atomname) != AllowedTerminalAtoms.end()) continue;

            _extras.push_back((*pos)->pdb_atmnam());
       }
}

void Residue::find_atom_type_mismatch(const int& Mol_ID)
{
       if (_ResName == "UNL" || _ResName == "N" || _ResName == "UNX") return;

       try {
            const ConnectFormat& drug = _ccDic->find_drug(_ResName);

            char buffer[200];
            std::string dict_atomtype, atomtype, atomName;
            for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
                 if ((*pos)->atom_type().empty()) continue;

                 dict_atomtype.clear();
                 try {
                      const AtomFormat& atom = drug.find_atom((*pos)->pdb_atmnam());
                      dict_atomtype = atom.atomtype();
     
                 } catch (const std::exception& exc) {
                      if ((*pos)->atom_type() == "D") {
                           atomName = (*pos)->pdb_atmnam();
                           atomName[0] = 'H';
                           try {
                                const AtomFormat& atom = drug.find_atom(atomName);
                                dict_atomtype = atom.atomtype();
                           } catch (const std::exception& exc) {}
                      }
                 }

                 bool atom_type_match = true;

                 String::UpperCase((*pos)->pdb_atmnam(), atomName);
                 String::UpperCase((*pos)->atom_type(), atomtype);
                 if (atomName.substr(0, atomtype.size()) != atomtype) atom_type_match = false;

                 if (!dict_atomtype.empty()) {
                      String::UpperCase(dict_atomtype);
                      if (!((dict_atomtype == atomtype) || ((dict_atomtype == "H") && (atomtype == "D")))) atom_type_match = false;
                 }
                 if (atom_type_match) continue;
                 
                 std::string Chain_ID = (*pos)->pdb_chnid();
                 if (Chain_ID.empty()) Chain_ID = "  ";
                 if (Mol_ID < 0)
                      sprintf(buffer, FORMAT1, (*pos)->pdb_resnam().c_str(), Chain_ID.c_str(), (*pos)->pdb_resnum().c_str(),
                              (*pos)->pdb_atmnam().c_str(), (*pos)->atom_type().c_str(), dict_atomtype.c_str());
                 else sprintf(buffer, FORMAT2, (*pos)->pdb_resnam().c_str(), Mol_ID, Chain_ID.c_str(), (*pos)->pdb_resnum().c_str(),
                              (*pos)->pdb_atmnam().c_str(), (*pos)->atom_type().c_str(), dict_atomtype.c_str());
                 _messageIo->insertMessage("atom_type", "label", "true");
                 _messageIo->insertMessage("atom_type_print", "error", buffer);
            }
       } catch (const std::exception& exc) {}
}

void Residue::check_hydrogen_bond_distance(const int& Mol_ID)
{
       if (_ResName == "UNK" || _ResName == "UNL" || _ResName == "N") return;

       try {
            const ConnectFormat& drug = _ccDic->find_drug(_ResName);

            bool is_cys = false;
            if (_ResName == "CYS") is_cys = true;

            std::vector<RCSB::Atom*> atom_list;
            for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
                 if (!(*pos)->is_hydrogen()) continue;

                 bool is_hg = false;
                 if (is_cys && ((*pos)->pdb_atmnam() == "HG" || ((*pos)->atom_type() == "D" || (*pos)->pdb_atmnam() == "DG"))) is_hg = true;
                 bool is_deuterium = false;
                 if ((*pos)->atom_type() == "D") is_deuterium = true;

                 int ipos = drug.find_atom_pos((*pos)->pdb_atmnam());
                 if (ipos < 0 && (*pos)->atom_type() == "D") {
                      std::string atomname = (*pos)->pdb_atmnam();
                      for (unsigned int j = 0; j < atomname.size(); ++j) {
                           if (atomname[j] == 'D') {
                                atomname[j] = 'H';
                                break;
                           }
                      }
                      ipos = drug.find_atom_pos(atomname);
                 }
                 if (ipos < 0) continue;

                 const std::vector<std::vector<std::string> >& bonds = drug.find_atom(ipos).bonds();
                 if (bonds.empty()) continue;

                 find_atom(bonds[0][0], atom_list);
                 if (atom_list.empty()) {
                      _wrong_connectivity = true;
                      std::string cs = "Atom '" + (*pos)->atmtype() + " " + (*pos)->pdb_resnam() + " " + (*pos)->pdb_chnid() + " " + (*pos)->pdb_resnum()
                                     + (*pos)->ins_code() + "'";
                      if (Mol_ID >= 0) cs += " of Model " + String::IntToString(Mol_ID);
                      cs += " is not connected.";
                      _messageIo->insertMessage("covalent_h_bond", "label", "true");
                      _messageIo->insertMessage("covalent_h_bond_print", "warning", cs);
                      continue;
                 }
                 for (unsigned int k = 0; k < atom_list.size(); ++k) {
                      if (!(*pos)->alt_loc().empty() && !atom_list[k]->alt_loc().empty() && (*pos)->alt_loc() != atom_list[k]->alt_loc()) continue;

                      double distance = cal_distance((*pos), atom_list[k]);
                      if (distance >= 0.6 && distance <= 1.2) continue;
                      if (is_hg && distance >= 0.6 && distance <= 1.4) continue;
                      if (is_deuterium && distance >= 0.6 && distance <= 1.6) continue;
                      if (BondUtil::is_a_bond(*pos, atom_list[k])) continue;

                      _wrong_connectivity = true;
                      std::string cs = "The bond distance is " + FloatToString(distance, 0, 3) + " between '" + (*pos)->atmtype() + " "
                                     + (*pos)->pdb_resnam() + " " + (*pos)->pdb_chnid() + " " + (*pos)->pdb_resnum() + (*pos)->ins_code() + "' and '"
                                     + atom_list[k]->atmtype() + " " + atom_list[k]->pdb_resnam() + " " + atom_list[k]->pdb_chnid() + " "
                                     + atom_list[k]->pdb_resnum() + atom_list[k]->ins_code() + "'";
                      if (Mol_ID >= 0) cs += " of Model " + String::IntToString(Mol_ID);
                      cs += ".";
                      _messageIo->insertMessage("covalent_h_bond", "label", "true");
                      _messageIo->insertMessage("covalent_h_bond_print", "warning", cs);
                 }
            }
       } catch (const std::exception& exc) {}
}

void Residue::find_wrong_alt_id(const int& Mol_ID, std::string& c_9999_records)
{
       // adding check '*' in atom name
       char buffer[200];
       std::string record;

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if ((*pos)->is_9999_flag()) {
                 _writeAtomRecord(*pos, record);
                 c_9999_records += record;
            }
            if ((*pos)->alt_loc().empty()) continue;
            if (((*pos)->alt_loc()[0] >= 'A' && (*pos)->alt_loc()[0] <= 'Z') || ((*pos)->alt_loc()[0] >= 'a' && (*pos)->alt_loc()[0] <= 'z')) continue;
            std::string Chain_ID = (*pos)->pdb_chnid();
            if (Chain_ID.empty()) Chain_ID = " ";
            if (Mol_ID < 0)
                 sprintf(buffer, FORMAT3, (*pos)->pdb_atmnam().c_str(), _ResName.c_str(), Chain_ID.c_str(), (*pos)->pdb_resnum().c_str(),
                        (*pos)->x(), (*pos)->y(), (*pos)->z(), (*pos)->alt_loc().c_str());
            else sprintf(buffer, FORMAT4, Mol_ID, (*pos)->pdb_atmnam().c_str(), _ResName.c_str(), Chain_ID.c_str(), (*pos)->pdb_resnum().c_str(),
                        (*pos)->x(), (*pos)->y(), (*pos)->z(), (*pos)->alt_loc().c_str());
            _messageIo->insertMessage("alt_code", "label", "true");
            _messageIo->insertMessage("alt_code_print", "warning", buffer);
       }
}

void Residue::check_same_coordinates(const int& Mol_ID)
{
       for (unsigned i = 0; i < _atoms.size() - 1; ++i) {
            for (unsigned j = i + 1; j < _atoms.size(); ++j) {
                 if (_atoms[i]->alt_loc().empty() || _atoms[j]->alt_loc().empty() || (_atoms[i]->alt_loc() == _atoms[j]->alt_loc()) /* ||
                    (_atoms[i]->alt_loc() != _atoms[j]->alt_loc() && _atoms[i]->pdb_atmnam() == _atoms[j]->pdb_atmnam()) */ ) {
                      float dx = fabs(_atoms[i]->x() - _atoms[j]->x());
                      float dy = fabs(_atoms[i]->y() - _atoms[j]->y());
                      float dz = fabs(_atoms[i]->z() - _atoms[j]->z());
                      float db = fabs(atof(_atoms[i]->t_fct().c_str()) - atof(_atoms[j]->t_fct().c_str()));
                      if (dx < 0.001 && dy < 0.001 && dz < 0.001 && db < 0.001)
                           _writeErrorMessage(_atoms[i], _atoms[j], "same_coor", "Atoms with same coordinates:\n", Mol_ID);
                 } 
            }
       }
}

bool Residue::is_zero_occupancy_residue(std::vector<std::vector<std::string> > &atom_list)
{
       atom_list.clear();

       std::vector<std::string> data;

       int total_count = 0;
       int zero_count = 0;
       std::set<std::string> UniIndex;
       UniIndex.clear();
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            total_count++;
            if (fabs(atof((*pos)->occ().c_str())) < ZERO_OCCUPANCY) {
                 zero_count++;
                 if (!(*pos)->is_hydrogen()) {
                      std::string cs = (*pos)->pdb_atmnam() + "_" + (*pos)->alt_loc();
                      if (UniIndex.find(cs) == UniIndex.end()) {
                           UniIndex.insert(cs);
                           data.clear();
                           data.push_back((*pos)->pdb_atmnam());
                           data.push_back((*pos)->alt_loc());
                           atom_list.push_back(data);
                      }
                 }
            }
       }

       if (zero_count && zero_count == total_count)
            return true;
       else return false;
}

bool Residue::is_zero_occupancy_ligand()
{
       int total_count = 0;
       int zero_count = 0;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            total_count++;
            if (fabs(atof((*pos)->occ().c_str())) < ZERO_OCCUPANCY_LIGAND) zero_count++;
       }

       if (zero_count && ((zero_count * 10) >= (total_count * 8)))
            return true;
       else return false;
}

bool Residue::changed_zero_occupancy()
{
       bool changed = false;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (fabs(atof((*pos)->occ().c_str())) < ZERO_OCCUPANCY) {
                 (*pos)->setValue("1.0", 11);
                 changed = true;
            }
       }

       return changed;
}

float Residue::get_partial_occupancy()
{
       float partial_occupancy = 1.0;

       std::set<std::string> unique_occ_set;
       unique_occ_set.clear();
       for (std::multimap<std::string, int>::const_iterator mpos = _atomNames.begin(); mpos != _atomNames.end(); ++mpos) {
            std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNames.equal_range(mpos->first);
            if (range.first == range.second) continue;

            float occ = 0.0;
            for (std::multimap<std::string, int>::iterator pos = range.first; pos != range.second; ++pos) {
                 occ += atof(_atoms[pos->second]->occ().c_str());
            }
            unique_occ_set.insert(FloatToString(occ, 0, 2));
       }
       if (unique_occ_set.size() == 1) partial_occupancy = atof(unique_occ_set.begin()->c_str());

       return partial_occupancy;
}

float Residue::get_q_score()
{
       float q_score_sum = 0;
       int count = 0;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            std::string atom_q_socre = (*pos)->getValue("Q-score");
            if (atom_q_socre.empty() || !String::IsNumber(atom_q_socre)) continue;
            q_score_sum += atof(atom_q_socre.c_str());
            count += 1;
       }
       if (count > 0)
            return (q_score_sum / (float) count);
       else throw NotFoundException();
}

bool Residue::is_pdb_format_compatible()
{
       // Five letter residue name
       if (_ResName.size() > 3) return false;

       for (std::vector<RCSB::Atom*>::const_iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            // B-factor > 1000
            if (atof((*pos)->t_fct().c_str()) > 999.995) return false;
       }

       return true;
}

void Residue::check_occupancy(const int& Mol_ID, std::list<RCSB::Atom*>& zero_occupancy_atom_list)
{
       if (_atoms.empty()) return;

       float occupancy = 0;
       std::string name = "XXXX";
       for (std::vector<RCSB::Atom*>::const_iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            float occ = atof((*pos)->occ().c_str());
            if (name != (*pos)->pdb_atmnam()) {
                 if (occupancy > 1.005) _writeOccupancyMessage(name, occupancy, Mol_ID);
                 name = (*pos)->pdb_atmnam();
                 occupancy = 0;
            }
            occupancy += occ;
            if (fabs(occ) < ZERO_OCCUPANCY) zero_occupancy_atom_list.push_back(*pos);
            if (occ < 0) _writeOccupancyMessage((*pos)->pdb_atmnam(), occ, Mol_ID);

       }
       if (occupancy > 1.005) _writeOccupancyMessage(name, occupancy, Mol_ID);
}

bool Residue::check_occupancy_and_altloc(const int& Mol_ID)
{
       bool error_free_flag = true;

       std::map<std::string, int> altloc_count_mapping;
       altloc_count_mapping.clear();

       float occupancy = 0;
       std::string name = "____XXXX";
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            std::string idx = (*pos)->pdb_atmnam() + "_" + (*pos)->alt_loc();
            std::map<std::string, int>::iterator
                mpos = altloc_count_mapping.find(idx);
            if (mpos != altloc_count_mapping.end())
                 mpos->second++;
            else altloc_count_mapping.insert(std::make_pair(idx, 1));

            float occ = atof((*pos)->occ().c_str());
            if (name != (*pos)->pdb_atmnam()) {
                 if (occupancy > 1.005) {
                      _writeOccupancyError(name, occupancy, Mol_ID);
                      error_free_flag = false;
                 }
                 name = (*pos)->pdb_atmnam();
                 occupancy = 0;
            }
            occupancy += occ;

       }
       if (occupancy > 1.005) {
            _writeOccupancyError(name, occupancy, Mol_ID);
            error_free_flag = false;
       }

       std::vector<std::string> data;
       std::string error;
       for (std::map<std::string, int>::const_iterator mpos = altloc_count_mapping.begin(); mpos != altloc_count_mapping.end(); ++mpos) {
            if (mpos->second < 2) continue; 
            get_wordarray_delimit_by_string(data, mpos->first, "_");

            error.clear();
            if (Mol_ID >= 0) error += "In model " + String::IntToString(Mol_ID) + ",";

            error += "Atom '" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + " " + data[0];
            if (!data[1].empty()) error += "[InsCode:" + data[1] + "]";
            error += "' repeats " + String::IntToString(mpos->second) + " times.\n";
            _logIo->message(error.c_str());
            error_free_flag = false;
       }

       return error_free_flag;
}

std::string Residue::check_b_factor(int& atom_count, int& high_b_factor_count, int& low_b_factor_count, int& zero_b_factor_count,
                                    std::list<RCSB::Atom*>& negative_b_factor_atoms)
{
       if (_atoms.empty()) return "";

       bool has_same_b_factor = true;
       float t_fct = -1000.0;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            atom_count++;

            float fct = atof((*pos)->t_fct().c_str());
            if (fct > HIGH_B_FACTOR) high_b_factor_count++;
            if (fct < LOW_B_FACTOR) low_b_factor_count++;
            if (fct < LESS_ZERO_B_FACTOR) negative_b_factor_atoms.push_back(*pos);
            else if (fct < ZERO_B_FACTOR) zero_b_factor_count++;
            if (t_fct < -999.0) t_fct = fct;
            if (fabs(fct - t_fct) > 0.01) has_same_b_factor = false;

       }

       if (has_same_b_factor) return FloatToString(t_fct, 0, 2);

       return "";
}

void Residue::check_occupancy_and_b_factor(std::vector<std::string>& occ_list, std::vector<std::string>& fct_list)
{
       if (_atoms.empty()) return;

       float occupancy = 0;
       std::string name = "XXXX";
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            std::string atomName = (*pos)->pdb_atmnam();
            if (!(*pos)->alt_loc().empty()) atomName = (*pos)->pdb_atmnam() + "(" + (*pos)->alt_loc() + ")";

            float occ = atof((*pos)->occ().c_str());
            if (occ < 0) {
                 occ_list.push_back("Atom '" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + " " + atomName
                                  + "' has occupancy " + FloatToString(occ, 0, 2) + ".");
            }
            if (name != (*pos)->pdb_atmnam()) {
                 if (occupancy > 1.005) {
                      occ_list.push_back("Atom '" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + " " + name
                                       + "' has occupancy " + FloatToString(occupancy, 0, 2) + ".");
                 }
                 name = (*pos)->pdb_atmnam();
                 occupancy = 0;
            }
            occupancy += occ;

            float fct = atof((*pos)->t_fct().c_str());
            if (fct < ZERO_B_FACTOR) {
                 fct_list.push_back("Atom '" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + " " + atomName
                                  + "' has temperature factor " + FloatToString(fct, 0, 2) + ".");
            }
       }
       if (occupancy > 1.005) {
            occ_list.push_back("Atom '" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + " " + name
                             + "' has occupancy " + FloatToString(occupancy, 0, 2) + ".");
       }
}

void Residue::update_occupancy(const std::string& new_occ, const bool& check_full_occ_flag)
{
       if (_atoms.empty()) return;

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (!check_full_occ_flag || (check_full_occ_flag && (fabs(atof((*pos)->occ().c_str()) - 1.0) < 0.001))) (*pos)->setValue(new_occ, 11);
       }
}

void Residue::CorresInfoChecking(bool& high_b_factor, std::list<RCSB::Atom*>& zero_occupancy_atom_list, std::list<RCSB::Atom*>& large_occupancy_atom_list,
                                 std::list<RCSB::Atom*>& zero_b_factor_atom_list)
{
       high_b_factor = false; 
       if (_atoms.empty()) return;

       int atom_count = 0;
       int high_b_factor_count = 0;
       float occupancy = 0;
       std::string name = "XXXX";

       std::list<RCSB::Atom*> atom_list;
       atom_list.clear();

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            atom_count++;
            float fct = atof((*pos)->t_fct().c_str());
            if (fct > HIGH_B_FACTOR) high_b_factor_count++;
            else if (fct < ZERO_B_FACTOR) zero_b_factor_atom_list.push_back(*pos);

            float occ = atof((*pos)->occ().c_str());
            if (fabs(occ) < ZERO_OCCUPANCY) zero_occupancy_atom_list.push_back(*pos);

            if (name != (*pos)->pdb_atmnam()) {
                 if (occupancy > 1.005 && !atom_list.empty()) {
                      for (std::list<RCSB::Atom*>::const_iterator lpos = atom_list.begin(); lpos != atom_list.end(); ++lpos) {
                           large_occupancy_atom_list.push_back(*lpos);
                      }
                 }
                 name = (*pos)->pdb_atmnam();
                 occupancy = 0;
                 atom_list.clear();
            }
            occupancy += occ;
            atom_list.push_back(*pos);
       }
       if (occupancy > 1.005 && !atom_list.empty()) {
            for (std::list<RCSB::Atom*>::const_iterator lpos = atom_list.begin(); lpos != atom_list.end(); ++lpos) {
                 large_occupancy_atom_list.push_back(*lpos);
            }
       }

       if (atom_count == high_b_factor_count) high_b_factor = true;
}

bool Residue::change_water_nomenclature()
{
       bool changed = false;
       bool changed_residue = false;
       if (_ResName == "H2O" || _ResName == "WAT" || _ResName == "TIP" || _ResName == "WTR" || _ResName == "OH2") {
            _ResName = "HOH";
            changed_residue = true;
            changed = true;
       }
       if (_ResName != "HOH" && _ResName != "DOD") return  changed_residue;

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (changed) {
                 (*pos)->set_restype(_ResName);
                 (*pos)->set_pdb_resnam(_ResName);
            }

            std::string p = (*pos)->atmtype();
            if (p == "OW" || p == "WAT" || p == "OH2" || p == "WTR" || p == "H2O" || p == "W" || p[0] == 'O') p = "O";
            else if (p == "1H") p = "H1";
            else if (p == "2H") p = "H2";
            else if (p == "1D") p = "D1";
            else if (p == "2D") p = "D2";
            if (p != (*pos)->atmtype()) {
                 (*pos)->set_atmtype(p);
                 (*pos)->set_pdb_atmnam(p);
                 changed = true;
            }
       }

       if (changed) {
            if (_atomNameAlts_Orininal.empty()) _atomNameAlts_Orininal = _atomNameAlts;
            _UpdateIndices();
       }
       return changed_residue;
}

void Residue::AtomNameUpperCase()
{
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            (*pos)->UpperCase();
       }
       _atomNameAlts_Orininal.clear();
       _UpdateIndices();
       reorder_atoms();
}

void Residue::update_nomenclature(const std::string &type, const std::string& res_num)
{
       _res_no = res_num;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (!type.empty()) (*pos)->set_type(type);
            (*pos)->set_chnid(_chnid);
            (*pos)->set_resnum(res_num);
       }
}

void Residue::update_nomenclature(const std::string &type, const std::string &chain_id, const std::string& res_num)
{
       _chnid = chain_id;
       _res_no = res_num;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (!type.empty()) (*pos)->set_type(type);
            (*pos)->set_chnid(chain_id);
            (*pos)->set_resnum(res_num);
       }
}

void Residue::update_nomenclature(const std::string &type, const std::string &chain_id, const std::string& res_num, const std::string &pdb_chainid,
                                  const std::string& pdb_res_num, const std::string& ins_code)
{
       _chnid = chain_id;
       _res_no = res_num;
       _pdb_chnid = pdb_chainid;
       _pdb_res_no = pdb_res_num;
       _ins_code = ins_code;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (!type.empty()) (*pos)->set_type(type);
            (*pos)->set_chnid(chain_id);
            (*pos)->set_resnum(res_num);
            (*pos)->set_pdb_chnid(pdb_chainid);
            (*pos)->set_pdb_resnum(pdb_res_num);
            (*pos)->set_ins_code(ins_code);
       }
}

void Residue::update_nomenclature(const std::map<std::string, std::string>& mapping)
{
       if (mapping.empty()) return;

       bool chainid_flag = false;
       bool number_flag = false;
       bool inscode_flag = false;
       std::map<std::string, std::string>::const_iterator mpos = mapping.find("chainid");
       if (mpos != mapping.end()) {
            _pdb_chnid = mpos->second;
            chainid_flag = true;
       }

       mpos = mapping.find("number");
       if (mpos != mapping.end()) {
            _pdb_res_no = mpos->second;
            number_flag = true;
       }

       mpos = mapping.find("inscode");
       if (mpos != mapping.end()) {
            _ins_code = mpos->second;
            inscode_flag = true;
       }

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (chainid_flag) (*pos)->set_pdb_chnid(_pdb_chnid);
            if (number_flag) (*pos)->set_pdb_resnum(_pdb_res_no);
            if (inscode_flag) (*pos)->set_ins_code(_ins_code);
       }
}

void Residue::update_pdb_chnid(const std::string &pdb_chainid)
{
       _pdb_chnid = pdb_chainid;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            (*pos)->set_pdb_chnid(pdb_chainid);
       }
}

void Residue::correction_name(const std::string &resname)
{
       if (_ResName == resname) return;

       _ResName = resname;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            (*pos)->set_restype(_ResName);
            (*pos)->set_pdb_resnam(_ResName);
       }
       _atomNameAlts_Orininal.clear();
       _UpdateIndices();
       reorder_atoms();

}

void Residue::convert_between_orthogonal_and_fractional(CrySymmetry& crySymm, const int& iflag)
{
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            COORD coord = (*pos)->orig();
            crySymm.ndb_trans_coord_3by3(coord, iflag);
            (*pos)->set_orig(coord);
       }
}

void Residue::symmetry_operations(CrySymmetry& crySymm, const int& lx, const int& ly, const int& lz, const NDBSYMMETRY& op, const bool& not_moving)
{
       bool exist_anisou = false;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if ((*pos)->has_anisou()) {
                 exist_anisou = true;
                 break;
            }
       }
       if (exist_anisou && not_moving) return;

       COORD coord;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            crySymm.do_symmetry_operation(coord, (*pos)->orig(), op);
            coord.x += (float) lx;
            coord.y += (float) ly;
            coord.z += (float) lz;
            (*pos)->set_orig(coord);
       }
}

void Residue::transformation(double m[4][4])
{
       COORD coord;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            matrix_rotate_translate(coord, (*pos)->orig(), m);
            (*pos)->set_orig(coord);
       }
}

bool Residue::rename(const std::string& sidechain_pattern, const std::string& NewName)
{
       if (!SeqCodeUtil::is_standard_aa_residue_plus_MSE(NewName) && (NewName != "UNK")) return false;

       std::string key_atom = "";

       std::set<std::string> delete_atom_set;

       SideChainPattern::get_key_atoms(sidechain_pattern, delete_atom_set);
       if (!delete_atom_set.empty()) {
            std::vector<std::string> key_atom_list;
            key_atom_list.clear();
            for (std::set<std::string>::const_iterator spos = delete_atom_set.begin(); spos != delete_atom_set.end(); ++spos) {
                 if (find_atom(*spos)) key_atom_list.push_back(*spos);
            }
            if (key_atom_list.size() != 1) return false;
            key_atom = key_atom_list[0];
       }

       std::map<std::string, std::pair<std::string, std::string> > rename_atom_map;

       if (!SideChainPattern::get_change_info(key_atom, sidechain_pattern, _ResName, NewName, delete_atom_set, rename_atom_map)) return false;

       for (std::map<std::string, std::pair<std::string, std::string> >::const_iterator mpos = rename_atom_map.begin(); mpos != rename_atom_map.end(); ++mpos) {
            std::string cs = CompositeIndex::getIndex(_ResName, mpos->first);
            std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNames.equal_range(cs);
            if (range.first == range.second) continue;

            for (std::multimap<std::string, int>::const_iterator pos = range.first; pos != range.second; ++pos) {
                 _atoms[pos->second]->set_atmtype(mpos->second.first);
                 _atoms[pos->second]->set_pdb_atmnam(mpos->second.first);
                 if (!mpos->second.second.empty()) _atoms[pos->second]->set_atom_type(mpos->second.second);
            }
       }

       if (!delete_atom_set.empty()) {
            std::set<int> delete_index;
            delete_index.clear();
            for (std::set<std::string>::const_iterator spos = delete_atom_set.begin(); spos != delete_atom_set.end(); ++spos) {
                 std::string cs = CompositeIndex::getIndex(_ResName, *spos);
                 std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNames.equal_range(cs);
                 if (range.first == range.second) continue;

                 for (std::multimap<std::string, int>::const_iterator pos = range.first; pos != range.second; ++pos) {
                      delete_index.insert(pos->second);
                 }
            }

            if (!delete_index.empty()) {
                 std::set<long> atom_set;
                 atom_set.clear();
                 _delete_atoms(delete_index, atom_set);
            }
       }

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            (*pos)->set_restype(NewName);
            (*pos)->set_pdb_resnam(NewName);
       }

       _OrigResName = _ResName;
       _ResName = NewName;
       _atomNameAlts_Orininal.clear();
       if (_UpdateIndices()) _change_hydrogen_names();
       reorder_atoms();

       return true;
}

void Residue::fix_n_terminal_hydrogen(const bool& first_residue, const std::vector<std::map<std::string, std::vector<std::string> > >& atom_name_mapping_list)
{
       const std::string& N_atom = _ccDic->find_terminal_atom(_ResName, "N");
       if (N_atom.empty()) return;

       std::vector<RCSB::Atom*> terminal_h_atoms = _get_terminal_hydrogens(N_atom);
       if (terminal_h_atoms.empty()) return;

       std::vector<RCSB::Atom*> atom_list;
       std::map<std::string, std::vector<RCSB::Atom*> > terminal_h_atom_name_mapping;
       terminal_h_atom_name_mapping.clear();

       for (std::vector<RCSB::Atom*>::const_iterator pos = terminal_h_atoms.begin(); pos != terminal_h_atoms.end(); ++pos) {
            std::map<std::string, std::vector<RCSB::Atom*> >::iterator mpos = terminal_h_atom_name_mapping.find((*pos)->pub_atmnam());
            if (mpos != terminal_h_atom_name_mapping.end()) mpos->second.push_back(*pos);
            else {
                 atom_list.clear();
                 atom_list.push_back(*pos);
                 terminal_h_atom_name_mapping.insert(std::make_pair((*pos)->pub_atmnam(), atom_list));
            }
       }

       std::vector<std::pair<std::vector<std::string>, std::vector<RCSB::Atom*> > > final_mapped_list, tmp_mapped_list;
       final_mapped_list.clear();

       for (std::vector<std::map<std::string, std::vector<std::string> > >::const_iterator vpos = atom_name_mapping_list.begin();
            vpos != atom_name_mapping_list.end(); ++vpos) {
            tmp_mapped_list.clear();
            for (std::map<std::string, std::vector<std::string> >::const_iterator name_mpos = vpos->begin(); name_mpos != vpos->end(); ++name_mpos) {
                 std::map<std::string, std::vector<RCSB::Atom*> >::const_iterator atom_mpos = terminal_h_atom_name_mapping.find(name_mpos->first);
                 if (atom_mpos != terminal_h_atom_name_mapping.end()) tmp_mapped_list.push_back(std::make_pair(name_mpos->second, atom_mpos->second));
            }
            if (tmp_mapped_list.empty()) continue;
            if (tmp_mapped_list.size() > final_mapped_list.size()) final_mapped_list = tmp_mapped_list;
            if (final_mapped_list.size() == terminal_h_atom_name_mapping.size()) break;
       }
       if (final_mapped_list.empty()) return;

       for (std::vector<std::pair<std::vector<std::string>, std::vector<RCSB::Atom*> > >::iterator vpos = final_mapped_list.begin();
            vpos != final_mapped_list.end(); ++vpos) {
            for (std::vector<RCSB::Atom*>::iterator pos = vpos->second.begin(); pos != vpos->second.end(); ++pos) {
                 if (final_mapped_list.size() == 1) {
                      if (first_residue) {
                           (*pos)->set_atmtype(vpos->first[1]);
                           (*pos)->set_pdb_atmnam(vpos->first[1]);
                      } else {
                           (*pos)->set_atmtype(vpos->first[2]);
                           (*pos)->set_pdb_atmnam(vpos->first[2]);
                      }
                 } else {
                      (*pos)->set_atmtype(vpos->first[0]);
                      (*pos)->set_pdb_atmnam(vpos->first[0]);
                 }
            }
       }

       _atomNameAlts_Orininal.clear();
       _UpdateIndices();
       reorder_atoms();
}

void Residue::fix_n_terminal_hydrogen_based_on_mapping()
{
       int changed = 0;

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (_ResName == "PRO") {
                 if ((*pos)->atmtype() == "HT1" || (*pos)->atmtype() == "1H") {
                      (*pos)->set_atmtype("H2");
                      (*pos)->set_pdb_atmnam("H2");
                      changed++;
                 } else if ((*pos)->atmtype() == "HT2" || (*pos)->atmtype() == "2H") { 
                      (*pos)->set_atmtype("H3");
                      (*pos)->set_pdb_atmnam("H3");
                      changed++;
                 }
                 continue;
            }

            if ((*pos)->atmtype() == "1H") {
                 (*pos)->set_atmtype("H1");
                 (*pos)->set_pdb_atmnam("H1");
                 changed++;
            } else if ((*pos)->atmtype() == "2H") {
                 (*pos)->set_atmtype("H2");
                 (*pos)->set_pdb_atmnam("H2");
                 changed++;
            } else if ((*pos)->atmtype() == "3H") {
                 (*pos)->set_atmtype("H3");
                 (*pos)->set_pdb_atmnam("H3");
                 changed++;
            } else if ((*pos)->atmtype() == "1D") {
                 (*pos)->set_atmtype("D1");
                 (*pos)->set_pdb_atmnam("D1");
                 changed++;
            } else if ((*pos)->atmtype() == "2D") {
                 (*pos)->set_atmtype("D2");
                 (*pos)->set_pdb_atmnam("D2");
                 changed++;
            } else if ((*pos)->atmtype() == "3D") {
                 (*pos)->set_atmtype("D3");
                 (*pos)->set_pdb_atmnam("D3");
                 changed++;
            }
       }

       if (changed) {
            _atomNameAlts_Orininal.clear();
            _UpdateIndices();
            reorder_atoms();
       }
}

void Residue::fix_5_terminal_hydrogen()
{
       std::vector<RCSB::Atom*> terminal_h_atoms = _get_terminal_hydrogens("O5'");
       if (terminal_h_atoms.empty()) return;

       std::string H_atom = "HO5'";

       int changed = 0;
       for (std::vector<RCSB::Atom*>::iterator pos = terminal_h_atoms.begin(); pos != terminal_h_atoms.end(); ++pos) {
            if ((*pos)->pub_atmnam() == "H5T" || (*pos)->pub_atmnam() == "HO5*" || (*pos)->pub_atmnam() == "H5*" || (*pos)->pub_atmnam() == "HO5'") {
                 (*pos)->set_atmtype(H_atom);
                 (*pos)->set_pdb_atmnam(H_atom);
                 changed++;
            }
       }
       if (changed) {
            _atomNameAlts_Orininal.clear();
            _UpdateIndices();
            reorder_atoms();
       }
}

void Residue::fix_5_terminal_hydrogen_based_on_mapping()
{
       int changed = 0;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if ((*pos)->atmtype() == "H5T") {
                 (*pos)->set_atmtype("HO5'");
                 (*pos)->set_pdb_atmnam("HO5'");
                 changed++;
            }
       }
       if (changed) {
            _atomNameAlts_Orininal.clear();
            _UpdateIndices();
            reorder_atoms();
       }
}

void Residue::check_connectivity(const int& Mol_ID)
{
       if (_uniqueAtomNames.size() < 2) return;

       try {
            const ConnectFormat& drug = _ccDic->find_drug(_ResName);
            const std::vector<std::vector<std::string> >& bonds = drug.bonds();
            // if (bonds.size() < 2) return;

            std::string ResString = _ResName + " " + _pdb_chnid + " " + _pdb_res_no + _ins_code;
            if (Mol_ID > 0) ResString += " ( Model " + String::IntToString(Mol_ID) + ")";
            bool break_flag = false;
            std::vector<RCSB::Atom*> atom1, atom2;
            for (unsigned int i = 0; i < bonds.size(); i++) {
                 find_atom(bonds[i][0], atom1);
                 find_atom(bonds[i][1], atom2);
                 for (std::vector<RCSB::Atom*>::iterator pos1 = atom1.begin(); pos1 != atom1.end(); ++pos1) {
                      for (std::vector<RCSB::Atom*>::iterator pos2 = atom2.begin(); pos2 != atom2.end(); ++pos2) {
                           if (!(*pos1)->alt_loc().empty() && !(*pos2)->alt_loc().empty() && (*pos1)->alt_loc() != (*pos2)->alt_loc()) continue;
                           double dist = cal_distance(*pos1, *pos2);
                           if (dist <= 3.2) continue;
     
                           _messageIo->insertMessage("dissociated_residue", "warning", ResString);
                           _messageIo->insertMessage("dissociated_residue_print", "warning", "Dissociated residue: " + ResString);
                           break_flag = true;
                           break;
                      }
                      if (break_flag) break;
                 }
                 if (break_flag) break;
            }

            std::set<std::string> bond_set, extra_atom_set;
            bond_set.clear();
            for (unsigned int i = 0; i < bonds.size(); i++) {
                 bond_set.insert(bonds[i][0] + "_" + bonds[i][1]);
                 bond_set.insert(bonds[i][1] + "_" + bonds[i][0]);
                 std::string name1 = bonds[i][0];
                 std::string name2 = bonds[i][1];
                 if (name1[0] == 'H' || name2[0] == 'H') {
                      if (name1[0] == 'H') name1[0] = 'D';
                      if (name2[0] == 'H') name2[0] = 'D';
                      bond_set.insert(name1 + "_" + name2);
                      bond_set.insert(name2 + "_" + name1);
                 }
            }
            extra_atom_set.clear();
            for (std::vector<std::string>::const_iterator pos = _extras.begin(); pos != _extras.end(); ++pos) extra_atom_set.insert(*pos);

            double distance = 0.0, lower_limit = 0.0, upper_limit = 0.0;
            for (unsigned int i = 0; i < _atoms.size() - 1; i++) {
                 if (extra_atom_set.find(_atoms[i]->pdb_atmnam()) != extra_atom_set.end()) continue;
                 for (unsigned int j = i + 1; j < _atoms.size(); j++) {
                      if (extra_atom_set.find(_atoms[j]->pdb_atmnam()) != extra_atom_set.end()) continue;
                      if (!_atoms[i]->alt_loc().empty() && !_atoms[j]->alt_loc().empty() && _atoms[i]->alt_loc() != _atoms[j]->alt_loc()) continue;
                      if ((_atoms[i]->atom_type() == "H" || _atoms[i]->atom_type() == "D") &&
                          (_atoms[j]->atom_type() == "H" || _atoms[j]->atom_type() == "D")) continue;
                      if (!BondUtil::is_a_bond_with_range(_atoms[i], _atoms[j], distance, lower_limit, upper_limit)) continue;
                      if (bond_set.find(_atoms[i]->pdb_atmnam() + "_" + _atoms[j]->pdb_atmnam()) != bond_set.end()) continue;
                      if ((_atoms[i]->pdb_atmnam() == "N" && (_atoms[j]->pdb_atmnam() == "H1" || _atoms[j]->pdb_atmnam() == "H2" ||
                           _atoms[j]->pdb_atmnam() == "H3" || _atoms[j]->pdb_atmnam() == "D1" || _atoms[j]->pdb_atmnam() == "D2" ||
                           _atoms[j]->pdb_atmnam() == "D3")) || (_atoms[j]->pdb_atmnam() == "N" && (_atoms[i]->pdb_atmnam() == "H1" ||
                           _atoms[i]->pdb_atmnam() == "H2" || _atoms[i]->pdb_atmnam() == "H3" ||_atoms[i]->pdb_atmnam() == "D1" ||
                           _atoms[i]->pdb_atmnam() == "D2" || _atoms[i]->pdb_atmnam() == "D3")) || (_atoms[i]->pdb_atmnam() == "O5'" &&
                          (_atoms[j]->pdb_atmnam() == "HO5'" || _atoms[j]->pdb_atmnam() == "DO5'")) || ((_atoms[i]->pdb_atmnam() == "HO5'" ||
                           _atoms[i]->pdb_atmnam() == "DO5'") && _atoms[j]->pdb_atmnam() == "O5'")) continue;

                      std::string atom_i = _atoms[i]->pdb_atmnam();
                      if (!_atoms[i]->alt_loc().empty()) atom_i += "(" + _atoms[i]->alt_loc() + ")";
                      std::string atom_j = _atoms[j]->pdb_atmnam();
                      if (!_atoms[j]->alt_loc().empty()) atom_j += "(" + _atoms[j]->alt_loc() + ")";
                      // std::string cs = "In '" + ResString + "' residue: incorrect bonding between atom '" + atom_i + "' and atom '" + atom_j + "'.";
                      std::string cs = "In '" + ResString + "' residue: the distance ( " + FloatToString(distance, 0, 2) + " ) between atom '"
                                     + atom_i + "' and atom '" + atom_j + "' is in bond distance range [ " + FloatToString(lower_limit, 0, 2)
                                     + ", " + FloatToString(upper_limit, 0, 2) + " ].";
                      _messageIo->insertMessage("residue_match_print", "warning", cs);
                 }
            }
       } catch (const std::exception& exc) {}
}

bool Residue::update_atom_name(const std::map<std::string, std::string>& atom_mapping, const int& Mol_ID, const bool& index_flag)
{
       bool _successful = true;
       if (atom_mapping.empty()) return _successful;

       std::map<std::string, std::string> new_atom_mapping, uppercase_atom_mapping;
       new_atom_mapping.clear();

       uppercase_atom_mapping.clear();
       for (std::map<std::string, std::string>::const_iterator mpos = atom_mapping.begin(); mpos != atom_mapping.end(); ++mpos) {
            std::string cs = mpos->first;
            String::UpperCase(cs);
            uppercase_atom_mapping.insert(std::make_pair(cs, mpos->second));
       }

       std::set<std::string> extra_atom_set;
       extra_atom_set.clear();
       std::vector<std::string> messages;
       messages.clear();
       for (std::vector<RCSB::Atom*>::const_iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            std::string atomname = (*pos)->atmtype();
            String::UpperCase(atomname);
            std::map<std::string, std::string>::const_iterator mpos = uppercase_atom_mapping.find(atomname);
            if (mpos != uppercase_atom_mapping.end()) {
                 new_atom_mapping.insert(std::make_pair((*pos)->atmtype(), mpos->second));
                 continue;
            }
            if (((*pos)->atom_type() == "D") || ((*pos)->atom_type() == "d")) {
                 atomname[0] = 'H';
                 mpos = uppercase_atom_mapping.find(atomname);
                 if (mpos != uppercase_atom_mapping.end()) {
                      std::string mapped_name = mpos->second;
                      mapped_name[0] = 'D';
                      new_atom_mapping.insert(std::make_pair((*pos)->atmtype(), mapped_name));
                      continue;
                 }
            } else if (((*pos)->atom_type() == "H") || ((*pos)->atom_type() == "h")) {
                 atomname[0] = 'D';
                 mpos = uppercase_atom_mapping.find(atomname);
                 if (mpos != uppercase_atom_mapping.end()) {
                      std::string mapped_name = mpos->second;
                      mapped_name[0] = 'H';
                      new_atom_mapping.insert(std::make_pair((*pos)->atmtype(), mapped_name));
                      continue;
                 }
            }
            extra_atom_set.insert((*pos)->atmtype());
            std::string warning = "Residue (" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + ") in Model " + String::IntToString(Mol_ID)
                                + ": Atom " + (*pos)->atmtype() + " is not mapped.\n";
            messages.push_back(warning);
            if ((*pos)->atom_type() != "H" && (*pos)->atom_type() != "D") _successful = false;
       }

       if (!extra_atom_set.empty()) {
            std::map<std::string, std::string> conflict_atom_mapping;
            conflict_atom_mapping.clear();

            for (std::map<std::string, std::string>::const_iterator pos = new_atom_mapping.begin(); pos != new_atom_mapping.end(); ++pos) {
                 // pos = atom_mapping.begin(); pos != atom_mapping.end(); ++pos) {
                 if (extra_atom_set.find(pos->second) == extra_atom_set.end()) {
                      // new_atom_mapping.insert(std::make_pair(pos->first, pos->second));
                      continue;
                 }
                 conflict_atom_mapping.insert(std::make_pair(pos->first, pos->second));
                 std::string error = "Residue (" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + ") in Model "
                                   + String::IntToString(Mol_ID) + ": mapping Atom " + pos->first + " to " + pos->second + " which already existed.\n";
                 messages.push_back(error);
            }

            if (!_successful) {
                 if (!messages.empty()) {
                      for (std::vector<std::string>::const_iterator
                           pos = messages.begin(); pos != messages.end(); ++pos) {
                           _logIo->message(pos->c_str());
                      }
                 }
                 return _successful;
            }

            if (!conflict_atom_mapping.empty()) {
                 for (std::map<std::string, std::string>::const_iterator pos = conflict_atom_mapping.begin(); pos != conflict_atom_mapping.end(); ++pos) {
                      if (_connected_to_same_atom(pos->first, pos->second)) {
                           new_atom_mapping.insert(std::make_pair(pos->second, pos->second));
                           extra_atom_set.erase(pos->second);
                           extra_atom_set.insert(pos->first);
                      } else new_atom_mapping.insert(std::make_pair(pos->first, pos->second));
                 }
            }

            std::set<std::string> used_atom_names;
            used_atom_names.clear();
            for (std::map<std::string, std::string>::const_iterator pos = new_atom_mapping.begin(); pos != new_atom_mapping.end(); ++pos) {
                 used_atom_names.insert(pos->second);
            }

            int count = 0;
            for (std::set<std::string>::const_iterator pos = extra_atom_set.begin(); pos != extra_atom_set.end(); ++pos) {
                 if (used_atom_names.find(*pos) == used_atom_names.end()) {
                      new_atom_mapping.insert(std::make_pair(*pos, *pos));
                      used_atom_names.insert(*pos);
                      continue;
                 }
                 std::string type = "H";
                 if ((*pos)[0] == 'D') type = "D";
                 while (true) {
                      count++;
                      std::string name = type + String::IntToString(count);
                      if (used_atom_names.find(name) == used_atom_names.end()) {
                           used_atom_names.insert(name);
                           new_atom_mapping.insert(std::make_pair(*pos, name));
                           break;
                      }
                 }
            }
       }

       if (!_successful) {
            if (!messages.empty()) {
                 for (std::vector<std::string>::const_iterator pos = messages.begin(); pos != messages.end(); ++pos) {
                      _logIo->message(pos->c_str());
                 }
            }
            return _successful;
       }

       for (std::map<std::string, std::string>::const_iterator pos = new_atom_mapping.begin(); pos != new_atom_mapping.end(); ++pos) {
            std::string cs = CompositeIndex::getIndex(_ResName, pos->first);
            std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNames.equal_range(cs);
            for (std::multimap<std::string, int>::iterator mpos = range.first; mpos != range.second; ++mpos) {
                 std::string mapped_name = pos->second;
                 if ((_atoms[mpos->second]->atom_type() == "D") && (mapped_name[0] == 'H')) mapped_name[0] = 'D';
                 _atoms[mpos->second]->set_atmtype(mapped_name);
                 _atoms[mpos->second]->set_pdb_atmnam(mapped_name);
            }
       }
       if (index_flag) {
            _atomNameAlts_Orininal.clear();
            _UpdateIndices();
            reorder_atoms();
       }

       return _successful;
}

void Residue::update_partial_atom_name(const std::map<std::string, std::string>& atom_mapping, const int& Mol_ID, const std::string& New_ResName)
{
       if (atom_mapping.empty()) return;

       std::set<std::string> rename_atom_set;
       rename_atom_set.clear();
       for (std::map<std::string, std::string>::const_iterator mpos = atom_mapping.begin(); mpos != atom_mapping.end(); ++mpos) {
            rename_atom_set.insert(mpos->first);
       }

       bool found_old_atom = false;
       bool found_new_atom = false;
       for (std::map<std::string, std::string>::const_iterator mpos = atom_mapping.begin(); mpos != atom_mapping.end(); ++mpos) {
            std::string cs = CompositeIndex::getIndex(_ResName, mpos->first);
            std::multimap<std::string, int>::const_iterator pos = _atomNames.find(cs);
            if (pos == _atomNames.end()) continue;

            found_old_atom = true;

            cs = CompositeIndex::getIndex(_ResName, mpos->second);
            pos = _atomNames.find(cs);
            if (pos == _atomNames.end()) continue;
            if (rename_atom_set.find(mpos->second) != rename_atom_set.end()) continue;

            found_new_atom = true;

            std::string error = "Residue (" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + ") in Model " + String::IntToString(Mol_ID)
                              + ": Map atom '" + mpos->first + "' to a existing atom '" + mpos->second + "'.\n";
            _logIo->message(error.c_str());
       }

       if (!found_old_atom || found_new_atom) return;
            
       if (!New_ResName.empty()) {
            if (_ResName != New_ResName) _OrigResName = _ResName;
            _ResName = New_ResName;
       }

       for (std::vector<RCSB::Atom*>::const_iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (!New_ResName.empty()) {
                 (*pos)->set_restype(_ResName);
                 (*pos)->set_pdb_resnam(_ResName);
            }

            std::map<std::string, std::string>::const_iterator mpos = atom_mapping.find((*pos)->atmtype());
            if (mpos == atom_mapping.end()) continue;

            (*pos)->set_atmtype(mpos->second);
            (*pos)->set_pdb_atmnam(mpos->second);
       }

       _atomNameAlts_Orininal.clear();
       _UpdateIndices();
       reorder_atoms();
}

bool Residue::update_atom_mapping(const int& Mol_ID, const std::string& New_ResName, const std::map<std::string, std::string>& atom_mapping)
{
       std::map<std::string, std::string> new_atom_mapping;
       if (_is_OP123_permutation(New_ResName, atom_mapping, new_atom_mapping)) {
            if (!update_atom_name(new_atom_mapping, Mol_ID, false)) return false;
       } else if (!update_atom_name(atom_mapping, Mol_ID, false)) return false;

       if (_ResName != New_ResName) _OrigResName = _ResName;
       _ResName = New_ResName;

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            (*pos)->set_restype(_ResName);
            (*pos)->set_pdb_resnam(_ResName);
       }

       if (_atomNameAlts_Orininal.empty()) _atomNameAlts_Orininal = _atomNameAlts;
       _UpdateIndices();
       reorder_atoms();

       return true;
}

bool Residue::is_connect(RCSB::Residue* res)
{
       for (std::vector<RCSB::Atom*>::const_iterator pos1 = _atoms.begin(); pos1 != _atoms.end(); ++pos1) {
            for (std::vector<RCSB::Atom*>::const_iterator pos2 = res->_atoms.begin(); pos2 != res->_atoms.end(); ++pos2) {
                 if (!(*pos1)->alt_loc().empty() && !(*pos2)->alt_loc().empty() && (*pos1)->alt_loc() != (*pos2)->alt_loc()) continue;
                 if (BondUtil::is_a_bond(*pos1, *pos2)) return true;
            }
       }
       return false;
}

bool Residue::is_overlap(RCSB::Residue* res)
{
       bool has_multiple_conformation = false;
       if ((_alt_loc_list.size() > 1) && _full_disorder_flag && (res->_alt_loc_list.size() > 1) && res->_full_disorder_flag) has_multiple_conformation = true;

       bool is_single_atom = false;
       if ((_atoms.size() == 1) || (res->_atoms.size() == 1)) is_single_atom = true;

       unsigned int overlap_cutoff = NUM_OVERLAP_CUTOFF + 1;
       if (_atoms.size() < overlap_cutoff) overlap_cutoff = _atoms.size();
       if (res->_atoms.size() < overlap_cutoff) overlap_cutoff = res->_atoms.size();

       unsigned int num_of_overlap = 0;
       for (std::vector<RCSB::Atom*>::const_iterator pos1 = _atoms.begin(); pos1 != _atoms.end(); ++pos1) {
            for (std::vector<RCSB::Atom*>::const_iterator pos2 = res->_atoms.begin(); pos2 != res->_atoms.end(); ++pos2) {
                 // added for entry D_1000019974/1Q6F
                 if (has_multiple_conformation && !(*pos1)->alt_loc().empty() && !(*pos2)->alt_loc().empty() &&
                    (*pos1)->alt_loc() != (*pos2)->alt_loc()) continue;
                 double dist = cal_distance(*pos1, *pos2);
                 if ((dist < 1.3) || (!is_single_atom && BondUtil::is_a_link((*pos1)->atom_type(), (*pos2)->atom_type(), dist))) {
                      num_of_overlap++;
                      if (num_of_overlap >= overlap_cutoff) return true;
                 }
            }
       }

       if (num_of_overlap >= overlap_cutoff) return true;
       return false;
}

bool Residue::_is_OP123_permutation(const std::string& New_ResName, const std::map<std::string, std::string>& atom_mapping,
                                    std::map<std::string, std::string>& new_atom_mapping)
{
       new_atom_mapping.clear();

       std::set<std::string> OP123_set;
       OP123_set.clear();
       OP123_set.insert("OP1");
       OP123_set.insert("OP2");
       OP123_set.insert("OP3");

       std::vector<std::pair<std::string, std::string> > mapping_list;
       mapping_list.clear();
       bool found_mapped_OP3 = false;
       bool found_mismatch = false;
       for (std::set<std::string>::const_iterator spos = OP123_set.begin(); spos != OP123_set.end(); ++spos) {
            std::map<std::string, std::string>::const_iterator mpos = atom_mapping.find(*spos);
            if (mpos == atom_mapping.end() || OP123_set.find(mpos->second) == OP123_set.end()) continue;
            if (mpos->second == "OP3") found_mapped_OP3 = true;
            if (mpos->first != mpos->second) found_mismatch = true;
            mapping_list.push_back(std::make_pair(mpos->first, mpos->second));
       }
       if (mapping_list.empty() || (!found_mapped_OP3 && !found_mismatch)) return false;
       if (found_mapped_OP3 && !found_mismatch && (mapping_list.size() == 3)) return false;

       std::vector<unsigned int> index;
       index.clear();
       std::map<std::string, unsigned int> indexMap;
       indexMap.clear();
       mapping_list.clear();

       for (std::map<std::string, std::string>::const_iterator mpos = atom_mapping.begin(); mpos != atom_mapping.end(); ++mpos) {
            if (!find_atom(mpos->first)) continue;
            if (OP123_set.find(mpos->first) != OP123_set.end() && OP123_set.find(mpos->second) != OP123_set.end()) index.push_back(mapping_list.size());
            indexMap.insert(std::make_pair(mpos->second, mapping_list.size()));
            mapping_list.push_back(std::make_pair(mpos->first, mpos->second));
       }
       if (index.empty()) return false;

       try {
            const ConnectFormat& drug = _ccDic->find_drug(New_ResName);
            const std::map<std::string, std::vector<std::string> >& linked_atom_mapping = drug.getLinkedAtoms();
            std::vector<std::vector<std::string> > linked_atom_list;
            linked_atom_list.clear();
            for (std::vector<unsigned int>::const_iterator ipos = index.begin(); ipos != index.end(); ++ipos) {
                 std::map<std::string, std::vector<std::string> >::const_iterator mpos = linked_atom_mapping.find(mapping_list[*ipos].second);
                 if (mpos != linked_atom_mapping.end()) linked_atom_list.push_back(mpos->second);
            }
            if (index.size() != linked_atom_list.size()) return false;

            std::string P_atom = "";
            for (std::vector<std::string>::const_iterator pos0 = linked_atom_list[0].begin(); pos0 != linked_atom_list[0].end(); ++pos0) {
                 bool found = true;
                 for (unsigned int i = 1; i < linked_atom_list.size(); ++i) {
                      bool found_same_atom = false;
                      for (std::vector<std::string>::const_iterator posi = linked_atom_list[i].begin(); posi != linked_atom_list[i].end(); ++posi) {
                           if (*pos0 == *posi) { found_same_atom = true; break; }
                      }
                      if (!found_same_atom) { found = false; break; }
                 }
                 if (!found) continue;

                 try {
                      const AtomFormat& atom = drug.find_atom(*pos0);
                      if (atom.atomtype() == "P") { P_atom = *pos0; break; }
                 } catch (const std::exception& exc) { }
            }

            if (P_atom.empty()) return false;

            std::map<std::string, unsigned int>::const_iterator ipos = indexMap.find(P_atom);
            if (ipos == indexMap.end()) return false;

            std::map<std::string, std::vector<std::string> >::const_iterator lpos = linked_atom_mapping.find(P_atom);
            if (lpos == linked_atom_mapping.end()) return false;
      
            std::vector<std::pair<std::string, std::string> > atomlist, bondlist, ref_a_list, ref_b_list;
            atomlist.clear();
            bondlist.clear();
            ref_a_list.clear();
            ref_b_list.clear();

            atomlist.push_back(std::make_pair("P", mapping_list[ipos->second].first));
            ref_a_list.push_back(std::make_pair("P", P_atom));

            std::set<std::string> linked_OP_atom_set;
            linked_OP_atom_set.clear();
            for (std::vector<std::string>::const_iterator pos = lpos->second.begin(); pos != lpos->second.end(); ++pos) {
                 if (OP123_set.find(*pos) == OP123_set.end()) continue;
                 linked_OP_atom_set.insert(*pos);
                 std::map<std::string, std::vector<std::string> >::const_iterator opos = linked_atom_mapping.find(*pos);
                 if (opos == linked_atom_mapping.end()) return false;
                 std::string H_atom = "";
                 bool found_other_atom = false;
                 for (std::vector<std::string>::const_iterator pos1 = opos->second.begin(); pos1 != opos->second.end(); ++pos1) {
                      if (*pos1 == P_atom) continue;
                      try {
                           const AtomFormat& atom = drug.find_atom(*pos1);
                           if (atom.atomtype() == "H") H_atom = *pos1;
                           else found_other_atom = true;
                      } catch (const std::exception& exc) { return false; }
                 }
                 if (found_other_atom) continue;

                 ref_a_list.push_back(std::make_pair("O", *pos));
                 ref_b_list.push_back(std::make_pair(P_atom, *pos));
                 if (!H_atom.empty()) {
                      ref_a_list.push_back(std::make_pair("H", H_atom));
                      ref_b_list.push_back(std::make_pair(*pos, H_atom));
                 }
                 std::map<std::string, unsigned int>::const_iterator oipos = indexMap.find(*pos);
                 if (oipos == indexMap.end()) continue;
                 atomlist.push_back(std::make_pair("O", mapping_list[oipos->second].first));
                 bondlist.push_back(std::make_pair(mapping_list[ipos->second].first, mapping_list[oipos->second].first));
                 if (H_atom.empty()) continue;
                 std::map<std::string, unsigned int>::const_iterator hipos = indexMap.find(H_atom);
                 if (hipos == indexMap.end()) continue;
                 if (mapping_list[hipos->second].first[0] == 'D')
                      atomlist.push_back(std::make_pair("D", mapping_list[hipos->second].first));
                 else atomlist.push_back(std::make_pair("H", mapping_list[hipos->second].first));
                 bondlist.push_back(std::make_pair(mapping_list[oipos->second].first, mapping_list[hipos->second].first));
            }
            for (std::vector<unsigned int>::const_iterator ipos = index.begin(); ipos != index.end(); ++ipos) {
                 if (linked_OP_atom_set.find(mapping_list[*ipos].second) == linked_OP_atom_set.end()) return false;
            }
            if (bondlist.empty()) return false;

            bool is_substructure_match = false;
            std::map<std::string, std::string> OP_list_mapping;
            GraphMatch::GetMatch(ref_a_list, ref_b_list, atomlist, bondlist, OP_list_mapping, is_substructure_match, TARGET_TO_REF, false, true);
            if (OP_list_mapping.empty()) return false;

            for (std::vector<std::pair<std::string, std::string> >::const_iterator pos = mapping_list.begin(); pos != mapping_list.end(); ++pos) {
                 std::map<std::string, std::string>::const_iterator mpos = OP_list_mapping.find(pos->first);
                 if (mpos != OP_list_mapping.end())
                      new_atom_mapping.insert(std::make_pair(pos->first, mpos->second));
                 else new_atom_mapping.insert(std::make_pair(pos->first, pos->second));
            }
       } catch (const std::exception& exc) { return false; }

       return true;
}

bool Residue::_connected_to_same_atom(const std::string& first_atom, const std::string& second_atom)
{
       RCSB::Atom* fst = find_atom(first_atom);
       if (!fst) return false;
       RCSB::Atom* fst_connected = _find_connected_atom(fst);
       if (!fst_connected) return false; 
       RCSB::Atom* snd = find_atom(second_atom);
       if (!snd) return false;
       RCSB::Atom* snd_connected = _find_connected_atom(snd);
       if (!snd_connected) return false;
       if (fst_connected == snd_connected) return true;
       return false;
}

RCSB::Atom* Residue::_find_connected_atom(RCSB::Atom* atom)
{
       RCSB::Atom* found_atom = NULL;
       double distance = 100000.0;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (*pos == atom || (*pos)->atom_type() == "H") continue;
            if (atom->alt_loc().empty() || (*pos)->alt_loc().empty() || atom->alt_loc() == (*pos)->alt_loc()) {
                 double dist = cal_distance(atom, *pos);
                 if (dist > 1.4) continue;
                 if (dist < distance) {
                      distance = dist;
                      found_atom = *pos;
                 }
            }
       }
       return found_atom;
}

std::vector<RCSB::Atom*> Residue::_get_terminal_hydrogens(const std::string& atom_name)
{
       std::vector<RCSB::Atom*> heavy_atoms, terminal_h_atoms;
       terminal_h_atoms.clear();

       find_atom(atom_name, heavy_atoms);
       if (!heavy_atoms.empty()) {
            for (unsigned int i = 0; i < heavy_atoms.size(); i++) {
                 for (unsigned int j = 0; j < _atoms.size(); j++) {
                      if (_atoms[j]->atom_type() != "H" && _atoms[j]->atom_type() != "D") continue;
                      if (!heavy_atoms[i]->alt_loc().empty() && !_atoms[j]->alt_loc().empty() && heavy_atoms[i]->alt_loc() != _atoms[j]->alt_loc()) continue;
                      if (!BondUtil::is_a_bond(heavy_atoms[i], _atoms[j])) continue;
                      terminal_h_atoms.push_back(_atoms[j]);
                 }
            }
       }

       return terminal_h_atoms;
}

void Residue::_delete_atoms(const std::set<int>& index_set, std::set<long>& atom_set)
{
       if (index_set.empty()) return;

       std::vector<RCSB::Atom*> atoms;
       atoms.clear();
       int i = 0;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (index_set.find(i) != index_set.end()) {
                 atom_set.insert((long) (*pos));
                 std::string msg = "Deleting atom '" + _pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + " " + (*pos)->pdb_atmnam() + "'.\n";
                 _logIo->message(msg.c_str());
                 delete *pos;
            } else atoms.push_back(*pos);
            i++;
       }
       _atoms = atoms;
}

bool is_connect(const std::string &atm_name_a, RCSB::Residue* res_a, const std::string &atm_name_b, RCSB::Residue* res_b)
{
       std::vector<RCSB::Atom*> atom_a, atom_b;
       res_a->find_atom(atm_name_a, atom_a);
       res_b->find_atom(atm_name_b, atom_b);
       for (std::vector<RCSB::Atom*>::iterator pos1 = atom_a.begin(); pos1 != atom_a.end(); ++pos1) {
            for (std::vector<RCSB::Atom*>::iterator pos2 = atom_b.begin(); pos2 != atom_b.end(); ++pos2) {
                 if (!(*pos1)->alt_loc().empty() && !(*pos2)->alt_loc().empty() && (*pos1)->alt_loc() != (*pos2)->alt_loc()) continue;
                 if (BondUtil::is_a_bond(*pos1, *pos2)) return true;
            }
       }
       return false;
}

bool is_connect(const std::string &atm_name_a, RCSB::Residue* res_a, const std::string &atm_name_b, RCSB::Residue* res_b, const double cut_off)
{
       std::vector<RCSB::Atom*> atom_a, atom_b;
       res_a->find_atom(atm_name_a, atom_a);
       res_b->find_atom(atm_name_b, atom_b);
       for (std::vector<RCSB::Atom*>::iterator pos1 = atom_a.begin(); pos1 != atom_a.end(); ++pos1) {
            for (std::vector<RCSB::Atom*>::iterator pos2 = atom_b.begin(); pos2 != atom_b.end(); ++pos2) {
                 if (!(*pos1)->alt_loc().empty() && !(*pos2)->alt_loc().empty() && (*pos1)->alt_loc() != (*pos2)->alt_loc()) continue;
                 double dist = cal_distance(*pos1, *pos2);
                 if (dist <= cut_off)  return true;
            }
       }
       return false;
}

bool is_connect(RCSB::Residue* res_a, const std::string &atm_type_a, RCSB::Residue* res_b, const std::string &atm_type_b)
{
       std::vector<RCSB::Atom*> atom_a, atom_b;
       res_a->find_atom_by_type(atm_type_a, atom_a);
       res_b->find_atom_by_type(atm_type_b, atom_b);
       for (std::vector<RCSB::Atom*>::iterator pos1 = atom_a.begin(); pos1 != atom_a.end(); ++pos1) {
            for (std::vector<RCSB::Atom*>::iterator pos2 = atom_b.begin(); pos2 != atom_b.end(); ++pos2) {
                 if (!(*pos1)->alt_loc().empty() && !(*pos2)->alt_loc().empty() && (*pos1)->alt_loc() != (*pos2)->alt_loc()) continue;
                 if (BondUtil::is_a_bond(*pos1, *pos2)) return true;
            }
       }
       return false;
}

bool is_glycosidic_link(RCSB::Residue* res_a, RCSB::Residue* res_b)
{
       std::vector<std::pair<std::string, std::string> > atom_type_pair_list;
       atom_type_pair_list.clear();
       atom_type_pair_list.push_back(std::make_pair("C", "O"));
       atom_type_pair_list.push_back(std::make_pair("O", "C"));

       std::vector<RCSB::Atom*> atom_a, atom_b;
       for (std::vector<std::pair<std::string, std::string> >::const_iterator vpos = atom_type_pair_list.begin(); vpos != atom_type_pair_list.end(); ++vpos) {
            res_a->find_atom_by_type(vpos->first, atom_a);
            res_b->find_atom_by_type(vpos->second, atom_b);
            for (std::vector<RCSB::Atom*>::iterator pos1 = atom_a.begin(); pos1 != atom_a.end(); ++pos1) {
                 for (std::vector<RCSB::Atom*>::iterator pos2 = atom_b.begin(); pos2 != atom_b.end(); ++pos2) {
                      if (!(*pos1)->alt_loc().empty() && !(*pos2)->alt_loc().empty() && (*pos1)->alt_loc() != (*pos2)->alt_loc()) continue;
                      if (BondUtil::is_glycosidic_bond(*pos1, *pos2)) return true;
                 }
            }
      }
      return false;
}

bool is_connect(RCSB::Residue* res_a, const std::string &atm_type_a, RCSB::Residue* res_b, const std::string &atm_type_b, std::vector<std::pair<RCSB::Atom*, RCSB::Atom*> >& pair_list)
{
       bool found = false;
       std::vector<RCSB::Atom*> atom_a, atom_b;
       res_a->find_atom_by_type(atm_type_a, atom_a);
       res_b->find_atom_by_type(atm_type_b, atom_b);
       for (std::vector<RCSB::Atom*>::iterator pos1 = atom_a.begin(); pos1 != atom_a.end(); ++pos1) {
            for (std::vector<RCSB::Atom*>::iterator pos2 = atom_b.begin(); pos2 != atom_b.end(); ++pos2) {
                 if (!(*pos1)->alt_loc().empty() && !(*pos2)->alt_loc().empty() && (*pos1)->alt_loc() != (*pos2)->alt_loc()) continue;
                 if (BondUtil::is_a_bond(*pos1, *pos2)) {
                      pair_list.push_back(std::make_pair(*pos1, *pos2));
                      found = true;
                 }
            }
       }
       return found;
}

bool is_polymer_connect(const std::string& chain_type, RCSB::Residue* res1, RCSB::Residue* res2)
{
       std::string F_atom = "C";
       std::string F_atom_alt = "";
       std::string S_atom = "N";
       double dist_cutoff = 1.76;
       if (chain_type == "ATOMN") {
            F_atom = "O3'";
            F_atom_alt = "O2'";
            S_atom = "P";
            dist_cutoff = 1.90;
       }

       if (is_connect(F_atom, res1, S_atom, res2, dist_cutoff)) return true;
       if (!F_atom_alt.empty() && is_connect(F_atom_alt, res1, S_atom, res2, dist_cutoff)) return true;

       return false;
}

void get_proline_n_terminal_hydrogen_mapping_list(std::vector<std::map<std::string, std::vector<std::string> > >& atom_name_mapping_list)
{
       atom_name_mapping_list.clear();

       std::map<std::string, std::vector<std::string> > atom_name_mapping;
       std::vector<std::string> name_list;

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("H2", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("H3", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HT1", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HT2", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HN2", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HN1", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("1HT", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("2HT", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("1H", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H2");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("2H", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("D2", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("D3", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DT1", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DT2", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DN2", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DN1", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("1DT", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("2DT", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("1D", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D2");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("2D", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);
}

void get_other_n_terminal_hydrogen_mapping_list(std::vector<std::map<std::string, std::vector<std::string> > >& atom_name_mapping_list)
{
       atom_name_mapping_list.clear();

       std::map<std::string, std::vector<std::string> > atom_name_mapping;
       std::vector<std::string> name_list;

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H1");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("H1", name_list));
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("H2", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("H3", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H1");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("H", name_list));
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HN", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H1");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HT1", name_list));
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HT2", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HT3", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H1");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("H", name_list));
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("H1", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("H2", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H1");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HN1", name_list));
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HN2", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("HN3", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H1");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("1H", name_list));
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("2H", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("3H", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H1");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("1HT", name_list));
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("2HT", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("3HT", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("H1");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("1HN", name_list));
       name_list.clear();
       name_list.push_back("H2");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("2HN", name_list));
       name_list.clear();
       name_list.push_back("H3");
       name_list.push_back("H1");
       name_list.push_back("H");
       atom_name_mapping.insert(std::make_pair("3HN", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D1");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("D1", name_list));
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("D2", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("D3", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D1");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("D", name_list));
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DN", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D1");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DT1", name_list));
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DT2", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DT3", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D1");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("D", name_list));
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("D1", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("D2", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D1");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DN1", name_list));
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DN2", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("DN3", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D1");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("1D", name_list));
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("2D", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("3D", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D1");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("1DT", name_list));
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("2DT", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("3DT", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);

       atom_name_mapping.clear();
       name_list.clear();
       name_list.push_back("D1");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("1DN", name_list));
       name_list.clear();
       name_list.push_back("D2");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("2DN", name_list));
       name_list.clear();
       name_list.push_back("D3");
       name_list.push_back("D1");
       name_list.push_back("D");
       atom_name_mapping.insert(std::make_pair("3DN", name_list));
       atom_name_mapping_list.push_back(atom_name_mapping);
}
