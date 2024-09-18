/*
FILE:     Molecule.C
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

#include "algorithm-util.h"
#include "algorithm-util.C"
#include "CompositeIndex.h"
#include "Molecule.h"
#include "SeqCodeUtil.h"
#include "SplitUtil.h"
#include "TypeDef.h"
#include "utillib.h"

using namespace RCSB;

#define __NO_STATUS       0
#define __NOT_SURE_STATUS 1
#define __YES_STATUS      2

Molecule::Molecule()
{
       _ccDic = NULL;
       _logIo = NULL;
       _messageIo = NULL;
       clear();
}

void Molecule::clear()
{
       _format_checking = false;
       _residue_reorder = true;
       _is_pdb_input_format = false;
       _found_sugar_flag = false;
       _rename_residue_flag = false;
       _find_sugar_chain_flag = true;
       _Mol_ID = -1;
       _index = 0;
       _chain_index = 0;
       _type = 0;
       _chains.clear();
       _residues.clear();
       _atoms.clear();
       _resIndex.clear();
       _chainOrder.clear();
       _chainIndex.clear();
       _pdbIndex.clear();
       _origPdbIndex.clear();
       _origPdbNumIndex.clear();
       // _cifIndex.clear();
       _origCifIndex.clear();
       _chainIdIndex.clear();
       _asymIdIndex.clear();
       _sugarMergedAsmyIDMap.clear();
       _PDBchainID_changed_map.clear();
       _Branch_Seq_Scheme_Mapping.clear();
       _residue_re_name_mapping.clear();
       _removed_residue_index_set.clear();
       _removed_chain_index_set.clear();
       _order_pos = _chainOrder.end();
}

void Molecule::Reset()
{
       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            if (*pos) delete *pos;
       }
       for (std::vector<RCSB::Residue*>::iterator pos = _residues.begin(); pos != _residues.end(); ++pos) {
            if (*pos) delete *pos;
       }
       clear();
}

const int& Molecule::Mol_ID() const             { return _Mol_ID; }
const int& Molecule::index() const              { return _index; }
const int  Molecule::Num_Chain() const          { return ((int) _chainOrder.size()); }
const int& Molecule::type() const               { return _type; }
// const bool& Molecule::has_sugar_entity() const  { return _found_sugar_flag; }
const std::vector<RCSB::Residue*>& Molecule::Residues() const { return _residues; }
const std::list<RCSB::Atom*>& Molecule::atoms() const  { return _atoms; }
void Molecule::setCCDic(ConnectDic *ccdic)      { _ccDic = ccdic; }
void Molecule::setLog(LogUtil *logPt)           { _logIo = logPt; }
void Molecule::setMessage(MessageUtil *message) { _messageIo = message; }
void Molecule::setFormatChecking()              { _format_checking = true; }
void Molecule::unsetResidueReroder()            { _residue_reorder = false; }
void Molecule::setPDBInputFormat()              { _is_pdb_input_format = true; }
void Molecule::setRenameResidueFlag(const bool& flag) { _rename_residue_flag = flag; }
void Molecule::setFindSugarChainFlag(const bool& flag) { _find_sugar_chain_flag = flag; }
void Molecule::set_Mol_ID(const int& mol_id)    { _Mol_ID = mol_id; }
void Molecule::set_Carbohydrate_Annotation_Info(const std::map<std::string, BRANCH_INFO>& mapping) { _Branch_Seq_Scheme_Mapping = mapping; }
void Molecule::set_index(const int& index_id)   { _index = index_id; }

const bool Molecule::has_sugar_entity()
{
       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMS") return true;
            chain = GetNextChain();
       }
       return false;
}

const bool Molecule::need_sugar_entity_update()
{
       std::set<std::string> used_chain_ids;
       used_chain_ids.clear();
       bool has_sugar_chain = false;

       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if ((chain->chain_type() == "ATOMP") || (chain->chain_type() == "ATOMN")) used_chain_ids.insert(chain->PDB_ChainID());
            else if (chain->chain_type() == "ATOMS") has_sugar_chain = true;
            chain = GetNextChain();
       }
       if (!has_sugar_chain) return false;

       chain = GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMS") {
                 if ((used_chain_ids.find(chain->PDB_ChainID()) != used_chain_ids.end()) || chain->need_carbohydrate_annotation()) return true;
            }
            chain = GetNextChain();
       }

       return false;
}

const int Molecule::MaxNumEntity()
{
       std::map<int, int> entitycount;
       entitycount.clear(); 

       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if (chain->int_entity_id() > 0 && (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN")) {
                 std::map<int, int>::iterator mpos = entitycount.find(chain->int_entity_id());
                 if (mpos != entitycount.end())
                      mpos->second++;
                 else entitycount.insert(std::make_pair(chain->int_entity_id(), 1));
            }
            chain = GetNextChain();
       }
       int _MaxNumEntity = 0;
       for (std::map<int, int>::const_iterator mpos = entitycount.begin(); mpos != entitycount.end(); ++mpos) {
            if (mpos->second > _MaxNumEntity)  _MaxNumEntity = mpos->second;
       }
       return _MaxNumEntity;
}

void Molecule::update_residue_indices()
{
       _pdbIndex.clear();
       // _cifIndex.clear();
       int i = 0;
       for (std::vector<RCSB::Residue*>::iterator pos = _residues.begin(); pos != _residues.end(); ++pos) {
            if ((*pos) != NULL) {
                 std::string cs = CompositeIndex::getIndex((*pos)->pdb_chnid(), (*pos)->ResName(), (*pos)->pdb_res_no(), (*pos)->ins_code());
                 _pdbIndex.insert(std::make_pair(cs, i));
/*
                 cs = CompositeIndex::getIndex((*pos)->chnid(), (*pos)->ResName(), (*pos)->res_no(), "");
                 _cifIndex.insert(std::make_pair(cs, i));
*/
            }
            i++;
       }
}

void  Molecule::insert_a_atom(RCSB::Atom* atom)
{
       _atoms.push_back(atom);
}

bool Molecule::insert_a_atom(const std::vector<std::string>& vals, const bool& convert_flag, const int& lineno)
{
       if (/* vals[0] == "SIGATM" || vals[0] == "SIGUIJ" || */ vals[0] == "TLSU") return true;
       if (vals[0] == "SIGATM" || vals[0] == "ANISOU" || vals[0] == "SIGUIJ" || vals[0] == "TER") {
            if (!_atoms.empty()) {
                 RCSB::Atom* atom = _atoms.back();
                 if (vals[0] == "TER")
                      atom->set_ter_flag(lineno);
                 else {
                      const std::vector<std::string>& record = atom->getValue();
                      for (int i = 2; i < 8; i++) {
                           if (record[i] != vals[i]) return false;
                      }
                      atom->setAuxiliaryValue(vals, convert_flag);
                 }
            } else if (vals[0] != "TER") return false;
            return true;
       }

       RCSB::Atom* atom = new Atom;
       atom->setValue(vals, lineno);
       _atoms.push_back(atom);

       return true;
}

void Molecule::insert_a_residue(RCSB::Residue* residue)
{
       residue->set_index(_residues.size());
       _resIndex.insert(std::make_pair(_residues.size(), _residues.size()));
/*
       std::string cs = CompositeIndex::getIndex(residue->chnid(), residue->ResName(), residue->res_no(), "");
       _cifIndex.insert(std::make_pair(cs, _residues.size()));
*/
       RCSB::Atom* atom = residue->GetFirstAtom();
       // _origPdbIndex only created when the Molecule was biult and never updated afterwards.
       std::string cs = CompositeIndex::getIndex(atom->getOrigValue(17), atom->getOrigValue(16), atom->getOrigValue(15), atom->getOrigValue(7));
       _origPdbIndex.insert(std::make_pair(cs, _residues.size()));
       cs = CompositeIndex::getIndex(atom->getOrigValue(17), atom->getOrigValue(15), atom->getOrigValue(7));
       _origPdbNumIndex.insert(std::make_pair(cs, _residues.size()));

       // _origCifIndex only created when the Molecule was biult and never updated afterwards.
       cs = CompositeIndex::getIndex(atom->getOrigValue(5), atom->getOrigValue(4), atom->getOrigValue(6), "");
       _origCifIndex.insert(std::make_pair(cs, _residues.size()));

       int i;
       if (!IsInteger(residue->pdb_res_no(), i) && _messageIo) {
            std::string error = "Residue ( " + residue->ResName() + " " + residue->pdb_chnid() + " " + residue->pdb_res_no() 
                              + " ) has non-integer residue number.";
            if (_format_checking) _messageIo->insertMessage("error", "model", error, true);
            else _messageIo->insertMessage("residue_number_print", "error", error);
       }
       cs = CompositeIndex::getIndex(residue->pdb_chnid(), residue->ResName(), residue->pdb_res_no(), residue->ins_code());
       _pdbIndex.insert(std::make_pair(cs, _residues.size()));
       _residues.push_back(residue);
}

void Molecule::changed_zero_occupancy()
{
       for (std::list<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (atof((*pos)->occ().c_str()) < 0.001) (*pos)->setValue("1.00", 11);
       }
}

bool Molecule::find_residues(const bool& UNX_flag)
{
       // key: Instance ID
       // pair.first: Reference Comp_ID
       // pair.second: Instance vs. Reference atom mapping;
       std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > > atom_mapping;
       atom_mapping.clear();

       std::set<std::string> instanceid_set;
       instanceid_set.clear();

       return find_residues(false, atom_mapping, instanceid_set, UNX_flag);
}

bool Molecule::find_residues(const bool& split_flag, const std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > >&
                             atom_mapping, std::set<std::string>& instanceid_set, const bool& UNX_flag)
{
       bool _successful = true;

       if (_atoms.empty()) return true;

       std::map<std::string, RCSB::Residue*> residue_mapping;
       residue_mapping.clear();

       std::map<std::string, int> last_number_mapping;
       last_number_mapping.clear();

       std::list<RCSB::Residue*> residue_list;
       residue_list.clear();

       for (std::list<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            _insert_a_atom(*pos, residue_mapping, last_number_mapping, residue_list, split_flag, "");
       }

       residue_mapping.clear();
       _atoms.clear();

       std::vector<std::vector<unsigned int> > special_na_lists;
       special_na_lists.clear();

       std::vector<unsigned int> special_na_list;
       special_na_list.clear();

       _residues.reserve(residue_list.size());
       for (std::list<RCSB::Residue*>::iterator pos = residue_list.begin(); pos != residue_list.end(); ++pos) {
            if (!atom_mapping.empty()) {
                 std::string index = CompositeIndex::getIndex((*pos)->pdb_chnid(), (*pos)->ResName(), (*pos)->pdb_res_no(), (*pos)->ins_code());
                 std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > >::const_iterator mpos = atom_mapping.find(index);
                 if (mpos != atom_mapping.end()) {
                      if ((*pos)->ResName() != mpos->second.first) {
                           // std::string cs = CompositeIndex::getIndex((*pos)->pdb_chnid(), (*pos)->ResName(), (*pos)->pdb_res_no(), (*pos)->ins_code());
                           _pdbIndex.insert(std::make_pair(index, _residues.size()));
                           _residue_re_name_mapping.insert(std::make_pair(index, mpos->second.first));
                      }
      
                      if (!(*pos)->update_atom_mapping(_Mol_ID, mpos->second.first, mpos->second.second)) _successful = false;
                      
                 }
            }
            if (split_flag && _is_split_residue(*pos, last_number_mapping))
                 delete *pos;
            else {
                 if (SeqCodeUtil::is_special_na_residue((*pos)->ResName())) {
                      if (!special_na_list.empty()) {
                           unsigned int last = special_na_list.back();
                           if (_residues[last]->pdb_chnid() != (*pos)->pdb_chnid() || _residues[last]->chnid() != (*pos)->chnid()) {
                                special_na_lists.push_back(special_na_list);
                                special_na_list.clear();
                           }
                      }
                      special_na_list.push_back(_residues.size());
                      if ((*pos)->ter_flag()) {
                           special_na_lists.push_back(special_na_list);
                           special_na_list.clear();
                      }
                 } else (*pos)->change_water_nomenclature();
                 insert_a_residue(*pos);
            }
       }
       residue_list.clear();

       if (!special_na_list.empty()) special_na_lists.push_back(special_na_list);

       if (!special_na_lists.empty()) _update_special_na_residues(special_na_lists);

       for (std::vector<RCSB::Residue*>::iterator rpos = _residues.begin(); rpos != _residues.end(); ++rpos) {
            (*rpos)->set_alt_loc();
            (*rpos)->get_atom_type();
            if (_residue_reorder) (*rpos)->reorder_atoms();
            if (UNX_flag && (*rpos)->ResName() == "UNX") (*rpos)->correction_atom_name(_Mol_ID);
       }

       return _successful;
}

void Molecule::_update_special_na_residues(const std::vector<std::vector<unsigned int> >& idx_lists)
{
       bool flag = false;
       for (std::vector<std::vector<unsigned int> >::const_iterator pos = idx_lists.begin(); pos != idx_lists.end(); ++pos) {
            if (pos->size() < 2) continue;

            bool first = true;
            for (std::vector<unsigned int>::const_iterator ppos = pos->begin(); ppos != pos->end(); ++ppos) {
                 _residues[*ppos]->correction_special_na_name(first, _Mol_ID);
                 first = false;
            }
            flag = true;
       }

       if (flag) update_residue_indices();
}

void Molecule::_insert_a_atom(RCSB::Atom* atom, std::map<std::string, RCSB::Residue*>& residue_mapping, std::map<std::string, int>& last_number_mapping,
                              std::list<RCSB::Residue*>& residue_list, const bool& split_flag, const std::string& polymer_type)
{
       std::string cs = CompositeIndex::getIndex(atom->pdb_chnid(), atom->pdb_resnam(), atom->pdb_resnum(), atom->ins_code());
       RCSB::Residue* residue = NULL;
       std::map<std::string, RCSB::Residue*>::iterator
           mpos = residue_mapping.find(cs);
       if (mpos != residue_mapping.end())
            residue = mpos->second;
       else {
            residue = new Residue;
            residue->setCCDic(_ccDic);
            residue->setLog(_logIo);
            residue->setMessage(_messageIo);
            residue->set_token(atom->type());
            residue->set_ResName(atom->pdb_resnam());
            residue->set_chnid(atom->chnid());
            residue->set_res_no(atom->resnum());
            residue->set_pdb_chnid(atom->pdb_chnid());
            residue->set_pdb_res_no(atom->pdb_resnum());
            residue->set_ins_code(atom->ins_code());
            if (split_flag) {
                 int number = atoi(atom->pdb_resnum().c_str());
                 std::map<std::string, int>::iterator impos = last_number_mapping.find(atom->pdb_chnid());
                 if (impos == last_number_mapping.end())
                      last_number_mapping.insert(std::make_pair(atom->pdb_chnid(), number));
                 else if (number > impos->second)
                      impos->second = number;
            }
            if (SeqCodeUtil::is_water_residue(atom->pdb_resnam()))
                 residue->set_chain_type("HETAS");
            else {
                 std::string type = _ccDic->find_residue_type(atom->pdb_resnam());
                 if ((type == "ATOMP" || type == "ATOMN") && polymer_type == "non-polymer") type = "HETAIN";
                 residue->set_chain_type(type);
            }
            residue_mapping.insert(std::make_pair(cs, residue));
            residue_list.push_back(residue);
       }
       residue->insert_a_atom(atom, _Mol_ID);
}

bool Molecule::_is_split_residue(RCSB::Residue* residue, std::map<std::string, int>& last_number_mapping)
{
       _SPLIT_RESIDUE_LIST *exist = SplitUtil::split_residue(residue->ResName());
       if (!exist) return false;

       std::map<std::string, int>::iterator mpos = last_number_mapping.find(residue->pdb_chnid());
       if (mpos == last_number_mapping.end()) return false;

       int start_number = mpos->second;

       std::vector<RCSB::Residue*> split_list;
       split_list.clear();
       split_list.reserve(exist->num_residue);

       // check if need to run atom name matching first for split residue
 
       unsigned int atom_count = 0;
       std::vector<RCSB::Atom*> atom_list;
       for (int i = 0; i < exist->num_residue; ++i) {
            if (i) start_number++;
            RCSB::Residue* res = new Residue;
            res->setCCDic(_ccDic);
            res->setLog(_logIo);
            res->setMessage(_messageIo);
            res->set_ResName(exist->split_map[i].new_residue);
            res->set_chnid(residue->chnid());
            res->set_pdb_chnid(residue->pdb_chnid());
            if (i) {
                 res->set_res_no(String::IntToString(start_number));
                 res->set_pdb_res_no(String::IntToString(start_number));
                 res->set_ins_code("");
            } else {
                 res->set_res_no(residue->res_no());
                 res->set_pdb_res_no(residue->pdb_res_no());
                 res->set_ins_code(residue->ins_code());
            }
            for (int j = 0; j < exist->split_map[i].num_atom; ++j) {
                 residue->find_atom(exist->split_map[i].atom_map[j].old_atom, atom_list);
                 if (atom_list.empty()) continue;
                 for (std::vector<RCSB::Atom*>::const_iterator
                      apos = atom_list.begin(); apos != atom_list.end(); ++apos) {
                      atom_count++;
                      RCSB::Atom* atom = new Atom;
                      *atom = *(*apos);
                      atom->set_atmtype(exist->split_map[i].atom_map[j].new_atom);
                      atom->set_pdb_atmnam(exist->split_map[i].atom_map[j].new_atom);
                      atom->set_restype(exist->split_map[i].new_residue);
                      atom->set_pdb_resnam(exist->split_map[i].new_residue);
                      if (i) {
                           atom->set_resnum(String::IntToString(start_number));
                           atom->set_pdb_resnum(String::IntToString(start_number));
                      }
                      res->insert_a_atom(atom, _Mol_ID);
                 }
            }
            if (res->NumAtoms())
                 split_list.push_back(res);
            else delete res;
       }

       if (atom_count != residue->NumAtoms()) {
            for (std::vector<RCSB::Residue*>::iterator rpos = split_list.begin(); rpos != split_list.end(); ++rpos) {
                 if (*rpos) delete *rpos;
            }
            split_list.clear();
            return false;
       }

       for (std::vector<RCSB::Residue*>::iterator rpos = split_list.begin(); rpos != split_list.end(); ++rpos) {
            insert_a_residue(*rpos);
       }

       mpos->second = start_number;

       return true;
}

void Molecule::find_chains(const std::set<std::string>& link_residue_set, const bool& find_sugar_link_flag)
{
       std::map<std::string, std::vector<std::vector<std::string> > > scheme_mapping;
       scheme_mapping.clear();

       std::map<std::string, SEQ> seqs;
       seqs.clear();

       find_chains(scheme_mapping, seqs, link_residue_set, false, find_sugar_link_flag);
}

void Molecule::find_chains(const std::map<std::string, std::vector<std::vector<std::string> > >& scheme_mapping, std::map<std::string, SEQ>& seqs,
                           const std::set<std::string>& input_link_residue_set, const bool& keep_solvent_position_flag, const bool& find_sugar_link_flag)
{
       if (_residues.empty()) return;

       std::set<std::string> link_residue_set = input_link_residue_set;
       std::vector<std::set<std::string> > sugar_cluster_list;
       sugar_cluster_list.clear();

       std::vector<std::string> data;

       if (!_Branch_Seq_Scheme_Mapping.empty()) {
            link_residue_set.clear();
            std::vector<std::pair<std::string, std::string> > link_pair_list;
            link_pair_list.clear();
            for (std::set<std::string>::const_iterator lpos = input_link_residue_set.begin(); lpos != input_link_residue_set.end(); ++lpos) {
                 get_wordarray_with_space(data, *lpos, "_");
                 if ((_ccDic->find_residue_type(data[1]) != "ATOMS") || (_ccDic->find_residue_type(data[5]) != "ATOMS")) continue;
                 std::string idx1 = CompositeIndex::getIndex(data[0], data[1], data[2], data[3]);
                 std::map<std::string, std::string>::const_iterator mpos = _residue_re_name_mapping.find(idx1);
                 if (mpos != _residue_re_name_mapping.end()) idx1 = CompositeIndex::getIndex(data[0], mpos->second, data[2], data[3]);
                 std::string idx2 = CompositeIndex::getIndex(data[4], data[5], data[6], data[7]);
                 mpos = _residue_re_name_mapping.find(idx2);
                 if (mpos != _residue_re_name_mapping.end()) idx2 = CompositeIndex::getIndex(data[4], mpos->second, data[6], data[7]);
                 link_pair_list.push_back(std::make_pair(idx1, idx2));
                 link_residue_set.insert(idx1 + "_" + idx2);
            }
            clustering_with_insertion(sugar_cluster_list, link_pair_list);
       }

       std::set<unsigned int> chain_breaks;
       chain_breaks.clear();
       bool has_ter_card = false;
       for (unsigned int i = 0; i < _residues.size(); ++i) {
            if (_residues[i]->ter_flag()) {
                 chain_breaks.insert(i);
                 has_ter_card = true;
            }
            if (!i) continue;

            if (_residues[i-1]->pdb_chnid() != _residues[i]->pdb_chnid() ||
                _residues[i-1]->chnid() != _residues[i]->chnid()) chain_breaks.insert(i-1);

            if ((!SeqCodeUtil::is_water_residue(_residues[i-1]->ResName()) &&
                  SeqCodeUtil::is_water_residue(_residues[i]->ResName())) ||
                 (SeqCodeUtil::is_water_residue(_residues[i-1]->ResName()) &&
                 !SeqCodeUtil::is_water_residue(_residues[i]->ResName()))) chain_breaks.insert(i-1);
       }

       std::vector<std::pair<unsigned int, unsigned int> > pair_array, tmp_array;
       get_chain_start_and_end(chain_breaks, _residues.size(), pair_array);

       std::set<std::string> polymer_index, branch_polymer_chain_id_set;
       polymer_index.clear();

       branch_polymer_chain_id_set.clear();
       for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator pos = pair_array.begin(); pos != pair_array.end(); ++pos) {
            std::string pdb_chnid = _residues[pos->first]->pdb_chnid();
            std::map<std::string, BRANCH_INFO>::const_iterator bmpos = _Branch_Seq_Scheme_Mapping.find(pdb_chnid);
            if (bmpos == _Branch_Seq_Scheme_Mapping.end()) continue;

            unsigned int residue_number = pos->second - pos->first + 1;
            if (bmpos->second.seqs.size() != residue_number) continue;

            bool found_residue_match = true;
            data.clear();
            for (unsigned int i = 0; i < bmpos->second.seqs.size(); ++i) {
                 if ((_residues[pos->first + i]->ResName() != bmpos->second.seqs[i][1]) ||
                     (_residues[pos->first + i]->pdb_res_no() != bmpos->second.seqs[i][2])) {
                      found_residue_match = false;
                 }
                 data.push_back(CompositeIndex::getIndex(_residues[pos->first + i]->pdb_chnid(), _residues[pos->first + i]->ResName(),
                                                         _residues[pos->first + i]->pdb_res_no(), _residues[pos->first + i]->ins_code()));
            }
            if (!found_residue_match) continue;

            found_residue_match = false;
            for (std::vector<std::set<std::string> >::const_iterator vpos = sugar_cluster_list.begin(); vpos != sugar_cluster_list.end(); ++vpos) {
                 if (vpos->size() != data.size()) continue;

                 bool found_residue_match1 = true;
                 for (std::vector<std::string>::const_iterator dpos = data.begin(); dpos != data.end(); ++dpos) {
                      if (vpos->find(*dpos) == vpos->end()) {
                           found_residue_match1 = false;
                           break;
                      }
                 }
                 if (found_residue_match1) {
                      found_residue_match = true;
                      break;
                 }
            }
            if (!found_residue_match) continue;

            std::string idx = pdb_chnid + "_" + String::IntToString(pos->first) + "_" + String::IntToString(pos->second);
            branch_polymer_chain_id_set.insert(idx);
       }

       std::map<std::string, std::string> chain_type_mapping;
       chain_type_mapping.clear();

       // key: PDB Chain ID
       // value: pair array
       std::map<std::string, std::vector<std::pair<unsigned int, unsigned int> > > pair_array_mapping;
       pair_array_mapping.clear();
       for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator pos = pair_array.begin(); pos != pair_array.end(); ++pos) {
            std::string pdb_chnid = _residues[pos->first]->pdb_chnid();
            std::string idx = pdb_chnid + "_" + String::IntToString(pos->first) + "_" + String::IntToString(pos->second);
            if (branch_polymer_chain_id_set.find(idx) != branch_polymer_chain_id_set.end()) continue;

            std::map<std::string, std::vector<std::pair<unsigned int, unsigned int> > >::iterator mpos = pair_array_mapping.find(pdb_chnid);
            if (mpos != pair_array_mapping.end())
                 mpos->second.push_back(*pos);
            else {
                 tmp_array.clear();
                 tmp_array.push_back(*pos);
                 pair_array_mapping.insert(std::make_pair(pdb_chnid, tmp_array));
            }
       }

       std::vector<std::vector<std::string> > mapping;
       for (std::map<std::string, std::vector<std::pair<unsigned int, unsigned int> > >::const_iterator
            mpos = pair_array_mapping.begin(); mpos != pair_array_mapping.end(); ++mpos) {
            mapping.clear();
            std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator mspos = scheme_mapping.find(mpos->first);
            if (mspos != scheme_mapping.end()) {
                 for (std::vector<std::vector<std::string> >::const_iterator vpos = mspos->second.begin(); vpos != mspos->second.end(); ++vpos) {
                      if ((*vpos)[2].empty()) continue;
                      mapping.push_back(*vpos);
                 }
            }
            std::string chain_type_from_seq = "";
            int polymer_status = __NOT_SURE_STATUS;
            std::map<std::string, SEQ>::const_iterator qpos = seqs.find(mpos->first);
            if (qpos != seqs.end()) {
                 polymer_status = __YES_STATUS;
                 chain_type_from_seq = qpos->second.chain_type;
            }
/*
            if (!seqs.empty()) {
                 if (seqs.find(mpos->first) != seqs.end())
                      polymer_status = __YES_STATUS;
                 // else polymer_status = __NO_STATUS;
            }
*/
            _analysis_chain_group(mpos->first, mpos->second, has_ter_card, polymer_status, chain_type_from_seq, mapping,
                                  branch_polymer_chain_id_set, chain_breaks, polymer_index, chain_type_mapping);
       }

       get_chain_start_and_end(chain_breaks, _residues.size(), pair_array);

       std::vector<unsigned int> sugar_residue_set;
       sugar_residue_set.clear();
       if ((_Branch_Seq_Scheme_Mapping.empty() || (_Branch_Seq_Scheme_Mapping.size() != branch_polymer_chain_id_set.size())) && _find_sugar_chain_flag) {
            sugar_residue_set.reserve(pair_array.size());
            for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator pos = pair_array.begin(); pos != pair_array.end(); ++pos) {
                 if ((pos->first == pos->second) && (_residues[pos->first]->chain_type() == "ATOMS")) sugar_residue_set.push_back(pos->first);
            }
       }

       // key: first residue index in the linked cluster sugar residues
       // value: all residue indices in the linked cluster sugar residues
       std::map<unsigned int, std::list<std::list<unsigned int> > > sugar_chains;
       _find_sugar_chains(sugar_chains, sugar_residue_set, link_residue_set, find_sugar_link_flag);

       // key: "_pdbx_poly_seq_scheme.seq_id"_"_pdbx_poly_seq_scheme.mon_id" with "_pdbx_poly_seq_scheme.seq_id" starting from 0 instead of 1
       // value.first: _pdbx_poly_seq_scheme.pdb_seq_num
       // value.second: _pdbx_poly_seq_scheme.pdb_ins_code
       std::map<std::string, std::pair<std::string, std::string> > missing_residue_numbering_mapping;

       std::set<std::string> aligned_with_coordinate_polymer_pdb_chain_ids;
       aligned_with_coordinate_polymer_pdb_chain_ids.clear();

       std::vector<RCSB::Residue*> residue_list;
       std::vector<RCSB::Chain*> possible_with_misplace_residue_chains;
       possible_with_misplace_residue_chains.clear();

       for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator pos = pair_array.begin(); pos != pair_array.end(); ++pos) {
            bool sugar_residue_flag = false;
            if ((pos->first == pos->second) && (_residues[pos->first]->chain_type() == "ATOMS"))
                 sugar_residue_flag = true;

            std::map<unsigned int, std::list<std::list<unsigned int> > >::const_iterator smpos = sugar_chains.find(pos->first);
            if (sugar_residue_flag && !sugar_chains.empty() && smpos == sugar_chains.end()) continue; // only 

            unsigned int first = pos->first;
            // if (sugar_residue_flag && smpos != sugar_chains.end()) first = smpos->second.front();

            RCSB::Chain* chain = new Chain;
            chain->setCCDic(_ccDic);
            chain->setLog(_logIo);
            chain->setMessage(_messageIo);

            RCSB::Atom* atom = _residues[first]->GetFirstAtom();
            chain->set_ChainID(atom->chnid());
            chain->set_PDB_ChainID(atom->pdb_chnid());
            chain->set_PUB_ChainID(atom->pub_chnid());
            insert_a_chain(chain);

            bool found_prev_branch_info_flag = false;
            std::string pdb_chnid = _residues[pos->first]->pdb_chnid();
            std::string idx = pdb_chnid + "_" + String::IntToString(pos->first) + "_" + String::IntToString(pos->second);
            if (branch_polymer_chain_id_set.find(idx) != branch_polymer_chain_id_set.end()) {
                 std::map<std::string, BRANCH_INFO>::const_iterator bmpos = _Branch_Seq_Scheme_Mapping.find(pdb_chnid);
                 if (bmpos != _Branch_Seq_Scheme_Mapping.end()) {
                      if (!bmpos->second.linkages.empty()) chain->set_branch_links(bmpos->second.linkages);
                      if (!bmpos->second.descriptors.empty()) chain->set_linear_descriptors(bmpos->second.descriptors);
                      found_prev_branch_info_flag = true;
                 }
            }

            if (sugar_residue_flag && smpos != sugar_chains.end()) {
                 _found_sugar_flag = true;
                 bool is_polymer_flag = false;
                 if (smpos->second.size() > 1) is_polymer_flag = true;

                 chain->setReserve(smpos->second.size());
                 for (std::list<std::list<unsigned int> >::const_iterator l2pos = smpos->second.begin(); l2pos != smpos->second.end(); ++l2pos) {
                      residue_list.clear();
                      bool insert_as_vector = false;
                      if (is_polymer_flag && l2pos->size() > 1) insert_as_vector = true;
                      for (std::list<unsigned int>::const_iterator lpos = l2pos->begin(); lpos != l2pos->end(); ++lpos) {
                           if (_residues[*lpos]->chnid() != atom->chnid()) {
                                _sugarMergedAsmyIDMap.insert(std::make_pair(_residues[*lpos]->chnid(), atom->chnid()));
                           }
                           if (insert_as_vector) residue_list.push_back(_residues[*lpos]);
                           else chain->insert_a_residue(_residues[*lpos]);
                      }
                      if (!residue_list.empty()) chain->insert_hetero_residues(residue_list);
                 }
                 chain->get_prev_entity_id();
                 if (chain->SeqLen() > 1)
                      chain->set_chain_type("ATOMS");
                 else chain->set_chain_type("HETAIN");
                 continue;
            }

            chain->setReserve(pos->second - pos->first + 1);
            for (unsigned int i = pos->first; i <= pos->second; ++i) {
                 chain->insert_a_residue(_residues[i]);
            }
            chain->get_prev_entity_id();
            if (found_prev_branch_info_flag) {
                 chain->set_chain_type("ATOMS");
                 continue;
            }

            bool is_polymer_chain = false;
            std::string frag_idx = String::IntToString(pos->first) + "_" + String::IntToString(pos->second);
            if (polymer_index.find(frag_idx) != polymer_index.end()) is_polymer_chain = true;

            std::map<std::string, std::string>::const_iterator ctmpos = chain_type_mapping.find(frag_idx);

            std::string chain_type = "";
            if (is_polymer_chain) {
                 std::map<std::string, SEQ>::iterator mpos = seqs.find(chain->PDB_ChainID());
                 if (mpos != seqs.end()) {
                      chain_type = mpos->second.chain_type;
                      chain->set_PolyType(mpos->second.PolyType);
                      chain->set_EvidenceCode(mpos->second.EvidenceCode);
                      mpos->second.never_used = false;
                      bool found_alignment = false;
                      missing_residue_numbering_mapping.clear();
                      std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator mspos = scheme_mapping.find(chain->PDB_ChainID());
                      if (mspos != scheme_mapping.end()) {
                           found_alignment = chain->seq_alignment(mspos->second, _rename_residue_flag);
                           if (keep_solvent_position_flag && !found_alignment) {
                                int i = 0;
                                for (std::vector<std::vector<std::string> >::const_iterator spos = mspos->second.begin(); spos != mspos->second.end(); ++spos) {
                                     if ((*spos)[2].empty()) {
                                          std::string key = String::IntToString(i) + "_" + (*spos)[1];
                                          missing_residue_numbering_mapping.insert(std::make_pair(key, std::make_pair((*spos)[3], (*spos)[4])));
                                     }
                                     i++;
                                }
                           }
                      }
                      if (!found_alignment) chain->seq_alignment(mpos->second, missing_residue_numbering_mapping);

                      if (chain->found_missing_single_residues()) possible_with_misplace_residue_chains.push_back(chain);

                      aligned_with_coordinate_polymer_pdb_chain_ids.insert(chain->PDB_ChainID());
                 } else if (pos->first == pos->second) chain_type = "HETAIN";
                 if (chain_type.empty() && ctmpos != chain_type_mapping.end()) chain_type = ctmpos->second;
            } else {
                 if (ctmpos != chain_type_mapping.end()) chain_type = ctmpos->second;
                 if (chain_type.empty()) {
                      bool is_water_chain = false, has_ATOM_token = false, is_100_percent = false;
                      chain_type = _analysis_chain_type(pos->first, pos->second, has_ATOM_token, is_water_chain, is_100_percent);
                      if (is_water_chain) chain_type = "HETAS";
                      else if ((chain_type == "ATOMP") || (chain_type == "ATOMN") || (chain_type == "ATOMS")) chain_type = "HETAIN";
                 }
            }

            if ((chain->SeqLen() == 1) && (chain_type == "ATOMS")) chain_type = "HETAIN";
            chain->set_chain_type(chain_type);
            chain->set_chain_type_to_residues();
            if (chain_type == "ATOMP" || chain_type == "ATOMN")
                 chain->check_ca_or_p_atom_only();
            if (chain_type == "ATOMN") chain->update_na_type();
       }

       // if (individual_flag) find_sugar_links();
       
       // Insert sequence only empty chains
       for (std::map<std::string, SEQ>::const_iterator mpos = seqs.begin(); mpos != seqs.end(); ++mpos) {
            if (aligned_with_coordinate_polymer_pdb_chain_ids.find(mpos->first) != aligned_with_coordinate_polymer_pdb_chain_ids.end()) continue;
            std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator mspos = scheme_mapping.find(mpos->first);

            RCSB::Chain* chain = new Chain;
            chain->setCCDic(_ccDic);
            chain->setLog(_logIo);
            chain->setMessage(_messageIo);

            chain->set_ChainID(mpos->second.ChainId);
            chain->set_PDB_ChainID(mpos->first);
            // chain->set_entity_id(mpos->second.entity_id);
            chain->set_prev_entity_id(mpos->second.entity_id);
            chain->set_PolyType(mpos->second.PolyType);
            chain->set_chain_type(mpos->second.chain_type);

            if ((mspos != scheme_mapping.end()) && (mspos->second.size() == mpos->second.res.size()))
                 chain->insert_sequence(mspos->second);
            else chain->insert_sequence(mpos->second);
              
            insert_a_chain(chain);
       }

       if (!possible_with_misplace_residue_chains.empty()) {
            // first key: PDB Chain ID
            // second key: Residue Name
            std::map<std::string, std::map<std::string, std::vector<RCSB::Chain*> > > single_residue_chain_mapping;
            single_residue_chain_mapping.clear();

            std::vector<RCSB::Chain*> tmp_vec;
            std::map<std::string, std::vector<RCSB::Chain*> > tmp_map;
            for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
                 if ((*pos)->SeqLen() > 1) continue;
                 RCSB::Residue* res = (*pos)->GetFirstResidue();
                 if (!res) continue;
                 std::map<std::string, std::map<std::string, std::vector<RCSB::Chain*> > >::iterator
                     mpos = single_residue_chain_mapping.find((*pos)->PDB_ChainID());
                 if (mpos != single_residue_chain_mapping.end()) {
                      std::map<std::string, std::vector<RCSB::Chain*> >::iterator mmpos = mpos->second.find(res->ResName());
                      if (mmpos != mpos->second.end()) mmpos->second.push_back(*pos);
                      else {
                           tmp_vec.clear();
                           tmp_vec.push_back(*pos);
                           mpos->second.insert(std::make_pair(res->ResName(), tmp_vec));
                      }
                 } else {
                      tmp_vec.clear();
                      tmp_vec.push_back(*pos);
                      tmp_map.clear();
                      tmp_map.insert(std::make_pair(res->ResName(), tmp_vec));
                      single_residue_chain_mapping.insert(std::make_pair((*pos)->PDB_ChainID(), tmp_map));
                 }
            }

            std::set<int> removed_chain_index_set;
            removed_chain_index_set.clear();
            for (std::vector<RCSB::Chain*>::iterator pos = possible_with_misplace_residue_chains.begin();
                 pos != possible_with_misplace_residue_chains.end(); ++pos) {
                 std::map<std::string, std::map<std::string, std::vector<RCSB::Chain*> > >::const_iterator
                     mpos = single_residue_chain_mapping.find((*pos)->PDB_ChainID());
                 if (mpos == single_residue_chain_mapping.end()) continue;
                 (*pos)->check_missing_misplace_residues(mpos->second, removed_chain_index_set);
            }

            if (!removed_chain_index_set.empty()) {
                 tmp_vec.clear();
                 for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
                      if (*pos == NULL) continue;
                      if (removed_chain_index_set.find((*pos)->index()) != removed_chain_index_set.end()) {
                           delete *pos;
                           continue;
                      }
                      tmp_vec.push_back(*pos);
                 }
                 _chains = tmp_vec;
            }
       }

       InternalOrder();
}

void Molecule::InternalOrder(const bool& numbering_based_order_flag)
{
       std::vector<std::string> chainIds;
       chainIds.clear();

       std::set<std::string> chainId_set;
       chainId_set.clear();

       int pdb_order = 0;

       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            if (*pos == NULL) continue;
            if ((*pos)->chain_type() == "ATOMN" || (*pos)->chain_type() == "ATOMP") {
                 if (chainId_set.find((*pos)->PDB_ChainID()) == chainId_set.end()) {
                      chainIds.push_back((*pos)->PDB_ChainID());
                      chainId_set.insert((*pos)->PDB_ChainID());
                 }
                 (*pos)->set_order(++pdb_order);
            }
       }

       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            if (*pos == NULL) continue;
            if ((*pos)->chain_type() == "ATOMS") {
                 if (chainId_set.find((*pos)->PDB_ChainID()) == chainId_set.end()) {
                      chainIds.push_back((*pos)->PDB_ChainID());
                      chainId_set.insert((*pos)->PDB_ChainID());
                 }
                 (*pos)->set_order(++pdb_order);
            }
       }

       // first key: first residue numner
       // second key: first residue insertion code
       std::map<int, std::map<std::string, std::vector<RCSB::Chain*> > > numbering_based_order_map;
       std::vector<RCSB::Chain*> tmp_chain_list;
       std::map<std::string, std::vector<RCSB::Chain*> > tmp_map;
       for (std::vector<std::string>::const_iterator vpos = chainIds.begin(); vpos != chainIds.end(); ++vpos) {
            numbering_based_order_map.clear();
            for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
                 if (*pos == NULL) continue;
                 if (*vpos == (*pos)->PDB_ChainID() && (((*pos)->chain_type() != "ATOMN" && (*pos)->chain_type() != "ATOMP" &&
                    (*pos)->chain_type() != "ATOMS" && (*pos)->chain_type() != "HETAS") ||
                    ((*pos)->chain_type() == "HETAS" && (*pos)->SeqRes(0)->Field[0] != "HOH" && (*pos)->SeqRes(0)->Field[0] != "DOD"))) {
                      if (numbering_based_order_flag) {
                           RCSB::Residue* res = (*pos)->GetFirstResidue();
                           int number = atoi(res->pdb_res_no().c_str());
                           std::map<int, std::map<std::string, std::vector<RCSB::Chain*> > >::iterator mpos = numbering_based_order_map.find(number);
                           if (mpos != numbering_based_order_map.end()) {
                                std::map<std::string, std::vector<RCSB::Chain*> >::iterator mpos1 = mpos->second.find(res->ins_code());
                                if (mpos1 != mpos->second.end()) mpos1->second.push_back(*pos);
                                else {
                                     tmp_chain_list.clear();
                                     tmp_chain_list.push_back(*pos);
                                     mpos->second.insert(std::make_pair(res->ins_code(), tmp_chain_list));
                                }
                           } else {
                                tmp_chain_list.clear();
                                tmp_chain_list.push_back(*pos);
                                tmp_map.clear();
                                tmp_map.insert(std::make_pair(res->ins_code(), tmp_chain_list));
                                numbering_based_order_map.insert(std::make_pair(number, tmp_map));
                           }
                      } else (*pos)->set_order(++pdb_order);
                 }
            }
            for (std::map<int, std::map<std::string, std::vector<RCSB::Chain*> > >::const_iterator
                 mpos = numbering_based_order_map.begin(); mpos != numbering_based_order_map.end(); ++mpos) {
                 for (std::map<std::string, std::vector<RCSB::Chain*> >::const_iterator mpos1 = mpos->second.begin(); mpos1 != mpos->second.end(); ++mpos1) {
                      for (std::vector<RCSB::Chain*>::const_iterator cpos = mpos1->second.begin(); cpos != mpos1->second.end(); ++cpos) {
                           (*cpos)->set_order(++pdb_order);
                      }
                 }
            }
       }

       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            if (*pos == NULL) continue;
            if (chainId_set.find((*pos)->PDB_ChainID()) != chainId_set.end()) continue;

            if (((*pos)->chain_type() != "ATOMN" && (*pos)->chain_type() != "ATOMP" && (*pos)->chain_type() != "ATOMS" && (*pos)->chain_type() != "HETAS") ||
                ((*pos)->chain_type() == "HETAS" &&  (*pos)->SeqRes(0)->Field[0] != "HOH" && (*pos)->SeqRes(0)->Field[0] != "DOD")) {
                 (*pos)->set_order(++pdb_order);
            }
       }

       for (std::vector<std::string>::const_iterator vpos = chainIds.begin(); vpos != chainIds.end(); ++vpos) {
            for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
                 if (*pos == NULL) continue;
                 if (*vpos == (*pos)->PDB_ChainID() && (*pos)->chain_type() == "HETAS" &&
                    ((*pos)->SeqRes(0)->Field[0] == "HOH" || (*pos)->SeqRes(0)->Field[0] == "DOD")) {
                      (*pos)->set_order(++pdb_order);
                 }
            }
       }

       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            if (*pos == NULL) continue;
            if (chainId_set.find((*pos)->PDB_ChainID()) != chainId_set.end()) continue;

            if ((*pos)->chain_type() == "HETAS" && ((*pos)->SeqRes(0)->Field[0] == "HOH" || (*pos)->SeqRes(0)->Field[0] == "DOD")) {
                 (*pos)->set_order(++pdb_order);
            }
       }

       GenerateInternalOrderIndex();
}

void Molecule::GenerateInternalOrderIndex()
{
       _chainOrder.clear();
       _chainIndex.clear();
       _chainIdIndex.clear();
       _asymIdIndex.clear();
       for (unsigned int i = 0; i < _chains.size(); ++i) {
            if (_chains[i] == NULL) continue;
            _chainOrder.insert(std::make_pair(_chains[i]->order(), i));
            _chainIndex.insert(std::make_pair(_chains[i]->index(), i));
            _chainIdIndex.insert(std::make_pair(_chains[i]->PDB_ChainID(), i));
            _asymIdIndex.insert(std::make_pair(_chains[i]->ChainID(), i));
       }
}

void Molecule::AssignAsymId()
{
       GenerateInternalOrderIndex();

       std::set<std::string> used_ids;
       used_ids.clear();

       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            std::string asym_id = get_next_available_asym_id(used_ids);
            chain->set_ChainID(asym_id);
            chain->update_asymId();
            chain = GetNextChain();
       }

       _asymIdIndex.clear();
       for (unsigned int i = 0; i < _chains.size(); ++i) {
            if (_chains[i] == NULL) continue;
            _asymIdIndex.insert(std::make_pair(_chains[i]->ChainID(), i));
       }
}

void Molecule::CheckUniquePDBNumbering()
{
       std::vector<RCSB::Residue*> residue_list;
       std::map<std::string, std::string> numbering;
       std::set<std::string> tmp_set;

       // Make sure the individual residue's PDB numbering are unique DAOTHER-7590
       // key: PDB ChainID
       // value: ResNum_PDBNumbering_InsCode
       std::map<std::string, std::set<std::string> > unique_numbering;
       unique_numbering.clear();

       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            if (*pos == NULL) continue;
            if (((*pos)->chain_type() == "ATOMN") || ((*pos)->chain_type() == "ATOMP") || ((*pos)->chain_type() == "ATOMS")) {
                 (*pos)->GetFirstResidueList(residue_list);
                 while (!residue_list.empty()) {
                      for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                           std::string res_index = CompositeIndex::getIndex((*rpos)->ResName(), (*rpos)->pdb_res_no(), (*rpos)->ins_code());
                           std::map<std::string, std::set<std::string> >::iterator mpos = unique_numbering.find((*rpos)->pdb_chnid());
                           if (mpos != unique_numbering.end()) mpos->second.insert(res_index);
                           else {
                                tmp_set.clear();
                                tmp_set.insert(res_index);
                                unique_numbering.insert(std::make_pair((*rpos)->pdb_chnid(), tmp_set));
                           }
                      }
                      (*pos)->GetNextResidueList(residue_list);
                 }
            }
       }

       bool found_changed = false;
       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            if (*pos == NULL) continue;
            if (((*pos)->chain_type() == "ATOMN") || ((*pos)->chain_type() == "ATOMP") || ((*pos)->chain_type() == "ATOMS") || 
                ((*pos)->chain_type() == "HETAS" && (((*pos)->SeqRes(0)->Field[0] == "HOH") || ((*pos)->SeqRes(0)->Field[0] == "DOD")))) continue;

            (*pos)->GetFirstResidueList(residue_list);
            while (!residue_list.empty()) {
                 bool found_redundant_number = false;
                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                      std::map<std::string, std::set<std::string> >::iterator mpos = unique_numbering.find((*rpos)->pdb_chnid());
                      if (mpos == unique_numbering.end()) continue;
                      if (mpos->second.find(CompositeIndex::getIndex((*rpos)->ResName(), (*rpos)->pdb_res_no(), (*rpos)->ins_code())) != mpos->second.end()) {
                           found_redundant_number = true;
                           break;
                      }
                 }

                 if (found_redundant_number) {
                      int res_number = 1;
                      while (true) {
                           bool found_unique_number = true;
                           for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                                std::map<std::string, std::set<std::string> >::iterator mpos = unique_numbering.find((*rpos)->pdb_chnid());
                                if (mpos == unique_numbering.end()) continue;
                                if (mpos->second.find(CompositeIndex::getIndex((*rpos)->ResName(), res_number, (*rpos)->ins_code())) != mpos->second.end()) {
                                     found_unique_number = false;
                                     break;
                                }
                           }
                           if (found_unique_number) break;

                           res_number += 1;
                      }

                      numbering.clear();
                      numbering.insert(std::make_pair("number", String::IntToString(res_number)));
                      for (std::vector<RCSB::Residue*>::iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                           (*rpos)->update_nomenclature(numbering);
                      }
                      found_changed = true;
                 }

                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                      std::string res_index = CompositeIndex::getIndex((*rpos)->ResName(), (*rpos)->pdb_res_no(), (*rpos)->ins_code());
                      std::map<std::string, std::set<std::string> >::iterator mpos = unique_numbering.find((*rpos)->pdb_chnid());
                      if (mpos != unique_numbering.end()) mpos->second.insert(res_index);
                      else {
                           tmp_set.clear();
                           tmp_set.insert(res_index);
                           unique_numbering.insert(std::make_pair((*rpos)->pdb_chnid(), tmp_set));
                      }
                 }
                 (*pos)->GetNextResidueList(residue_list);
            }
       }

       if (found_changed) update_residue_indices();
}

int Molecule::ClearNonAsymAndWaterChains()
{
       _chainOrder.clear();
       _chainIndex.clear();
       _chainIdIndex.clear();
       _asymIdIndex.clear();
       for (unsigned int i = 0; i < _chains.size(); ++i) {
            if (_chains[i] == NULL) continue;
            RCSB::Residue* residue = _chains[i]->GetFirstResidue();
            if (!residue) {
                 continue;
            }

            if (residue->ResName() == "HOH" || residue->ResName() == "DOD") {
                 _removed_chain_index_set.insert(_chains[i]->index());
                 delete _chains[i];
                 _chains[i] = NULL;
                 continue;
            }

            _chainOrder.insert(std::make_pair(_chains[i]->order(), i));
            _chainIndex.insert(std::make_pair(_chains[i]->index(), i));
            _chainIdIndex.insert(std::make_pair(_chains[i]->PDB_ChainID(), i));
            _asymIdIndex.insert(std::make_pair(_chains[i]->ChainID(), i));
       }
       return ((int) _chainOrder.size());
}

RCSB::Chain* Molecule::GetFirstChain()
{
       _order_pos = _chainOrder.begin();
       if (_order_pos != _chainOrder.end())
            return (_chains[_order_pos->second]);
       else return NULL;
}

RCSB::Chain* Molecule::GetNextChain()
{
       ++_order_pos;
       if (_order_pos != _chainOrder.end())
            return (_chains[_order_pos->second]);
       else return NULL;
}

RCSB::Chain* Molecule::GetIndexChain(const int& index, bool& is_removed)
{
       is_removed = false;
       if (_removed_chain_index_set.find(index) != _removed_chain_index_set.end()) is_removed = true;
       std::multimap<int, int>::const_iterator pos = _chainIndex.find(index);
       if (pos != _chainIndex.end()) return _chains[pos->second];
       else return NULL;
}

RCSB::Chain* Molecule::GetPolyChain(const std::string& pdb_chnid)
{
       std::multimap<std::string, int>::const_iterator pos = _chainIdIndex.find(pdb_chnid);
       if (pos == _chainIdIndex.end()) return NULL;

       std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator>
       range = _chainIdIndex.equal_range(pdb_chnid);
       for (pos = range.first; pos != range.second; ++pos) {
            if ((_chains[pos->second]->chain_type() == "ATOMP") || (_chains[pos->second]->chain_type() == "ATOMN") ||
                (_chains[pos->second]->chain_type() == "ATOMS"))
                 return _chains[pos->second];
       }

       return NULL;
}

void Molecule::GetPdbChains(const std::string& pdb_chnid, std::vector<RCSB::Chain*>& chain_list)
{
       chain_list.clear();

       std::multimap<std::string, int>::const_iterator pos = _chainIdIndex.find(pdb_chnid);
       if (pos == _chainIdIndex.end()) return;

       std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator>
       range = _chainIdIndex.equal_range(pdb_chnid);
       for (pos = range.first; pos != range.second; ++pos) {
            chain_list.push_back(_chains[pos->second]);
       }
}

RCSB::Chain* Molecule::GetAsymChain(const std::string& asymId)
{
       std::multimap<std::string, int>::const_iterator pos = _asymIdIndex.find(asymId);
       if (pos != _asymIdIndex.end())
            return _chains[pos->second];
       else return NULL;
}

RCSB::Chain* Molecule::GetOrigAsymSugarChain(const std::string& asymId)
{
       std::map<std::string, std::string>::const_iterator mpos = _sugarMergedAsmyIDMap.find(asymId);
       if (mpos != _sugarMergedAsmyIDMap.end()) return GetAsymChain(mpos->second);
       return GetAsymChain(asymId);
}

RCSB::Residue* Molecule::find_pdb_residue(const std::string& chnid, const std::string& resname, const int& resnum, const std::string& ins_code)
{
       std::string cs = CompositeIndex::getIndex(chnid, resname, resnum, ins_code);
       std::multimap<std::string, int>::const_iterator pos = _pdbIndex.find(cs);
       if (pos != _pdbIndex.end())
            return _residues[pos->second];
       else return NULL;
}

RCSB::Residue* Molecule::find_pdb_residue(const std::string& chnid, const std::string& resname, const std::string& resnum, const std::string& ins_code)
{
       std::string cs = CompositeIndex::getIndex(chnid, resname, resnum, ins_code);
       std::multimap<std::string, int>::const_iterator pos = _pdbIndex.find(cs);
       if (pos != _pdbIndex.end())
            return _residues[pos->second];
       else return NULL;
}

RCSB::Residue* Molecule::find_pdb_residue(const std::string& pdb_res_index)
{
       std::multimap<std::string, int>::const_iterator pos = _pdbIndex.find(pdb_res_index);
       if (pos != _pdbIndex.end())
            return _residues[pos->second];
       else return NULL;
}

RCSB::Residue* Molecule::find_orig_pdb_residue(const std::string& chnid, const std::string& resname, const std::string& resnum, const std::string& ins_code)
{
       std::string cs = CompositeIndex::getIndex(chnid, resname, resnum, ins_code);
       std::multimap<std::string, int>::const_iterator pos = _origPdbIndex.find(cs);
       if (pos != _origPdbIndex.end())
            return _residues[pos->second];
       else return NULL;
}

RCSB::Residue* Molecule::find_orig_pdb_residue(const std::string& chnid, const std::string& resnum, const std::string& ins_code)
{
       std::string cs = CompositeIndex::getIndex(chnid, resnum, ins_code);
       std::multimap<std::string, int>::const_iterator pos = _origPdbNumIndex.find(cs);
       if (pos != _origPdbNumIndex.end())
            return _residues[pos->second];
       else return NULL;
}

/*
RCSB::Residue* Molecule::find_cif_residue(const std::string& chnid, const std::string& resname, const int resnum)
{
       std::string cs = CompositeIndex::getIndex(chnid, resname, resnum, "");
       std::multimap<std::string, int>::iterator pos = _cifIndex.find(cs);
       if (pos != _cifIndex.end())
            return _residues[pos->second];
       else return NULL;
}

RCSB::Residue* Molecule::find_cif_residue(const std::string& chnid, const std::string& resname, const std::string& resnum)
{
       std::string cs = CompositeIndex::getIndex(chnid, resname, resnum, "");
       std::multimap<std::string, int>::iterator pos = _cifIndex.find(cs);
       if (pos != _cifIndex.end())
            return _residues[pos->second];
       else return NULL;
}
*/

RCSB::Residue* Molecule::find_orig_cif_residue(const std::string& chnid, const std::string& resname, const std::string& resnum)
{
       std::string cs = CompositeIndex::getIndex(chnid, resname, resnum, "");
       std::multimap<std::string, int>::iterator pos = _origCifIndex.find(cs);
       if (pos != _origCifIndex.end())
            return _residues[pos->second];
       else return NULL;
}

RCSB::Residue* Molecule::find_residue(const int& res_index, bool& is_removed)
{
       is_removed = false;
       if (_removed_residue_index_set.find(res_index) != _removed_residue_index_set.end()) is_removed = true;
       std::map<int, int>::const_iterator mpos = _resIndex.find(res_index);
       if (mpos != _resIndex.end()) return _residues[mpos->second];
       else return NULL;
}

RCSB::Atom* Molecule::find_atom(RCSB::Atom* atom)
{
       if (atom == NULL) return NULL;

       RCSB::Residue* res = find_pdb_residue(atom->pdb_chnid(), atom->pdb_resnam(), atom->pdb_resnum(), atom->ins_code());
       if (res == NULL) return NULL;

       return (res->find_atom(atom->atmtype(), atom->alt_loc()));
}

#if 0
void Molecule::GenResidueIndex()
{
       // _cifIndex.clear();
       _pdbIndex.clear();
       for (std::vector<RCSB::Residue*>::iterator pos = _residues.begin(); pos != _residues.end(); ++pos) {
            RCSB::Residue* residue = *pos;
            RCSB::Atom* atom = residue->GetFirstAtom();
            if (atom == NULL) continue;
/*
            std::string cs = CompositeIndex::getIndex(atom->chnid(), atom->restype(), atom->resnum(), "");
            _cifIndex.insert(std::make_pair(cs, residue->index()));
*/
            std::string cs = CompositeIndex::getIndex(atom->pdb_chnid(), atom->pdb_resnam(), atom->pdb_resnum(), atom->ins_code());
            _pdbIndex.insert(std::make_pair(cs, residue->index()));
       }      
}
#endif
void Molecule::insert_a_chain(RCSB::Chain* chain, const bool& set_index_flag)
{
       int num = _chains.size();
       chain->setCCDic(_ccDic);
       chain->setLog(_logIo);
       chain->setMessage(_messageIo);
       if (set_index_flag) chain->set_index(_chain_index);
       chain->set_order(_chain_index);
       _chains.push_back(chain);
       _chainOrder.insert(std::make_pair(_chain_index, num));
       _chainIndex.insert(std::make_pair(chain->index(), num));
       _chainIdIndex.insert(std::make_pair(chain->PDB_ChainID(), num));
       _asymIdIndex.insert(std::make_pair(chain->ChainID(), num));
       _chain_index++;
}

void Molecule::insert_a_chain(const std::list<RCSB::Residue*>& residue_list)
{
       if (residue_list.empty()) return;

       RCSB::Chain* chain = insert_a_chain(residue_list.front());

       chain->setReserve(residue_list.size());
       for (std::list<RCSB::Residue*>::const_iterator pos = residue_list.begin(); pos != residue_list.end(); ++pos) {
            chain->insert_a_residue(*pos);
       }

       chain->get_prev_entity_id();
}

RCSB::Chain* Molecule::insert_a_chain(RCSB::Residue* res, const bool& set_prev_entity_id_flag)
{
       RCSB::Chain* chain = new Chain;
       chain->setCCDic(_ccDic);
       chain->setLog(_logIo);
       chain->setMessage(_messageIo);

       RCSB::Atom* atom = res->GetFirstAtom();
       chain->set_ChainID(atom->chnid());
       chain->set_PDB_ChainID(atom->pdb_chnid());
       chain->set_PUB_ChainID(atom->pub_chnid());
       chain->set_chain_type(res->chain_type());
       if (set_prev_entity_id_flag) chain->set_prev_entity_id(atom->getValue("label_entity_id"));
       insert_a_chain(chain);

       return chain;
}

void Molecule::build_residue_index()
{
       std::vector<RCSB::Residue*> residue_list;
       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            chain->GetFirstResidueList(residue_list);
            while (!residue_list.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator pos = residue_list.begin(); pos != residue_list.end(); ++pos) {
                      _resIndex.insert(std::make_pair((*pos)->index(), _residues.size()));
                      std::string cs = CompositeIndex::getIndex((*pos)->pdb_chnid(), (*pos)->ResName(), (*pos)->pdb_res_no(), (*pos)->ins_code());
                      _pdbIndex.insert(std::make_pair(cs, _residues.size()));
                      _residues.push_back(*pos);
                 }
                 chain->GetNextResidueList(residue_list);
            }
            chain = GetNextChain();
       }
}

int Molecule::count_atom_no(const std::string& card_id)
{
       int atom_no = 0;
       std::vector<RCSB::Residue*> residue_list;
       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if (chain->chain_type() == card_id || card_id.empty()) {
                 chain->GetFirstResidueList(residue_list);
                 while (!residue_list.empty()) {
                      for (unsigned int i = 0; i < residue_list.size(); ++i) {
                           atom_no += residue_list[i]->NonHydrogenAtomNumbers();
                      }
                      chain->GetNextResidueList(residue_list);
                 }
            }
            chain = GetNextChain();
       }
       return atom_no;
}

int Molecule::count_atom_no(const std::set<std::string>& token_set)
{
       int atom_no = 0;
       std::vector<RCSB::Residue*> residue_list;
       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if (token_set.empty() || token_set.find(chain->chain_type()) != token_set.end()) {
                 chain->GetFirstResidueList(residue_list);
                 while (!residue_list.empty()) {
                      for (unsigned int i = 0; i < residue_list.size(); ++i) {
                           atom_no += residue_list[i]->NonHydrogenAtomNumbers();
                      }
                      chain->GetNextResidueList(residue_list);
                 }
            }
            chain = GetNextChain();
       }
       return atom_no;
}

int Molecule::count_total_atom_number()
{
       int atom_no = 0;
       int chain_number_count = 0;
       std::vector<RCSB::Residue*> residue_list;
       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") chain_number_count++;
            chain->GetFirstResidueList(residue_list);
            while (!residue_list.empty()) {
                 for (unsigned int i = 0; i < residue_list.size(); ++i) {
                      atom_no += residue_list[i]->AtomNumbers();
                 }
                 chain->GetNextResidueList(residue_list);
            }
            chain = GetNextChain();
       }
       return (atom_no + chain_number_count);
}

void Molecule::update_cif_nomenclature()
{
       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            int serial_no = 1;
            (*pos)->update_residues_nomenclature(true, false, serial_no, false);
       }
}

void Molecule::assign_carbohydrate_chain_id(std::map<std::string, std::string>& sugar_chain_id_mapping)
{
       sugar_chain_id_mapping.clear();

       std::set<std::string> used_pdb_chain_ids;
       used_pdb_chain_ids.clear();

       bool has_carbohydrate_chain = false;
       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") used_pdb_chain_ids.insert(chain->PDB_ChainID());
            if (chain->chain_type() == "ATOMS") has_carbohydrate_chain = true;
            chain = GetNextChain();
       }
       if (!has_carbohydrate_chain) return;

       bool changed = false;
       chain = GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMS") {
                 std::string old_chain_id = chain->PDB_ChainID();
                 std::string new_chain_id = chain->PDB_ChainID();
                 if (used_pdb_chain_ids.find(chain->PDB_ChainID()) == used_pdb_chain_ids.end()) used_pdb_chain_ids.insert(chain->PDB_ChainID());
                 else {
                      new_chain_id = get_next_available_pdb_chain_id(used_pdb_chain_ids);
                      // chain->set_PDB_ChainID_to_residues(new_chain_id);
                      chain->set_PDB_ChainID(new_chain_id);
                 }

                 std::map<std::string, std::set<std::string> >::const_iterator mpos = _PDBchainID_changed_map.find(old_chain_id);
                 if (mpos != _PDBchainID_changed_map.end()) {
                      for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                           sugar_chain_id_mapping.insert(std::make_pair(new_chain_id, *spos));
                      }
                 } else if (old_chain_id != new_chain_id) sugar_chain_id_mapping.insert(std::make_pair(new_chain_id, old_chain_id));

                 chain->renumbering_polymer_chain_starting_with_one();
                 chain->update_indices();
                 changed = true;
            }
            chain = GetNextChain();
       }

       if (changed) {
            update_residue_indices();
            GenerateInternalOrderIndex();
       }
}

void Molecule::merged_linked_sugar_chains(const std::set<std::string>& residue_set, const std::set<std::string>& link_set)
{
       if (_chains.empty() || residue_set.empty() || link_set.empty()) return;

       std::set<std::string> tmp_set;
       std::map<int, std::set<std::string> > residue_idx_map;
       residue_idx_map.clear();
       std::vector<int> index;
       index.clear();
       for (unsigned int i = 0; i < _chains.size(); ++i) {
            if ((_chains[i]->chain_type() == "ATOMP") || (_chains[i]->chain_type() == "ATOMN")) continue;
            RCSB::Residue* res = _chains[i]->GetFirstResidue();
            if ((res->ResName() == "HOH") || (res->ResName() == "DOD")) continue;

            tmp_set.clear();
            bool found = false;
            while (res) {
                 std::string idx = CompositeIndex::getIndex(res->pdb_chnid(), res->ResName(), res->pdb_res_no(), res->ins_code());
                 if (residue_set.find(idx) != residue_set.end()) found = true;
                 tmp_set.insert(idx);
                 res = _chains[i]->GetNextResidue();
            }
            if (!found) continue;

            index.push_back(i);
            residue_idx_map.insert(std::make_pair(i, tmp_set));
       }
       if (index.size() < 2) return;

       std::vector<RCSB::Residue*> residue_list;
       std::vector<std::map<std::string, std::string> > empty_vec;
       empty_vec.clear();
       std::set<int> removed_chains;
       removed_chains.clear();

       for (unsigned int i = 0; i < index.size() - 1; ++i) {
            if (index[i] < 0) continue;
            std::map<int, std::set<std::string> >::const_iterator impos = residue_idx_map.find(index[i]);
            if (impos == residue_idx_map.end()) continue;

            for (unsigned int j = i + 1; j < index.size(); ++j) {
                 if (index[j] < 0) continue;
                 std::map<int, std::set<std::string> >::const_iterator jmpos = residue_idx_map.find(index[j]);
                 if (jmpos == residue_idx_map.end()) continue;

                 bool found = false;
                 for (std::set<std::string>::const_iterator ispos = impos->second.begin(); ispos != impos->second.end(); ++ispos) {
                      for (std::set<std::string>::const_iterator jspos = jmpos->second.begin(); jspos != jmpos->second.end(); ++jspos) {
                           if (link_set.find(*ispos + "_" + *jspos) != link_set.end()) {
                                found = true;
                                break;
                           }
                      }
                      if (found) break;
                 }
                 if (!found) continue;

                 _chains[index[j]]->GetFirstResidueList(residue_list);
                 while (!residue_list.empty()) {
                      for (unsigned int k = 0; k < residue_list.size(); ++k) {
                           _chains[index[i]]->insert_a_residue(residue_list[k]);
                      }
                      _chains[index[j]]->GetNextResidueList(residue_list);
                 }
                 removed_chains.insert(index[j]);
                 index[j] = -1;

                 _chains[index[i]]->set_chain_type("ATOMS");
                 _chains[index[i]]->get_prev_entity_id();
                 _chains[index[i]]->set_linear_descriptors(empty_vec);
                 _chains[index[i]]->set_branch_links(empty_vec);
            }
       }

       _remove_empty_chains(removed_chains);
}

void Molecule::_find_sugar_chains(std::map<unsigned int, std::list<std::list<unsigned int> > >& sugar_chains, const std::vector<unsigned int>&
                                  sugar_residue_set, const std::set<std::string>& link_residue_set, const bool& find_sugar_link_flag)
{
       sugar_chains.clear();
       if (sugar_residue_set.empty()) return;

       std::set<unsigned int> tmp_set, include_set;

       std::vector<std::set<unsigned int> > groups;
       groups.clear();
       groups.reserve(sugar_residue_set.size());
       for (unsigned int i = 0; i < sugar_residue_set.size(); ++i) {
            tmp_set.clear();
            tmp_set.insert(sugar_residue_set[i]);
            groups.push_back(tmp_set);
       }
       std::vector<std::pair<unsigned int, unsigned int> > links;
       links.clear();

       std::map<unsigned int, float> residue_occupancy_map;
       residue_occupancy_map.clear();
       for (unsigned int i = 0; i < sugar_residue_set.size(); ++i) {
            residue_occupancy_map.insert(std::make_pair(i, _residues[sugar_residue_set[i]]->get_partial_occupancy()));
       }

       // key: residue index
       // value: overlaped residue index set
       std::map<unsigned int, std::set<unsigned int> > overlap_index_map, tmp_index_mapping;
       overlap_index_map.clear();

       for (unsigned int i = 0; i < sugar_residue_set.size() - 1; ++i) {
            if (overlap_index_map.find(_residues[sugar_residue_set[i]]->index()) != overlap_index_map.end()) continue;

            for (unsigned int j = i + 1; j < sugar_residue_set.size(); ++j) {
                 if (overlap_index_map.find(_residues[sugar_residue_set[j]]->index()) != overlap_index_map.end()) continue;

                 if (((residue_occupancy_map[i] + residue_occupancy_map[j]) < 1.02) &&
                      _residues[sugar_residue_set[i]]->is_overlap(_residues[sugar_residue_set[j]])) {
                      _insert_overlap_residue_index_map(sugar_residue_set[i], sugar_residue_set[j], overlap_index_map);
                      _insert_overlap_residue_index_map(sugar_residue_set[j], sugar_residue_set[i], overlap_index_map);
                      continue;
                 }
                 std::string idx = _residues[sugar_residue_set[i]]->pdb_chnid() + "_" + _residues[sugar_residue_set[i]]->ResName() + "_"
                                 + _residues[sugar_residue_set[i]]->pdb_res_no() + "_" + _residues[sugar_residue_set[i]]->ins_code() + "_"
                                 + _residues[sugar_residue_set[j]]->pdb_chnid() + "_" + _residues[sugar_residue_set[j]]->ResName() + "_"
                                 + _residues[sugar_residue_set[j]]->pdb_res_no() + "_" + _residues[sugar_residue_set[j]]->ins_code();
                 if (link_residue_set.find(idx) != link_residue_set.end())
                      links.push_back(std::make_pair(sugar_residue_set[i], sugar_residue_set[j]));
                 else if (find_sugar_link_flag && is_glycosidic_link(_residues[sugar_residue_set[i]], _residues[sugar_residue_set[j]])) {
                      links.push_back(std::make_pair(sugar_residue_set[i], sugar_residue_set[j]));
                 }
            }
       }
       if (!links.empty()) clustering_with_merging(groups, links);

       std::vector<unsigned int> two_residue_vec;

       std::list<std::list<unsigned int> > tmp_list_list;
       std::list<unsigned int> tmp_list;
       for (std::vector<std::set<unsigned int> >::const_iterator vpos = groups.begin(); vpos != groups.end(); ++vpos) {
            tmp_index_mapping.clear();
            include_set.clear();

            bool found_overlaps = false;
            std::set<unsigned int>::const_iterator spos = vpos->begin();
            unsigned int key = *spos;
            for ( ; spos != vpos->end(); ++spos) {
                 if (include_set.find(*spos) != include_set.end()) continue;

                 tmp_set.clear();
                 tmp_set.insert(*spos);

                 // merge overlap residues as single residue list
                 std::map<unsigned int, std::set<unsigned int> >::const_iterator mpos = overlap_index_map.find(*spos);
                 if (mpos != overlap_index_map.end()) {
                      for (std::set<unsigned int>::const_iterator smpos = mpos->second.begin(); smpos != mpos->second.end(); ++smpos) {
                           if (vpos->find(*smpos) == vpos->end()) continue;
                           tmp_set.insert(*smpos);
                           include_set.insert(*smpos);
                           found_overlaps = true;
                      }
                 }
                 tmp_index_mapping.insert(std::make_pair(*spos, tmp_set));
            }

            tmp_list_list.clear();
            two_residue_vec.clear();
            for (std::map<unsigned int, std::set<unsigned int> >::const_iterator mpos = tmp_index_mapping.begin(); mpos != tmp_index_mapping.end(); ++mpos) {
                 tmp_list.clear();
                 for (std::set<unsigned int>::const_iterator spos = mpos->second.begin(); spos !=  mpos->second.end(); ++spos) {
                      tmp_list.push_back(*spos);
                      if (!found_overlaps) two_residue_vec.push_back(*spos);
                 }
                 tmp_list_list.push_back(tmp_list);
            }
            if (two_residue_vec.size() == 2) {
                 if (_residues[two_residue_vec[0]]->is_overlap(_residues[two_residue_vec[1]])) {
                      for (std::vector<unsigned int>::const_iterator trvpos = two_residue_vec.begin(); trvpos != two_residue_vec.end(); ++trvpos) {
                           tmp_list.clear();
                           tmp_list.push_back(*trvpos);
                           tmp_list_list.clear();
                           tmp_list_list.push_back(tmp_list);
                           sugar_chains.insert(std::make_pair(*trvpos, tmp_list_list));
                      }
                 } else sugar_chains.insert(std::make_pair(key, tmp_list_list));
            } else sugar_chains.insert(std::make_pair(key, tmp_list_list));
       }
}

/*
void Molecule::find_sugar_links()
{
       if (_chains.empty()) return;
       std::vector<int> index;
       index.clear();
       int count = 0; 
       for (std::vector<RCSB::Chain*>::const_iterator
            pos = _chains.begin(); pos != _chains.end(); ++pos) {
            if ((*pos)->chain_type() == "ATOMS") {
                 index.push_back(count);
            }
            count++;
       }
       if (index.size() < 2) return;

       std::set<int> removed_chains;
       removed_chains.clear();
       for (unsigned int i = 0; i < index.size() - 1; ++i) {
            if (index[i] < 0) continue;
            for (unsigned int j = i + 1; j < index.size(); ++j) {
                 if (index[j] < 0) continue;
                 bool found = false;
                 RCSB::Residue* residue1 = _chains[index[i]]->GetFirstResidue();
                 while (residue1) {
                      RCSB::Residue* residue2 = _chains[index[j]]->GetFirstResidue();
                      while (residue2) {
                           if (is_connect(residue1, "O", residue2, "C")) {
                                found = true; break;
                           }
                           residue2 = _chains[index[j]]->GetNextResidue();
                      }
                      if (found) break;
                      residue1 = _chains[index[i]]->GetNextResidue();
                 }
                 if (found) {
                      std::vector<RCSB::Residue*> residue_list;
                      _chains[index[j]]->GetFirstResidueList(residue_list);
                      while (!residue_list.empty()) {
                           for (unsigned int k = 0; k < residue_list.size(); ++k) {
                                _chains[index[i]]->insert_a_residue(residue_list[k]);
                           }
                           _chains[index[j]]->GetNextResidueList(residue_list);
                      }
                      removed_chains.insert(index[j]);
                      index[j] = -1;

                      _chains[index[i]]->get_prev_entity_id();
                 }
            }
       }

       _remove_empty_chains(removed_chains);
}
*/

void Molecule::_remove_empty_chains(const std::set<int>& chain_set)
{
       if (chain_set.empty()) return;

       std::vector<RCSB::Chain*> tmp_chain;
       tmp_chain.clear();
       int count = 0; 
       for (std::vector<RCSB::Chain*>::const_iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            RCSB::Chain* chain = *pos;
            if (chain_set.find(count) != chain_set.end()) delete chain;
            else tmp_chain.push_back(chain);
            count++;
       }

       _chains = tmp_chain;
       GenerateInternalOrderIndex();
}

void Molecule::_analysis_chain_group(const std::string PDB_Chain_ID, const std::vector<std::pair<unsigned int, unsigned int> >& pair_array,
                                     const bool& has_ter_card, const int& polymer_status, const std::string& chain_type_from_seq,
                                     const std::vector<std::vector<std::string> >& mapping, const std::set<std::string>& branch_chain_id_set,
                                     std::set<unsigned int>& chain_breaks, std::set<std::string>& polymer_index,
                                     std::map<std::string, std::string>& chain_type_mapping)
{
       std::set<int> ATOM_chain_set, polymer_to_nonpolymer_set;
       ATOM_chain_set.clear();
       polymer_to_nonpolymer_set.clear();

       std::map<int, std::string> type_mapping;
       type_mapping.clear();

       std::vector<int> polymer_chains, nonpolymer_chains;
       polymer_chains.clear();
       polymer_chains.reserve(pair_array.size());
       nonpolymer_chains.clear();
       nonpolymer_chains.reserve(pair_array.size());
       std::set<int> polymer_chain_100_percent_set;
       polymer_chain_100_percent_set.clear();

       bool has_ATOM_token = false, is_water_chain = false, is_100_percent = false;
       unsigned int idx = 0;
       for (std::vector<std::pair<unsigned int, unsigned int> >::const_iterator pos = pair_array.begin(); pos != pair_array.end(); ++pos) {
            std::string branch_idx = _residues[pos->first]->pdb_chnid() + "_" + String::IntToString(pos->first) + "_" + String::IntToString(pos->second);
            std::string chain_type = "";
            if (branch_chain_id_set.find(branch_idx) != branch_chain_id_set.end()) chain_type = "ATOMS";
            else chain_type = _analysis_chain_type(pos->first, pos->second, has_ATOM_token, is_water_chain, is_100_percent);
            if (has_ATOM_token) ATOM_chain_set.insert(idx);

            bool found_match_pdbx_poly_seq_scheme = false;
            if (mapping.size() == (pos->second - pos->first + 1)) {
                 unsigned int count = 0;
                 for (unsigned int i = 0; i < mapping.size(); ++i) {
                      if (mapping[i][2] != _residues[pos->first + i]->ResName() || mapping[i][3] != _residues[pos->first + i]->pdb_res_no() ||
                          mapping[i][4] != _residues[pos->first + i]->ins_code()) break;
                      count++;
                 }
                 if (count == mapping.size()) found_match_pdbx_poly_seq_scheme = true;
            } 

            if (found_match_pdbx_poly_seq_scheme) {
                 polymer_chains.push_back(idx);
                 if (!chain_type_from_seq.empty()) type_mapping.insert(std::make_pair(idx, chain_type_from_seq));
            } else if (is_water_chain) {
                 std::string frag_idx = String::IntToString(pos->first) + "_" + String::IntToString(pos->second);
                 chain_type = "HETAS";
                 chain_type_mapping.insert(std::make_pair(frag_idx, chain_type));
            } else {
                 type_mapping.insert(std::make_pair(idx, chain_type));
                 if (branch_chain_id_set.find(branch_idx) == branch_chain_id_set.end()) {
                      if (polymer_status == __NO_STATUS) nonpolymer_chains.push_back(idx);
                      else {
                           if (chain_type == "ATOMP" || chain_type == "ATOMN") {
                                polymer_chains.push_back(idx);
                                if (is_100_percent) polymer_chain_100_percent_set.insert(idx);
                           } else nonpolymer_chains.push_back(idx);
                      }
                 }
            }
            idx++;
       }

       if (!polymer_chains.empty()) {
            // currently select polymer based on chain length
            // may use sequence information in future
            idx = 0;
            int max_length = pair_array[polymer_chains[0]].second - pair_array[polymer_chains[0]].first;
            for (unsigned int i = 1; i < polymer_chains.size(); ++i) {
                 int length = pair_array[polymer_chains[i]].second - pair_array[polymer_chains[i]].first;
                 if (length > max_length) {
                      max_length = length;
                      idx = i;
                 }
            }

            std::set<int> ter_set;
            ter_set.clear();
            for (unsigned int i = 0; i < polymer_chains.size(); ++i) {
                 if (_residues[pair_array[polymer_chains[i]].second]->ter_flag()) ter_set.insert(_residues[pair_array[polymer_chains[i]].second]->ter_flag());
                 if (i == idx) continue;
                 nonpolymer_chains.push_back(polymer_chains[i]);
                 polymer_to_nonpolymer_set.insert(polymer_chains[i]);
            }
            if (_format_checking && _messageIo) {
                 if (_is_pdb_input_format) {
                      if (ter_set.size() > 1) {
                           std::string lines = "";
                           for (std::set<int>::const_iterator spos = ter_set.begin(); spos != ter_set.end(); ++spos) {
                                if (!lines.empty()) lines += ", ";
                                lines += String::IntToString(*spos);
                           }
                           std::string error = "There are multiple TER cards for polymer chain '" + PDB_Chain_ID + "' at lines ("
                                   + lines + "). Only polymer chain should have a TER card at the end. No TER cards should be";
                           error += " included at the end of non-polymer residues. Please fix and upload the file again.";
                           _messageIo->insertMessage("error", "model", error, true);
                      } else if (ter_set.empty() && max_length > 20) {
                           std::string error = "There is no TER card at the end of polymer chain '" + PDB_Chain_ID
                                + "'. All polymer chains should have a TER card at the end of each chain. Please add TER card and restart again.";
                           _messageIo->insertMessage("error", "model", error, true);
                      }
                 } else {
                      int last_i = (int) polymer_chains.size() - 1;
                      for (int i = (int) polymer_chains.size() - 1; i >= 0; --i) {
                           last_i = i;
                           int length = pair_array[polymer_chains[i]].second - pair_array[polymer_chains[i]].first + 1;
                           if ((i == ((int) polymer_chains.size() - 1) && length > 5) || (i != ((int) polymer_chains.size() - 1) && length > 1)) break;
                      }

                      std::set<std::string> tmp_set;
                      tmp_set.clear();
                      std::vector<std::string> label_asym_ids;
                      label_asym_ids.clear();
                      for (int i = 0; i <= last_i; ++i) {
                           if (tmp_set.find(_residues[pair_array[polymer_chains[i]].first]->chnid()) != tmp_set.end()) continue;

                           if ((i == (int) idx) || (polymer_chain_100_percent_set.find(i) != polymer_chain_100_percent_set.end()))  {
                                tmp_set.insert(_residues[pair_array[polymer_chains[i]].first]->chnid());
                                label_asym_ids.push_back(_residues[pair_array[polymer_chains[i]].first]->chnid());
                           }
                      }

                      if (label_asym_ids.size() > 1) {
                           std::string error = "Polymer chain '" + PDB_Chain_ID + "' is separated into multiple fragments with '_atom_site.label_asym_id' values [ ";
                           std::string asym_id = "";
                           for (std::vector<std::string>::const_iterator pos = label_asym_ids.begin(); pos != label_asym_ids.end(); ++pos) {
                                if (!asym_id.empty()) asym_id += ", ";
                                asym_id += "'" + *pos + "'";
                           }
                           error += asym_id + " ]. Please put them together with same '_atom_site.label_asym_id' value.";
                           _messageIo->insertMessage("error", "model", error, true);
                      }
                 }
            }

            unsigned int first = pair_array[polymer_chains[idx]].first;
            unsigned int second = pair_array[polymer_chains[idx]].second;
            bool found_polymer = true;

            if (polymer_status == __NOT_SURE_STATUS) {
                 // for case like D_123812
                 if (has_ter_card && !_residues[second]->ter_flag() &&
                     ATOM_chain_set.find(polymer_chains[idx]) == ATOM_chain_set.end())
                      found_polymer = false;
                 if (found_polymer) {
                      std::map<int, std::string>::const_iterator mpos = type_mapping.find(polymer_chains[idx]);
                      if (mpos != type_mapping.end() && ((second - first) < 15) && !_residues[second]->ter_flag())
                           found_polymer = _check_linkage_between_residues(mpos->second, first, second);
                  }
            }

            if (found_polymer) {
                 std::string frag_idx = String::IntToString(first) + "_" + String::IntToString(second);
                 std::map<int, std::string>::const_iterator mpos = type_mapping.find(polymer_chains[idx]);
                 if (mpos != type_mapping.end()) chain_type_mapping.insert(std::make_pair(frag_idx, mpos->second));
                 polymer_index.insert(frag_idx);
            } else {
                 nonpolymer_chains.push_back(polymer_chains[idx]);
                 polymer_to_nonpolymer_set.insert(polymer_chains[idx]);
            }
       }

       if (nonpolymer_chains.empty()) return;

       for (std::vector<int>::const_iterator pos = nonpolymer_chains.begin(); pos != nonpolymer_chains.end(); ++pos) {
            if (_format_checking && _messageIo && _is_pdb_input_format && _residues[pair_array[*pos].second]->ter_flag() &&
                polymer_to_nonpolymer_set.find(*pos) == polymer_to_nonpolymer_set.end()) {
                 // std::string error = "All polymer chains should have a TER card at the end. No TER cards should be ";
                 std::string error = "Please remove 'TER' card after single residue '" + _residues[pair_array[*pos].second]->ResName() + " "
                        + _residues[pair_array[*pos].second]->pdb_chnid() + " " + _residues[pair_array[*pos].second]->pdb_res_no()
                        + "' in line " + String::IntToString(_residues[pair_array[*pos].second]->ter_flag())
                        + " and upload the file again.";
                 _messageIo->insertMessage("error", "model", error);
            }
            if (pair_array[*pos].first == pair_array[*pos].second) {
                 std::map<int, std::string>::const_iterator mpos = type_mapping.find(*pos);
                 if (mpos != type_mapping.end()) {
                      std::string chain_type = mpos->second;
                      std::string frag_idx = String::IntToString(pair_array[*pos].first) + "_"
                                           + String::IntToString(pair_array[*pos].second);
                      if (chain_type == "ATOMP" || chain_type == "ATOMN") chain_type = "HETAIN";
                      chain_type_mapping.insert(std::make_pair(frag_idx, chain_type));
                 }
            } else {
                 for (unsigned int i = pair_array[*pos].first; i < pair_array[*pos].second; ++i) {
                      if (!SeqCodeUtil::is_water_residue(_residues[i]->ResName()) ||
                          !SeqCodeUtil::is_water_residue(_residues[i+1]->ResName()))
                           chain_breaks.insert(i);
                 }
            }
       }
}

std::string Molecule::_analysis_chain_type(const unsigned int& start, const unsigned int& end, bool& has_ATOM_token, bool& is_water_chain, bool& is_100_percent)
{
       has_ATOM_token = false;
       is_water_chain = true;
       is_100_percent = false;

       bool is_short_chain = false;
       if ((end > start) && ((end - start) < 30)) is_short_chain = true;
       bool has_ATOMP = false;
       bool has_ATOMN = false;
       int count_P_O3 = 0;
       std::map<std::string, int> all_types;
       all_types.clear();
       for (unsigned int i = start; i <= end; ++i) {
            if (_residues[i]->token() == "ATOM") has_ATOM_token = true;
            if (_residues[i]->chain_type() == "ATOMP") has_ATOMP = true;
            if (_residues[i]->chain_type() == "ATOMN") has_ATOMN = true;
            if (is_short_chain && _residues[i]->find_atom("P") && (_residues[i]->find_atom("O3'") || _residues[i]->find_atom("O3*"))) count_P_O3++;
            add_type(all_types, _residues[i]->chain_type());
            if (!SeqCodeUtil::is_water_residue(_residues[i]->ResName())) is_water_chain = false;
       }
       if (all_types.size() == 1) is_100_percent = true;
       std::string chain_type = get_type(all_types);
       if ((chain_type != "ATOMP") && (chain_type != "ATOMN") && has_ATOMP && (end > start) && ((end - start) < 11)) chain_type = "ATOMP";
       if ((chain_type != "ATOMP") && (chain_type != "ATOMN")) {
            if (is_short_chain && has_ATOMN) chain_type = "ATOMN";
            // else if (count_P_O3 && ((count_P_O3 * 2) >= (int) (end - start))) chain_type = "ATOMN";
       }
       return chain_type;
}

bool Molecule::_check_linkage_between_residues(const std::string& chain_type, const unsigned int& start,
                                               const unsigned int& end)
{
       for (unsigned int i = start; i < end; ++i) {
            if (is_polymer_connect(chain_type, _residues[i], _residues[i + 1])) return true;
       }

       return false;
}

void Molecule::_insert_overlap_residue_index_map(const unsigned int& r_index_1, const unsigned int& r_index_2,
                                        std::map<unsigned int, std::set<unsigned int> >& overlap_index_map)
{
       std::map<unsigned int, std::set<unsigned int> >::iterator mpos = overlap_index_map.find(r_index_1);
       if (mpos != overlap_index_map.end()) mpos->second.insert(r_index_2);
       else {
            std::set<unsigned int> t_set;
            t_set.clear();
            t_set.insert(r_index_2);
            overlap_index_map.insert(std::make_pair(r_index_1, t_set));
       }
}

#if 0
void get_residue_indices(const std::vector<Molecule*>& mols, const int& format, int& mol_index,
                     std::vector<int>& ResIndex, const std::vector<std::vector<std::string> >& data,
                     std::vector<std::string>& errors)
{
       mol_index = 0;
       ResIndex.clear();
       errors.clear();

       RCSB::Residue* residue;
       bool found = false;
       for (unsigned int i = 0; i < data.size(); ++i) {
            for (unsigned int j = 0; j < mols.size(); ++j) {
                 // if (format == NDB_FILE_FORMAT_PDB)
                      residue = mols[j]->find_pdb_residue(data[i][0],
                                     data[i][1], data[i][2], data[i][3]);
/*
                 else residue = mols[j]->find_cif_residue(data[i][0],
                                     data[i][1], data[i][2]);
*/
                 if (residue) {
                      mol_index = j;
                      found = true;
                      break;
                 }
            }
            if (found) break;
       }

       if (!found) {
            for (unsigned int i = 0; i < data.size(); ++i) {
                 std::string cs = "Residue (" + data[i][0] + " " + data[i][1] + " "
                           + data[i][2] + data[i][3] + ") can not be found.";
                 errors.push_back(cs);
            }
            return;
       }

       for (unsigned int i = 0; i < data.size(); ++i) {
            // if (format == NDB_FILE_FORMAT_PDB)
                 residue = mols[mol_index]->find_pdb_residue(data[i][0],
                                     data[i][1], data[i][2], data[i][3]);
/*
            else residue = mols[mol_index]->find_cif_residue(data[i][0],
                                     data[i][1], data[i][2]);
*/
            if (residue)
                 ResIndex.push_back(residue->index());
            else {
                 std::string cs = "Residue (" + data[i][0] + " " + data[i][1] + " "
                           + data[i][2] + data[i][3] + ") can not be found.";
                 errors.push_back(cs);
            }
       }
}

void get_atom_indices(const std::vector<Molecule*>& mols, const int& format, int& mol_index,
                     std::vector<RCSB::Atom*>& atoms, const std::vector<std::vector<std::string> >& data,
                     std::vector<std::string>& errors)
{
       mol_index = 0;
       atoms.clear();
       errors.clear();

       RCSB::Residue* residue;
       bool found = false;
       for (unsigned int i = 0; i < data.size(); ++i) {
            for (unsigned int j = 0; j < mols.size(); ++j) {
                 // if (format == NDB_FILE_FORMAT_PDB)
                      residue = mols[j]->find_pdb_residue(data[i][0],
                                     data[i][1], data[i][2], data[i][3]);
/*
                 else residue = mols[j]->find_cif_residue(data[i][0],
                                     data[i][1], data[i][2]);
*/
                 if (residue) {
                      mol_index = j;
                      found = true;
                      break;
                 }
            }
            if (found) break;
       }

       if (!found) {
            for (unsigned int i = 0; i < data.size(); ++i) {
                 std::string cs = "Atom (" + data[i][4] + data[i][5] + " "
                           + data[i][0] + " " + data[i][1] + " " + data[i][2]
                           + data[i][3] + ") can not be found.";
                 errors.push_back(cs);
            }
            return;
       }

       for (unsigned int i = 0; i < data.size(); ++i) {
            // if (format == NDB_FILE_FORMAT_PDB)
                 residue = mols[mol_index]->find_pdb_residue(data[i][0],
                                     data[i][1], data[i][2], data[i][3]);
/*
            else residue = mols[mol_index]->find_cif_residue(data[i][0],
                                     data[i][1], data[i][2]);
*/
            RCSB::Atom* atom = NULL;
            if (residue) atom = residue->find_atom(data[i][4], data[i][5]);
            if (atom) atoms.push_back(atom);
            else {
                 std::string cs = "Atom (" + data[i][4] + data[i][5] + " "
                           + data[i][0] + " " + data[i][1] + " " + data[i][2]
                           + data[i][3] + ") can not be found.";
                 errors.push_back(cs);
            }
       }
}
#endif
