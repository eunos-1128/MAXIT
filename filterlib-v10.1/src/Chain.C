/*
FILE:     Chain.C
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


#include "AlignUtil.h"
#include "Chain.h"
#include "CompositeIndex.h"
#include "SeqCodeUtil.h"
#include "SideChainPattern.h"
#include "TypeDef.h"
#include "utillib.h"

using namespace RCSB;

static double score(const int& i, const int& j, void* data);

void Chain::clear()
{
       _ccDic = NULL;
       _logIo = NULL;
       _messageIo = NULL;
       _ChainID.clear();
       _PDB_ChainID.clear();
       _PUB_ChainID.clear();
       _PDB_ChainID_Flag.clear();
       _details.clear();
       _PolyType.clear();
       _EvidenceCode.clear();
       _chain_type.clear();
       _prev_entity_id.clear();
       _descriptor.clear();
       _entity_id.clear();
       _entity_key.clear();
       _has_sequence = false;
       _empty_chain = false;
       _missing_sequence = false;
       _ca_or_p_atom_only = 0;
       _isModified = false;
       _index = 0;
       _order = 0;
       _na_type = 0;
       _ResidueNumbers = 0;
       _residues.clear();
       _SeqRes.clear();
       _numIndex.clear();
       _pdbIndex.clear();
       // _cifIndex.clear();
       _linear_descriptors.clear();
       _branch_links.clear();
       _firstresidueIndex.clear();
       _first_pos = _firstresidueIndex.begin();
}

void Chain::Reset()
{
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if (*pos) delete *pos;
       }
       clear();
}

const bool Chain::empty() const            { return _SeqRes.empty(); }
const std::string& Chain::ChainID() const          { return _ChainID; }
const std::string& Chain::PDB_ChainID() const      { return _PDB_ChainID; }
const std::string& Chain::PUB_ChainID() const      { return _PUB_ChainID; }
const std::string& Chain::PDB_ChainID_Flag() const { return _PDB_ChainID_Flag; }
const std::string& Chain::details() const          { return _details; }
const std::string& Chain::PolyType() const         { return _PolyType; }
const std::string& Chain::EvidenceCode() const     { return _EvidenceCode; }
const std::string& Chain::chain_type() const       { return _chain_type; }
const int&  Chain::ResidueNumbers() const          { return _ResidueNumbers; }
const bool& Chain::has_sequence() const            { return _has_sequence; }
const bool& Chain::empty_chain() const             { return _empty_chain; }
const bool& Chain::missing_sequence() const        { return _missing_sequence; }
const int& Chain::ca_or_p_atom_only() const        { return _ca_or_p_atom_only; }
const bool& Chain::isModified() const              { return _isModified; }
const int&  Chain::index() const                   { return _index; }
const int&  Chain::order() const                   { return _order; }
const int&  Chain::na_type() const                 { return _na_type; }
const std::vector<std::map<std::string, std::string> >& Chain::get_linear_descriptors() const { return _linear_descriptors; }
const std::vector<std::map<std::string, std::string> >& Chain::get_branch_links() const { return _branch_links; }
const std::string& Chain::prev_entity_id() const   { return _prev_entity_id; }
const std::string& Chain::entity_id() const        { return _entity_id; }
const int Chain::int_prev_entity_id() const        { return atoi(_prev_entity_id.c_str()); }
const int Chain::int_entity_id() const             { return atoi(_entity_id.c_str()); }
const std::string& Chain::entity_key()
{
       if (_entity_key.empty()) _get_entity_key();
       return _entity_key;
}
const unsigned int Chain::SeqLen() const           { return _SeqRes.size(); }
const bool Chain::IsLastResidue() const { return (_first_pos == _firstresidueIndex.end()); }

bool Chain::need_carbohydrate_annotation()
{
       if (_linear_descriptors.empty()) return true;

       for (std::vector<std::map<std::string, std::string> >::const_iterator vpos = _linear_descriptors.begin(); vpos != _linear_descriptors.end(); ++vpos) {
            std::map<std::string, std::string>::const_iterator mpos = vpos->find("type");
            if ((mpos != vpos->end()) && (mpos->second == "WURCS")) return false;
       }

       return true;
}

void Chain::setCCDic(ConnectDic *ccdic)                { _ccDic = ccdic; }
void Chain::setLog(LogUtil *logPt)                     { _logIo = logPt; }
void Chain::setMessage(MessageUtil *message)           { _messageIo = message; }
void Chain::set_ChainID(const std::string& a)          { _ChainID = a; }
void Chain::set_PDB_ChainID(const std::string& a)      { _PDB_ChainID = a; }
void Chain::set_PUB_ChainID(const std::string& a)      { _PUB_ChainID = a; }
void Chain::set_PDB_ChainID_Flag(const std::string& a) { _PDB_ChainID_Flag = a; }
void Chain::set_details(const std::string& a)          { _details = a; }
void Chain::set_PolyType(const std::string& a)         { _PolyType = a; }
void Chain::set_EvidenceCode(const std::string& a)     { _EvidenceCode = a; }
void Chain::set_chain_type(const std::string& a)       { _chain_type = a; }
void Chain::set_missing_sequence()                     { _missing_sequence = true; } 
void Chain::Chain::set_index(const int& a)             { _index = a; __set_chain_index(); }
void Chain::set_order(const int& a)                    { _order = a; }
void Chain::set_na_type(const int& a)                  { _na_type = a; }
void Chain::set_isModified(const bool& a)              { _isModified = a; }
void Chain::set_entity_id(const int& a)                { _entity_id = String::IntToString(a); }
void Chain::set_entity_id(const std::string& a)        { _entity_id = a; }
void Chain::set_prev_entity_id(const std::string& a)   { _prev_entity_id = a; }
void Chain::set_descriptor(const std::string& a)       { _descriptor = a; }
void Chain::set_entity_key(const std::string& a)       { _entity_key = a; }
void Chain::set_linear_descriptors(const std::vector<std::map<std::string, std::string> >& linear_descriptors) { _linear_descriptors = linear_descriptors; }
void Chain::set_branch_links(const std::vector<std::map<std::string, std::string> >& branch_links) { _branch_links = branch_links; }
void Chain::set_has_sequence()                         { _has_sequence = true; }
void Chain::setReserve(const int& number)
{
       _residues.reserve(number);
       _SeqRes.reserve(number);
}

void Chain::_get_entity_key()
{
       if ((_chain_type == "ATOMS") && !_branch_links.empty()) {
            std::vector<std::string> items, reverse_items, links;
            items.clear();
            items.push_back("comp_id_1");
            items.push_back("entity_branch_list_num_1");
            items.push_back("atom_id_1");
            items.push_back("comp_id_2");
            items.push_back("entity_branch_list_num_2");
            items.push_back("atom_id_2");

            reverse_items.clear();
            reverse_items.push_back("comp_id_2");
            reverse_items.push_back("entity_branch_list_num_2");
            reverse_items.push_back("atom_id_2");
            reverse_items.push_back("comp_id_1");
            reverse_items.push_back("entity_branch_list_num_1");
            reverse_items.push_back("atom_id_1");

            bool found_all_items = true;
            for (std::vector<std::map<std::string, std::string> >::const_iterator vpos = _branch_links.begin(); vpos != _branch_links.end(); ++vpos) {
                 for (std::vector<std::string>::const_iterator ivpos = items.begin(); ivpos != items.end(); ++ivpos) {
                      std::map<std::string, std::string>::const_iterator mpos = vpos->find(*ivpos);
                      if (mpos == vpos->end()) {
                           found_all_items = false;
                           break;
                      }
                 }
                 if (!found_all_items) break;
            }

            // first key: smaller index number
            // second key: larger index number
            // value[0]: residue name with smaller index number
            // value[1]: residue number (index+1) with smaller index number
            // value[2]: linked atom name with smaller index number
            // value[3]: residue name with larger index number
            // value[4]: residue number (index+1) with larger index number
            // value[5]: linked atom name with larger index number
            std::multimap<int, std::multimap<int, std::vector<std::string> > > residue_link_mapping;
            residue_link_mapping.clear();
     
            if (found_all_items) {
                 std::multimap<int, std::vector<std::string> > tmp_link_mapping;
                 for (std::vector<std::map<std::string, std::string> >::const_iterator vpos = _branch_links.begin(); vpos != _branch_links.end(); ++vpos) {
                      links.clear();
                      std::map<std::string, std::string>::const_iterator mpos1 = vpos->find("entity_branch_list_num_1");
                      std::map<std::string, std::string>::const_iterator mpos2 = vpos->find("entity_branch_list_num_2");
                      int small_index = atoi(mpos1->second.c_str());
                      int large_index = atoi(mpos2->second.c_str());
                      if (small_index > large_index) {
                           small_index = atoi(mpos2->second.c_str());
                           large_index = atoi(mpos1->second.c_str());
                           for (std::vector<std::string>::const_iterator ivpos = reverse_items.begin(); ivpos != reverse_items.end(); ++ivpos) {
                                std::map<std::string, std::string>::const_iterator mpos = vpos->find(*ivpos);
                                links.push_back(mpos->second);
                           }
                      } else {
                           for (std::vector<std::string>::const_iterator ivpos = items.begin(); ivpos != items.end(); ++ivpos) {
                                std::map<std::string, std::string>::const_iterator mpos = vpos->find(*ivpos);
                                links.push_back(mpos->second);
                           }
                      }

                      std::multimap<int, std::multimap<int, std::vector<std::string> > >::iterator rlmpos = residue_link_mapping.find(small_index);
                      if (rlmpos != residue_link_mapping.end()) rlmpos->second.insert(std::make_pair(large_index, links));
                      else {
                           tmp_link_mapping.clear();
                           tmp_link_mapping.insert(std::make_pair(large_index, links));
                           residue_link_mapping.insert(std::make_pair(small_index, tmp_link_mapping));
                      }
                 }
            }

            if (!residue_link_mapping.empty()) {
                 // Create entity key with residue name, residue number and linked atom names
                 _entity_key = "";
                 for (std::multimap<int, std::multimap<int, std::vector<std::string> > >::const_iterator
                      mpos = residue_link_mapping.begin(); mpos != residue_link_mapping.end(); ++mpos) {
                      for (std::multimap<int, std::vector<std::string> >::const_iterator mpos1 = mpos->second.begin(); mpos1 != mpos->second.end(); ++mpos1) {
                           if (!_entity_key.empty()) _entity_key += "_";
                           _entity_key += CompositeIndex::getIndex(mpos1->second);
                      }
                 }
                 return;                 
            }
       }

       std::vector<std::string> seqs;
       get_seq(seqs);
       if (seqs.size() == 1)
            _entity_key = seqs[0];
       else _entity_key = CompositeIndex::getIndex(seqs);
}

void Chain::__set_chain_index()
{
       for (std::vector<std::vector<RCSB::Residue*> >::iterator pos = _residues.begin(); pos != _residues.end(); ++pos) {
            for (std::vector<RCSB::Residue*>::iterator pos1 = pos->begin(); pos1 != pos->end(); ++pos1) {
                 (*pos1)->set_chain_index(_index);
            }
       }
}

const unsigned int Chain::NumAtoms()
{
       int num = 0;
       RCSB::Residue* res = GetFirstResidue();
       while (res) {
            num += res->NumAtoms();
            res = GetNextResidue();
       }
       return num;
}

RCSB::Residue* Chain::GetFirstResidue()
{
       _first_pos = _firstresidueIndex.begin();
       if (_first_pos != _firstresidueIndex.end())
            return _residues[*_first_pos][0];
       else return NULL;
}

RCSB::Residue* Chain::GetNextResidue()
{
       ++_first_pos;
       if (_first_pos != _firstresidueIndex.end())
            return _residues[*_first_pos][0];
       else return NULL;
}

void Chain::GetFirstResidueList(std::vector<RCSB::Residue*>& residue_list)
{
       residue_list.clear();
       _first_pos = _firstresidueIndex.begin();
       if (_first_pos != _firstresidueIndex.end())
            residue_list = _residues[*_first_pos];
}

void Chain::GetNextResidueList(std::vector<RCSB::Residue*>& residue_list)
{
       residue_list.clear();
       ++_first_pos;
       if (_first_pos != _firstresidueIndex.end())
            residue_list = _residues[*_first_pos];
}

void Chain::GetLastResidueList(std::vector<RCSB::Residue*>& residue_list)
{
       residue_list.clear();
       if (!_firstresidueIndex.empty()) {
            int index = _firstresidueIndex.back();
            residue_list = _residues[index];
       }
}

void Chain::GetResidueList(const int& pos, std::vector<RCSB::Residue*>& residue_list)
{
       residue_list.clear();
       if (pos >= 0 && pos < (int) _SeqRes.size() && _SeqRes[pos]->ResIndex >= 0) {
                 residue_list = _residues[_SeqRes[pos]->ResIndex];
       }
}

void Chain::GetResidueListByIndex(const int& index, std::vector<RCSB::Residue*>& residue_list)
{
       residue_list.clear();
       if (index >= 0 && index < (int) _residues.size()) {
            residue_list = _residues[index];
       }
}

_FIELD* Chain::SeqRes(const int& pos)
{
       if (pos >= 0 && pos < (int) _SeqRes.size())
            return (_SeqRes[pos]);
       else return NULL;
}

void Chain::merge_chain(Chain* chn)
{
       _has_sequence = true;

       for (std::vector<_FIELD*>::iterator fpos = chn->_SeqRes.begin(); fpos != chn->_SeqRes.end(); ++fpos) {
            int ResIndex = (*fpos)->ResIndex;
            if (ResIndex >= 0) {
                 (*fpos)->ResIndex = _residues.size();
                 for (std::vector<RCSB::Residue*>::iterator rpos = chn->_residues[ResIndex].begin(); rpos != chn->_residues[ResIndex].end(); ++rpos) {
                      (*rpos)->set_chain_index(_index);
                      (*rpos)->set_position(_residues.size());
                 }
                 _residues.push_back(chn->_residues[ResIndex]);
            }
            _SeqRes.push_back(*fpos);
            *fpos = NULL;
       }
}

void Chain::merge_residue(RCSB::Residue *residue)
{
       _has_sequence = true;

       residue->set_chain_index(_index);
       residue->set_position(_residues.size());

       _FIELD *seqres = new _FIELD;
       seqres->ResIndex = _residues.size();

       RCSB::Atom *atom = residue->GetFirstAtom();
       seqres->InsCode  = atom->ins_code();
       seqres->Field[0] = atom->restype();
       seqres->Field[1] = atom->restype();
       seqres->Field[2] = atom->resnum();
       seqres->Field[3] = atom->pdb_resnam();
       seqres->Field[4] = atom->pdb_resnum();
       seqres->Field[5] = atom->pub_resnam();
       seqres->Field[6] = atom->pub_resnum();

       _SeqRes.push_back(seqres);

       std::vector<RCSB::Residue*> residue_list;
       residue_list.clear();
       residue_list.push_back(residue);
       _residues.push_back(residue_list);
}

void Chain::insert_a_residue(RCSB::Residue *residue, const bool& skip_alt_flag)
{
       residue->set_chain_index(_index);

       std::string cs = CompositeIndex::getIndex(residue->pdb_res_no(), residue->ins_code());
       std::map<std::string, int>::iterator pos = _numIndex.find(cs);
       if (!skip_alt_flag && pos != _numIndex.end() && residue->has_alt_loc()) {
            int index = _SeqRes[pos->second]->ResIndex;
            residue->set_position(index);
            cs = CompositeIndex::getIndex(residue->pdb_chnid(), residue->ResName(), residue->pdb_res_no(), residue->ins_code());
            _pdbIndex.insert(std::make_pair(cs, std::make_pair(pos->second, _residues[index].size())));
            _residues[index].push_back(residue);
#if 0
            /*
             * for heterogeneity residues, reorder residues based alternate code
             */
            std::multimap<std::string, RCSB::Residue*> tmp_map;
            tmp_map.clear();
            for (std::vector<RCSB::Residue*>::iterator vpos = _residues[index].begin(); vpos != _residues[index].end(); ++vpos) {
                 cs = CompositeIndex::getIndex((*vpos)->pdb_chnid(), (*vpos)->ResName(), (*vpos)->pdb_res_no(), (*vpos)->ins_code());
                 _pdbIndex.erase(cs);
/*
                 cs = CompositeIndex::getIndex((*vpos)->ResName(), (*vpos)->res_no(), "");
                 _cifIndex.erase(cs);
*/
                 tmp_map.insert(std::make_pair((*vpos)->alt_loc(), *vpos));
            }
            tmp_map.insert(std::make_pair(residue->alt_loc(), residue));

            _residues[index].clear();
            for (std::multimap<std::string, RCSB::Residue*>::iterator tpos = tmp_map.begin(); tpos != tmp_map.end(); ++tpos) {
                 cs = CompositeIndex::getIndex(tpos->second->pdb_chnid(), tpos->second->ResName(), tpos->second->pdb_res_no(), tpos->second->ins_code());
                 _pdbIndex.insert(std::make_pair(cs, std::make_pair(pos->second, _residues[index].size())));
/*
                 cs = CompositeIndex::getIndex(tpos->second->ResName(), tpos->second->res_no(), "");
                 _cifIndex.insert(std::make_pair(cs, std::make_pair(pos->second, _residues[index].size())));
*/
                 _residues[index].push_back(tpos->second);
            }
#endif
            return;
       }

       residue->set_position(_residues.size());
       std::vector<RCSB::Residue*> residue_list;
       residue_list.clear();
       residue_list.push_back(residue);
       _residues.push_back(residue_list);

       _FIELD *seqres = new _FIELD;
       seqres->ResIndex = _residues.size() - 1;

       RCSB::Atom *atom = residue->GetFirstAtom();
       seqres->InsCode  = atom->ins_code();
       seqres->Field[0] = atom->restype();
       seqres->Field[1] = atom->restype();
       seqres->Field[2] = atom->resnum();
       seqres->Field[3] = atom->pdb_resnam();
       seqres->Field[4] = atom->pdb_resnum();
       seqres->Field[5] = atom->pub_resnam();
       seqres->Field[6] = atom->pub_resnum();

       cs = CompositeIndex::getIndex(seqres->Field[4], seqres->InsCode);
       _numIndex.insert(std::make_pair(cs, _SeqRes.size()));

       cs = CompositeIndex::getIndex(atom->pdb_chnid(), seqres->Field[3], seqres->Field[4], seqres->InsCode);
       _pdbIndex.insert(std::make_pair(cs, std::make_pair(_SeqRes.size(), 0)));

       // cs = CompositeIndex::getIndex(seqres->Field[1], seqres->Field[2], "");
       // _cifIndex.insert(std::make_pair(cs, std::make_pair(_SeqRes.size(), 0)));

       _SeqRes.push_back(seqres);

       _firstresidueIndex.push_back(seqres->ResIndex);
       _ResidueNumbers++;
}

void Chain::insert_hetero_residues(const std::vector<RCSB::Residue*>& residue_list)
{
       for (unsigned int i = 0; i < residue_list.size(); ++i) {
            residue_list[i]->set_chain_index(_index);
            residue_list[i]->set_position(_residues.size());
            std::string cs = CompositeIndex::getIndex(residue_list[i]->pdb_res_no(), residue_list[i]->ins_code());
            _numIndex.insert(std::make_pair(cs, _SeqRes.size()));
            cs = CompositeIndex::getIndex(residue_list[i]->pdb_chnid(), residue_list[i]->ResName(), residue_list[i]->pdb_res_no(), residue_list[i]->ins_code());
            _pdbIndex.insert(std::make_pair(cs, std::make_pair(_SeqRes.size(), i)));
       }

       _residues.push_back(residue_list);

       RCSB::Atom *atom = residue_list[0]->GetFirstAtom();

       _FIELD *seqres = new _FIELD;
       seqres->ResIndex = _residues.size() - 1;
       seqres->InsCode  = atom->ins_code();
       seqres->Field[0] = atom->restype();
       seqres->Field[1] = atom->restype();
       seqres->Field[2] = atom->resnum();
       seqres->Field[3] = atom->pdb_resnam();
       seqres->Field[4] = atom->pdb_resnum();
       seqres->Field[5] = atom->pub_resnam();
       seqres->Field[6] = atom->pub_resnum();
       _SeqRes.push_back(seqres);

       _firstresidueIndex.push_back(seqres->ResIndex);
       _ResidueNumbers++;
}

void Chain::insert_OXT2N_residue(const std::map<unsigned int, RCSB::Residue*>& residue_mapping)
{
       if (residue_mapping.empty()) return; 

       std::vector<std::vector<RCSB::Residue*> > _tmp_residues;
       _tmp_residues.clear();
       _tmp_residues.reserve(_residues.size() + residue_mapping.size());

       std::vector<RCSB::Residue*> residue_list;

       for (unsigned int i = 0; i < _SeqRes.size(); ++i) {
            if (_SeqRes[i]->ResIndex >= 0) {
                 for (unsigned int j = 0; j < _residues[_SeqRes[i]->ResIndex].size(); ++j) {
                      _residues[_SeqRes[i]->ResIndex][j]->set_position(_tmp_residues.size());
                 }
                 _tmp_residues.push_back(_residues[_SeqRes[i]->ResIndex]);
                 _SeqRes[i]->ResIndex = _tmp_residues.size() - 1;
            } else {
                 std::map<unsigned int, RCSB::Residue*>::const_iterator mpos = residue_mapping.find(i);
                 if (mpos != residue_mapping.end()) {
                      RCSB::Residue *residue = mpos->second;
                      residue->set_chain_index(_index);
                      residue->set_position(_tmp_residues.size());
                      _SeqRes[i]->ResIndex = _tmp_residues.size();
                      residue_list.clear();
                      residue_list.push_back(residue);
                      _tmp_residues.push_back(residue_list);
                 }
            }
       }

       _residues = _tmp_residues;

       update_indices();
}

void Chain::remove_residue(const std::set<int>& residue_set)
{
       std::vector<RCSB::Residue*> residue_list;

       for (std::vector<_FIELD*>::iterator fpos = _SeqRes.begin(); fpos != _SeqRes.end(); ++fpos) {
            int ResIndex = (*fpos)->ResIndex;
            if (ResIndex < 0) continue;

            residue_list.clear();
            for (std::vector<RCSB::Residue*>::const_iterator rpos = _residues[ResIndex].begin(); rpos != _residues[ResIndex].end(); ++rpos) {
                 if (residue_set.find((*rpos)->index()) != residue_set.end()) continue;
                 residue_list.push_back(*rpos);
            }
            _residues[ResIndex] = residue_list;
            if (residue_list.empty()) {
                 delete *fpos;
                 *fpos = NULL;
            }
       }

       std::vector<_FIELD*> tmp_seqres;
       tmp_seqres.clear();
       for (std::vector<_FIELD*>::iterator fpos = _SeqRes.begin(); fpos != _SeqRes.end(); ++fpos) {
            if (*fpos == NULL) continue;
            tmp_seqres.push_back(*fpos);
       }
       _SeqRes = tmp_seqres;

       update_indices();
}

int Chain::find_field_index(const int& resnum, const std::string& ins_code)
{
       return find_field_index(String::IntToString(resnum), ins_code);
}

int Chain::find_field_index(const std::string& resnum, const std::string& ins_code)
{
       if (resnum.empty()) return (-1);
       std::string cs = CompositeIndex::getIndex(resnum, ins_code);
       std::map<std::string, int>::iterator pos = _numIndex.find(cs);
       if (pos != _numIndex.end())
            return pos->second;
       else return (-1);

}
/*
RCSB::Residue* Chain::find_cif_residue(const std::string& resname, const int resnum)
{
       return find_cif_residue(resname, String::IntToString(resnum));
}

RCSB::Residue* Chain::find_cif_residue(const std::string& resname, const std::string& resnum)
{
       std::string cs = CompositeIndex::getIndex(resname, resnum, "");
       std::map<std::string, std::pair<int, int> >::iterator pos = _cifIndex.find(cs);
       if (pos != _cifIndex.end()) {
            if (pos->second.first >= 0 && pos->second.second >= 0)
                 return _residues[pos->second.first][pos->second.second];
       }
       return NULL;
}
*/
RCSB::Residue* Chain::find_pdb_residue(const std::string& resname, const int resnum, const std::string& ins_code)
{
       return find_pdb_residue(resname, String::IntToString(resnum), ins_code); 
}

RCSB::Residue* Chain::find_pdb_residue(const std::string& resname, const std::string& resnum, const std::string& ins_code)
{
       std::string cs = CompositeIndex::getIndex(_PDB_ChainID, resname, resnum, ins_code);
       std::map<std::string, std::pair<int, int> >::iterator pos = _pdbIndex.find(cs);
       if (pos == _pdbIndex.end()) return NULL;

       if (pos->second.first < 0 || pos->second.second < 0) return NULL;

       if (_SeqRes[pos->second.first]->ResIndex < 0) return NULL;

       return _residues[_SeqRes[pos->second.first]->ResIndex][pos->second.second];
}

std::pair<int, RCSB::Residue*> Chain::find_pdb_residue_with_sort_index(const std::string& pdb_chnid, const std::string& resname, const std::string& resnum,
                                                                       const std::string& ins_code)
{
       RCSB::Residue* res = NULL;
       std::string cs = CompositeIndex::getIndex(pdb_chnid, resname, resnum, ins_code);
       std::map<std::string, std::pair<int, int> >::iterator pos = _pdbIndex.find(cs);
       if (pos == _pdbIndex.end()) return std::make_pair(-1, res);

       if (pos->second.first < 0 || pos->second.second < 0) return std::make_pair(-1, res);

       if (_SeqRes[pos->second.first]->ResIndex < 0) return std::make_pair(-1, res);

       // int sort_index = _SeqRes[pos->second.first]->ResIndex * 10 + pos->second.second;
       // return std::make_pair(sort_index, _residues[_SeqRes[pos->second.first]->ResIndex][pos->second.second]);
       return std::make_pair(_SeqRes[pos->second.first]->ResIndex, _residues[_SeqRes[pos->second.first]->ResIndex][pos->second.second]);
}

RCSB::Residue* Chain::find_prev_residue(const std::string& resname, const int resnum, const std::string& ins_code)
{
       return find_prev_residue(resname, String::IntToString(resnum), ins_code); 
}

RCSB::Residue* Chain::find_prev_residue(const std::string& resname, const std::string& resnum, const std::string& ins_code)
{
       std::string cs = CompositeIndex::getIndex(_PDB_ChainID, resname, resnum, ins_code);
       std::map<std::string, std::pair<int, int> >::iterator pos = _pdbIndex.find(cs);
       if (pos != _pdbIndex.end()) {
            if (pos->second.first >= 0 && pos->second.second >= 0) {
                 for (int i = pos->second.first - 1; i >= 0; --i) {
                      if (_SeqRes[i]->ResIndex < 0) continue;
                      return _residues[_SeqRes[i]->ResIndex][0];
                 }
            }
       }
       return NULL;
}

RCSB::Residue* Chain::find_next_residue(const std::string& resname, const std::string& resnum, const std::string& ins_code)
{
       std::string cs = CompositeIndex::getIndex(_PDB_ChainID, resname, resnum, ins_code);
       std::map<std::string, std::pair<int, int> >::iterator pos = _pdbIndex.find(cs);
       if (pos != _pdbIndex.end()) {
            if (pos->second.first >= 0 && pos->second.second >= 0) {
                 for (int i = pos->second.first + 1; i < (int) _SeqRes.size(); ++i) {
                      if (_SeqRes[i]->ResIndex < 0) continue;
                      return _residues[_SeqRes[i]->ResIndex][0];
                 }
            }
       }
       return NULL;
}

void Chain::insert_sequence(const std::vector<std::vector<std::string> >& mapping)
{
       _has_sequence = true;
       _empty_chain = true;

       for (unsigned int i = 0; i < mapping.size(); ++i) {
            _FIELD *seqres = new _FIELD;
            seqres->ResIndex = -1;

            seqres->InsCode.clear();
            seqres->Field[0] = mapping[i][1];
            seqres->Field[1] = mapping[i][1];
            seqres->Field[2] = String::IntToString(i + 1);
            seqres->Field[3].clear();
            seqres->Field[4] = mapping[i][3];
            seqres->Field[5].clear();
            seqres->Field[6].clear();

            _SeqRes.push_back(seqres);
       }
       update_indices();
}

void Chain::insert_sequence(const SEQ& seq)
{
       _has_sequence = true;
       _empty_chain = true;

       int i = 0;
       for (std::vector<std::string>::const_iterator pos = seq.res.begin(); pos != seq.res.end(); ++pos) {
            i++;
            _FIELD *seqres = new _FIELD;
            seqres->ResIndex = -1;

            seqres->InsCode.clear();
            seqres->Field[0] = *pos;
            seqres->Field[1] = *pos;
            seqres->Field[2] = String::IntToString(i);
            seqres->Field[3].clear();
            seqres->Field[4] = String::IntToString(i);
            seqres->Field[5].clear();
            seqres->Field[6].clear();

            _SeqRes.push_back(seqres);
       }
       update_indices();
}

bool Chain::seq_alignment(const std::vector<std::vector<std::string> >& mapping, const bool& rename_residue_flag)
{
       unsigned int count = 0;
       for (unsigned int i = 0; i < mapping.size(); ++i) {
            if (mapping[i][2].empty()) continue;
            if (count >= _SeqRes.size()) return false;

            if (mapping[i][2] != _SeqRes[count]->Field[3] ||
                mapping[i][3] != _SeqRes[count]->Field[4] ||
                mapping[i][4] != _SeqRes[count]->InsCode) return false;
            count++;
       }
       if (count != _SeqRes.size()) return false;

       _has_sequence = true;

       // Assign previous entity id if it does not exist in coordinate
       if (_prev_entity_id.empty()) {
            _prev_entity_id = mapping[0][5];
            _entity_id = _prev_entity_id;
       }

       if (mapping.size() > _SeqRes.size()) {
            std::vector<_FIELD*> tmp = _SeqRes;
            _SeqRes.clear();
            _SeqRes.reserve(mapping.size());
            count = 0;
            for (unsigned int i = 0; i < mapping.size(); ++i) {
                 if (mapping[i][2].empty()) {
                      _FIELD *seqres = new _FIELD;
                      seqres->ResIndex = -1;
                      for (int k = 0; k < NUM_FIELD; ++k) {
                           seqres->Field[k].clear();
                      }
                      seqres->Field[4] = mapping[i][3];
                      seqres->InsCode = mapping[i][4];
                      _SeqRes.push_back(seqres);
                 } else {
                      _SeqRes.push_back(tmp[count]);
                      count++;
                 }
            }
            tmp.clear();
       }

       bool need_renumber = false;
       std::set<std::string> index_set;
       index_set.clear();
       for (unsigned int i = 0; i < mapping.size(); ++i) {
            _SeqRes[i]->Field[0] = mapping[i][1];
            if (!_SeqRes[i]->Field[0].empty() && !_SeqRes[i]->Field[1].empty() && _SeqRes[i]->Field[0] != _SeqRes[i]->Field[1] && _SeqRes[i]->ResIndex >= 0) {
                 RCSB::Residue *res = _residues[_SeqRes[i]->ResIndex][0];
                 if (res->OrigResName() == _SeqRes[i]->Field[0] && res->ResName() == _SeqRes[i]->Field[1]) {
                      _SeqRes[i]->Field[0] = res->ResName();
                 }
            }
            // if (mapping[i][2].empty()) need_renumber = true;
            std::string index = _SeqRes[i]->Field[4] + "_" + _SeqRes[i]->InsCode;
            if (index_set.find(index) != index_set.end())
                 need_renumber = true;
            else index_set.insert(index);
       }
       if (need_renumber) check_unique_number();

       if (rename_residue_flag) {
            for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
                 if (((*pos)->ResIndex < 0) || ((*pos)->Field[0] == (*pos)->Field[1])) continue;
                 if (((*pos)->Field[1] == "GLY") || ((*pos)->Field[1] == "ALA")) {
                      for (unsigned int i = 0; i < _residues[(*pos)->ResIndex].size(); ++i) {
                           if (_residues[(*pos)->ResIndex][i]->ResName() == (*pos)->Field[1]) {
                                if (_residues[(*pos)->ResIndex][i]->rename("ala-gly-like", (*pos)->Field[0])) {
                                     (*pos)->Field[1] = (*pos)->Field[0];
                                     (*pos)->Field[3] = (*pos)->Field[0];
                                }
                           }
                      }
                 }
            }
       }

       update_indices();
       // if (_chain_type == "ATOMN" && !original_flag) rename_dna_residues(_na_type);

       return true;
}

bool Chain::Mapping_SeqTool_Alignment(std::vector<std::vector<std::string> >& mapping)
{
       // maping[][0]: auth_mon_id
       // maping[][1]: auth_mon_num
       // maping[][2]: pdb_mon_id
       // maping[][3]: pdb_mon_num
       // maping[][4]: insertion code

       // find all mapped residues with coordinates
       std::vector<std::vector<std::string> > tmp_mapping;
       tmp_mapping.clear();
       tmp_mapping.reserve(mapping.size());
       for (std::vector<std::vector<std::string> >::const_iterator pos = mapping.begin(); pos != mapping.end(); ++pos) {
            if ((*pos)[2].empty() || (*pos)[3].empty()) continue;
            tmp_mapping.push_back(*pos);
       }

       // find all residues with coordinates
       std::vector<_FIELD*> tmp_seqres;
       tmp_seqres.clear();
       tmp_seqres.reserve(_SeqRes.size());
 
       for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->ResIndex < 0) continue;
            tmp_seqres.push_back(*pos);
       }

       // check if have same number of residues with coordinates
       if (tmp_mapping.size() != tmp_seqres.size()) return false;

       int count = 0;
       for (unsigned int i = 0; i < tmp_mapping.size(); ++i) {
            // residue numbering must be same
            if (tmp_mapping[i][3] != tmp_seqres[i]->Field[4]) return false;
            if ((tmp_mapping[i][2] == tmp_seqres[i]->Field[3]) || (tmp_seqres[i]->Field[3] == "ALA") ||
                (tmp_seqres[i]->Field[3] == "GLY") || (tmp_seqres[i]->Field[3] == "UNK")) count++;
       }

       // must have more than 90% residues with same residue name
       double ratio = (double) count / (double) tmp_mapping.size();
       if (ratio < 0.9) return false;

       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->ResIndex < 0) delete *pos;
       }

       _has_sequence = true;

       _SeqRes.clear();
       _SeqRes.reserve(mapping.size());
       count = 0;
       bool need_renumber = false;
       for (std::vector<std::vector<std::string> >::const_iterator pos = mapping.begin(); pos != mapping.end(); ++pos) {
            if ((*pos)[2].empty() || (*pos)[3].empty()) {
                 need_renumber = true;
                 _FIELD *seqres = new _FIELD;
                 seqres->ResIndex = -1;
                 seqres->InsCode.clear();
                 for (int k = 0; k < NUM_FIELD; ++k) {
                      seqres->Field[k].clear();
                 }
                 seqres->Field[0] = (*pos)[0];
                 seqres->Field[1] = (*pos)[0];
                 _SeqRes.push_back(seqres);
            } else {
                 if ((*pos)[2] != tmp_seqres[count]->Field[3]) {
                      bool changed = false;
                      std::string sidechain_pattern = "";
                      if ((*pos)[2] == "UNK")
                           changed = true;
                      else if (tmp_seqres[count]->Field[3] == "GLY" || tmp_seqres[count]->Field[3] == "ALA") {
                           changed = true;
                           sidechain_pattern = "ala-gly-like";
                      } else if (/* (tmp_seqres[count]->Field[3] == "ASP" && (*pos)[2] == "ASN") || */
                               (tmp_seqres[count]->Field[3] == "ASN" && (*pos)[2] == "ASP"))
                           changed = true;
                      else if (/* (tmp_seqres[count]->Field[3] == "GLU" && (*pos)[2] == "GLN") || */
                               (tmp_seqres[count]->Field[3] == "GLN" && (*pos)[2] == "GLU"))
                           changed = true;
                      else {
                           sidechain_pattern = _residues[tmp_seqres[count]->ResIndex][0]->get_side_chain_pattern();
                           if (SideChainPattern::is_a_valid_rename(sidechain_pattern, tmp_seqres[count]->Field[3], (*pos)[2])) changed = true;
                      }
                      if (changed) {
                           int i = tmp_seqres[count]->ResIndex;
                           for (unsigned int j = 0; j < _residues[i].size(); ++j) {
                                if (_residues[i][j]->ResName() == tmp_seqres[count]->Field[3]) {
                                     if (_residues[i][j]->rename(sidechain_pattern, (*pos)[2])) {
                                          tmp_seqres[count]->Field[1] = (*pos)[2];
                                          tmp_seqres[count]->Field[3] = (*pos)[2];
                                     }
                                }
                           }
                      }
                 }
                 tmp_seqres[count]->Field[0] = (*pos)[0];
                 _SeqRes.push_back(tmp_seqres[count]);
                 count++;
            }
       }
       tmp_seqres.clear();

       if (need_renumber) check_unique_number();

       for (unsigned int i = 0; i < mapping.size(); ++i) {
            mapping[i][3] = _SeqRes[i]->Field[4];
            mapping[i][4] = _SeqRes[i]->InsCode;
       }

       update_indices();

       _get_entity_key();

       // if (_chain_type == ATOMN_TOKEN && !original_flag) rename_dna_residues(_na_type);

       return true;
}

void Chain::seq_alignment(const SEQ& seq, const std::map<std::string, std::pair<std::string, std::string> >& missing_residue_numbering_mapping)
{
       _has_sequence = true;

       // Assign previous entity id if it does not exist in coordinate
       if (_prev_entity_id.empty()) {
            _prev_entity_id = seq.entity_id;
            _entity_id = _prev_entity_id;
       }

       std::vector<std::vector<int> > sa, sb, ss;
       AlignUtil::initAssignment(sa, seq.res.size());
       AlignUtil::initAssignment(sb, _SeqRes.size());

       std::vector<int> relative_numbering;
       relative_numbering.clear();

       _DataContainer cdata;
       cdata._seqa = seq.res;
       cdata._seqb.clear();
       cdata._seqb.reserve(_SeqRes.size());

       int start_number = 0;
       std::string prev_number = "";
       std::string prev_inscode = "";
       for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            cdata._seqb.push_back((*pos)->Field[1]);
            if ((*pos)->ResIndex < 0) start_number++;
            else {
                 if (prev_number.empty()) start_number++;
                 else {
                      if (prev_number == (*pos)->Field[4]) {
                           if (!(*pos)->InsCode.empty()) {
                                char prev_ins = '@';
                                if (!prev_inscode.empty()) prev_ins = prev_inscode[0];
                                char curr_ins = (*pos)->InsCode[0];
                                start_number += abs(int(curr_ins) - int(prev_ins)); 
                           } else start_number++;
                      } else {
                           start_number += abs(atoi((*pos)->Field[4].c_str()) - atoi(prev_number.c_str()));
                      }
                 }
                 prev_number = (*pos)->Field[4];
                 prev_inscode = (*pos)->InsCode;
            }
            relative_numbering.push_back(start_number);
       }

       std::vector<bool> linkage;
       linkage.clear();

       for (unsigned int i = 1; i < _SeqRes.size(); ++i) {
            if (_SeqRes[i - 1]->ResIndex < 0 || _SeqRes[i]->ResIndex < 0)
                 linkage.push_back(false);
            else if (_IsConnect(_residues[_SeqRes[i - 1]->ResIndex][0], _residues[_SeqRes[i]->ResIndex][0]))
                 linkage.push_back(true);
            else linkage.push_back(false);
       }
       if (_SeqRes[_SeqRes.size()-1]->ResIndex < 0)
            linkage.push_back(false);
       else linkage.push_back(true);

       AlignUtil::Alignment(ss, sa, sb, 2, 0, score, &cdata, linkage);

       bool found_mismatch = false;
       bool found_missing = false;
       for (std::vector<std::vector<int> >::iterator pos = ss.begin(); pos != ss.end(); ++pos) {
            if ((*pos)[0] >= 0 && (*pos)[1] < 0) found_missing = true;
            if ((*pos)[0] < 0 || (*pos)[1] < 0) continue;
            if (seq.res[(*pos)[0]] != _SeqRes[(*pos)[1]]->Field[1] && _SeqRes[(*pos)[1]]->Field[1] != "ALA" && _SeqRes[(*pos)[1]]->Field[1] != "GLY") {
                 found_mismatch = true;
                 break;
            }
       }
       if (found_mismatch) {
            ss.clear();
            linkage.clear();
            AlignUtil::Alignment(ss, sa, sb, 2, 0, score, &cdata, linkage);
       }

       if (found_missing) check_alignment(cdata, relative_numbering, ss);

       std::vector<_FIELD*> tmp = _SeqRes;
       _SeqRes.clear();
       _SeqRes.reserve(ss.size());

       bool found_gap = false;
       for (std::vector<std::vector<int> >::iterator pos = ss.begin(); pos != ss.end(); ++pos) {
            if ((*pos)[1] >= 0) {
                 _SeqRes.push_back(tmp[(*pos)[1]]);
            } else {
                 _FIELD *seqres = new _FIELD;
                 seqres->ResIndex = -1;
                 seqres->InsCode.clear();
                 for (int k = 0; k < NUM_FIELD; ++k) {
                      seqres->Field[k].clear();
                 }
                 _SeqRes.push_back(seqres);
                 found_gap = true;
            }
            if ((*pos)[0] >= 0)
                 _SeqRes[_SeqRes.size()-1]->Field[0] = seq.res[(*pos)[0]];
            else _SeqRes[_SeqRes.size()-1]->Field[0].clear();
       }
       tmp.clear();

       std::string cs;
       if (found_gap) {
            if (!missing_residue_numbering_mapping.empty()) {
                 found_gap = false;
                 int i = 0;
                 for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
                      if ((*pos)->ResIndex < 0 && (*pos)->Field[4].empty()) {
                           cs = String::IntToString(i) + "_" + (*pos)->Field[0];
                           std::map<std::string, std::pair<std::string, std::string> >::const_iterator mpos = missing_residue_numbering_mapping.find(cs);
                           if (mpos != missing_residue_numbering_mapping.end()) {
                                (*pos)->Field[4] = mpos->second.first;
                                (*pos)->InsCode = mpos->second.second;
                           } else found_gap = true;
                      }
                      i++;
                 }
            }
            if (found_gap) check_unique_number();
       }

       bool need_renumber = false;
       std::set<std::string> UniqueIndex;
       UniqueIndex.clear();

       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if (!(*pos)->Field[0].empty() && !(*pos)->Field[1].empty() && (*pos)->Field[0] != (*pos)->Field[1] && (*pos)->ResIndex >= 0) {
                 RCSB::Residue *res = _residues[(*pos)->ResIndex][0];
                 if (res->OrigResName() == (*pos)->Field[0] && res->ResName() == (*pos)->Field[1]) {
                      (*pos)->Field[0] = res->ResName();
                 }
            }
            if ((*pos)->ResIndex < 0 && (*pos)->Field[4].empty()) {
                 need_renumber = true;
                 // break;
            }
            std::string cs = (*pos)->Field[4] + (*pos)->InsCode;
            if (UniqueIndex.find(cs) != UniqueIndex.end()) {
                 need_renumber = true;
                 // break;
            }
            UniqueIndex.insert(cs);
       }
       if (need_renumber) check_unique_number();

       update_indices();
       // if (_chain_type == "ATOMN" && !original_flag) rename_dna_residues(_na_type);
}

bool Chain::check_seq_alignment(const std::vector<std::string>& authSeqs, std::vector<std::string>& seq_coord_labels)
{
       unsigned int seq_count = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->Field[0].empty()) continue;
            seq_count++;
       }
       if (seq_count != authSeqs.size()) return false;

       seq_count = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->Field[0].empty()) continue;
            if ((*pos)->ResIndex >= 0) seq_coord_labels[seq_count] = "1";
            seq_count++;
       }
       
       return true;
}

std::string Chain::checking_missing_residue_in_coordinates()
{
       unsigned int count_missing_residue = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->ResIndex >= 0) continue;
            count_missing_residue++;
       }
       if ((count_missing_residue * 10) < (_SeqRes.size() * 3)) return "";

       double percent = (((double) count_missing_residue) * 100.0) / ((double) _SeqRes.size());
       return FloatToString(percent, 0, 1) + "% residues of chain '" + _PDB_ChainID + "' are missing in coordinates.";
}

unsigned int Chain::count_missing_residue_in_coordinates()
{
       unsigned int count_missing_residue = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->ResIndex >= 0) continue;
            count_missing_residue++;
       }
       return count_missing_residue;
}

void Chain::update_missing_residues(const std::vector<std::vector<std::string> >& missinglist)
{
       int index = -1;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->ResIndex >= 0) continue;
            index++;
            if (index >= (int) missinglist.size()) continue;
            if (missinglist[index][0] == (*pos)->Field[0]) {
                 (*pos)->Field[4] = missinglist[index][1];
                 (*pos)->InsCode = missinglist[index][2];
            }
       }
}

void Chain::check_alignment(const _DataContainer& cdata, const std::vector<int>& relative_numbering, std::vector<std::vector<int> >& ss)
{
       // adjust a single residue at the N-/C- terminals of loop range
       int begin = -1;
       int end = -1;
       for (unsigned int i = 1; i < ss.size() - 1; ++i) {
            if (ss[i][0] >= 0 && ss[i][1] >= 0 && ss[i+1][0] >= 0 && ss[i+1][1] < 0) {
                 begin = i; end = -1;
            }
            if (ss[i-1][0] >= 0 && ss[i-1][1] < 0 && ss[i][0] >= 0 && ss[i][1] >= 0) {
                 if (begin >= 0 && begin < (int) i)
                      end = i;
                 else {
                      begin = -1;
                      end = -1;
                 }
            }
            if (begin < 0 || end < 0) continue;
      
            int first = ss[begin][1];
            int second = ss[end][1]; 
            RCSB::Residue *rfirst = _residues[_SeqRes[first]->ResIndex][0];
            RCSB::Residue *rsecond = _residues[_SeqRes[second]->ResIndex][0];
            RCSB::Residue *prev = NULL;
            RCSB::Residue *next = NULL;
            if ((first - 1) >= 0 && _SeqRes[first-1]->ResIndex >= 0)
                 prev = _residues[_SeqRes[first-1]->ResIndex][0];
            if ((second + 1) < (int) _SeqRes.size() && _SeqRes[second + 1]->ResIndex >= 0)
                 next = _residues[_SeqRes[second + 1]->ResIndex][0];

            if (prev && next && _IsConnect(rfirst, rsecond)) {
                 bool moved = false;
                 int dif1 = abs(atoi(rfirst->pdb_res_no().c_str()) - atoi(prev->pdb_res_no().c_str()));
                 int dif2 = abs(atoi(rsecond->pdb_res_no().c_str()) - atoi(rfirst->pdb_res_no().c_str()));
                 int dif3 = abs(atoi(next->pdb_res_no().c_str()) - atoi(rsecond->pdb_res_no().c_str()));
                 if (cdata._seqa[ss[begin+1][0]] == cdata._seqa[ss[end][0]] && dif2 < dif3) {
                      // Move to N-terminal of the loop range
                      ss[begin+1][1] = ss[end][1];
                      ss[end][1] = -1;
                      moved = true;
                 }
                 if (!moved && dif1 > dif2 && cdata._seqa[ss[begin][0]] == cdata._seqa[ss[end-1][0]]) {
                      // Mobe to C-terminal of the loop range
                      ss[end-1][1] = ss[begin][1];
                      ss[begin][1] = -1;
                 }
            }
            i = end;
            begin = -1;
            end = -1;
       }

       std::vector<int> tmp_vec;
       tmp_vec.clear();
       for (int i = 0; i < 4; ++i) tmp_vec.push_back(0);

       // vector[0]: 1 for block, 0 for loop
       // vector[1]: length of the range
       // vector[2]: begin index of the range
       // vector[3]: end index of the range
       std::vector<std::vector<int> > block_loop_ranges;
       block_loop_ranges.clear();
       int loop_begin = -1;
       int loop_end = -1;
       int block_begin = -1;
       int block_end = -1;
       for (unsigned int i = 0; i < ss.size(); ++i) {
            if (ss[i][0] < 0) continue;
            if (ss[i][1] < 0) {
                 if (block_begin >= 0) {
                      tmp_vec[0] = 1;
                      tmp_vec[1] = block_end - block_begin + 1;
                      tmp_vec[2] = block_begin;
                      tmp_vec[3] = block_end;
                      block_loop_ranges.push_back(tmp_vec);
                 }
                 block_begin = -1;
                 block_end = -1;
                 if (loop_begin < 0) loop_begin = i;
                 loop_end = i;
            } else {
                 if (loop_begin >= 0) {
                      tmp_vec[0] = 0;
                      tmp_vec[1] = loop_end - loop_begin + 1;
                      tmp_vec[2] = loop_begin;
                      tmp_vec[3] = loop_end;
                      block_loop_ranges.push_back(tmp_vec);
                 }
                 loop_begin = -1;
                 loop_end = -1;
                 if (block_begin < 0) block_begin = i;
                 block_end = i;
            }
       }
       if (block_begin >= 0) {
            tmp_vec[0] = 1;
            tmp_vec[1] = block_end - block_begin + 1;
            tmp_vec[2] = block_begin;
            tmp_vec[3] = block_end;
            block_loop_ranges.push_back(tmp_vec);
       }
       if (loop_begin >= 0) {
            tmp_vec[0] = 1;
            tmp_vec[1] = loop_end - loop_begin + 1;
            tmp_vec[2] = loop_begin;
            tmp_vec[3] = loop_end;
            block_loop_ranges.push_back(tmp_vec);
       }

       std::vector<int> match_list, prev_block, next_block;
       for (unsigned int i = 0; i < block_loop_ranges.size(); ++i) {
            if (!block_loop_ranges[i][0]) continue;

            prev_block.clear();
            if ((i >= 2) && block_loop_ranges[i-2][0]) prev_block = block_loop_ranges[i-2];
            next_block.clear();
            if ((i < (block_loop_ranges.size() - 2)) && block_loop_ranges[i+2][0]) next_block = block_loop_ranges[i+2];
            if (prev_block.empty() && next_block.empty()) continue;

            // check with begore loop range
            if ((i >= 1) && !block_loop_ranges[i-1][0] && (block_loop_ranges[i-1][1] > (block_loop_ranges[i][1] + 1)) &&
                find_match_list(cdata, relative_numbering, ss, block_loop_ranges[i], block_loop_ranges[i-1], prev_block, next_block, match_list)) {
                 for (unsigned int j = 0; j < match_list.size(); ++j) {
                      ss[match_list[j]][1] = ss[block_loop_ranges[i][2] + j][1];
                      ss[block_loop_ranges[i][2] + j][1] = -1;
                 }
                 block_loop_ranges[i-1][3] = match_list[0] - 1;
                 block_loop_ranges[i-1][1] = block_loop_ranges[i-1][3] - block_loop_ranges[i-1][2] + 1;
                 block_loop_ranges[i][2] = match_list[0];
                 block_loop_ranges[i][3] = match_list[match_list.size() - 1];
                 if ((i < (block_loop_ranges.size() - 1)) && !block_loop_ranges[i+1][0]) {
                      block_loop_ranges[i+1][2] = match_list[match_list.size() - 1] + 1;
                      block_loop_ranges[i+1][1] = block_loop_ranges[i+1][3] - block_loop_ranges[i+1][2] + 1;
                 }
                 continue;
            }

            // check with after loop range
            if ((i < (block_loop_ranges.size() - 1)) && !block_loop_ranges[i+1][0] && (block_loop_ranges[i+1][1] > (block_loop_ranges[i][1] + 1)) &&
                find_match_list(cdata, relative_numbering, ss, block_loop_ranges[i], block_loop_ranges[i+1], prev_block, next_block, match_list)) {
                 for (unsigned int j = 0; j < match_list.size(); ++j) {
                      ss[match_list[j]][1] = ss[block_loop_ranges[i][2] + j][1];
                      ss[block_loop_ranges[i][2] + j][1] = -1;
                 }
                 if ((i >= 1) && !block_loop_ranges[i-1][0]) {
                      block_loop_ranges[i-1][3] = match_list[0] - 1;
                      block_loop_ranges[i-1][1] = block_loop_ranges[i-1][3] - block_loop_ranges[i-1][2] + 1;
                 }
                 block_loop_ranges[i][2] = match_list[0];
                 block_loop_ranges[i][3] = match_list[match_list.size() - 1];
                 block_loop_ranges[i+1][2] = match_list[match_list.size() - 1] + 1;
                 block_loop_ranges[i+1][1] = block_loop_ranges[i+1][3] - block_loop_ranges[i+1][2] + 1;
            }
       }
}

bool Chain::find_match_list(const _DataContainer& cdata, const std::vector<int>& relative_numbering, const std::vector<std::vector<int> >& ss,
                            const std::vector<int>& block_range, const std::vector<int>& loop_range, const std::vector<int>& prev_block,
                            const std::vector<int>& next_block, std::vector<int>& match_list)
{
       match_list.clear();

       int diff_value = 0;
       if (!prev_block.empty()) {
            diff_value += abs(abs(relative_numbering[ss[block_range[2]][1]] - relative_numbering[ss[prev_block[3]][1]])
                        - abs(ss[block_range[2]][0] - ss[prev_block[3]][0]));
       }
       if (!next_block.empty()) {
            diff_value += abs(abs(relative_numbering[ss[next_block[2]][1]] - relative_numbering[ss[block_range[3]][1]])
                        - abs(ss[next_block[2]][0] - ss[block_range[3]][0]));
       }

       std::vector<int> tmp_match_list;
       for (int i = loop_range[2] + 1; i <= loop_range[3] - block_range[1]; ++i) {
            tmp_match_list.clear();
            for (int j = 0; j < block_range[1]; ++j) {
                 if (cdata._seqa[ss[block_range[2] + j][0]] == cdata._seqa[ss[i + j][0]]) { 
                      tmp_match_list.push_back(i + j); 
                 }
            }
            if ((int) tmp_match_list.size() == block_range[1]) {
                 int value = 0;
                 if (!prev_block.empty()) {
                      value += abs(relative_numbering[ss[block_range[2]][1]] - relative_numbering[ss[prev_block[3]][1]])
                             - abs(ss[tmp_match_list[0]][0] - ss[prev_block[3]][0]);
                 }
                 if (!next_block.empty()) {
                      value += abs(relative_numbering[ss[next_block[2]][1]] - relative_numbering[ss[block_range[3]][1]])
                             - abs(ss[next_block[2]][0] - ss[tmp_match_list[tmp_match_list.size()-1]][0]);
                 }
                 if (value < diff_value) {
                      diff_value = value;
                      match_list = tmp_match_list;
                 }
            }
       }
       return !match_list.empty();
}

void Chain::check_unique_number()
{
       bool is_increased = false;
       bool is_decreased = false;
       int last_number = -999999;
       for (int i = 0; i < (int) _SeqRes.size(); ++i) {
            _SeqRes[i]->Field[2] = String::IntToString(i + 1);
            if (_SeqRes[i]->ResIndex < 0) continue;
       
            if (last_number != -999999) {
                 if (atoi(_SeqRes[i]->Field[4].c_str()) < last_number)
                      is_decreased = true;
                 else is_increased = true;
            }
            last_number = atoi(_SeqRes[i]->Field[4].c_str());
       }
       bool has_reverse_number = false;
       if (is_decreased && !is_increased) has_reverse_number = true;

       char alphabet[27] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
       int last = -1;
       last_number = -999999;
       int last_insert_code = 0;
       for (int i = 0; i < (int) _SeqRes.size(); ++i) {
            if (_SeqRes[i]->ResIndex < 0) continue;
            if (i > last + 1) {
                 int seq_num = atoi(_SeqRes[i]->Field[4].c_str());
                 int current_insert_code = _SeqRes[i]->InsCode[0];
                 if ((seq_num - last_number) >= (i - last)) {
                      for (int j = i - 1; j > last; j--) {
                           if (has_reverse_number)
                                seq_num++;
                           else seq_num--;
                           _SeqRes[j]->Field[4]= String::IntToString(seq_num);
                      }
                 } else if ((i - last - 1) <= 26 * (seq_num - last_number)) {
                      seq_num = 0;
                      if (last_insert_code > 0) {
                           for (int j = 0; j < 26; j++) {
                                if (alphabet[j] == last_insert_code) {
                                     seq_num = j + 1;
                                     if (seq_num >= 26) seq_num = 0;
                                     break;
                                }
                           }
                      }
                      for (int j = last + 1; j < i; j++) {
                           _SeqRes[j]->Field[4] = String::IntToString(last_number);
                           _SeqRes[j]->InsCode = alphabet[seq_num];
                           seq_num++;
                           if (seq_num >= 26) {
                                seq_num = 0;
                                last_number++;
                           }
                      }
                 } else if (seq_num == last_number && current_insert_code > 0 && last_insert_code > 0 && (current_insert_code-last_insert_code) >= (i-last)) {
                      for (int j = last + 1; j < i; j++) {
                           _SeqRes[j]->Field[4] = String::IntToString(last_number);
                           last_insert_code++;
                           _SeqRes[j]->InsCode = (char) last_insert_code;
                      }
                 } else if (seq_num == last_number && current_insert_code > 0 && last_insert_code == 0 && (current_insert_code - 64) >= (i - last)) {
                      for (int j = i - 1; j > last; j--) {
                           _SeqRes[j]->Field[4] = String::IntToString(last_number);
                           current_insert_code--;
                           _SeqRes[j]->InsCode = (char) current_insert_code;
                      }
                 }
            }
            last_number = atoi(_SeqRes[i]->Field[4]. c_str());
            last_insert_code = _SeqRes[i]->InsCode[0];
            last = i;
       }
       if (last >= 0) {
            int seq_num = atoi(_SeqRes[last]->Field[4].c_str());
            for (int i = last + 1; i < (int) _SeqRes.size(); ++i) {
                 if (_SeqRes[i]->ResIndex >= 0) break;
                 if (has_reverse_number)
                      seq_num--;
                 else seq_num++;
                 _SeqRes[i]->Field[4] = String::IntToString(seq_num);
            }
       }
}

void Chain::update_indices()
{
       _pdbIndex.clear();
       _numIndex.clear();
       // _cifIndex.clear();
       _firstresidueIndex.clear();
       _ResidueNumbers = 0;
       int i = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            std::string cs = CompositeIndex::getIndex((*pos)->Field[4], (*pos)->InsCode);
            _numIndex.insert(std::make_pair(cs, i));
     
            if ((*pos)->ResIndex < 0) {
                 cs = CompositeIndex::getIndex(_PDB_ChainID, (*pos)->Field[0], (*pos)->Field[4], (*pos)->InsCode);
                 _pdbIndex.insert(std::make_pair(cs, std::make_pair(i, -1)));
          
                 // cs = CompositeIndex::getIndex((*pos)->Field[1], (*pos)->Field[2], "");
                 // _cifIndex.insert(std::make_pair(cs, std::make_pair(i, -1)));
            } else {
                 _ResidueNumbers++;
                 _firstresidueIndex.push_back((*pos)->ResIndex);
                 int j = 0;
                 for (std::vector<RCSB::Residue*>::iterator vpos = _residues[(*pos)->ResIndex].begin(); vpos != _residues[(*pos)->ResIndex].end(); ++vpos) {
                      cs = CompositeIndex::getIndex((*vpos)->pdb_chnid(), (*vpos)->ResName(), (*pos)->Field[4], (*pos)->InsCode);
                      _pdbIndex.insert(std::make_pair(cs, std::make_pair(i, j)));
                      //  add extra index if the pdb number or insertion code are different
                      if (((*pos)->Field[4] != (*vpos)->pdb_res_no()) || ((*pos)->InsCode != (*vpos)->ins_code())) {
                           cs = CompositeIndex::getIndex((*vpos)->pdb_chnid(), (*vpos)->ResName(), (*vpos)->pdb_res_no(), (*vpos)->ins_code());
                           _pdbIndex.insert(std::make_pair(cs, std::make_pair(i, j)));
                      }
/*
                      cs = CompositeIndex::getIndex((*vpos)->ResName(), (*pos)->Field[2], "");
                      _cifIndex.insert(std::make_pair(cs, std::make_pair(i, j)));
*/
                      j++;

                 }
            }
            i++;
       }
}

void Chain::update_number()
{
       update_indices();
       int i = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            i++;
            if ((*pos)->ResIndex < 0)
                 (*pos)->Field[4] = String::IntToString(i);
            else (*pos)->Field[4] = _residues[(*pos)->ResIndex][0]->pdb_res_no();
            (*pos)->InsCode.clear();
       }
}

bool Chain::has_sequence_mismatch()
{
       if (_chain_type != "ATOMN" && _chain_type != "ATOMP") return false;

       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->Field[1].empty()) continue;
            if (SeqCodeUtil::isSameResName((*pos)->Field[0], (*pos)->Field[1])) continue;
            if (String::IsEqual((*pos)->Field[1], "ALA", Char::eCASE_INSENSITIVE) ||
                String::IsEqual((*pos)->Field[1], "GLY", Char::eCASE_INSENSITIVE)) continue;

            return true;
       }
       return false;
}

void Chain::extract_alignment(std::vector<std::vector<std::string> >& alignments)
{
       // [0]: PDB Chain ID
       // [1]: Entity ID
       // [2]: three letter code sequence
       // [3]: one letter code sequence
       // [4]: sequence number
       // [5]: three letter code residue
       // [6]: one letter code residue
       // [7]: residue number + insertion code
       // [8]: mismatch
       alignments.clear();

       std::vector<std::string> data;
       int count = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            data.clear();
            data.push_back(_PDB_ChainID);
            data.push_back(_entity_id);
            data.push_back((*pos)->Field[0]);
            if (!(*pos)->Field[0].empty()) {
                 data.push_back(SeqCodeUtil::GetOneLetterCodeWithMSE((*pos)->Field[0]));
                 count++;
                 data.push_back(String::IntToString(count));
            } else {
                 data.push_back("");
                 data.push_back("");
            }
            data.push_back((*pos)->Field[1]);
            if (!(*pos)->Field[1].empty())
                 data.push_back(SeqCodeUtil::GetOneLetterCodeWithMSE((*pos)->Field[1]));
            else data.push_back("");
            data.push_back((*pos)->Field[4] + (*pos)->InsCode);

            std::string mismatch = "";
            if (!(*pos)->Field[1].empty() && !SeqCodeUtil::isSameResName((*pos)->Field[0], (*pos)->Field[1])) {
                 if ((String::IsEqual((*pos)->Field[1], "ALA", Char::eCASE_INSENSITIVE) ||
                      String::IsEqual((*pos)->Field[1], "GLY", Char::eCASE_INSENSITIVE)) && !(*pos)->Field[0].empty())
                      mismatch = "ala_gly_mismatch";
                 else mismatch = "conflict";
            }
            data.push_back(mismatch);

            alignments.push_back(data);
       }
}

void Chain::check_polymer_sequence(const int& Mol_ID)
{
       if (_chain_type != "ATOMN" && _chain_type != "ATOMP") return;

       int met_count = 0;
       int mse_count = 0;
       bool has_mismatch = false;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->Field[0] == "MET") met_count++;
            if ((*pos)->Field[0] == "MSE") mse_count++;
            if ((*pos)->Field[1].empty()) continue;
            if (SeqCodeUtil::isSameResName((*pos)->Field[0], (*pos)->Field[1])) continue;

            has_mismatch = true;

            std::string error = "SEQRES/COORD MISMATCH:";
            if (Mol_ID >= 0) error += " Model=" + String::IntToString(Mol_ID);
            error += " ChainID=" + _PDB_ChainID + " ResNum=" + (*pos)->Field[4] + " (" + (*pos)->Field[0] + " <---> " + (*pos)->Field[1] + ")";
            _messageIo->insertMessage("mismatch_print", "error", error);
       }
       if (met_count > 0 && mse_count > 0) {
            std::string message = "Both MET and MSE residues are present in the sequence. ";
            message += "Verify with the authors.";
            _messageIo->insertMessage("met_mse_mixture_print", "warning", message);
       }

       if (!has_mismatch) return;

       std::string alignment = print_alignment(40);
       _messageIo->insertMessage("mismatch_print", "error", alignment);
}

std::string Chain::print_alignment(const int& num_per_line)
{
       std::string name, message;

       message.clear();

       message += "PDB Chain_ID: " + _PDB_ChainID + "\n\n";

       int l = _SeqRes.size() / num_per_line;
       int x = _SeqRes.size() % num_per_line;
       int m = l;
       if (x == 0) m = l - 1;
       int num = 0, n, y;
       int first = -1, first_num = 0, last = -1, last_num = 0;
       for (int i = 0; i <= m; ++i) {
            n = num_per_line;
            if (i == l) n = x;

            message += "           ";
            first = -1;
            last = -1;
            for (y = 0; y < n; y++) {
                 if (!_SeqRes[i * num_per_line + y]->Field[0].empty()) {
                      num++;
                      if (first == -1) {
                           first = y;
                           first_num = num;
                      }
                      last = y;
                      last_num = num;
                 }
            }
            for (y = 0; y < n; y++) {
                 if (y == first)
                      message += FloatToString((double) first_num, 4, 0, true, false);
                 else if (y == last)
                      message += FloatToString((double) last_num, 4, 0, false, false);
                 else if (y == last-1)
                      message += "   ";
                 else message += "    ";
            }

            message += "\n   SEQRES: ";
            for (y = 0; y < n; y++) {
                 name = _SeqRes[i * num_per_line + y]->Field[0];
                 if (name.empty()) name = "?";
                 message += FormattedString(name, 3) + " ";
            }

            message += "\n   COORDS: ";
            bool error = false;
            std::string tmp_line = "           ";
            for (y = 0; y < n; y++) {
                 name = _SeqRes[i * num_per_line + y]->Field[1];
                 if (name.empty()) name = "?";

                 if (!SeqCodeUtil::isSameResName(_SeqRes[i*num_per_line+y]->Field[0],
                      _SeqRes[i * num_per_line + y]->Field[1]) &&
                     !_SeqRes[i * num_per_line + y]->Field[1].empty()) {
                      error = true;
                      tmp_line += FormattedString("^^^", 3) + " "; 
                      String::LowerCase(name);
                 } else tmp_line += "    ";

                 message += FormattedString(name, 3) + " ";
            }
            message += "\n           ";
            first = -1; last = -1;
            for (y = 0; y < n; y++) {
                 if (!_SeqRes[i * num_per_line + y]->Field[1].empty()) {
                      if (first == -1) first = y;
                      last = y;
                 }
            }
            for (y = 0; y < n; y++) {
                 if (y == first)
                      message += FormattedString(_SeqRes[i * num_per_line + y]->Field[4], 4, true, false);
                 else if (y == last)
                      message += FormattedString(_SeqRes[i * num_per_line + y]->Field[4], 4);
                 else if (y == last - 1)
                      message += "   ";
                 else message += "    ";
            }
            message += "\n";
            if (error) { message += tmp_line; message += "\n"; }
            message += "\n";
       }

       return message;
}

void Chain::update_asymId()
{
       int i = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            (*pos)->Field[2] = String::IntToString(i + 1);
            if ((*pos)->ResIndex >= 0) {
                 for (std::vector<RCSB::Residue*>::iterator vpos = _residues[(*pos)->ResIndex].begin(); vpos != _residues[(*pos)->ResIndex].end(); ++vpos) {
                      std::string resnum = (*pos)->Field[2];
                      if (_chain_type != "ATOMN" && _chain_type != "ATOMP") resnum.clear();
                      /* if (_chain_type == "ATOMS") (*vpos)->update_nomenclature(_chain_type, resnum);
                      else */ (*vpos)->update_nomenclature(_chain_type, _ChainID, resnum);
                 }
            }
            i++;
       }
}

void Chain::update_residues_nomenclature(const bool& cif_only_flag, const bool& serial_flag, int &serial_no, const bool& insertion_flag)
{
       bool cif_only = cif_only_flag;
       if (serial_flag) cif_only = false;

       // _cifIndex.clear();
       if (!cif_only) {
            _pdbIndex.clear();
            _numIndex.clear();
       }
       std::string cs;
       int i = 0;
       for (std::vector<_FIELD*>::iterator
            pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            (*pos)->Field[2] = String::IntToString(i + 1);
            if (serial_flag) (*pos)->Field[4] = String::IntToString(serial_no);
            serial_no++;
            if (insertion_flag) (*pos)->InsCode.clear();

            if (!cif_only) {
                 cs = (*pos)->Field[4] + (*pos)->InsCode;
                 _numIndex.insert(std::make_pair(cs, i));
            }
     
            if ((*pos)->ResIndex < 0) {
                 if (!cif_only) {
                      cs = CompositeIndex::getIndex(_PDB_ChainID, (*pos)->Field[3], (*pos)->Field[4], (*pos)->InsCode);
                      _pdbIndex.insert(std::make_pair(cs, std::make_pair(i, -1)));
                 }
/*
                 cs = CompositeIndex::getIndex((*pos)->Field[1], (*pos)->Field[2], "");
                 _cifIndex.insert(std::make_pair(cs, std::make_pair(i, -1)));
*/
            } else {
                 int j = 0;
                 for (std::vector<RCSB::Residue*>::iterator vpos = _residues[(*pos)->ResIndex].begin();
                      vpos != _residues[(*pos)->ResIndex].end(); ++vpos) {
                      std::string resnum = (*pos)->Field[2];
                      if (_chain_type != "ATOMN" && _chain_type != "ATOMP") resnum.clear();
                      if (cif_only)
                           (*vpos)->update_nomenclature(_chain_type, _ChainID, resnum);
                      else (*vpos)->update_nomenclature(_chain_type, _ChainID, resnum, _PDB_ChainID, (*pos)->Field[4], (*pos)->InsCode);
                      if (!cif_only) {
                           cs = CompositeIndex::getIndex((*vpos)->pdb_chnid(), (*vpos)->ResName(), (*pos)->Field[4], (*pos)->InsCode);
                           _pdbIndex.insert(std::make_pair(cs, std::make_pair(i, j)));
                      }
/*
                      cs = CompositeIndex::getIndex((*vpos)->ResName(), (*pos)->Field[2], "");
                      _cifIndex.insert(std::make_pair(cs, std::make_pair(i, j)));
*/
                      j++;
                 }
            }
            i++;
       }
}

void Chain::update_residues_nomenclature(const std::vector<std::string>& numbering, const bool& insertion_flag)
{
       if (numbering.size() != _SeqRes.size()) return;

       // _cifIndex.clear();
       _pdbIndex.clear();
       _numIndex.clear();

       std::string cs;
       int i = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            (*pos)->Field[2] = String::IntToString(i + 1);
            (*pos)->Field[4] = numbering[i];
            if (insertion_flag) (*pos)->InsCode.clear();

            cs = (*pos)->Field[4] + (*pos)->InsCode;
            _numIndex.insert(std::make_pair(cs, i));
     
            if ((*pos)->ResIndex < 0) {
                 cs = CompositeIndex::getIndex(_PDB_ChainID, (*pos)->Field[3], (*pos)->Field[4], (*pos)->InsCode);
                 _pdbIndex.insert(std::make_pair(cs, std::make_pair(i, -1)));
/*
                 cs = CompositeIndex::getIndex((*pos)->Field[1], (*pos)->Field[2], "");
                 _cifIndex.insert(std::make_pair(cs, std::make_pair(i, -1)));
*/
            } else {
                 int j = 0;
                 for (std::vector<RCSB::Residue*>::iterator vpos = _residues[(*pos)->ResIndex].begin();
                      vpos != _residues[(*pos)->ResIndex].end(); ++vpos) {
                      std::string resnum = (*pos)->Field[2];
                      if (_chain_type != "ATOMN" && _chain_type != "ATOMP") resnum.clear();
                      (*vpos)->update_nomenclature(_chain_type, _ChainID, resnum, _PDB_ChainID, (*pos)->Field[4], (*pos)->InsCode);
                      cs = CompositeIndex::getIndex((*vpos)->pdb_chnid(), (*vpos)->ResName(), (*pos)->Field[4], (*pos)->InsCode);
                      _pdbIndex.insert(std::make_pair(cs, std::make_pair(i, j)));
/*
                      cs = CompositeIndex::getIndex((*vpos)->ResName(), (*pos)->Field[2], "");
                      _cifIndex.insert(std::make_pair(cs, std::make_pair(i, j)));
*/
                      j++;
                 }
            }
            i++;
       }
}

void Chain::set_PDB_ChainID_to_residues(const std::string& chain_id)
{
       _PDB_ChainID = chain_id;
       for (std::vector<std::vector<RCSB::Residue*> >::iterator pos = _residues.begin(); pos != _residues.end(); ++pos) {
            for (std::vector<RCSB::Residue*>::iterator vpos = pos->begin(); vpos != pos->end(); ++vpos) {
                 (*vpos)->set_pdb_chnid(chain_id);
                 RCSB::Atom *atom = (*vpos)->GetFirstAtom();
                 while (atom) {
                      atom->set_pdb_chnid(chain_id);
                      atom = (*vpos)->GetNextAtom();
                 }
            }
       }
}

void Chain::renumbering_polymer_chain(const bool& check_uniqueness_flag)
{
       int start = -1;
       bool is_pdb_number_unique = true;
       std::set<std::string> unique_set;
       unique_set.clear();
       for (unsigned int i = 0; i < _SeqRes.size(); ++i) {
            std::string idx = _SeqRes[i]->Field[4] + "_" + _SeqRes[i]->InsCode;
            if (_SeqRes[i]->Field[4].empty() || (unique_set.find(idx) != unique_set.end())) {
                 is_pdb_number_unique = false;
            }
            unique_set.insert(idx);
            if (!_SeqRes[i]->Field[3].empty() && (start < 0)) start = i;
       }
       if (is_pdb_number_unique && check_uniqueness_flag) return;

       int serial_no = atoi(_SeqRes[0]->Field[4].c_str());
       if ((serial_no < 0) && (start > 0)) serial_no = atoi(_SeqRes[start]->Field[4].c_str()) - start;
       update_residues_nomenclature(false, true, serial_no, true);
}

void Chain::renumbering_polymer_chain(Chain *template_chain)
{
       if (_SeqRes.size() != template_chain->_SeqRes.size())
            renumbering_polymer_chain();
       else {
            for (unsigned int i = 0; i < _SeqRes.size(); ++i) {
                 _SeqRes[i]->Field[2] = template_chain->_SeqRes[i]->Field[2];
                 _SeqRes[i]->Field[4] = template_chain->_SeqRes[i]->Field[4];
                 _SeqRes[i]->InsCode  = template_chain->_SeqRes[i]->InsCode;
                 if (_SeqRes[i]->ResIndex < 0) continue;
     
                 for (std::vector<RCSB::Residue*>::iterator vpos = _residues[_SeqRes[i]->ResIndex].begin(); vpos != _residues[_SeqRes[i]->ResIndex].end(); ++vpos) {
                      (*vpos)->update_nomenclature(_chain_type, _ChainID, _SeqRes[i]->Field[2], _PDB_ChainID, _SeqRes[i]->Field[4], _SeqRes[i]->InsCode);
                 }
            }
            update_indices();
       }
}

void Chain::renumbering_polymer_chain_starting_with_one()
{
       for (unsigned int i = 0; i < _SeqRes.size(); ++i) {
            _SeqRes[i]->Field[4] = String::IntToString(i + 1);
            _SeqRes[i]->InsCode.clear(); // remove insertion code
            if (_SeqRes[i]->ResIndex < 0) continue;

            for (std::vector<RCSB::Residue*>::iterator vpos = _residues[_SeqRes[i]->ResIndex].begin(); vpos != _residues[_SeqRes[i]->ResIndex].end(); ++vpos) {
                 (*vpos)->update_nomenclature(_chain_type, _ChainID, _SeqRes[i]->Field[2], _PDB_ChainID, _SeqRes[i]->Field[4], _SeqRes[i]->InsCode);
            }
       }

}

void Chain::correction_residues_name()
{
       int i = 0;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->Field[0] != (*pos)->Field[1] || (*pos)->Field[0] != (*pos)->Field[3]) {
                 std::string cs = CompositeIndex::getIndex(_PDB_ChainID, (*pos)->Field[3], (*pos)->Field[4], (*pos)->InsCode);
                 _pdbIndex.erase(cs);
                 (*pos)->Field[3] = (*pos)->Field[0];
                 cs = CompositeIndex::getIndex(_PDB_ChainID, (*pos)->Field[3], (*pos)->Field[4], (*pos)->InsCode);
                 if ((*pos)->ResIndex >= 0)
                      _pdbIndex.insert(std::make_pair(cs, std::make_pair(i, 0)));
                 else _pdbIndex.insert(std::make_pair(cs, std::make_pair(i, -1)));
/*
                 cs = CompositeIndex::getIndex((*pos)->Field[1], (*pos)->Field[2], "");
                 _cifIndex.erase(cs);
                 (*pos)->Field[1] = (*pos)->Field[0];
                 cs = CompositeIndex::getIndex((*pos)->Field[1], (*pos)->Field[2], "");
                 if ((*pos)->ResIndex >= 0)
                      _cifIndex.insert(std::make_pair(cs, std::make_pair(i, 0)));
                 else _cifIndex.insert(std::make_pair(cs, std::make_pair(i, -1)));
*/
                 if ((*pos)->ResIndex < 0) continue;

                 _residues[(*pos)->ResIndex][0]->correction_name((*pos)->Field[0]);
            }
            i++;
       }
}

float Chain::weight()
{
       if (_chain_type == "HETAS" && (_SeqRes[0]->Field[0] == "HOH" || _SeqRes[0]->Field[0] == "DOD")) return 18.015;
       else {
            float _weight = 0.0, Water = 18.0152, PO2 = 62.97256;;
            for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
                 _weight += _ccDic->get_molecule_weight((*pos)->Field[0]);
            }
            if (_SeqRes.size() > 1) {
                 _weight -= (double) (_SeqRes.size() - 1) * Water;
                 if (_chain_type == "ATOMN") _weight -= PO2;
            }
            return _weight;
       }
}

void Chain::get_seq(std::vector<std::string>& seqs)
{
       seqs.clear();
       if (_chain_type == "HETAS" && (_SeqRes[0]->Field[0] == "HOH" || _SeqRes[0]->Field[0] == "DOD")) seqs.push_back(_SeqRes[0]->Field[0]);
       else {
            seqs.reserve(_SeqRes.size());
            for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
                 if (!(*pos)->Field[0].empty()) seqs.push_back((*pos)->Field[0]);
            }
       }
}

void Chain::get_seq(std::list<std::string>& seqs)
{
       seqs.clear();

       if (!_has_sequence && (_chain_type != "ATOMN") && (_chain_type != "ATOMP") && (_chain_type != "ATOMS")) return;

       for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if (!(*pos)->Field[0].empty()) seqs.push_back((*pos)->Field[0]);
       }
}

void Chain::get_seq(std::list<std::vector<std::string> >& seqs)
{
       seqs.clear();

       if (!_has_sequence && (_chain_type != "ATOMN") && (_chain_type != "ATOMP") && (_chain_type != "ATOMS")) return;

       std::vector<std::string> data;
       for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            data.clear();
            data.push_back((*pos)->Field[0]);
            if ((*pos)->ResIndex >= 0) {
                 std::string index = _residues[(*pos)->ResIndex][0]->pdb_chnid() + "_" + _residues[(*pos)->ResIndex][0]->ResName() + "_"
                                   + _residues[(*pos)->ResIndex][0]->pdb_res_no() + "_" + _residues[(*pos)->ResIndex][0]->ins_code();
                 data.push_back(index);
            } else data.push_back("");
            seqs.push_back(data);
       }
}

void Chain::get_seqres(const int& format, const int& NumField, std::list<std::vector<std::string> >& seqres)
{
       if (!_has_sequence && _chain_type != "ATOMN" && _chain_type != "ATOMP") return;
/*
       if (_chain_type == "ATOMP" && _SeqRes.size() < 3 || _chain_type == "ATOMN" && _SeqRes.size() < 2) return;
*/

       std::vector<std::string> data;
       for (int i = 0; i < NumField; ++i) data.push_back("");
       data[0] = "SEQRES";
       if (format == NDB_FILE_FORMAT_PDB)
            data[2] = _PDB_ChainID;
       else data[2] = _ChainID;
       data[3] = String::IntToString(_SeqRes.size());
       int seqField = 4;
       for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if (seqField == NumField) {
                 if (!data[4].empty()) seqres.push_back(data);
                 for (int i = 4; i < NumField; ++i) data[i].clear();
                 seqField = 4;
            }
            data[seqField] = (*pos)->Field[0];
            seqField++;
       }
       if (!data[4].empty()) seqres.push_back(data);
}

std::string Chain::get_poly_type()
{
       if (_SeqRes.empty() || (_chain_type != "ATOMP" && _chain_type != "ATOMN")) return "";

       if (_chain_type == "ATOMP") {
            int dtype = 0;
            int ltype = 0;
            for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
                 if (SeqCodeUtil::is_standard_aa_residue((*pos)->Field[0])) ltype++;
                 else if (SeqCodeUtil::is_standard_Daa_residue((*pos)->Field[0])) dtype++;
            }
            if (dtype > ltype)
                 return "polypeptide(D)";
            else return "polypeptide(L)";
       } else {
            int dtype = 0;
            int rtype = 0;
            for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
                 if (SeqCodeUtil::is_rna_residue((*pos)->Field[0])) rtype++;
                 else if (SeqCodeUtil::is_dna_residue((*pos)->Field[0])) dtype++;
            }
            if (dtype && rtype)
                 return "polydeoxyribonucleotide/polyribonucleotide hybrid";
            else if (rtype && !dtype)
                 return "polyribonucleotide";
            else return "polydeoxyribonucleotide";
       }
}

void Chain::set_chain_type_to_residues()
{
       for (std::vector<std::vector<RCSB::Residue*> >::iterator pos = _residues.begin(); pos != _residues.end(); ++pos) {
            for (std::vector<RCSB::Residue*>::iterator ppos = pos->begin(); ppos != pos->end(); ++ppos) {
                 (*ppos)->set_chain_type(_chain_type);
            }
       }
}

void Chain::check_ca_or_p_atom_only()
{
       bool is_ca_only_chain = true;
       bool is_p_only_chain = true;
       for (std::vector<std::vector<RCSB::Residue*> >::const_iterator pos = _residues.begin(); pos != _residues.end(); ++pos) {
            if ((*pos)[0]->ca_or_p_atom_only() != CAAtom_ONLY) is_ca_only_chain = false;
            if ((*pos)[0]->ca_or_p_atom_only() != PAtom_ONLY) is_p_only_chain = false;
       }
       if (is_ca_only_chain) _ca_or_p_atom_only = CAAtom_ONLY;
       if (is_p_only_chain) _ca_or_p_atom_only = PAtom_ONLY;
}

void Chain::update_na_type()
{
       if (_PolyType == "polydeoxyribonucleotide/polyribonucleotide hybrid")
            _na_type = ATOMN_TYPE_DNA_RNA;
       else if (_PolyType == "polyribonucleotide")
            _na_type = ATOMN_TYPE_RNA_ONLY;
       else if (_PolyType == "polydeoxyribonucleotide")
            _na_type = ATOMN_TYPE_DNA_ONLY;

       if (_na_type) return;

       bool exist_DNA = false;
       bool exist_RNA = false; 
       for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            int type = _ccDic->find_monomer_type((*pos)->Field[0]);
            if (type == MONOMER_TYPE_DNA)
                 exist_DNA = true;
            else if (type == MONOMER_TYPE_RNA)
                 exist_RNA = true;
            if ((*pos)->ResIndex < 0) continue;

            for (std::vector<RCSB::Residue*>::iterator rpos = _residues[(*pos)->ResIndex].begin(); rpos != _residues[(*pos)->ResIndex].end(); ++rpos) {
                 (*rpos)->set_type(type);
            }
       }

       if (exist_DNA && exist_RNA)
            _na_type = ATOMN_TYPE_DNA_RNA;
       else if (exist_RNA)
            _na_type = ATOMN_TYPE_RNA_ONLY;
       else if (exist_DNA)
            _na_type = ATOMN_TYPE_DNA_ONLY;
}

void Chain::get_descriptor(std::string &descriptor)
{
       if (!_descriptor.empty()) {
            descriptor = _descriptor;
            return;
       }

       std::string sugar, prev_sugar;

       descriptor.clear();
       if (_chain_type == "ATOMP") {
            if (_SeqRes.size() > 1 && _SeqRes.size() <= 24) {
                 for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
                      if ((*pos)->Field[0].empty() || (*pos)->Field[0] == "?") continue;
                      if (!descriptor.empty()) descriptor += "-";
                      descriptor += (*pos)->Field[0];
                 }
            }
       } else if (_chain_type == "ATOMS" && _SeqRes.size() > 1) {
            // descriptor = String::IntToString(_SeqRes.size()) + "-MER";
       } else if (_chain_type == "ATOMN") {
            if (_SeqRes.size() > 24) {
                 descriptor = String::IntToString(_SeqRes.size()) + "-MER";
                 return;
            }

            descriptor = "5\'-";

            std::string prev_sugar = "XXX";
            bool first = true;
            for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
                 std::string sugar = "";
                 std::string phosphate = "";
                 if ((*pos)->ResIndex < 0) {
                      if (!first) phosphate = "P";
                      // int type = _ccDic->find_monomer_type((*pos)->Field[0]);
                      // if (type == MONOMER_TYPE_RNA) sugar = "R";
                 } else {
                      RCSB::Residue *residue = _residues[(*pos)->ResIndex][0];
                      // if (residue->type() == MONOMER_TYPE_RNA) sugar = "R";
                      RCSB::Atom *atom = residue->find_atom("P");
                      if (atom) phosphate = "P";
                 }
                 first = false;

                 int type = _ccDic->find_monomer_type((*pos)->Field[0]);
                 if (type == MONOMER_TYPE_RNA)
                      sugar = "R";
                 else if (type == MONOMER_TYPE_DNA)
                      sugar = "D";
                 else if (prev_sugar != "XXX")
                      sugar = prev_sugar;
                 else sugar = "D";

                 if (sugar != prev_sugar) {
                      if (prev_sugar != "XXX") descriptor += ")-";
                      descriptor += sugar + "(";
                 }
                 prev_sugar = sugar;

                 descriptor += phosphate + "*";
                 std::string code = SeqCodeUtil::GetOneLetterCode((*pos)->Field[0]);
                 if (code.size() == 1)
                      descriptor += code;
                 else descriptor += "(" + code + ")";
            }
            descriptor += ")-3\'";
       } else if (_SeqRes[0]->Field[0] == "UNL")
            descriptor = "UNKNOWN LIGAND";
       else {
            try {
                 const ConnectFormat& drug = _ccDic->find_drug(_SeqRes[0]->Field[0]);
                 descriptor = drug.chemical_name();
            } catch (const std::exception& exc) {} 
       }
}

void Chain::check_inscode(const int& Mol_ID)
{
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->InsCode.empty()) continue;
            if ((*pos)->InsCode[0] < 'A' || ((*pos)->InsCode[0] > 'Z' && (*pos)->InsCode[0] < 'a') || (*pos)->InsCode[0] > 'z') {
                 std::string cs = "";
                 if (Mol_ID >= 0) cs += "In model " + String::IntToString(Mol_ID) + ",";
                 cs += "Residue " + (*pos)->Field[3] + " " + _PDB_ChainID + " " + (*pos)->Field[4] + " has wrong insertion code '" + (*pos)->InsCode + "'\n";
                 _messageIo->insertMessage("ins_code", "label", "true");
                 _messageIo->insertMessage("ins_code_print", "warning", cs);
            }
       }
}
/*
void Chain::rename_dna_residues(const int& type)
{
       if (type != ATOMN_TYPE_DNA_ONLY && type != ATOMN_TYPE_DNA_RNA) return;

       bool changed = false;
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->ResIndex < 0) {
                 if (type != ATOMN_TYPE_DNA_ONLY) continue;
                 if ((*pos)->Field[0].size() > 1) continue;
                 std::string name = SeqCodeUtil::GetThreeLetterCode((*pos)->Field[0][0], "POLYDEOXYRIBONUCLEOTIDE");
                 if (name != (*pos)->Field[0]) {
                      (*pos)->Field[0] = name;
                      changed = true;
                 }
                 continue;
            }

            for (unsigned int j = 0; j < _residues[(*pos)->ResIndex].size(); ++j) {
                 RCSB::Residue *residue = _residues[(*pos)->ResIndex][j];
                 if (residue->ResName().size() > 1) continue;
                 int residue_type = ATOMN_TYPE_DNA_ONLY;
                 if (type == ATOMN_TYPE_DNA_RNA && (residue->find_atom("O2'") || residue->find_atom("O2*"))) residue_type = ATOMN_TYPE_RNA_ONLY;
                 if (residue_type != ATOMN_TYPE_DNA_ONLY) continue;
                 std::string name = SeqCodeUtil::GetThreeLetterCode(residue->ResName()[0], "POLYDEOXYRIBONUCLEOTIDE");
                 if (name != residue->ResName() || name == "DT") {
                      if (j == 0) {
                           (*pos)->Field[0] = name;
                           (*pos)->Field[1] = name;
                           (*pos)->Field[3] = name;
                           changed = true;
                      }
                      residue->correction_name(name);
                 }
            }
       }
       if (changed) update_indices();
}
*/
void Chain::get_prev_entity_id()
{
       std::map<std::string, int> all_types;
       all_types.clear();

       for (std::vector<std::vector<RCSB::Residue*> >::const_iterator pos = _residues.begin(); pos != _residues.end(); ++pos) {
            for (std::vector<RCSB::Residue*>::const_iterator ptr = pos->begin(); ptr != pos->end(); ++ptr) {
                 RCSB::Atom *atom = (*ptr)->GetFirstAtom();
                 const std::string& entity_id = atom->getValue("label_entity_id");
                 if (!entity_id.empty()) add_type(all_types, entity_id);
            }
       }

       if (all_types.empty()) return;

       std::string entity_id = get_type(all_types, "");
       if (!entity_id.empty()) {
            _prev_entity_id = entity_id;
            _entity_id = _prev_entity_id;
       }
}

void Chain::get_split_chains(const std::vector<std::vector<std::string> >& splits, const std::vector<std::vector<std::string> >& chain_info,
                             int& chain_index, std::vector<Chain*>& split_chains)
{
       split_chains.clear();

       std::map<std::string, std::pair<std::string, int> > number_mapping;
       number_mapping.clear();

       int count = 0;
       int idx = -1;
       for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            idx++;
            if ((*pos)->Field[0].empty()) continue;
            count++; 
            number_mapping.insert(std::make_pair(String::IntToString(count), std::make_pair((*pos)->Field[0], idx)));
       }

       std::vector<int> end_points;
       end_points.clear();

       bool found_error = false;
       for (std::vector<std::vector<std::string> >::const_iterator pos = splits.begin(); pos != splits.end(); ++pos) {
            std::map<std::string, std::pair<std::string, int> >::const_iterator mpos = number_mapping.find((*pos)[1]);
            if (mpos == number_mapping.end()) {
                 std::string error = "Residue " + (*pos)[0] + " " + (*pos)[1] + " does not exist chain " + _PDB_ChainID + ".\n";
                 _logIo->message(error.c_str());
                 found_error = true;
            } else if (!String::IsEqual((*pos)[0], mpos->second.first, Char::eCASE_INSENSITIVE)) {
                 std::string error = "Residue " + (*pos)[0] + " " + (*pos)[1] + " does not exist chain " + _PDB_ChainID + ".\n";
                 _logIo->message(error.c_str());
                 found_error = true;
            }

            if (mpos != number_mapping.end()) end_points.push_back(mpos->second.second);

            mpos = number_mapping.find((*pos)[3]);
            if (mpos == number_mapping.end()) {
                 std::string error = "Residue " + (*pos)[2] + " " + (*pos)[3] + " does not exist chain " + _PDB_ChainID + ".\n";
                 _logIo->message(error.c_str());
                 found_error = true;
            } else if (!String::IsEqual((*pos)[2], mpos->second.first, Char::eCASE_INSENSITIVE)) {
                 std::string error = "Residue " + (*pos)[2] + " " + (*pos)[3] + " does not exist chain " + _PDB_ChainID + ".\n";
                 _logIo->message(error.c_str());
                 found_error = true;
            }
       }
       if (found_error) return;

       std::vector<std::pair<int, int> > split_ranges, non_empty_split_ranges;
       split_ranges.clear();

       for (unsigned int i = 0; i < end_points.size(); ++i) {
            if (i == 0)
                 split_ranges.push_back(std::make_pair(0, end_points[i]));
            else split_ranges.push_back(std::make_pair(end_points[i - 1] + 1, end_points[i]));
       }

       split_ranges.push_back(std::make_pair(end_points[end_points.size() - 1] + 1, _SeqRes.size() - 1));

       if (split_ranges.size() != chain_info.size()) {
            std::string error = "Spliting chain " + _PDB_ChainID + "failed.\n";
            _logIo->message(error.c_str());
            return;
       }

       non_empty_split_ranges.clear();
       for (unsigned int i = 0; i < split_ranges.size(); ++i) {
            bool found_coord_residue_flag = false;
            for (int j = split_ranges[i].first; j <= split_ranges[i].second; ++j) {
                 if (_SeqRes[j]->ResIndex >= 0) {
                      found_coord_residue_flag = true;
                      break;
                 }
            }
            if (!found_coord_residue_flag) continue;
            non_empty_split_ranges.push_back(std::make_pair(split_ranges[i].first, split_ranges[i].second));
       }
       if (non_empty_split_ranges.empty()) {
            std::string error = "Spliting chain " + _PDB_ChainID + "failed.\n";
            _logIo->message(error.c_str());
            return;
       }

       split_ranges = non_empty_split_ranges;

       std::map<int, std::set<unsigned int> > index_map;
       index_map.clear();
       std::vector<int> index;
       index.clear();
       std::set<unsigned int> t_set;

       for (unsigned int i = 0; i < split_ranges.size(); ++i) {
            index.push_back(i);
            int size = split_ranges[i].second - split_ranges[i].first + 1;
            std::map<int, std::set<unsigned int> >::iterator mpos = index_map.find(size);
            if (mpos != index_map.end()) mpos->second.insert(i);
            else { t_set.clear(); t_set.insert(i); index_map.insert(std::make_pair(size, t_set)); }
       }

       idx = 0;
       for (std::map<int, std::set<unsigned int> >::const_reverse_iterator mpos = index_map.rbegin(); mpos != index_map.rend(); ++mpos) {
            for (std::set<unsigned int>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                 index[*spos] = idx;
                 idx++;
            }
       }

       for (unsigned int i = 0; i < split_ranges.size(); ++i) {
            int size = split_ranges[i].second - split_ranges[i].first + 1;
/*
            if ((size == 1) && (_SeqRes[split_ranges[i].first]->ResIndex < 0)) {
                 delete _SeqRes[split_ranges[i].first];
                 _SeqRes[split_ranges[i].first] = NULL;
                 continue;
            }
*/
            Chain* chain = new Chain;

            chain->setCCDic(_ccDic);
            chain->setLog(_logIo);
            chain->setMessage(_messageIo);

            chain->set_ChainID(chain_info[index[i]][1]);
            chain->set_PDB_ChainID(chain_info[index[i]][2]);
            chain->set_PUB_ChainID(chain_info[index[i]][2]);
            chain->set_PDB_ChainID_Flag(_PDB_ChainID_Flag);
            chain->set_details(_details);
            chain->set_prev_entity_id(chain_info[index[i]][0]);
            chain->set_entity_id(chain_info[index[i]][0]);
            chain->set_chain_type(_chain_type);
            chain->set_index(chain_index);
            chain->set_order(chain_index);
            chain_index++;
            chain->set_na_type(_na_type);
            chain->_has_sequence = true;

            chain->_SeqRes.clear();
            chain->_SeqRes.reserve(size);
            chain->_ResidueNumbers = 0;
            std::string residue_type = "HETAIN";
            for (int j = split_ranges[i].first; j <= split_ranges[i].second; ++j) {
                 chain->_SeqRes.push_back(_SeqRes[j]);
                 if (size == 1) residue_type = _ccDic->find_residue_type(_SeqRes[j]->Field[0]);
                 if (_SeqRes[j]->ResIndex >= 0) chain->_ResidueNumbers++;
                 _SeqRes[j] = NULL;
            }

            if (size == 1) {
                 if (residue_type == "ATOMP" || residue_type == "ATOMN") residue_type = "HETAIN";
                 chain->set_chain_type(residue_type);
                 chain->set_na_type(0);
                 chain->_has_sequence = false;
            }

            chain->_residues.clear();
            chain->_residues.reserve(chain->_ResidueNumbers);
            for (std::vector<_FIELD*>::iterator pos = chain->_SeqRes.begin(); pos != chain->_SeqRes.end(); ++pos) {
                 if ((*pos)->ResIndex < 0) continue;
                 int j = (*pos)->ResIndex;
                 (*pos)->ResIndex = chain->_residues.size();
                 for (std::vector<RCSB::Residue*>::iterator vpos = _residues[j].begin(); vpos != _residues[j].end(); ++vpos) {
                      (*vpos)->update_pdb_chnid(chain_info[index[i]][2]);
                 }
                 chain->_residues.push_back(_residues[j]);
            }
            
            chain->renumbering_polymer_chain(true);
            chain->update_indices();
            split_chains.push_back(chain);
       }
}

void Chain::get_nonpolymer_chains(const std::vector<std::vector<std::string> >& deletes, int& chain_index, std::vector<Chain*>& nonpolymer_chains)
{
       nonpolymer_chains.clear();

       std::map<int, std::vector<std::string> > delete_mapping;
       delete_mapping.clear();

       for (std::vector<std::vector<std::string> >::const_iterator vpos = deletes.begin(); vpos != deletes.end(); ++vpos) {
            int idx = atoi((*vpos)[1].c_str()) - 1;
            delete_mapping.insert(std::make_pair(idx, *vpos));
       }

       bool found_error = false;
       for (std::map<int, std::vector<std::string> >::const_iterator mpos = delete_mapping.begin(); mpos != delete_mapping.end(); ++mpos) {
            if ((mpos->first >= 0) && (mpos->first < (int) _SeqRes.size())) {
                 if (mpos->second[0] == _SeqRes[mpos->first]->Field[0]) continue;
            }

            std::string error = "Residue " + mpos->second[0] + " " + mpos->second[1] + " does not exist chain " + _PDB_ChainID + ".\n";
            _logIo->message(error.c_str());
            found_error = true;
       }
       if (found_error) return;

       std::vector<_FIELD*> tmp_seqres;
       tmp_seqres.clear();
       int idx = -1;
       for (std::vector<_FIELD*>::const_iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            idx++;
            std::map<int, std::vector<std::string> >::const_iterator mpos = delete_mapping.find(idx);
            if (mpos != delete_mapping.end()) {
                 if ((*pos)->ResIndex < 0) continue;

                 Chain* chain = new Chain;

                 chain->setCCDic(_ccDic);
                 chain->setLog(_logIo);
                 chain->setMessage(_messageIo);

                 chain->set_ChainID(mpos->second[2]);
                 chain->set_PDB_ChainID(mpos->second[3]);
                 chain->set_PUB_ChainID(mpos->second[3]);
                 chain->set_PDB_ChainID_Flag(_PDB_ChainID_Flag);

                 chain->set_prev_entity_id(mpos->second[4]);
                 chain->set_entity_id(mpos->second[4]);

                 chain->set_chain_type("HETAIN");

                 chain->set_index(chain_index);
                 chain->set_order(chain_index);
                 chain_index++;

                 chain->_residues.clear();
                 chain->_residues.push_back(_residues[(*pos)->ResIndex]);

                 chain->_has_sequence = false;
                 chain->_SeqRes.clear();
                 chain->_SeqRes.push_back(*pos);
                 chain->_SeqRes[0]->ResIndex = 0;
                 chain->_ResidueNumbers = 1;

                 chain->update_indices();

                 nonpolymer_chains.push_back(chain);
            } else tmp_seqres.push_back(*pos);
       }
       _SeqRes = tmp_seqres;

       update_indices();
       _get_entity_key();
}

bool Chain::Update_SeqRes(const std::string& res_name, const std::string& res_num, const std::string& ins_code, const std::string& new_name)
{
       std::string cs = CompositeIndex::getIndex(_PDB_ChainID, res_name, res_num, ins_code);
       std::map<std::string, std::pair<int, int> >::iterator pos = _pdbIndex.find(cs);
       if (pos == _pdbIndex.end() || pos->second.first < 0 || pos->second.second < 0 || _SeqRes[pos->second.first]->ResIndex < 0) {
            std::string error = "Residue (" + res_name + " " + res_num + ins_code + ") does not eixt in chain " + _PDB_ChainID + ".\n";
            _logIo->message(error.c_str());
            return false;
       }

       // if (_SeqRes[pos->second.first]->Field[0] == res_name) {
            _SeqRes[pos->second.first]->Field[0] = new_name;
            _SeqRes[pos->second.first]->Field[1] = new_name;
            _SeqRes[pos->second.first]->Field[3] = new_name;
       // }

       _entity_key.clear();

       return true;
}

bool Chain::Update_SeqRes(const std::string& pdb_chnid, const std::string& res_name, const std::string& res_num, const std::string& ins_code,
                          const std::string& new_name)
{
       std::string cs = CompositeIndex::getIndex(pdb_chnid, res_name, res_num, ins_code);
       std::map<std::string, std::pair<int, int> >::iterator pos = _pdbIndex.find(cs);
       if (pos == _pdbIndex.end() || pos->second.first < 0 || pos->second.second < 0 || _SeqRes[pos->second.first]->ResIndex < 0) {
            std::string error = "Residue (" + res_name + " " + res_num + ins_code + ") does not eixt in chain " + _PDB_ChainID + ".\n";
            _logIo->message(error.c_str());
            return false;
       }

       _SeqRes[pos->second.first]->Field[0] = new_name;
       _SeqRes[pos->second.first]->Field[1] = new_name;
       _SeqRes[pos->second.first]->Field[3] = new_name;

       _entity_key.clear();

       return true;
}

bool Chain::Update_SeqRes(const std::string& res_name, const std::string& res_num, const std::string& ins_code, const std::map<std::string, std::string>& mapping)
{
       std::string cs = CompositeIndex::getIndex(_PDB_ChainID, res_name, res_num, ins_code);
       std::map<std::string, std::pair<int, int> >::iterator pos = _pdbIndex.find(cs);
       if (pos == _pdbIndex.end() || pos->second.first < 0 /* || _SeqRes[pos->second.first]->ResIndex >= 0*/ ) return false;

       std::map<std::string, std::string>::const_iterator mpos = mapping.find("number");
       if (mpos != mapping.end()) _SeqRes[pos->second.first]->Field[4] = mpos->second;

       mpos = mapping.find("inscode");
       if (mpos != mapping.end()) _SeqRes[pos->second.first]->InsCode = mpos->second;

       return true;
}

bool Chain::_IsConnect(RCSB::Residue *res1, RCSB::Residue *res2)
{
       RCSB::Atom *atom1 = res1->GetFirstAtom();
       RCSB::Atom *atom2 = res2->GetFirstAtom();
       if (res1->is_connect(res2))
            return true;
       else if (res1->AtomNumbers() == 1 && res2->AtomNumbers() == 1) {
            if (cal_distance(atom1, atom2) < 4.5)
                 return true;
            else if (abs(atoi(res1->pdb_res_no().c_str()) - atoi(res2->pdb_res_no().c_str())) <= 1)
                 return true;
       }
       return false;
}

void Chain::_update_seqres()
{
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if (*pos) delete *pos;
       }
       _SeqRes.clear();
       _SeqRes.reserve(_residues.size());

       int idx = 0;
       for (std::vector<std::vector<RCSB::Residue*> >::const_iterator vpos = _residues.begin(); vpos != _residues.end(); ++vpos) {
            RCSB::Residue* res = (*vpos)[0];
            _FIELD *seqres = new _FIELD;
            seqres->ResIndex = idx;

            RCSB::Atom *atom = res->GetFirstAtom();
            seqres->InsCode  = atom->ins_code();
            seqres->Field[0] = atom->restype();
            seqres->Field[1] = atom->restype();
            seqres->Field[2] = atom->resnum();
            seqres->Field[3] = atom->pdb_resnam();
            seqres->Field[4] = atom->pdb_resnum();
            seqres->Field[5] = atom->pub_resnam();
            seqres->Field[6] = atom->pub_resnum();

            _SeqRes.push_back(seqres);
            idx++;
       }
       update_indices();
}

void Chain::update_nomenclature(const std::list<std::pair<int, std::string> >& idx_list, std::list<std::list<RCSB::Residue*> >& res_list)
{
       res_list.clear();

       // simply changing numbering
       if (idx_list.empty()) {
            update_indices();
            return;
       }

       std::set<int> remove_set;
       remove_set.clear();
       std::string new_chain_id = idx_list.front().second;
       bool same_chain_id = true;
       for (std::list<std::pair<int, std::string> >::const_iterator pos = idx_list.begin(); pos != idx_list.end(); ++pos) {
            if (pos->second != new_chain_id) same_chain_id = false;
            remove_set.insert(pos->first);
       }

       // simply changing chain ID
       if (idx_list.size() == _SeqRes.size() && _residues[0][0]->ResName() != "HOH" && _residues[0][0]->ResName() != "DOD" && same_chain_id) {
            update_indices();
            _PDB_ChainID = new_chain_id;
            return;
       }

       std::map<std::string, std::list<RCSB::Residue*> > residue_mapping;
       residue_mapping.clear();

       std::vector<RCSB::Residue*> tmp_vec;
       std::list<RCSB::Residue*> tmp_list;

       for (std::vector<_FIELD*>::iterator fpos = _SeqRes.begin(); fpos != _SeqRes.end(); ++fpos) {
            int ResIndex = (*fpos)->ResIndex;
            if (ResIndex < 0) continue;

            tmp_vec.clear();
            for (std::vector<RCSB::Residue*>::const_iterator rpos = _residues[ResIndex].begin(); rpos != _residues[ResIndex].end(); ++rpos) {
                 if (remove_set.find((*rpos)->index()) != remove_set.end()) {
                      std::map<std::string, std::list<RCSB::Residue*> >::iterator mpos = residue_mapping.find((*rpos)->pdb_chnid());
                      if (mpos != residue_mapping.end())
                           mpos->second.push_back(*rpos);
                      else {
                           tmp_list.clear();
                           tmp_list.push_back(*rpos);
                           residue_mapping.insert(std::make_pair((*rpos)->pdb_chnid(), tmp_list));
                      }
                 } else tmp_vec.push_back(*rpos);
            }
            _residues[ResIndex] = tmp_vec;
            if (tmp_vec.empty()) {
                 delete *fpos;
                 *fpos = NULL;
            }
       }

       std::vector<_FIELD*> tmp_seqres;
       tmp_seqres.clear();
       for (std::vector<_FIELD*>::iterator fpos = _SeqRes.begin(); fpos != _SeqRes.end(); ++fpos) {
            if (*fpos == NULL) continue;
            tmp_seqres.push_back(*fpos);
       }
       _SeqRes = tmp_seqres;

       update_indices();

       for (std::map<std::string, std::list<RCSB::Residue*> >::const_iterator mpos = residue_mapping.begin(); mpos != residue_mapping.end(); ++mpos) {
            res_list.push_back(mpos->second);
       }
}

void Chain::merge_waters(const std::list<RCSB::Residue*>& res_list)
{
       if (res_list.empty()) return;

       std::multimap<int, std::vector<RCSB::Residue*> > index_mapping;
       index_mapping.clear();

       for (std::vector<std::vector<RCSB::Residue*> >::const_iterator vpos = _residues.begin(); vpos != _residues.end(); ++vpos) {
            RCSB::Residue* res = (*vpos)[0];
            index_mapping.insert(std::make_pair(atoi(res->pdb_res_no().c_str()), *vpos));
       }

       std::vector<RCSB::Residue*> t_vec;
       for (std::list<RCSB::Residue*>::const_iterator lpos = res_list.begin(); lpos != res_list.end(); ++lpos) {
            t_vec.clear();
            t_vec.push_back(*lpos);
            index_mapping.insert(std::make_pair(atoi((*lpos)->pdb_res_no().c_str()), t_vec));
       }

       _residues.clear();
       _residues.reserve(index_mapping.size());
       for (std::multimap<int, std::vector<RCSB::Residue*> >::const_iterator mpos = index_mapping.begin(); mpos != index_mapping.end(); ++mpos) {
            _residues.push_back(mpos->second);
       }

       _update_seqres();
}

bool Chain::Update(const std::string& resname, const std::string& resnum, const std::string& ins_code, const std::vector<RCSB::Residue*>& res_list)
{
       if (res_list.empty()) return true;

       std::string cs = CompositeIndex::getIndex(_PDB_ChainID, resname, resnum, ins_code);
       std::map<std::string, std::pair<int, int> >::const_iterator mpos = _pdbIndex.find(cs);
       if (mpos == _pdbIndex.end() || mpos->second.first < 0 || mpos->second.second < 0 || _SeqRes[mpos->second.first]->ResIndex < 0) {
            std::string error = "Residue " + resname + " " + resnum + ins_code + " does not exist chain " + _PDB_ChainID + ".\n";
            _logIo->message(error.c_str());
            return false;
       }

       int res_idx = mpos->second.first;
       int pos_idx = mpos->second.second;
       int cif_number = mpos->second.first + 1;
       int pdb_number = atoi(resnum.c_str());

       // vector[0]: _ChainID
       // vector[1]: _PDB_ChainID
       // vector[2]: _res_no
       // vector[3]: _pdb_res_no
       // vector[4]: _ins_code
       std::vector<std::vector<std::string> > nomenclature;
       nomenclature.clear();

       std::vector<std::string> data;
       std::string insCode;

       if (res_idx == 0) {
            pdb_number -= (int) res_list.size();
            for (unsigned int i = 0; i < res_list.size() - 1; ++i) {
                 cif_number++;
                 std::string cs = CompositeIndex::getIndex(String::IntToString(pdb_number + 1), "");
                 if (_numIndex.find(cs) == _numIndex.end()) {
                      pdb_number++;
                      data.clear();
                      data.push_back(_ChainID);
                      data.push_back(_PDB_ChainID);
                      data.push_back(String::IntToString(cif_number));
                      data.push_back(String::IntToString(pdb_number));
                      data.push_back("");
                      nomenclature.push_back(data);
                 } else {
                      char icode = 'A';
                      while (true) {
                           insCode.clear();
                           insCode += icode;
                           std::string cs = CompositeIndex::getIndex(String::IntToString(pdb_number), insCode);
                           if ((_numIndex.find(cs) == _numIndex.end()) || (icode == 'Z')) {
                                data.clear();
                                data.push_back(_ChainID);
                                data.push_back(_PDB_ChainID);
                                data.push_back(String::IntToString(cif_number));
                                data.push_back(String::IntToString(pdb_number));
                                data.push_back(insCode);
                                nomenclature.push_back(data);
                                break;
                           }
                           icode++;
                      }
                 }
            }

            data.clear();
            data.push_back(_ChainID);
            data.push_back(_PDB_ChainID);
            data.push_back(String::IntToString(cif_number));
            data.push_back(resnum);
            data.push_back(ins_code);
            nomenclature.push_back(data);
       } else {
            data.clear();
            data.push_back(_ChainID);
            data.push_back(_PDB_ChainID);
            data.push_back(String::IntToString(cif_number));
            data.push_back(resnum);
            data.push_back(ins_code);
            nomenclature.push_back(data);

            for (unsigned int i = 1; i < res_list.size(); ++i) {
                 cif_number++;
                 std::string cs = CompositeIndex::getIndex(String::IntToString(pdb_number + 1), "");
                 if (_numIndex.find(cs) == _numIndex.end()) {
                      pdb_number++;
                      data.clear();
                      data.push_back(_ChainID);
                      data.push_back(_PDB_ChainID);
                      data.push_back(String::IntToString(cif_number));
                      data.push_back(String::IntToString(pdb_number));
                      data.push_back("");
                      nomenclature.push_back(data);
                 } else {
                      char icode = 'A';
                      while (true) {
                           insCode.clear();
                           insCode += icode;
                           std::string cs = CompositeIndex::getIndex(String::IntToString(pdb_number), insCode);
                           if ((_numIndex.find(cs) == _numIndex.end()) || (icode == 'Z')) {
                                data.clear();
                                data.push_back(_ChainID);
                                data.push_back(_PDB_ChainID);
                                data.push_back(String::IntToString(cif_number));
                                data.push_back(String::IntToString(pdb_number));
                                data.push_back(insCode);
                                nomenclature.push_back(data);
                                break;
                           }
                           icode++;
                      }
                 }
            }
       }

       std::list<RCSB::Atom*> atom_H_list, atom_H2_list;
       for (unsigned int i = 0; i < res_list.size(); ++i) {
            atom_H_list.clear();
            atom_H2_list.clear();
            RCSB::Atom *atom = res_list[i]->GetFirstAtom();
            while (atom) {
                 if (SeqCodeUtil::is_standard_aa_residue(res_list[i]->ResName()) || SeqCodeUtil::is_standard_Daa_residue(res_list[i]->ResName())) { 
                      if (atom->atmtype() == "H") atom_H_list.push_back(atom);
                      if (atom->atmtype() == "H2") atom_H2_list.push_back(atom);
                 }
                 atom->set_chnid(nomenclature[i][0]);
                 atom->set_pdb_chnid(nomenclature[i][1]);
                 atom->set_resnum(nomenclature[i][2]);
                 atom->set_pdb_resnum(nomenclature[i][3]);
                 atom->set_ins_code(nomenclature[i][4]);
                 atom = res_list[i]->GetNextAtom();
            }
            if (atom_H_list.empty() && !atom_H2_list.empty()) {
                 for (std::list<RCSB::Atom*>::iterator lpos = atom_H2_list.begin(); lpos != atom_H2_list.end(); ++lpos) {
                      (*lpos)->set_atmtype("H");
                      (*lpos)->set_pdb_atmnam("H");
                 }
            }
            res_list[i]->set_chnid(nomenclature[i][0]);
            res_list[i]->set_pdb_chnid(nomenclature[i][0]);
            res_list[i]->set_res_no(nomenclature[i][2]);
            res_list[i]->set_pdb_res_no(nomenclature[i][3]);
            res_list[i]->set_ins_code(nomenclature[i][4]);
            res_list[i]->set_chain_index(_index);
       }

       _SeqRes[res_idx]->Field[0] = res_list[0]->ResName();
       _SeqRes[res_idx]->Field[1] = res_list[0]->ResName();
       _SeqRes[res_idx]->Field[3] = res_list[0]->ResName(); 
       _SeqRes[res_idx]->Field[4] = res_list[0]->pdb_res_no();
       _SeqRes[res_idx]->InsCode = res_list[0]->ins_code();
       _residues[_SeqRes[res_idx]->ResIndex][pos_idx] = res_list[0];
       if (res_list.size() == 1) {
            cs = CompositeIndex::getIndex(_PDB_ChainID, resname, resnum, ins_code);
            _pdbIndex.erase(cs);
            cs = CompositeIndex::getIndex(_PDB_ChainID, res_list[0]->ResName(), resnum, ins_code);
            _pdbIndex.insert(std::make_pair(cs, std::make_pair(res_idx, pos_idx)));
            return true;
       }

       std::vector<_FIELD*> tmp_seqres = _SeqRes;
       _SeqRes.clear();
       _SeqRes.reserve(tmp_seqres.size() + res_list.size());

       bool need_renumber = false;
       std::set<std::string> index_set;
       index_set.clear();

       for (int i = 0; i <= res_idx; ++i) {
            _SeqRes.push_back(tmp_seqres[i]);

            cs = tmp_seqres[i]->Field[4] + "_" + tmp_seqres[i]->InsCode;
            if (index_set.find(cs) != index_set.end()) need_renumber = true;
            index_set.insert(cs);
       }
       for (unsigned int i = 1; i < res_list.size(); ++i) {
            insert_a_residue(res_list[i], true);

            cs = res_list[i]->pdb_res_no() + "_" + res_list[i]->ins_code();
            if (index_set.find(cs) != index_set.end()) need_renumber = true;
            index_set.insert(cs);
       }
       for (unsigned int i = res_idx + 1; i < tmp_seqres.size(); ++i) {
            _SeqRes.push_back(tmp_seqres[i]);

            cs = tmp_seqres[i]->Field[4] + "_" + tmp_seqres[i]->InsCode;
            if (index_set.find(cs) != index_set.end()) need_renumber = true;
            index_set.insert(cs);
       }
       tmp_seqres.clear();

       if (need_renumber) check_unique_number();

       update_indices();
 
       return true;
}

void Chain::reorder_residues(const std::map<int, std::string>& index_mapping)
{
       std::map<std::string, int> pos_index_mapping;
       pos_index_mapping.clear();
       for (unsigned int i = 0; i < _residues.size(); ++i) {
            std::string idx = CompositeIndex::getIndex(_residues[i][0]->pdb_chnid(), _residues[i][0]->ResName(), _residues[i][0]->pdb_res_no(),
                                                       _residues[i][0]->ins_code());
            pos_index_mapping.insert(std::make_pair(idx, (int) i));
       }

       std::vector<std::vector<RCSB::Residue*> > new_residues;
       new_residues.clear();
       new_residues.reserve(_residues.size());
       for (std::map<int, std::string>::const_iterator impos = index_mapping.begin(); impos != index_mapping.end(); ++impos) {
            std::map<std::string, int>::const_iterator pmpos = pos_index_mapping.find(impos->second);
            if (pmpos != pos_index_mapping.end()) new_residues.push_back(_residues[pmpos->second]);
       }
       if (new_residues.size() != _residues.size()) return;

       _residues = new_residues;
       _update_seqres();
}

bool Chain::found_missing_single_residues()
{
       for (unsigned int i = 0; i < _SeqRes.size(); ++i) {
            if (_SeqRes[i]->ResIndex >= 0) continue;

            if (i == 0) {
                 if (_SeqRes[1]->ResIndex >= 0) return true;
            } else if (i == (_SeqRes.size() - 1)) {
                 if (_SeqRes[i - 1]->ResIndex >= 0) return true;
            } else {
                 if ((_SeqRes[i - 1]->ResIndex >= 0) && (_SeqRes[i + 1]->ResIndex >= 0)) return true;
            }
       }

       return false;
}

void Chain::check_missing_misplace_residues(const std::map<std::string, std::vector<RCSB::Chain*> >& res_name_chain_mapping,
                                            std::set<int>& removed_chain_index_set)
{
       std::vector<std::vector<RCSB::Residue*> > tmp_residues;
       tmp_residues.clear();

       std::vector<_FIELD*> tmp_SeqRes;
       tmp_SeqRes.clear();

       bool found_misplace_residue = false;
       bool need_renumber = _check_need_renumber();

       for (unsigned int i = 0; i < _SeqRes.size(); ++i) {
            _FIELD *seqres = new _FIELD;
            for (int j = 0; j < NUM_FIELD; ++j) seqres->Field[j] = _SeqRes[i]->Field[j];
            seqres->InsCode = _SeqRes[i]->InsCode;
            seqres->ResIndex = _SeqRes[i]->ResIndex;

            if ((_SeqRes[i]->ResIndex >= 0) || _SeqRes[i]->Field[0].empty()) {
                 if (_SeqRes[i]->ResIndex >= 0) {
                      int new_idx = tmp_residues.size();
                      tmp_residues.push_back(_residues[_SeqRes[i]->ResIndex]);
                      seqres->ResIndex = new_idx;
                 }
                 tmp_SeqRes.push_back(seqres);
                 continue;
            }

            int found = -1;
            if (i == 0) {
                 if (_SeqRes[1]->ResIndex >= 0) found = 0;
            } else if (i == (_SeqRes.size() - 1)) {
                 if (_SeqRes[i - 1]->ResIndex >= 0) found = 1;
            } else {
                 if ((_SeqRes[i - 1]->ResIndex >= 0) && (_SeqRes[i + 1]->ResIndex >= 0)) found = 2;
            }
            if (found < 0) {
                 tmp_SeqRes.push_back(seqres);
                 continue;
            }

            std::map<std::string, std::vector<RCSB::Chain*> >::const_iterator mpos = res_name_chain_mapping.find(_SeqRes[i]->Field[0]);
            if (mpos == res_name_chain_mapping.end()) {
                 tmp_SeqRes.push_back(seqres);
                 continue;
            }

            for (std::vector<RCSB::Chain*>::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 RCSB::Chain* chain = *vpos;
                 if ((chain->_SeqRes.size() > 1) || (chain->_residues.size() != 1)) continue;
                 bool is_linked = false;
                 if (found == 0) {
                     if (is_polymer_connect(_chain_type, chain->_residues[0][0], _residues[_SeqRes[i + 1]->ResIndex][0])) is_linked = true;
                 } else if (found == 1) {
                     if (is_polymer_connect(_chain_type, _residues[_SeqRes[i - 1]->ResIndex][0], chain->_residues[0][0])) is_linked = true;
                 } else if (is_polymer_connect(_chain_type, _residues[_SeqRes[i - 1]->ResIndex][0], chain->_residues[0][0]) &&
                            is_polymer_connect(_chain_type, chain->_residues[0][0], _residues[_SeqRes[i + 1]->ResIndex][0])) is_linked = true;
                 if (!is_linked) continue;

                 found_misplace_residue = true;
                 removed_chain_index_set.insert(chain->_index);

                 seqres->ResIndex = tmp_residues.size();
                 if (!need_renumber) {
                      for (std::vector<RCSB::Residue*>::iterator rpos = chain->_residues[0].begin(); rpos != chain->_residues[0].end(); ++rpos) {
                           (*rpos)->update_nomenclature(_chain_type, _ChainID, _SeqRes[i]->Field[2], _PDB_ChainID, _SeqRes[i]->Field[4], _SeqRes[i]->InsCode);
                      }
                 }
                 tmp_residues.push_back(chain->_residues[0]);

                 RCSB::Atom *atom = chain->_residues[0][0]->GetFirstAtom();
                 seqres->InsCode  = atom->ins_code();
                 seqres->Field[1] = atom->restype();
                 seqres->Field[2] = atom->resnum();
                 seqres->Field[3] = atom->pdb_resnam();
                 seqres->Field[4] = atom->pdb_resnum();
                 seqres->Field[5] = atom->pub_resnam();
                 seqres->Field[6] = atom->pub_resnum();
                 
                 break;
            }
            tmp_SeqRes.push_back(seqres);
       }

       if (!found_misplace_residue) {
            for (std::vector<_FIELD*>::iterator pos = tmp_SeqRes.begin(); pos != tmp_SeqRes.end(); ++pos) {
                 if (*pos) delete *pos;
            }
            tmp_residues.clear();
            tmp_SeqRes.clear();
            return;
       }

       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if (*pos) delete *pos;
       }
    
       _residues = tmp_residues;
       _SeqRes = tmp_SeqRes;

       if (_check_need_renumber()) check_unique_number();

       update_indices();
}

bool Chain::_check_need_renumber()
{
       bool need_renumber = false;
       std::set<std::string> UniqueIndex;
       UniqueIndex.clear();
       for (std::vector<_FIELD*>::iterator pos = _SeqRes.begin(); pos != _SeqRes.end(); ++pos) {
            if ((*pos)->ResIndex < 0 && (*pos)->Field[4].empty()) {
                 need_renumber = true;
            }
            std::string cs = (*pos)->Field[4] + (*pos)->InsCode;
            if (UniqueIndex.find(cs) != UniqueIndex.end()) {
                 need_renumber = true;
            }
            UniqueIndex.insert(cs);
       }
       return need_renumber;
}

static double score(const int& i, const int& j, void* data)
{
       _DataContainer *cdata = (_DataContainer*) data;
       if (!SeqCodeUtil::isSameResName(cdata->_seqa[i], cdata->_seqb[j]))
            return 0;
       else if ((cdata->_seqa[i] == "ALA") && (cdata->_seqb[j] == "ALA"))
            return 2;
       else return 5;
}
