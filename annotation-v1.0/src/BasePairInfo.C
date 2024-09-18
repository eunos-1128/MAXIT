/*
FILE:     BasePairInfo.C
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

#include "BasePairInfo.h"
#include "BaseParameter.h"
#include "BasePairParameter.h"
#include "CompositeIndex.h"
#include "MatrixUtil.h"
#include "PiaAreaUtil.h"
#include "SeqCodeUtil.h"
#include "VectorUtil.h"
#include "utillib.h"
#include "xtal.h"

#define MIN_CUTOFF 2.2
#define MAX_CUTOFF 3.8

#define MAX_DORG  10.0
#define MAX_DPRJ  2.5
#define MAX_ANGLE 65.0
#define MIN_N91   4.5
#define XBIG      1.0e+18
#define OVERLAP   0.01

#define NUM_BASE_PAIR_TYPE      30
#define NUM_HYDROGEN_BOND_ATOMS  9

BasePairInfo::BasePairInfo()
{
       _clear();
}

BasePairInfo::~BasePairInfo()
{
       _clear();
}

void BasePairInfo::_clear()
{
       _cell = NULL;
       _ccDic = NULL;
       _mol = NULL;
       _classification = 0;
       _paired_chain_id_map.clear();
       _basepair_contacts.clear();
       _basebase_params.clear();
       _interbase_params.clear();
       _baseframes.clear();
       _basepairs.clear();
}

void BasePairInfo::setCell(CrySymmetry* cell)
{
       _cell = cell;
}

void BasePairInfo::setCCDic(ConnectDic* ccdic)
{
       _ccDic = ccdic;
}

void BasePairInfo::setMolecule(RCSB::Molecule* mol)
{
       _mol = mol;
}

void BasePairInfo::calculateBasePairInfo()
{
       if (_mol == NULL || _cell == NULL) return;

       std::set<std::string> chain_types, residue_set;
       chain_types.clear();
       chain_types.insert("ATOMN");
       residue_set.clear();

       xtal mycrys;
       if (!set_mycrys(mycrys, *_cell, _mol, chain_types, residue_set, true, false, false)) return;

       BasePairParameter::initialize();

       std::vector<CONTACT> contact;
       get_contact(mycrys, contact, MIN_CUTOFF, MAX_CUTOFF, "asym", "all_asym", "all_asym", 0, 1);
       int atomn_residues = _getBaseFrames(contact, "", 0, 0, 0, 0);
       int residue_count = _getBasePairInformation(contact, "", true);
       _getBestBasePair();
       _re_ordering();

       if (_cell->is_artifical() || ((residue_count * 2) > atomn_residues)) return;

       std::vector<std::pair<std::string, std::vector<int> > > chnid_symop_list;
       get_contact(mycrys, contact, MIN_CUTOFF, MAX_CUTOFF, "sym", "all_asym", "all_sym", 0, 0);
       _getSymOperator(contact, chnid_symop_list);
       for (std::vector<std::pair<std::string, std::vector<int> > >::const_iterator vpos = chnid_symop_list.begin(); vpos != chnid_symop_list.end(); ++vpos) {
            std::map<std::string, std::set<int> >::const_iterator mpos = _paired_chain_id_map.find(vpos->first);
            if ((mpos != _paired_chain_id_map.end()) && ((int) mpos->second.size() > vpos->second[4])) continue;
            if (vpos->second[0] || vpos->second[1] || vpos->second[2] || vpos->second[3]) {
                 _getBaseFrames(contact, vpos->first, vpos->second[0], vpos->second[1], vpos->second[2], vpos->second[3]);
                 _getBasePairInformation(contact, vpos->first, false);
                 _getBestBasePair();
                 _re_ordering();
            }
       }
}

void BasePairInfo::getBasePairs(std::list<_BSPAIR>& bspairs)
{
       bspairs.clear();
       if (_basepair_contacts.empty()) return;

       const char *_base_pair_type[NUM_BASE_PAIR_TYPE] = {
              "",
              "TYPE_1_PAIR",
              "TYPE_2_PAIR",
              "TYPE_3_PAIR",
              "TYPE_4_PAIR",
              "TYPE_5_PAIR",
              "TYPE_6_PAIR",
              "TYPE_7_PAIR",
              "TYPE_8_PAIR",
              "TYPE_9_PAIR",
              "TYPE_10_PAIR",
              "TYPE_11_PAIR",
              "TYPE_12_PAIR",
              "TYPE_13_PAIR",
              "TYPE_14_PAIR",
              "TYPE_15_PAIR",
              "TYPE_16_PAIR",
              "TYPE_17_PAIR",
              "TYPE_18_PAIR",
              "WATSON-CRICK",
              "WATSON-CRICK",
              "REVERSED WATSON-CRICK",
              "REVERSED WATSON-CRICK",
              "HOOGSTEEN",
              "REVERSED HOOGSTEEN",
              "TYPE_25_PAIR",
              "TYPE_26_PAIR",
              "TYPE_27_PAIR",
              "TYPE_28_PAIR",
              "TYPE_29_PAIR"
       };

       _BSPAIR bspair;
       bspair.SymOP_1.clear();
       bspair.type = "hydrog";
       bspair.dist.clear();
       bspair.bondtype.clear();
       bspair.leaving_flag.clear();
       for (std::multimap<BPIndex, CONTACT>::const_iterator mpos = _basepair_contacts.begin(); mpos != _basepair_contacts.end(); ++mpos) {
            bspair.mol_index = mpos->second.mol_index;
            bspair.fstAtom   = mpos->second.a_atm;
            bspair.sndAtom   = mpos->second.b_atm;
            bspair.SymOP_2.clear();
            if (mpos->second.sym > 1 || mpos->second.lx || mpos->second.ly || mpos->second.lz) {
                 bspair.SymOP_2 = String::IntToString(mpos->second.sym) + "_" + String::IntToString(mpos->second.lx + 5)
                                + String::IntToString(mpos->second.ly + 5) + String::IntToString(mpos->second.lz + 5);
            }
            bspair.details.clear();
            if (mpos->second.type > 0 && mpos->second.type < NUM_BASE_PAIR_TYPE)
                 bspair.details = _base_pair_type[mpos->second.type];
            else if ((SeqCodeUtil::isSameResName(mpos->second.a_atm->restype(), "A") && (SeqCodeUtil::isSameResName(mpos->second.b_atm->restype(), "T") ||
                      SeqCodeUtil::isSameResName(mpos->second.b_atm->restype(), "U"))) || ((SeqCodeUtil::isSameResName(mpos->second.a_atm->restype(), "T") ||
                      SeqCodeUtil::isSameResName(mpos->second.a_atm->restype(), "U")) && SeqCodeUtil::isSameResName(mpos->second.b_atm->restype(), "A")) ||
                     (SeqCodeUtil::isSameResName(mpos->second.a_atm->restype(), "G") && SeqCodeUtil::isSameResName(mpos->second.b_atm->restype(), "C")) ||
                     (SeqCodeUtil::isSameResName(mpos->second.a_atm->restype(), "C") && SeqCodeUtil::isSameResName(mpos->second.b_atm->restype(), "G")))
                 bspair.details = mpos->second.a_atm->restype() + "-" + mpos->second.b_atm->restype() + " PAIR";
            else bspair.details = mpos->second.a_atm->restype() + "-" + mpos->second.b_atm->restype() + " MISPAIR";
            bspairs.push_back(bspair);
       }
}

const std::list<_BASEBASE_PARAMS>& BasePairInfo::getBaseBaseParams() const
{
       return _basebase_params;
}

const std::list<_INTERBASE_PARAMS>& BasePairInfo::getInterBaseParams() const
{
       return _interbase_params;
}

int BasePairInfo::_getBaseFrames(const std::vector<CONTACT>& contact, const std::string& selected_chainid, const int& op, const int& lx,
                                 const int& ly, const int& lz)
{
       if (contact.empty()) return 0;

       _baseframes.clear();

       _BASE _base;
       BaseParameter base_parameter;
       base_parameter.setCell(_cell);
       base_parameter.setCCDic(_ccDic);

       int _atomn_residues = 0;
       RCSB::Chain* chain = _mol->GetFirstChain();
       while (chain) {
            if ((chain->chain_type() == "ATOMN") && (selected_chainid.empty() || (selected_chainid == chain->PDB_ChainID()))) {
                 RCSB::Residue* res = chain->GetFirstResidue();
                 while (res) {
                      _atomn_residues++;
                      if (op || lx || ly || lz) {
                           base_parameter.getBaseFrame(_base, res, 0, 0, 0, 0);
                           _baseframes.push_back(_base);
                      }
                      base_parameter.getBaseFrame(_base, res, op, lx, ly, lz);
                      _baseframes.push_back(_base);
                      res = chain->GetNextResidue();
                 }
            }
            chain = _mol->GetNextChain();
       }
       return _atomn_residues;
}

int BasePairInfo::_getBasePairInformation(const std::vector<CONTACT>& contact, const std::string& selected_chainid, const bool& count_paired_flag)
{
       if (_baseframes.empty()) return 0;

       std::vector<std::string> data;
       std::vector<CONTACT> contact_vec;
       std::map<std::string, vector<CONTACT> > contact_index;
       contact_index.clear();
       for (std::vector<CONTACT>::const_iterator pos = contact.begin(); pos != contact.end(); ++pos) {
            data.clear();
            pos->a_atm->getAtomIndex(data, false, true);
            data.push_back(String::IntToString(pos->sym));
            data.push_back(String::IntToString(pos->lx));
            data.push_back(String::IntToString(pos->ly));
            data.push_back(String::IntToString(pos->lz));
            pos->b_atm->getAtomIndex(data, false, true);
            std::string index = CompositeIndex::getIndex(data);
            std::map<std::string, vector<CONTACT> >::iterator mpos = contact_index.find(index);
            if (mpos != contact_index.end()) mpos->second.push_back(*pos);
            else {
                 contact_vec.clear();
                 contact_vec.push_back(*pos);
                 contact_index.insert(std::make_pair(index, contact_vec));
            }
       }

       _basepairs.clear();
       _basepairs.reserve(_baseframes.size());
       std::vector<_BASEPAIR> t_vec;
       t_vec.clear();
       for (unsigned int i = 0; i < _baseframes.size(); ++i) _basepairs.push_back(t_vec);

       _BASEPAIR a_pair;
       int count = 0;
       for (unsigned int i = 0; i < _baseframes.size() - 1; ++i) {
            if (!_baseframes[i].has_base) continue;
            for (unsigned int j = i + 1; j < _baseframes.size(); ++j) {
                 if (!_baseframes[j].has_base) continue;
                 int i1 = i;
                 int j1 = j;
                 if (_baseframes[i].op || _baseframes[i].lx || _baseframes[i].ly || _baseframes[i].lz) {
                      i1 = j; j1 = i;
                 }
                 if (_check_base_pair(_baseframes[i1], _baseframes[j1], a_pair)) {
                      count++;
                      a_pair.i = i1;
                      a_pair.j = j1;
                      BasePairParameter::getParameter(_baseframes[i1], _baseframes[j1], contact_index, a_pair.type, a_pair.type_i, a_pair.type_j,
                                                      a_pair.cis_or_trans, a_pair.bp_contact);
                      if (a_pair.bp_contact.empty()) continue;
                      if (a_pair.type > 0 && a_pair.type < 30) {
                           a_pair.dsum -= 1.0;
                           if (a_pair.type == 19 || a_pair.type == 20) a_pair.dsum -= 1.0;
                      }
                      for (std::vector<CONTACT>::const_iterator pos = a_pair.bp_contact.begin(); pos != a_pair.bp_contact.end(); ++pos) {
                           if (pos->dist < 3.5) a_pair.dsum -= 1.0;
                      }
                      _basepairs[i1].push_back(a_pair);
                 }
            }
       }

       if (!selected_chainid.empty()) {
            std::map<std::string, std::set<int> >::const_iterator mpos = _paired_chain_id_map.find(selected_chainid);
            if (mpos != _paired_chain_id_map.end() && ((int) mpos->second.size() > count)) {
                 _basepairs.clear();
                 return 0;
            }
       }

       return _refineBasePairs(count_paired_flag);
}

int BasePairInfo::_refineBasePairs(const bool& count_paired_flag)
{
       if (_basepairs.empty()) return 0;

       std::set<std::string> residue_count;
       residue_count.clear();

       BPIndex bpidx;
       for (std::vector<std::vector<_BASEPAIR> >::iterator pos = _basepairs.begin(); pos != _basepairs.end(); ++pos) {
            if (pos->empty()) continue;
            _getBasePairs(*pos);
            for (std::vector<_BASEPAIR>::const_iterator ppos = pos->begin(); ppos != pos->end(); ++ppos) {
                 if (count_paired_flag) {
                      _updatePairedChainIdMap(ppos->i);
                      _updatePairedChainIdMap(ppos->j);
                 }
                 for (std::vector<CONTACT>::const_iterator cpos = ppos->bp_contact.begin(); cpos != ppos->bp_contact.end(); ++cpos) {
                      residue_count.insert(cpos->a_atm->pdb_chnid() + "_" + cpos->a_atm->pdb_resnum() + "_" + cpos->a_atm->pdb_resnam());
                      residue_count.insert(cpos->b_atm->pdb_chnid() + "_" + cpos->b_atm->pdb_resnum() + "_" + cpos->b_atm->pdb_resnam());
                      bpidx.set_a_chnid(cpos->a_atm->pdb_chnid());
                      bpidx.set_a_resnum(cpos->a_atm->pdb_resnum());
                      bpidx.set_a_resnam(cpos->a_atm->pdb_resnam());
                      bpidx.set_a_atomtyp(cpos->a_atm->pdb_atmnam());
                      bpidx.set_b_chnid(cpos->b_atm->pdb_chnid());
                      bpidx.set_b_resnum(cpos->b_atm->pdb_resnum());
                      bpidx.set_b_resnam(cpos->b_atm->pdb_resnam());
                      bpidx.set_b_atomtyp(cpos->b_atm->pdb_atmnam());
                      _basepair_contacts.insert(std::make_pair(bpidx, *cpos));
                 }
            }
       }
       return residue_count.size();
}

void BasePairInfo::_updatePairedChainIdMap(const int& idx)
{
       std::map<std::string, std::set<int> >::iterator mpos = _paired_chain_id_map.find(_baseframes[idx].res->pdb_chnid());
       if (mpos != _paired_chain_id_map.end())
            mpos->second.insert(idx);
       else {
            std::set<int> index_count;
            index_count.clear();
            index_count.insert(idx);
            _paired_chain_id_map.insert(std::make_pair(_baseframes[idx].res->pdb_chnid(), index_count));
       }
}

void BasePairInfo::_getBasePairs(std::vector<_BASEPAIR>& base_pairs)
{
       if (base_pairs.size() < 2) return; 

       std::vector<std::vector<_BASEPAIR> > base_pair_clusters;
       base_pair_clusters.clear();

       std::vector<_BASEPAIR> base_pair_cluster;
       base_pair_cluster.clear();
       base_pair_cluster.push_back(base_pairs[0]);
       base_pair_clusters.push_back(base_pair_cluster);

       for (unsigned int i = 1; i < base_pairs.size(); ++i) {
            bool found = false;
            for (std::vector<std::vector<_BASEPAIR> >::iterator pos = base_pair_clusters.begin(); pos != base_pair_clusters.end(); ++pos) {
                 for (std::vector<_BASEPAIR>::const_iterator ppos = pos->begin(); ppos != pos->end(); ++ppos) {
                      if (_is_overlap(_baseframes[ppos->j], _baseframes[base_pairs[i].j])) {
                           found = true;
                           break;
                      }
                 }
                 if (found) {
                      pos->push_back(base_pairs[i]);
                      break;
                 }
            }
            if (found) continue;
            base_pair_cluster.clear();
            base_pair_cluster.push_back(base_pairs[i]);
            base_pair_clusters.push_back(base_pair_cluster);
       }

       std::vector<_BASEPAIR> return_base_pairs;
       return_base_pairs.clear();
       std::set<unsigned int> found_index;
       found_index.clear();
       for (unsigned int i = 0; i < base_pair_clusters.size(); ++i) {
            for (std::vector<_BASEPAIR>::const_iterator pos = base_pair_clusters[i].begin(); pos != base_pair_clusters[i].end(); ++pos) {
                 // remove bad base pairs from clusters if there are any good base pairs
                 if (pos->type > 0 && pos->type < 30) {
                      return_base_pairs.push_back(*pos);
                      found_index.insert(i);
                 }
            }
       }
       if (!return_base_pairs.empty()) {
            for (unsigned int i = 0; i < base_pair_clusters.size(); ++i) {
                 if (found_index.find(i) != found_index.end()) continue; 
                 // remove other bad base pairs from clusters if there are any good base pairs
                 for (std::vector<_BASEPAIR>::const_iterator pos = base_pair_clusters[i].begin(); pos != base_pair_clusters[i].end(); ++pos) {
                      if (pos->bp_contact[0].dist < 3.5) return_base_pairs.push_back(*pos);
                 }
            }
            base_pairs = return_base_pairs;
       }
}

void BasePairInfo::_getBestBasePair()
{
       if (_basepairs.empty()) return;

       for (unsigned int i = 0; i < _basepairs.size() - 1; ++i) {
            for (unsigned int k = 0; k < _basepairs[i].size(); ++k) {
                 if (_basepairs[i][k].type < 0) continue;
                 int j = _basepairs[i][k].j;
                 bool found = false;
                 for (unsigned int l = 0; l < _basepairs[j].size(); ++l) {
                      if (_basepairs[j][l].j == (int) i) {
                           found = true;
                           break;
                      }
                 }
                 if (found) continue;
                 _BASEPAIR a_pair = _basepairs[i][k];
                 a_pair.i = _basepairs[i][k].j;
                 a_pair.j = _basepairs[i][k].i;
                 a_pair.type_i = _basepairs[i][k].type_j;
                 a_pair.type_j = _basepairs[i][k].type_i;
                 a_pair.zave.x = -_basepairs[i][k].zave.x;
                 a_pair.zave.y = -_basepairs[i][k].zave.y;
                 a_pair.zave.z = -_basepairs[i][k].zave.z;
                 _basepairs[j].push_back(a_pair);
            }
       }

       std::vector<int> index;
       index.clear();
       for (unsigned int i = 0; i < _basepairs.size(); ++i) index.push_back(0);

       _BASEPAIR bp_pair_i, bp_pair_j;
       _best_basepairs.clear();
       for (unsigned int i = 0; i < _basepairs.size(); ++i) {
            if (index[i]) continue;
            int j = _getBestPair(_basepairs[i], bp_pair_i, index);
            if (j >= 0) {
                 int k = _getBestPair(_basepairs[j], bp_pair_j, index);
                 if (k == (int) i) {
                      index[i] = 1;
                      index[j] = 1;
                      _best_basepairs.push_back(bp_pair_i);
                 }
            }
       }
}

void BasePairInfo::_getSymOperator(const std::vector<CONTACT>& contact, std::vector<std::pair<std::string, std::vector<int> > >& chnid_symop_list)
{
       chnid_symop_list.clear();

       const char *_hydrogen_bond_atoms[NUM_HYDROGEN_BOND_ATOMS] = { "N1", "N2", "N3", "N4", "N6", "N7", "O2", "O4", "O6" };

       std::set<std::string> hydrogen_bond_atom_set;
       hydrogen_bond_atom_set.clear();
       for (int i = 0; i < NUM_HYDROGEN_BOND_ATOMS; ++i) hydrogen_bond_atom_set.insert(_hydrogen_bond_atoms[i]);

       std::multimap<std::string, std::vector<int> > sym_op_map;
       sym_op_map.clear();
       std::set<std::string> sym_op_set;
       sym_op_set.clear();
       std::vector<int> sym_op;
       for (std::vector<CONTACT>::const_iterator pos = contact.begin(); pos != contact.end(); ++pos) {
            if (pos->sym <= 1 && pos->lx == 0 && pos->ly == 0 && pos->lz == 0) continue;
            if (pos->a_atm->pdb_chnid() != pos->b_atm->pdb_chnid()) continue;
            if (hydrogen_bond_atom_set.find(pos->a_atm->pdb_atmnam()) == hydrogen_bond_atom_set.end() ||
                hydrogen_bond_atom_set.find(pos->b_atm->pdb_atmnam()) == hydrogen_bond_atom_set.end()) continue;
            std::string cs = pos->a_atm->pdb_chnid() + "-" + String::IntToString(pos->sym) + "_" + String::IntToString(pos->lx) + "_"
                           + String::IntToString(pos->ly) + "_" + String::IntToString(pos->lz);
            sym_op.clear();
            sym_op.push_back(pos->sym - 1);
            sym_op.push_back(pos->lx);
            sym_op.push_back(pos->ly);
            sym_op.push_back(pos->lz);
            sym_op_map.insert(std::make_pair(cs, sym_op));
            sym_op_set.insert(cs);
       }
       if (sym_op_map.empty()) return;

       std::vector<std::string> data;
       std::multimap<int, std::vector<int> > tmp_map;

       std::map<std::string, std::multimap<int, std::vector<int> > > chain_id_sym_op_map;
       chain_id_sym_op_map.clear();

       for (std::set<std::string>::const_iterator spos = sym_op_set.begin(); spos != sym_op_set.end(); ++spos) {
            std::multimap<std::string, std::vector<int> >::const_iterator mpos = sym_op_map.find(*spos);
            if (mpos == sym_op_map.end()) continue;

            get_wordarray(data, *spos, "-");
            int count = sym_op_map.count(*spos);
            if (count < 9) continue;

            sym_op.clear();
            for (int i = 0; i < 4; ++i) sym_op.push_back(mpos->second[i]);
            sym_op.push_back(count);

            std::map<std::string, std::multimap<int, std::vector<int> > >::iterator mmpos = chain_id_sym_op_map.find(data[0]);
            if (mmpos != chain_id_sym_op_map.end())
                 mmpos->second.insert(std::make_pair(count, sym_op));
            else {
                 tmp_map.clear();
                 tmp_map.insert(std::make_pair(count, sym_op));
                 chain_id_sym_op_map.insert(std::make_pair(data[0], tmp_map));
            }
       }

       if (chain_id_sym_op_map.empty()) return;

       for (std::map<std::string, std::multimap<int, std::vector<int> > >::const_iterator mpos = chain_id_sym_op_map.begin();
            mpos != chain_id_sym_op_map.end(); ++mpos) {
            std::multimap<int, std::vector<int> >::const_reverse_iterator rpos = mpos->second.rbegin();
            int max_count = rpos->first * 9 / 10;
            for (std::multimap<int, std::vector<int> >::const_iterator mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                 if (mmpos->first < max_count) continue;
                 chnid_symop_list.push_back(std::make_pair(mpos->first, mmpos->second));
            }
       }
}

bool BasePairInfo::_check_base_pair(const _BASE& base_i, const _BASE& base_j, _BASEPAIR& bs_pair)
{
       COORD coord, dorg, zave;
       bs_pair.dir.x = vector_dot_product(base_j.orien[0], base_i.orien[0]); // relative x direction
       bs_pair.dir.y = vector_dot_product(base_j.orien[1], base_i.orien[1]); // relative y direction
       bs_pair.dir.z = vector_dot_product(base_j.orien[2], base_i.orien[2]); // relative z direction
       double angle = 90.0 - fabs(dot2ang(bs_pair.dir.z) - 90.0); // angle between base normals
       double dn = MIN_N91 + 1.0; // RN9-YN1 distance
       std::map<std::string, COORD>::const_iterator mpos_i = base_i.coord_xyz.find("N");
       std::map<std::string, COORD>::const_iterator mpos_j = base_j.coord_xyz.find("N");
       if (mpos_i != base_i.coord_xyz.end() && mpos_j != base_j.coord_xyz.end()) {
            vector_difference(coord, mpos_j->second, mpos_i->second);
            dn = vector_length(coord); // RN9-YN1 distance
       }
       vector_difference(dorg, base_j.org, base_i.org);
       double d = vector_length(dorg); // distance between origins
       if (bs_pair.dir.z <= 0.0)
            vector_difference(bs_pair.zave, base_i.orien[2], base_j.orien[2]);
       else vector_sum(bs_pair.zave, base_i.orien[2], base_j.orien[2]);
       vector_normalize(bs_pair.zave);
       zave.x = -bs_pair.zave.x;
       zave.y = -bs_pair.zave.y;
       zave.z = -bs_pair.zave.z;
       bs_pair.dv = fabs(vector_dot_product(dorg, zave)); // projection onto mean normal
       bs_pair.dsum  = d + 2.0 * bs_pair.dv + angle / 20.0;
       // bs_pair.dsum  = bs_pair.d + 2.0 * bs_pair.dv; // version 1.5
       if (d <= MAX_DORG && bs_pair.dv <= MAX_DPRJ && angle <= MAX_ANGLE && dn >= MIN_N91) {
            vector_middle(bs_pair.oave, base_j.org, base_i.org);
            if (!base_i.polygons.empty() && !base_j.polygons.empty()) {
                 double overlap_area = _get_overlap_area(base_i, base_j, bs_pair.oave, zave);
                 if (overlap_area < OVERLAP) return true;
            } else return true;
       }
       return false;
}

bool BasePairInfo::_is_overlap(const _BASE& base_i, const _BASE& base_j)
{
       if (!base_i.polygons.empty() && !base_j.polygons.empty()) {
            COORD oave, zave;
            vector_middle(oave, base_i.org, base_j.org);
            double z = vector_dot_product(base_i.orien[2], base_j.orien[2]);
            if (z <= 0.0)
                 vector_difference(zave, base_i.orien[2], base_j.orien[2]);
            else vector_sum(zave, base_i.orien[2], base_j.orien[2]);
            if (_get_overlap_area(base_i, base_j, oave, zave) > OVERLAP) return true;
       }
       return false;
}

double BasePairInfo::_get_overlap_area(const _BASE& base_i, const _BASE& base_j, const COORD& oave, const COORD& zave)
{
       COORD z, hinge;
       z.x = 0.0;
       z.y = 0.0;
       z.z = 1.0;
       vector_cross_product(hinge, zave, z);
       double ang_deg = magang(zave, z);

       COORD rot_mtx[3], rot_mtxT[3];
       arb_rotation(rot_mtx, hinge, ang_deg);
       transpose_matrix(rot_mtxT, rot_mtx);

       std::vector<COORD> polygon_i = _get_polygon(base_i.polygons, oave, rot_mtxT);
       std::vector<COORD> polygon_j = _get_polygon(base_j.polygons, oave, rot_mtxT);
       return pia_area(polygon_i, polygon_j);
}

std::vector<COORD> BasePairInfo::_get_polygon(const std::vector<COORD>& polygon_in, const COORD& oave, const COORD rot_mtx[3])
{
       std::vector<COORD> polygon_out;
       polygon_out.clear();

       COORD coord, coordT;
       for (std::vector<COORD>::const_iterator pos = polygon_in.begin(); pos != polygon_in.end(); ++pos) {
            coord = *pos;
            vector_subtract(coord, oave);
            multi_matrix_vector(coordT, rot_mtx, coord);
            polygon_out.push_back(coordT);
       }

       return polygon_out;
}

int BasePairInfo::_getBestPair(const std::vector<_BASEPAIR>&  base_pair, _BASEPAIR& bp_pair, const std::vector<int>& skip)
{
       int j = -1;
       double pair_score = XBIG;
       for (unsigned int i = 0; i < base_pair.size(); ++i) {
            if (skip[base_pair[i].j]) continue;
            if (base_pair[i].type < 0) continue;
            if (base_pair[i].dsum < pair_score) {
                 bp_pair = base_pair[i];
                 j = base_pair[i].j;
                 pair_score = base_pair[i].dsum;
            }
       }
       return j;
}
