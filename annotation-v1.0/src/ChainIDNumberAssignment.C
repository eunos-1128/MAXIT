/*
FILE:     ChainIDNumberAssignment.C
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

#include "CompositeIndex.h"
#include "ChainIDNumberAssignment.h"
#include "utillib.h"

ChainIDNumberAssignment::ChainIDNumberAssignment()
{
       _init();
}

ChainIDNumberAssignment::~ChainIDNumberAssignment()
{
       _init();
}

void ChainIDNumberAssignment::_init()
{
       _water_cutoff_distance = 5.805;
       _skip_flag = false;
       _multi_model_chain_flag = false;
       _cell = NULL;
       _links = NULL;
       _moving_waters = NULL;
       _molecules = NULL;
       _mol_nonpolymer_chains.clear();
       _mol_water_residues.clear();
       _nonpolymer_group_mapping.clear();
       _minumum_number = -9999;
       _non_polymer_numbering.clear();
       _unique_index.clear();
       _last_numbers.clear();
}

void ChainIDNumberAssignment::setWaterCutoffDistance(const double& cutoff_distance)
{
       _water_cutoff_distance = cutoff_distance;
}

void ChainIDNumberAssignment::setSkipFlag()
{
       _skip_flag = true;
}

void ChainIDNumberAssignment::setLink(std::list<_LINK>* links)
{
       _links = links;
}

void ChainIDNumberAssignment::setCrySymmetry(CrySymmetry* cell)
{
       _cell = cell;
}

void ChainIDNumberAssignment::setMovingWaters(std::list<_MOVING_ATOM>* moving_waters)
{
       _moving_waters = moving_waters;
}

void ChainIDNumberAssignment::setMolecule(std::vector<RCSB::Molecule*>* mols)
{
       _molecules = mols;
}

void ChainIDNumberAssignment::doAssignment()
{
       bool skip = _checking_unique_numbering();
       if (skip) return;

       _get_nonpolymer_group_mapping();

       _checking_multi_model_chain_flag();

       for (std::vector<RCSB::Molecule*>::iterator pos = _molecules->begin(); pos != _molecules->end(); ++pos) {
            _reposition_waters(*pos);
       }

       if (_multi_model_chain_flag) {
            _assign_multi_model_nonpolymer_chainIDs();
            _assign_multi_model_water_chainIDs();
       }

       for (std::vector<RCSB::Molecule*>::iterator pos = _molecules->begin(); pos != _molecules->end(); ++pos) {
            (*pos)->InternalOrder();
            _assign_numering(*pos);
            (*pos)->update_residue_indices();
       }
}

void ChainIDNumberAssignment::_get_nonpolymer_group_mapping()
{
       _nonpolymer_group_mapping.clear();

       if (!_molecules || _molecules->empty() || !_links || _links->empty()) return;

       unsigned int size = _nonpolymer_group_mapping.size();
       while (true) {
            for (std::list<_LINK>::const_iterator pos = _links->begin(); pos != _links->end(); ++pos) {
                 if ((!pos->SymOP_1.empty() && pos->SymOP_1 != "1_555") || (!pos->SymOP_2.empty() && pos->SymOP_2 != "1_555") ||
                     (pos->type != "covale")) continue;
                 RCSB::Residue* res1 = (*_molecules)[0]->find_pdb_residue(pos->fstAtom->pdb_chnid(),
                    pos->fstAtom->pdb_resnam(), pos->fstAtom->pdb_resnum(), pos->fstAtom->ins_code());
                 if (!res1 || (res1->chain_type() == "ATOMS")) continue;
                 RCSB::Residue* res2 = (*_molecules)[0]->find_pdb_residue(pos->sndAtom->pdb_chnid(),
                    pos->sndAtom->pdb_resnam(), pos->sndAtom->pdb_resnum(), pos->sndAtom->ins_code());
                 if (!res2 || (res2->chain_type() == "ATOMS")) continue;
                 if (((res1->chain_type() == "ATOMP") || (res1->chain_type() == "ATOMN")) &&
                     ((res2->chain_type() == "ATOMP") || (res2->chain_type() == "ATOMN"))) continue;

                 std::string idx1 = CompositeIndex::getIndex(res1->pdb_chnid(), res1->ResName(), res1->pdb_res_no(), res1->ins_code());
                 std::string idx2 = CompositeIndex::getIndex(res2->pdb_chnid(), res2->ResName(), res2->pdb_res_no(), res2->ins_code());
                 if ((res1->chain_type() == "ATOMP") || (res1->chain_type() == "ATOMN")) {
                      if (_nonpolymer_group_mapping.find(idx2) == _nonpolymer_group_mapping.end())
                           _nonpolymer_group_mapping.insert(std::make_pair(idx2, res1->pdb_chnid()));
                 } else if ((res2->chain_type() == "ATOMP") || (res2->chain_type() == "ATOMN")) {
                      if (_nonpolymer_group_mapping.find(idx1) == _nonpolymer_group_mapping.end())
                           _nonpolymer_group_mapping.insert(std::make_pair(idx1, res2->pdb_chnid()));
                 } else {
                      std::map<std::string, std::string>::const_iterator mpos1 = _nonpolymer_group_mapping.find(idx1);
                      std::map<std::string, std::string>::const_iterator mpos2 = _nonpolymer_group_mapping.find(idx2);
                      if (mpos1 == _nonpolymer_group_mapping.end() && mpos2 != _nonpolymer_group_mapping.end())
                           _nonpolymer_group_mapping.insert(std::make_pair(idx1, mpos2->second));
                      else if (mpos1 != _nonpolymer_group_mapping.end() && mpos2 == _nonpolymer_group_mapping.end())
                           _nonpolymer_group_mapping.insert(std::make_pair(idx2, mpos1->second));
                 }
            }
            if (_nonpolymer_group_mapping.size() == size) break;
            size = _nonpolymer_group_mapping.size();
       }
}

void ChainIDNumberAssignment::_checking_multi_model_chain_flag()
{
       if (_molecules->size() < 2) return;

       std::set<std::string> polymer_chain_id_set_1, polymer_chain_id_set_2;
       std::map<std::string, std::vector<RCSB::Chain*> > polymer_chains, other_chains;

       std::vector<RCSB::Molecule*>::const_iterator pos = _molecules->begin();
       (*pos)->GetChains(polymer_chain_id_set_1, polymer_chains, other_chains);
       if (polymer_chain_id_set_1.size() > 1) {
            ++pos;
            bool has_multi_chain_flag = true;
            for (; pos != _molecules->end(); ++pos) {
                 (*pos)->GetChains(polymer_chain_id_set_2, polymer_chains, other_chains);
                 if (polymer_chain_id_set_2 != polymer_chain_id_set_1) {
                      has_multi_chain_flag = false;
                      break;
                 }
            }
            if (has_multi_chain_flag) _multi_model_chain_flag = true;
       }
}

void ChainIDNumberAssignment::_reposition_waters(RCSB::Molecule* mol)
{
       std::vector<RCSB::Chain*> polymer_chains, carbohydrate_chains, nonpolymer_chains;
       std::list<RCSB::Residue*> water_lists;
       mol->GetChainsAndWaters(polymer_chains, carbohydrate_chains, nonpolymer_chains, water_lists);
       if (!carbohydrate_chains.empty()) {
            std::vector<RCSB::Residue*> residues;
            std::set<std::string> used_chain_ids;
            used_chain_ids.clear();
            for (std::vector<RCSB::Chain*>::const_iterator pos = polymer_chains.begin(); pos != polymer_chains.end(); ++pos) {
                 used_chain_ids.insert((*pos)->PDB_ChainID());
            }
            for (std::vector<RCSB::Chain*>::const_iterator pos = carbohydrate_chains.begin(); pos != carbohydrate_chains.end(); ++pos) {
                 if (used_chain_ids.find((*pos)->PDB_ChainID()) == used_chain_ids.end()) {
                      used_chain_ids.insert((*pos)->PDB_ChainID());
                      continue;
                 }
                 std::string pdb_chain_id = get_next_available_pdb_chain_id(used_chain_ids);
                 (*pos)->set_PDB_ChainID(pdb_chain_id);

                 (*pos)->GetFirstResidueList(residues);
                 while (!residues.empty()) {
                      for (std::vector<RCSB::Residue*>::const_iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                           (*rpos)->set_pdb_chnid(pdb_chain_id);
                           RCSB::Atom* atom = (*rpos)->GetFirstAtom();
                           while (atom) {
                                atom->set_pdb_chnid(pdb_chain_id);
                                atom = (*rpos)->GetNextAtom();
                           }
                      }
                      (*pos)->GetNextResidueList(residues);
                 }
            }
       }
       if (!nonpolymer_chains.empty() || !water_lists.empty()) {
            xtal mycrys;
            bool found = false;
            if ((polymer_chains.size() > 1) || !water_lists.empty()) {
                 std::set<std::string> chain_types;
                 chain_types.clear();
                 std::set<std::string> res_names;
                 res_names.clear();
                 found = set_mycrys(mycrys, *_cell, mol, chain_types, res_names, true, true, true);
            }
            if (polymer_chains.size() == 1)
                 _assign_single_chainID(polymer_chains[0]->PDB_ChainID(), mol);
            else if ((polymer_chains.size() > 1) && found)
                 _assign_multiple_chainIDs(mycrys, nonpolymer_chains);
            __reposition_waters(mycrys, mol, water_lists);
       }
}

void ChainIDNumberAssignment::_assign_multi_model_nonpolymer_chainIDs()
{
       if (_mol_nonpolymer_chains.empty()) return;

       // key: residue_id - chnid_resname_resnum_ins_code
       // value.key: assigned pdb_chainid
       // value.value: count
       std::map<std::string, std::map<std::string, int> > count_mapping;
       count_mapping.clear();
       std::map<std::string, int> t_map;
       std::vector<RCSB::Residue*> residues;

       for (std::vector<std::vector<std::pair<std::string, RCSB::Chain*> > >::const_iterator vcpos = _mol_nonpolymer_chains.begin();
            vcpos != _mol_nonpolymer_chains.end(); ++vcpos) {
            for (std::vector<std::pair<std::string, RCSB::Chain*> >::const_iterator cpos = vcpos->begin(); cpos != vcpos->end(); ++cpos) {
                 if (cpos->first.empty()) continue;

                 cpos->second->GetFirstResidueList(residues);
                 while (!residues.empty()) {
                      for (std::vector<RCSB::Residue*>::const_iterator vpos = residues.begin(); vpos != residues.end(); ++vpos) {
                           std::string idx = CompositeIndex::getIndex((*vpos)->pdb_chnid(), (*vpos)->ResName(), (*vpos)->pdb_res_no(), (*vpos)->ins_code());
                           std::map<std::string, std::map<std::string, int> >::iterator mpos = count_mapping.find(idx);
                           if (mpos == count_mapping.end()) {
                                t_map.clear();
                                t_map.insert(std::make_pair(cpos->first, 1));
                                count_mapping.insert(std::make_pair(idx, t_map));
                           } else {
                                std::map<std::string, int>::iterator mmpos = mpos->second.find(cpos->first);
                                if (mmpos != mpos->second.end())
                                     mmpos->second++;
                                else mpos->second.insert(std::make_pair(cpos->first, 1));
                           } 
                      }
                      cpos->second->GetNextResidueList(residues);
                 }
            }
       }

       // key: residue_id - chnid_resname_resnum_ins_code
       // value: assigned pdb_chainid
       std::map<std::string, std::string> chain_id_mapping;
       chain_id_mapping.clear();
       for (std::map<std::string, std::map<std::string, int> >::const_iterator mpos = count_mapping.begin(); mpos != count_mapping.end(); ++mpos) {
            std::string pdb_chainid = "";
            int count = 0;
            for (std::map<std::string, int>::const_iterator mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                 if (mmpos->second > count) {
                      count = mmpos->second;
                      pdb_chainid = mmpos->first;
                 }
            }
            if (!pdb_chainid.empty()) chain_id_mapping.insert(std::make_pair(mpos->first, pdb_chainid));
       }

       // update chain IDs for nonpolymer
       for (std::vector<std::vector<std::pair<std::string, RCSB::Chain*> > >::iterator vcpos = _mol_nonpolymer_chains.begin();
            vcpos != _mol_nonpolymer_chains.end(); ++vcpos) {
            for (std::vector<std::pair<std::string, RCSB::Chain*> >::const_iterator cpos = vcpos->begin(); cpos != vcpos->end(); ++cpos) {
                 std::string pdb_chainid = "";
                 cpos->second->GetFirstResidueList(residues);
                 std::string idx = CompositeIndex::getIndex(residues[0]->pdb_chnid(), residues[0]->ResName(), residues[0]->pdb_res_no(),
                                                            residues[0]->ins_code());
                 std::map<std::string, std::string>::const_iterator mpos = chain_id_mapping.find(idx);
                 if (mpos != chain_id_mapping.end()) pdb_chainid = mpos->second;
                 else if (!cpos->first.empty()) pdb_chainid = cpos->first;
                 if (!pdb_chainid.empty()) cpos->second->set_PDB_ChainID_to_residues(pdb_chainid);
            }
       }
}

void ChainIDNumberAssignment::_assign_multi_model_water_chainIDs()
{
       if (_mol_water_residues.empty()) return;


       // key: residue_id - chnid_resname_resnum_ins_code
       // value.key: assigned pdb_chainid
       // value.value: dist set
       std::map<std::string, std::map<std::string, std::multiset<double> > > count_mapping;
       count_mapping.clear();
       std::map<std::string, std::multiset<double> > t_map;
       std::multiset<double> t_set;

       for (std::vector<std::pair<RCSB::Molecule*, std::vector<_WATER_RESIDUE> > >::iterator vvpos = _mol_water_residues.begin();
            vvpos != _mol_water_residues.end(); ++vvpos) {
            for (std::vector<_WATER_RESIDUE>::iterator vpos = vvpos->second.begin(); vpos != vvpos->second.end(); ++vpos) {
                 if (vpos->ref_atom) vpos->pdb_chnid = vpos->ref_atom->pdb_chnid();
                 std::string idx = CompositeIndex::getIndex(vpos->residue->pdb_chnid(), vpos->residue->ResName(), vpos->residue->pdb_res_no(),
                                                            vpos->residue->ins_code());
                 std::map<std::string, std::map<std::string, std::multiset<double> > >::iterator mpos = count_mapping.find(idx);
                 if (mpos == count_mapping.end()) {
                      t_set.clear();
                      t_set.insert(vpos->dist);
                      t_map.clear();
                      t_map.insert(std::make_pair(vpos->pdb_chnid, t_set));
                      count_mapping.insert(std::make_pair(idx, t_map));
                 } else {
                      std::map<std::string, std::multiset<double> >::iterator mmpos = mpos->second.find(vpos->pdb_chnid);
                      if (mmpos != mpos->second.end())
                           mmpos->second.insert(vpos->dist);
                      else {
                           t_set.clear();
                           t_set.insert(vpos->dist);
                           mpos->second.insert(std::make_pair(vpos->pdb_chnid, t_set));
                      }
                 }
            }
       }

       // key: residue_id - chnid_resname_resnum_ins_code
       // value.first: assigned pdb_chainid
       // value.second: dist
       std::map<std::string, std::pair<std::string, double> > chain_id_mapping;
       chain_id_mapping.clear();

       for (std::map<std::string, std::map<std::string, std::multiset<double> > >::const_iterator
            mpos = count_mapping.begin(); mpos != count_mapping.end(); ++mpos) {
            std::string pdb_chainid = "";
            t_set.clear();
            for (std::map<std::string, std::multiset<double> >::const_iterator mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                 if (mmpos->second.size() > t_set.size()) {
                      t_set = mmpos->second;
                      pdb_chainid = mmpos->first;
                 }
            }
            if (!pdb_chainid.empty() && !t_set.empty()) chain_id_mapping.insert(std::make_pair(mpos->first, std::make_pair(pdb_chainid, *(t_set.begin()))));
       }

       // update chain IDs for water
       for (std::vector<std::pair<RCSB::Molecule*, std::vector<_WATER_RESIDUE> > >::iterator vvpos = _mol_water_residues.begin();
            vvpos != _mol_water_residues.end(); ++vvpos) {
            for (std::vector<_WATER_RESIDUE>::iterator vpos = vvpos->second.begin(); vpos != vvpos->second.end(); ++vpos) {
                 std::string idx = CompositeIndex::getIndex(vpos->residue->pdb_chnid(), vpos->residue->ResName(), vpos->residue->pdb_res_no(),
                                                            vpos->residue->ins_code());
                 std::map<std::string, std::pair<std::string, double> >::const_iterator mpos = chain_id_mapping.find(idx);
                 if (mpos != chain_id_mapping.end()) {
                      vpos->pdb_chnid = mpos->second.first;
                      vpos->dist = mpos->second.second;
                 }
            }
            __reposition_waters(vvpos->first, vvpos->second);
       }
}
