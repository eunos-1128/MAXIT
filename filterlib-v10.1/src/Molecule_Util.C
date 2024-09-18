/*
FILE:     Molecule_Util.C
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
#include "Molecule.h"
#include "utillib.h"

using namespace RCSB;

void Molecule::GetChainsAndWaters(std::vector<RCSB::Chain*>& polymer_chains, std::vector<RCSB::Chain*>& carbohydrate_chains,
                                  std::vector<RCSB::Chain*>& nonpolymer_chains, std::list<RCSB::Residue*>& water_lists)
{
       polymer_chains.clear();
       carbohydrate_chains.clear();
       nonpolymer_chains.clear();
       water_lists.clear();

       std::vector<RCSB::Residue*> residues;
       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN")
                 polymer_chains.push_back(chain);
            else if (chain->chain_type() == "ATOMS")
                 carbohydrate_chains.push_back(chain);
            else if (chain->GetFirstResidue()->ResName() == "HOH" || chain->GetFirstResidue()->ResName() == "DOD") {
                 chain->GetFirstResidueList(residues);
                 while (!residues.empty()) {
                      for (std::vector<RCSB::Residue*>::const_iterator vpos = residues.begin(); vpos != residues.end(); ++vpos) {
                           water_lists.push_back(*vpos);
                      }
                      chain->GetNextResidueList(residues);
                 }
            } else nonpolymer_chains.push_back(chain);

            chain = GetNextChain();
       }
}

void Molecule::GetChains(const std::string& pdb_chainid, std::vector<RCSB::Chain*>& polymer_chains, std::vector<RCSB::Chain*>& nonpolymer_chains,
                         std::vector<RCSB::Chain*>& waters)
{
       polymer_chains.clear();
       nonpolymer_chains.clear();
       waters.clear();
       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if (chain->PDB_ChainID() != pdb_chainid) {
                 chain = GetNextChain();
                 continue;
            }

            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMS")
                 polymer_chains.push_back(chain);
            else if (chain->GetFirstResidue()->ResName() == "HOH" ||
                     chain->GetFirstResidue()->ResName() == "DOD")
                 waters.push_back(chain);
            else nonpolymer_chains.push_back(chain);

            chain = GetNextChain();
       }
}

void Molecule::GetChains(std::set<std::string>& pdb_chainids, std::map<std::string, std::vector<RCSB::Chain*> >& polymer_chains,
                         std::map<std::string, std::vector<RCSB::Chain*> >& other_chains)
{
       pdb_chainids.clear();
       polymer_chains.clear();
       other_chains.clear();

       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            pdb_chainids.insert(chain->PDB_ChainID());

            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMNS")
                 _insert_chains(polymer_chains, chain);
            else _insert_chains(other_chains, chain);

            chain = GetNextChain();
       }
}

void Molecule::GetChains(std::vector<std::string>& pdb_chain_ids, std::map<std::string, std::vector<RCSB::Chain*> >& polymer_chains,
                         std::map<std::string, std::vector<RCSB::Chain*> >& non_polymer_chains, std::map<std::string, std::vector<RCSB::Chain*> >& water_chains)
{
       pdb_chain_ids.clear();
       polymer_chains.clear();
       non_polymer_chains.clear();
       water_chains.clear();

       std::set<std::string> chain_id_set;
       chain_id_set.clear();

       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            if (chain_id_set.find(chain->PDB_ChainID()) == chain_id_set.end()) {
                 chain_id_set.insert(chain->PDB_ChainID());
                 pdb_chain_ids.push_back(chain->PDB_ChainID());
            }

            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMS")
                 _insert_chains(polymer_chains, chain);
            else if (chain->GetFirstResidue()->ResName() == "HOH" || chain->GetFirstResidue()->ResName() == "DOD")
                 _insert_chains(water_chains, chain);
            else _insert_chains(non_polymer_chains, chain);

            chain = GetNextChain();
       }
}

void Molecule::GetMolInfo(std::vector<std::string>& pdb_chain_ids, std::map<std::string, RCSB::Chain*>& polymers, std::map<std::string,
                          std::list<RCSB::Residue*> >& nonpolymers, std::map<std::string, std::list<RCSB::Residue*> >& waters)
{
       pdb_chain_ids.clear();
       polymers.clear();
       nonpolymers.clear();
       waters.clear();

       std::set<std::string> chain_id_set;
       chain_id_set.clear();

       std::list<RCSB::Residue*> r_list;

       std::vector<RCSB::Residue*> residues;
       RCSB::Chain* chain = GetFirstChain();
       while (chain) {
            std::string chainID = chain->PDB_ChainID();
            if (chain_id_set.find(chainID) == chain_id_set.end()) {
                 chain_id_set.insert(chainID);
                 pdb_chain_ids.push_back(chainID);
            }

            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN" /* || chain->chain_type() == "ATOMS" */ )
                 polymers.insert(std::make_pair(chainID, chain));
            else if (chain->GetFirstResidue()->ResName() == "HOH" || chain->GetFirstResidue()->ResName() == "DOD") {
                 chain->GetFirstResidueList(residues);
                 while (!residues.empty()) {
                      for (std::vector<RCSB::Residue*>::const_iterator vpos = residues.begin(); vpos != residues.end(); ++vpos) {
                           std::map<std::string, std::list<RCSB::Residue*> >::iterator mpos = waters.find(chainID);
                           if (mpos != waters.end()) mpos->second.push_back(*vpos);
                           else {
                                r_list.clear();
                                r_list.push_back(*vpos);
                                waters.insert(std::make_pair(chainID, r_list));
                           }
                      }
                      chain->GetNextResidueList(residues);
                 }
            } else {
                 chain->GetFirstResidueList(residues);
                 while (!residues.empty()) {
                      for (std::vector<RCSB::Residue*>::const_iterator vpos = residues.begin(); vpos != residues.end(); ++vpos) {
                           std::map<std::string, std::list<RCSB::Residue*> >::iterator mpos = nonpolymers.find(chainID);
                           if (mpos != nonpolymers.end()) mpos->second.push_back(*vpos);
                           else {
                                r_list.clear();
                                r_list.push_back(*vpos);
                                nonpolymers.insert(std::make_pair(chainID, r_list));
                           }
                      }
                      chain->GetNextResidueList(residues);
                 }
            }

            chain = GetNextChain();
       }
}

bool Molecule::PolymerUpdate(const std::vector<std::vector<std::string> >& mapping, const std::vector<std::vector<std::string> >& in_links,
                             std::vector<std::pair<std::string, std::vector<std::vector<RCSB::Atom*> > > >& out_links)
{
       bool _successful = true;

       if (mapping.empty()) return _successful;

       // map key: residue ID (pdbchnid_name_number_inscode)
       // value: residue index
       std::multimap<std::string, int> residue_mapping;
       residue_mapping.clear();

       std::vector<int> res_idx;
       res_idx.clear();

       std::set<int> involved_chain_set;
       involved_chain_set.clear();

       for (std::vector<std::vector<std::string> >::const_iterator vpos = mapping.begin(); vpos != mapping.end(); ++vpos) {
            std::string cs = (*vpos)[1] + "_" + (*vpos)[2] + "_" + (*vpos)[3] + "_" + (*vpos)[4];
            if (residue_mapping.find(cs) != residue_mapping.end()) continue;
            std::multimap<std::string, int>::iterator mpos = _pdbIndex.find(cs);
            if (mpos == _pdbIndex.end()) {
                 std::string error = "Residue (" + (*vpos)[1] + " " + (*vpos)[2] + " " + (*vpos)[3] + (*vpos)[4] + ") does not eixt in Model "
                                   + String::IntToString(_Mol_ID) + ".\n";
                 _logIo->message(error.c_str());
                 _successful = false;
            }
            residue_mapping.insert(std::make_pair(cs, mpos->second));
            res_idx.push_back(mpos->second);
            involved_chain_set.insert(_residues[mpos->second]->chain_index());
       }
       if (!_successful || res_idx.empty()) return _successful;

       int chn_idx = _residues[res_idx[0]]->chain_index();

       std::multimap<std::string, RCSB::Atom*> found_atoms;
       found_atoms.clear();
       for (std::vector<int>::const_iterator rpos = res_idx.begin(); rpos != res_idx.end(); ++rpos) {
            _get_involved_atoms(_residues[*rpos], found_atoms);
       }
       if (found_atoms.empty()) return _successful;

       std::list<RCSB::Residue*> residue_list;
       std::map<std::string, std::string> polymer_type_mapping;
       _successful = _get_new_residue_list(mapping, found_atoms, in_links, residue_list, polymer_type_mapping, out_links);
       if (!_successful) return _successful;

       std::vector<RCSB::Residue*> polymer_residues, non_polymer_residues;
       polymer_residues.clear();
       non_polymer_residues.clear();
       for (std::list<RCSB::Residue*>::const_iterator lpos = residue_list.begin(); lpos != residue_list.end(); ++lpos) {
            std::string idx = (*lpos)->pdb_chnid() + "_" + (*lpos)->ResName() + "_" + (*lpos)->pdb_res_no();
            std::map<std::string, std::string>::const_iterator mpos = polymer_type_mapping.find(idx);
            if (mpos != polymer_type_mapping.end() && mpos->second == "polymer")
                 polymer_residues.push_back(*lpos);
            else non_polymer_residues.push_back(*lpos);
       }
       if (polymer_residues.empty()) {
            _logIo->message("No polymer residue defined in chopper output file.\n");
            return false;
       }

       std::vector<RCSB::Chain*> polymer_chains, nonpolymer_chains, waters;
       GetChains(mapping[0][1], polymer_chains, nonpolymer_chains, waters);

       for (std::vector<RCSB::Chain*>::const_iterator pos = polymer_chains.begin(); pos != polymer_chains.end(); ++pos) {
            if ((*pos)->index() == chn_idx) {
                 _successful = (*pos)->Update(mapping[0][2], mapping[0][3], mapping[0][4], polymer_residues);
                 involved_chain_set.erase(chn_idx);
            }
       }
       if (!_successful) return _successful;

       for (std::vector<RCSB::Residue*>::const_iterator pos = polymer_residues.begin(); pos != polymer_residues.end(); ++pos) {
            insert_a_residue(*pos);
       }

       _remove_old_residues_and_chains(residue_mapping, involved_chain_set);
       if (!involved_chain_set.empty()) {
            std::vector<RCSB::Chain*> tmp_chain;
            tmp_chain.clear();
            for (std::vector<RCSB::Chain*>::const_iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
                 RCSB::Chain* chain = *pos;
                 if (chain == NULL) continue;
                 tmp_chain.push_back(chain);
            }
            _chains = tmp_chain;
       }

       if (!non_polymer_residues.empty()) {
            int last_number = 0;
            _get_last_number(last_number, polymer_chains);
            _get_last_number(last_number, nonpolymer_chains);

            for (std::vector<RCSB::Residue*>::const_iterator pos = non_polymer_residues.begin(); pos != non_polymer_residues.end(); ++pos) {
                 last_number++;
                 RCSB::Atom* atom = (*pos)->GetFirstAtom();
                 while (atom) {
                      atom->set_pdb_resnum(String::IntToString(last_number));
                      atom->set_pdb_chnid(mapping[0][1]);
                      atom = (*pos)->GetNextAtom();
                 }

                 (*pos)->set_pdb_res_no(String::IntToString(last_number));
                 (*pos)->set_pdb_chnid(mapping[0][1]);
                 if ((*pos)->chain_type() == "ATOMP" || (*pos)->chain_type() == "ATOMN") (*pos)->set_chain_type("HETAIN");
                 (*pos)->set_alt_loc();
                 (*pos)->get_atom_type();
                 (*pos)->reorder_atoms();
                 insert_a_residue(*pos);
                 RCSB::Chain* chain = insert_a_chain(*pos, true);
                 chain->insert_a_residue(*pos);
            }
       }

       if (!involved_chain_set.empty() || !non_polymer_residues.empty()) InternalOrder();

       return _successful;
}

bool Molecule::PrdUpdate(const std::vector<std::vector<std::string> >& mapping, const std::map<std::string, SEQ>& seqs,
                         const std::vector<std::vector<std::string> >& in_links, std::vector<std::pair<std::string,
                         std::vector<std::vector<RCSB::Atom*> > > >& out_links)
{
       bool _successful = true;

       if (mapping.empty()) return _successful;

       bool branched_flag = true;
       if (!seqs.empty()) branched_flag = false;

       std::set<std::string> old_unique_set, new_unique_set;
       old_unique_set.clear();
       new_unique_set.clear();

       for (std::vector<std::vector<std::string> >::const_iterator vpos = mapping.begin(); vpos != mapping.end(); ++vpos) {
            if ((*vpos)[11] != "branched") branched_flag = false;
            std::string cs = CompositeIndex::getIndex((*vpos)[1], (*vpos)[2], (*vpos)[3], (*vpos)[4]);
            old_unique_set.insert(cs);
            new_unique_set.insert((*vpos)[6] + "_" + (*vpos)[7]);
       }

       // integer is chain index
       std::set<int> prd_involved_chain_set;
       _get_PRD_involved_chains(mapping, prd_involved_chain_set, _successful);

       if (!_successful) return _successful;
       if (prd_involved_chain_set.empty()) return _successful;

       if (branched_flag && (old_unique_set.size() == 1) && (prd_involved_chain_set.size() == 1) && (new_unique_set.size() == 1)) {
            RCSB::Residue* res = this->find_pdb_residue(mapping[0][1], mapping[0][2], mapping[0][3], mapping[0][4]);
            if (res) {
                 std::map<std::string, std::string> atom_mapping;
                 atom_mapping.clear();
                 for (std::vector<std::vector<std::string> >::const_iterator vpos = mapping.begin(); vpos != mapping.end(); ++vpos) {
                      atom_mapping.insert(std::make_pair((*vpos)[5], (*vpos)[8]));
                 }
                 res->update_atom_mapping(_Mol_ID, mapping[0][6], atom_mapping);
                 bool is_removed = false;
                 RCSB::Chain* chain = GetIndexChain((*(prd_involved_chain_set.begin())), is_removed);
                 if (chain) {
                      chain->Update_SeqRes(mapping[0][2], mapping[0][3], mapping[0][4], mapping[0][6]);
                      chain->set_chain_type("HETAIN");
                 }
                 return _successful;
            }
       }

       // map key: residue ID (pdbchnid_name_number_inscode)
       // value: residue index
       std::multimap<std::string, int> residue_mapping;

       // map key: atom ID (pdbchnid_name_number_inscode_atomname)
       // value: atom pointer
       std::multimap<std::string, RCSB::Atom*> found_atoms;

       _get_PRD_involved_residue_and_atoms(prd_involved_chain_set, residue_mapping, found_atoms, _successful);

       if (!_successful) return _successful;
       if (found_atoms.empty()) return _successful;

       std::list<RCSB::Residue*> residue_list, t_residue_list;
       std::map<std::string, std::string> polymer_type_mapping;
       _successful = _get_new_residue_list(mapping, found_atoms, in_links, residue_list, polymer_type_mapping, out_links);
       if (!_successful) return _successful;

       _remove_old_residues_and_chains(residue_mapping, prd_involved_chain_set);
       if (branched_flag) {
            if ((old_unique_set.size() == 1) && (prd_involved_chain_set.size() == 1) && in_links.empty()) {
                 char ins_code = '@';
                 std::string ins_code_str;
                 bool first = true;
                 for (std::list<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                      if (first) {
                           ins_code_str = (*rpos)->ins_code();
                           if (!ins_code_str.empty()) ins_code = ins_code_str[0];
                      } else {
                           ins_code++;
                           ins_code_str.clear();
                           ins_code_str += ins_code;
                           (*rpos)->set_ins_code(ins_code_str);
                           RCSB::Atom* atom = (*rpos)->GetFirstAtom();
                           while (atom) {
                                atom->set_ins_code(ins_code_str);
                                atom = (*rpos)->GetNextAtom();
                           }
                      }

                      t_residue_list.clear();
                      t_residue_list.push_back(*rpos);
                      _add_new_residues_and_chains(t_residue_list, seqs, "HETAIN");
                      first = false;
                 }
            } else _add_new_residues_and_chains(residue_list, seqs, "ATOMS");
       } else _add_new_residues_and_chains(residue_list, seqs, "");

       std::vector<std::map<std::string, std::vector<std::string> > > pro_atom_name_mapping_list, other_atom_name_mapping_list;
       get_proline_n_terminal_hydrogen_mapping_list(pro_atom_name_mapping_list);
       get_other_n_terminal_hydrogen_mapping_list(other_atom_name_mapping_list);

       std::vector<RCSB::Residue*> residues;
       std::vector<RCSB::Chain*> tmp_chain;
       tmp_chain.clear();
       for (std::vector<RCSB::Chain*>::const_iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            RCSB::Chain* chain = *pos;
            if (chain == NULL) continue;
            if (chain->chain_type() == "ATOMP") {
                 for (unsigned int i = 0; i < chain->SeqLen(); ++i) {
                      chain->GetResidueList(i, residues);
                      if (residues.empty()) continue;
                      for (std::vector<RCSB::Residue*>::iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                           if ((*rpos)->ResName() == "PRO")
                                (*rpos)->fix_n_terminal_hydrogen((i == 0), pro_atom_name_mapping_list);
                           else (*rpos)->fix_n_terminal_hydrogen((i == 0), other_atom_name_mapping_list);
                      }
                      break;
                 }
            }
            tmp_chain.push_back(chain);
       }
       _chains = tmp_chain;

       InternalOrder();

       return _successful;
}

bool Molecule::PrdUpdateAtom(const std::map<std::string, std::map<std::string, std::string> >& atom_mapping)
{
       bool _successful = true;
       if (atom_mapping.empty()) return _successful;

       std::vector<std::string> data;
       for (std::map<std::string, std::map<std::string,std::string> >::const_iterator mpos = atom_mapping.begin(); mpos != atom_mapping.end(); ++mpos) {
            std::multimap<std::string, int>::iterator pos = _pdbIndex.find(mpos->first);
            if (pos == _pdbIndex.end()) {
                 get_wordarray(data, mpos->first, "_");
                 if (data.size() == 3) data.push_back("");
                 std::string error = "Residue (" + data[0] + " " + data[1] + " " + data[2] + data[3] + ") does not eixt in Model "
                                   + String::IntToString(_Mol_ID) + ".\n";
                 _logIo->message(error.c_str());
                 _successful = false;
                 continue;
            }
            if (!_residues[pos->second]->update_atom_name(mpos->second, _Mol_ID, true))
                 _successful = false;
       }

       return _successful;
}

bool Molecule::PrdUpdateNumbering(const std::map<std::string, std::string>& num_mapping)
{
       bool _successful = true;
       if (num_mapping.empty()) return _successful;

       std::set<int> chain_index_set;
       chain_index_set.clear();

       std::vector<std::string> data;
       for (std::map<std::string, std::string>::const_iterator mpos = num_mapping.begin(); mpos != num_mapping.end(); ++mpos) {
            std::multimap<std::string, int>::iterator pos = _pdbIndex.find(mpos->first);
            if (pos == _pdbIndex.end()) {
                 get_wordarray(data, mpos->first, "_");
                 if (data.size() == 3) data.push_back("");
                 std::string error = "Residue (" + data[0] + " " + data[1] + " " + data[2] + data[3] + ") does not eixt in Model "
                                   + String::IntToString(_Mol_ID) + ".\n";
                 _logIo->message(error.c_str());
                 _successful = false;
                 continue;
            }

            int index = pos->second;
            _residues[index]->set_pdb_res_no(mpos->second);
            chain_index_set.insert(_residues[index]->chain_index());
            RCSB::Atom* atom = _residues[index]->GetFirstAtom();
            while (atom) {
                 atom->set_pdb_resnum(mpos->second);
                 atom = _residues[index]->GetNextAtom();
            }

            _pdbIndex.erase(mpos->first);
            std::string cs = CompositeIndex::getIndex(_residues[index]->pdb_chnid(), _residues[index]->ResName(),
                                     _residues[index]->pdb_res_no(), _residues[index]->ins_code());
            _pdbIndex.insert(std::make_pair(cs, index));
       }

       if (_successful && !chain_index_set.empty()) {
            for (std::set<int>::const_iterator spos = chain_index_set.begin(); spos != chain_index_set.end(); ++spos) {
                 bool is_removed = false;
                 RCSB::Chain* chain = GetIndexChain(*spos, is_removed);
                 if (!chain) continue;
                 chain->update_number();
            }
       }

       return _successful;
}

bool Molecule::MergePolymer(const int& new_entity_id, const std::string& new_asym_id, const std::string& new_chain_id,
                            const std::vector<std::vector<std::string> >& group, std::vector<std::string>& seqs,
                            std::vector<std::vector<int> >& merged_entityids)
{
       seqs.clear();
       merged_entityids.clear();

       bool _successful = true;
       if (group.empty()) return _successful;

       RCSB::Chain* chain = new Chain;
       chain->setCCDic(_ccDic);
       chain->setLog(_logIo);
       chain->setMessage(_messageIo);

       chain->set_ChainID(new_asym_id);
       chain->set_PDB_ChainID(new_chain_id);
       chain->set_PUB_ChainID(new_chain_id);
       chain->set_prev_entity_id(String::IntToString(new_entity_id));
       chain->set_index(_chain_index);
       chain->set_order(_chain_index);
       _chain_index++;

       std::vector<int> data;
       std::set<int> removed_chain_set, removed_residue_chain_set;
       std::set<int>  removed_residue_set;
       removed_chain_set.clear();
       removed_residue_chain_set.clear();
       removed_residue_set.clear();

       std::map<std::string, int> type_count;
       type_count.clear();

       int residue_count = 0;
       std::string residue_type_info = "";
       std::string chain_type = "";
       for (std::vector<std::vector<std::string> >::const_iterator pos = group.begin(); pos != group.end(); ++pos) {
            if (pos->size() == 1) {
                 RCSB::Chain* chn = GetPolyChain((*pos)[0]);
                 if (!chn) {
                      std::string error = "Chain " + (*pos)[0] + " does not exist in Model " + String::IntToString(_Mol_ID) + ".\n";
                      _logIo->message(error.c_str());
                      _successful = false;
                      continue;
                 }

                 if (chn->chain_type() == "ATOMP")
                      chain_type = chn->chain_type();
                 else if (chn->chain_type() == "ATOMN" && chain_type != "ATOMP")
                      chain_type = chn->chain_type();
                 else if (chn->chain_type() == "ATOMS" && chain_type != "ATOMP" && chain_type != "ATOMN")
                      chain_type = chn->chain_type();

                 data.clear();
                 data.push_back(chn->int_prev_entity_id());
                 data.push_back(1);
                 data.push_back(chn->SeqLen());
                 merged_entityids.push_back(data);
                 removed_chain_set.insert(chn->index());
                 chain->merge_chain(chn);
            } else if (pos->size() == 4) {
                 residue_count++;
                 RCSB::Residue* res = find_pdb_residue((*pos)[0], (*pos)[1], (*pos)[2], (*pos)[3]);
                 if (!res) {
                      std::string error = "Residue (" + (*pos)[0] + " " + (*pos)[1] + " " + (*pos)[2] + (*pos)[3] + ") does not eixt in Model "
                                        + String::IntToString(_Mol_ID) + ".\n";
                      _logIo->message(error.c_str());
                      _successful = false;
                      continue;
                 }

                 add_type(type_count, _ccDic->find_residue_type(res->ResName()));

                 std::string monomer_type = _ccDic->find_cc_metadata(res->ResName(), "type");
                 if (monomer_type.empty()) monomer_type = "unknown";
                 if (!residue_type_info.empty()) residue_type_info += "\n";
                 residue_type_info += "Residue (" + (*pos)[0] + " " + (*pos)[1] + " " + (*pos)[2] + (*pos)[3] + ") has '" + monomer_type + "' type.";

                 removed_residue_chain_set.insert(res->chain_index());
                 removed_residue_set.insert(res->index());
                 chain->merge_residue(res);
            }
       }

       if (!_successful) {
            delete chain;
            return _successful;
       }

       if (chain_type.empty()) {
            int atomp_type_count = 0;
            std::map<std::string, int>::const_iterator mpos = type_count.find("ATOMP"); if (mpos != type_count.end()) atomp_type_count = mpos->second;
            int atomn_type_count = 0; mpos = type_count.find("ATOMN"); if (mpos != type_count.end()) atomn_type_count = mpos->second;
            int atoms_type_count = 0; mpos = type_count.find("ATOMS"); if (mpos != type_count.end()) atoms_type_count = mpos->second;

            if (atoms_type_count == residue_count) chain_type = "ATOMS";
            else if (atomp_type_count > atomn_type_count) {
                 if ((atomp_type_count * 3) >= residue_count) chain_type = "ATOMP";
            } else if (atomn_type_count > 0) {
                 if ((atomn_type_count * 3) >= residue_count) chain_type = "ATOMN";
            }

            if (chain_type.empty()) {
                 std::string error = "Can not define the group's polymer type based on the following residue type:\n" + residue_type_info + "\n";
                 _logIo->message(error.c_str());
                 _successful = false;
                 return _successful;
            }
       }

       if (merged_entityids.size() == 1) {
            merged_entityids[0][2] = chain->SeqLen();
       } else if (merged_entityids.size() > 1) {
            for (unsigned int i = 1; i < merged_entityids.size(); ++i) {
                 merged_entityids[i][1] += merged_entityids[i - 1][2];
                 merged_entityids[i][2] += merged_entityids[i - 1][2];
            }
       }

       for (std::set<int>::const_iterator spos = removed_residue_chain_set.begin(); spos != removed_residue_chain_set.end(); ++spos) {
            bool is_removed = false;
            RCSB::Chain* chn = GetIndexChain(*spos, is_removed);
            if (!chn) continue;
            chn->remove_residue(removed_residue_set);
            if (chn->SeqLen() == 0) removed_chain_set.insert(*spos);
       }

       chain->set_chain_type(chain_type);
       chain->renumbering_polymer_chain();
       chain->update_indices();

       bool first = true;
       std::vector<RCSB::Chain*> tmp_chains;
       tmp_chains.clear();
       for (std::vector<RCSB::Chain*>::const_iterator pos = _chains.begin(); pos != _chains.end(); ++pos) { 
            if (removed_chain_set.find((*pos)->index()) == removed_chain_set.end()) {
                 tmp_chains.push_back(*pos);
                 continue;
            }
            if (first) tmp_chains.push_back(chain);
            first = false;
            _removed_chain_index_set.insert((*pos)->index());
            delete *pos;
       }
       _chains = tmp_chains;

       GenerateInternalOrderIndex();

       chain->get_seq(seqs);

       return _successful;
}

bool Molecule::SplitPolymer(const std::vector<std::vector<std::string> >& splits, const std::vector<std::vector<std::string> >& chain_info,
                            std::map<int, std::pair<std::vector<std::string>, std::vector<std::string> > >& splited_entityids)
{
       bool _successful = true;
       splited_entityids.clear();
       if (splits.empty() || chain_info.empty()) return _successful;

       RCSB::Chain* chain = GetPolyChain(chain_info[0][2]);
       if (!chain) {
            std::string error = "Chain " + chain_info[0][2] + " does not exist in Model " + String::IntToString(_Mol_ID) + ".\n";
            _logIo->message(error.c_str());
            _successful = false;
            return _successful;
       }

       std::vector<RCSB::Chain*> split_chains;
       chain->get_split_chains(splits, chain_info, _chain_index, split_chains);
       if (split_chains.empty()) {
            _successful = false;
            return _successful;
       }

       std::vector<std::string> seqs, data;

       std::vector<RCSB::Chain*> tmp_chains;
       tmp_chains.clear();
       for (std::vector<RCSB::Chain*>::const_iterator pos = _chains.begin(); pos != _chains.end(); ++pos) { 
            if (*pos == chain) {
                 _removed_chain_index_set.insert((*pos)->index());
                 delete *pos;
                 for (std::vector<RCSB::Chain*>::const_iterator spos = split_chains.begin(); spos != split_chains.end(); ++spos) {
                      if (((*spos)->chain_type() == "ATOMP") || ((*spos)->chain_type() == "ATOMN")) {
                           (*spos)->get_seq(seqs);
                           data.clear();
                           data.push_back("1");
                           data.push_back(String::IntToString((*spos)->SeqLen()));
                           data.push_back((*spos)->PDB_ChainID());
                           splited_entityids.insert(std::make_pair((*spos)->int_prev_entity_id(), std::make_pair(seqs, data)));
                      }
                      tmp_chains.push_back(*spos);
                 }
                 continue;
            }
            tmp_chains.push_back(*pos);
       }
       _chains = tmp_chains;

       GenerateInternalOrderIndex();

       return _successful;
}

bool Molecule::EditPolymer(const std::vector<std::vector<std::string> >& deletes)
{
       bool _successful = true;
       if (deletes.empty()) return _successful;

       // deletes[0]: Residue Name
       // deletes[1]: Residue position in sequence starting with 1
       // deletes[2]: New Asym ID
       // deletes[3]: Current PDB Chain ID
       // deletes[4]: New entity ID

       RCSB::Chain* chain = GetPolyChain(deletes[0][3]);
       if (!chain) {
            std::string error = "Chain " + deletes[0][3] + " does not exist in Model " + String::IntToString(_Mol_ID) + ".\n";
            _logIo->message(error.c_str());
            _successful = false;
            return _successful;
       }

       std::vector<RCSB::Chain*> nonpolymer_chains;
       chain->get_nonpolymer_chains(deletes, _chain_index, nonpolymer_chains);

       if (nonpolymer_chains.empty()) return _successful;

       for (std::vector<RCSB::Chain*>::const_iterator vpos = nonpolymer_chains.begin(); vpos != nonpolymer_chains.end(); ++vpos) {
            _chains.push_back(*vpos);
       }

       GenerateInternalOrderIndex();

       return _successful;
}

void Molecule::MergeResidue(const std::vector<std::pair<std::string, std::vector<std::vector<std::string> > > >& merge_residue_list)
{
       if (merge_residue_list.empty()) return;

       bool merge_flag = false;
       std::vector<int> index_list;
       std::set<std::string> name_set;
       std::map<std::string, std::string> name_mapping;
       for (std::vector<std::pair<std::string, std::vector<std::vector<std::string> > > >::const_iterator
            mpos = merge_residue_list.begin(); mpos != merge_residue_list.end(); ++mpos) {
            index_list.clear();
            for (std::vector<std::vector<std::string> >::const_iterator rpos = mpos->second.begin(); rpos != mpos->second.end(); ++rpos) {
                 std::string cs = CompositeIndex::getIndex(*rpos);
                 std::multimap<std::string, int>::const_iterator ipos = _pdbIndex.find(cs);
                 if (ipos != _pdbIndex.end()) index_list.push_back(ipos->second);
            }
            if (index_list.size() != mpos->second.size()) continue;

            name_set.clear();
            RCSB::Residue* newRes = _residues[index_list[0]];
            newRes->set_ResName(mpos->first);
            RCSB::Atom* atom = newRes->GetFirstAtom();
            while (atom) {
                 name_set.insert(atom->pdb_atmnam());
                 atom->set_restype(mpos->first);
                 atom->set_pdb_resnam(mpos->first);
                 atom = newRes->GetNextAtom();
            }
            for (unsigned int i = 1; i < index_list.size(); ++i) {
                 name_mapping.clear();
                 RCSB::Residue* removeRes = _residues[index_list[i]];
                 atom = removeRes->GetFirstAtom();
                 while (atom) {
                      std::map<std::string, std::string>::const_iterator npos = name_mapping.find(atom->pdb_atmnam());
                      if (npos != name_mapping.end()) {
                           atom->set_atmtype(npos->second);
                           atom->set_pdb_atmnam(npos->second);
                      } else {
                           std::string new_name = get_unique_atom_name(name_set, atom->pdb_atmnam(), atom->atom_type());
                           atom->set_atmtype(new_name);
                           atom->set_pdb_atmnam(new_name);
                           name_mapping.insert(std::make_pair(atom->pdb_atmnam(), new_name));
                      }
                      atom->set_chnid(newRes->chnid());
                      atom->set_restype(mpos->first);
                      atom->set_resnum(newRes->res_no());
                      atom->set_pdb_chnid(newRes->pdb_chnid());
                      atom->set_pdb_resnam(mpos->first);
                      atom->set_pdb_resnum(newRes->pdb_res_no());
                      newRes->insert_a_atom(atom, _Mol_ID);

                      atom = removeRes->GetNextAtom();
                 }
                 // clear first before delete the residue so that atom records will not be deleted.
                 _removed_residue_index_set.insert(index_list[i]);
                 removeRes->clear();
                 delete removeRes;
                 _residues[index_list[i]] = NULL;
            }
            newRes->UpdateIndices();
            merge_flag = true;
       }
       if (!merge_flag) return;

       std::vector<RCSB::Residue*> tmp_residues;
       tmp_residues.clear();
       for (std::vector<RCSB::Residue*>::iterator rpos = _residues.begin(); rpos != _residues.end(); ++rpos) {
            if ((*rpos) == NULL) continue;
            tmp_residues.push_back(*rpos);
       }
       _residues = tmp_residues;
       update_residue_indices();
}

void Molecule::_get_PRD_involved_chains(const std::vector<std::vector<std::string> >& mapping, std::set<int>& chain_set, bool& _successful)
{
       // get unique chain set
       chain_set.clear();

       // value: pdbchnid_name_number_inscode
       std::set<std::string> unique_residue_set;
       unique_residue_set.clear();

       // Index for mapping array
       // vector[0]: model ID
       // vector[1]: PDB Chain ID
       // vector[2]: Residue Name
       // vector[3]: Residue Number
       // vector[4]: Insertion Code
       // vector[5]: Atom Name
       // vector[6]: New Residue Name
       // vector[7]: New Residue Number
       // vector[8]: New Atom Name
       // vector[9]: New polymer type
       // vector[10]: Component ID
       // vector[11]: Old polymer type
       // vector[12]: Old Asym ID
       // vector[13]: New PDB Chain ID
       // vector[14]: New Asym ID
       for (std::vector<std::vector<std::string> >::const_iterator pos = mapping.begin(); pos != mapping.end(); ++pos) {
            std::string cs = (*pos)[1] + "_" + (*pos)[2] + "_" + (*pos)[3] + "_" + (*pos)[4];
            if (unique_residue_set.find(cs) != unique_residue_set.end()) continue;

            unique_residue_set.insert(cs);

            std::multimap<std::string, int>::iterator mpos = _pdbIndex.find(cs);
            if (mpos == _pdbIndex.end()) {
                 std::string error = "Residue (" + (*pos)[1] + " " + (*pos)[2] + " " + (*pos)[3] + (*pos)[4] + ") does not eixt in Model "
                                   + String::IntToString(_Mol_ID) + ".\n";
                 _logIo->message(error.c_str());
                 _successful = false;
                 continue;
            }
            chain_set.insert(_residues[mpos->second]->chain_index());
       }
}

void Molecule::_get_PRD_involved_residue_and_atoms(const std::set<int>& chain_set, std::multimap<std::string, int>& residue_mapping,
                                                   std::multimap<std::string, RCSB::Atom*>& found_atoms, bool& _successful)
{
       residue_mapping.clear();
       found_atoms.clear();

       std::vector<RCSB::Residue*> residue_list;
       for (std::set<int>::const_iterator spos = chain_set.begin(); spos != chain_set.end(); ++spos) {
            bool is_removed = false;
            RCSB::Chain* chain = GetIndexChain(*spos, is_removed);
            if (chain == NULL) {
                 if (!is_removed) _logIo->message("Chain index out of range.\n");
                 _successful = false;
                 continue;
            }

            chain->GetFirstResidueList(residue_list);
            while (!residue_list.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator vpos = residue_list.begin(); vpos != residue_list.end(); ++vpos) {
                      std::string cs = CompositeIndex::getIndex((*vpos)->pdb_chnid(), (*vpos)->ResName(), (*vpos)->pdb_res_no(), (*vpos)->ins_code());
                      residue_mapping.insert(std::make_pair(cs, (*vpos)->index()));
                      _get_involved_atoms(*vpos, found_atoms);
                 }
                 chain->GetNextResidueList(residue_list);
            }
       }
}

void Molecule::_get_involved_atoms(RCSB::Residue* res, std::multimap<std::string, RCSB::Atom*>& found_atoms)
{
       std::vector<std::string> data;

       RCSB::Atom* atom = res->GetFirstAtom();
       while (atom) {
            data.clear();
            data.push_back(atom->pdb_chnid());
            data.push_back(atom->pdb_resnam());
            data.push_back(atom->pdb_resnum());
            data.push_back(atom->ins_code());
            data.push_back(atom->pdb_atmnam());
            std::string cs = CompositeIndex::getIndex(data);
            found_atoms.insert(std::make_pair(cs, atom));

            atom = res->GetNextAtom();
       }
}

bool Molecule::_get_new_residue_list(const std::vector<std::vector<std::string> >& mapping, const std::multimap<std::string, RCSB::Atom*>& found_atoms, const
                                     std::vector<std::vector<std::string> >& in_links, std::list<RCSB::Residue*>& residue_list, std::map<std::string, std::string>&
                                     polymer_type_mapping, std::vector<std::pair<std::string, std::vector<std::vector<RCSB::Atom*> > > >& out_links)
{
       residue_list.clear();
       polymer_type_mapping.clear();

       std::vector<std::vector<std::string> > missing_atoms;
       missing_atoms.clear();

       std::set<std::string> found_mapped_atom_set, t_set;
       found_mapped_atom_set.clear();
       for (std::vector<std::vector<std::string> >::const_iterator pos = mapping.begin(); pos != mapping.end(); ++pos) {
            std::string cs = (*pos)[1] + "_" + (*pos)[2] + "_" + (*pos)[3] + "_" + (*pos)[4] + "_" + (*pos)[5];
            std::pair<std::multimap<std::string, RCSB::Atom*>::const_iterator, std::multimap<std::string, RCSB::Atom*>::const_iterator> range = found_atoms.equal_range(cs);
            if (range.first == range.second) {
                 missing_atoms.push_back(*pos);
                 continue;
            }
            found_mapped_atom_set.insert(cs);
       }

       bool extra_atom_flag = false;
       std::vector<RCSB::Atom*> extra_atom_list, atoms;
       extra_atom_list.clear();

       std::string extra_atoms, record;
       extra_atoms.clear();
       for (std::multimap<std::string, RCSB::Atom*>::const_iterator apos = found_atoms.begin(); apos != found_atoms.end(); ++apos) {
            if (found_mapped_atom_set.find(apos->first) != found_mapped_atom_set.end()) continue;
            if (!apos->second->is_hydrogen()) extra_atom_flag = true;
            apos->second->writeAtomRecord(record);
            extra_atoms += record;
            extra_atom_list.push_back(apos->second);
       }

       std::map<std::string, std::vector<RCSB::Atom*> > connected_extra_h_atom_map;
       connected_extra_h_atom_map.clear();

       // extra hydrogen(s) only
       if (!extra_atom_flag && !extra_atom_list.empty()) {
            for (std::vector<RCSB::Atom*>::const_iterator lpos = extra_atom_list.begin(); lpos != extra_atom_list.end(); ++lpos) {
                 bool found_connected_atom = false;
                 for (std::multimap<std::string, RCSB::Atom*>::const_iterator apos = found_atoms.begin(); apos != found_atoms.end(); ++apos) {
                      if (apos->second->is_hydrogen()) continue;
                      if (!(*lpos)->alt_loc().empty() && !apos->second->alt_loc().empty() && ((*lpos)->alt_loc() != apos->second->alt_loc())) continue;
                      if (!BondUtil::is_a_bond(*lpos, apos->second)) continue;

                      std::map<std::string, std::vector<RCSB::Atom*> >::iterator mpos = connected_extra_h_atom_map.find(apos->first);
                      if (mpos != connected_extra_h_atom_map.end()) mpos->second.push_back(*lpos);
                      else {
                           atoms.clear();
                           atoms.push_back(*lpos);
                           connected_extra_h_atom_map.insert(std::make_pair(apos->first, atoms));
                      }

                      found_connected_atom = true;
                      break;
                 }
                 if (!found_connected_atom) {
                      extra_atom_flag = true;
                      break;
                 }
            }
       }

       record.clear();
       if (!missing_atoms.empty()) {
            if (missing_atoms.size() > 1)
                 record = "The following atoms can not be found:\n\n";
            else record = "The following atom can not be found:\n\n";
            for (std::vector<std::vector<std::string> >::const_iterator pos = missing_atoms.begin(); pos != missing_atoms.end(); ++pos) {
                 record += (*pos)[1] + " " + (*pos)[2] + " " + (*pos)[3] + (*pos)[4] + " " + (*pos)[5] + "\n";
            }
       }
       if (extra_atom_flag && !extra_atoms.empty()) {
            if (!record.empty()) record += "\n";
            if (extra_atom_list.size() > 1)
                 record += "The following atoms are not mapped:\n\n";
            else record += "The following atom is not mapped:\n\n";
            record += extra_atoms;
       }
       if (!record.empty()) {
            extra_atom_list.clear();
            connected_extra_h_atom_map.clear();
            _logIo->message(record.c_str());
            return false;
       }

       std::map<std::string, RCSB::Residue*> residue_mapping;
       residue_mapping.clear();

       std::map<std::string, int> last_number_mapping;
       last_number_mapping.clear();

       std::vector<std::vector<RCSB::Atom*> > atom_lists;
       std::map<std::string, std::vector<RCSB::Atom*> > comp_atom_mapping;
       comp_atom_mapping.clear();

       // Index for mapping array
       // vector[0]: model ID
       // vector[1]: PDB Chain ID
       // vector[2]: Residue Name
       // vector[3]: Residue Number
       // vector[4]: Insertion Code
       // vector[5]: Atom Name
       // vector[6]: New Residue Name
       // vector[7]: New Residue Number
       // vector[8]: New Atom Name
       // vector[9]: New polymer type
       // vector[10]: Component ID
       // vector[11]: Old polymer type
       // vector[12]: Old Asym ID
       // vector[13]: New PDB Chain ID
       // vector[14]: New Asym ID
       // vector[15]: New Entity ID

       // insert mapped atoms
       for (std::vector<std::vector<std::string> >::const_iterator pos = mapping.begin(); pos != mapping.end(); ++pos) {
            std::string cs = (*pos)[1] + "_" + (*pos)[2] + "_" + (*pos)[3] + "_" + (*pos)[4] + "_" + (*pos)[5];
            std::pair<std::multimap<std::string, RCSB::Atom*>::const_iterator, std::multimap<std::string, RCSB::Atom*>::const_iterator> range = found_atoms.equal_range(cs);
            if (range.first == range.second) continue;

            if ((!(*pos)[1].empty()) && (!(*pos)[13].empty()) && ((*pos)[1] != (*pos)[13])) {
                 std::map<std::string, std::set<std::string> >::iterator mpos = _PDBchainID_changed_map.find((*pos)[13]);
                 if (mpos != _PDBchainID_changed_map.end()) mpos->second.insert((*pos)[1]);
                 else {
                      t_set.clear();
                      t_set.insert((*pos)[1]);
                      _PDBchainID_changed_map.insert(std::make_pair((*pos)[13], t_set));
                 }
            }

            for (std::multimap<std::string, RCSB::Atom*>::const_iterator apos = range.first; apos != range.second; ++apos) {
                 apos->second->set_chnid((*pos)[14]);
                 apos->second->set_pdb_chnid((*pos)[13]);
                 apos->second->set_restype((*pos)[6]);
                 apos->second->set_pdb_resnam((*pos)[6]);
                 apos->second->set_resnum((*pos)[7]);
                 apos->second->set_pdb_resnum((*pos)[7]);
                 apos->second->set_ins_code("");
                 apos->second->set_atmtype((*pos)[8]);
                 apos->second->set_pdb_atmnam((*pos)[8]);
                 if (!(*pos)[15].empty()) apos->second->setValue((*pos)[15], 23);
                 _insert_a_atom(apos->second, residue_mapping, last_number_mapping, residue_list, false, (*pos)[9]);
                 std::string idx = (*pos)[13] + "_" + (*pos)[6] + "_" + (*pos)[7];
                 if (polymer_type_mapping.find(idx) == polymer_type_mapping.end()) polymer_type_mapping.insert(std::make_pair(idx, (*pos)[9]));

                 if (in_links.empty()) continue;

                 std::string comp_atom_id = (*pos)[6] + "_" + (*pos)[7] + "_" + (*pos)[8] + "_" + (*pos)[10];
                 std::map<std::string, std::vector<RCSB::Atom*> >::iterator mpos = comp_atom_mapping.find(comp_atom_id);
                 if (mpos != comp_atom_mapping.end()) mpos->second.push_back(apos->second);
                 else {
                      atoms.clear();
                      atoms.push_back(apos->second);
                      comp_atom_mapping.insert(std::make_pair(comp_atom_id, atoms));
                 }
            }
       }

       // insert extra hydrogen atom(s) into proper residue(s)
       for (std::map<std::string, std::vector<RCSB::Atom*> >::const_iterator cpos = connected_extra_h_atom_map.begin(); cpos != connected_extra_h_atom_map.end(); ++cpos) {
            std::multimap<std::string, RCSB::Atom*>::const_iterator apos = found_atoms.find(cpos->first);
            if (apos == found_atoms.end()) continue;

            std::string cs = CompositeIndex::getIndex(apos->second->pdb_chnid(), apos->second->pdb_resnam(), apos->second->pdb_resnum(), apos->second->ins_code());
            std::map<std::string, RCSB::Residue*>::const_iterator rpos = residue_mapping.find(cs);
            if (rpos == residue_mapping.end()) continue;

            for (std::vector<RCSB::Atom*>::const_iterator vpos = cpos->second.begin(); vpos != cpos->second.end(); ++vpos) {
                 // not rename hydrogen atom name
                 (*vpos)->set_chnid(apos->second->chnid());
                 (*vpos)->set_pdb_chnid(apos->second->pdb_chnid());
                 (*vpos)->set_restype(apos->second->restype());
                 (*vpos)->set_pdb_resnam(apos->second->pdb_resnam());
                 (*vpos)->set_resnum(apos->second->resnum());
                 (*vpos)->set_pdb_resnum(apos->second->pdb_resnum());
                 (*vpos)->set_ins_code(apos->second->ins_code());
                 rpos->second->insert_a_atom(*vpos, _Mol_ID);
            }
       }

       extra_atom_list.clear();
       connected_extra_h_atom_map.clear();
       residue_mapping.clear();
       last_number_mapping.clear();

       if (comp_atom_mapping.empty()) return true;

       for (std::vector<std::vector<std::string> >::const_iterator pos = in_links.begin(); pos != in_links.end(); ++pos) {
            std::string cs = (*pos)[0] + "_" + (*pos)[1] + "_" + (*pos)[2] + "_" + (*pos)[3];
            std::map<std::string, std::vector<RCSB::Atom*> >::const_iterator mpos1 = comp_atom_mapping.find(cs);
            if (mpos1 == comp_atom_mapping.end()) continue;
            cs = (*pos)[4] + "_" + (*pos)[5] + "_" + (*pos)[6] + "_" + (*pos)[7];
            std::map<std::string, std::vector<RCSB::Atom*> >::const_iterator mpos2 = comp_atom_mapping.find(cs);
            if (mpos2 == comp_atom_mapping.end()) continue;

            atom_lists.clear();
            atom_lists.push_back(mpos1->second);
            atom_lists.push_back(mpos2->second);
            out_links.push_back(std::make_pair((*pos)[8], atom_lists));
       }

       return true;
}

void Molecule::_remove_old_residues_and_chains(const std::multimap<std::string, int>& residue_mapping, const std::set<int>& chain_set)
{
       // removed old residues and indices 
       for (std::multimap<std::string, int>::const_iterator mpos = residue_mapping.begin(); mpos != residue_mapping.end(); ++mpos) {
            _removed_residue_index_set.insert(_residues[mpos->second]->index());
            _residues[mpos->second]->clear();
            delete _residues[mpos->second];
            _residues[mpos->second] = NULL;
            _pdbIndex.erase(mpos->first);
       }

       // removed old chains
       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            if (chain_set.find((*pos)->index()) == chain_set.end()) continue;
            _removed_chain_index_set.insert((*pos)->index());
            delete *pos;
            *pos = NULL;
       }
}

void Molecule::_add_new_residues_and_chains(std::list<RCSB::Residue*>& residue_list, const std::map<std::string, SEQ>& seqs,
                                            const std::string& assigned_chain_type)
{
       std::map<std::string, RCSB::Chain*> chain_mapping;
       chain_mapping.clear();

       for (std::list<RCSB::Residue*>::iterator pos = residue_list.begin(); pos != residue_list.end(); ++pos) {
            (*pos)->set_alt_loc();
            (*pos)->get_atom_type();
            (*pos)->reorder_atoms();
            insert_a_residue(*pos);

            RCSB::Chain* chain = NULL;
            RCSB::Atom* atom = (*pos)->GetFirstAtom();
            std::map<std::string, RCSB::Chain*>::const_iterator mpos = chain_mapping.find(atom->chnid());
            if (mpos == chain_mapping.end()) {
                 chain = insert_a_chain(*pos, true);
                 chain_mapping.insert(std::make_pair(atom->chnid(), chain));
            } else chain = mpos->second;

            chain->insert_a_residue(*pos);
       }
       residue_list.clear();

       std::map<std::string, std::pair<std::string, std::string> > missing_residue_numbering_mapping;
       missing_residue_numbering_mapping.clear();

       std::string chain_type;
       std::map<std::string, int> all_types;
       for (std::map<std::string, RCSB::Chain*>::iterator mpos = chain_mapping.begin(); mpos != chain_mapping.end(); ++mpos) {
            chain_type.clear();
            std::map<std::string, SEQ>::const_iterator spos = seqs.find(mpos->second->ChainID());
            if (spos != seqs.end()) {
                 mpos->second->seq_alignment(spos->second, missing_residue_numbering_mapping);
                 chain_type = spos->second.chain_type;
            } else if (!assigned_chain_type.empty()) chain_type = assigned_chain_type;
            if (chain_type.empty()) {
                 all_types.clear();
                 bool has_ATOMP = false;
                 int count = 0;
                 RCSB::Residue* residue = mpos->second->GetFirstResidue();
                 while (residue) {
                      count++;
                      add_type(all_types, residue->chain_type());
                      if (residue->chain_type() == "ATOMP") has_ATOMP = true;
                      residue = mpos->second->GetNextResidue();
                 }
                 chain_type = get_type(all_types);
                 if ((chain_type != "ATOMP") && has_ATOMP && (count > 1) && (count < 11)) chain_type = "ATOMP";
            }
            mpos->second->set_chain_type(chain_type);
            if (chain_type == "ATOMP" || chain_type == "ATOMN")
                 mpos->second->set_has_sequence();
       }
       chain_mapping.clear();
}

bool Molecule::UpdateInstance(const std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > >& atom_mapping)
{
       if (atom_mapping.empty()) return true;

       std::vector<std::string> data;

       std::set<std::string> unique_index;
       unique_index.clear();

       // bool _successful = true;

       for (std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > >::const_iterator
            mpos = atom_mapping.begin(); mpos != atom_mapping.end(); ++mpos) {
            get_wordarray(data, mpos->first, "_");
            std::string pdb_chnid = data[1];
            std::string res_name  = data[2];
            std::string res_num   = data[3];
            std::string ins_code  = "";
            if (data.size() == 5) ins_code = data[4];
            std::string index = CompositeIndex::getIndex(pdb_chnid, res_name, res_num, ins_code);
            if (unique_index.find(index) != unique_index.end()) continue;

            unique_index.insert(index);

            RCSB::Residue* res = find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!res) {
                 std::string error = "Residue (" + pdb_chnid + " " + res_name + " " + res_num + ins_code + ") does not eixt in Model "
                                   + String::IntToString(_Mol_ID) + ".\n";
                 _logIo->message(error.c_str());
                 // _successful = false;
                 continue;
            }

            if (!res->update_atom_mapping(_Mol_ID, mpos->second.first, mpos->second.second)) {
                 // _successful = false;
                 continue;
            }
       
            if (res_name != mpos->second.first) {
                 bool is_removed = false;
                 RCSB::Chain* chain = GetIndexChain(res->chain_index(), is_removed);
                 if (!chain) {
                      if (!!is_removed) {
                           std::string error = "Search chain index " + String::IntToString(res->chain_index()) + " in Model "
                                             + String::IntToString(_Mol_ID) + " failed.\n";
                           _logIo->message(error.c_str());
                      }
                      // _successful = false;
                      continue;
                 }
                 if (!chain->Update_SeqRes(res_name, res_num, ins_code, mpos->second.first)) {
                      // _successful = false;
                      continue;
                 }
            }
       }

       return true;
}

bool Molecule::update_residue_nomenclature(const std::map<std::string, std::map<std::string, std::string> >& mapping)
{
       // key: chain_index
       // pair.first: residue_index
       // pair.second: new pdb chainid
       std::map<int, std::list<std::pair<int, std::string> > > chain_index_mapping;
       chain_index_mapping.clear();

       std::list<std::pair<int, std::string> > t_vec;

       bool _successful = true;
       std::vector<std::string> data;
       std::string new_chainid;

       for (std::map<std::string, std::map<std::string, std::string> >::const_iterator mpos = mapping.begin(); mpos != mapping.end(); ++mpos) {
            get_wordarray_delimit_by_string(data, mpos->first, "_");

            RCSB::Residue* res = find_pdb_residue(data[0], data[1], data[2], data[3]);
            if (res) {
                 res->update_nomenclature(mpos->second);
                 bool is_removed = false;
                 RCSB::Chain* chain = GetIndexChain(res->chain_index(), is_removed);
                 if (chain) chain->Update_SeqRes(data[1], data[2], data[3], mpos->second);

                 new_chainid.clear();
                 std::map<std::string, std::string>::const_iterator mmpos = mpos->second.find("chainid");
                 if (mmpos != mpos->second.end()) new_chainid = mmpos->second;

                 std::map<int, std::list<std::pair<int, std::string> > >::iterator cpos = chain_index_mapping.find(res->chain_index());
                 if (cpos != chain_index_mapping.end()) {
                      if (!new_chainid.empty()) cpos->second.push_back(std::make_pair(res->index(), new_chainid));
                 } else {
                      t_vec.clear();
                      if (!new_chainid.empty()) t_vec.push_back(std::make_pair(res->index(), new_chainid));
                      chain_index_mapping.insert(std::make_pair(res->chain_index(), t_vec));
                 }
                 continue;
            } else if (mpos->second.find("chainid") == mpos->second.end()) {
                 RCSB::Chain* chain = GetPolyChain(data[0]);
                 if (chain && chain->Update_SeqRes(data[1], data[2], data[3], mpos->second)) {
                      std::map<int, std::list<std::pair<int, std::string> > >::iterator cpos = chain_index_mapping.find(chain->index());
                      if (cpos == chain_index_mapping.end()) {
                           t_vec.clear();
                           chain_index_mapping.insert(std::make_pair(chain->index(), t_vec));
                      }
                      continue;
                 }
            }

            std::string error = "Residue (" + data[0] + " " + data[1] + " " + data[2] + data[3] + ") does not eixt in Model "
                              + String::IntToString(_Mol_ID) + ".\n";
            _logIo->message(error.c_str());
            _successful = false;
       }

       if (_successful) _successful = _update_chain_nomenclature(chain_index_mapping);

       return _successful;
}

bool Molecule::_update_chain_nomenclature(const std::map<int, std::list<std::pair<int, std::string> > >& mapping)
{
       if (mapping.empty()) return true;

       std::set<int> removed_chains;
       removed_chains.clear();

       std::list<std::list<RCSB::Residue*> > relocate_residues;
       std::vector<RCSB::Chain*> chain_list;

       bool reorder_flag = false;

       for (std::map<int, std::list<std::pair<int, std::string> > >::const_iterator mpos = mapping.begin(); mpos != mapping.end(); ++mpos) {
            bool is_removed = false;
            RCSB::Chain* chn = GetIndexChain(mpos->first, is_removed);
            if (!chn) {
                 if (!is_removed) _logIo->message("Get index chain failed.\n");
                 return false;
            }

            chn->update_nomenclature(mpos->second, relocate_residues);
            for (std::list<std::list<RCSB::Residue*> >::const_iterator vpos = relocate_residues.begin(); vpos != relocate_residues.end(); ++vpos) {
                 if ((*vpos).front()->ResName() == "HOH" || (*vpos).front()->ResName() == "DOD") {
                      GetPdbChains((*vpos).front()->pdb_chnid(), chain_list);
                      bool found = false;
                      for (std::vector<RCSB::Chain*>::const_iterator cpos = chain_list.begin(); cpos != chain_list.end(); ++cpos) {
                           if ((*cpos)->empty()) continue;
                           RCSB::Residue* res = (*cpos)->GetFirstResidue();
                           if (!res) continue;
                           if (res->ResName() == "HOH" || res->ResName() == "DOD") {
                                found = true; 
                                (*cpos)->merge_waters(*vpos);
                                break;
                           }
                      }
                      if (found) continue;

                      insert_a_chain(*vpos);
                      reorder_flag = true;
                 }
            }

            if (chn->empty()) removed_chains.insert(chn->index());
       }

       if (!removed_chains.empty()) {
            std::vector<RCSB::Chain*> tmp_chain;
            tmp_chain.clear();
            for (std::vector<RCSB::Chain*>::const_iterator pos = _chains.begin(); pos != _chains.end(); ++pos) { 
                 RCSB::Chain* chn = *pos;
                 if (removed_chains.find(chn->index()) != removed_chains.end()) {
                      _removed_chain_index_set.insert(chn->index());
                      delete chn;
                 } else tmp_chain.push_back(chn);
            }
            _chains = tmp_chain;

            GenerateInternalOrderIndex();
       }

       if (reorder_flag) InternalOrder();

       return true;
}

void Molecule::update_numbering(const std::map<std::string, std::map<std::string, std::string> >& mapping)
{
       std::vector<RCSB::Chain*> polymer_chains, nonpolymer_chains, water_chains;
       std::map<std::string, std::vector<RCSB::Chain*> > chain_mapping;
       for (std::map<std::string, std::map<std::string, std::string> >::const_iterator mpos = mapping.begin(); mpos != mapping.end(); ++mpos) {
            GetChains(mpos->first, polymer_chains, nonpolymer_chains, water_chains);
            chain_mapping.clear();
            if (!polymer_chains.empty()) chain_mapping.insert(std::make_pair("Polymer", polymer_chains));
            if (!nonpolymer_chains.empty()) chain_mapping.insert(std::make_pair("Nonpolymer", nonpolymer_chains));
            if (!water_chains.empty()) chain_mapping.insert(std::make_pair("Water", water_chains));
            for (std::map<std::string, std::string>::const_iterator mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                 std::map<std::string, std::vector<RCSB::Chain*> >::const_iterator cpos = chain_mapping.find(mmpos->first);
                 if (cpos == chain_mapping.end()) continue;

                 int serial_no = atoi(mmpos->second.c_str()); 
                 for (std::vector<RCSB::Chain*>::const_iterator vpos = cpos->second.begin(); vpos != cpos->second.end(); ++vpos) {
                      (*vpos)->update_residues_nomenclature(false, true, serial_no, true);
                 }
            }
       }
}

void Molecule::update_numbering(const std::map<std::string, std::vector<std::string> >& mapping)
{
       for (std::map<std::string, std::vector<std::string> >::const_iterator mpos = mapping.begin(); mpos != mapping.end(); ++mpos) {
            RCSB::Chain* chain = GetPolyChain(mpos->first);
            if (!chain) continue;
            chain->update_residues_nomenclature(mpos->second, true);
       }
}

void Molecule::update_chainid(const std::map<std::string, std::string>& mapping)
{
       for (std::vector<RCSB::Chain*>::iterator pos = _chains.begin(); pos != _chains.end(); ++pos) {
            std::map<std::string, std::string>::const_iterator mpos = mapping.find((*pos)->PDB_ChainID());
            if (mpos != mapping.end()) (*pos)->set_PDB_ChainID_to_residues(mpos->second);
       }
}

void Molecule::_insert_chains(std::map<std::string, std::vector<RCSB::Chain*> >& mapping_chains, RCSB::Chain* chain)
{
       std::vector<RCSB::Chain*> t_vec;
       std::map<std::string, std::vector<RCSB::Chain*> >::iterator mpos = mapping_chains.find(chain->PDB_ChainID());
       if (mpos != mapping_chains.end()) mpos->second.push_back(chain);
       else {
            t_vec.clear();
            t_vec.push_back(chain);
            mapping_chains.insert(std::make_pair(chain->PDB_ChainID(), t_vec));
       }
}

void Molecule::_get_last_number(int& last_number, const std::vector<RCSB::Chain*>& chains)
{
       if (chains.empty()) return;

       for (std::vector<RCSB::Chain*>::const_iterator pos = chains.begin(); pos != chains.end(); ++pos) {
            for (unsigned int i = 0; i < (*pos)->SeqLen(); ++i) {
                 _FIELD* SeqRes = (*pos)->SeqRes(i);
                 int number = atoi(SeqRes->Field[4].c_str());
                 if (number > last_number) last_number = number;
            }
       }
}
