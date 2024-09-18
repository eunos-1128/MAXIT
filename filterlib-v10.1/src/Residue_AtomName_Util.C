/*
FILE:     Residue_AtomName_Util.C
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
#include <ctype.h>

#include "algorithm-util.h"
#include "algorithm-util.C"
#include "BondUtil.h"
#include "CompositeIndex.h"
#include "Element.h"
#include "GetPairList.h"
#include "GetPairList.C"
#include "GraphMatch.h"
#include "Residue.h"
#include "SeqCodeUtil.h"
#include "utillib.h"

#define NUM_NA_RESIDUES  6

static const char *__na_residues[NUM_NA_RESIDUES] = { "DA", "DC", "DG", "DI", "DT", "DU" };

#define NUM_AA_RESIDUES  22

static const char *__aa_residues[NUM_AA_RESIDUES] = { "TRP", "TYR", "PYL", "ARG", "LYS", "PHE", "HIS", "MSE", "MET", "GLN", "GLU",
                                                      "LEU", "ASN", "ASP", "PRO", "ILE", "THR", "VAL", "CYS", "SEC", "SER", "ALA" };

using namespace RCSB;

void Residue::correction_special_na_name(const bool& first_flag, const int& Mol_ID)
{
       std::string resname = SeqCodeUtil::GetStandardNaCode(_ResName);
       if (resname.empty()) {
            if (find_atom("O2'") || find_atom("O2*") || ((_uniqueAtomNames.size() == 1) && find_atom("P"))) {
                 if (_ResName == "ADE")
                      resname = "A";
                 else if (_ResName == "CYT")
                      resname = "C";
                 else if (_ResName == "GUA")
                      resname = "G";
            } else {
                 if (_ResName == "ADE")
                      resname = "DA";
                 else if (_ResName == "CYT")
                      resname = "DC";
                 else if (_ResName == "GUA")
                      resname = "DG";
            }
       }
       if (resname.empty()) return;

       if (_uniqueAtomNames.size() == 1) {
            if (find_atom("P")) {
                 _chain_type = "ATOMN";
                 if (_ResName != resname) {
                      _ResName = resname;
                      for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
                           (*pos)->set_restype(_ResName);
                           (*pos)->set_pdb_resnam(_ResName);
                      }
                      if (_atomNameAlts_Orininal.empty()) _atomNameAlts_Orininal = _atomNameAlts;
                      _UpdateIndices();
                 }
            }
            return;
       }

       bool has_checking_atoms = false;
       bool has_sugar_atom = false;
       if ((find_atom("C1'") || find_atom("C1*")) && (find_atom("C2'") || find_atom("C2*")) && (find_atom("C3'") || find_atom("C3*")) &&
           (find_atom("C4'") || find_atom("C4*")) && (find_atom("C5'") || find_atom("C5*")) && (find_atom("O3'") || find_atom("O3*")) &&
           (find_atom("O4'") || find_atom("O4*")) && (find_atom("O5'") || find_atom("O5*"))) has_sugar_atom = true;
       if ((first_flag || find_atom("P")) && has_sugar_atom) has_checking_atoms = true;

       int ret = _correction_atom_name(resname, Mol_ID, true, has_checking_atoms);
       if (ret == 0) return;

       _chain_type = "ATOMN";
       bool change_indices = false;
       if (ret == 1) change_indices = true;
       if (_ResName != resname) {
            _ResName = resname;
            for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
                 (*pos)->set_restype(_ResName);
                 (*pos)->set_pdb_resnam(_ResName);
            }
            change_indices = true;
       }
       if (!change_indices) return;

       if (_atomNameAlts_Orininal.empty()) _atomNameAlts_Orininal = _atomNameAlts;
       _UpdateIndices();
}

void Residue::correction_atom_name(const int& Mol_ID)
{
       if (_ResName == "UNL") return;

       if (_uniqueAtomNames.size() == 1) {
            try {
                 const ConnectFormat& drug = _ccDic->find_drug(_ResName);
                 const std::vector<AtomFormat>& atoms = drug.atoms();
                 if (atoms.size() == 1) {
                      for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
                           (*pos)->set_atmtype(atoms[0].atomname());
                           (*pos)->set_pdb_atmnam(atoms[0].atomname());
                           (*pos)->set_atom_type(atoms[0].atomtype());
                      }
                      if (_atomNameAlts_Orininal.empty()) _atomNameAlts_Orininal = _atomNameAlts;
                      _UpdateIndices();
                 }
            } catch (const std::exception& exc) {}
            return;
       }

       if (_correction_atom_name(_ResName, Mol_ID) != 1) return;

       if (_atomNameAlts_Orininal.empty()) _atomNameAlts_Orininal = _atomNameAlts;
       _UpdateIndices();
       reorder_atoms();
}

bool Residue::is_matched_with_CCD(std::vector<std::string>& extra_atom_list, const bool& check_extra_atom_flag)
{
       extra_atom_list.clear();
       try {
            const ConnectFormat& drug = _ccDic->find_drug(_ResName);

            std::vector<std::pair<std::string, std::string> > atomlist, bondlist;
            // key: equavalent hydrogen name
            // value.first: hydrogen type
            // value.second: matched hydrogen name
            std::map<std::string, std::pair<std::string, std::string> > hydrogenMap;
            _createScratchGraph(-1, atomlist, bondlist, hydrogenMap);
            if (atomlist.empty() || bondlist.empty()) return false;

            bool is_substructure_match = false;
            std::map<std::string, std::string> list;
            list.clear();
            std::vector<std::pair<std::string, std::string> > ref_a_list = drug.getAtomList();
            std::vector<std::pair<std::string, std::string> > ref_b_list = drug.getBondList();
            GraphMatch::GetMatch(ref_a_list, ref_b_list, atomlist, bondlist, list, is_substructure_match);

            if (!list.empty()) return true;

            if (check_extra_atom_flag) {
                 std::set<std::string> terminal_atom_set;
                 terminal_atom_set.clear();
                 // Remove hydrogen atoms
                 _get_heavy_atom_only_list(atomlist, bondlist, terminal_atom_set);
                 // Get terminal atoms
                 for (std::vector<std::pair<std::string, std::string> >::const_iterator vpos = ref_a_list.begin(); vpos != ref_a_list.end(); ++vpos) {
                      if (drug.isTerminalAtom(vpos->second)) terminal_atom_set.insert(vpos->second);
                 }
                 // Remove hydrogen and terminal atoms
                 _get_heavy_atom_only_list(ref_a_list, ref_b_list, terminal_atom_set);

                 is_substructure_match = false;
                 list.clear();
                 GraphMatch::GetMatch(atomlist, bondlist, ref_a_list, ref_b_list, list, is_substructure_match, REF_TO_TARGET, true, false);
                 if (is_substructure_match && !list.empty()) {
                      for (std::vector<std::pair<std::string, std::string> >::const_iterator apos = atomlist.begin(); apos != atomlist.end(); ++apos) {
                           if ((apos->first == "H") || (apos->first == "D") || (list.find(apos->second) != list.end())) continue;
                           extra_atom_list.push_back(apos->second);
                      }
                      if (extra_atom_list.empty()) return true;
                 }
            }
       } catch (const std::exception& exc) {}

       return false;
}

bool Residue::has_modification(const std::string& chain_type)
{
       std::vector<std::pair<std::string, std::string> > atomlist, bondlist;
       std::map<std::string, std::pair<std::string, std::string> > hydrogenMap;
       // Create the graph with heavy atoms only
       _createScratchGraph(-1, atomlist, bondlist, hydrogenMap, true);
       // Creating graph failed.
       if (atomlist.empty() || bondlist.empty()) return false;

       bool is_substructure_match = false;
       std::map<std::string, std::string> list;

       try {
            const ConnectFormat& drug = _ccDic->find_drug(_ResName);
            if (drug.getMetaData("pdbx_release_status") != "REF_ONLY") {
                 std::vector<std::pair<std::string, std::string> > ref_a_list = drug.getAtomList();
                 std::vector<std::pair<std::string, std::string> > ref_b_list = drug.getBondList();
                 list.clear();
                 GraphMatch::GetMatch(ref_a_list, ref_b_list, atomlist, bondlist, list, is_substructure_match);
                 // Match with CCD definition
                 if (!list.empty()) return false;
            }
       } catch (const std::exception& exc) {}

       int num_of_residues = 0;
       const char** residue_names = NULL;
       if (chain_type == "ATOMN") {
            num_of_residues = NUM_NA_RESIDUES;
            residue_names = __na_residues;
       } else if (chain_type == "ATOMP") {
            num_of_residues = NUM_AA_RESIDUES;
            residue_names = __aa_residues;
       }

       std::vector<std::set<std::string> > atom_sets;
       for (int i = 0; i < num_of_residues; ++i) {
            try {
                 const ConnectFormat& drug = _ccDic->find_drug(residue_names[i]);
                 std::vector<std::pair<std::string, std::string> > ref_a_list = drug.getAtomList();
                 std::vector<std::pair<std::string, std::string> > ref_b_list = drug.getBondList();

                 std::set<std::string> terminal_atom_set;
                 terminal_atom_set.clear();
                 // Get terminal atoms
                 for (std::vector<std::pair<std::string, std::string> >::const_iterator vpos = ref_a_list.begin(); vpos != ref_a_list.end(); ++vpos) {
                      if (drug.isTerminalAtom(vpos->second)) terminal_atom_set.insert(vpos->second);
                 }
                 // Remove hydrogen and terminal atoms
                 _get_heavy_atom_only_list(ref_a_list, ref_b_list, terminal_atom_set);
                 if ((atomlist.size() + 3) < ref_a_list.size()) continue;
                 // Run sub graph match
                 is_substructure_match = false;
                 list.clear();
                 GraphMatch::GetMatch(atomlist, bondlist, ref_a_list, ref_b_list, list, is_substructure_match, REF_TO_TARGET, true, false);
                 if (is_substructure_match && !list.empty()) {
                      atom_sets.clear();
                      for (std::vector<std::pair<std::string, std::string> >::const_iterator apos = atomlist.begin(); apos != atomlist.end(); ++apos) {
                           if (list.find(apos->second) != list.end()) continue;
                           terminal_atom_set.clear();
                           terminal_atom_set.insert(apos->second);
                           atom_sets.push_back(terminal_atom_set);
                      }
                      if (atom_sets.empty()) continue;

                      clustering_with_merging(atom_sets, bondlist);
                      for (std::vector<std::set<std::string> >::const_iterator vpos = atom_sets.begin(); vpos != atom_sets.end(); ++vpos) {
                           if (vpos->size() > 2) return true;
                      }
                 }
            } catch (const std::exception& exc) {}
       }
       return false;
}

void Residue::check_linkage(const ConnectFormat& drug, std::list<std::string>& error_list, const bool& reset_flag)
{
       if (reset_flag) error_list.clear();

       std::map<std::string, std::vector<std::string> > bond_info_map;
       bond_info_map.clear();

       const std::vector<std::vector<std::string> >& bondList = drug.bonds();
       for (std::vector<std::vector<std::string> >::const_iterator vpos = bondList.begin(); vpos != bondList.end(); ++vpos) {
            RCSB::Atom* atom = find_atom((*vpos)[0]);
            if (!atom) continue;
            atom = find_atom((*vpos)[1]);
            if (!atom) continue;

            bond_info_map.insert(std::make_pair((*vpos)[0] + "_" + (*vpos)[1], *vpos));
       }

       std::map<std::string, std::string> atom_name_type_map;
       atom_name_type_map.clear();

       const std::vector<std::vector<std::string> >& atomList = drug.atomArray();
       for (std::vector<std::vector<std::string> >::const_iterator vpos = atomList.begin(); vpos != atomList.end(); ++vpos) {
            atom_name_type_map.insert(std::make_pair((*vpos)[0], (*vpos)[2]));
       }

       std::string error = _check_existing_linkage(bond_info_map);
       error += _check_missing_linkage(bond_info_map, atom_name_type_map);

       if (error.empty()) return;

       error_list.push_back(_pdb_chnid + " " + _ResName + " " + _pdb_res_no + _ins_code + ":\n" + error);
}

int Residue::_correction_atom_name(const std::string& resname, const int& Mol_ID, const bool& special_na_flag, const bool& has_checking_atoms)
{
       std::vector<ConnectFormat> drug_list;
       _ccDic->find_drugs(resname, drug_list);
       if (drug_list.empty()) return 0;

       std::vector<std::pair<std::string, std::string> > atomlist, bondlist;
       // key: equavalent hydrogen name
       // value.first: hydrogen type
       // value.second: matched hydrogen name
       std::map<std::string, std::pair<std::string, std::string> > hydrogenMap;
       _createScratchGraph(Mol_ID, atomlist, bondlist, hydrogenMap);

       if (atomlist.empty() || bondlist.empty()) return 0;

       std::map<std::string, std::string> atomtype;
       atomtype.clear();
       for (std::vector<std::pair<std::string, std::string> >::iterator pos = atomlist.begin(); pos != atomlist.end(); ++pos) {
            atomtype.insert(std::make_pair(pos->second, pos->first));
       }

       std::multimap<unsigned int, std::pair<int, std::map<std::string, std::string> > > all_mapping_list;
       all_mapping_list.clear();
       unsigned int max_mapped_vertex = 0;
       // record all mapping from component dictionary (including variants)
       // sort by matched vertex numbers

       bool is_substructure_match = 0;
       std::map<std::string, std::string> list, list_permutation;
       for (std::vector<ConnectFormat>::const_iterator dpos = drug_list.begin(); dpos != drug_list.end(); ++dpos) {
            std::vector<std::pair<std::string, std::string> > ref_a_list = dpos->getAtomList();
            std::vector<std::pair<std::string, std::string> > ref_b_list = dpos->getBondList();
            GraphMatch::GetMatch(ref_a_list, ref_b_list, atomlist, bondlist, list, is_substructure_match);
            if (list.empty()) continue;
            if (_is_OP123_permutation(resname, list, list_permutation)) list = list_permutation;

            if (is_substructure_match && special_na_flag) {
                 if (!has_checking_atoms) continue;
            }

            bool bad_mapping = 0;
            int count = 0;
            for (std::map<std::string, std::string>::const_iterator mpos = list.begin(); mpos != list.end(); ++mpos) {
                 std::string aName = mpos->first;
                 std::string eType = ""; 
                 std::map<std::string, std::string>::const_iterator mpos1 = atomtype.find(aName);
                 if (mpos1 != atomtype.end()) eType = mpos1->second;
                 if (eType == "D") {
                      for (unsigned int k = 0; k < aName.size(); k++) {
                           if (aName[k] == 'D') {
                                aName[k] = 'H';
                                break;
                           }
                      }

                 }

                 std::string cs = CompositeIndex::getIndex(_ResName, mpos->second);
                 if ((list.find(mpos->second) == list.end()) && (hydrogenMap.find(mpos->second) == hydrogenMap.end()) &&
                     (eType != "D") && (_atomNames.find(cs) != _atomNames.end())) {
                      bad_mapping = true;
                 }
                 if (aName != mpos->second) count++; 
            }
            if (bad_mapping) continue;

            if (list.size() == atomlist.size() && count == 0) {
                 if (hydrogenMap.empty()) return 2;
                 for (std::map<std::string, std::pair<std::string, std::string> >::const_iterator
                      mpos = hydrogenMap.begin(); mpos != hydrogenMap.end(); ++mpos) {
                      _update_atom_nomenclature(mpos->first, mpos->second.second, mpos->second.first);
                 }
                 return 1;
            }

            if (list.size() > max_mapped_vertex) max_mapped_vertex = list.size();
            all_mapping_list.insert(std::make_pair(list.size(), std::make_pair(count, list)));
       }

       std::multimap<unsigned int, std::map<std::string, std::string> > selected_mapping_list;
       selected_mapping_list.clear();
       // select mapping with maximum matched vertex numbers
       // sort by the number of atom name difference
       if (!all_mapping_list.empty()) {
            std::pair<std::multimap<unsigned int, std::pair<int, std::map<std::string, std::string> > >::iterator,
                      std::multimap<unsigned int, std::pair<int, std::map<std::string, std::string> > >::iterator>
                  range = all_mapping_list.equal_range(max_mapped_vertex);
            for (std::multimap<unsigned int, std::pair<int, std::map<std::string, std::string> > >::iterator
                 mmpos = range.first; mmpos != range.second; ++mmpos) {
                 selected_mapping_list.insert(std::make_pair(mmpos->second.first, mmpos->second.second));
            }
       }
       if (selected_mapping_list.empty() && (_ResName == resname)) {
            std::string cs = "In model " + String::IntToString(Mol_ID) + ", Can't find match for residue (";
            if (Mol_ID < 0) cs = "Can't find match for residue (";
            cs += GetFirstAtom()->pdb_chnid() + " " + GetFirstAtom()->pdb_resnam() + " " + GetFirstAtom()->pdb_resnum() + ").";
            _messageIo->insertMessage("residue_match_print", "warning", cs);
            cs += "\n";
            _logIo->message(cs.c_str());
            return 0;
       }

       std::multimap<unsigned int, std::map<std::string, std::string> >::iterator mmpos = selected_mapping_list.begin();
       if (mmpos->first == 0) return 0;

       for (std::map<std::string, std::string>::const_iterator mpos = mmpos->second.begin(); mpos != mmpos->second.end(); ++mpos) {
            std::string eType = "";
            std::map<std::string, std::string>::const_iterator mpos1 = atomtype.find(mpos->first);
            if (mpos1 != atomtype.end()) eType = mpos1->second;
            _update_atom_nomenclature(mpos->first, mpos->second, eType);
       }

       for (std::map<std::string, std::pair<std::string, std::string> >::const_iterator mpos = hydrogenMap.begin(); mpos != hydrogenMap.end(); ++mpos) {
            std::map<std::string, std::string>::const_iterator mpos1 = mmpos->second.find(mpos->second.second);
            if (mpos1 == mmpos->second.end()) continue;
            _update_atom_nomenclature(mpos->first, mpos1->second, mpos->second.first);
       }

       return 1;
}

void Residue::_update_atom_nomenclature(const std::string& old_name, const std::string& new_name, const std::string& atom_type)
{
       std::string cs = CompositeIndex::getIndex(_ResName, old_name);
       std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNames.equal_range(cs);
       if (range.first == range.second) return;

       _update_atom_nomenclature(range, new_name, atom_type);
}

void Residue::_update_atom_nomenclature(const std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator>& range,
                                        const std::string& new_name, const std::string& atom_type)
{
       std::string name = new_name;
       if ((atom_type == "D") && (name[0] == 'H')) name[0] = 'D';
       else if ((atom_type == "H") && (name[0] == 'D')) name[0] = 'H';

       for (std::multimap<std::string, int>::const_iterator mpos = range.first; mpos != range.second; ++mpos) {
            _atoms[mpos->second]->set_atmtype(name);
            _atoms[mpos->second]->set_pdb_atmnam(name);
       }
}

void Residue::_createScratchGraph(const int& Mol_ID, std::vector<std::pair<std::string, std::string> >& atomlist, std::vector<std::pair<std::string,
                                  std::string> >& bondlist, std::map<std::string, std::pair<std::string, std::string> >& hydrogenMap,
                                  const bool& heavy_atom_only_flag)
{
       atomlist.clear();
       bondlist.clear();
       hydrogenMap.clear();

       std::map<std::string, std::string> name_type;
       _getAtomNameTypeMapping(name_type);

       // key: atom name
       // value.first: alt_loc
       // value.second.first: atom_type
       std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > > heavyAtoms;

       // first: alt_loc
       // second.first: atom_type
       std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > hydrogenAtoms;

       _getAtoms(name_type, heavyAtoms, hydrogenAtoms);
       if ((heavyAtoms.size() < 2) && hydrogenAtoms.empty()) return;

       std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > > allConformerLists;

       _getAllConformerLists(heavyAtoms, allConformerLists);

       _getHeavyAtomAndBondList(Mol_ID, allConformerLists, atomlist, bondlist);
       if (atomlist.empty() || hydrogenAtoms.empty() || heavy_atom_only_flag) return;

       // key: heavy atom name
       // value.first: alt_loc
       // value.second.first:   hydrogen atom type
       // value.second.second:  hydrogen atom name
       std::map<std::string, std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > > hydrogenBondingInfo;
       _getHydrogenBondingInfo(heavyAtoms, hydrogenAtoms, hydrogenBondingInfo);

       std::vector<std::pair<std::string, std::string> > linkedHydrogens;
       for (std::map<std::string, std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > >::const_iterator
            mpos = hydrogenBondingInfo.begin(); mpos != hydrogenBondingInfo.end(); ++mpos) {
            _getHydrogenBondList(mpos->second, linkedHydrogens, hydrogenMap);
            for (std::vector<std::pair<std::string, std::string> >::const_iterator vpos = linkedHydrogens.begin(); vpos != linkedHydrogens.end(); ++vpos) {
                 atomlist.push_back(std::make_pair(vpos->first, vpos->second));
                 bondlist.push_back(std::make_pair(mpos->first, vpos->second));
            }
       }
}

void Residue::_change_hydrogen_names()
{
       std::map<std::string, std::string> name_type;
       _getAtomNameTypeMapping(name_type);

       // key: atom name
       // value.first: alt_loc
       // value.second.first: atom_type
       std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > > heavyAtoms;

       // first: alt_loc
       // second.first: atom_type
       std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > hydrogenAtoms;

       _getAtoms(name_type, heavyAtoms, hydrogenAtoms);
       if ((heavyAtoms.size() < 2) && hydrogenAtoms.empty()) return;


       std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > > allConformerLists;
       _getAllConformerLists(heavyAtoms, allConformerLists);

       // key: heavy atom name
       // value.first: alt_loc
       // value.second.first:   hydrogen atom type
       // value.second.second:  hydrogen atom name
       std::map<std::string, std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > > hydrogenBondingInfo;
       _getHydrogenBondingInfo(heavyAtoms, hydrogenAtoms, hydrogenBondingInfo);

       if (hydrogenBondingInfo.empty()) return;

       std::map<std::string, std::vector<std::string> > dic_heavy_hydrogen_atom_mapping;
       dic_heavy_hydrogen_atom_mapping.clear();

       try {
            std::set<std::string> hydrogen_set;
            std::vector<std::string> hydrogen_list;

            const ConnectFormat& drug = _ccDic->find_drug(_ResName);
            const std::map<std::string, std::vector<std::string> >& linkage_mapping = drug.getLinkedAtoms();

            for (std::map<std::string, std::vector<std::string> >::const_iterator mpos = linkage_mapping.begin(); mpos != linkage_mapping.end(); ++mpos) {
                 const AtomFormat& refAtom = drug.find_atom(mpos->first);
                 if (refAtom.ishydrogen()) continue;

                 hydrogen_set.clear();
                 for (std::vector<std::string>::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                      const AtomFormat& linkAtom = drug.find_atom(*vpos);
                      if (linkAtom.ishydrogen()) hydrogen_set.insert(*vpos);
                 }
                 if (hydrogen_set.empty()) continue;

                 hydrogen_list.clear();
                 for (std::set<std::string>::const_iterator spos = hydrogen_set.begin(); spos != hydrogen_set.end(); ++spos) {
                      hydrogen_list.push_back(*spos);
                 }
                 dic_heavy_hydrogen_atom_mapping.insert(std::make_pair(mpos->first, hydrogen_list));
            }
       } catch (const std::exception& exc) { return; }

       // vector[0]: old name
       // vector[1]: atom type
       // vector[2]: alt_loc
       // vector[3]: new name
       std::vector<std::vector<std::string> > rename_hydrogen_list;
       rename_hydrogen_list.clear();

       // vector[0]: atom name
       // vector[1]: alt_loc
       std::vector<std::vector<std::string> > remove_hydrogen_list;
       remove_hydrogen_list.clear();

       std::map<std::string, std::string> name_type_mapping;
       std::vector<std::vector<std::string> > name_type_list;
       std::vector<std::string> data;

       for (std::map<std::string, std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > >::const_iterator
            mpos = hydrogenBondingInfo.begin(); mpos != hydrogenBondingInfo.end(); ++mpos) {
            // skip backbone atoms (except for "CA")
            if ((mpos->first == "N") || (mpos->first == "C") || (mpos->first == "O") || (mpos->first == "OXT")) continue;

            for (std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >::const_iterator
                 vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 // re-order hydrogen names
                 name_type_mapping.clear();
                 for (std::vector<std::pair<std::string, std::string> >::const_iterator vvpos = vpos->second.begin(); vvpos != vpos->second.end(); ++vvpos) {
                      name_type_mapping.insert(std::make_pair(vvpos->second, vvpos->first));
                 }
                 name_type_list.clear();
                 for (std::map<std::string, std::string>::const_iterator mmpos = name_type_mapping.begin(); mmpos != name_type_mapping.end(); ++mmpos) {
                      data.clear();
                      data.push_back(mmpos->first);
                      data.push_back(mmpos->second);
                      name_type_list.push_back(data);
                 }

                 std::map<std::string, std::vector<std::string> >::const_iterator dpos = dic_heavy_hydrogen_atom_mapping.find(mpos->first);
                 if (dpos == dic_heavy_hydrogen_atom_mapping.end()) {
                      // insert remove hydrogen list
                      for (unsigned int i = 0; i < name_type_list.size(); ++i) {
                           data.clear();
                           data.push_back(name_type_list[i][0]);
                           data.push_back(vpos->first);
                           remove_hydrogen_list.push_back(data);
                      }
                 } else {
                      unsigned int len = name_type_list.size();
                      if (dpos->second.size() < len) len = dpos->second.size();

                      // insert rename hydrogen list
                      for (unsigned int i = 0; i < len; ++i) {
                           if (name_type_list[i][0] == dpos->second[i]) continue;

                           data = name_type_list[i];
                           data.push_back(vpos->first);
                           data.push_back(dpos->second[i]);
                           rename_hydrogen_list.push_back(data);
                      }

                      // insert remove hydrogen list
                      for (unsigned int i = len; i < name_type_list.size(); ++i) {
                           data.clear();
                           data.push_back(name_type_list[i][0]);
                           data.push_back(vpos->first);
                           remove_hydrogen_list.push_back(data);
                      }
                 }
            }
       }
       if (rename_hydrogen_list.empty() && remove_hydrogen_list.empty()) return;

       for (std::vector<std::vector<std::string> >::const_iterator vpos = rename_hydrogen_list.begin(); vpos != rename_hydrogen_list.end(); ++vpos) {
            std::string cs = CompositeIndex::getIndex(_ResName, (*vpos)[0], (*vpos)[2]);
            std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNameAlts.equal_range(cs);
            if (range.first == range.second) continue;

            _update_atom_nomenclature(range, (*vpos)[3], (*vpos)[1]);
       }

       std::set<int> delete_index;
       delete_index.clear();
       for (std::vector<std::vector<std::string> >::const_iterator vpos = remove_hydrogen_list.begin(); vpos != remove_hydrogen_list.end(); ++vpos) {
            std::string cs = CompositeIndex::getIndex(_ResName, (*vpos)[0], (*vpos)[1]);
            std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range = _atomNameAlts.equal_range(cs);
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
       _UpdateIndices();
}

void Residue::_getAtomNameTypeMapping(std::map<std::string, std::string>& name_type)
{
       name_type.clear();
       bool not_found_atom_type = false;
       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            name_type.insert(std::make_pair((*pos)->atmtype(), ""));
            if ((*pos)->atom_type().empty() || (*pos)->atom_type() == ".")
                 not_found_atom_type = true;
       }
       if (not_found_atom_type) {
            Element::getAtomSymbols(name_type);
       } else name_type.clear();
}

void Residue::_getAtoms(const std::map<std::string, std::string>& name_type, std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string,
                        RCSB::Atom*> > > >& heavyAtoms, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >& hydrogenAtoms)
{
       // separate heavy vs. hydrogen atoms to handle weird hydrogen naming from SHELXL program (DAOTHER-2810)
       heavyAtoms.clear();
       hydrogenAtoms.clear();

       // first: alt_loc
       // second.first: atom_type
       std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > tmpAtoms;

       for (std::vector<RCSB::Atom*>::iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            std::string type = (*pos)->atom_type();
            if (type.empty() || type == ".") {
                 std::map<std::string, std::string>::const_iterator mpos = name_type.find((*pos)->atmtype());
                 if (mpos != name_type.end()) type = mpos->second;
            }
            if (type == "Q") continue;

            if ((type == "H") || (type == "D")) {
                 hydrogenAtoms.push_back(std::make_pair((*pos)->alt_loc(), std::make_pair(type, *pos)));
            } else {
                 std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >::iterator
                     mpos = heavyAtoms.find((*pos)->atmtype());
                 if (mpos != heavyAtoms.end()) {
                      mpos->second.push_back(std::make_pair((*pos)->alt_loc(), std::make_pair(type, *pos)));
                 } else {
                      tmpAtoms.clear();
                      tmpAtoms.push_back(std::make_pair((*pos)->alt_loc(), std::make_pair(type, *pos)));
                      heavyAtoms.insert(std::make_pair((*pos)->atmtype(), tmpAtoms));
                 }
            }
       }
}

void Residue::_getAllConformerLists(const std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& heavyAtoms,
                                    std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& allConformerLists)
{      // list all possible alternative conformers. It should take care of weird case reported in DAOTHER-2000 ticket.
       allConformerLists.clear();

       std::set<std::string> unique_alts;
       // pair.first: alt_loc
       // pair.second.first: atom_type
       std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > regAtoms, altAtoms;
       std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > > allAtomLists;
       allAtomLists.clear();
       for (std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >::const_iterator
            mpos = heavyAtoms.begin(); mpos != heavyAtoms.end(); ++mpos) {
            if (mpos->second.size() < 2) {
                 allAtomLists.push_back(mpos->second);
                 continue;
            }
            regAtoms.clear();
            altAtoms.clear();
            unique_alts.clear();
            for (std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >::const_iterator
                 vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 if (vpos->first.empty()) {
                      regAtoms.push_back(*vpos);
                      break;
                 } else if (unique_alts.find(vpos->first) == unique_alts.end()) {
                      altAtoms.push_back(*vpos);
                      unique_alts.insert(vpos->first);
                 }
            }
            if (!regAtoms.empty()) {
                 // pick the atom without alt_loc ID
                 allAtomLists.push_back(regAtoms);
            } else if (!altAtoms.empty()) {
                 // pick the atom(s) with alt_loc ID
                 if (altAtoms.size() == 1) altAtoms[0].first = "";
                 allAtomLists.push_back(altAtoms);
            }
       }

       _getPairLists(0, allAtomLists, allConformerLists);
}

void Residue::_getPairLists(const unsigned int& index, const std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& a_list,
                            std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& p_list)
{
       if (index == 0) p_list.clear();
       if (index == a_list.size()) return;

       std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > > tmp_list;
       tmp_list.clear();

       std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > data;

       if (index == 0) {
            for (std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >::const_iterator
                 pos = a_list[index].begin(); pos != a_list[index].end(); ++pos) {
                 data.clear();
                 data.push_back(*pos);
                 tmp_list.push_back(data);
            } 
       } else {
            for (std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >::const_iterator
                 ppos = p_list.begin(); ppos != p_list.end(); ++ppos) {
                 std::set<std::string> existing_alts = _get_alt_loc_set(*ppos);
                 std::set<std::string> new_alts = _get_alt_loc_set(a_list[index]);
                 std::set<std::string> common_alts = _get_set_common_element(existing_alts, new_alts);
                 for (std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >::const_iterator
                      pos = a_list[index].begin(); pos != a_list[index].end(); ++pos) {
                      if (!pos->first.empty() && !common_alts.empty() && common_alts.find(pos->first) == common_alts.end()) continue;
                      data.clear();
                      for (std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >::const_iterator
                           pos1 = ppos->begin(); pos1 != ppos->end(); ++pos1) {
                           data.push_back(*pos1);
                      }
                      data.push_back(*pos);
                      tmp_list.push_back(data);
                 }
            }
       }
       p_list = tmp_list;
       _getPairLists(index + 1, a_list, p_list);
}

std::set<std::string> Residue::_get_alt_loc_set(const std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >& a_list)
{
       std::set<std::string> loc_set;
       loc_set.clear();
       for (std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >::const_iterator pos = a_list.begin(); pos != a_list.end(); ++pos) {
            if (!pos->first.empty()) loc_set.insert(pos->first);
       }
       return loc_set;
}

std::set<std::string> Residue::_get_set_common_element(const std::set<std::string>& set1, const std::set<std::string>& set2)
{
       std::set<std::string> common_set;
       common_set.clear();
       if (set1.empty() || set2.empty()) return common_set;
       for (std::set<std::string>::const_iterator pos = set1.begin(); pos != set1.end(); ++pos) {
            if (set2.find(*pos) == set2.end()) continue;
            common_set.insert(*pos);
       }
       return common_set;
}

void Residue::_getHeavyAtomAndBondList(const int& Mol_ID, const std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >&
                                       allConformerLists, std::vector<std::pair<std::string, std::string> >& atomlist,
                                       std::vector<std::pair<std::string, std::string> >& bondlist)
{
       atomlist.clear();
       bondlist.clear();

       // key: atom name pair - atomName1_atomName2
       // value: first.first  - first atom type
       // value: first.second - second atom type
       // value: second - distance
       std::map<std::string, std::pair<std::pair<std::string, std::string>, double> > dist_map;
       std::vector<std::map<std::string, std::pair<std::pair<std::string, std::string>, double> > > dist_map_array;
       dist_map_array.clear();

       std::vector<std::vector<std::pair<std::string, std::string> > > atomlist_array, bondlist_array;
       atomlist_array.clear();
       bondlist_array.clear();

       std::vector<std::set<std::string> > atom_sets;
       std::vector<std::vector<std::set<std::string> > > atom_sets_array;
       atom_sets_array.clear();

       for (std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >::const_iterator
            pos = allConformerLists.begin(); pos != allConformerLists.end(); ++pos) {
            _calDistMap(*pos, atomlist, bondlist, dist_map, atom_sets);
            clustering_with_merging(atom_sets, bondlist);
            if (atom_sets.size() < 2) return;
            
            dist_map_array.push_back(dist_map);
            atomlist_array.push_back(atomlist);
            bondlist_array.push_back(bondlist);
            atom_sets_array.push_back(atom_sets);
       }

       for (unsigned int i = 0; i < dist_map_array.size(); ++i) {
            _findLinkBetweenSubGraph(dist_map_array[i], atomlist_array[i], atom_sets_array[i], bondlist_array[i]);
            if (atom_sets_array[i].size() < 2) {
                 atomlist = atomlist_array[i];
                 bondlist = bondlist_array[i];
                 return;
            }
       }

       atomlist.clear();
       bondlist.clear();

       std::string cs = "";
       if (!atom_sets_array.empty()) {
            for (std::vector<std::set<std::string> >::iterator pos = atom_sets_array[0].begin(); pos != atom_sets_array[0].end(); ++pos) {
                 for (std::set<std::string>::iterator s1 = pos->begin(); s1 != pos->end(); ++s1) {
                      cs += *s1 + " ";
                 }
                 cs += "\n";
            }
       }

       if (Mol_ID < 0)
            cs += "\nResidue (";
       else cs += "\nIn model " + String::IntToString(Mol_ID) + ", Residue (";
       cs += GetFirstAtom()->pdb_chnid() + " " + GetFirstAtom()->pdb_resnam() + " " + GetFirstAtom()->pdb_resnum() + ") is not connected.";
       _messageIo->insertMessage("residue_match_print", "warning", cs);
       cs += "\n";
       _logIo->message(cs.c_str());
}

void Residue::_calDistMap(const std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >& atomList, std::vector<std::pair<std::string,
                          std::string> >& atomlist, std::vector<std::pair<std::string, std::string> >& bondlist, std::map<std::string,
                          std::pair<std::pair<std::string, std::string>, double> >& dist_map, std::vector<std::set<std::string> >& atom_sets)
{
       atomlist.clear();
       bondlist.clear();
       dist_map.clear();
       atom_sets.clear();

       std::set<std::string> tmp_set;
       for (unsigned int i = 0; i < atomList.size() - 1; ++i) {
            atomlist.push_back(std::make_pair(atomList[i].second.first, atomList[i].second.second->atmtype()));
            tmp_set.clear();
            tmp_set.insert(atomList[i].second.second->atmtype());
            atom_sets.push_back(tmp_set);
            for (unsigned int j = i + 1; j < atomList.size(); ++j) {
                 double dist = cal_distance(atomList[i].second.second, atomList[j].second.second);
                 if (BondUtil::is_a_bond_with_extension(atomList[i].second.first, atomList[j].second.first, dist)) {
                      bondlist.push_back(std::make_pair(atomList[i].second.second->atmtype(), atomList[j].second.second->atmtype()));
                 }
                 dist_map.insert(std::make_pair(atomList[i].second.second->atmtype() + "_" + atomList[j].second.second->atmtype(),
                                 std::make_pair(std::make_pair(atomList[i].second.first, atomList[j].second.first), dist)));
                 dist_map.insert(std::make_pair(atomList[j].second.second->atmtype() + "_" + atomList[i].second.second->atmtype(),
                                 std::make_pair(std::make_pair(atomList[i].second.first, atomList[j].second.first), dist)));
            }
       }
       atomlist.push_back(std::make_pair(atomList[atomList.size() - 1].second.first, atomList[atomList.size() - 1].second.second->atmtype()));
       tmp_set.clear();
       tmp_set.insert(atomList[atomList.size() - 1].second.second->atmtype());
       atom_sets.push_back(tmp_set);
}

void Residue::_findLinkBetweenSubGraph(const std::map<std::string, std::pair<std::pair<std::string, std::string>, double> >& dist_map, const
                                       std::vector<std::pair<std::string, std::string> >& atomlist, std::vector<std::set<std::string> >&
                                       atom_sets, std::vector<std::pair<std::string, std::string> >& bondlist)
{
       std::vector<std::pair<std::string, std::string> > links;
       links.clear();

       for (unsigned int i = 0; i < atom_sets.size() - 1; ++i) {
            for (unsigned int j = i + 1; j < atom_sets.size(); ++j) {
                 std::string atom1 = "";
                 std::string atom2 = "";
                 std::string type_1 = "";
                 std::string type_2 = "";
                 double minimum_dist = 100;
                 for (std::set<std::string>::const_iterator s1 = atom_sets[i].begin(); s1 != atom_sets[i].end(); ++s1) {
                      for (std::set<std::string>::const_iterator s2 = atom_sets[j].begin(); s2 != atom_sets[j].end(); ++s2) {
                           std::string key = *s1 + "_" + *s2;
                           std::map<std::string, std::pair<std::pair<std::string, std::string>, double> >::const_iterator mpos = dist_map.find(key);
                           if (mpos == dist_map.end()) continue;
                           if (mpos->second.second > minimum_dist) continue;
                           atom1 = *s1;
                           atom2 = *s2;
                           type_1 = mpos->second.first.first;
                           type_2 = mpos->second.first.second;
                           minimum_dist = mpos->second.second;
                      }
                 }
                 if (atom1.empty() || atom2.empty() || type_1.empty() || type_2.empty()) continue;
                 if (!BondUtil::is_a_link_loose(type_1, type_2, minimum_dist)) continue;
                 links.push_back(std::make_pair(atom1, atom2));
            }
       }
       if (links.empty()) return;

       clustering_with_merging(atom_sets, links);
       for (std::vector<std::pair<std::string, std::string> >::const_iterator pos = links.begin(); pos != links.end(); ++pos) {
            bondlist.push_back(*pos);
       }
}

void Residue::_getHydrogenBondingInfo(const std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& heavyAtoms,
                                      const std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >& hydrogenAtoms, std::map<std::string,
                                      std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > >& hydrogenBondingInfo)
{
       hydrogenBondingInfo.clear();

       std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > t_vec1;
       std::vector<std::pair<std::string, std::string> > t_vec2;
       for (std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >::const_iterator
            hpos = hydrogenAtoms.begin(); hpos != hydrogenAtoms.end(); ++hpos) {
            std::string atom = "";
            double minimum_dist = 100;
            for (std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >::const_iterator
                 mpos = heavyAtoms.begin(); mpos != heavyAtoms.end(); ++mpos) {
                 for (std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >::const_iterator
                      vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                      if (hpos->first != vpos->first && !hpos->first.empty() && !vpos->first.empty()) continue;

                      double dist = cal_distance(hpos->second.second, vpos->second.second);
                      if (!BondUtil::is_a_hydrogen_bonding(hpos->second.first, vpos->second.first, dist)) continue;

                      if (dist > minimum_dist) continue;

                      minimum_dist = dist;
                      atom = mpos->first;
                 }
            }
            if (atom.empty()) continue;

            std::map<std::string, std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > >::iterator
                mpos = hydrogenBondingInfo.find(atom);
            if (mpos == hydrogenBondingInfo.end()) {
                 t_vec2.clear();
                 t_vec2.push_back(std::make_pair(hpos->second.first, hpos->second.second->atmtype()));
                 t_vec1.clear();
                 t_vec1.push_back(std::make_pair(hpos->first, t_vec2));
                 hydrogenBondingInfo.insert(std::make_pair(atom, t_vec1));
            } else {
                 bool found = false;
                 for (std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >::iterator
                      vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                      if (hpos->first == vpos->first) {
                           found = true;
                           bool exist = false;
                           for (std::vector<std::pair<std::string, std::string> >::const_iterator
                                vvpos = vpos->second.begin(); vvpos != vpos->second.end(); ++vvpos) {
                                if (vvpos->second == hpos->second.second->atmtype()) {
                                     exist = true;
                                     break;
                                }
                           }
                           if (!exist) vpos->second.push_back(std::make_pair(hpos->second.first, hpos->second.second->atmtype()));
                           break;
                      }
                 }
                 if (found) continue;
                 t_vec2.clear();
                 t_vec2.push_back(std::make_pair(hpos->second.first, hpos->second.second->atmtype()));
                 mpos->second.push_back(std::make_pair(hpos->first, t_vec2));
            }
       }

       // key: heavy atom name
       // value.first: alt_loc
       // value.second.first:   hydrogen atom type
       // value.second.second:  hydrogen atom name
       // re-assign non-alt-loc hydrogen(s) to alt-loc list
       std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > second;
       for (std::map<std::string, std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > >::iterator
            mpos = hydrogenBondingInfo.begin(); mpos != hydrogenBondingInfo.end(); ++mpos) {
            if (mpos->second.size() < 2) continue;
            bool has_alt_loc = false;
            bool has_non_alt_loc = false;
            for (std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >::const_iterator
                 pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                 if (pos->first.empty()) has_non_alt_loc = true;
                 if (!pos->first.empty()) has_alt_loc = true;
            }
            if (!has_alt_loc || !has_non_alt_loc) continue;

            second.clear();
            for (std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >::const_iterator
                 pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                 if (pos->first.empty()) continue;
                 second.push_back(std::make_pair(pos->first, pos->second));
            }
            for (std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >::const_iterator
                 pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                 if (!pos->first.empty()) continue;
                 for (std::vector<std::pair<std::string, std::string> >::const_iterator ppos = pos->second.begin(); ppos != pos->second.end(); ++ppos) {
                      for (unsigned int i = 0; i < second.size(); ++i) {
                           second[i].second.push_back(std::make_pair(ppos->first, ppos->second));
                      }
                 }
            }
            mpos->second = second;
       }

       // re-order the hydrogen atom based on number-part of atom name
       std::vector<std::pair<std::string, std::string> > order;
       std::multimap<int, unsigned int> order_map;
       std::string number = "";
       for (std::map<std::string, std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > >::iterator
            mpos = hydrogenBondingInfo.begin(); mpos != hydrogenBondingInfo.end(); ++mpos) {
            for (std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >::iterator
                 pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                 if (pos->second.size() < 2) continue;

                 order_map.clear();
                 for (unsigned int i = 0; i < pos->second.size(); ++i) {
                      number.clear();
                      for (unsigned int j = 0; j < pos->second[i].second.size(); ++j) {
                           if (isdigit(pos->second[i].second[j])) number += pos->second[i].second[j];
                      }
                      if (number.empty()) {
                           order_map.clear();
                           break;
                      }
                      order_map.insert(std::make_pair(atoi(number.c_str()), i));
                 }
                 if (order_map.empty()) continue;

                 order.clear();
                 for (std::multimap<int, unsigned int>::const_iterator ompos = order_map.begin(); ompos != order_map.end(); ++ompos) {
                      order.push_back(pos->second[ompos->second]);
                 }
                 pos->second = order;
            }
       }
}

void Residue::_getHydrogenBondList(const std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >& hydrogens,
                                   std::vector<std::pair<std::string, std::string> >& linkedHydrogens,
                                   std::map<std::string, std::pair<std::string, std::string> >& hydrogenMap)
{
       linkedHydrogens.clear();
       if (hydrogens.size() == 1) {
            linkedHydrogens = hydrogens[0].second;
            return;
       }

       // finding the conformer with most hydrogens linked.
       unsigned int len = hydrogens[0].second.size();
       unsigned int idx = 0;
       for (unsigned int i = 1; i < hydrogens.size(); ++i) {
            if (hydrogens[i].second.size() > len) {
                 len = hydrogens[i].second.size();
                 idx = i;
            }
       }

       linkedHydrogens = hydrogens[idx].second;
       // key: hydrogen atom name
       // valye: H/D mixture
       std::map<std::string, std::string> atom_set;
       atom_set.clear();
       std::set<std::string> used_atom_set;
       for (std::vector<std::pair<std::string, std::string> >::const_iterator pos = linkedHydrogens.begin(); pos != linkedHydrogens.end(); ++pos) {
            std::string name = pos->second;
            if ((pos->first == "D") && (name[0] == 'D')) name[0] = 'H';
            atom_set.insert(std::make_pair(name, pos->second));
       }

       // mapping other conformers' hydrogen(s)
       for (unsigned int i = 0; i < hydrogens.size(); ++i) {
            if (i == idx) continue;
            used_atom_set.clear();
            for (std::vector<std::pair<std::string, std::string> >::const_iterator pos = hydrogens[i].second.begin(); pos != hydrogens[i].second.end(); ++pos) {
                 std::string name = pos->second;
                 if ((pos->first == "D") && (name[0] == 'D')) name[0] = 'H';
                 std::map<std::string, std::string>::const_iterator mpos = atom_set.find(name);

                 // finding the existing equivalent hydrogen
                 if (mpos != atom_set.end()) {
                      used_atom_set.insert(name);
                      // inserting equivalent hydrogen with different atom name (for H/D mixture)
                      if (pos->second != mpos->second) hydrogenMap.insert(std::make_pair(pos->second, std::make_pair(pos->first, mpos->second)));
                      continue;
                 }
                 // assigning equivalent hydrogen based on atom order
                 for (mpos = atom_set.begin(); mpos != atom_set.end(); ++mpos) {
                      if (used_atom_set.find(mpos->first) != used_atom_set.end()) continue;
                      used_atom_set.insert(name);
                      hydrogenMap.insert(std::make_pair(pos->second, std::make_pair(pos->first, mpos->second)));
                      break;
                 }
            }
       }
}

std::string Residue::_check_existing_linkage(std::map<std::string, std::vector<std::string> >& bond_info_map)
{
       std::string error = "";
       double dist, lower_limit, upper_limit;
       std::set<std::string> validated_keys;
       validated_keys.clear();
       for (unsigned int i = 0; i < _atoms.size() - 1; i++) {
            for (unsigned int j = i + 1; j < _atoms.size(); j++) {
                 if (!_atoms[i]->alt_loc().empty() && !_atoms[j]->alt_loc().empty() && (_atoms[i]->alt_loc() != _atoms[j]->alt_loc())) continue;
                 int type = BondUtil::is_a_bond_with_range(_atoms[i], _atoms[j], dist, lower_limit, upper_limit);
                 if (type == 0) continue;

                 std::string key1 = _atoms[i]->pdb_atmnam() + "_" + _atoms[j]->pdb_atmnam();
                 std::string key2 = _atoms[j]->pdb_atmnam() + "_" + _atoms[i]->pdb_atmnam();
                 if ((validated_keys.find(key1) != validated_keys.end()) || (validated_keys.find(key2) != validated_keys.end())) continue;

                 if (bond_info_map.find(key1) != bond_info_map.end()) {
                      bond_info_map.erase(key1);
                      validated_keys.insert(key1);
                      validated_keys.insert(key2);
                      continue;
                 } else if (bond_info_map.find(key2) != bond_info_map.end()) {
                      bond_info_map.erase(key2);
                      validated_keys.insert(key1);
                      validated_keys.insert(key2);
                      continue;
                 } else if (type != 1) continue;

                 std::string atom_i = _atoms[i]->pdb_atmnam();
                 if (!_atoms[i]->alt_loc().empty()) atom_i += "(" + _atoms[i]->alt_loc() + ")";
                 std::string atom_j = _atoms[j]->pdb_atmnam();
                 if (!_atoms[j]->alt_loc().empty()) atom_j += "(" + _atoms[j]->alt_loc() + ")";
                 error += "The distance ( " + FloatToString(dist, 0, 2) + " ) between atom '" + atom_i + "' and atom '" + atom_j
                        + "' is in bond distance range [ " + FloatToString(lower_limit, 0, 2) + ", " + FloatToString(upper_limit, 0, 2) + " ].\n";
            }
       }
       return error;
}

std::string Residue::_check_missing_linkage(const std::map<std::string, std::vector<std::string> >& bond_info_map,
                                            const std::map<std::string, std::string>& atom_name_type_map)
{
       std::string error = "";
       if (bond_info_map.empty()) return error;

       std::vector<RCSB::Atom*> atoms;
       std::vector<std::vector<RCSB::Atom*> > atom_lists, pair_lists;

       for (std::map<std::string, std::vector<std::string> >::const_iterator mpos = bond_info_map.begin(); mpos != bond_info_map.end(); ++mpos) {
            double lower_limit = -1.0;
            double upper_limit = -1.0;
            if ((mpos->second.size() == 7) && !mpos->second[5].empty() && String::IsNumber(mpos->second[5]) && !mpos->second[5].empty() &&
                 String::IsNumber(mpos->second[5])) {
                 lower_limit = atof(mpos->second[5].c_str()) - 10.0 * atof(mpos->second[6].c_str());
                 upper_limit = atof(mpos->second[5].c_str()) + 10.0 * atof(mpos->second[6].c_str());
            }
            std::map<std::string, std::string>::const_iterator apos0 = atom_name_type_map.find(mpos->second[0]);
            std::map<std::string, std::string>::const_iterator apos1 = atom_name_type_map.find(mpos->second[1]);
            if ((apos0 != atom_name_type_map.end()) && (apos1 != atom_name_type_map.end())) {
                 BondUtil::get_bond_range(apos0->second, apos1->second, lower_limit, upper_limit);
            }
            if (upper_limit < 0) continue;

            atom_lists.clear();
            find_atom(mpos->second[0], atoms);
            if (atoms.empty()) continue;
            atom_lists.push_back(atoms);

            find_atom(mpos->second[1], atoms);
            if (atoms.empty()) continue;
            atom_lists.push_back(atoms);

            GetPairList(atom_lists, pair_lists);
            if (pair_lists.empty()) continue;

            for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = pair_lists.begin(); apos != pair_lists.end(); ++apos) {
                 double dist = cal_distance((*apos)[0], (*apos)[1]);
                 if ((dist >= lower_limit) && (dist <= upper_limit)) continue;

                 std::string atom_i = (*apos)[0]->pdb_atmnam();
                 if (!(*apos)[0]->alt_loc().empty()) atom_i += "(" + (*apos)[0]->alt_loc() + ")";
                 std::string atom_j = (*apos)[1]->pdb_atmnam();
                 if (!(*apos)[1]->alt_loc().empty()) atom_j += "(" + (*apos)[1]->alt_loc() + ")";
                 error += "The distance ( " + FloatToString(dist, 0, 2) + " ) between atom '" + atom_i + "' and atom '" + atom_j
                        + "' is out of bond distance range [ " + FloatToString(lower_limit, 0, 2) + ", " + FloatToString(upper_limit, 0, 2) + " ].\n";
            }
       }

       return error;
}

void Residue::_get_heavy_atom_only_list(std::vector<std::pair<std::string, std::string> >& atomlist, std::vector<std::pair<std::string, std::string> >& bondlist,
                                        const std::set<std::string>& terminal_atom_set)
{
       std::vector<std::pair<std::string, std::string> > heavy_atomlist, heavy_bondlist;
       std::set<std::string> removed_atom_set;

       heavy_atomlist.clear();
       removed_atom_set.clear();
       for (std::vector<std::pair<std::string, std::string> >::const_iterator vpos = atomlist.begin(); vpos != atomlist.end(); ++vpos) {
            if ((vpos->first == "H") || (vpos->first == "D") || (terminal_atom_set.find(vpos->second) != terminal_atom_set.end()))
                 removed_atom_set.insert(vpos->second);
            else heavy_atomlist.push_back(*vpos);
       }
       atomlist = heavy_atomlist;

       heavy_bondlist.clear();
       for (std::vector<std::pair<std::string, std::string> >::const_iterator vpos = bondlist.begin(); vpos != bondlist.end(); ++vpos) {
            if ((removed_atom_set.find(vpos->first) != removed_atom_set.end()) || (removed_atom_set.find(vpos->second)!= removed_atom_set.end())) continue;
            heavy_bondlist.push_back(*vpos);
       }
       bondlist = heavy_bondlist;
}
