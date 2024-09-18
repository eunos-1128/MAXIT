/*
FILE:     ValidateObj.C
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
#include "GetMatthewsDistributionValues.h"
#include "GetPairList.h"
#include "GetPairList.C"
#include "SeqCodeUtil.h"
#include "ValidateObj.h"
#include "utillib.h"

#define NUM_A_A_B_A_LIST  4

const static char *_amino_acid_backbone_atom_list[NUM_A_A_B_A_LIST] = { "N", "CA", "C", "O" };

ValidateObj::ValidateObj(): FileObj()
{
       clear();
       for (int i = 0; i < NUM_A_A_B_A_LIST; ++i) _amino_acid_backbone_atom_set.insert(_amino_acid_backbone_atom_list[i]);
}

ValidateObj::~ValidateObj()
{
       clear();
}

void ValidateObj::clear()
{
       _exclude_extra_hydrogen_flag = false;
       _clear_contact();
       _clear_geometry();
       _clear_chirality_and_planarity();
       _cal_out_range_water_flag = false;
       _out_range_waters.clear();
       _unobs_or_zero_occ_residues.clear();
       _unobs_or_zero_occ_atoms.clear();
       _extra_atoms.clear();
       _links.clear();
       _ssbonds.clear();
       _sltbrgs.clear();
       _bspairs.clear();
       _defined_links.clear();
       _covale_link_mapping.clear();
       _amino_acid_backbone_atom_set.clear();
       _distance_checking_upper_limit_only_flag = false;
}

void ValidateObj::_clear_contact()
{
       _cal_a_contact_flag = false;
       _cal_s_contact_flag = false;
       _cal_special_position_flag = false;
       _a_contact.clear();
       _s_contact.clear();
       _special_position_atoms.clear();
}

void ValidateObj::_clear_geometry()
{
       _rmsbond = 0;
       _bond_count = 0;
       _rmsangle = 0;
       _angle_count = 0;
       _cis_trans_count = 0;
       _main_planarity_count = 0;
       _phi_psi_count = 0;
       _cal_bond_flag = false;
       _cal_angle_flag = false;
       _cal_cis_trans_flag = false;
       _cal_main_planarity_flag = false;
       _cal_phi_psi_flag = false;
       _cal_polymer_linkage_flag = false;
       _violated_bond.clear();
       _violated_angle.clear();
       _violated_cis_trans.clear();
       _violated_main_planarity.clear();
       _violated_phi_psi.clear();
       _violated_polymer_linkage.clear();
}

void ValidateObj::_clear_chirality_and_planarity()
{
       _cal_side_planarity_flag = false;
       _cal_chiral_center_flag = false;
       _violated_chiral_center.clear();
       _violated_side_planarity.clear();
       _side_planarity_count = 0;
}

std::string ValidateObj::_reformat_symmetry(const std::string& old_symmetry)
{
       std::string new_symmetry;
       new_symmetry.clear();
       if (old_symmetry.empty()) return new_symmetry;

       if (old_symmetry.find("_") != std::string::npos)
            new_symmetry = old_symmetry;
       else if (old_symmetry.size() > 3) {
            int len = old_symmetry.size() - 3;
            new_symmetry = old_symmetry.substr(0, len) + "_" + old_symmetry.substr(len);
       }
       if (new_symmetry == "1_555") new_symmetry.clear();
       return new_symmetry;
}

void ValidateObj::_insert_a_sltbrg(const _SLTBRG& sltbrg)
{
       for (std::list<_SLTBRG>::iterator
            pos = _sltbrgs.begin(); pos != _sltbrgs.end(); ++pos) {
            if (pos->mol_index == sltbrg.mol_index &&
               ((pos->fstAtom == sltbrg.fstAtom && pos->sndAtom == sltbrg.sndAtom &&
                 pos->SymOP_1 == sltbrg.SymOP_1 && pos->SymOP_2 == sltbrg.SymOP_2) ||
                (pos->fstAtom == sltbrg.sndAtom && pos->sndAtom == sltbrg.fstAtom &&
                 pos->SymOP_1 == sltbrg.SymOP_2 && pos->SymOP_2 == sltbrg.SymOP_1))) return;
       }

       _sltbrgs.push_back(sltbrg);
}

void ValidateObj::_insert_a_bspair(const _BSPAIR& bspair)
{
       for (std::list<_BSPAIR>::iterator
            pos = _bspairs.begin(); pos != _bspairs.end(); ++pos) {
            if (pos->mol_index == bspair.mol_index &&
               ((pos->fstAtom == bspair.fstAtom && pos->sndAtom == bspair.sndAtom &&
                 pos->SymOP_1 == bspair.SymOP_1 && pos->SymOP_2 == bspair.SymOP_2) ||
                (pos->fstAtom == bspair.sndAtom && pos->sndAtom == bspair.fstAtom &&
                 pos->SymOP_1 == bspair.SymOP_2 && pos->SymOP_2 == bspair.SymOP_1)) &&
                pos->details == bspair.details) return;
       }

       _bspairs.push_back(bspair);
}

void ValidateObj::_create_link_mapping()
{
       _defined_links.clear();
       _insert_defined_links(_links);
       _insert_defined_links(_ssbonds);
}

void ValidateObj::_insert_defined_links(const std::list<_LINK>& links)
{
       if (links.empty()) return;
       for (std::list<_LINK>::const_iterator
            pos = links.begin(); pos != links.end(); ++pos) {
            _insert_links(pos->fstAtom, pos->sndAtom, pos->dist);
            _insert_links(pos->sndAtom, pos->fstAtom, pos->dist);
       }
}

void ValidateObj::_insert_links(RCSB::Atom* atom1, RCSB::Atom* atom2, const std::string& dist)
{
       std::vector<std::string> data1, data2;
       data1.clear();
       data2.clear();
       atom1->getAtomIndex(data1);
       atom2->getAtomIndex(data2);
       std::string index1 = CompositeIndex::getIndex(data1);
       std::string index2 = CompositeIndex::getIndex(data2);
       std::string index = index1 + "_" + index2;
       std::map<std::string, std::set<double> >::iterator mpos = _defined_links.find(index);
       if (mpos == _defined_links.end()) {
            std::set<double> tset;
            tset.clear();
            tset.insert(atof(dist.c_str()));
            _defined_links.insert(std::make_pair(index, tset));
       } else mpos->second.insert(atof(dist.c_str()));

       index = index2 + "_" + index1;
       mpos = _defined_links.find(index);
       if (mpos == _defined_links.end()) {
            std::set<double> tset;
            tset.clear();
            tset.insert(atof(dist.c_str()));
            _defined_links.insert(std::make_pair(index, tset));
       } else mpos->second.insert(atof(dist.c_str()));
}

bool ValidateObj::_is_a_link(RCSB::Atom* atom1, RCSB::Atom* atom2, const double& dist)
{
       std::vector<std::string> data;
       data.clear();
       atom1->getAtomIndex(data);
       atom2->getAtomIndex(data);
       std::string index = CompositeIndex::getIndex(data);
       std::map<std::string, std::set<double> >::const_iterator
           mpos = _defined_links.find(index);
       if (mpos == _defined_links.end()) return false;
       return true;
/*
       for (std::set<double>::const_iterator
            spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
            if (fabs(dist - *spos) < 0.001) return true;
       }
       return false;
*/
}

int ValidateObj::_cal_z_value()
{
       if (_molecules.empty()) return 0;

       if (_cell.is_artifical()) return 1;

       unsigned int num_ncs_operations = 1;
       if (_CifObj) {
            Block& block = _CifObj->GetBlock(_firstBlockName);
            ISTable *t = getTablePtr(block, "struct_ncs_oper");
            if (t) {
                 std::string cs;
                 for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                      get_value_clean_lower(cs, t, i, "code");
                      if (cs == "generate") num_ncs_operations++;
                 }
            }
       }

       return (_cell.symops().size() * num_ncs_operations * _molecules[0]->MaxNumEntity());
}

bool ValidateObj::_cal_matthew_and_solvent(std::vector<double>& return_values)
{
       return_values.clear();

       if (_molecules.empty() || _cell.is_artifical() ||
         !(_experiment_type & EXPERIMENT_TYPE_XRAY)) return false;

       double weight = 0.0;

       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") {
                 weight += chain->weight();
            }
            chain = _molecules[0]->GetNextChain();
       }

       if (weight < 10) return false;

       double resolution = 200.0;
       if (_CifObj) {
            Block& block = _CifObj->GetBlock(_firstBlockName);
            ISTable *t = getTablePtr(block, "reflns");
            if (t) {
                 std::string cs;
                 get_value_clean(cs, t, 0, "d_resolution_high");
                 if (!cs.empty() && String::IsNumber(cs)) resolution = atof(cs.c_str());
            }
       }

       std::vector<double> distrib_values = get_matthews_coefficient_and_solvent_content_distrib_range_values(weight, resolution);

       double rad = acos(-1.0)/180.0;
       double alpha = atof(_cell.alpha().c_str());
       double beta = atof(_cell.beta().c_str());
       double gamma = atof(_cell.gamma().c_str());
       double a = atof(_cell.a().c_str());
       double b = atof(_cell.b().c_str());
       double c = atof(_cell.c().c_str());
       int nop = _cell.symops().size();

       double cosa = cos(alpha * rad);
       if (fabs(alpha - 90.0) < 0.001) cosa = 0.0;
       double cosb = cos(beta * rad);
       if (fabs(beta - 90.0) < 0.001) cosb = 0.0;
       double cosc = cos(gamma * rad);
       if (fabs(gamma - 90.0) < 0.001) cosc = 0.0;
       double vm = sqrt(1.00 + 2.0 * cosa * cosb * cosc - cosa * cosa
                 - cosb * cosb - cosc * cosc) * a * b * c;
       double matthew = vm / (weight * (double) nop);
       double solvent = 1 - 1.23 / matthew;

       return_values.push_back(matthew);
       return_values.push_back(solvent);
       for (std::vector<double>::const_iterator pos = distrib_values.begin(); pos != distrib_values.end(); ++pos) {
            return_values.push_back(*pos);
       }

       return true;
}

void ValidateObj::_check_polymer_linkage(RCSB::Chain* chain, const int& Mol_ID, const bool& check_gap_distance)
{
       if (chain->chain_type() != "ATOMP" && chain->chain_type() != "ATOMN") return;

       std::vector<std::string> warning_messages;
       _run_polymer_linkage_checking(chain, Mol_ID, check_gap_distance, warning_messages);

       if (chain->chain_type() != "ATOMP") return;

       std::vector<RCSB::Residue*> residues;
       chain->GetLastResidueList(residues);
       if (residues.empty()) return;

       RCSB::Atom* atom1 = residues[0]->find_atom("C");
       RCSB::Atom* atom2 = residues[0]->find_atom("OXT");
       if (!atom1 || !atom2) return;

       double dist = cal_distance(atom1, atom2);
       if (dist <= 1.6) return;

       std::string cs = "";
       if (Mol_ID >= 0) cs = "In model " + String::IntToString(Mol_ID) + ", ";

       cs += "Distance of C-OXT bond for residue (" + chain->PDB_ChainID() + " "
           + residues[0]->ResName() + " " + residues[0]->pdb_res_no()
           + residues[0]->ins_code() + ") is " + FloatToString(dist, 0, 2) + ".";
       _error_messages.insertMessage("oxt_distance_print", "warning", cs);
}

void ValidateObj::_run_polymer_linkage_checking(RCSB::Chain* chain, const int& Mol_ID, const bool& check_gap_distance, std::vector<std::string>& warning_messages)
{
       warning_messages.clear();

       std::vector<std::vector<RCSB::Residue*> > r_pair_list, residue_lists;
       std::vector<RCSB::Residue*> prev, curr;
       prev.clear();
       bool gap_flag = false;
       int gap_count = 0;

       for (unsigned int i = 0; i < chain->SeqLen(); ++i) {
            _FIELD* Seq = chain->SeqRes(i);
            if (!Seq) continue;

            if (Seq->ResIndex < 0) {
                 gap_flag = true;
                 gap_count++;
                 continue;
            }

            chain->GetResidueListByIndex(Seq->ResIndex, curr);
            if (curr.empty()) continue;

            if (!prev.empty()) {
                 residue_lists.clear();
                 residue_lists.push_back(prev);
                 residue_lists.push_back(curr);
                 GetPairList(residue_lists, r_pair_list);
                 if (!r_pair_list.empty()) {
                      for (std::vector<std::vector<RCSB::Residue*> >::const_iterator rpos = r_pair_list.begin(); rpos != r_pair_list.end(); ++rpos) {
                           _run_adjacent_linkage_checking(Mol_ID, chain->chain_type(), (*rpos)[0], (*rpos)[1], gap_flag);
                      }
                 }
            }

            if ((chain->chain_type() == "ATOMP") && check_gap_distance && gap_count) {
                 std::string dist_warning = _run_gap_distance_checking(Mol_ID, prev, curr, gap_count);
                 if (!dist_warning.empty()) {
                      std::string gap_warning = "The gap distance " + dist_warning + " is too large.\nThere ";
                      if (gap_count > 1) gap_warning += "are";
                      else gap_warning += "is";
                      gap_warning += " only " + String::IntToString(gap_count) + " residues between them. It is not enough to cover this gap.\n";
                      warning_messages.push_back(gap_warning);
                 }
            }

            prev = curr;
            gap_flag = false;
            gap_count = 0;
       }
}

void ValidateObj::_run_adjacent_linkage_checking(const int& Mol_ID, const std::string& chain_type, RCSB::Residue* prev, RCSB::Residue* curr, const bool& gap_flag)
{
       std::set<std::string> tailTerminalAtoms, headTerminalAtoms;
       std::string standardTailTerminalAtom = _ccDic->getTailTerminalAtoms(chain_type, prev->ResName(), tailTerminalAtoms);
       std::string standardHeadTerminalAtom = _ccDic->getHeadTerminalAtoms(chain_type, curr->ResName(), headTerminalAtoms);

       // key: distance
       // pair.first: is_a_link flag
       // pair.second.first: tail terminal atom
       // pair.second.second: head terminal atom
       std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > > terminal_atom_pairs, linked_atom_pairs;

       std::set<std::string> checked_pair_set;
       checked_pair_set.clear();

       _cal_adjacent_linkage_distances(prev, curr, tailTerminalAtoms, headTerminalAtoms, terminal_atom_pairs, checked_pair_set);
       _cal_adjacent_link_record_distances(prev, curr, linked_atom_pairs, checked_pair_set);

       if (checked_pair_set.find(standardTailTerminalAtom + "_" + standardHeadTerminalAtom) == checked_pair_set.end() &&
           !gap_flag && terminal_atom_pairs.empty() && linked_atom_pairs.empty())
           _cal_adjacent_linkage_distance(prev, curr, standardTailTerminalAtom, standardHeadTerminalAtom, terminal_atom_pairs, checked_pair_set);

       if (terminal_atom_pairs.empty()) {
            if (!gap_flag && !prev->is_connect(curr)) {
                 std::string error = "";
                 if (Mol_ID >= 0) error = "In model " + String::IntToString(Mol_ID) + ", ";
                 error += "Residues " + prev->ResName() + " " + prev->pdb_chnid() + prev->pdb_res_no() + prev->ins_code() + " and "
                        + curr->ResName() + " " + curr->pdb_chnid() + curr->pdb_res_no() + curr->ins_code()
                        + " that are next to each other in the sample sequence are not properly linked.";
                 _error_messages.insertMessage("mislinkage_print", "error", error, true);
            }
            return;
       }

       std::string inproper_linked_error = "";
       Value value;

       for (std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > >::const_iterator
            mpos = terminal_atom_pairs.begin(); mpos != terminal_atom_pairs.end(); ++mpos) {
            if ((gap_flag && !mpos->second.first) || (!gap_flag && mpos->second.first)) continue;

            std::string error = "";
            if (Mol_ID >= 0) error = "In model " + String::IntToString(Mol_ID) + ", ";
            error += "Residues " + mpos->second.second.first->pdb_resnam() + " " + mpos->second.second.first->pdb_chnid() + " "
                   + mpos->second.second.first->pdb_resnum() + mpos->second.second.first->ins_code() + " and "
                   + mpos->second.second.second->pdb_resnam() + " " + mpos->second.second.second->pdb_chnid() + " "
                   + mpos->second.second.second->pdb_resnum() + mpos->second.second.second->ins_code();
         
            if (gap_flag && mpos->second.first) {
                 if (inproper_linked_error.empty()) {
                      inproper_linked_error = error
                        + " are linked together in the model; however, there are residues between them in\nthe deposited polymeric sequence.";
                 }
            } else if (!gap_flag && !mpos->second.first) {
                 if (!linked_atom_pairs.empty() && (mpos->first > 3.0)) continue;

                 error += " that are next to each other in the sample sequence are not properly linked: ";
                 std::string atom_name_1 = mpos->second.second.first->pdb_atmnam();
                 if (!mpos->second.second.first->alt_loc().empty()) atom_name_1 += " (" + mpos->second.second.first->alt_loc() + " conformer)";
                 std::string atom_name_2 = mpos->second.second.second->pdb_atmnam();
                 if (!mpos->second.second.second->alt_loc().empty()) atom_name_2 += " (" + mpos->second.second.second->alt_loc() + " conformer)";
                 error += "distance between " + atom_name_1 + " and " + atom_name_2 + " is " + FloatToString(mpos->first, 0, 2) + ".";
                 _error_messages.insertMessage("mislinkage_print", "error", error, true);

                 // remove pdbx_validate_polymer_linkage records for case like 7KZN/D_1000253131
                 // if ((mpos->second.second.first->pdb_resnam() != "UNK") && (mpos->second.second.second->pdb_resnam() != "UNK")) {
                      value.clear();
                      value.set_Mol_ID(Mol_ID);
                      value.set_val(mpos->first);
                      value.insert_atom(mpos->second.second.first);
                      value.insert_atom(mpos->second.second.second);
                      _violated_polymer_linkage.push_back(value);
                 // }
            }
       }
       if (!inproper_linked_error.empty()) _error_messages.insertMessage("mislinkage_print", "error", inproper_linked_error);
}

void ValidateObj::_cal_adjacent_linkage_distances(RCSB::Residue* prev, RCSB::Residue* curr, const std::set<std::string>& tailAtoms, const std::set<std::string>& headAtoms,
                                                  std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > >& global_atom_pair_map,
                                                  std::set<std::string>& checked_pair_set)
{
       global_atom_pair_map.clear();
       if (tailAtoms.empty() || headAtoms.empty()) return;

       double global_distance = 100000000.0;
       std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > > local_atom_pair_map;

       for (std::set<std::string>::const_iterator tpos = tailAtoms.begin(); tpos != tailAtoms.end(); ++tpos) {
            for (std::set<std::string>::const_iterator hpos = headAtoms.begin(); hpos != headAtoms.end(); ++hpos) {
                 double local_distance = _cal_adjacent_linkage_distance(prev, curr, *tpos, *hpos, local_atom_pair_map, checked_pair_set);
                 if (local_atom_pair_map.empty()) continue;
                 if (local_distance < global_distance) {
                      global_distance = local_distance;
                      global_atom_pair_map = local_atom_pair_map;
                 }
            }
       }
}

void ValidateObj::_cal_adjacent_link_record_distances(RCSB::Residue* prev, RCSB::Residue* curr, std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > >&
                                                      global_atom_pair_map, std::set<std::string>& checked_pair_set)
{
       global_atom_pair_map.clear();

       std::string idx = prev->pdb_chnid() + "_" + prev->ResName() + "_" + prev->pdb_res_no() + "_" + prev->ins_code() + "_"
                       + curr->pdb_chnid() + "_" + curr->ResName() + "_" + curr->pdb_res_no() + "_" + curr->ins_code();
       std::map<std::string, std::vector<std::pair<std::string, std::string> > >::const_iterator mpos = _covale_link_mapping.find(idx);
       if (mpos == _covale_link_mapping.end()) return;

       double global_distance = 100000000.0;
       std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > > local_atom_pair_map;

       for (std::vector<std::pair<std::string, std::string> >::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
            double local_distance = _cal_adjacent_linkage_distance(prev, curr, vpos->first, vpos->second, local_atom_pair_map, checked_pair_set);
            if (local_atom_pair_map.empty()) continue;
            if (local_distance < global_distance) {
                 global_distance = local_distance;
                 global_atom_pair_map = local_atom_pair_map;
            }
       }
}

double ValidateObj::_cal_adjacent_linkage_distance(RCSB::Residue* prev, RCSB::Residue* curr, const std::string& tailAtom, const std::string& headAtom,
                                                   std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > >& atom_pair_map,
                                                   std::set<std::string>& checked_pair_set)
{
       atom_pair_map.clear();
       double global_distance = 100000000.0;

       if (prev->find_atom(tailAtom) == NULL) return global_distance;
       if (curr->find_atom(headAtom) == NULL) return global_distance;

       checked_pair_set.insert(tailAtom + "_" + headAtom);

       std::vector<std::pair<RCSB::Residue*, std::string> > res_atom_pair_list;
       std::vector<std::vector<RCSB::Atom*> > a_pair_list;

       res_atom_pair_list.clear();
       res_atom_pair_list.push_back(std::make_pair(prev, tailAtom));
       res_atom_pair_list.push_back(std::make_pair(curr, headAtom));
       _get_atom_pair_list(res_atom_pair_list, a_pair_list);
       if (a_pair_list.empty()) return global_distance;

       for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = a_pair_list.begin(); apos != a_pair_list.end(); ++apos) {
            double dist = cal_distance((*apos)[0], (*apos)[1]);
            bool flag = false;
            if (_distance_checking_upper_limit_only_flag) {
                 if (BondUtil::is_a_link_with_upper_limit((*apos)[0]->atom_type(), (*apos)[1]->atom_type(), dist)) flag = true;
            } else if (BondUtil::is_a_link((*apos)[0]->atom_type(), (*apos)[1]->atom_type(), dist)) flag = true;
            if (dist < global_distance) global_distance = dist;
            atom_pair_map.insert(std::make_pair(dist, std::make_pair(flag, std::make_pair((*apos)[0], (*apos)[1]))));
       }
       return global_distance;
}

std::string ValidateObj::_run_gap_distance_checking(const int& Mol_ID, const std::vector<RCSB::Residue*>& prev, const std::vector<RCSB::Residue*>& curr,
                                                    const int& gap_count)
{
       if (prev.empty()) return "";

       std::vector<std::vector<RCSB::Residue*> > r_pair_list, residue_lists;
       residue_lists.clear();
       residue_lists.push_back(prev);
       residue_lists.push_back(curr);
       GetPairList(residue_lists, r_pair_list);
       if (r_pair_list.empty()) return "";

       double cutoff_dist = ((double) (gap_count + 1)) * 3.8;
       if (cutoff_dist < 8.0) cutoff_dist = 8.0;

       for (std::vector<std::vector<RCSB::Residue*> >::const_iterator rpos = r_pair_list.begin(); rpos != r_pair_list.end(); ++rpos) {
            if (((*rpos)[0]->ResName() == "UNK") || ((*rpos)[1]->ResName() == "UNK")) continue;
            RCSB::Atom* atom1 = (*rpos)[0]->find_atom("CA");
            RCSB::Atom* atom2 = (*rpos)[1]->find_atom("CA");
            if (!atom1 || !atom2) {
                 atom1 = (*rpos)[0]->GetFirstAtom();
                 atom2 = (*rpos)[1]->GetFirstAtom();
            }
            if (!atom1 || !atom2) continue;

            double dist = cal_distance(atom1, atom2);
            if (dist < cutoff_dist) continue;

            std::string error = "The ";
            if (Mol_ID >= 0) error = "In model " + String::IntToString(Mol_ID) + ", the ";
            error += "distance between Residue " + atom1->pdb_chnid() + " " + atom1->pdb_resnam() + " "
                   + atom1->pdb_resnum() + atom1->ins_code() + " and Residue " + atom2->pdb_chnid() + " "
                   + atom2->pdb_resnam() + " " + atom2->pdb_resnum() + atom2->ins_code()
                   + " is " + FloatToString(dist, 0, 2) + " Angstrom. But there ";
            if (gap_count > 1)
                 error += "are only " + String::IntToString(gap_count) + " residues"; 
            else error += "is only " + String::IntToString(gap_count) + " residue";
            error += " (not enough sequence) to cover the gap region.";
            _error_messages.insertMessage("mislinkage_print", "warning", error);

            error = "( " + FloatToString(dist, 0, 2) + " Angstrom ) between residues ( " + atom1->pdb_chnid() + " " + atom1->pdb_resnam()
                  + " " + atom1->pdb_resnum() + atom1->ins_code() + " ) and " + "( " + atom2->pdb_chnid() + " " + atom2->pdb_resnam()
                  + " " + atom2->pdb_resnum() + atom2->ins_code() + " )"; 
            return error ;
       }
       return "";
}

void ValidateObj::_get_software_ordinal(std::set<std::string>& software_ordinal_set)
{
       software_ordinal_set.clear();

       if (!_CifObj) return;

       Block& block = _CifObj->GetBlock(_firstBlockName);
       ISTable *t = getTablePtr(block, "pdbx_nmr_refine");
       if (!t) return;

       std::string cs;
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(cs, t, i, "software_ordinal");
            if (!cs.empty()) software_ordinal_set.insert(cs);
       }
}

std::string ValidateObj::_get_struct_descriptor(Block& block)
{
       std::string descriptor, cs;
       std::string na_type = "";
       std::string na_descriptor = "";
       std::string prot_descriptor = "";
       std::string prot_ec_descriptor = "";
       std::set<std::string> descriptor_set;
       descriptor_set.clear();
       for (std::map<int, Entity>::const_iterator
            mpos = _entities.begin(); mpos != _entities.end(); ++mpos) {
            if (mpos->second.chain_type() != "ATOMN" && mpos->second.chain_type() != "ATOMP") continue;
            if (mpos->second.chain_type() == "ATOMN" && na_type.empty()) {
                 cs = mpos->second.getValue("poly_type");
                 if (cs == "polyribonucleotide")
                      na_type = "RNA";
                 else na_type = "DNA";
            }

            cs = mpos->second.getValue("pdbx_description");
            if (cs.empty()) continue;
            if (descriptor_set.find(cs) != descriptor_set.end()) continue;
            descriptor_set.insert(cs);

            if (mpos->second.chain_type() == "ATOMN") {
                 if (!na_descriptor.empty()) na_descriptor += ", ";
                 na_descriptor += cs;
            } else {
                 if (!prot_descriptor.empty()) prot_descriptor += ", ";
                 prot_descriptor += cs;
                 if (!prot_ec_descriptor.empty()) prot_ec_descriptor += ", ";
                 prot_ec_descriptor += cs;
                 cs = mpos->second.getValue("pdbx_ec");
                 if (!cs.empty()) prot_ec_descriptor += " (E.C." + cs + ")";
            }
       }

       descriptor.clear();
       if (!na_type.empty()) {
            if (!prot_descriptor.empty())
                 descriptor = prot_descriptor + "/" + na_type + " Complex";
            else descriptor = na_descriptor;
       } else if (!prot_ec_descriptor.empty()) {
            descriptor = prot_ec_descriptor;
       } else {
            ISTable *t = getTablePtr(block, "struct_keywords");
            if (t) get_value_clean(descriptor, t, 0, "pdbx_keywords");
       }

       return descriptor;
}

void ValidateObj::_get_oneletter_code(std::string& code, std::string& code_can)
{
       code_can = "X";
       if (code.empty()) return;

       try {
            const ConnectFormat& drug = _ccDic->find_drug(code);
            std::string one_letter_code = drug.getMetaData("one_letter_code");
            if (!one_letter_code.empty()) code_can = one_letter_code;
       } catch (const std::exception& exc) {}

       std::string buffer = "(" + code + ")";

       if (code.size() == 1) {
            if (code_can == "X") code_can = code;
            if (!SeqCodeUtil::is_rna_residue(code)) code = buffer;
            return;
       }

       if (code_can != "X" && SeqCodeUtil::is_standard_aa_residue(code))
            code = code_can;
       else code = buffer;
}

bool ValidateObj::_check_struct_title(const std::string& title)
{
       if ((title.size() < 20) || (title.size() > 200)) return false;

       std::string cs;
       String::UpperCase(title, cs);
       std::string stop_word_list = _rcsbroot + "/data/ascii/stop_word_list";
       std::set<std::string> stop_word_set = get_stop_words(stop_word_list);
       std::multimap<unsigned int, std::string> stop_word_map;
       stop_word_map.clear();
       for (std::set<std::string>::const_iterator spos = stop_word_set.begin(); spos != stop_word_set.end(); ++spos) {
            stop_word_map.insert(std::make_pair(spos->size(), *spos));
       }
       for (std::multimap<unsigned int, std::string>::const_reverse_iterator mpos = stop_word_map.rbegin(); mpos != stop_word_map.rend(); ++mpos) {
            std::string::size_type pos = cs.find(mpos->second);
            if (pos != std::string::npos) cs.replace(pos, mpos->first, "");
       }
       String::StripAndCompressWs(cs);
       if (cs.size() < 4) return false;
       return true;
}

void ValidateObj::_get_missing_backbone_atoms(RCSB::Residue* res, std::vector<std::string>& atoms)
{
       atoms.clear();
       if (!SeqCodeUtil::is_standard_aa_residue_plus_MSE(res->ResName())) return;
       if (res->missing().empty()) return;
       for (std::vector<std::string>::const_iterator vpos = res->missing().begin(); vpos != res->missing().end(); ++vpos) {
            if (_amino_acid_backbone_atom_set.find((*vpos)) != _amino_acid_backbone_atom_set.end()) atoms.push_back(*vpos);
       }
}

void ValidateObj::_check_missing_extra_atoms(const std::string& mol_ID, RCSB::Residue* res, const std::vector<std::string>& atoms, std::string& message)
{
       if (atoms.empty()) return;

       message += " " + FormattedString(res->ResName(), 7) + "(" + FormattedString(mol_ID, 3) + " " + FormattedString(res->pdb_chnid(), 2, false, true)
                + FormattedString(res->pdb_res_no(), 4) + ")        ";
       int len = 27;
       for (std::vector<std::string>::const_iterator pos = atoms.begin(); pos != atoms.end(); ++pos) {
            len += 5;
            if (len > 81) {
                 message += "\n                           ";
                 len = 27;
            }
            message += FormattedString(*pos, 4) + " ";
       }
       message += "\n";
}

void ValidateObj::_checking_pdbx_database_related(Block& block, const bool& release_flag)
{
       // DAOTHER-5197
       ISTable *t = getTablePtr(block, "pdbx_database_related");
       if (!t) {
            if (!release_flag && (_experiment_type & EXPERIMENT_TYPE_EM))
                 _error_messages.insertMessage("database_print", "error", "No pdbx_database_related category.");
            return;
       }

       std::map<std::string, std::set<std::string> > content_type_value_mapping;
       content_type_value_mapping.clear();

       std::set<std::string> t_set;
       t_set.clear();
       t_set.insert("unspecified");
       t_set.insert("split");
       t_set.insert("re-refinement");
       t_set.insert("complete structure");
       t_set.insert("ensemble");
       t_set.insert("representative structure");
       content_type_value_mapping.insert(std::make_pair("PDB", t_set));

       t_set.clear();
       t_set.insert("associated EM volume");
       t_set.insert("other EM volume");
       t_set.insert("focused EM volume");
       t_set.insert("consensus EM volume");
       content_type_value_mapping.insert(std::make_pair("EMDB", t_set));

       std::map<std::string, std::string> id_pattern_mapping, db_name_enumeration;
       id_pattern_mapping.clear();
       id_pattern_mapping.insert(std::make_pair("PDB", "[1-9][A-Za-z0-9]{3}"));
       id_pattern_mapping.insert(std::make_pair("EMDB", "[Ee][Mm][Dd][-][0-9]{4,5}"));
       id_pattern_mapping.insert(std::make_pair("BMRB", "[1-9][0-9]{3,4}"));

       db_name_enumeration.clear();
       if (_dictUtil.Read()) db_name_enumeration = _dictUtil.getEnumeration("pdbx_database_related", "db_name");
       std::string cs, error;
       std::vector<std::string> items, data;
       items.clear();
       items.push_back("db_name");
       items.push_back("db_id");
       items.push_back("content_type");

       bool found_emdb_id = false;
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            data.clear();
            for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
                 get_value_clean(cs, t, i, *pos);
                 data.push_back(cs);
            }
            cs.clear();
            int count = 0;
            for (unsigned int j = 0; j < 2; ++j) {
                 if (!data[j].empty()) continue;
                 count++;
                 if (!cs.empty()) cs += " & ";
                 cs += "'_pdbx_database_related." + items[j] + "'";
            }
            if (!cs.empty()) {
                 if (count > 1) error = "Missing values in " + cs + " fields at row " + String::IntToString(i + 1) + ".";
                 else error = "Missing value in " + cs + " field at row " + String::IntToString(i + 1) + ".";
                 _error_messages.insertMessage("database_print", "error", error);
                 continue;
            }

            String::UpperCase(data[0], cs);
            if (!db_name_enumeration.empty() && (db_name_enumeration.find(cs) == db_name_enumeration.end())) {
                 error = "Unknown value '" + data[0] + "' in '_pdbx_database_related.db_name' field at row " + String::IntToString(i + 1) + ".";
                 _error_messages.insertMessage("database_print", "error", error);
                 continue;
            }

            std::map<std::string, std::string>::const_iterator mpos = id_pattern_mapping.find(cs);
            if (mpos != id_pattern_mapping.end()) {
                 if (!checking_regular_expression_match_ok(mpos->second, data[1])) {
                      error = "Incorrect "  + data[0] + " ID '" + data[1] + "' in 'pdbx_database_related.db_id' field at row " + String::IntToString(i + 1) + ".";
                      _error_messages.insertMessage("database_print", "error", error);
                 }
            }

            if ((cs == "EMDB") && !data[1].empty() && (data[2] == "associated EM volume")) {
                 found_emdb_id = true;
            }

            if (!release_flag && (cs == "PDB") && !data[1].empty() && (data[2] == "re-refinement")) {
                 error = "This entry is a re-refinement of the original data from " + data[1] + ". Please verify.";
                 _error_messages.insertMessage("re_refinement_print", "warning", error);
            }

            std::map<std::string, std::set<std::string> >::const_iterator cpos = content_type_value_mapping.find(cs);
            if (cpos == content_type_value_mapping.end()) continue;
            
            if (data[2].empty()) {
                 error = "Missing value in '_pdbx_database_related.content_type' field at row " + String::IntToString(i + 1) + ".";
                 _error_messages.insertMessage("database_print", "error", error);
            } else if (cpos->second.find(data[2]) == cpos->second.end()) { 
                 error = "Incorrect '" + data[2] + "' value in '_pdbx_database_related.content_type' field at row " + String::IntToString(i + 1) + ".";
                 _error_messages.insertMessage("database_print", "error", error);
            }
       }
       if (!found_emdb_id && !release_flag && (_experiment_type & EXPERIMENT_TYPE_EM)) {
            _error_messages.insertMessage("database_print", "error", "No related EMDB in pdbx_database_related category.");
       }
}

void ValidateObj::_get_struct_ncs_dom_lim_item_names(ISTable* t, std::vector<std::string>& all_beg_items, std::vector<std::vector<std::string> >& beg_items,
                                                     std::vector<std::string>& all_end_items, std::vector<std::vector<std::string> >& end_items)
{
       std::vector<std::string> data;

       beg_items.clear();
       data.clear();
       data.push_back("beg_label_asym_id");
       data.push_back("beg_label_comp_id");
       data.push_back("beg_label_seq_id");
       beg_items.push_back(data);
       data.clear();
       data.push_back("beg_auth_asym_id");
       if (!t->IsColumnPresent("beg_auth_comp_id") && t->IsColumnPresent("beg_label_comp_id"))
            data.push_back("beg_label_comp_id");
       else data.push_back("beg_auth_comp_id");
       data.push_back("beg_auth_seq_id");
       data.push_back("pdbx_beg_PDB_ins_code"); // not exist in struct_ncs_dom_lim
       beg_items.push_back(data);

       all_beg_items.clear();
       all_beg_items.push_back("beg_label_asym_id");
       all_beg_items.push_back("beg_label_comp_id");
       all_beg_items.push_back("beg_label_seq_id");
       all_beg_items.push_back("beg_auth_asym_id");
       all_beg_items.push_back("beg_auth_comp_id");
       all_beg_items.push_back("beg_auth_seq_id");

       end_items.clear();
       data.clear();
       data.push_back("end_label_asym_id");
       data.push_back("end_label_comp_id");
       data.push_back("end_label_seq_id");
       end_items.push_back(data);
       data.clear();
       data.push_back("end_auth_asym_id");
       if (!t->IsColumnPresent("end_auth_comp_id") && t->IsColumnPresent("end_label_comp_id"))
            data.push_back("end_label_comp_id");
       else data.push_back("end_auth_comp_id");
       data.push_back("end_auth_seq_id");
       data.push_back("pdbx_end_PDB_ins_code"); // not exist in struct_ncs_dom_lim
       end_items.push_back(data);

       all_end_items.clear();
       all_end_items.push_back("end_label_asym_id");
       all_end_items.push_back("end_label_comp_id");
       all_end_items.push_back("end_label_seq_id");
       all_end_items.push_back("end_auth_asym_id");
       all_end_items.push_back("end_auth_comp_id");
       all_end_items.push_back("end_auth_seq_id");
}

RCSB::Residue* ValidateObj::_find_struct_ncs_dom_lim_residue(ISTable* t, const unsigned int& row, const std::vector<std::vector<std::string> >& items,
                                                             const std::vector<std::string>& all_items, std::string& msg)
{
       msg.clear();

       std::vector<std::vector<std::string> > values, tmp_items;
       get_values_clean(t, row, items, values);

       RCSB::Residue* res = NULL;
       if (!_molecules.empty()) {
            if (!values[0][0].empty() && !values[0][1].empty() && !values[0][2].empty()) {
                 res = _molecules[0]->find_orig_cif_residue(values[0][0], values[0][1], values[0][2]);
                 if (!res && !values[1][0].empty() && !values[1][1].empty() && !values[1][2].empty()) {
                      res = _molecules[0]->find_orig_pdb_residue(values[1][0], values[1][1], values[1][2], values[1][3]);
                 }
            } else if (!values[1][0].empty() && !values[1][1].empty() && !values[1][2].empty()) {
                 res = _molecules[0]->find_orig_pdb_residue(values[1][0], values[1][1], values[1][2], values[1][3]);
            } else {
                 res = _molecules[0]->find_orig_cif_residue(values[0][0], values[0][1], values[0][2]);
                 if (!res && !values[1][0].empty() && values[1][1].empty() && !values[1][2].empty()) {
                      res = _molecules[0]->find_orig_pdb_residue(values[1][0], values[1][2], values[1][3]);
                 }
            }
       }

       if (!res) {
            tmp_items.clear();
            tmp_items.push_back(all_items);
            get_values_clean(t, row, tmp_items, values);
 
            for (unsigned int i = 0; i < all_items.size(); ++i) {
                 if (values[0][i].empty()) continue;
                 if (!msg.empty()) msg += ", ";
                 msg += "_" + t->GetName() + "." + all_items[i] + "=" + values[0][i];
            }
       }

       return res;
}
