/*
FILE:     SwitchLabeling.C
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
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include "ChiralVolumeUtil.h"
#include "GetPairList.h"
#include "GetPairList.C"
#include "SeqCodeUtil.h"
#include "SwitchLabelingUtil.h"
#include "ValidateObj.h"
#include "utillib.h"

static bool switch_atom_labeling(ConnectDic* ccDic, const CHIRAL_LABELINGS& labelings, std::vector<RCSB::Residue*>& prev, RCSB::Residue* res);
static bool switch_atom_labeling(const bool& check_linkage, const CHIRAL_LABELINGS& labelings, RCSB::Residue* prev, RCSB::Residue* curr);
static void find_atom_list(RCSB::Residue* res, const std::string& atom_name, const int& isHydrogen, std::vector<RCSB::Atom*>& atom_list);
static bool switch_both_atom_labeling(const bool& check_linkage, const CHIRAL_LABELING* labelings, std::vector<std::vector<RCSB::Atom*> >& atom_lists,
                                      RCSB::Residue* res);
static bool switch_single_atom_labeling(const bool& check_linkage, const CHIRAL_LABELING& labeling, std::vector<std::vector<RCSB::Atom*> >& atom_lists,
                                        RCSB::Residue* res);
static void update_atom_name(RCSB::Atom* atom, const std::string& atom_name, const int& isHydrogen);
static void switch_hydrogen_labeling(const CHIRAL_LABELING& labeling, RCSB::Atom* atom, RCSB::Residue* res);
static void check_atom_labeling(const bool& check_linkage, std::string& error, const std::string& mol_id, const CHIRAL_LABELING& chiral,
                                RCSB::Residue* prev, RCSB::Residue* curr);
static std::string get_atom_name_with_alt_loc(RCSB::Atom* atom);

void ValidateObj::SwitchLabeling()
{
       if (_molecules.empty()) return;

       SwitchLabelingUtil::initialize();

       std::vector<RCSB::Residue*> prev, curr;

       for (std::vector<RCSB::Molecule*>::const_iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            RCSB::Chain* chain = (*mpos)->GetFirstChain();
            while (chain) {
                 if (chain->chain_type() != "ATOMP" && chain->chain_type() != "ATOMN") {
                      chain = (*mpos)->GetNextChain();
                      continue;
                 }

                 prev.clear();
                 chain->GetFirstResidueList(curr);
                 while (!curr.empty()) {
                      for (std::vector<RCSB::Residue*>::iterator rpos = curr.begin(); rpos != curr.end(); ++rpos) {
                           CHIRAL_LABELING_POOL* ext = SwitchLabelingUtil::find_chiral_labeling((*rpos)->ResName());
                           bool non_standard_flag = false;
                           if (!ext && chain->chain_type() == "ATOMN") {
                                non_standard_flag = true;
                                ext = SwitchLabelingUtil::find_chiral_labeling("DA");
                           }
                           if (!ext) continue;

                           bool relabel_flag = false;
                           for (int i = 0; i < ext->nums; ++i) {
                                if (non_standard_flag && ext->lists[i].isPhosphorus == 0) continue;
                                if (switch_atom_labeling(_ccDic, ext->lists[i], prev, *rpos)) relabel_flag = true;
                           }
                           if (relabel_flag) {
                                (*rpos)->UpdateIndices();
                                (*rpos)->reorder_atoms();
                           }
                      }
                      prev = curr;
                      chain->GetNextResidueList(curr);
                 }
                 chain = (*mpos)->GetNextChain();
            }
       }
}

void ValidateObj::CheckLabeling()
{
       if (_molecules.empty()) return;

       SwitchLabelingUtil::initialize();

       std::vector<RCSB::Residue*> prev, curr;

       std::string error = "";

       for (std::vector<RCSB::Molecule*>::const_iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            std::string mol_id = "";
            if (_molecules.size() > 1) mol_id = String::IntToString((*mpos)->Mol_ID());
            RCSB::Chain* chain = (*mpos)->GetFirstChain();
            while (chain) {
                 if (chain->chain_type() != "ATOMP" && chain->chain_type() != "ATOMN") {
                      chain = (*mpos)->GetNextChain();
                      continue;
                 }

                 prev.clear();
                 chain->GetFirstResidueList(curr);
                 while (!curr.empty()) {
                      for (std::vector<RCSB::Residue*>::iterator rpos = curr.begin(); rpos != curr.end(); ++rpos) {
                           CHIRAL_LABELING_POOL* ext = SwitchLabelingUtil::find_chiral_labeling((*rpos)->ResName());
                           if (!ext) continue;

                           for (int i = 0; i < ext->nums; ++i) {
                                for (int j = 0; j < ext->lists[i].num; ++j) {
                                     if (ext->lists[i].isPhosphorus) {
                                          if (!prev.empty()) {
                                               for (std::vector<RCSB::Residue*>::iterator ppos = prev.begin(); ppos != prev.end(); ++ppos) {
                                                    check_atom_labeling(true, error, mol_id, ext->lists[i].chiral[j], *ppos, *rpos);
                                               }
                                          }
                                     } else check_atom_labeling(false, error, mol_id, ext->lists[i].chiral[j], *rpos, *rpos);
                                }
                           }
                      }
                      prev = curr;
                      chain->GetNextResidueList(curr);
                 }
                 chain = (*mpos)->GetNextChain();
            }
       }

       if (error.empty()) return;

       std::string message = "Some of O1P/O2P and/or hydrogen atoms do not follow the convention defined in the\n";
       message += "standard IUBMB nomenclature (Liebecq, C. Compendium of Biochemical Nomenclature and\n";
       message += "Related Documents, 2nd ed.; Portland Press: London and Chapel Hill, 1992). If you do\n";
       message += "not indicate otherwise, we will switch the labels of these atoms as shown below.\n\n";
       if (_molecules.size() > 1) {
            message += "Model  Chain  Residue  Residue  Chiral       Chiral Neighbors      Original\n";
            message += "               Name    Number   Center                             Atom Name\n";
            message += "-----  -----  -------  ------- --------  ------------------------  ---------\n";
       } else {
            message += "Chain  Residue  Residue  Chiral       Chiral Neighbors      Original\n";
            message += "        Name    Number   Center                             Atom Name\n";
            message += "-----  -------  ------- --------  ------------------------  ---------\n";
       }
       message += error;

       _error_messages.insertMessage("label_print", "warning", message);
}

static bool switch_atom_labeling(ConnectDic* ccDic, const CHIRAL_LABELINGS& labelings, std::vector<RCSB::Residue*>& prev, RCSB::Residue* res)
{
       bool relabel_flag = false;

       if (labelings.isPhosphorus) {
            if (!prev.empty() && (SeqCodeUtil::is_na_residue(res->ResName()) || ccDic->hasValidPhosphorylGroup(res->ResName()))) {
                 for (std::vector<RCSB::Residue*>::iterator rpos = prev.begin(); rpos != prev.end(); ++rpos) {
                      if (!(*rpos)->alt_loc().empty() && !res->alt_loc().empty() && (*rpos)->alt_loc() != res->alt_loc()) continue;
                      if (switch_atom_labeling(true, labelings, *rpos, res)) relabel_flag = true;
                 }
            }
       } else if (switch_atom_labeling(false, labelings, res, res)) relabel_flag = true;

       return relabel_flag;
}

static bool switch_atom_labeling(const bool& check_linkage, const CHIRAL_LABELINGS& labelings, RCSB::Residue* prev, RCSB::Residue* curr)
{
       std::vector<std::vector<RCSB::Atom*> > atom_lists;
       std::vector<RCSB::Atom*> first_list, second_list;

       atom_lists.clear();

       curr->find_atom(labelings.chiral[0].center_atom, first_list);
       if (first_list.empty()) return false;
       atom_lists.push_back(first_list);

       curr->find_atom(labelings.chiral[0].first_atom, first_list);
       if (first_list.empty()) return false;
       atom_lists.push_back(first_list);

       prev->find_atom(labelings.chiral[0].second_atom, first_list);
       if (first_list.empty()) return false;
       atom_lists.push_back(first_list);

       find_atom_list(curr, labelings.chiral[0].check_atom, labelings.chiral[0].isHydrogen, first_list);
       find_atom_list(curr, labelings.chiral[1].check_atom, labelings.chiral[1].isHydrogen, second_list);
       if (first_list.empty() && second_list.empty()) return false;

       bool relabel_flag = false;
       if (!first_list.empty() && !second_list.empty()) {
            atom_lists.push_back(first_list);
            atom_lists.push_back(second_list);
            relabel_flag = switch_both_atom_labeling(check_linkage, labelings.chiral, atom_lists, curr);
       } else if (!first_list.empty()) {
            atom_lists.push_back(first_list);
            relabel_flag = switch_single_atom_labeling(check_linkage, labelings.chiral[0], atom_lists, curr);
       } else if (!second_list.empty()) {
            atom_lists.push_back(second_list);
            relabel_flag = switch_single_atom_labeling(check_linkage, labelings.chiral[1], atom_lists, curr);
       }

       return relabel_flag;
}

static void find_atom_list(RCSB::Residue* res, const std::string& atom_name, const int& isHydrogen, std::vector<RCSB::Atom*>& atom_list)
{
       atom_list.clear();
       res->find_atom(atom_name, atom_list);
       if (isHydrogen) {
            std::string name = "D" + atom_name.substr(1);
            std::vector<RCSB::Atom*> d_atom_list;
            res->find_atom(name, d_atom_list);
            for (std::vector<RCSB::Atom*>::const_iterator vpos = d_atom_list.begin(); vpos != d_atom_list.end(); ++vpos) {
                 atom_list.push_back(*vpos);
            }
       }
}

static bool switch_both_atom_labeling(const bool& check_linkage, const CHIRAL_LABELING* labelings, std::vector<std::vector<RCSB::Atom*> >& atom_lists,
                                      RCSB::Residue* res)
{
       bool relabel_flag = false;

       std::vector<std::vector<RCSB::Atom*> > pair_lists;
       std::vector<RCSB::Atom*> atom_list0, atom_list1;

       GetPairList(atom_lists, pair_lists);
       if (pair_lists.empty()) return relabel_flag;

       double ang;
       if ((atom_lists[0].size() == 1) && atom_lists[0][0]->alt_loc().empty() && (atom_lists[1].size() == 1) && atom_lists[1][0]->alt_loc().empty() &&
           (atom_lists[2].size() == 1) && atom_lists[2][0]->alt_loc().empty() && (atom_lists[3].size() == 1) && atom_lists[3][0]->alt_loc().empty() &&
           (atom_lists[4].size() > 1)) {
            if (check_linkage) {
                 double dist = cal_distance(atom_lists[0][0], atom_lists[2][0]);
                 if (dist > 2.0) return relabel_flag;
            }

            atom_list0.clear();
            atom_list0.push_back(atom_lists[0][0]);
            atom_list0.push_back(atom_lists[1][0]);
            atom_list0.push_back(atom_lists[2][0]);
            atom_list0.push_back(atom_lists[3][0]);
            double value0 = calc_chiral_volume(atom_list0, ang);

            unsigned int switch_count = 0;
            for (std::vector<RCSB::Atom*>::const_iterator apos = atom_lists[4].begin(); apos != atom_lists[4].end(); ++apos) {
                 atom_list1.clear();
                 atom_list1.push_back(atom_lists[0][0]);
                 atom_list1.push_back(atom_lists[1][0]);
                 if ((res->ResName() == "LEU" || res->ResName() == "VAL") && !labelings[0].isHydrogen)
                      atom_list1.push_back(atom_lists[3][0]);
                 else atom_list1.push_back(atom_lists[2][0]);
                 atom_list1.push_back(*apos);
                 double value1 = calc_chiral_volume(atom_list1, ang);

                 if (value0 * labelings[0].correct_value >= 0 || value1 * labelings[1].correct_value >= 0) continue;
                 switch_count++;
            }

            if (switch_count == atom_lists[4].size()) {
                 relabel_flag = true;
                 update_atom_name(atom_lists[3][0], labelings[0].change_to, labelings[0].isHydrogen);
                 switch_hydrogen_labeling(labelings[0], atom_lists[3][0], res);
                 for (std::vector<RCSB::Atom*>::iterator apos = atom_lists[4].begin(); apos != atom_lists[4].end(); ++apos) {
                      update_atom_name(*apos, labelings[1].change_to, labelings[1].isHydrogen);
                      switch_hydrogen_labeling(labelings[1], *apos, res);
                 }
            }
       } else if ((atom_lists[0].size() == 1) && atom_lists[0][0]->alt_loc().empty() && (atom_lists[1].size() == 1) && atom_lists[1][0]->alt_loc().empty() &&
                  (atom_lists[2].size() == 1) && atom_lists[2][0]->alt_loc().empty() && (atom_lists[4].size() == 1) && atom_lists[4][0]->alt_loc().empty() &&
                  (atom_lists[3].size() > 1)) {
            if (check_linkage) {
                 double dist = cal_distance(atom_lists[0][0], atom_lists[2][0]);
                 if (dist > 2.0) return relabel_flag;
            }

            unsigned int switch_count = 0;
            for (std::vector<RCSB::Atom*>::const_iterator apos = atom_lists[3].begin(); apos != atom_lists[3].end(); ++apos) {
                 atom_list0.clear();
                 atom_list0.push_back(atom_lists[0][0]);
                 atom_list0.push_back(atom_lists[1][0]);
                 atom_list0.push_back(atom_lists[2][0]);
                 atom_list0.push_back(*apos);
                 double value0 = calc_chiral_volume(atom_list0, ang);

                 atom_list1.clear();
                 atom_list1.push_back(atom_lists[0][0]);
                 atom_list1.push_back(atom_lists[1][0]);
                 if ((res->ResName() == "LEU" || res->ResName() == "VAL") && !labelings[0].isHydrogen)
                      atom_list1.push_back(*apos);
                 else atom_list1.push_back(atom_lists[2][0]);
                 atom_list1.push_back(atom_lists[4][0]);
                 double value1 = calc_chiral_volume(atom_list1, ang);

                 if (value0 * labelings[0].correct_value >= 0 || value1 * labelings[1].correct_value >= 0) continue;
                 switch_count++;
            }

            if (switch_count == atom_lists[3].size()) {
                 relabel_flag = true;
                 for (std::vector<RCSB::Atom*>::iterator apos = atom_lists[3].begin(); apos != atom_lists[3].end(); ++apos) {
                      update_atom_name(*apos, labelings[0].change_to, labelings[0].isHydrogen);
                      switch_hydrogen_labeling(labelings[0], *apos, res);
                 }
                 update_atom_name(atom_lists[4][0], labelings[1].change_to, labelings[1].isHydrogen);
                 switch_hydrogen_labeling(labelings[1], atom_lists[4][0], res);
            }
       } else {
            for (std::vector<std::vector<RCSB::Atom*> >::iterator pos = pair_lists.begin(); pos != pair_lists.end(); ++pos) {
                 if (check_linkage) {
                      double dist = cal_distance((*pos)[0], (*pos)[2]);
                      if (dist > 2.0) continue;
                 }
                 atom_list0.clear();
                 atom_list0.push_back((*pos)[0]);
                 atom_list0.push_back((*pos)[1]);
                 atom_list0.push_back((*pos)[2]);
                 atom_list0.push_back((*pos)[3]);
                 double value0 = calc_chiral_volume(atom_list0, ang);

                 atom_list1.clear();
                 atom_list1.push_back((*pos)[0]);
                 atom_list1.push_back((*pos)[1]);
                 if ((res->ResName() == "LEU" || res->ResName() == "VAL") && !labelings[0].isHydrogen)
                      atom_list1.push_back((*pos)[3]);
                 else atom_list1.push_back((*pos)[2]);
                 atom_list1.push_back((*pos)[4]);
                 double value1 = calc_chiral_volume(atom_list1, ang);

                 if (value0 * labelings[0].correct_value >= 0 || value1 * labelings[1].correct_value >= 0) continue;

                 if ((*pos)[3]->alt_loc() != (*pos)[4]->alt_loc()) continue;

                 relabel_flag = true;
                 update_atom_name((*pos)[3], labelings[0].change_to, labelings[0].isHydrogen);
                 switch_hydrogen_labeling(labelings[0], (*pos)[3], res);
                 update_atom_name((*pos)[4], labelings[1].change_to, labelings[1].isHydrogen);
                 switch_hydrogen_labeling(labelings[1], (*pos)[4], res);
            }
       }

       return relabel_flag;
}

static bool switch_single_atom_labeling(const bool& check_linkage, const CHIRAL_LABELING& labeling, std::vector<std::vector<RCSB::Atom*> >& atom_lists,
                                        RCSB::Residue* res)
{
       std::vector<std::vector<RCSB::Atom*> > pair_lists;

       GetPairList(atom_lists, pair_lists);
       if (pair_lists.empty()) return false;

       bool relabel_flag = false;
       double ang;
       for (std::vector<std::vector<RCSB::Atom*> >::iterator pos = pair_lists.begin(); pos != pair_lists.end(); ++pos) {
            if (check_linkage) {
                 double dist = cal_distance((*pos)[0], (*pos)[2]);
                 if (dist > 2.0) continue;
            }

            double value_ref = calc_chiral_volume(*pos, ang);
            if (value_ref * labeling.correct_value >= 0) continue;

            relabel_flag = true;
            update_atom_name((*pos)[3], labeling.change_to, labeling.isHydrogen);
            switch_hydrogen_labeling(labeling, (*pos)[3], res);
       }

       return relabel_flag;
}

static void update_atom_name(RCSB::Atom* atom, const std::string& atom_name, const int& isHydrogen)
{
       if (isHydrogen && (atom->atom_type() == "D")) {
            atom->set_atmtype("D" + atom_name.substr(1));
            atom->set_pdb_atmnam("D" + atom_name.substr(1));
       } else {
            atom->set_atmtype(atom_name);
            atom->set_pdb_atmnam(atom_name);
       }
}

static void switch_hydrogen_labeling(const CHIRAL_LABELING& labeling, RCSB::Atom* atom, RCSB::Residue* res)
{
       if (labeling.n_change == 0) return;

       std::vector<RCSB::Atom*> atom_list;
       for (int i = 0; i < labeling.n_change; ++i) {
            find_atom_list(res, labeling.change_list[i][0], 1, atom_list);
            if (atom_list.empty()) continue;
            for (unsigned int j = 0; j < atom_list.size(); ++j) {
                 if (!atom->alt_loc().empty() && !atom_list[j]->alt_loc().empty() && (atom->alt_loc() != atom_list[j]->alt_loc())) continue;
                 update_atom_name(atom_list[j], labeling.change_list[i][1], 1);
            }
       }
}

static void check_atom_labeling(const bool& check_linkage, std::string& error, const std::string& mol_id, const CHIRAL_LABELING& chiral,
                                RCSB::Residue* prev, RCSB::Residue* curr)
{
       std::vector<std::vector<RCSB::Atom*> > atom_lists, pair_lists;
       std::vector<RCSB::Atom*> atom_list;

       atom_lists.clear();

       curr->find_atom(chiral.center_atom, atom_list);
       if (atom_list.empty()) return;
       atom_lists.push_back(atom_list);

       curr->find_atom(chiral.first_atom, atom_list);
       if (atom_list.empty()) return;
       atom_lists.push_back(atom_list);

       prev->find_atom(chiral.second_atom, atom_list);
       if (atom_list.empty()) return;
       atom_lists.push_back(atom_list);

       find_atom_list(curr, chiral.check_atom, chiral.isHydrogen, atom_list);
       if (atom_list.empty()) return;
       atom_lists.push_back(atom_list);

       GetPairList(atom_lists, pair_lists);
       if (pair_lists.empty()) return;

       double ang;
       std::vector<std::string> names;

       for (std::vector<std::vector<RCSB::Atom*> >::const_iterator pos = pair_lists.begin(); pos != pair_lists.end(); ++pos) {
            if (check_linkage) {
                 double dist = cal_distance((*pos)[0], (*pos)[2]);
                 if (dist > 2.0) continue;
            }

            double value_ref = calc_chiral_volume(*pos, ang);
            if (value_ref * chiral.correct_value >= 0) continue;

            names.clear();
            for (std::vector<RCSB::Atom*>::const_iterator apos = pos->begin(); apos != pos->end(); ++apos) {
                 names.push_back(get_atom_name_with_alt_loc(*apos));
            }
            if (!mol_id.empty()) error += FormattedString(mol_id, 3, false, false) + "    ";
            error += (*pos)[0]->pdb_chnid() + "      " + FormattedString((*pos)[0]->pdb_resnam(), 4, false, false) + "    "
                   + FormattedString((*pos)[0]->pdb_resnum(), 5, false, false) + "    " + FormattedString(names[0], 8, true, false) + "  "
                   + FormattedString(names[3], 24, true, false) + "  ";
            if ((*pos)[3]->pdb_atmnam() != (*pos)[3]->pub_atmnam())
                 error += names[3] + "\n";
            else error += "\n";
       }
}

static std::string get_atom_name_with_alt_loc(RCSB::Atom* atom)
{
       std::string name = atom->pdb_atmnam();
       if (!atom->alt_loc().empty()) name += "(" + atom->alt_loc() + ")";
       return name;
}
