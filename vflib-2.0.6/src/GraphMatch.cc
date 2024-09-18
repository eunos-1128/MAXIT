/*
FILE:     GraphMatch.cc
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

/*!
** \file GraphMatch.C
**
** \brief Implementation file for GraphMatch class.
*/
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <stdexcept>
#include <ctype.h>

#include "GraphMatch.h"
#include "argedit.h"
#include "vf2_sub_state.h"
#include "match.h"

typedef struct {
       int order;
       std::vector<std::string> target;
       std::vector<std::string> reference;
       std::vector<std::vector<std::pair<std::string, std::string> > > matchlists;      
} _CallBackMatch;

static bool matchCallBack(int n,  node_id ni1[], node_id ni2[], void *usr_data)
{
       _CallBackMatch *cbm = (_CallBackMatch *) usr_data;

       std::vector<std::pair<std::string, std::string> > matchlist;
       matchlist.clear();
       for (int i = 0; i < n; ++i) {
            if (cbm->order == TARGET_TO_REF)
                 matchlist.push_back(std::make_pair(cbm->target[ni1[i]],
                                               cbm->reference[ni2[i]]));
            else matchlist.push_back(std::make_pair(cbm->reference[ni2[i]],
                                                    cbm->target[ni1[i]]));
       }
       cbm->matchlists.push_back(matchlist);

       return false;
}

void GraphMatch::VF_GetMatch(const std::vector<std::pair<std::string, std::string> >& ref_atomlist,
                             const std::vector<std::pair<std::string, std::string> >& ref_bondlist,
                             const std::vector<std::pair<std::string, std::string> >& tgt_atomlist,
                             const std::vector<std::pair<std::string, std::string> >& tgt_bondlist,
                             std::map<std::string, std::string>& matchlist, const int& order)
{
       matchlist.clear();

       std::vector<std::vector<std::pair<std::string, std::string> > > match_pair_lists;
       DoMatch(ref_atomlist, ref_bondlist, tgt_atomlist, tgt_bondlist, match_pair_lists, order);
       if (match_pair_lists.empty()) return;

       std::map<std::string, std::set<std::string> > ref_linked_H_atoms, tgt_linked_H_atoms;
       ref_linked_H_atoms.clear();
       tgt_linked_H_atoms.clear();
       GetBestMatch(matchlist, ref_linked_H_atoms, tgt_linked_H_atoms, match_pair_lists);
}

void GraphMatch::GetMatch(const std::vector<std::pair<std::string, std::string> >& ref_atomlist,
                          const std::vector<std::pair<std::string, std::string> >& ref_bondlist,
                          const std::vector<std::pair<std::string, std::string> >& tgt_atomlist,
                          const std::vector<std::pair<std::string, std::string> >& tgt_bondlist,
                          std::map<std::string, std::string>& matchlist, bool& is_substructure_match,
                          const int& order, const bool& heavy_atom_only, const bool& exclude_OP3_flag)
{
       matchlist.clear();
       is_substructure_match = false;

       std::vector<std::pair<std::string, std::string> > ref_atom_no_H;
       std::vector<std::pair<std::string, std::string> > ref_bond_no_H;
       removeHydrogenAtoms(ref_atomlist, ref_atom_no_H, ref_bondlist, ref_bond_no_H);

       std::vector<std::pair<std::string, std::string> > tgt_atom_no_H;
       std::vector<std::pair<std::string, std::string> > tgt_bond_no_H;
       removeHydrogenAtoms(tgt_atomlist, tgt_atom_no_H, tgt_bondlist, tgt_bond_no_H);
/*
       if (tgt_atom_no_H.size() == tgt_atomlist.size()) {
            VF_GetMatch(ref_atom_no_H, ref_bond_no_H, tgt_atom_no_H, tgt_bond_no_H,
                        matchlist, order);
            if (matchlist.empty()) return;
            if (matchlist.size() < ref_atom_no_H.size()) is_substructure_match = true;
       } else {
*/
            std::vector<std::vector<std::pair<std::string, std::string> > > match_pair_lists;
            DoMatch(ref_atom_no_H, ref_bond_no_H, tgt_atom_no_H, tgt_bond_no_H, match_pair_lists, REF_TO_TARGET);
            if (match_pair_lists.empty()) return;

            std::map<std::string, std::set<std::string> > ref_linked_H_atoms;
            getLinkedHydrogenAtoms(ref_linked_H_atoms, ref_atomlist, ref_bondlist);

            std::map<std::string, std::set<std::string> > tgt_linked_H_atoms;
            getLinkedHydrogenAtoms(tgt_linked_H_atoms, tgt_atomlist, tgt_bondlist);

            GetBestMatch(matchlist, ref_linked_H_atoms, tgt_linked_H_atoms, match_pair_lists, exclude_OP3_flag);
            if (matchlist.size() < ref_atom_no_H.size()) is_substructure_match = true;

            if (!heavy_atom_only) AddHydrogenMapping(matchlist, ref_linked_H_atoms, tgt_linked_H_atoms);

            if (order == TARGET_TO_REF) {
                 std::map<std::string, std::string> reverse_mapping_list;
                 reverse_mapping_list.clear();
                 for (std::map<std::string, std::string>::const_iterator
                      pos = matchlist.begin(); pos != matchlist.end(); ++pos) {
                      reverse_mapping_list.insert(std::make_pair(pos->second, pos->first));
                 }
                 matchlist = reverse_mapping_list;
            }
/*
       }
*/
}

ARGraph<_Node_*, _Edge_*>* GraphMatch::getGraph(const std::vector<std::pair<std::string, std::string> >& atomlist,
                                              const std::vector<std::pair<std::string, std::string> >& bondlist)
{
       if (atomlist.size() < 2) return NULL;
       if (bondlist.empty()) return NULL;

       ARGEdit ed;
       std::map<std::string, unsigned int> atomindex;
       atomindex.clear();
       for (unsigned int i = 0; i < atomlist.size(); ++i) {
            ed.InsertNode(new _Node_(atomlist[i].first, atomlist[i].second));
            atomindex.insert(std::make_pair(atomlist[i].second, i));
       }
       for (unsigned int i = 0; i < bondlist.size(); ++i) {
            std::map<std::string, unsigned int>::iterator
                ptr1 = atomindex.find(bondlist[i].first);
            if (ptr1 == atomindex.end()) continue;
            std::map<std::string, unsigned int>::iterator
                ptr2 = atomindex.find(bondlist[i].second);
            if (ptr2 == atomindex.end()) continue;
            ed.InsertEdge(ptr1->second, ptr2->second, new _Edge_(1));
            ed.InsertEdge(ptr2->second, ptr1->second, new _Edge_(1));
       }

       ARGraph<_Node_*, _Edge_*> *graph = new ARGraph<_Node_*, _Edge_*>(&ed);
       graph->SetNodeComparator(new _Node_Comparator());
       graph->SetEdgeComparator(new _Edge_Comparator());
       graph->SetNodeDestroyer(new _Node_Destroyer());
       graph->SetEdgeDestroyer(new _Edge_Destroyer());
       return graph;
}

void GraphMatch::DoMatch(const std::vector<std::pair<std::string, std::string> >& ref_atomlist,
                         const std::vector<std::pair<std::string, std::string> >& ref_bondlist,
                         const std::vector<std::pair<std::string, std::string> >& tgt_atomlist,
                         const std::vector<std::pair<std::string, std::string> >& tgt_bondlist,
                         std::vector<std::vector<std::pair<std::string, std::string> > >& 
                         match_pair_lists, const int& order)
{
       match_pair_lists.clear();

       if (ref_atomlist.size() == 1 && tgt_atomlist.size() == 1) {
            if (ref_atomlist[0].first == tgt_atomlist[0].first) {
                 std::vector<std::pair<std::string, std::string> > tmp_vec;
                 tmp_vec.clear();
                 if (order == TARGET_TO_REF)
                      tmp_vec.push_back(std::make_pair(tgt_atomlist[0].second, ref_atomlist[0].second));
                 else tmp_vec.push_back(std::make_pair(ref_atomlist[0].second, tgt_atomlist[0].second));
                 match_pair_lists.push_back(tmp_vec);
            }
            return;
       }

       try {
            _CallBackMatch cbm;
            cbm.order = order;
            cbm.target.clear();
            cbm.reference.clear();
            cbm.matchlists.clear();

            ARGraph<_Node_*, _Edge_*> *graph1 = getGraph(ref_atomlist, ref_bondlist);
            if (graph1 == NULL) return;
            ARGraph<_Node_*, _Edge_*> *graph2 = getGraph(tgt_atomlist, tgt_bondlist);
            if (graph2 == NULL) {
                 delete graph1;
                 return;
            }

            for (std::vector<std::pair<std::string, std::string> >::const_iterator
                 pos = ref_atomlist.begin(); pos != ref_atomlist.end(); ++pos) {
                 cbm.reference.push_back(pos->second);
            }

            for (std::vector<std::pair<std::string, std::string> >::const_iterator
                 pos = tgt_atomlist.begin(); pos != tgt_atomlist.end(); ++pos) {
                 cbm.target.push_back(pos->second);
            }

            VF2SubState s0(graph2, graph1);
            int nMatch = match(&s0, matchCallBack, &cbm);

            delete graph1;
            delete graph2;

            if (nMatch == 0 || cbm.matchlists.empty()) return;

            match_pair_lists = cbm.matchlists;
       } catch (const std::exception& exc) {}
}

void GraphMatch::GetBestMatch(std::map<std::string, std::string>& matchlist,
                     const std::map<std::string, std::set<std::string> >& linked_H_atoms1,
                     const std::map<std::string, std::set<std::string> >& linked_H_atoms2,
                     const std::vector<std::vector<std::pair<std::string, std::string> > >&
                     match_pair_lists, const bool& exclude_OP3_flag)
{
       int index = 0;
       if (match_pair_lists.size() > 1)
            index = GetBestMatchIndex(linked_H_atoms1, linked_H_atoms2, match_pair_lists, exclude_OP3_flag);

       for (std::vector<std::pair<std::string, std::string> >::const_iterator pos =
            match_pair_lists[index].begin(); pos != match_pair_lists[index].end(); ++pos) {
            matchlist.insert(std::make_pair(pos->first, pos->second));
       }
}

int GraphMatch::GetBestMatchIndex(const std::map<std::string, std::set<std::string> >& linked_H_atoms1,
                                  const std::map<std::string, std::set<std::string> >& linked_H_atoms2,
                                  const std::vector<std::vector<std::pair<std::string, std::string> > >& matchlists,
                                  const bool& exclude_OP3_flag)
{
       int index = 0;
       int max_score = -10;
#if 0
       for (std::vector<std::pair<std::string, std::string> >::const_iterator pos = matchlists[0].begin(); pos != matchlists[0].end(); ++pos) {
            if (exclude_OP3_flag && pos->first == "OP3") {}
            else max_score += getAtomNameScore(pos->first, pos->second);
/*
            if (pos->first == pos->second) max_score += 2;
            else if (isCloseName(pos->first, pos->second))
                 max_score += 1;
*/
            if (linked_H_atoms2.empty() && exclude_OP3_flag && pos->first != "OP3" && pos->second == "OP3") max_score += 1;
            if (linked_H_atoms1.empty() || linked_H_atoms2.empty()) continue;

            std::map<std::string, std::set<std::string> >::const_iterator mpos1 = linked_H_atoms1.find(pos->first);
            std::map<std::string, std::set<std::string> >::const_iterator mpos2 = linked_H_atoms2.find(pos->second);
            if (mpos1 != linked_H_atoms1.end() && mpos2 != linked_H_atoms2.end()) {
                 int s = abs((int) mpos1->second.size() - (int) mpos2->second.size());
                 max_score += 5 - s;
                 if (exclude_OP3_flag && pos->first != "OP3" && pos->second == "OP3") max_score += 1;
            } else if (mpos1 == linked_H_atoms1.end() && mpos2 == linked_H_atoms2.end()) {
                 max_score += 1;
                 if (exclude_OP3_flag && pos->first != "OP3" && pos->second == "OP3") max_score += 1;
            }
       }

       for (unsigned int i = 1; i < matchlists.size(); ++i) {
#endif
       for (unsigned int i = 0; i < matchlists.size(); ++i) {
            int score = 0;
            for (std::vector<std::pair<std::string, std::string> >::const_iterator pos = matchlists[i].begin(); pos != matchlists[i].end(); ++pos) {
                 if (exclude_OP3_flag && pos->first == "OP3") {}
                 else score += getAtomNameScore(pos->first, pos->second);
/*
                 if (pos->first == pos->second) score += 2;
                 else if (isCloseName(pos->first, pos->second))
                      score += 1;
                 if (linked_H_atoms1.empty() || linked_H_atoms2.empty()) continue;
*/

                 std::map<std::string, std::set<std::string> >::const_iterator mpos1 = linked_H_atoms1.find(pos->first);
                 std::map<std::string, std::set<std::string> >::const_iterator mpos2 = linked_H_atoms2.find(pos->second);
                 if (mpos1 != linked_H_atoms1.end() && mpos2 != linked_H_atoms2.end()) {
                      int s = abs((int) mpos1->second.size() - (int) mpos2->second.size());
                      score += 5 - s;
                      if (exclude_OP3_flag && pos->first != "OP3" && pos->second == "OP3") score += 1;
                 } else if (mpos1 == linked_H_atoms1.end() && mpos2 == linked_H_atoms2.end()) {
                      score += 1;
                      if (exclude_OP3_flag && pos->first != "OP3" && pos->second == "OP3") score += 1;
                 } else if (mpos1 != linked_H_atoms1.end() && mpos2 == linked_H_atoms2.end()) {
                      if (exclude_OP3_flag && pos->first != "OP3" && pos->second == "OP3") score += 1;
                 } else if (mpos1 == linked_H_atoms1.end() && mpos2 != linked_H_atoms2.end()) {
                      score -= 5;
                 }
            }
            if (score > max_score) {
                 max_score = score;
                 index = i;
            }
       }

       return index;
}

void GraphMatch::reorderName(std::string& name)
{
       std::set<std::string> charlist;
       charlist.clear();
       for (unsigned int i = 0; i < name.size(); ++i) {
            charlist.insert(name.substr(i, 1));
       }
       name.clear();
       for (std::set<std::string>::iterator
            pos = charlist.begin(); pos != charlist.end(); ++pos) {
            name += *pos;
       }
}

bool GraphMatch::isCloseName(const std::string& name1, const std::string& name2)
{
       if ((name1 == "H" && (name2 == "H1" || name2 == "2H")) ||
           (name2 == "H" && (name1 == "H1" || name1 == "2H")))
            return true;

       std::string _name1 = name1;
       reorderName(_name1);
       std::string _name2 = name2;
       reorderName(_name2);
       if (_name1 == _name2) return true;

       return false;
}

int GraphMatch::getAtomNameScore(const std::string& name1, const std::string& name2)
{
       if (name1 == name2) return 2;

       std::string _name2 = name2;
       if (name1.substr(0, 1) == "H") {
            for (unsigned int i = 0; i < _name2.size(); ++i) {
                 if (isalpha(_name2[i])) {
                      if (_name2[i] == 'D') _name2[i] = 'H';
                      break;
                 }
            }
       }

       if (name1 == _name2) return 2;

       if (isCloseName(name1, _name2)) return 1;

       return 0;
}

void GraphMatch::removeHydrogenAtoms(const std::vector<std::pair<std::string, std::string> >& alist,
                                     std::vector<std::pair<std::string, std::string> >& alist_no_H,
                                     const std::vector<std::pair<std::string, std::string> >& blist,
                                     std::vector<std::pair<std::string, std::string> >& blist_no_H)
{
       alist_no_H.clear();
       blist_no_H.clear();

       std::set<std::string> atom_set;
       atom_set.clear();
       for (std::vector<std::pair<std::string, std::string> >::const_iterator
            pos = alist.begin(); pos != alist.end(); ++pos) {
            if (pos->first == "H" || pos->first == "D") continue;
            atom_set.insert(pos->second);
            alist_no_H.push_back(*pos);
       }
       if (alist_no_H.size() == alist.size()) {
            blist_no_H = blist; 
            return;
       }

       for (std::vector<std::pair<std::string, std::string> >::const_iterator
            pos = blist.begin(); pos != blist.end(); ++pos) {
            if (atom_set.find(pos->first) == atom_set.end()) continue;
            if (atom_set.find(pos->second) == atom_set.end()) continue;
            blist_no_H.push_back(*pos);
       }
}

void GraphMatch::AddHydrogenMapping(std::map<std::string, std::string>& atom_mapping_list,
                       const std::map<std::string, std::set<std::string> >& linked_H_atoms1,
                       const std::map<std::string, std::set<std::string> >& linked_H_atoms2)
{
       std::vector<std::pair<std::string, std::string> > atom_pair, h_atom_pair;
       atom_pair.clear();
       for (std::map<std::string, std::string>::iterator pos = atom_mapping_list.begin(); pos != atom_mapping_list.end(); ++pos) {
            atom_pair.push_back(std::make_pair(pos->first, pos->second));
       }

       for (std::vector<std::pair<std::string, std::string> >::const_iterator pos = atom_pair.begin(); pos != atom_pair.end(); ++pos) {
            std::map<std::string, std::set<std::string> >::const_iterator mpos1 = linked_H_atoms1.find(pos->first);
            if (mpos1 == linked_H_atoms1.end()) continue;
            std::map<std::string, std::set<std::string> >::const_iterator mpos2 = linked_H_atoms2.find(pos->second);
            if (mpos2 == linked_H_atoms2.end()) continue;
            getHydrogenAtomPair(h_atom_pair, mpos1->second, mpos2->second);
            for (std::vector<std::pair<std::string, std::string> >::const_iterator ptr = h_atom_pair.begin(); ptr != h_atom_pair.end(); ++ptr) {
                 atom_mapping_list.insert(std::make_pair(ptr->first, ptr->second));
            }
       }
}
 
void GraphMatch::getLinkedHydrogenAtoms(std::map<std::string, std::set<std::string> >& atoms,
                             const std::vector<std::pair<std::string, std::string> >& a_list,
                             const std::vector<std::pair<std::string, std::string> >& b_list)
{
       atoms.clear();

       std::map<std::string, unsigned int> index;
       index.clear();
       unsigned int count = 0;
       for (std::vector<std::pair<std::string, std::string> >::const_iterator
            pos = a_list.begin(); pos != a_list.end(); ++pos) {
            index.insert(std::make_pair(pos->second, count));
            count++;
       }

       std::set<std::string> data;
       for (std::vector<std::pair<std::string, std::string> >::const_iterator
            pos = b_list.begin(); pos != b_list.end(); ++pos) {
            std::map<std::string, unsigned int>::const_iterator
                ipos1 = index.find(pos->first);
            if (ipos1 == index.end()) continue;
            std::map<std::string, unsigned int>::const_iterator
                ipos2 = index.find(pos->second);
            if (ipos2 == index.end()) continue;
            if (a_list[ipos1->second].first != "H" &&
                a_list[ipos1->second].first != "D" &&
                a_list[ipos2->second].first != "H" &&
                a_list[ipos2->second].first != "D") continue;

            unsigned int n_index = ipos1->second;
            unsigned int h_index = ipos2->second;
            if (a_list[ipos2->second].first != "H" &&
                a_list[ipos2->second].first != "D") {
                 n_index = ipos2->second;
                 h_index = ipos1->second;
            }
            std::map<std::string, std::set<std::string> >::iterator
                apos = atoms.find(a_list[n_index].second);
            if (apos == atoms.end()) {
                 data.clear();
                 data.insert(a_list[h_index].second);
                 atoms.insert(std::make_pair(a_list[n_index].second, data));
            } else apos->second.insert(a_list[h_index].second);
       }
}

void GraphMatch::getHydrogenAtomPair(std::vector<std::pair<std::string, std::string> >& atom_pair,
                             const std::set<std::string>& set1, const std::set<std::string>& set2)
{
       atom_pair.clear();
       if (set1.empty() || set2.empty()) return;

       std::set<std::string> used_atom_set1, used_atom_set2;
       used_atom_set1.clear();
       used_atom_set2.clear();
       for (std::set<std::string>::const_iterator spos = set1.begin(); spos != set1.end(); ++spos) {
            if (set2.find(*spos) != set2.end()) {
                 atom_pair.push_back(std::make_pair(*spos, *spos));
                 used_atom_set1.insert(*spos);
                 used_atom_set2.insert(*spos);
            } else {
                 std::string cs = *spos;
                 if (cs[0] == 'D') cs[0] = 'H';
                 else if (cs[0] == 'H') cs[0] = 'D';
                 if (set2.find(cs) != set2.end()) {
                      atom_pair.push_back(std::make_pair(*spos, cs));
                      used_atom_set1.insert(*spos);
                      used_atom_set2.insert(cs);
                 }
            }
       }

       if ((set1.size() == used_atom_set1.size()) || (set2.size() == used_atom_set2.size())) return;

       std::vector<std::string> list1, list2;
       if (getHydrogenList(set1, list1) && getHydrogenList(set2, list2)) {
            atom_pair.clear();
       } else {
            // convert set into vactor
            list1.clear();
            for (std::set<std::string>::const_iterator spos = set1.begin(); spos != set1.end(); ++spos) {
                 if (used_atom_set1.find(*spos) != used_atom_set1.end()) continue;
                 list1.push_back(*spos);
            }
            if (list1.empty()) return;

            list2.clear();
            for (std::set<std::string>::const_iterator spos = set2.begin(); spos != set2.end(); ++spos) {
                 if (used_atom_set2.find(*spos) != used_atom_set2.end()) continue;
                 list2.push_back(*spos);
            }
            if (list2.empty()) return;
       }

       if (list1.size() == list2.size()) {
            for (unsigned int i = 0; i < list1.size(); ++i) {
                 atom_pair.push_back(std::make_pair(list1[i], list2[i]));
            }
            return;
       }

       std::vector<std::vector<std::pair<std::string, std::string> > > matchlists;
       matchlists.clear();
       std::vector<std::pair<std::string, std::string> > pairlist;

       if (list1.size() > list2.size()) {
            for (unsigned int i = 0; i < (list1.size() - list2.size()); ++i) {
                 pairlist.clear();
                 for (unsigned int j = 0; j < list2.size(); ++j) {
                      pairlist.push_back(std::make_pair(list1[i+j], list2[j]));
                 }
                 matchlists.push_back(pairlist);
            }
       } else {
            for (unsigned int i = 0; i < (list2.size() - list1.size()); ++i) {
                 pairlist.clear();
                 for (unsigned int j = 0; j < list1.size(); ++j) {
                      pairlist.push_back(std::make_pair(list1[j], list2[i+j]));
                 }
                 matchlists.push_back(pairlist);
            }
       }

       std::map<std::string, std::set<std::string> > a_atoms, b_atoms;
       a_atoms.clear();
       b_atoms.clear();
       int index = GetBestMatchIndex(a_atoms, b_atoms, matchlists);
       for (unsigned int i = 0; i < matchlists[index].size(); ++i) {
            atom_pair.push_back(matchlists[index][i]);
       }
}

bool GraphMatch::getHydrogenList(const std::set<std::string>& atom_set, std::vector<std::string>& atom_list)
{
       atom_list.clear();

       std::map<int, std::string> num_name_mapping;
       num_name_mapping.clear();
       for (std::set<std::string>::const_iterator spos = atom_set.begin(); spos != atom_set.end(); ++spos) {
            std::string num = "";
            for (unsigned int i = 0; i < spos->size(); ++i) {
                 if (isdigit((*spos)[i])) num += (*spos)[i];
            }
            if (num.empty()) return false;

            num_name_mapping.insert(std::make_pair(atoi(num.c_str()), *spos));
       }
       if (num_name_mapping.size() != atom_set.size()) return false;

       for (std::map<int, std::string>::const_iterator mpos = num_name_mapping.begin(); mpos != num_name_mapping.end(); ++mpos) {
            atom_list.push_back(mpos->second);
       }
       return true;
}
