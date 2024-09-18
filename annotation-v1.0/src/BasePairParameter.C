/*
FILE:     BasePairParameter.C
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

#include "BasePairParameter.h"
#include "BasePairParameter_global.h"
#include "CompositeIndex.h"
#include "VectorUtil.h"
#include "utillib.h"

void BasePairParameter::initialize()
{
       _base_type_index.clear();
       for (int i = 0; i < NUM_BASE_TYPE; ++i) {
            _base_type_index.insert(std::make_pair(_base_type[i].type, i));
       }
       _base_hbond_type_index.clear();
       for (int i = 0; i < NUM_BASE_HBOND_TYPE; ++i) {
            _base_hbond_type_index.insert(std::make_pair(_base_hbond_type[i].type, i));
       }
}

void BasePairParameter::getParameter(const _BASE& base_i, const _BASE& base_j, const std::map<std::string, std::vector<CONTACT> >& found_contacts,
                                     int& type, int& type_i, int& type_j, int& cis_or_trans, std::vector<CONTACT>& contacts)
{
       contacts.clear();
       type = -1;
       type_i = -1;
       type_j = -1;
       cis_or_trans = -1;

       std::string name_i = base_i.code;
       if (name_i == "P") name_i = "U";
       std::map<std::string, int>::iterator pos = _base_type_index.find(name_i);
       if (pos == _base_type_index.end()) return;
       _DONOR_ACCEPTOR& exist_i = _base_type[pos->second];

       std::string name_j = base_j.code;
       if (name_j == "P") name_j = "U";
       pos = _base_type_index.find(name_j);
       if (pos == _base_type_index.end()) return;
       _DONOR_ACCEPTOR& exist_j = _base_type[pos->second];

       std::vector<std::vector<CONTACT> > contact_list;
       contact_list.clear();
       std::vector<std::pair<std::string, std::string> > atom_pair;
       atom_pair.clear();
       std::vector<int> atom_pair_type;
       atom_pair_type.clear();

       std::vector<std::string> donor, acceptor;
       std::vector<int> donor_type, acceptor_type;

       donor.clear();
       donor_type.clear();
       for (int i = 0; i < exist_i.da->n_donor; ++i) {
            donor.push_back(exist_i.da->donor[i]);
            donor_type.push_back(exist_i.da->donor_type[i]);
       }

       acceptor.clear();
       acceptor_type.clear();
       for (int i = 0; i < exist_j.da->n_acceptor; ++i) {
            acceptor.push_back(exist_j.da->acceptor[i]);
            acceptor_type.push_back(exist_j.da->acceptor_type[i]);
       }
       
       _getAtomPairInfo(base_i, base_j, found_contacts, donor, donor_type, acceptor, acceptor_type, atom_pair, atom_pair_type, contact_list);

       donor.clear();
       donor_type.clear();
       for (int i = 0; i < exist_i.da->n_acceptor; ++i) {
            donor.push_back(exist_i.da->acceptor[i]);
            donor_type.push_back(exist_i.da->acceptor_type[i]);
       }

       acceptor.clear();
       acceptor_type.clear();
       for (int i = 0; i < exist_j.da->n_donor; ++i) {
            acceptor.push_back(exist_j.da->donor[i]);
            acceptor_type.push_back(exist_j.da->donor_type[i]);
       }

       _getAtomPairInfo(base_i, base_j, found_contacts, donor, donor_type, acceptor, acceptor_type, atom_pair, atom_pair_type, contact_list);

       if (contact_list.empty()) return;

       bool found = false;
       for (std::vector<int>::const_iterator pos = atom_pair_type.begin(); pos != atom_pair_type.end(); ++pos) {
            if (*pos == 2) {
                 found = true;
                 break;
            }
       }
       if (!found) return;

       _getBasePairAndType(name_i, name_j, atom_pair, atom_pair_type, contact_list, type, contacts);
       type_i = 0;
       type_j = 0;
       cis_or_trans = 0;
       if (type < 0 || type == HBOND_TYPE_19 || type == HBOND_TYPE_20) return;

       std::map<std::string, int> edge_i, edge_j;
       edge_i.clear();
       edge_j.clear();
       for (int i = 0; i < exist_i.be->num; ++i) {
            edge_i.insert(std::make_pair(exist_i.be->_atoms[i].atom, exist_i.be->_atoms[i].edge));
       }
       for (int i = 0; i < exist_j.be->num; ++i) {
            edge_j.insert(std::make_pair(exist_j.be->_atoms[i].atom, exist_j.be->_atoms[i].edge));
       }

       _getBaseEdgeType(atom_pair, edge_i, edge_j, type_i, type_j);

       cis_or_trans = _getCisTrans(base_i, base_j);
}

void BasePairParameter::_getAtomPairInfo(const _BASE& base_i, const _BASE& base_j, const std::map<std::string, std::vector<CONTACT> >& contacts,
                                         const std::vector<std::string>& donor, const std::vector<int>& donor_type, const std::vector<std::string>& acceptor,
                                         const std::vector<int>& acceptor_type, std::vector<std::pair<std::string, std::string> >& atom_pair,
                                         std::vector<int>& atom_pair_type, std::vector<std::vector<CONTACT> >& contact_list)
{
       RCSB::Atom* atom_i = base_i.res->GetFirstAtom();
       if (!atom_i) return;
       RCSB::Atom* atom_j = base_j.res->GetFirstAtom();
       if (!atom_j) return;

       std::vector<std::string> data;
       for (unsigned int i = 0; i < donor.size(); ++i) {
            for (unsigned int j = 0; j < acceptor.size(); ++j) {
                 data.clear();
                 atom_i->getResidueIndex(data);
                 data.push_back(donor[i]);
                 data.push_back(String::IntToString(base_j.op + 1));
                 data.push_back(String::IntToString(base_j.lx));
                 data.push_back(String::IntToString(base_j.ly));
                 data.push_back(String::IntToString(base_j.lz));
                 atom_j->getResidueIndex(data);
                 data.push_back(acceptor[j]);
                 std::string index = CompositeIndex::getIndex(data);
                 std::map<std::string, std::vector<CONTACT> >::const_iterator mpos = contacts.find(index);
                 if (mpos == contacts.end()) continue;

                 atom_pair.push_back(std::make_pair(donor[i], acceptor[j]));
                 int type = 0;
                 if ((donor_type[i] & 2) & (acceptor_type[j] & 2)) type = 2;
                 else if ((donor_type[i] & 1) & (acceptor_type[j] & 1)) type = 1;
                 atom_pair_type.push_back(type);
                 contact_list.push_back(mpos->second);
            }
       }
}

void BasePairParameter::_getBasePairAndType(const std::string& name_i, const std::string& name_j, const std::vector<std::pair<std::string, std::string> >&
                                            atom_pair, const std::vector<int>& atom_pair_type, const std::vector<std::vector<CONTACT> >& bs_contacts,
                                            int& type, std::vector<CONTACT>& return_contacts)
{
       std::string cs = name_i + name_j;
       bool reverse = false;
       if (name_i > name_j) {
            cs = name_j + name_i;
            reverse = true;
       }

       std::map<std::string, int>::const_iterator pos = _base_hbond_type_index.find(cs);
       if (pos == _base_hbond_type_index.end()) return;
       BASE_HBOND_TYPE& exist = _base_hbond_type[pos->second];

       std::map<std::string, unsigned int> _AtomPairIndex;
       _AtomPairIndex.clear();
       for (unsigned int i = 0; i < atom_pair.size(); ++i) {
             cs = atom_pair[i].first + "_" + atom_pair[i].second;
             if (reverse) cs = atom_pair[i].second + "_" + atom_pair[i].first;
            _AtomPairIndex.insert(std::make_pair(cs, i));
       }

       std::vector<CONTACT> all_contacts;
       std::vector<std::vector<CONTACT> > local_all_contacts, contacts;
       local_all_contacts.clear();
       std::vector<int> local_type_list;
       local_type_list.clear();
       std::vector<double> local_dist_list;
       local_dist_list.clear();
       for (int l = 0; l < exist.btype->num_base_type; ++l) {
            contacts.clear();
            for (int m = 0; m < exist.btype->base_type[l].num; ++m) {
                 cs = exist.btype->base_type[l].atom_pair[m].a_atomtyp;
                 cs += "_";
                 cs += exist.btype->base_type[l].atom_pair[m].b_atomtyp;
                 std::map<std::string, unsigned int>::const_iterator
                     mpos = _AtomPairIndex.find(cs);
                 if (mpos == _AtomPairIndex.end()) break;
                 contacts.push_back(bs_contacts[mpos->second]);
            }
            if (exist.btype->base_type[l].num == (int) contacts.size()) {
                 double dist = _getBasePair(contacts, exist.btype->base_type[l].type, all_contacts);
                 if (!all_contacts.empty()) {
                      local_all_contacts.push_back(all_contacts);
                      local_type_list.push_back(exist.btype->base_type[l].type);
                      local_dist_list.push_back(dist);        
                 }
            }
       }
       if (!local_all_contacts.empty()) {
            int index = _getBestIndex(local_type_list, local_dist_list);
            return_contacts = local_all_contacts[index];
            type = local_type_list[index];
            return;
       }

       contacts.clear();
       for (unsigned int l = 0; l < atom_pair_type.size(); ++l) {
            if (atom_pair_type[l] == 2) contacts.push_back(bs_contacts[l]);
       }
       if (contacts.empty()) return;

       double min_dist = 100000.0;
       int index = -1;
       for (unsigned int i = 0; i < contacts.size(); ++i) {
            double dist = contacts[i][0].dist;
            for (std::vector<CONTACT>::iterator cpos = contacts[i].begin(); cpos != contacts[i].end(); ++cpos) {
                 cpos->type = 0;
                 if (cpos->dist < dist) dist = cpos->dist;
            }
            if (dist < min_dist) {
                 min_dist = dist;
                 index = i;
            }
       }
       if (index < 0) return;
       return_contacts = contacts[index];
       type = 0;
}

int BasePairParameter::_getBestIndex(const std::vector<int>& type_list, const std::vector<double>& dist_list)
{
       int index = 0;
       if (type_list.size() > 1) {
            std::vector<int> wc_type, reverse_wc_type, other_type, no_type;
            wc_type.clear();
            reverse_wc_type.clear();
            other_type.clear();
            no_type.clear();
            for (unsigned int i = 0; i < type_list.size(); ++i) {
                 if (type_list[i] == HBOND_TYPE_19 || type_list[i] == HBOND_TYPE_20) wc_type.push_back(i);
                 else if (type_list[i] == HBOND_TYPE_21 || type_list[i] == HBOND_TYPE_22 || type_list[i] == HBOND_TYPE_23 || type_list[i] == HBOND_TYPE_24)
                     reverse_wc_type.push_back(i);
                 else if (type_list[i] > 0 && type_list[i] < 30) other_type.push_back(i);
                 no_type.push_back(i);
            }
            if (!wc_type.empty())
                 index = _getBestIndexByDistance(wc_type, dist_list);
            else if (!reverse_wc_type.empty())
                 index = _getBestIndexByDistance(reverse_wc_type, dist_list);
            else if (!other_type.empty())
                 index = _getBestIndexByDistance(other_type, dist_list);
            else index = _getBestIndexByDistance(no_type, dist_list);
       }
       return index;
}

int BasePairParameter::_getBestIndexByDistance(const std::vector<int>& index_list, const std::vector<double>& dist_list)
{
       if (index_list.size() == 1) return index_list[0];

       int index = index_list[0];
       double d_min = dist_list[index];
       for (unsigned int i = 1; i < index_list.size(); ++i) {
            if (dist_list[index_list[i]] < d_min) {
                 d_min = dist_list[index_list[i]];
                 index = index_list[i];
            }
       }
       return index;
}

void BasePairParameter::_getBaseEdgeType(const std::vector<std::pair<std::string, std::string> >& atom_pair, const std::map<std::string, int>& edge_i,
                                         const std::map<std::string, int>& edge_j, int& type_i, int& type_j) 
{
       int watson_i = 0;
       int hoogsteen_i = 0;
       int sugar_i = 0;
       int watson_j = 0;
       int hoogsteen_j = 0;
       int sugar_j = 0;
       for (unsigned int i = 0; i < atom_pair.size(); ++i) {
            std::map<std::string, int>::const_iterator
                mpos_i = edge_i.find(atom_pair[i].first);
            if (mpos_i == edge_i.end()) continue;
            std::map<std::string, int>::const_iterator
                mpos_j = edge_j.find(atom_pair[i].second);
            if (mpos_j == edge_j.end()) continue;

            if (mpos_i->second & 1) watson_i++;
            if (mpos_i->second & 2) hoogsteen_i++;
            if (mpos_i->second & 4) sugar_i++;
            if (mpos_j->second & 1) watson_j++;
            if (mpos_j->second & 2) hoogsteen_j++;
            if (mpos_j->second & 4) sugar_j++;
       }
       if ((watson_i > 1 || hoogsteen_i > 1 || sugar_i > 1) &&
           (watson_j > 1 || hoogsteen_j > 1 || sugar_j > 1)) {
            int max1 = (watson_i >= hoogsteen_i) ? watson_i : hoogsteen_i;
            int max2 = (hoogsteen_i >= sugar_i) ? hoogsteen_i : sugar_i;
            int max = (max1 >= max2) ? max1 : max2;
            if (max == watson_i)
                 type_i = WATSON;
            else if (max == hoogsteen_i)
                 type_i = HOOGSTEEN;
            else if (max == sugar_i)
                 type_i = SUGAR;
            max1 = (watson_j >= hoogsteen_j) ? watson_j : hoogsteen_j;
            max2 = (hoogsteen_j >= sugar_j) ? hoogsteen_j : sugar_j;
            max = (max1 >= max2) ? max1 : max2;
            if (max == watson_j)
                 type_j = WATSON;
            else if (max == hoogsteen_j)
                 type_j = HOOGSTEEN;
            else if (max == sugar_j)
                 type_j = SUGAR;
       }
}

int BasePairParameter::_getCisTrans(const _BASE& base_i, const _BASE& base_j)
{
       std::map<std::string, COORD>::const_iterator
           mpos_i_N = base_i.coord_xyz.find("N");
       if (mpos_i_N == base_i.coord_xyz.end()) return 0;
       std::map<std::string, COORD>::const_iterator
           mpos_i_C = base_i.coord_xyz.find("C");
       if (mpos_i_C == base_i.coord_xyz.end()) return 0;
       std::map<std::string, COORD>::const_iterator
           mpos_j_N = base_j.coord_xyz.find("N");
       if (mpos_j_N == base_j.coord_xyz.end()) return 0;
       std::map<std::string, COORD>::const_iterator
           mpos_j_C = base_j.coord_xyz.find("C");
       if (mpos_j_C == base_j.coord_xyz.end()) return 0;

       COORD coord1, coord2, coord, cross_coord1, cross_coord2;
       vector_difference(coord1, mpos_i_C->second, mpos_i_N->second);
       vector_difference(coord2, mpos_j_C->second, mpos_j_N->second);
       vector_difference(coord, mpos_j_N->second, mpos_i_N->second);
       vector_cross_product(cross_coord1, coord1, coord);
       vector_cross_product(cross_coord2, coord2, coord);
       if (vector_dot_product(cross_coord1, cross_coord2) > 0)
            return 1;
       else return 2;
}

double BasePairParameter::_getBasePair(const std::vector<std::vector<CONTACT> >& contacts, const int& type, std::vector<CONTACT>& final_contacts)
{
       final_contacts.clear();

       double min_dist = 100000.0;
       std::vector<std::vector<CONTACT> > pair_list;
       _getPairList(contacts, pair_list, 0);
       if (pair_list.empty()) return min_dist;

       for (std::vector<std::vector<CONTACT> >::iterator pos = pair_list.begin(); pos != pair_list.end(); ++pos) {
            if (pos->size() != contacts.size()) continue;
            for (std::vector<CONTACT>::iterator ppos = pos->begin(); ppos != pos->end(); ++ppos) {
                 ppos->type = type;
                 if (ppos->dist < min_dist) min_dist = ppos->dist;
                 final_contacts.push_back(*ppos);
            }
       }
       return min_dist;
}

void BasePairParameter::_getPairList(const std::vector<std::vector<CONTACT> >& in_list, std::vector<std::vector<CONTACT> >& out_list, const unsigned int& index)
{
       if (index == 0) out_list.clear();
       if (index == in_list.size()) return;

       std::vector<std::vector<CONTACT> > tmp_list;
       std::vector<CONTACT> data;

       tmp_list.clear();

       if (index == 0) {
            for (std::vector<CONTACT>::const_iterator pos = in_list[index].begin(); pos != in_list[index].end(); ++pos) {
                 data.clear();
                 data.push_back(*pos);
                 tmp_list.push_back(data);
            }
       } else {
            for (std::vector<std::vector<CONTACT> >::const_iterator ppos = out_list.begin(); ppos != out_list.end(); ++ppos) {
                 std::set<std::string> existing_alts = _getAltLocSet(*ppos);
                 std::set<std::string> new_alts = _getAltLocSet(in_list[index]);
                 std::set<std::string> common_alts = get_set_common_element(existing_alts, new_alts);
                 for (std::vector<CONTACT>::const_iterator pos = in_list[index].begin(); pos != in_list[index].end(); ++pos) {
                      std::string alt_loc = pos->a_atm->alt_loc();
                      if (!(pos->b_atm->alt_loc().empty())) alt_loc = pos->b_atm->alt_loc();
                      if (!alt_loc.empty() && !common_alts.empty() && common_alts.find(alt_loc) == common_alts.end()) continue;
                      data = *ppos;
                      data.push_back(*pos);
                      tmp_list.push_back(data);
                 }
            }
       }
       out_list = tmp_list;
       _getPairList(in_list, out_list, index + 1);
}

std::set<std::string> BasePairParameter::_getAltLocSet(const std::vector<CONTACT>& in_list)
{
       std::set<std::string> loc_set;
       loc_set.clear();
       for (std::vector<CONTACT>::const_iterator pos = in_list.begin(); pos != in_list.end(); ++pos) {
            if (!(pos->a_atm->alt_loc().empty())) loc_set.insert(pos->a_atm->alt_loc());
            if (!(pos->b_atm->alt_loc().empty())) loc_set.insert(pos->b_atm->alt_loc());
       }
       return loc_set;
}
