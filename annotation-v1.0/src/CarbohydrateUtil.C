/*
FILE:     CarbohydrateUtil.C
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
#include "CarbohydrateUtil.h"
#include "CompositeIndex.h"
#include "GetPairList.h"
#include "GetPairList.C"
#include "utillib.h"

CarbohydrateUtil::CarbohydrateUtil()
{
       _clear();
}

CarbohydrateUtil::~CarbohydrateUtil()
{
       _clear();
}

void CarbohydrateUtil::setLog(LogUtil* logPt)
{
       _logIo = logPt;
}

void CarbohydrateUtil::setCCDic(ConnectDic* ccdic)
{
       _ccDic = ccdic;
}

void CarbohydrateUtil::setPDBID(const std::string& id)
{
       _pdb_id = id;
}

void CarbohydrateUtil::addResidue(const std::vector<RCSB::Residue*>& residues, const std::string& type)
{
       if (residues.empty()) return;

       if (type == "ATOMP") _protein_residue = residues[0];
       else {
            std::string id = residues[0]->pdb_chnid() + "_" + residues[0]->ResName() + "_" + residues[0]->pdb_res_no() + "_" +  residues[0]->ins_code();
            _sugar_residue_index_map.insert(std::make_pair(id, _sugar_residue_list.size()));

/* only index first residue, ignore the rest residues.
            for (unsigned int i = 1; i < residues.size(); ++i) {
                 std::string tid = residues[i]->pdb_chnid() + "_" + residues[i]->ResName() + "_" + residues[i]->pdb_res_no() + "_" +  residues[i]->ins_code();
                 _sugar_residue_index_map.insert(std::make_pair(tid, _sugar_residue_list.size()));
            }
*/

            SugarResidue sugar_residue;
            _clear_sugar_residue(sugar_residue);
            sugar_residue._residues = residues;
            sugar_residue._id = id;
            sugar_residue._name = residues[0]->ResName();
            try {
                 const ConnectFormat& drug = _ccDic->find_drug(residues[0]->ResName());
                 sugar_residue._chemical_name = drug.chemical_name();
                 sugar_residue._condensed_symbol = drug.sugar_condensed_iupac_symbol();
                 sugar_residue._snfg_symbol = drug.sugar_snfg_symbol();
                 sugar_residue._primary_carbonyl_group = drug.sugar_primary_carbonyl_group();
                 sugar_residue._anomeric_carbon = drug.anomeric_carbon();
                 sugar_residue._anomeric_oxygen = drug.anomeric_oxygen();
/*
                 if (sugar_residue._primary_carbonyl_group == "aldose") {
                      sugar_residue._anomeric_carbon = "C1";
                      sugar_residue._anomeric_oxygen = "O1";
                 } else if (sugar_residue._primary_carbonyl_group == "ketose") {
                      sugar_residue._anomeric_carbon = "C2";
                      sugar_residue._anomeric_oxygen = "O2";
                 }

                 if ((sugar_residue._condensed_symbol.empty() || sugar_residue._snfg_symbol.empty() || sugar_residue._primary_carbonyl_group.empty()) &&
                     (_missing_info_ccd_id_set.find(residues[0]->ResName()) == _missing_info_ccd_id_set.end())) {
                      _missing_info_ccd_id_set.insert(residues[0]->ResName());
                      std::string error = "";
                      if (sugar_residue._primary_carbonyl_group.empty()) {
                           error = "Missing 'CARBOHYDRATE PRIMARY CARBONYL GROUP' feature";
                      }
                      int count = 0;
                      if (sugar_residue._condensed_symbol.empty()) {
                           count++;
                           if (error.empty())
                                error = "Missing 'CONDENSED IUPAC CARBOHYDRATE SYMBOL'";
                           else error += ", 'CONDENSED IUPAC CARBOHYDRATE SYMBOL'";
                      }
                      if (sugar_residue._snfg_symbol.empty()) {
                           count++;
                           if (error.empty())
                                error = "Missing 'SNFG CARBOHYDRATE SYMBOL' identifier";
                           else if (count)
                                error += " and 'SNFG CARBOHYDRATE SYMBOL' identifiers";
                      } else if (count) error += " identifier";
                      error += " in chemical component definition for '" + residues[0]->ResName() + "'.\n";
                      _logIo->message(error.c_str());
                 }
*/
            } catch (const std::exception& exc) {
                 if (_missing_info_ccd_id_set.find(residues[0]->ResName()) == _missing_info_ccd_id_set.end()) {
                      _missing_info_ccd_id_set.insert(residues[0]->ResName());
                      std::string error = "Can not find chemical component definition for '" + residues[0]->ResName() + "'.\n";
                      _logIo->message(error.c_str());
                      _error_flag = true;
                 }
            }
            _sugar_residue_list.push_back(sugar_residue);
       }
}

void CarbohydrateUtil::addLinks(const std::list<_LINK>& links)
{
       if (_error_flag) return;

       std::vector<std::pair<unsigned int, unsigned int> > reduced_links;
       reduced_links.clear();
       std::set<std::string> link_set, residue_link_set;
       link_set.clear();
       residue_link_set.clear();

       for (std::list<_LINK>::const_iterator lpos = links.begin(); lpos != links.end(); ++lpos) {
            if ((!lpos->SymOP_1.empty() && lpos->SymOP_1 != "1_555") || (!lpos->SymOP_2.empty() && lpos->SymOP_2 != "1_555")) continue;

            if (link_set.find(lpos->fstAtom->getAtomAllIndex(false) + "-" + lpos->sndAtom->getAtomAllIndex(false)) != link_set.end()) continue;

            link_set.insert(lpos->fstAtom->getAtomAllIndex(false) + "-" + lpos->sndAtom->getAtomAllIndex(false));
            link_set.insert(lpos->sndAtom->getAtomAllIndex(false) + "-" + lpos->fstAtom->getAtomAllIndex(false));

            std::string id1 = lpos->fstAtom->pdb_chnid() + "_" + lpos->fstAtom->pdb_resnam() + "_"
                             + lpos->fstAtom->pdb_resnum() + "_" + lpos->fstAtom->ins_code();
            std::map<std::string, unsigned int>::const_iterator mpos1 = _sugar_residue_index_map.find(id1);

            std::string id2 = lpos->sndAtom->pdb_chnid() + "_" + lpos->sndAtom->pdb_resnam() + "_"
                             + lpos->sndAtom->pdb_resnum() + "_" + lpos->sndAtom->ins_code();
            std::map<std::string, unsigned int>::const_iterator mpos2 = _sugar_residue_index_map.find(id2);

            if (_protein_residue != NULL) {
                 std::string id = _protein_residue->pdb_chnid() + "_" + _protein_residue->ResName() + "_"
                                 + _protein_residue->pdb_res_no() + "_" +  _protein_residue->ins_code();
                 if ((id == id1) && (mpos2 != _sugar_residue_index_map.end())) {
                      _sugar_residue_list[mpos2->second]._attached_protein_residue = _protein_residue;
                      continue;
                 } else if ((id == id2) && (mpos1 != _sugar_residue_index_map.end())) {
                      _sugar_residue_list[mpos1->second]._attached_protein_residue = _protein_residue;
                      continue;
                 }
            }
            if ((mpos1 == _sugar_residue_index_map.end()) || (mpos2 == _sugar_residue_index_map.end())) continue;

            // skip the link(s) if the other alternative atom was already transfered to other residue. 
            // if (mpos1->second == mpos2->second) continue;

            // DO proper O1/O2 anomeric oxygen transfer for non overlap residues
            if ((_sugar_residue_list[mpos1->second]._residues.size() == 1) && (_sugar_residue_list[mpos2->second]._residues.size() == 1) &&
                !_sugar_residue_list[mpos1->second]._anomeric_oxygen.empty() && !_sugar_residue_list[mpos1->second]._anomeric_carbon.empty() &&
                !_sugar_residue_list[mpos2->second]._anomeric_oxygen.empty() && !_sugar_residue_list[mpos2->second]._anomeric_carbon.empty()) {
                 if ((lpos->fstAtom->pdb_atmnam() == _sugar_residue_list[mpos1->second]._anomeric_oxygen) && (lpos->fstAtom->atom_type() == "O") &&
                     (lpos->sndAtom->atom_type() == "C") && (lpos->sndAtom->pdb_atmnam() != _sugar_residue_list[mpos2->second]._anomeric_carbon)) {
                      std::string target_oxygen = lpos->sndAtom->pdb_atmnam();
                      target_oxygen[0] = 'O';
                      if (_transfer_oxygen_atom(_sugar_residue_list[mpos1->second]._residues[0], _sugar_residue_list[mpos1->second]._anomeric_carbon,
                          lpos->fstAtom->pdb_atmnam(), _sugar_residue_list[mpos2->second]._residues[0], lpos->sndAtom->pdb_atmnam(), target_oxygen)) {
                           _link_list.push_back(std::make_pair(std::make_pair(mpos1->second, _sugar_residue_list[mpos1->second]._anomeric_carbon),
                                                               std::make_pair(mpos2->second, target_oxygen)));
                           if (residue_link_set.find(String::IntToString(mpos1->second) + "-" + String::IntToString(mpos2->second)) == residue_link_set.end()) {
                                reduced_links.push_back(std::make_pair(mpos1->second, mpos2->second));
                                residue_link_set.insert(String::IntToString(mpos1->second) + "-" + String::IntToString(mpos2->second));
                                residue_link_set.insert(String::IntToString(mpos2->second) + "-" + String::IntToString(mpos1->second));
                           }
                           continue;
                      }
                 } else if ((lpos->sndAtom->pdb_atmnam() == _sugar_residue_list[mpos2->second]._anomeric_oxygen) && (lpos->sndAtom->atom_type() == "O") &&
                            (lpos->fstAtom->atom_type() == "C") && (lpos->fstAtom->pdb_atmnam() != _sugar_residue_list[mpos1->second]._anomeric_carbon)) {
                      std::string target_oxygen = lpos->fstAtom->pdb_atmnam();
                      target_oxygen[0] = 'O';
                      if (_transfer_oxygen_atom(_sugar_residue_list[mpos2->second]._residues[0], _sugar_residue_list[mpos2->second]._anomeric_carbon,
                          lpos->sndAtom->pdb_atmnam(), _sugar_residue_list[mpos1->second]._residues[0], lpos->fstAtom->pdb_atmnam(), target_oxygen)) {
                           _link_list.push_back(std::make_pair(std::make_pair(mpos2->second, _sugar_residue_list[mpos2->second]._anomeric_carbon),
                                                               std::make_pair(mpos1->second, target_oxygen)));
                           if (residue_link_set.find(String::IntToString(mpos1->second) + "-" + String::IntToString(mpos2->second)) == residue_link_set.end()) {
                                reduced_links.push_back(std::make_pair(mpos1->second, mpos2->second));
                                residue_link_set.insert(String::IntToString(mpos1->second) + "-" + String::IntToString(mpos2->second));
                                residue_link_set.insert(String::IntToString(mpos2->second) + "-" + String::IntToString(mpos1->second));
                           }
                           continue;
                      }
                 }
            }

            _link_list.push_back(std::make_pair(std::make_pair(mpos1->second, lpos->fstAtom->pdb_atmnam()),
                                                std::make_pair(mpos2->second, lpos->sndAtom->pdb_atmnam())));
            if (residue_link_set.find(String::IntToString(mpos1->second) + "-" + String::IntToString(mpos2->second)) == residue_link_set.end()) {
                 reduced_links.push_back(std::make_pair(mpos1->second, mpos2->second));
                 residue_link_set.insert(String::IntToString(mpos1->second) + "-" + String::IntToString(mpos2->second));
                 residue_link_set.insert(String::IntToString(mpos2->second) + "-" + String::IntToString(mpos1->second));
            }
       }

       std::vector<std::set<unsigned int> > groups;
       clustering_with_insertion(groups, reduced_links);
       if ((groups.size() == 1) && (groups[0].size() == _sugar_residue_list.size())) return;

       for (unsigned int i = 0; i < _sugar_residue_list.size() - 1; ++i) {
            for (unsigned int j = i + 1; j < _sugar_residue_list.size(); ++j) {
                 if (residue_link_set.find(String::IntToString(i) + "-" + String::IntToString(j)) != residue_link_set.end()) continue;
                 if (_find_additional_links(i, j, "O", "C")) {
                      reduced_links.push_back(std::make_pair(i, j));
                      continue;
                 }
                 if (_find_additional_links(i, j, "C", "O")) reduced_links.push_back(std::make_pair(i, j));
            }
       }

       clustering_with_insertion(groups, reduced_links);
       if ((groups.size() != 1) || (groups[0].size() != _sugar_residue_list.size())) {
            _error_flag = true;
            std::string error = "The sugar polymer [";
            bool first = true;
            for (std::map<std::string, unsigned int>::const_iterator mpos = _sugar_residue_index_map.begin(); mpos != _sugar_residue_index_map.end(); ++mpos) {
                 if (!first) error += " ";
                 error += mpos->first;
                 first = false;
            }
            error +="] is not connected.\n";
            _logIo->message(error.c_str());
       }
}

void CarbohydrateUtil::buildOligoSaccharide()
{
       if (_error_flag) return;

       _update_linked_residue_list();

       std::vector<unsigned int> reducing_root, anomeric_root;
       reducing_root.clear();
       anomeric_root.clear();

       for (unsigned int i = 0; i < _sugar_residue_list.size(); ++i) {
            if ((_sugar_residue_list[i]._attached_protein_residue != NULL) || (_sugar_residue_list[i]._reducing_linkage_list.empty() &&
                _sugar_residue_list[i]._anomeric_anomeric_linkage_list.empty() && !_sugar_residue_list[i]._non_reducing_linkage_list.empty())) {
                 reducing_root.push_back(i);
                 break;
            }
       }
       if (reducing_root.empty()) {
            for (unsigned int i = 0; i < _sugar_residue_list.size(); ++i) {
                 if (!_sugar_residue_list[i]._reducing_linkage_list.empty() || _sugar_residue_list[i]._anomeric_oxygen.empty()) continue;
                 if ((_sugar_residue_list[i]._anomeric_anomeric_linkage_list.size() + _sugar_residue_list[i]._non_reducing_linkage_list.size()) != 1) continue;
                 bool found_linked_anomeric_oxygen = false;
                 for (unsigned int j = 0; j < _sugar_residue_list[i]._anomeric_anomeric_linkage_list.size(); ++j) {
                      if (_sugar_residue_list[i]._anomeric_oxygen == _sugar_residue_list[i]._anomeric_anomeric_linkage_list[j]._current_atom) {
                           found_linked_anomeric_oxygen = true;
                           break;
                      }
                 }
                 for (unsigned int j = 0; j < _sugar_residue_list[i]._non_reducing_linkage_list.size(); ++j) {
                      if (_sugar_residue_list[i]._anomeric_oxygen == _sugar_residue_list[i]._non_reducing_linkage_list[j]._current_atom) {
                           found_linked_anomeric_oxygen = true;
                           break;
                      }
                 }
                 if (found_linked_anomeric_oxygen) continue;

                 RCSB::Atom* atom = _sugar_residue_list[i]._residues[0]->find_atom(_sugar_residue_list[i]._anomeric_oxygen);
                 if (atom) {
                      reducing_root.push_back(i);
                      break;
                 }
            }
       }
       if (reducing_root.empty()) {
            for (unsigned int i = 0; i < _sugar_residue_list.size(); ++i) {
                 if (!_sugar_residue_list[i]._anomeric_anomeric_linkage_list.empty() || !_sugar_residue_list[i]._reducing_linkage_list.empty()) {
                      if (anomeric_root.empty()) anomeric_root.push_back(i);
                      else if (anomeric_root[0] > i) anomeric_root[0] = i;
                 }
            }
       }

       unsigned int root = 0;
       if (!reducing_root.empty()) root = reducing_root[0];
       else {
            std::vector<std::vector<unsigned int> > ringLists;
            ringLists.clear();
            if (_sugar_residue_list.size() > 3) find_rings(_sugar_residue_list.size(), _linked_residue_index_mapping, ringLists);

            if (ringLists.size() == 1) {
                 _is_cyclic = true;
                 if (ringLists[0].size() == _sugar_residue_list.size()) {
                      if (_update_cyclic_saccharide_info(root)) return;
                 } else {
                      for (std::vector<unsigned int>::const_iterator vpos = ringLists[0].begin(); vpos != ringLists[0].end(); ++vpos) {
                           unsigned size = _sugar_residue_list[*vpos]._reducing_linkage_list.size()
                                         + _sugar_residue_list[*vpos]._non_reducing_linkage_list.size()
                                         + _sugar_residue_list[*vpos]._anomeric_anomeric_linkage_list.size();
                           if (size > 2) {
                                root = *vpos;
                                break;
                           }
                      }
                 }
            } else if (!anomeric_root.empty()) {
                 root = anomeric_root[0];
            } else {
                 _error_flag = true;
                 std::string error = "Can not find the reducing end for sugar polymer [";
                 bool first = true;
                 for (std::map<std::string, unsigned int>::const_iterator mpos = _sugar_residue_index_map.begin();
                      mpos != _sugar_residue_index_map.end(); ++mpos) {
                      if (!first) error += " ";
                      error += mpos->first;
                      first = false;
                 }
                 error +="].\n";
                 _logIo->message(error.c_str());
                 return;
            }
       }

       _sugar_residue_list[root]._is_root = true;
       _traverse_graph(root, _descriptor_in_full_form, _descriptor_in_condensed_form, _max_branch_length);

       if (_descriptor_id_set.size() != _sugar_residue_index_map.size()) {
            std::string missing_residues = "";
            for (std::map<std::string, unsigned int>::const_iterator mpos = _sugar_residue_index_map.begin(); mpos != _sugar_residue_index_map.end(); ++mpos) {
                 if (_descriptor_id_set.find(mpos->first) != _descriptor_id_set.end()) continue;
                 if (!missing_residues.empty()) missing_residues += ", ";
                 missing_residues += mpos->first;
            }
            if (!missing_residues.empty()) {
                 std::string error = "Missing the following residue(s) ( " + missing_residues + " ) in the descriptor.\n";
                 _logIo->message(error.c_str());
                 _error_flag = true;
            }
       }

       if (!_error_flag && !_descriptor_in_condensed_form.empty()) {
            _reindex_residues(_descriptor_in_condensed_form);
            _generate_new_index_map();
            _transfer_oxygen_atom_list();
       }

       if (_is_cyclic) {
            _cyclic_descriptor.clear();
            std::string descriptor = getEntityDescriptor();
            if (!descriptor.empty()) _cyclic_descriptor = "Cyclic " + descriptor;
            _descriptor_in_full_form.clear();
            _descriptor_in_condensed_form.clear();
       }
}

bool CarbohydrateUtil::getStatus()
{
       return !_error_flag;
}

std::string CarbohydrateUtil::getEntityDescriptor()
{
       if (_error_flag) return "";
       if (!_cyclic_descriptor.empty()) return _cyclic_descriptor;
       if (_descriptor_in_full_form.empty()) return "";

       std::string descriptor = _descriptor_in_full_form;
       for (std::vector<SugarResidue>::const_iterator vpos = _sugar_residue_list.begin(); vpos != _sugar_residue_list.end(); ++vpos) {
            if (vpos->_chemical_name.empty()) {
                 _error_flag = true;
                 break;
            }
            replace_string(descriptor, "{" + vpos->_id + "}", vpos->_chemical_name);
       }
       if (_error_flag) descriptor.clear();

       return descriptor;
}

std::string CarbohydrateUtil::getCondensedDescriptor()
{
       if (_error_flag || _descriptor_in_condensed_form.empty()) return "";

       std::string descriptor = _descriptor_in_condensed_form;
       for (std::vector<SugarResidue>::const_iterator vpos = _sugar_residue_list.begin(); vpos != _sugar_residue_list.end(); ++vpos) {
            std::string symbol = vpos->_condensed_symbol;
            if (symbol.empty()) symbol = vpos->_snfg_symbol;
            if (symbol.empty()) {
                 _error_flag = true;
                 break;
            }
            replace_string(descriptor, "{" + vpos->_id + "}", symbol);
       }
       if (_error_flag) descriptor.clear();

       return descriptor;
}

const std::map<int, std::string>& CarbohydrateUtil::getReorderIndexMapping() const
{
       return _new_index_map;
}

const std::list<_LINK>& CarbohydrateUtil::getAdditionalLinks() const
{
       return _additional_links;
}

const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& CarbohydrateUtil::getBadLinks() const
{
       return _bad_links;
}

std::string CarbohydrateUtil::getEntityKey()
{
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

       std::vector<std::string> links;
       links.clear();
       for (int i = 0; i < 6; ++i) links.push_back("");

       std::multimap<int, std::vector<std::string> > tmp_link_mapping;

       for (std::vector<std::pair<std::pair<unsigned int, std::string>, std::pair<unsigned int, std::string> > >::const_iterator
            vpos = _link_list.begin(); vpos != _link_list.end(); ++vpos) {
            int small_index = _sugar_residue_list[vpos->first.first]._new_index;
            int large_index = _sugar_residue_list[vpos->second.first]._new_index;
            if (small_index > large_index) {
                 small_index = _sugar_residue_list[vpos->second.first]._new_index;
                 large_index = _sugar_residue_list[vpos->first.first]._new_index;
                 links[0] = _sugar_residue_list[vpos->second.first]._name;
                 links[1] = String::IntToString(small_index);
                 links[2] = vpos->second.second;
                 links[3] = _sugar_residue_list[vpos->first.first]._name;
                 links[4] = String::IntToString(large_index);
                 links[5] = vpos->first.second;
            } else {
                 links[0] = _sugar_residue_list[vpos->first.first]._name;
                 links[1] = String::IntToString(small_index);
                 links[2] = vpos->first.second;
                 links[3] = _sugar_residue_list[vpos->second.first]._name;
                 links[4] = String::IntToString(large_index);
                 links[5] = vpos->second.second;
            }

            std::multimap<int, std::multimap<int, std::vector<std::string> > >::iterator rlmpos = residue_link_mapping.find(small_index);
            if (rlmpos != residue_link_mapping.end()) rlmpos->second.insert(std::make_pair(large_index, links));
            else {
                 tmp_link_mapping.clear();
                 tmp_link_mapping.insert(std::make_pair(large_index, links));
                 residue_link_mapping.insert(std::make_pair(small_index, tmp_link_mapping));
            }
       }

       // Create entity key with residue name, residue number and linked atom names
       std::string entity_key = "";
       for (std::multimap<int, std::multimap<int, std::vector<std::string> > >::const_iterator
            mpos = residue_link_mapping.begin(); mpos != residue_link_mapping.end(); ++mpos) {
            for (std::multimap<int, std::vector<std::string> >::const_iterator mpos1 = mpos->second.begin(); mpos1 != mpos->second.end(); ++mpos1) {
                 if (!entity_key.empty()) entity_key += "_";
                 entity_key += CompositeIndex::getIndex(mpos1->second);
            }
       }
       return entity_key;
}

void CarbohydrateUtil::_clear()
{
       _error_flag = false;
       _is_cyclic = false;
       _pdb_id.clear();
       _logIo = NULL;
       _ccDic = NULL;
       _protein_residue = NULL;
       _sugar_residue_list.clear();
       _sugar_residue_index_map.clear();
       _missing_info_ccd_id_set.clear();
       _link_list.clear();
       _additional_links.clear();
       _bad_links.clear();
       _cyclic_descriptor.clear();
       _descriptor_in_full_form.clear();
       _descriptor_in_condensed_form.clear();
       _max_branch_length = 0;
       _descriptor_id_set.clear();
       _global_index = 1;
       _new_index_map.clear();
       _linked_residue_index_mapping.clear();
       _oxygen_transfer_list.clear();
}

void CarbohydrateUtil::_clear_sugar_residue(SugarResidue& sugar_residue)
{
       sugar_residue._residues.clear();
       sugar_residue._attached_protein_residue = NULL;
       sugar_residue._id.clear();
       sugar_residue._name.clear();
       sugar_residue._chemical_name.clear();
       sugar_residue._condensed_symbol.clear();
       sugar_residue._snfg_symbol.clear();
       sugar_residue._primary_carbonyl_group.clear();
       sugar_residue._anomeric_carbon.clear();
       sugar_residue._anomeric_oxygen.clear();
       sugar_residue._is_visited = false;
       sugar_residue._is_root = false;
       sugar_residue._new_index = 0;
       sugar_residue._reducing_linkage_list.clear();
       sugar_residue._non_reducing_linkage_list.clear();
       sugar_residue._anomeric_anomeric_linkage_list.clear();
}

void CarbohydrateUtil::_clear_sugar_linkage(SugarLinkage& sugar_linkage)
{
       sugar_linkage._current_atom.clear();
       sugar_linkage._current_number = 0;
       sugar_linkage._linked_residue_index = 0;
       sugar_linkage._linked_atom.clear();
       sugar_linkage._linked_number = 0;
       sugar_linkage._linkage_type.clear();
       sugar_linkage._inverse_linkage_type.clear();
}

bool CarbohydrateUtil::_transfer_oxygen_atom(RCSB::Residue* src_residue, const std::string& src_carbon, const std::string& src_oxygen,
                                             RCSB::Residue* tgt_residue, const std::string& tgt_carbon, const std::string& tgt_oxygen)
{
       // verify if tgt_oxygen exists
       RCSB::Atom* atom = tgt_residue->find_atom(tgt_oxygen);
       if (atom) return false;

       // verify if tgt_carbon exists
       atom = tgt_residue->find_atom(tgt_carbon);
       if (!atom) return false;

       try {
            // verify if tgt_oxygen & tgt_carbon are linked together in CCD denifition.
            bool found_bond = false;

            const ConnectFormat& drug = _ccDic->find_drug(tgt_residue->ResName());

            const std::map<std::string, std::string>& name_mapping = drug.getAtomIdToStandardAtomIdMapping();
            const std::vector<std::vector<std::string> >& bondlist = drug.bonds();
            for (std::vector<std::vector<std::string> >::const_iterator bpos = bondlist.begin(); bpos != bondlist.end(); ++bpos) {
                 std::string atom1 = (*bpos)[0];
                 std::map<std::string, std::string>::const_iterator nmpos = name_mapping.find((*bpos)[0]);
                 if (nmpos != name_mapping.end()) atom1 = nmpos->second;

                 std::string atom2 = (*bpos)[1];
                 nmpos = name_mapping.find((*bpos)[1]);
                 if (nmpos != name_mapping.end()) atom2 = nmpos->second;

                 if (((atom1 == tgt_oxygen) && (atom2 == tgt_carbon)) || ((atom1 == tgt_carbon) && (atom2 == tgt_oxygen))) {
                      found_bond = true;
                      break;
                 }
            }
            if (!found_bond) return false;
       } catch (const std::exception& exc) { return false; }

       std::vector<RCSB::Atom*> src_O_atom_list, attached_C_atom_list;
       src_residue->find_atom(src_oxygen, src_O_atom_list);
       src_residue->find_atom(src_carbon, attached_C_atom_list);
       if (!src_O_atom_list.empty() && !attached_C_atom_list.empty()) {
            OxygenTransfer ot;
            ot.src_residue = src_residue;
            ot.src_carbon = src_carbon;
            ot.src_oxygen = src_oxygen;
            ot.tgt_residue = tgt_residue;
            ot.tgt_carbon = tgt_carbon;
            ot.tgt_oxygen = tgt_oxygen;
            _oxygen_transfer_list.push_back(ot);

            return true;
       }

       return false;
}

void CarbohydrateUtil::_transfer_oxygen_atom_list()
{
       if (_oxygen_transfer_list.empty()) return;

       std::vector<RCSB::Atom*> src_O_atom_list, attached_C_atom_list;
       std::vector<std::vector<RCSB::Atom*> > atom_lists, pair_lists;

       _LINK link;
       link.mol_index = 0;
       link.SymOP_1.clear();
       link.SymOP_2.clear();
       link.details.clear();
       link.type = "covale";
       link.bondtype = "sing";
       link.leaving_flag.clear();
       link.pdbx_role.clear();

       for (std::list<OxygenTransfer>::const_iterator lpos = _oxygen_transfer_list.begin(); lpos != _oxygen_transfer_list.end(); ++lpos) {
            RCSB::Atom* atom = lpos->tgt_residue->find_atom(lpos->tgt_carbon);
            if (!atom) continue;

            lpos->src_residue->find_atom(lpos->src_carbon, attached_C_atom_list);
            if (attached_C_atom_list.empty()) continue;

            lpos->src_residue->remove_a_atom(lpos->src_oxygen, src_O_atom_list);
            if (src_O_atom_list.empty()) continue;

            for (unsigned int i = 0; i < src_O_atom_list.size(); ++i) {
                 src_O_atom_list[i]->set_chnid(atom->chnid());
                 src_O_atom_list[i]->set_pdb_chnid(atom->pdb_chnid());
                 src_O_atom_list[i]->set_restype(atom->restype());
                 src_O_atom_list[i]->set_pdb_resnam(atom->pdb_resnam());
                 src_O_atom_list[i]->set_resnum(atom->resnum());
                 src_O_atom_list[i]->set_pdb_resnum(atom->pdb_resnum());
                 src_O_atom_list[i]->set_ins_code(atom->ins_code());
                 src_O_atom_list[i]->set_atmtype(lpos->tgt_oxygen);
                 src_O_atom_list[i]->set_pdb_atmnam(lpos->tgt_oxygen);
                 lpos->tgt_residue->insert_a_atom(src_O_atom_list[i], -1);
            }
            lpos->tgt_residue->reorder_atoms();

            atom_lists.clear();
            atom_lists.push_back(attached_C_atom_list);
            atom_lists.push_back(src_O_atom_list);
            GetPairList(atom_lists, pair_lists);

            if (!pair_lists.empty()) {
                 for (std::vector<std::vector<RCSB::Atom*> >::iterator pos = pair_lists.begin(); pos != pair_lists.end(); ++pos) {      
                      double dist = cal_distance((*pos)[0], (*pos)[1]);
                      link.fstAtom = (*pos)[0];
                      link.sndAtom = (*pos)[1];
                      link.dist = FloatToString(dist, 0, 3);
                      _additional_links.push_back(link);
                 }
            }
       }
}

bool CarbohydrateUtil::_find_additional_links(unsigned int& i, unsigned int& j, const std::string& type_i, const std::string& type_j)
{
       _LINK link;
       link.mol_index = 0;
       link.SymOP_1.clear();
       link.SymOP_2.clear();
       link.details.clear();
       link.type = "covale";
       link.bondtype = "sing";
       link.leaving_flag.clear();
       link.pdbx_role.clear();
       
       std::vector<RCSB::Atom*> atom_i, atom_j;
       _sugar_residue_list[i]._residues[0]->find_atom_by_type(type_i, atom_i);
       _sugar_residue_list[j]._residues[0]->find_atom_by_type(type_j, atom_j);
       for (std::vector<RCSB::Atom*>::const_iterator pos_i = atom_i.begin(); pos_i != atom_i.end(); ++pos_i) {
            for (std::vector<RCSB::Atom*>::const_iterator pos_j = atom_j.begin(); pos_j != atom_j.end(); ++pos_j) {
                 if (!(*pos_i)->alt_loc().empty() && !(*pos_j)->alt_loc().empty() && (*pos_i)->alt_loc() != (*pos_j)->alt_loc()) continue;

                 double dist = cal_distance((*pos_i), (*pos_j));
                 if (BondUtil::is_a_link((*pos_i)->atom_type(), (*pos_j)->atom_type(), dist)) {
                      link.fstAtom = (*pos_i);
                      link.sndAtom = (*pos_j);
                      link.dist = FloatToString(dist, 0, 3);
                      _additional_links.push_back(link);

                      _link_list.push_back(std::make_pair(std::make_pair(i, (*pos_i)->pdb_atmnam()), std::make_pair(j, (*pos_j)->pdb_atmnam())));
                      return true;
                 }
            }
       }
       return false;
}

void CarbohydrateUtil::_update_linked_residue_list()
{
       std::vector<std::pair<std::pair<unsigned int, SugarLinkage>, std::pair<unsigned int, SugarLinkage> > > non_standard_linkage_list;
       non_standard_linkage_list.clear();

       std::set<std::string> added_linkage_set;
       added_linkage_set.clear();

       SugarLinkage first_linkage, second_linkage ;
       for (std::vector<std::pair<std::pair<unsigned int, std::string>, std::pair<unsigned int, std::string> > >::const_iterator
            vpos = _link_list.begin(); vpos != _link_list.end(); ++vpos) {
            bool first_anomeric_flag = false;
            bool second_anomeric_flag = false;
            if (_sugar_residue_list[vpos->first.first]._anomeric_carbon == vpos->first.second) {
                 first_anomeric_flag = true;
/*
                 std::string anomeric_oxygen = _sugar_residue_list[vpos->second.first]._anomeric_carbon;
                 if (!anomeric_oxygen.empty()) anomeric_oxygen[0] = 'O';
                 if (anomeric_oxygen == vpos->second.second) second_anomeric_flag = true;
*/
                 if ((_sugar_residue_list[vpos->second.first]._anomeric_oxygen == vpos->second.second) ||
                     (_sugar_residue_list[vpos->first.first]._anomeric_oxygen == vpos->second.second)) second_anomeric_flag = true;
            } else if (_sugar_residue_list[vpos->second.first]._anomeric_carbon == vpos->second.second) {
                 second_anomeric_flag = true;
/*
                 std::string anomeric_oxygen = _sugar_residue_list[vpos->first.first]._anomeric_carbon;
                 if (!anomeric_oxygen.empty()) anomeric_oxygen[0] = 'O';
                 if (anomeric_oxygen == vpos->first.second) first_anomeric_flag = true;
*/
                 if ((_sugar_residue_list[vpos->first.first]._anomeric_oxygen == vpos->first.second) ||
                     (_sugar_residue_list[vpos->second.first]._anomeric_oxygen == vpos->first.second)) first_anomeric_flag = true;
            }

            if (!_sugar_residue_list[vpos->first.first]._anomeric_carbon.empty() && !_sugar_residue_list[vpos->first.first]._anomeric_oxygen.empty() &&
                !_sugar_residue_list[vpos->second.first]._anomeric_carbon.empty() && !_sugar_residue_list[vpos->second.first]._anomeric_oxygen.empty()) {
                 if (((vpos->first.second == _sugar_residue_list[vpos->first.first]._anomeric_oxygen) &&
                      (vpos->second.second != _sugar_residue_list[vpos->second.first]._anomeric_carbon)) ||
                     ((vpos->second.second == _sugar_residue_list[vpos->second.first]._anomeric_oxygen) &&
                      (vpos->first.second != _sugar_residue_list[vpos->first.first]._anomeric_carbon))) {
                 }
            }

            unsigned int first_number = _get_number_part_from_atom_name(_sugar_residue_list[vpos->first.first]._name, vpos->first.second);
            unsigned int second_number = _get_number_part_from_atom_name(_sugar_residue_list[vpos->second.first]._name, vpos->second.second);

            _clear_sugar_linkage(first_linkage);
            first_linkage._current_atom = vpos->first.second;
            first_linkage._current_number = first_number;
            first_linkage._linked_residue_index = vpos->second.first;
            first_linkage._linked_atom = vpos->second.second;
            first_linkage._linked_number = second_number;
            if ((first_number > 0) && (second_number > 0)) {
                 first_linkage._linkage_type = String::IntToString(first_number) + "-" + String::IntToString(second_number);
                 first_linkage._inverse_linkage_type = String::IntToString(second_number) + "-" + String::IntToString(first_number);
            }

            _clear_sugar_linkage(second_linkage);
            second_linkage._current_atom = vpos->second.second;
            second_linkage._current_number = second_number;
            second_linkage._linked_residue_index = vpos->first.first;
            second_linkage._linked_atom = vpos->first.second;
            second_linkage._linked_number = first_number;
            if ((first_number > 0) && (second_number > 0)) {
                 second_linkage._linkage_type = String::IntToString(second_number) + "-" + String::IntToString(first_number);
                 second_linkage._inverse_linkage_type = String::IntToString(first_number) + "-" + String::IntToString(second_number);
            }

            bool added_flag = false;
            if (first_anomeric_flag && second_anomeric_flag) {
                 _sugar_residue_list[vpos->first.first]._anomeric_anomeric_linkage_list.push_back(first_linkage);
                 _sugar_residue_list[vpos->second.first]._anomeric_anomeric_linkage_list.push_back(second_linkage);
                 added_flag = true;
            } else if (first_anomeric_flag && !second_anomeric_flag) {
                 _sugar_residue_list[vpos->first.first]._reducing_linkage_list.push_back(first_linkage);
                 _sugar_residue_list[vpos->second.first]._non_reducing_linkage_list.push_back(second_linkage);
                 added_flag = true;
            } else if (!first_anomeric_flag && second_anomeric_flag) {
                 _sugar_residue_list[vpos->first.first]._non_reducing_linkage_list.push_back(first_linkage);
                 _sugar_residue_list[vpos->second.first]._reducing_linkage_list.push_back(second_linkage);
                 added_flag = true;
            } else {
                 non_standard_linkage_list.push_back(std::make_pair(std::make_pair(vpos->first.first, first_linkage),
                                                     std::make_pair(vpos->second.first, second_linkage)));
            }

            if (!added_flag || (added_linkage_set.find(String::IntToString(vpos->first.first) + "-" + String::IntToString(vpos->second.first))
                != added_linkage_set.end())) continue;

            _insert_linked_residue_index_mapping(vpos->first.first, vpos->second.first);
            _insert_linked_residue_index_mapping(vpos->second.first, vpos->first.first);
            added_linkage_set.insert(String::IntToString(vpos->first.first) + "-" + String::IntToString(vpos->second.first));
            added_linkage_set.insert(String::IntToString(vpos->second.first) + "-" + String::IntToString(vpos->first.first));
       }

       for (std::vector<std::pair<std::pair<unsigned int, SugarLinkage>, std::pair<unsigned int, SugarLinkage> > >::const_iterator
            vpos = non_standard_linkage_list.begin(); vpos != non_standard_linkage_list.end(); ++vpos) {
            if (added_linkage_set.find(String::IntToString(vpos->first.first) + "-" + String::IntToString(vpos->second.first)) == added_linkage_set.end()) {
                 _insert_linked_residue_index_mapping(vpos->first.first, vpos->second.first);
                 _insert_linked_residue_index_mapping(vpos->second.first, vpos->first.first);
                 _sugar_residue_list[vpos->first.first]._anomeric_anomeric_linkage_list.push_back(vpos->first.second);
                 _sugar_residue_list[vpos->second.first]._anomeric_anomeric_linkage_list.push_back(vpos->second.second);
            } else {
                 std::string header = "Incorrect linkage";
                 RCSB::Residue* res1 = _sugar_residue_list[vpos->first.first]._residues[0];
                 RCSB::Residue* res2 = _sugar_residue_list[vpos->second.first]._residues[0];
                 RCSB::Atom* atom1 = res1->find_atom(vpos->first.second._current_atom);
                 RCSB::Atom* atom2 = res2->find_atom(vpos->second.second._current_atom);
                 if (atom1 && atom2) {
                      _bad_links.push_back(std::make_pair(atom1, atom2));
                      header += " (removed)";
                 }

                 std::string error = header + " : " + res1->pdb_chnid() + " " + res1->ResName() + " " + res1->pdb_res_no() +  res1->ins_code() + " "
                                   + vpos->first.second._current_atom + " -- " + res2->pdb_chnid() + " " + res2->ResName() + " " + res2->pdb_res_no() 
                                   + res2->ins_code() + " " + vpos->second.second._current_atom + ".\n";
                 _logIo->message(error.c_str());
            }
       }
}

int CarbohydrateUtil::_get_number_part_from_atom_name(const std::string& res_name, const std::string& atom_name)
{
       if (atom_name.empty()) return 0;

       if (isdigit(atom_name[atom_name.size() - 1]) != 0) {
            for (unsigned int i = atom_name.size() - 1; i >= 0; --i) {
                 if (isdigit(atom_name[i]) == 0) {
                      return atoi(atom_name.substr(i+1).c_str());
                 }
                 if (i == 0) break;
            }
       }

       try {
            const ConnectFormat& drug = _ccDic->find_drug(res_name);
            int number = drug.get_sugar_atom_locant_number(atom_name);
            if (number > 0) return number;
       } catch (const std::exception& exc) {}

       return 0;
}

void CarbohydrateUtil::_insert_linked_residue_index_mapping(const unsigned int& first_index, const unsigned int& second_index)
{
       std::map<unsigned int, std::vector<unsigned int> >::iterator mpos = _linked_residue_index_mapping.find(first_index);
       if (mpos != _linked_residue_index_mapping.end()) mpos->second.push_back(second_index);
       else {
            std::vector<unsigned int> t_vec;
            t_vec.clear();
            t_vec.push_back(second_index);
            _linked_residue_index_mapping.insert(std::make_pair(first_index, t_vec));
       }
}

bool CarbohydrateUtil::_update_cyclic_saccharide_info(unsigned int& root)
{
       root = 0;
       int min_number = atoi(_sugar_residue_list[root]._residues[0]->pdb_res_no().c_str());
       for (unsigned int i = 1; i < _sugar_residue_list.size(); ++i) {
            int number = atoi(_sugar_residue_list[i]._residues[0]->pdb_res_no().c_str());
            if (number < min_number) {
                 min_number = number;
                 root = i;
            }
       }

       std::vector<unsigned int> cycleList;
       cycleList.clear();
       cycleList.push_back(root);
       while (true) {
            unsigned last = cycleList[cycleList.size() - 1];
            if (!_sugar_residue_list[last]._non_reducing_linkage_list.empty()) {
                 unsigned int next = _sugar_residue_list[last]._non_reducing_linkage_list[0]._linked_residue_index;
                 if (next == root) break;
                 cycleList.push_back(next);
            } else if (!_sugar_residue_list[last]._anomeric_anomeric_linkage_list.empty()) {
                 unsigned int next = _sugar_residue_list[last]._anomeric_anomeric_linkage_list[0]._linked_residue_index;
                 if (next == root) break;
                 cycleList.push_back(next);
            } else break;
       }

       if (cycleList.size() != _sugar_residue_list.size()) return false;

       std::set<std::string> unique_name_set, unique_linkage_type_set;
       unique_name_set.clear();
       unique_linkage_type_set.clear();

       for (unsigned int i = 0; i < cycleList.size(); ++i) {
            unique_name_set.insert(_sugar_residue_list[cycleList[i]]._name);
            if (!_sugar_residue_list[cycleList[i]]._non_reducing_linkage_list.empty()) {
                 unique_linkage_type_set.insert(_sugar_residue_list[cycleList[i]]._non_reducing_linkage_list[0]._inverse_linkage_type);
            } else if (!_sugar_residue_list[cycleList[i]]._anomeric_anomeric_linkage_list.empty()) {
                 unique_linkage_type_set.insert(_sugar_residue_list[cycleList[i]]._anomeric_anomeric_linkage_list[0]._inverse_linkage_type);
            }
       }

       if ((unique_name_set.size() == 1) && (unique_linkage_type_set.size() == 1)) {
            std::string prefix = other_numerical_multiplier(cycleList.size());
            if (!prefix.empty()) {
                 for (unsigned int i = 0; i < cycleList.size(); ++i) {
                      _sugar_residue_list[cycleList[i]]._new_index = _global_index;
                      _global_index += _sugar_residue_list[cycleList[i]]._residues.size();
                 }
                 _generate_new_index_map();
                 _transfer_oxygen_atom_list();
                 _cyclic_descriptor = "Cyclo" + prefix + "(" + (*unique_linkage_type_set.begin()) + ")-(" + _sugar_residue_list[0]._chemical_name + ")";
                 return true;
            }
       }

       return false;
}

void CarbohydrateUtil::_traverse_graph(const unsigned int& idx, std::string& full_form_descriptor, std::string& short_form_descriptor, int& max_length)
{
       full_form_descriptor.clear();
       short_form_descriptor.clear();
       max_length = 0;

       if (_sugar_residue_list[idx]._is_visited) return;

       _sugar_residue_list[idx]._is_visited = true;

       _reassign_anomeric_anomeric_linkage_list(idx);
       if (!_is_cyclic) _reassign_reducing_linkage_list(idx);

       if (_sugar_residue_list[idx]._non_reducing_linkage_list.size() > 1) {
            _reorder_non_reducing_linkage_list(idx);
       }

       std::string sub_full_descriptor = "";
       std::string sub_short_descriptor = "";
       int sub_max_length = 0;

       //key: _non_reducing_linkage_list index
       //value.first: branch length
       //value.second.first: branch full descriptor 
       //value.second.second: branch short descriptor 
       std::map<unsigned int, std::pair<int, std::pair<std::string, std::string> > > descriptor_length_map;
       descriptor_length_map.clear();
       for (unsigned int i = 0; i < _sugar_residue_list[idx]._non_reducing_linkage_list.size(); ++i) {
            _traverse_graph(_sugar_residue_list[idx]._non_reducing_linkage_list[i]._linked_residue_index, sub_full_descriptor,
                            sub_short_descriptor, sub_max_length);
            if (sub_full_descriptor.empty() || sub_short_descriptor.empty() || (sub_max_length == 0)) continue;

            if (sub_max_length > max_length) max_length = sub_max_length;
            descriptor_length_map.insert(std::make_pair(i, std::make_pair(sub_max_length, std::make_pair(sub_full_descriptor, sub_short_descriptor))));
       }
       if (!descriptor_length_map.empty()) {
            std::map<unsigned int, std::pair<std::string, std::string> > tmp_map;

            // key: branch length
            // value.first: _non_reducing_linkage_list index
            // value.second.first branch full descriptor
            // value.second.second: branch short descriptor
            std::map<int, std::map<unsigned int, std::pair<std::string, std::string> > > index_map;
            index_map.clear();
            for (std::map<unsigned int, std::pair<int, std::pair<std::string, std::string> > >::const_iterator
                 dmpos = descriptor_length_map.begin(); dmpos != descriptor_length_map.end(); ++dmpos) {
                 std::map<int, std::map<unsigned int, std::pair<std::string, std::string> > >::iterator impos = index_map.find(dmpos->second.first);
                 if (impos != index_map.end()) {
                      impos->second.insert(std::make_pair(dmpos->first, std::make_pair(dmpos->second.second.first, dmpos->second.second.second)));
                 } else {
                      tmp_map.clear();
                      tmp_map.insert(std::make_pair(dmpos->first, std::make_pair(dmpos->second.second.first, dmpos->second.second.second)));
                      index_map.insert(std::make_pair(dmpos->second.first, tmp_map));
                 }
            }
            std::vector<std::pair<std::string, std::string> > sub_descriptor_list;
            sub_descriptor_list.clear();
            for (std::map<int, std::map<unsigned int, std::pair<std::string, std::string> > >::const_reverse_iterator
                 rmpos = index_map.rbegin(); rmpos != index_map.rend(); ++rmpos) {
                 for (std::map<unsigned int, std::pair<std::string, std::string> >::const_iterator
                      mpos = rmpos->second.begin(); mpos != rmpos->second.end(); ++mpos) {
                      sub_descriptor_list.push_back(mpos->second);
                 }
            }
            for (unsigned int i = 0; i < sub_descriptor_list.size(); ++i) {
                 if (i) {
                      unsigned int size = sub_descriptor_list[i].first.size() - 1;
                      if (sub_descriptor_list[i].first[size] == '-')
                           full_form_descriptor += "[" + sub_descriptor_list[i].first.substr(0, size) + "]";
                      else full_form_descriptor += "[" + sub_descriptor_list[i].first + "]";
                      short_form_descriptor += "[" + sub_descriptor_list[i].second + "]";
                 } else {
                      full_form_descriptor += sub_descriptor_list[i].first;
                      short_form_descriptor += sub_descriptor_list[i].second;
                 }
            }
       }

       _descriptor_id_set.insert(_sugar_residue_list[idx]._id);

       full_form_descriptor += "{" + _sugar_residue_list[idx]._id + "}";
       short_form_descriptor += "{" + _sugar_residue_list[idx]._id + "}";
       if (_sugar_residue_list[idx]._is_root) {
            unsigned int number = _get_number_part_from_atom_name(_sugar_residue_list[idx]._name, _sugar_residue_list[idx]._anomeric_carbon);
            if (number > 0) {
                 std::set<unsigned int> connected_numbers;
                 connected_numbers.clear();
                 for (std::vector<SugarLinkage>::const_iterator vpos = _sugar_residue_list[idx]._reducing_linkage_list.begin();
                      vpos != _sugar_residue_list[idx]._reducing_linkage_list.end(); ++vpos) {
                      connected_numbers.insert(vpos->_current_number);
                 }
                 for (std::vector<SugarLinkage>::const_iterator vpos = _sugar_residue_list[idx]._non_reducing_linkage_list.begin();
                      vpos != _sugar_residue_list[idx]._non_reducing_linkage_list.end(); ++vpos) {
                      connected_numbers.insert(vpos->_current_number);
                 }
                 for (std::vector<SugarLinkage>::const_iterator vpos = _sugar_residue_list[idx]._anomeric_anomeric_linkage_list.begin();
                      vpos != _sugar_residue_list[idx]._anomeric_anomeric_linkage_list.end(); ++vpos) {
                      connected_numbers.insert(vpos->_current_number);
                 }
                 if (connected_numbers.find(number) == connected_numbers.end()) {
                      short_form_descriptor += String::IntToString(number) + "-";
                      if (_sugar_residue_list[idx]._attached_protein_residue == NULL) {
                           std::string terminal = _get_reducing_terminal(_sugar_residue_list[idx]._name, _sugar_residue_list[idx]._anomeric_carbon);
                           if (!terminal.empty()) short_form_descriptor += terminal;
                      }
                 }
            }
       } else if (!_sugar_residue_list[idx]._reducing_linkage_list.empty()) {
            full_form_descriptor += "-(" + _sugar_residue_list[idx]._reducing_linkage_list[0]._linkage_type + ")-";
            short_form_descriptor += _sugar_residue_list[idx]._reducing_linkage_list[0]._linkage_type;
       }
       max_length += 1;
}

void CarbohydrateUtil::_reindex_residues(const std::string& descriptor)
{
       std::vector<std::string> sub_branch_list;
       sub_branch_list.clear();

       std::string main_branch = "";
       std::string sub_branch = "";
       int bracket = 0;

       for (unsigned int i = 0; i < descriptor.size(); ++i) {
            if (descriptor[i] == '[') {
                 bracket++;
                 if (bracket == 1) continue;
            } else if (descriptor[i] == ']') {
                 bracket--;
                 if ((bracket == 0) && !sub_branch.empty()) {
                      sub_branch_list.push_back(sub_branch);
                      sub_branch.clear();
                      continue;
                 }
            }

            if (bracket == 0)
                 main_branch += descriptor[i];
            else sub_branch += descriptor[i];
       }

       std::vector<std::string> residue_list;
       residue_list.clear();
       std::string residue_id = "";

       bracket = 0;
       for (unsigned int i = 0; i < main_branch.size(); ++i) {
            if (main_branch[i] == '{') {
                 bracket++;
            } else if (main_branch[i] == '}') {
                 bracket--;
                 if (!residue_id.empty()) {
                      residue_list.push_back(residue_id);
                      residue_id.clear();
                 }
            } else if (bracket) residue_id += main_branch[i];
       }

       for (std::vector<std::string>::const_reverse_iterator vpos = residue_list.rbegin(); vpos != residue_list.rend(); ++vpos) {
            std::map<std::string, unsigned int>::const_iterator mpos = _sugar_residue_index_map.find(*vpos);
            if (mpos != _sugar_residue_index_map.end()) {
                 _sugar_residue_list[mpos->second]._new_index = _global_index;
                 _global_index += _sugar_residue_list[mpos->second]._residues.size();
            }
       }

       for (std::vector<std::string>::const_iterator vpos = sub_branch_list.begin(); vpos != sub_branch_list.end(); ++vpos) {
            _reindex_residues(*vpos);
       }
}

void CarbohydrateUtil::_generate_new_index_map()
{
       _new_index_map.clear();
       std::string error = "";
       std::string missing_residues = "";
       for (std::vector<SugarResidue>::const_iterator vpos = _sugar_residue_list.begin(); vpos != _sugar_residue_list.end(); ++vpos) {
            if (vpos->_new_index == 0) {
                 if (!missing_residues.empty()) missing_residues += ", ";
                 missing_residues += vpos->_id;
                 continue;
            }
            for (unsigned int i = 0; i < vpos->_residues.size(); ++i) {
                 int index = vpos->_new_index + i;
                 std::string id = vpos->_residues[i]->pdb_chnid() + "_" + vpos->_residues[i]->ResName() + "_"
                                + vpos->_residues[i]->pdb_res_no() + "_" +  vpos->_residues[i]->ins_code();

                 std::map<int, std::string>::const_iterator mpos = _new_index_map.find(index);
                 if (mpos != _new_index_map.end()) {
                      if (!error.empty()) error += "\n";
                      error += "Residues '" + mpos->second + "' and '" + id + "' were assigned to same index number.";
                 }
                 _new_index_map.insert(std::make_pair(index, id));
            }
       }
       if (!error.empty() || !missing_residues.empty()) {
            _new_index_map.clear();

            if (!error.empty()) error += "\n";
            if (!missing_residues.empty()) {
                 error += "Missing the following residue(s) ( " + missing_residues + " ) in numbering index list.\n";
            }
            _logIo->message(error.c_str());
            _error_flag = true;
       }
}

void CarbohydrateUtil::_reassign_anomeric_anomeric_linkage_list(const unsigned int& idx)
{
       if (_sugar_residue_list[idx]._anomeric_anomeric_linkage_list.empty()) return;

       std::vector<SugarLinkage> anomeric_anomeric_linkage_list;
       for (std::vector<SugarLinkage>::const_iterator ipos = _sugar_residue_list[idx]._anomeric_anomeric_linkage_list.begin();
            ipos != _sugar_residue_list[idx]._anomeric_anomeric_linkage_list.end(); ++ipos) {
            // move current _anomeric_anomeric_linkage_list into non reducing linkage
            _sugar_residue_list[idx]._non_reducing_linkage_list.push_back(*ipos);

            unsigned int next_idx = ipos->_linked_residue_index;
            if (_sugar_residue_list[next_idx]._is_visited) continue;

            anomeric_anomeric_linkage_list.clear();
            for (std::vector<SugarLinkage>::const_iterator npos = _sugar_residue_list[next_idx]._anomeric_anomeric_linkage_list.begin();
                 npos != _sugar_residue_list[next_idx]._anomeric_anomeric_linkage_list.end(); ++npos) {          
                 if (npos->_linked_residue_index == idx)
                      // move linked _anomeric_anomeric_linkage_list into reducing linkage
                      _sugar_residue_list[next_idx]._reducing_linkage_list.push_back(*npos);
                 else anomeric_anomeric_linkage_list.push_back(*npos);
            }
            _sugar_residue_list[next_idx]._anomeric_anomeric_linkage_list = anomeric_anomeric_linkage_list;
       }
 
       _sugar_residue_list[idx]._anomeric_anomeric_linkage_list.clear();
}

void CarbohydrateUtil::_reassign_reducing_linkage_list(const unsigned int& idx)
{
       if (_sugar_residue_list[idx]._reducing_linkage_list.empty()) return;

       std::vector<SugarLinkage> reducing_linkage_list, non_reducing_linkage_list;
       reducing_linkage_list.clear();
       for (std::vector<SugarLinkage>::iterator ipos = _sugar_residue_list[idx]._reducing_linkage_list.begin();
            ipos != _sugar_residue_list[idx]._reducing_linkage_list.end(); ++ipos) {
            unsigned int next_idx = ipos->_linked_residue_index;
            if (_sugar_residue_list[next_idx]._is_visited) {
                 // keep the true reducing linkage
                 reducing_linkage_list.push_back(*ipos);;
            } else {
                 // move the false reducing linkage into non reducing linkage
                 std::string linkage_type = ipos->_linkage_type;
                 ipos->_linkage_type = ipos->_inverse_linkage_type;
                 ipos->_inverse_linkage_type = linkage_type;
                 _sugar_residue_list[idx]._non_reducing_linkage_list.push_back(*ipos);

                 non_reducing_linkage_list.clear();
                 for (std::vector<SugarLinkage>::iterator npos = _sugar_residue_list[next_idx]._non_reducing_linkage_list.begin();
                      npos != _sugar_residue_list[next_idx]._non_reducing_linkage_list.end(); ++npos) {          
                      if (npos->_linked_residue_index == idx) {
                           // move linked _non_reducing_linkage_list into reducing linkage
                           linkage_type = npos->_linkage_type;
                           npos->_linkage_type = npos->_inverse_linkage_type;
                           npos->_inverse_linkage_type = linkage_type;
                           _sugar_residue_list[next_idx]._reducing_linkage_list.push_back(*npos);
                      } else non_reducing_linkage_list.push_back(*npos);
                 }
                 _sugar_residue_list[next_idx]._non_reducing_linkage_list = non_reducing_linkage_list;
            }
       }
 
       _sugar_residue_list[idx]._reducing_linkage_list = reducing_linkage_list;
}

void CarbohydrateUtil::_reorder_non_reducing_linkage_list(const unsigned int& idx)
{
       std::map<unsigned int, SugarLinkage> sugarLinkageMap;
       sugarLinkageMap.clear();

       for (std::vector<SugarLinkage>::const_iterator ipos = _sugar_residue_list[idx]._non_reducing_linkage_list.begin();
            ipos != _sugar_residue_list[idx]._non_reducing_linkage_list.end(); ++ipos) {
            sugarLinkageMap.insert(std::make_pair(ipos->_current_number, *ipos));
       }

       _sugar_residue_list[idx]._non_reducing_linkage_list.clear();
       for (std::map<unsigned int, SugarLinkage>::const_iterator mpos = sugarLinkageMap.begin(); mpos != sugarLinkageMap.end(); ++mpos) {
            _sugar_residue_list[idx]._non_reducing_linkage_list.push_back(mpos->second);
       }
}

std::string CarbohydrateUtil::_get_reducing_terminal(const std::string& ccd_id, const std::string& anomeric_carbon)
{
       try {
            const ConnectFormat& drug = _ccDic->find_drug(ccd_id);

            std::string pattern = "";
            const std::vector<std::map<std::string, std::string> >& tree_list = drug.get_sugar_side_atom_tree_list(anomeric_carbon);
            for (std::vector<std::map<std::string, std::string> >::const_iterator vpos = tree_list.begin(); vpos != tree_list.end(); ++vpos) {
                 for (std::map<std::string, std::string>::const_iterator mpos = vpos->begin(); mpos != vpos->end(); ++mpos) {
                      if (mpos->second == "H") continue;
                      pattern += mpos->second;
                 }
            } 

            if ((pattern == "O") || (pattern == "OH"))
                 return "ROH";
            else if ((pattern == "OC") || (pattern == "OCHHH"))
                 return "OME";
            else if ((pattern == "OCCCC") || (pattern == "OCCCCHHHHHHHHH"))
                 return "TBT";
            else return "";
       } catch (const std::exception& exc) {}

       return "";
}
