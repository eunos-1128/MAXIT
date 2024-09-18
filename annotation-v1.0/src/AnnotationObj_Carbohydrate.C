/*
FILE:     AnnotationObj_Carbohydrate.export.C
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
#include <sys/stat.h>
#include <ctype.h>

// #include "algorithm-util.h"
// #include "algorithm-util.C"
#include "AnnotationObj.h"
#include "CarbohydrateUtil.h"
#include "CompositeIndex.h"
#include "SeqCodeUtil.h"
#include "utillib.h"

void AnnotationObj::_annotate_branch_polymers()
{
       if (_branch_annotated_flag || !_find_sugar_chain_flag) return;

       bool need_carbohydrate_annotation = false;
       for (unsigned int i = 0; i < _molecules.size(); ++i) {
            if (!_molecules[i]->has_sugar_entity()) continue;
            need_carbohydrate_annotation = true;
            break;
       }
       _branch_annotated_flag = true;
       if (!need_carbohydrate_annotation) return;

       std::map<unsigned int, std::list<_LINK> > _modelLinkMap;
       _modelLinkMap.clear();
       if (_molecules.size() > 1) {
            _LINK link;
            link.SymOP_1.clear();
            link.SymOP_2.clear();
            link.details.clear();
            link.type.clear();
            link.dist.clear();
            link.bondtype.clear();
            link.leaving_flag.clear();
            link.pdbx_role.clear();

            std::list<_LINK> model_links;
            for (unsigned int i = 1; i < _molecules.size(); ++i) {
                 if (!_molecules[i]->has_sugar_entity()) continue;

                 model_links.clear();
                 for (std::list<_LINK>::const_iterator ptr = _links.begin(); ptr != _links.end(); ++ptr) {
                      if ((!ptr->SymOP_1.empty() && ptr->SymOP_1 != "1_555") || (!ptr->SymOP_2.empty() && ptr->SymOP_2 != "1_555")) continue;

                      RCSB::Residue* res1 = _molecules[i]->find_pdb_residue(ptr->fstAtom->pdb_chnid(), ptr->fstAtom->pdb_resnam(),
                                                                            ptr->fstAtom->pdb_resnum(), ptr->fstAtom->ins_code());
                      if (!res1) continue;
                      RCSB::Atom* atom1 = res1->find_atom(ptr->fstAtom->pdb_atmnam());
                      if (!atom1) continue;

                      RCSB::Residue* res2 = _molecules[i]->find_pdb_residue(ptr->sndAtom->pdb_chnid(), ptr->sndAtom->pdb_resnam(),
                                                                            ptr->sndAtom->pdb_resnum(), ptr->sndAtom->ins_code());
                      if (!res2) continue;
                      RCSB::Atom* atom2 = res2->find_atom(ptr->sndAtom->pdb_atmnam());
                      if (!atom2) continue;
            
                      link.fstAtom = atom1;
                      link.sndAtom = atom2;
                      model_links.push_back(link);
                 }
                 _modelLinkMap.insert(std::make_pair(i, model_links));
            }
       }

       _get_glyco_site_linkage_info();

       int sequential_count = 0;
       for (unsigned int i = 0; i < _molecules.size(); ++i) {
            if (!_molecules[i]->has_sugar_entity()) continue;

            RCSB::Chain* chain = _molecules[i]->GetFirstChain();
            while (chain) {
                 if ((chain->chain_type() == "ATOMS") /* && chain->need_carbohydrate_annotation() */) {
                      RCSB::Residue* glycoRes = _get_glyco_site_residue(_molecules[i], chain);
                      if (i == 0) _get_carb_annotation(glycoRes, chain, _links, sequential_count);
                      else {
                           std::map<unsigned int, std::list<_LINK> >::iterator mpos = _modelLinkMap.find(i);
                           _get_carb_annotation(glycoRes, chain, mpos->second, sequential_count);
                      }
                 }
                 chain = _molecules[i]->GetNextChain();
            }
       }

       _run_assign_carbohydrate_chain_id();
}

void AnnotationObj::_get_glyco_site_linkage_info()
{
       _glyco_site_linkage_info.clear();

       if (_links.empty()) return;

       std::string idx;
       std::vector<std::string> data;
       for (std::list<_LINK>::const_iterator pos = _links.begin(); pos != _links.end(); ++pos) {
            if (!pos->SymOP_1.empty() || !pos->SymOP_2.empty()) continue;

            idx.clear();
            data.clear();
            std::string glycosylation_type = _find_glycosylation_type(pos->fstAtom, pos->sndAtom);
            if (!glycosylation_type.empty()) {
                 if (SeqCodeUtil::is_standard_aa_residue(pos->sndAtom->pdb_resnam())) {
                      idx = CompositeIndex::getIndex(pos->fstAtom->pdb_chnid(), pos->fstAtom->pdb_resnam(), pos->fstAtom->pdb_resnum(),
                                                     pos->fstAtom->ins_code());
                      data.push_back(pos->sndAtom->pdb_chnid());
                      data.push_back(pos->sndAtom->pdb_resnam());
                      data.push_back(pos->sndAtom->pdb_resnum());
                      data.push_back(pos->sndAtom->ins_code());
                 } else{
                      idx = CompositeIndex::getIndex(pos->sndAtom->pdb_chnid(), pos->sndAtom->pdb_resnam(), pos->sndAtom->pdb_resnum(),
                                                     pos->sndAtom->ins_code());
                      data.push_back(pos->fstAtom->pdb_chnid());
                      data.push_back(pos->fstAtom->pdb_resnam());
                      data.push_back(pos->fstAtom->pdb_resnum());
                      data.push_back(pos->fstAtom->ins_code());
                 }
            }
            if (!idx.empty() && !data.empty()) _glyco_site_linkage_info.insert(std::make_pair(idx, data));
       }
}

RCSB::Residue* AnnotationObj::_get_glyco_site_residue(RCSB::Molecule* mol, RCSB::Chain* chain)
{
       if (_glyco_site_linkage_info.empty()) return NULL;

       // Only select first conformer
       std::vector<RCSB::Residue*> resList;
       chain->GetFirstResidueList(resList);
       while (!resList.empty()) {
            std::string idx = CompositeIndex::getIndex(resList[0]->pdb_chnid(), resList[0]->ResName(), resList[0]->pdb_res_no(), resList[0]->ins_code());
            std::map<std::string, std::vector<std::string> >::const_iterator mpos = _glyco_site_linkage_info.find(idx);
            if (mpos != _glyco_site_linkage_info.end()) {
                 RCSB::Residue* glycoRes = mol->find_pdb_residue(mpos->second[0], mpos->second[1], mpos->second[2], mpos->second[3]);
                 if (glycoRes) return glycoRes;
            }
            chain->GetNextResidueList(resList);
       }

       return NULL;
}

void AnnotationObj::_get_carb_annotation(RCSB::Residue* glycoRes, RCSB::Chain* chain, std::list<_LINK>& model_related_links, int& sequential_count)
{
       // Re-name pdb chain ID to take care of 2-letter pdb chain ID case
       // key: residue's orginal pdb chain ID
       // value: re-assign unique one letter pdb chain ID
       std::map<std::string, std::string> pdb_chnid_map;
       pdb_chnid_map.clear();

       // re-assign unique one letter pdb chain ID
       std::set<std::string> used_pdb_chnid_set;
       used_pdb_chnid_set.clear();

       std::vector<std::pair<std::vector<RCSB::Residue*>, std::pair<std::string, std::string> > > sugar_residue_list;
       sugar_residue_list.clear();

       std::vector<RCSB::Residue*> tmp_res_vec;

       std::string removeName = "";
       std::string resName = "Unknown";
       std::string re_assign_pdb_chnid;
       if (glycoRes) {
            re_assign_pdb_chnid.clear();
            std::map<std::string, std::string>::const_iterator mpos = pdb_chnid_map.find(glycoRes->pdb_chnid());
            if (mpos != pdb_chnid_map.end()) re_assign_pdb_chnid = mpos->second;
            else {
                 re_assign_pdb_chnid = get_next_available_pdb_chain_id(used_pdb_chnid_set);
                 pdb_chnid_map.insert(std::make_pair(glycoRes->pdb_chnid(), re_assign_pdb_chnid));
            }
            tmp_res_vec.clear();
            tmp_res_vec.push_back(glycoRes);
            sugar_residue_list.push_back(std::make_pair(tmp_res_vec, std::make_pair("ATOMP", re_assign_pdb_chnid)));
            resName = glycoRes->ResName();
       }
       removeName = "[" + resName + "]";

       std::set<int> overlapExcludedSet;
       overlapExcludedSet.clear();

       // Only select first conformer
       unsigned int residue_count = 0;
       std::vector<RCSB::Residue*> resList;
       chain->GetFirstResidueList(resList);
       while (!resList.empty()) {
            residue_count++;

            re_assign_pdb_chnid.clear();
            std::map<std::string, std::string>::const_iterator pmpos = pdb_chnid_map.find(resList[0]->pdb_chnid());
            if (pmpos != pdb_chnid_map.end()) re_assign_pdb_chnid = pmpos->second;
            else {
                 re_assign_pdb_chnid = get_next_available_pdb_chain_id(used_pdb_chnid_set);
                 pdb_chnid_map.insert(std::make_pair(resList[0]->pdb_chnid(), re_assign_pdb_chnid));
            }

            sugar_residue_list.push_back(std::make_pair(resList, std::make_pair("HETAIN", re_assign_pdb_chnid)));

            chain->GetNextResidueList(resList);
       }

       if (residue_count < 2) return;

       CarbohydrateUtil carbUtil;
       carbUtil.setLog(_logIo);
       carbUtil.setCCDic(_ccDic);
       carbUtil.setPDBID(_LowercasePdbId);

       for (std::vector<std::pair<std::vector<RCSB::Residue*>, std::pair<std::string, std::string> > >::const_iterator
            vpos = sugar_residue_list.begin(); vpos != sugar_residue_list.end(); ++vpos) {
            if (vpos->second.first == "ATOMP")
                 carbUtil.addResidue(vpos->first, vpos->second.first);
            else carbUtil.addResidue(vpos->first);
       }

       carbUtil.addLinks(model_related_links);
       carbUtil.buildOligoSaccharide();

       std::vector<std::map<std::string, std::string> > linear_descriptors;
       linear_descriptors.clear();

       if (carbUtil.getStatus()) {
            std::string entityName = carbUtil.getEntityDescriptor();
            if (!entityName.empty()) chain->set_descriptor(entityName);

            std::string descriptor = carbUtil.getCondensedDescriptor();
            if (!descriptor.empty()) {
                 std::map<std::string, std::string> tmp_map;
                 tmp_map.clear();
                 tmp_map["type"] = "Glycam Condensed Sequence";
                 tmp_map["descriptor"] = descriptor;
                 tmp_map["program"] = "GMML";
                 tmp_map["program_version"] = "1.0";
                 linear_descriptors.push_back(tmp_map);
            }

            const std::map<int, std::string>& reorder_index_mapping = carbUtil.getReorderIndexMapping();
            if (!reorder_index_mapping.empty()) chain->reorder_residues(reorder_index_mapping);

            const std::list<_LINK>& additional_links = carbUtil.getAdditionalLinks();
            if (!additional_links.empty()) {
                 for (std::list<_LINK>::const_iterator lpos = additional_links.begin(); lpos != additional_links.end(); ++lpos) {
                      model_related_links.push_back(*lpos);
                 }
            }

            std::set<std::string> bad_link_set;
            bad_link_set.clear();

            const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& bad_links = carbUtil.getBadLinks();
            for (std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >::const_iterator lpos = bad_links.begin(); lpos != bad_links.end(); ++lpos) {
                 bad_link_set.insert(lpos->first->getAtomAllIndex(false) + "-" + lpos->second->getAtomAllIndex(false));
                 bad_link_set.insert(lpos->second->getAtomAllIndex(false) + "-" + lpos->first->getAtomAllIndex(false));
            }

            for (std::list<_LINK>::iterator ptr = model_related_links.begin(); ptr != model_related_links.end(); ++ptr) {
                 if ((!ptr->SymOP_1.empty() && ptr->SymOP_1 != "1_555") || (!ptr->SymOP_2.empty() && ptr->SymOP_2 != "1_555")) continue;

                 if (bad_link_set.find(ptr->fstAtom->getAtomAllIndex(false) + "-" + ptr->sndAtom->getAtomAllIndex(false)) != bad_link_set.end()) {
                      ptr = model_related_links.erase(ptr);
                      --ptr;
                 } else if ((ptr->fstAtom->pdb_chnid() == ptr->sndAtom->pdb_chnid()) && (ptr->fstAtom->pdb_resnam() == ptr->sndAtom->pdb_resnam()) &&
                            (ptr->fstAtom->pdb_resnum() == ptr->sndAtom->pdb_resnum()) && (ptr->fstAtom->ins_code() == ptr->sndAtom->ins_code())) {
                      ptr = model_related_links.erase(ptr);
                      --ptr;
                 }
            }

            std::string entity_key = carbUtil.getEntityKey();
            if (!entity_key.empty()) chain->set_entity_key(entity_key);
       }

       _get_branch_entity_links(chain, model_related_links, carbUtil.getBadLinks());

       if (!linear_descriptors.empty()) chain->set_linear_descriptors(linear_descriptors);
}

void AnnotationObj::_run_assign_carbohydrate_chain_id()
{
       std::map<std::string, std::string> sugar_chain_id_mapping;

       bool first = true;
       for (std::vector<RCSB::Molecule*>::const_iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            sugar_chain_id_mapping.clear();
            if ((*mpos)->has_sugar_entity()) (*mpos)->assign_carbohydrate_chain_id(sugar_chain_id_mapping);
            if (first && !sugar_chain_id_mapping.empty()) {
                 std::map<std::string, std::vector<std::string> > sugar_old_new_chain_id_mapping;
                 sugar_old_new_chain_id_mapping.clear();
                 for (std::map<std::string, std::string>::const_iterator mpos = sugar_chain_id_mapping.begin(); mpos != sugar_chain_id_mapping.end(); ++mpos) {
                      std::map<std::string, std::vector<std::string> >::iterator smpos = sugar_old_new_chain_id_mapping.find(mpos->second);
                      if (smpos != sugar_old_new_chain_id_mapping.end()) {
                           smpos->second.push_back(mpos->first);
                      } else {
                           std::vector<std::string> t_vec;
                           t_vec.clear();
                           t_vec.push_back(mpos->first);
                           sugar_old_new_chain_id_mapping.insert(std::make_pair(mpos->second, t_vec));
                      }
                 }
                 if (!sugar_old_new_chain_id_mapping.empty()) {
                      for (unsigned int i = 0; i < _assemblies.size(); ++i) {
                           _assemblies[i].AddSugarChains(sugar_old_new_chain_id_mapping);
                      }
                 }
            }
            first = false;
       }
}
