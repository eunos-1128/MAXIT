/*
FILE:     Ndb2Pdb_Remark_400_group.C
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

#include "Maxit.h"
#include "NdbToken.h"
#include "Remark_define.h"
#include "Remark_extern.h"
#include "utillib.h"

typedef struct {
       std::string id;
       std::string prd_id;
       std::string name;
       std::string chainids;
       std::vector<std::string> components;
       std::string description;
} _Remark_400_Group;

static void get_polymer_prd_ids(ISTable *t, const std::map<std::string, RCSB::Chain*>& asym_id_chain, std::set<std::string>& prd_ids);
static void get_group_info(ISTable *t, const std::set<std::string>& polymer_prds, std::vector<_Remark_400_Group>& remark_groups,
                           std::vector<std::string>& general_remarks);
static void get_prd_pdb_chnid_mapping(ISTable *t, const std::set<std::string>& polymer_prds, const std::map<std::string, RCSB::Chain*>& asym_id_chain,
                                      std::map<std::string, std::set<std::string> >& prd_pdb_chnid_mapping,
                                      std::map<std::string, std::vector<RCSB::Chain*> >& prd_instance_mapping);
static void get_prd_component_mapping(const std::map<std::string, std::vector<RCSB::Chain*> >& prd_instance_mapping,
                                      std::map<std::string, std::vector<std::string> >& prd_component_mapping);
static void update_remark_groups(const std::map<std::string, int>& polymer_chain_order, const std::map<std::string, std::set<std::string> >&
                                 prd_pdb_chnid_mapping, const std::map<std::string, std::vector<std::string> >& prd_component_mapping,
                                 std::vector<_Remark_400_Group>& remark_groups);

void Maxit::_ndb_to_pdb_get_remark_400_group()
{
       _non_polymer_prd_info.clear();

       if (!_CifObj || _molecules.empty()) return;

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);
       ISTable *t1 = getTablePtr(_cifblock, "pdbx_molecule");
       if (!t1) return;
       ISTable *t2 = getTablePtr(_cifblock, "pdbx_molecule_features");
       if (!t2) return;

       // key: PDB chain ID
       // value: order
       std::map<std::string, int> polymer_chain_order;

       // key: Asym ID
       // value: chain pointer
       std::map<std::string, RCSB::Chain*> asym_id_chain;

       _ndb_to_pdb_get_polymer_chain_info(polymer_chain_order, asym_id_chain);

       std::set<std::string> polymer_prds;
       get_polymer_prd_ids(t1, asym_id_chain, polymer_prds);
       _ndb_to_pdb_get_non_polymer_prd_info(t1, t2, asym_id_chain, polymer_prds);
       if (polymer_prds.empty()) return;

       std::vector<_Remark_400_Group> remark_groups;
       std::vector<std::string> general_remarks;

       get_group_info(t2, polymer_prds, remark_groups, general_remarks);
       if (remark_groups.empty()) return;

       // key: prd id
       // value: PDB chain IDs
       std::map<std::string, std::set<std::string> > prd_pdb_chnid_mapping;
       // key: prd_id
       // value: instance components (pick the instance has most components presented)
       std::map<std::string, std::vector<RCSB::Chain*> > prd_instance_mapping;
       get_prd_pdb_chnid_mapping(t1, polymer_prds, asym_id_chain, prd_pdb_chnid_mapping, prd_instance_mapping);

       // key: prd id
       // value: list of component description
       std::map<std::string, std::vector<std::string> > prd_component_mapping;
       get_prd_component_mapping(prd_instance_mapping, prd_component_mapping);

       update_remark_groups(polymer_chain_order, prd_pdb_chnid_mapping, prd_component_mapping, remark_groups);

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
       int width = ndbformat.FieldList[2].FieldWidth;

       std::vector<std::string> remark_array;

       bool found_400 = false;
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("REMARK");
       if (ppos != _pdb_records.end()) {
            const std::vector<std::string>& Field = ppos->second.back();
            if (Field[1] == "400") found_400 = true;
       }

       if (!found_400) _addNewRemark(400, "COMPOUND");
 
       for (std::vector<std::string>::const_iterator pos = general_remarks.begin(); pos != general_remarks.end(); ++pos) {
            _addNewRemark(400, "");
            get_text_array_from_block(remark_array, *pos, width);
            _ndb_to_pdb_add_remark(400, remark_array);
       }

       for (std::vector<_Remark_400_Group>::const_iterator gpos = remark_groups.begin(); gpos != remark_groups.end(); ++gpos) {
            _addNewRemark(400, "");
            _addNewRemark(400, " GROUP: " + gpos->id);

            if (gpos->name.empty())
                 _addNewRemark(400, "  NAME: NULL");
            else {
                 get_text_array_from_block_with_prefix(remark_array, gpos->name, "  NAME: ", width - 8);
                 _ndb_to_pdb_add_remark(400, remark_array);
            }

            if (gpos->chainids.empty())
                 _addNewRemark(400, "  CHAIN: NULL");
            else {
                 get_text_array_from_block_with_prefix(remark_array, gpos->chainids, "  CHAIN: ", width - 9);
                 _ndb_to_pdb_add_remark(400, remark_array);
            }

            int count = 0;
            for (std::vector<std::string>::const_iterator vpos = gpos->components.begin(); vpos != gpos->components.end(); ++vpos) { 
                 count++;
                 std::string prefix = "  COMPONENT_" + String::IntToString(count) + ": ";
                 if (vpos->empty())
                      _addNewRemark(400, prefix + "NULL");
                 else {
                      get_text_array_from_block_with_prefix(remark_array, *vpos, prefix, width - 15);
                      _ndb_to_pdb_add_remark(400, remark_array);
                 }
            }

            if (gpos->description.empty())
                 _addNewRemark(400, "  DESCRIPTION: NULL");
            else {
                 get_text_array_from_block_with_prefix(remark_array, gpos->description, "  DESCRIPTION: ", width - 15);
                 _ndb_to_pdb_add_remark(400, remark_array);
            }
       }
}

void Maxit::_ndb_to_pdb_get_polymer_chain_info(std::map<std::string, int>& chain_order, std::map<std::string, RCSB::Chain*>& asym_id_chain)
{
       chain_order.clear();
       asym_id_chain.clear();

       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            if (chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP") {
                 chain_order.insert(std::make_pair(chain->PDB_ChainID(), (int) chain_order.size()));
            }
            asym_id_chain.insert(std::make_pair(chain->ChainID(), chain));
            chain = _molecules[0]->GetNextChain();
       }
}

void Maxit::_ndb_to_pdb_get_non_polymer_prd_info(ISTable *instances, ISTable *features, const std::map<std::string, RCSB::Chain*>& asym_id_chain,
                                                 const std::set<std::string>& polymer_prds)
{
       std::set<std::string> polymer_instances;
       polymer_instances.clear();

       std::string cs;
       int rowNo = instances->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs, instances, i, "prd_id");
            if (cs.empty()) continue;
            if (polymer_prds.find(cs) != polymer_prds.end()) {
                 get_value_clean(cs, instances, i, "instance_id");
                 polymer_instances.insert(cs);
            }
       }

       std::vector<RCSB::Residue*> residues;
       std::map<std::string, std::string> het_prd_id_mapping;
       het_prd_id_mapping.clear();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs, instances, i, "instance_id");
            if (cs.empty()) continue;
            if (polymer_instances.find(cs) != polymer_instances.end()) continue;
            get_value_clean(cs, instances, i, "asym_id");
            if (cs.empty()) continue;
            std::map<std::string, RCSB::Chain*>::const_iterator mpos = asym_id_chain.find(cs);
            if (mpos == asym_id_chain.end()) continue;
            get_value_clean(cs, instances, i, "prd_id");
            if (cs.empty()) continue;
            mpos->second->GetFirstResidueList(residues);
            het_prd_id_mapping.insert(std::make_pair(residues[0]->ResName(), cs));
       }
       if (het_prd_id_mapping.empty()) return;

       std::string prd_id, type, clss, details; 
       std::vector<std::string> data;
       // key: prd_id
       // vector[0]: prd_id
       // vector[1]: type + class
       // vector[2]: details
       std::map<std::string, std::vector<std::string> > prd_info;
       prd_info.clear();
       rowNo = features->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(prd_id, features, i, "prd_id");
            if (prd_id.empty()) continue;

            data.clear();
            data.push_back(prd_id);
            get_value_clean(type, features, i, "type");
            get_value_clean(clss, features, i, "class");
            get_value_clean(details, features, i, "details");
            cs = type;
            if (!clss.empty()) {
                 if (!cs.empty()) cs += " ";
                 cs += clss;
            }
            data.push_back(cs);
            data.push_back(details);
            prd_info.insert(std::make_pair(prd_id, data));
       }

       if (prd_info.empty()) return;

       for (std::map<std::string, std::string>::const_iterator pos = het_prd_id_mapping.begin(); pos != het_prd_id_mapping.end(); ++pos) {
            std::map<std::string, std::vector<std::string> >::const_iterator mpos = prd_info.find(pos->second);
            if (mpos == prd_info.end()) continue;
            _non_polymer_prd_info.insert(std::make_pair(pos->first, mpos->second));
       }
}

static void get_polymer_prd_ids(ISTable *t, const std::map<std::string, RCSB::Chain*>& asym_id_chain, std::set<std::string>& prd_ids)
{
       prd_ids.clear();

       std::string cs, cs1;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs, t, i, "prd_id");
            if (cs.empty()) continue;
            get_value_clean(cs1, t, i, "asym_id");
            if (cs1.empty()) continue;
            std::map<std::string, RCSB::Chain*>::const_iterator mpos = asym_id_chain.find(cs1);
            if (mpos == asym_id_chain.end()) continue;
            if (mpos->second->chain_type() != "ATOMP" && mpos->second->chain_type() != "ATOMN") continue;
            prd_ids.insert(cs);
       }
}

static void get_group_info(ISTable *t, const std::set<std::string>& polymer_prds, std::vector<_Remark_400_Group>& remark_groups,
                           std::vector<std::string>& general_remarks)
{
       remark_groups.clear();
       general_remarks.clear();

       _Remark_400_Group group;
       group.id.clear();
       group.prd_id.clear();
       group.name.clear();
       group.chainids.clear();
       group.components.clear();
       group.description.clear();

       std::string cs, cs1;
       int count = 0;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(group.prd_id, t, i, "prd_id");
            if (polymer_prds.find(group.prd_id) == polymer_prds.end()) continue;

            count++;
            group.id = String::IntToString(count);
            get_value_clean_upper(group.name, t, i, "name");
            get_value_clean(group.description, t, i, "details");
            remark_groups.push_back(group);

            get_value_upper(cs, t, i, "type");
            get_value_upper(cs1, t, i, "class");
            std::string str = "THE " + group.name + " IS " + cs + ", A MEMBER OF " + cs1 + " CLASS.";;
            general_remarks.push_back(str);
       }
}

static void get_prd_pdb_chnid_mapping(ISTable *t, const std::set<std::string>& polymer_prds, const std::map<std::string, RCSB::Chain*>& asym_id_chain,
                                      std::map<std::string, std::set<std::string> >& prd_pdb_chnid_mapping,
                                      std::map<std::string, std::vector<RCSB::Chain*> >& prd_instance_mapping)
{
       prd_pdb_chnid_mapping.clear();
       prd_instance_mapping.clear();

       // first key:  prd id
       // second key: instance id
       // vector: Chain(s) related to the instance
       std::map<std::string, std::map<std::string, std::vector<RCSB::Chain*> > > prd_component_id_mapping;
       prd_component_id_mapping.clear();

       std::set<std::string> data_set;
       std::vector<RCSB::Chain*> i_vector;
       std::map<std::string, std::vector<RCSB::Chain*> > t_map;
       std::string prd_id, cs;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value(prd_id, t, i, "prd_id");
            if (polymer_prds.find(prd_id) == polymer_prds.end()) continue;

            get_value(cs, t, i, "asym_id");

            // converted Asym ID to PDB chain ID
            std::map<std::string, RCSB::Chain*>::const_iterator mpos = asym_id_chain.find(cs);
            if (mpos == asym_id_chain.end()) continue;
            cs = mpos->second->PDB_ChainID();

            std::map<std::string, std::set<std::string> >::iterator mpos1 = prd_pdb_chnid_mapping.find(prd_id);
            if (mpos1 == prd_pdb_chnid_mapping.end()) {
                 data_set.clear();
                 data_set.insert(cs);
                 prd_pdb_chnid_mapping.insert(std::make_pair(prd_id, data_set));
            } else mpos1->second.insert(cs);

            get_value(cs, t, i, "instance_id"); 
            std::map<std::string, std::map<std::string, std::vector<RCSB::Chain*> > >::iterator mpos2 = prd_component_id_mapping.find(prd_id);
            if (mpos2 != prd_component_id_mapping.end()) {
                 std::map<std::string, std::vector<RCSB::Chain*> >::iterator mpos3 = mpos2->second.find(cs);
                 if (mpos3 != mpos2->second.end()) mpos3->second.push_back(mpos->second);
                 else {
                      i_vector.clear();
                      i_vector.push_back(mpos->second);
                      mpos2->second.insert(std::make_pair(cs, i_vector));
                 }
            } else {
                 i_vector.clear();
                 i_vector.push_back(mpos->second);
                 t_map.clear();
                 t_map.insert(std::make_pair(cs, i_vector));
                 prd_component_id_mapping.insert(std::make_pair(prd_id, t_map));
            }
       }

       // find largest instance example (has most components presented)
       for (std::map<std::string, std::map<std::string, std::vector<RCSB::Chain*> > >::const_iterator
            mpos = prd_component_id_mapping.begin(); mpos != prd_component_id_mapping.end(); ++mpos) {
            for (std::map<std::string, std::vector<RCSB::Chain*> >::const_iterator mpos1 = mpos->second.begin(); mpos1 != mpos->second.end(); ++mpos1) {
                 std::map<std::string, std::vector<RCSB::Chain*> >::iterator mpos2 = prd_instance_mapping.find(mpos->first);
                 if (mpos2 == prd_instance_mapping.end())
                      prd_instance_mapping.insert(std::make_pair(mpos->first, mpos1->second));
                 else if (mpos1->second.size() > mpos2->second.size())
                      mpos2->second = mpos1->second;
            }
       }
}

static void get_prd_component_mapping(const std::map<std::string, std::vector<RCSB::Chain*> >& prd_instance_mapping, std::map<std::string,
                                      std::vector<std::string> >& prd_component_mapping)
{
       prd_component_mapping.clear();

       std::vector<std::string> data;
       std::map<std::string, int> non_polymer_component_count;

       for (std::map<std::string, std::vector<RCSB::Chain*> >::const_iterator mpos =
            prd_instance_mapping.begin(); mpos != prd_instance_mapping.end(); ++mpos) {
            data.clear();
            non_polymer_component_count.clear();
            for (std::vector<RCSB::Chain*>::const_iterator cpos = mpos->second.begin(); cpos != mpos->second.end(); ++cpos) {
                 RCSB::Chain* chain = *cpos;
                 if (chain->chain_type() == "ATOMP")
                      data.push_back("PEPTIDE LIKE POLYMER");
                 else if (chain->chain_type() == "ATOMN")
                      data.push_back("POLYMER");
                 else {
                      for (unsigned int i = 0; i < chain->SeqLen(); ++i) {
                           std::string cs = chain->SeqRes(i)->Field[0];
                           std::map<std::string, int>::iterator mpos1 = non_polymer_component_count.find(cs);
                           if (mpos1 != non_polymer_component_count.end())
                                mpos1->second++;
                           else non_polymer_component_count.insert(std::make_pair(cs, 1));
                      }
                 }
            }

            for (std::map<std::string, int>::const_iterator mpos1 = non_polymer_component_count.begin();
                 mpos1 != non_polymer_component_count.end(); ++mpos1) {
                 std::string str = "RESIDUE " + mpos1->first;
                 if (mpos1->second > 1) {
                      str += ", " + String::IntToString(mpos1->second) + " COPIES";
                 }
                 data.push_back(str);
            }
            prd_component_mapping.insert(std::make_pair(mpos->first, data));
       }
}

static void update_remark_groups(const std::map<std::string, int>& polymer_chain_order, const std::map<std::string, std::set<std::string> >&
                                 prd_pdb_chnid_mapping, const std::map<std::string, std::vector<std::string> >& prd_component_mapping,
                                 std::vector<_Remark_400_Group>& remark_groups)
{
       std::multimap<int, std::string> order_chains;

       for (std::vector<_Remark_400_Group>::iterator gpos = remark_groups.begin(); gpos != remark_groups.end(); ++gpos) {
            std::map<std::string, std::set<std::string> >::const_iterator mpos = prd_pdb_chnid_mapping.find(gpos->prd_id);
            if (mpos != prd_pdb_chnid_mapping.end()) { // update chainids
                 order_chains.clear();
                 int j = polymer_chain_order.size();
                 for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                      std::map<std::string, int>::const_iterator mpos1 = polymer_chain_order.find(*spos);
                      if (mpos1 != polymer_chain_order.end())
                           order_chains.insert(std::make_pair(mpos1->second, *spos));
                      else order_chains.insert(std::make_pair(j++, *spos));
                 }
                 for (std::multimap<int, std::string>::const_iterator mpos2 = order_chains.begin(); mpos2 != order_chains.end(); ++mpos2) {
                      if (!gpos->chainids.empty()) gpos->chainids += ", ";
                      gpos->chainids += mpos2->second;
                 }
            }
            // update components
            std::map<std::string, std::vector<std::string> >::const_iterator mpos3 = prd_component_mapping.find(gpos->prd_id);
            if (mpos3 != prd_component_mapping.end()) gpos->components = mpos3->second;
       }
}
