/*
FILE:     CNA_Assign_ChainID.C
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
#include "ChainIDNumberAssignment.h"
#include "CompositeIndex.h"

#define N_CYCLE  6

static double _lower_bound[N_CYCLE] = { 1.2, 3.5,  5.0, 10.0, 20.0,  50.0 };
static double _upper_bound[N_CYCLE] = { 3.5, 5.0, 10.0, 20.0, 50.0, 100.0 };

void ChainIDNumberAssignment::_assign_single_chainID(const std::string& pdb_chain_id, RCSB::Molecule* mol)
{
       RCSB::Chain* chain = mol->GetFirstChain();
       while (chain) {
            if (chain->chain_type() != "ATOMS") chain->set_PDB_ChainID_to_residues(pdb_chain_id);
            chain = mol->GetNextChain();
       }
}

void ChainIDNumberAssignment::_assign_multiple_chainIDs(xtal& mycrys, std::vector<RCSB::Chain*>& nonpolymer_chains)
{
       mycrys.markAtoms("polymer");

       std::vector<std::pair<std::string, RCSB::Chain*> > nonpolymer_chain_pairs;
       nonpolymer_chain_pairs.clear();

       std::list<cryatom*> atm_list;
       mycrys.selectAtoms(atm_list, "nonpolymer");
       if (!atm_list.empty()) {
            std::map<int, int> resIndexOrder_mapping;
            _get_index_mapping(resIndexOrder_mapping, nonpolymer_chains);

            // nonpolymer_chains vector index vs. assigned PDB chain ID mapping
            std::map<int, std::string> index_chain_id_mapping;
            index_chain_id_mapping.clear();

            std::list<CONTACT> contacts;
            for (int cycle = 0; cycle < N_CYCLE / 2; ++cycle) {
                 mycrys.getContactList(contacts, atm_list, _lower_bound[2 * cycle], _upper_bound[2 * cycle + 1], false);
                 if (contacts.empty()) continue;
                 for (int i = 0; i < 2; ++i) {
                      _process_contact_list(contacts, _lower_bound[2 * cycle + i], _upper_bound[2 * cycle + i], resIndexOrder_mapping, index_chain_id_mapping);
                      if (index_chain_id_mapping.size() == nonpolymer_chains.size()) break;
                 }
                 if (index_chain_id_mapping.size() == nonpolymer_chains.size()) break;
            }

            std::string new_pdb_chnid;
            std::set<std::string> pdb_chnids;
            std::vector<RCSB::Residue*> residues;
            int count = 0;
            for (std::vector<RCSB::Chain*>::iterator cpos = nonpolymer_chains.begin(); cpos != nonpolymer_chains.end(); ++cpos) {
                 std::map<int, std::string>::iterator pos = index_chain_id_mapping.find(count);
                 count++;

                 pdb_chnids.clear();
                 (*cpos)->GetFirstResidueList(residues);
                 while (!residues.empty()) {
                      for (std::vector<RCSB::Residue*>::iterator vpos = residues.begin(); vpos != residues.end(); ++vpos) {
                           std::string idx = CompositeIndex::getIndex((*vpos)->pdb_chnid(), (*vpos)->ResName(), (*vpos)->pdb_res_no(), (*vpos)->ins_code());
                           std::map<std::string, std::string>::const_iterator mpos1 = _nonpolymer_group_mapping.find(idx);
                           if (mpos1 != _nonpolymer_group_mapping.end()) pdb_chnids.insert(mpos1->second);
                      }
                      (*cpos)->GetNextResidueList(residues);
                 }

                 new_pdb_chnid.clear();
                 if (pdb_chnids.size() == 1) new_pdb_chnid = *(pdb_chnids.begin());
                 else if (pos != index_chain_id_mapping.end()) new_pdb_chnid = pos->second;
                 // else add warning info.

                 if (_multi_model_chain_flag) nonpolymer_chain_pairs.push_back(std::make_pair(new_pdb_chnid, *cpos));
                 else if (!new_pdb_chnid.empty()) (*cpos)->set_PDB_ChainID_to_residues(new_pdb_chnid);
            }
       }

       if (!nonpolymer_chain_pairs.empty()) _mol_nonpolymer_chains.push_back(nonpolymer_chain_pairs);

       mycrys.unmarkAtoms();
}

void ChainIDNumberAssignment::_get_index_mapping(std::map<int, int>& mapping, const std::vector<RCSB::Chain*>& chains)
// create residue index vs. chain vector index mapping
{
       mapping.clear();

       std::vector<RCSB::Residue*> residues;
       int count = 0;
       for (std::vector<RCSB::Chain*>::const_iterator cpos = chains.begin(); cpos != chains.end(); ++cpos) {
            (*cpos)->GetFirstResidueList(residues);
            while (!residues.empty()) {
                 for (std::vector<RCSB::Residue*>::iterator vpos = residues.begin(); vpos != residues.end(); ++vpos) {
                      mapping.insert(std::make_pair((*vpos)->index(), count));
                 }
                 (*cpos)->GetNextResidueList(residues);
            }
            count++;
       }
}

void ChainIDNumberAssignment::_process_contact_list(const std::list<CONTACT>& contacts, const double& lo, const double& hi, std::map<int, int>& index,
                                                    std::map<int, std::string>& chain_id_mapping)
{
       // first key: nonpolymer chain order index
       // second key: assigned PDB chain ID
       // multiset holds contact distances
       std::map<int, std::map<std::string, std::multiset<double> > > results;
       results.clear();

       std::multiset<double> tmp_set;
       std::map<std::string, std::multiset<double> > tmp_map;

       for (std::list<CONTACT>::const_iterator ptr = contacts.begin(); ptr != contacts.end(); ++ptr) {
            if (ptr->sym || ptr->lx || ptr->ly || ptr->lz) continue;
            if (ptr->dist < lo || ptr->dist > hi) continue;
            // find nonpolymer residue's chain order from residue's index
            std::map<int, int>::iterator pos = index.find(ptr->a_res_index);
            if (pos == index.end()) continue;

            // ignore if PDB chainID already assigned
            if (chain_id_mapping.find(pos->second) != chain_id_mapping.end())
                 continue;

            std::string assigned_chain_id = ptr->b_atm->pdb_chnid();

            std::map<int, std::map<std::string, std::multiset<double> > >::iterator mpos = results.find(pos->second);
            if (mpos == results.end()) {
                 tmp_set.clear();
                 tmp_set.insert(ptr->dist);
                 tmp_map.clear();
                 tmp_map.insert(std::make_pair(assigned_chain_id, tmp_set));
                 results.insert(std::make_pair(pos->second, tmp_map));
            } else {
                 std::map<std::string, std::multiset<double> >::iterator mmpos = mpos->second.find(assigned_chain_id);
                 if (mmpos == mpos->second.end()) {
                      tmp_set.clear();
                      tmp_set.insert(ptr->dist);
                      mpos->second.insert(std::make_pair(assigned_chain_id, tmp_set)); 
                 } else mmpos->second.insert(ptr->dist);
            }
       }

       if (results.empty()) return;

       for (std::map<int, std::map<std::string, std::multiset<double> > >::iterator mpos = results.begin(); mpos != results.end(); ++mpos) {
            if (mpos->second.size() == 1) {
                 std::map<std::string, std::multiset<double> >::iterator mmpos = mpos->second.begin();
                 chain_id_mapping.insert(std::make_pair(mpos->first, mmpos->first));
                 continue;
            }

            double min_dist = 3.0;
            unsigned int max_count = 0;
            std::string min_assigned_id = "";
            std::string max_assigned_id = "";
            for (std::map<std::string, std::multiset<double> >::iterator mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                 if (mmpos->second.size() > max_count) {
                      max_count = mmpos->second.size();
                      max_assigned_id = mmpos->first;
                 }
                 std::multiset<double>::iterator spos = mmpos->second.begin();
                 if ((*spos) < min_dist) {
                      min_dist = *spos;
                      min_assigned_id = mmpos->first;
                 } 
            }
            if (min_dist < 1.75 && !min_assigned_id.empty())
                 chain_id_mapping.insert(std::make_pair(mpos->first, min_assigned_id));
            else if (max_count > 0 && !max_assigned_id.empty())
                 chain_id_mapping.insert(std::make_pair(mpos->first, max_assigned_id));
       }
}
