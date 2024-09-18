/*
FILE:     getLinkSSBondUtil.C
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
#include <sys/stat.h>

#include "BondUtil.h"
#include "CategoryMapping.h"
#include "CompositeIndex.h"
#include "Element.h"
#include "SeqCodeUtil.h"
#include "TypeDef.h"
#include "UpdateAtomInfo.h"
#include "AnnotationObj.h"
#include "utillib.h"
#include "xtal.h"

#define NUM_GLYCO_LINK   13

const static char *__glyco_link_list[NUM_GLYCO_LINK][4] = {
       { "GHP",  "O4",   "BGC",   "C1"   },
       { "OMY",  "ODE",  "RER",   "C1"   },
       { "OMX",  "OC",   "ERE",   "C1"   },
       { "OMY",  "ODE",  "NAG",   "C1"   },
       { "3FG",  "OD1",  "BMA",   "C1"   },
       { "GHP",  "O4",   "GCS",   "C1"   },
       { "GHP",  "O4",   "N1L",   "C1"   },
       { "D4P",  "O4",   "MAN",   "C1"   },
       { "OMY",  "ODE",  "DVC",   "C1"   },
       { "02V",  "C24",  "RAM",   "O1"   },
       { "OMX",  "OC",   "NAG",   "C1"   },
       { "3FG",  "OD1",  "MAN",   "C1"   },
       { "ALA",  "N",    "MUB",   "C10"  }
};

static bool is_residue_in_list(RCSB::Residue* res, const std::vector<RCSB::Residue*>& res_list);

void AnnotationObj::Calculate_Link_and_SSBond(const int& type, const bool& force_flag)
{
       if (_molecules.empty()) return;

       if (_skip_task_set.find("link") != _skip_task_set.end() &&
           _skip_task_set.find("ssbond") != _skip_task_set.end()) return;

       if (!force_flag && (!_links.empty() || !_ssbonds.empty())) return;

       if ((type & FIND_LINK) && _skip_task_set.find("link") == _skip_task_set.end()) {
            _links.clear();
            _link_residue_set.clear();
            _covale_link_mapping.clear();
       }
       if ((type & FIND_SSBOND) && _skip_task_set.find("ssbond") == _skip_task_set.end()) _ssbonds.clear();

       std::set<std::string> short_protein_chain_id_set;
       short_protein_chain_id_set.clear();

       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            if ((chain->chain_type() ==  "ATOMP") && (chain->SeqLen() < 20)) short_protein_chain_id_set.insert(chain->PDB_ChainID());
            chain = _molecules[0]->GetNextChain();
       }

       if (!short_protein_chain_id_set.empty()) {
            for (int i = 0; i < NUM_GLYCO_LINK; ++i) {
                 _glyco_link_set.insert(CompositeIndex::getIndex(__glyco_link_list[i][0], __glyco_link_list[i][1],
                                                                 __glyco_link_list[i][2], __glyco_link_list[i][3]));
                 _glyco_link_set.insert(CompositeIndex::getIndex(__glyco_link_list[i][2], __glyco_link_list[i][3],
                                                                 __glyco_link_list[i][0], __glyco_link_list[i][1]));
            }
       }

       std::set<std::string> chain_types;
       chain_types.clear();
       std::set<std::string> residue_set;
       residue_set.clear();
       xtal mycrys;
       if (!set_mycrys(mycrys, _cell, _molecules[0], chain_types, residue_set, true, false, false)) return;

       std::vector<CONTACT> contact;
       get_contact(mycrys, contact, 0.85, 4.0, "asym", "all_asym", "all_asym", 1, 0, false);
       _find_Link_and_SSBond(type, contact, short_protein_chain_id_set);
       
       if (!_cell.is_artifical()) {
            get_contact(mycrys, contact, 0.85, 4.0, "sym", "all_asym", "all", 1, 0, false);
            _find_Link_and_SSBond(type, contact, short_protein_chain_id_set);
       }
}

bool AnnotationObj::Converting_Close_Contact_to_Link(const std::string& datafile)
{
       if (datafile.empty()) {
            _logIo->message("No input selected record file.\n");
            return false;
       }

       struct stat statbuf;
       if (stat(datafile.c_str(), &statbuf) != 0) {
            std::string error = "Input data file " + datafile + " does not exist.\n";
            _logIo->message(error.c_str());
            return false;
       }

       std::string cs;
       std::vector<std::string> data;

       int count = 0;
       _LINK link;
       link.SymOP_1.clear();
       link.SymOP_2.clear();
       link.details.clear();
       link.bondtype.clear();

       std::set<std::string> sugar_residue_set, sugar_residue_link_set;
       sugar_residue_set.clear();
       sugar_residue_link_set.clear();

       FILE *fp = fopen(datafile.c_str(), "r");
       while (!feof(fp)) {
            get_one_line(fp, cs);
            if (feof(fp)) break;
            if (cs.empty()) continue;

            get_wordarray_delimit_by_string(data, cs, "_");
            if (data.size() != 13) continue;

            RCSB::Residue* fstRes = _molecules[0]->find_pdb_residue(data[0], data[1], data[2], data[3]);
            if (!fstRes) {
                 std::string error = "Atom (" + data[0] + " " + data[1] + " " + data[2] + data[3] + " " + data[4] + data[5] 
                                   + ") does not eixt in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                 _logIo->message(error.c_str());
                 continue;
            }
            RCSB::Atom* fstAtom = fstRes->find_atom(data[4], data[5]);
            if (!fstAtom) {
                 std::string error = "Atom (" + data[0] + " " + data[1] + " " + data[2] + data[3] + " " + data[4] + data[5] 
                                   + ") does not eixt in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                 _logIo->message(error.c_str());
                 continue;
            }
            RCSB::Residue* sndRes = _molecules[0]->find_pdb_residue(data[6], data[7], data[8], data[9]);
            if (!sndRes) {
                 std::string error = "Atom (" + data[6] + " " + data[7] + " " + data[8] + data[9] + " " + data[10] + data[11] 
                                   + ") does not eixt in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                 _logIo->message(error.c_str());
                 continue;
            }
            RCSB::Atom* sndAtom = sndRes->find_atom(data[10], data[11]);
            if (!sndAtom) {
                 std::string error = "Atom (" + data[6] + " " + data[7] + " " + data[8] + data[9] + " " + data[10] + data[11] 
                                   + ") does not eixt in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                 _logIo->message(error.c_str());
                 continue;
            }

            try {
                 const ConnectFormat& drug_1 = _ccDic->find_drug(fstRes->ResName());
                 const ConnectFormat& drug_2 = _ccDic->find_drug(sndRes->ResName());
                 if ((drug_1.getMetaData("pdbx_type") == "ATOMS") && (drug_2.getMetaData("pdbx_type") == "ATOMS")) {
                      std::string idx1 = CompositeIndex::getIndex(data[0], data[1], data[2], data[3]);
                      sugar_residue_set.insert(idx1);
                      std::string idx2 = CompositeIndex::getIndex(data[6], data[7], data[8], data[9]);
                      sugar_residue_set.insert(idx2);
                      sugar_residue_link_set.insert(idx1 + "_" + idx2);
                      sugar_residue_link_set.insert(idx2 + "_" + idx1);
                 }
            } catch (const std::exception& exc) {}


            link.fstAtom = fstAtom;
            link.sndAtom = sndAtom;
            link.dist = FloatToString(cal_distance(fstAtom, sndAtom), 0, 3); 

            link.type.clear();
            link.leaving_flag.clear();
            int bond_type = BondUtil::get_bond_type(fstAtom->atom_type(), sndAtom->atom_type());
            if (bond_type == 2)
                 link.type = "metalc";
            else link.type = "covale";
            if (link.type == "covale") link.leaving_flag = _get_leaving_flag(fstAtom, sndAtom);
            _insert_a_link(link);

            count++;
       }
       fclose (fp);

       if (count == 0) return false;

       if (!sugar_residue_set.empty() && !sugar_residue_link_set.empty()) {
            for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
                 (*mpos)->merged_linked_sugar_chains(sugar_residue_set, sugar_residue_link_set);
            }
       }

       if (_CifObj) {
            Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
            std::set<std::string> allow_set;
            allow_set.clear();
            allow_set.insert("link");
            _cif_update_add_value_to_pdbx_data_processing_status(_cifblock, allow_set);
       }

       return true;
}

bool AnnotationObj::Removing_Covalent_Link(const std::string& datafile)
{
       if (datafile.empty()) {
            _logIo->message("No input selected record file.\n");
            return false;
       }

       struct stat statbuf;
       if (stat(datafile.c_str(), &statbuf) != 0) {
            std::string error = "Input data file " + datafile + " does not exist.\n";
            _logIo->message(error.c_str());
            return false;
       }

       std::string cs;
       std::vector<std::string> data;

       int count = 0;
       _LINK link;
       link.mol_index = 0;
       link.SymOP_1.clear();
       link.SymOP_2.clear();
       link.details.clear();
       link.bondtype.clear();
       link.type.clear();
       link.leaving_flag.clear();

       FILE *fp = fopen(datafile.c_str(), "r");
       while (!feof(fp)) {
            get_one_line(fp, cs);
            if (feof(fp)) break;
            if (cs.empty()) continue;

            get_wordarray_delimit_by_string(data, cs, "_");
            if (data.size() != 16) continue;

            RCSB::Residue* fstRes = _molecules[0]->find_pdb_residue(data[0], data[1], data[2], data[3]);
            if (!fstRes) {
                 std::string error = "Atom (" + data[0] + " " + data[1] + " " + data[2] + data[3] + " " + data[4] + data[5] 
                                   + ") does not eixt in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                 _logIo->message(error.c_str());
                 continue;
            }
            RCSB::Atom* fstAtom = fstRes->find_atom(data[4], data[5]);
            if (!fstAtom) {
                 std::string error = "Atom (" + data[0] + " " + data[1] + " " + data[2] + data[3] + " " + data[4] + data[5] 
                                   + ") does not eixt in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                 _logIo->message(error.c_str());
                 continue;
            }
            RCSB::Residue* sndRes = _molecules[0]->find_pdb_residue(data[7], data[8], data[9], data[10]);
            if (!sndRes) {
                 std::string error = "Atom (" + data[7] + " " + data[8] + " " + data[9] + data[10] + " " + data[11] + data[12] 
                                   + ") does not eixt in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                 _logIo->message(error.c_str());
                 continue;
            }
            RCSB::Atom* sndAtom = sndRes->find_atom(data[11], data[12]);
            if (!sndAtom) {
                 std::string error = "Atom (" + data[7] + " " + data[8] + " " + data[9] + data[10] + " " + data[11] + data[12] 
                                   + ") does not eixt in Model " + String::IntToString(_molecules[0]->Mol_ID()) + ".\n";
                 _logIo->message(error.c_str());
                 continue;
            }

            link.fstAtom = fstAtom;

            replace_string(data[6], "-", "_");
            if (!data[6].empty() && (data[6] != "1_555"))
                 link.SymOP_1 = data[6];
            else link.SymOP_1.clear();

            link.sndAtom = sndAtom;

            replace_string(data[13], "-", "_");
            if (!data[13].empty() && (data[13] != "1_555"))
                 link.SymOP_2 = data[13];
            else link.SymOP_2.clear();

            _remove_a_link(link);

            count++;
       }
       fclose (fp);

       if (count == 0) return false;

       if (_CifObj) {
            Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
            std::set<std::string> allow_set;
            allow_set.clear();
            allow_set.insert("link");
            _cif_update_add_value_to_pdbx_data_processing_status(_cifblock, allow_set);
       }

       return true;
}

void AnnotationObj::_find_Link_and_SSBond(const int& type, const std::vector<CONTACT>& contact, const std::set<std::string>& short_protein_chain_id_set)
{
       if (contact.empty()) return;

       _LINK link;
       //link.mol_index = _molecules[0]->index();
       link.mol_index = 0;
       link.SymOP_1.clear();
       link.details.clear();
       link.bondtype.clear();
       for (std::vector<CONTACT>::const_iterator pos = contact.begin(); pos != contact.end(); ++pos) {
            int bond_type = BondUtil::is_a_link(pos->a_atm->atom_type(), pos->b_atm->atom_type(), pos->dist);
            // loose distance restriction for S-S bond
            if (!bond_type && pos->a_atm->atmtype() == "SG" && pos->b_atm->atmtype() == "SG") {
                 if (pos->dist >= 1.8 && pos->dist <= 3.0) bond_type = 1;
            }
            if ((bond_type != 1) && (bond_type != 2))  continue;

            int assign_type = _assign_contact_type(bond_type, *pos, short_protein_chain_id_set, link.pdbx_role);
            if (!assign_type) continue;

            link.fstAtom = pos->a_atm;
            link.sndAtom = pos->b_atm;
            link.SymOP_2 = String::IntToString(pos->sym) + "_" + String::IntToString(pos->lx + 5)
                         + String::IntToString(pos->ly + 5) + String::IntToString(pos->lz + 5);
            if (link.SymOP_2 == "1_555") link.SymOP_2.clear();
            link.type.clear();
            link.dist = FloatToString(pos->dist, 0, 3);
            link.leaving_flag.clear();

            if ((type & FIND_LINK) && (assign_type & FIND_LINK)) {
                 if (bond_type == 2)
                      link.type = "metalc";
                 else link.type = "covale";
                 if (link.type == "covale") link.leaving_flag = _get_leaving_flag(link.fstAtom, link.sndAtom);
                 // skip non-sense links for standard residues (DAOTHER-3747)
                 if ((link.leaving_flag != "both") && SeqCodeUtil::is_a_standard_residue(link.fstAtom->pdb_resnam()) &&
                     SeqCodeUtil::is_a_standard_residue(link.sndAtom->pdb_resnam())) continue;
                 if (_skip_task_set.find("link") == _skip_task_set.end()) _insert_a_link(link);
            } else if ((type & FIND_SSBOND) && (assign_type & FIND_SSBOND)) {
                 link.type = "disulf";
                 // if (link.fstAtom->atmtype() == "SG" && link.sndAtom->atmtype() == "SG")
                 if (_skip_task_set.find("ssbond") == _skip_task_set.end()) _insert_a_ssbond(link);
            }
       }
}

int AnnotationObj::_assign_contact_type(const int& bond_type, const CONTACT& contact, const std::set<std::string>& short_protein_chain_id_set,
                                        std::string& glycosylation_type)
{
/*
       const std::set<std::string>& a_connected_atoms = _ccDic->get_connected_atoms(contact.a_atm->pdb_resnam());
       const std::set<std::string>& b_connected_atoms = _ccDic->get_connected_atoms(contact.b_atm->pdb_resnam());
*/

       glycosylation_type.clear();
       std::string res_type1 = _ccDic->find_residue_type(contact.a_atm->pdb_resnam());
       std::string res_type2 = _ccDic->find_residue_type(contact.b_atm->pdb_resnam());

       std::vector<RCSB::Residue*> firstResidueList, lastResidueList;
       int type = 0;
       if (contact.sym > 1 && contact.a_atm->pdb_resnam() == contact.b_atm->pdb_resnam() && contact.a_atm->pdb_chnid() == contact.b_atm->pdb_chnid() &&
           contact.a_atm->pdb_resnum() == contact.b_atm->pdb_resnum() && contact.a_atm->ins_code() == contact.b_atm->ins_code()) {
           if (res_type1 == "ATOMP" && res_type2 == "ATOMP" && contact.a_atm->atom_type() == "S" && contact.b_atm->atom_type() == "S") type = FIND_SSBOND;
       } else if (bond_type == 2)
            type = FIND_LINK;
       else if (res_type1 == "ATOMP" && res_type2 == "ATOMP" && contact.a_atm->atom_type() == "S" && contact.b_atm->atom_type() == "S")
            type = FIND_SSBOND;
       else if (contact.a_atm->atom_type() == "O" && contact.b_atm->atom_type() == "O") {
       } else if ((contact.a_atm->pdb_resnam() == "HOH") || (contact.a_atm->pdb_resnam() == "DOD") ||
                (contact.b_atm->pdb_resnam() == "HOH") || (contact.b_atm->pdb_resnam() == "DOD")) {
       } else if ((SeqCodeUtil::is_standard_aa_residue(contact.a_atm->pdb_resnam()) || contact.a_atm->pdb_resnam() == "UNK") &&
                  (SeqCodeUtil::is_standard_aa_residue(contact.b_atm->pdb_resnam()) || contact.b_atm->pdb_resnam() == "UNK")) {
            if ((contact.a_atm->atmtype() != "C" && contact.a_atm->atmtype() != "N") || (contact.b_atm->atmtype() != "C" && contact.b_atm->atmtype() != "N") ||
                (contact.a_atm->pdb_chnid() != contact.b_atm->pdb_chnid()) || (contact.a_atm->chnid() != contact.b_atm->chnid())) {
/*
                 if (a_connected_atoms.find(contact.a_atm->atmtype()) != a_connected_atoms.end() &&
                     b_connected_atoms.find(contact.b_atm->atmtype()) != b_connected_atoms.end())
*/
                      type = FIND_LINK;
            } else {
                 RCSB::Chain* chain = _molecules[0]->GetPolyChain(contact.a_atm->pdb_chnid());
                 RCSB::Residue* a_res = _molecules[0]->find_pdb_residue(contact.a_atm->pdb_chnid(), contact.a_atm->pdb_resnam(),
                                                                        contact.a_atm->pdb_resnum(), contact.a_atm->ins_code());
                 RCSB::Residue* b_res = _molecules[0]->find_pdb_residue(contact.b_atm->pdb_chnid(), contact.b_atm->pdb_resnam(),
                                                                        contact.b_atm->pdb_resnum(), contact.b_atm->ins_code());
                 if (chain && a_res && b_res && (chain->ResidueNumbers() > 2)) {
                      chain->GetFirstResidueList(firstResidueList);
                      chain->GetLastResidueList(lastResidueList);
                      if ((is_residue_in_list(a_res, firstResidueList) && is_residue_in_list(b_res, lastResidueList)) ||
                          (is_residue_in_list(a_res, lastResidueList) && is_residue_in_list(b_res, firstResidueList))) type = FIND_LINK;
                 }
            }
       } else if (SeqCodeUtil::is_na_residue(contact.a_atm->pdb_resnam()) && SeqCodeUtil::is_na_residue(contact.b_atm->pdb_resnam())) {
            if ((contact.a_atm->atmtype() != "O3'" && contact.a_atm->atmtype() != "P") ||
                (contact.b_atm->atmtype() != "O3'" && contact.b_atm->atmtype() != "P")) {
/*
                 if (a_connected_atoms.find(contact.a_atm->atmtype()) != a_connected_atoms.end() &&
                     b_connected_atoms.find(contact.b_atm->atmtype()) != b_connected_atoms.end())
*/
                      type = FIND_LINK;
            }
       } else if (((res_type1 == "ATOMP") && (res_type2 == "ATOMS")) || ((res_type1 == "ATOMS") && (res_type2 == "ATOMP"))) {
            glycosylation_type = _find_glycosylation_type(contact.a_atm, contact.b_atm);
            if (!glycosylation_type.empty()) {
                 type = FIND_LINK;
                 RCSB::Chain* chain1 = _molecules[0]->GetAsymChain(contact.a_atm->chnid());
                 RCSB::Chain* chain2 = _molecules[0]->GetAsymChain(contact.b_atm->chnid());
                 // remove Glycosylation type if the amino acid residue is not in a polymer. ( DAOTHER-9215 )
                 if (chain1 && (chain1->chain_type() != "ATOMP") && chain2 && (chain2->chain_type() != "ATOMP")) glycosylation_type.clear();
            } else if ((((res_type1 == "ATOMP") && (short_protein_chain_id_set.find(contact.a_atm->pdb_chnid()) != short_protein_chain_id_set.end())) ||
                      ((res_type2 == "ATOMP") && (short_protein_chain_id_set.find(contact.b_atm->pdb_chnid()) != short_protein_chain_id_set.end()))) &&
                        _find_glyco_like_linkage(contact.a_atm, contact.b_atm)) type = FIND_LINK;
#if 0
       } else if (res_type1 == "ATOMP" && res_type2 == "ATOMS") {
            if (contact.b_atm->atom_type() == "C" &&
               (contact.a_atm->pdb_resnam() == "ASN" && // contact.b_atm->pdb_resnam() == "NAG" &&
                contact.a_atm->atmtype() == "ND2" ||
                contact.a_atm->pdb_resnam() == "SER" && // contact.b_atm->pdb_resnam() == "NGA" &&
                contact.a_atm->atmtype() == "OG" ||
                contact.a_atm->pdb_resnam() == "THR" && // contact.b_atm->pdb_resnam() == "NGA" &&
                contact.a_atm->atmtype() == "OG1")) type = FIND_LINK;
       } else if (res_type1 == "ATOMS" && res_type2 == "ATOMP") {
            if (contact.a_atm->atom_type() == "C" &&
               (/* contact.a_atm->pdb_resnam() == "NAG" && */ contact.b_atm->pdb_resnam() == "ASN" &&
                contact.b_atm->atmtype() == "ND2" ||
                /* contact.a_atm->pdb_resnam() == "NGA" && */ contact.b_atm->pdb_resnam() == "SER" &&
                contact.b_atm->atmtype() == "OG" ||
                /* contact.a_atm->pdb_resnam() == "NGA" && */ contact.b_atm->pdb_resnam() == "THR" &&
                contact.b_atm->atmtype() == "OG1")) type = FIND_LINK;
#endif
       } else if (res_type1 == "ATOMS" && res_type2 == "ATOMS") {
            if ((contact.a_atm->atom_type() == "C" && contact.b_atm->atom_type() == "O") ||
                (contact.a_atm->atom_type() == "O" && contact.b_atm->atom_type() == "C"))
                 type = FIND_LINK;
       } else /* if (a_connected_atoms.find(contact.a_atm->atmtype()) != a_connected_atoms.end() &&
                  b_connected_atoms.find(contact.b_atm->atmtype()) != b_connected_atoms.end()) */
            type = FIND_LINK;

       return type;
}

bool AnnotationObj::_find_glyco_like_linkage(RCSB::Atom* atom1, RCSB::Atom* atom2)
{
       std::string idx = CompositeIndex::getIndex(atom1->pdb_resnam(), atom1->atmtype(), atom2->pdb_resnam(), atom2->atmtype());
       if (_glyco_link_set.find(idx) != _glyco_link_set.end()) return true;

       return false;
}

static bool is_residue_in_list(RCSB::Residue* res, const std::vector<RCSB::Residue*>& res_list)
{
       if (!res || res_list.empty()) return false;

       for (unsigned int i = 0; i < res_list.size(); ++i) {
            if (res == res_list[i]) return true;
       }
       return false;
}

#if 0
std::string AnnotationObj::_get_leaving_flag(const Atom* atom1, const Atom* atom2)
{
       const std::set<std::string>& a_connected_atoms = _ccDic->get_connected_atoms(atom1->pdb_resnam());
       const std::set<std::string>& b_connected_atoms = _ccDic->get_connected_atoms(atom2->pdb_resnam());

       int count = 0;
       if (a_connected_atoms.find(atom1->atmtype()) != a_connected_atoms.end()) count++;
       if (b_connected_atoms.find(atom2->atmtype()) != b_connected_atoms.end()) count++;

       if (count == 2)
            return "both";
       else if (count == 1)
            return "one";
       else return "none";
}

void AnnotationObj::_update_leaving_atom_flag()
{
       if (_links.empty()) return;

       for (std::list<_LINK>::iterator pos = _links.begin(); pos != _links.end(); ++pos) {
            if ((pos->type != "covale") || !pos->leaving_flag.empty()) continue;
            pos->leaving_flag = _get_leaving_flag(pos->fstAtom, pos->sndAtom);
       }
}
#endif
void AnnotationObj::Update_StructLink_and_StructConn()
{
       if (!_CifObj) return;

       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
       _update_pdbx_struct_link(_cifblock);
       _update_struct_conn(_cifblock);
}

void AnnotationObj::_update_pdbx_struct_link(Block& block)
{
       // if (_links.empty()) {
            deleteTable(block, "pdbx_struct_link");
            return;
       // }

       try {
            ISTable *t = _newTablePtr("pdbx_struct_link");
            std::multimap<double, _LINK> sort_link;
            sort_link.clear();
            for (std::list<_LINK>::const_iterator pos = _links.begin(); pos != _links.end(); ++pos) {
                 sort_link.insert(std::make_pair(atof(pos->dist.c_str()), *pos));
            }
            int i = 0;
            for (std::multimap<double, _LINK>::const_iterator pos = sort_link.begin(); pos != sort_link.end(); ++pos) {
                 t->AddRow();
                 std::string sym1 = pos->second.SymOP_1; if (sym1.empty()) sym1 = "1_555";
                 std::string sym2 = pos->second.SymOP_2; if (sym2.empty()) sym2 = "1_555";
                 t->UpdateCell(i, "id", String::IntToString(i + 1)); 
                 t->UpdateCell(i, "ptnr1_symmetry", sym1);
                 t->UpdateCell(i, "ptnr2_symmetry", sym2);
                 t->UpdateCell(i, "details", pos->second.details);
                 t->UpdateCell(i, "type", pos->second.type);
                 t->UpdateCell(i, "pdbx_dist_value", pos->second.dist);
                 UpdateAtomInfo::UpdateTable_2(t, i, pos->second.fstAtom, "ptnr1_", NUM_ATOM_ITEM);
                 UpdateAtomInfo::UpdateTable_2(t, i, pos->second.sndAtom, "ptnr2_", NUM_ATOM_ITEM);
                 i++;
            }
            block.WriteTable(t);
       } catch (const std::exception& exc) {
            _logIo->messageError(exc.what());
       }
}
#if 0
void AnnotationObj::_update_struct_conn(Block& block)
{
       if (_links.empty() && _ssbonds.empty() && _sltbrgs.empty() && _bspairs.empty()) {
            deleteTable(block, "struct_conn");
            deleteTable(block, "struct_conn_type");
            return;
       }

       try {
            _update_leaving_atom_flag();

            std::vector<std::string> type_array;
            std::set<std::string> type_set;
            type_array.clear();
            type_set.clear();

            ISTable *t = _newTablePtr("struct_conn");
            _update_struct_conn_table(t, type_array, type_set, _ssbonds);
            _update_struct_conn_table(t, type_array, type_set, _links);
            _update_struct_conn_table(t, type_array, type_set, _bspairs);
            _update_struct_conn_table(t, type_array, type_set, _sltbrgs);
            block.WriteTable(t);

            if (!type_array.empty()) {
                 t = _newTablePtr("struct_conn_type");
                 for (unsigned int i = 0; i < type_array.size(); ++i) {
                      t->AddRow();
                      t->UpdateCell(i, "id", type_array[i]);
                 }
                 block.WriteTable(t);
            }
       } catch (const std::exception& exc) {
            _logIo->messageError(exc.what());
       }
}

void AnnotationObj::_update_struct_conn_table(ISTable* t, std::vector<std::string>& type_array, std::set<std::string>& type_set, const std::list<_LINK>& links)
{
       if (links.empty()) return;
/*
       std::multimap<double, _LINK> sort_link;
       sort_link.clear();
       for (std::list<_LINK>::const_iterator
            pos = links.begin(); pos != links.end(); ++pos) {
            sort_link.insert(std::make_pair(atof(pos->dist.c_str()), *pos));
       }
*/
       std::map<std::string, int> id_mapping;
       id_mapping.clear();

       int i = t->GetNumRows();
/*
       for (std::multimap<double, _LINK>::const_iterator pos = sort_link.begin(); pos != sort_link.end(); ++pos) {
            t->AddRow();
            std::string sym1 = pos->second.SymOP_1; if (sym1.empty()) sym1 = "1_555";
            std::string sym2 = pos->second.SymOP_2; if (sym2.empty()) sym2 = "1_555";
            std::map<std::string, int>::iterator mpos = id_mapping.find(pos->second.type);
            if (mpos == id_mapping.end()) {
                id_mapping.insert(std::make_pair(pos->second.type, 1));
                mpos = id_mapping.find(pos->second.type);
            } else mpos->second++;
            t->UpdateCell(i, "id", pos->second.type + String::IntToString(mpos->second)); 
            if (type_set.find(pos->second.type) == type_set.end()) {
                 type_array.push_back(pos->second.type);
                 type_set.insert(pos->second.type);
            }
            t->UpdateCell(i, "conn_type_id", pos->second.type);
            t->UpdateCell(i, "ptnr1_symmetry", sym1);
            t->UpdateCell(i, "ptnr2_symmetry", sym2);
            t->UpdateCell(i, "details", pos->second.details);
            t->UpdateCell(i, "pdbx_dist_value", pos->second.dist);
            t->UpdateCell(i, "pdbx_value_order", pos->second.bondtype);
            t->UpdateCell(i, "pdbx_leaving_atom_flag", pos->second.leaving_flag);
            UpdateAtomInfo::UpdateTable_3(t, i, pos->second.fstAtom, "ptnr1_", NUM_ALL_ITEM);
            UpdateAtomInfo::UpdateTable_3(t, i, pos->second.sndAtom, "ptnr2_", NUM_ALL_ITEM);
            i++;
       }
*/
       for (std::list<_LINK>::const_iterator pos = links.begin(); pos != links.end(); ++pos) {
            t->AddRow();
            std::string sym1 = pos->SymOP_1; if (sym1.empty()) sym1 = "1_555";
            std::string sym2 = pos->SymOP_2; if (sym2.empty()) sym2 = "1_555";
            std::map<std::string, int>::iterator mpos = id_mapping.find(pos->type);
            if (mpos == id_mapping.end()) {
                id_mapping.insert(std::make_pair(pos->type, 1));
                mpos = id_mapping.find(pos->type);
            } else mpos->second++;
            t->UpdateCell(i, "id", pos->type + String::IntToString(mpos->second)); 
            if (type_set.find(pos->type) == type_set.end()) {
                 type_array.push_back(pos->type);
                 type_set.insert(pos->type);
            }
            t->UpdateCell(i, "conn_type_id", pos->type);
            t->UpdateCell(i, "ptnr1_symmetry", sym1);
            t->UpdateCell(i, "ptnr2_symmetry", sym2);
            t->UpdateCell(i, "details", pos->details);
            t->UpdateCell(i, "pdbx_dist_value", pos->dist);
            t->UpdateCell(i, "pdbx_value_order", pos->bondtype);
            t->UpdateCell(i, "pdbx_leaving_atom_flag", pos->leaving_flag);
            UpdateAtomInfo::UpdateTable_3(t, i, pos->fstAtom, "ptnr1_", NUM_ALL_ITEM);
            UpdateAtomInfo::UpdateTable_3(t, i, pos->sndAtom, "ptnr2_", NUM_ALL_ITEM);
            i++;
       }
}
#endif
