/*
FILE:     LinkSSBondUtil.C
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

#include "BondUtil.h"
#include "CategoryMapping.h"
#include "CompositeIndex.h"
#include "Element.h"
#include "GetPairList.h"
#include "GetPairList.C"
#include "SeqCodeUtil.h"
#include "TypeDef.h"
#include "UpdateAtomInfo.h"
#include "ValidateObj.h"
#include "utillib.h"
#include "xtal.h"

void ValidateObj::_insert_a_link(const _LINK& link)
{
       if ((link.SymOP_1.empty() || link.SymOP_1 == "1_555") && (link.SymOP_2.empty() || link.SymOP_2 == "1_555")) {
            _insert_link_mapping(link.fstAtom, link.sndAtom, link.type);
            _insert_link_mapping(link.sndAtom, link.fstAtom, link.type);
       }

       for (std::list<_LINK>::iterator pos = _links.begin(); pos != _links.end(); ++pos) {
            if (pos->mol_index == link.mol_index && (((pos->fstAtom == link.fstAtom) && (pos->sndAtom == link.sndAtom) && (pos->SymOP_1 == link.SymOP_1) &&
               (pos->SymOP_2 == link.SymOP_2)) || ((pos->fstAtom == link.sndAtom) && (pos->sndAtom == link.fstAtom) && (pos->SymOP_1 == link.SymOP_2) &&
               (pos->SymOP_2 == link.SymOP_1)))) {
                 if (pos->details.empty()) pos->details = link.details;
                 if (pos->type.empty()) pos->type = link.type;
                 if (pos->dist.empty()) pos->dist = link.dist;
                 if (pos->bondtype.empty()) pos->bondtype = link.bondtype;
                 if (pos->leaving_flag.empty()) pos->leaving_flag = link.leaving_flag;
                 if (pos->pdbx_role.empty()) pos->pdbx_role = link.pdbx_role;
                 return;
            }
       }

       _links.push_back(link);
}

void ValidateObj::_remove_a_link(const _LINK& link)
{
       for (std::list<_LINK>::iterator pos = _links.begin(); pos != _links.end(); ++pos) {
            if (pos->mol_index == link.mol_index && (((pos->fstAtom == link.fstAtom) && (pos->sndAtom == link.sndAtom) && (pos->SymOP_1 == link.SymOP_1) &&
               (pos->SymOP_2 == link.SymOP_2)) || ((pos->fstAtom == link.sndAtom) && (pos->sndAtom == link.fstAtom) && (pos->SymOP_1 == link.SymOP_2) &&
               (pos->SymOP_2 == link.SymOP_1)))) {
                 _links.erase(pos);
                 return;
            }
       }
}

void ValidateObj::_insert_link_mapping(RCSB::Atom* fstAtom, RCSB::Atom* sndAtom, const std::string& type)
{
       std::string idx = fstAtom->pdb_chnid() + "_" + fstAtom->pdb_resnam() + "_" + fstAtom->pdb_resnum() + "_" + fstAtom->ins_code() + "_"
                       + sndAtom->pdb_chnid() + "_" + sndAtom->pdb_resnam() + "_" + sndAtom->pdb_resnum() + "_" + sndAtom->ins_code();
       if (!_is_glycosylation_site(fstAtom, sndAtom)) _link_residue_set.insert(idx);

       if (type != "covale") return;

       std::map<std::string, std::vector<std::pair<std::string, std::string> > >::iterator mpos = _covale_link_mapping.find(idx);
       if (mpos == _covale_link_mapping.end()) {
            std::vector<std::pair<std::string, std::string> > p_vec;
            p_vec.clear();
            p_vec.push_back(std::make_pair(fstAtom->pdb_atmnam(), sndAtom->pdb_atmnam()));
            _covale_link_mapping.insert(std::make_pair(idx, p_vec));
            return;
       }

       for (std::vector<std::pair<std::string, std::string> >::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
            if (vpos->first == fstAtom->pdb_atmnam() && vpos->second == sndAtom->pdb_atmnam()) return;
       }

       mpos->second.push_back(std::make_pair(fstAtom->pdb_atmnam(), sndAtom->pdb_atmnam()));
}

void ValidateObj::_insert_a_ssbond(const _SSBOND& ssbond)
{
       for (std::list<_SSBOND>::iterator
            pos = _ssbonds.begin(); pos != _ssbonds.end(); ++pos) {
            if (pos->mol_index == ssbond.mol_index &&
               ((pos->fstAtom == ssbond.fstAtom && pos->sndAtom == ssbond.sndAtom && pos->SymOP_1 == ssbond.SymOP_1 && pos->SymOP_2 == ssbond.SymOP_2) ||
                (pos->fstAtom == ssbond.sndAtom && pos->sndAtom == ssbond.fstAtom && pos->SymOP_1 == ssbond.SymOP_2 && pos->SymOP_2 == ssbond.SymOP_1))) {
                 if (pos->details.empty()) pos->details = ssbond.details;
                 if (pos->type.empty()) pos->type = ssbond.type;
                 if (pos->dist.empty()) pos->dist = ssbond.dist;
                 if (pos->bondtype.empty()) pos->bondtype = ssbond.bondtype;
                 if (pos->leaving_flag.empty()) pos->leaving_flag = ssbond.leaving_flag;
                 if (pos->pdbx_role.empty()) pos->pdbx_role = ssbond.pdbx_role;
                 return;
            }
       }

       _ssbonds.push_back(ssbond);
}

void ValidateObj::InputPDBLINKRecord()
{
       if (_molecules.empty()) return;

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = _pdb_records.find("LINK");
       if (mpos == _pdb_records.end()) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > names;
       std::vector<RCSB::Atom*> atoms;
       _LINK link;

       for (std::list<std::vector<std::string> >::const_iterator pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
            names.clear();
            data.clear();
            data.push_back((*pos)[4]);
            data.push_back((*pos)[3]);
            data.push_back((*pos)[5]);
            data.push_back((*pos)[6]);
            data.push_back((*pos)[1]);
            data.push_back((*pos)[2]);
            names.push_back(data);
            data.clear();
            data.push_back((*pos)[10]);
            data.push_back((*pos)[9]);
            data.push_back((*pos)[11]);
            data.push_back((*pos)[12]);
            data.push_back((*pos)[7]);
            data.push_back((*pos)[8]);
            names.push_back(data);

            for (std::vector<RCSB::Molecule*>::iterator MolPos = _molecules.begin(); MolPos != _molecules.end(); ++MolPos) {
                 _getAtomsfromNames((*MolPos)->Mol_ID(), names, link.mol_index, atoms, data, true);
                 if (!atoms.empty()) break;
            }
            if (atoms.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 continue;
            }

            link.fstAtom = atoms[0];
            link.sndAtom = atoms[1];
            link.SymOP_1 = _reformat_symmetry((*pos)[13]);
            link.SymOP_2 = _reformat_symmetry((*pos)[14]);
            link.details.clear();
            link.type = "covale";
            if (Element::findMetalFlag(atoms[0]->atom_type()) || Element::findMetalFlag(atoms[1]->atom_type()))
                 link.type = "metalc";
            double dist = _cell.cal_symmetry_distance(atoms[0], atoms[1], link.SymOP_1, link.SymOP_2);
            link.dist = FloatToString(dist, 0, 3);
            link.bondtype.clear();
            link.leaving_flag.clear();
            link.pdbx_role.clear();

            _insert_a_link(link);
       }
}

void ValidateObj::InputPDBSSBONDRecord()
{
       if (_molecules.empty()) return;

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = _pdb_records.find("SSBOND");
       if (mpos == _pdb_records.end()) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > names;
       std::vector<RCSB::Residue*> residues;
       std::vector<RCSB::Atom*> atoms;
       std::vector<std::vector<RCSB::Atom*> > atom_lists, pair_lists;
       _SSBOND ssbond;
       ssbond.type = "disulf";
       ssbond.bondtype.clear();
       ssbond.leaving_flag.clear();
       ssbond.pdbx_role.clear();

       for (std::list<std::vector<std::string> >::const_iterator pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
            names.clear();
            data.clear();
            data.push_back((*pos)[3]);
            data.push_back((*pos)[2]);
            data.push_back((*pos)[4]);
            data.push_back((*pos)[5]);
            names.push_back(data);
            data.clear();
            data.push_back((*pos)[7]);
            data.push_back((*pos)[6]);
            data.push_back((*pos)[8]);
            data.push_back((*pos)[9]);
            names.push_back(data);

            _getResiduesfromNames(0, names, ssbond.mol_index, residues, data);
            if (residues.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 continue;
            }

            atom_lists.clear();
            for (unsigned int i = 0; i < residues.size(); ++i) {
                 residues[i]->find_atom("SG", atoms);
                 if (atoms.empty()) {
                      std::string error = "Atom ('SG' '" + names[i][0] + "' '" + names[i][1] + "' '" + names[i][2] + names[i][3] + "') can not be found.";
                      _logIo->message("%s\n", error.c_str());
                      continue;
                 }
                 atom_lists.push_back(atoms);
            }
            if (atom_lists.size() != residues.size()) continue;

            GetPairList(atom_lists, pair_lists);
            if (pair_lists.empty()) continue;
  
            ssbond.SymOP_1 = _reformat_symmetry((*pos)[10]);
            ssbond.SymOP_2 = _reformat_symmetry((*pos)[11]);

            for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = pair_lists.begin(); apos != pair_lists.end(); ++apos) {
                 double dist = _cell.cal_symmetry_distance((*apos)[0], (*apos)[1], ssbond.SymOP_1, ssbond.SymOP_2);
                 if (dist > 3.0) continue;
                 ssbond.fstAtom = (*apos)[0];
                 ssbond.sndAtom = (*apos)[1];
                 ssbond.dist = FloatToString(dist, 0, 3);
                 _insert_a_ssbond(ssbond);
            }
       }
}

void ValidateObj::_read_struct_conn_category(ISTable* t)
{
       if (_molecules.empty()) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();

       data.clear();
       data.push_back("ptnr1_auth_asym_id");
       data.push_back("ptnr1_auth_comp_id");
       data.push_back("ptnr1_auth_seq_id");
       data.push_back("pdbx_ptnr1_PDB_ins_code");
       data.push_back("ptnr1_label_atom_id");
       data.push_back("pdbx_ptnr1_label_alt_id");
       items.push_back(data);
       data.clear();
       data.push_back("ptnr2_auth_asym_id");
       data.push_back("ptnr2_auth_comp_id");
       data.push_back("ptnr2_auth_seq_id");
       data.push_back("pdbx_ptnr2_PDB_ins_code");
       data.push_back("ptnr2_label_atom_id");
       data.push_back("pdbx_ptnr2_label_alt_id");
       items.push_back(data);

       _LINK link;
       std::vector<RCSB::Atom*> atoms;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean_lower(link.type, t, i, "conn_type_id");
            if (link.type.empty()) continue;
            if (link.type != "link" && link.type != "covale" && link.type != "metalc" && link.type != "disulf" &&
             /* link.type != "sltbrg" && link.type != "saltbr" && */ link.type != "hydrog") continue;

            get_values_clean(t, i, items, names);
            for (std::vector<RCSB::Molecule*>::iterator MolPos = _molecules.begin(); MolPos != _molecules.end(); ++MolPos) {
                 _getAtomsfromNames((*MolPos)->Mol_ID(), names, link.mol_index, atoms, data, true);
                 if (!atoms.empty()) break;
            }
            if (atoms.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 continue;
            }

            link.fstAtom = atoms[0];
            link.sndAtom = atoms[1];
            get_value(link.SymOP_1, t, i, "ptnr1_symmetry");
            if (link.SymOP_1 == "1_555") link.SymOP_1.clear();
            get_value(link.SymOP_2, t, i, "ptnr2_symmetry");
            if (link.SymOP_2 == "1_555") link.SymOP_2.clear();
            get_value(link.details, t, i, "details");
            get_value(link.bondtype, t, i, "pdbx_value_order");
            get_value(link.leaving_flag, t, i, "pdbx_leaving_atom_flag");
            get_value(link.pdbx_role, t, i, "pdbx_role");

            link.dist.clear();
            get_value(link.dist, t, i, "pdbx_dist_value");
            if (link.type == "link" || link.type == "covale" || link.type == "metalc" || link.type == "disulf") {
                 if (link.type == "link" && link.fstAtom->atmtype() == "SG" && link.fstAtom->pdb_resnam() == "CYS" &&
                     link.sndAtom->atmtype() == "SG" && link.sndAtom->pdb_resnam() == "CYS") link.type = "disulf";
                 if (link.type != "disulf") {
                      link.type = "covale";
                      if (Element::findMetalFlag(atoms[0]->atom_type()) || Element::findMetalFlag(atoms[1]->atom_type()))
                           link.type = "metalc";
                 }
                 // if (link.dist.empty()) {
                      double d = _cell.cal_symmetry_distance(atoms[0], atoms[1], link.SymOP_1, link.SymOP_2);
                      link.dist = FloatToString(d, 0, 3);
                 // }
            }

            if (link.type == "covale" || link.type == "metalc") {
/*
                 if (link.type == "covale") {
                      try {
                           const ConnectFormat& drug_1 = _ccDic->find_drug(link.fstAtom->pdb_resnam());
                           const ConnectFormat& drug_2 = _ccDic->find_drug(link.sndAtom->pdb_resnam());
                           if ((drug_1.getMetaData("pdbx_type") == "ATOMS") && (drug_2.getMetaData("pdbx_type") == "ATOMS")) {
                                if (!BondUtil::is_a_bond_with_extension(link.fstAtom->atom_type(), link.sndAtom->atom_type(), atof(link.dist.c_str())))
                                     continue;
                           }
                      } catch (const std::exception& exc) {}
                 }
*/
                 if ((link.type == "covale") && link.pdbx_role.empty()) {
                      link.pdbx_role = _find_glycosylation_type(link.fstAtom, link.sndAtom);
                      if (!link.pdbx_role.empty()) {
                           RCSB::Chain* chain1 = _molecules[link.mol_index]->GetAsymChain(link.fstAtom->chnid());
                           RCSB::Chain* chain2 = _molecules[link.mol_index]->GetAsymChain(link.sndAtom->chnid());
                           // remove Glycosylation type if the amino acid residue is not in a polymer. ( DAOTHER-9215 )
                           if (chain1 && (chain1->chain_type() != "ATOMP") && chain2 && (chain2->chain_type() != "ATOMP")) link.pdbx_role.clear();
                      }
                 }
                 _insert_a_link(link);
            } else if (link.type == "disulf")
                 _insert_a_ssbond(link);
            else if (link.type == "sltbrg" || link.type == "saltbr") {
                 link.type = "sltbrg";
                 _insert_a_sltbrg(link);
            } else if (link.type == "hydrog")
                 _insert_a_bspair(link);
       }
}

void ValidateObj::_read_pdbx_struct_link_category(ISTable* t)
{
       if (_molecules.empty()) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > items, names;
       items.clear();

       data.clear();
       data.push_back("ptnr1_label_asym_id");
       data.push_back("ptnr1_label_comp_id");
       data.push_back("ptnr1_label_seq_id");
       data.push_back("ptnr1_label_ins_code");
       data.push_back("ptnr1_label_atom_id");
       data.push_back("ptnr1_label_alt_id");
       items.push_back(data);
       data.clear();
       data.push_back("ptnr2_label_asym_id");
       data.push_back("ptnr2_label_comp_id");
       data.push_back("ptnr2_label_seq_id");
       data.push_back("ptnr2_label_ins_code");
       data.push_back("ptnr2_label_atom_id");
       data.push_back("ptnr2_label_alt_id");
       items.push_back(data);

       _LINK link;
       std::vector<RCSB::Atom*> atoms;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_values(t, i, items, names);
            for (std::vector<RCSB::Molecule*>::iterator MolPos = _molecules.begin(); MolPos != _molecules.end(); ++MolPos) {
                 _getAtomsfromNames((*MolPos)->Mol_ID(), names, link.mol_index, atoms, data, true);
                 if (!atoms.empty()) break;
            }
            if (atoms.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 continue;
            }

            link.fstAtom = atoms[0];
            link.sndAtom = atoms[1];
            get_value(link.SymOP_1, t, i, "ptnr1_symmetry");
            if (link.SymOP_1 == "1_555") link.SymOP_1.clear();
            get_value(link.SymOP_2, t, i, "ptnr2_symmetry");
            if (link.SymOP_2 == "1_555") link.SymOP_2.clear();
            get_value(link.details, t, i, "details");
            get_value(link.bondtype, t, i, "pdbx_value_order");
            link.leaving_flag.clear();
            link.pdbx_role.clear();

            link.type = "covale";
            if (Element::findMetalFlag(atoms[0]->atom_type()) ||
                Element::findMetalFlag(atoms[1]->atom_type()))
                 link.type = "metalc";

            double d = _cell.cal_symmetry_distance(atoms[0], atoms[1], link.SymOP_1, link.SymOP_2);
            link.dist = FloatToString(d, 0, 3);

            _insert_a_link(link);
       }
}

void ValidateObj::_read_pdbx_struct_link_and_struct_conn(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_struct_link");
       if (t) _read_pdbx_struct_link_category(t);
       t = getTablePtr(block, "struct_conn");
       if (t) {
            if (_DepUI_Flag) {
                 vector<string> data;
                 std::vector<std::pair<std::string, std::string> > label_auth_item_pairs;
                 label_auth_item_pairs.clear();
                 label_auth_item_pairs.push_back(std::make_pair("ptnr1_label_comp_id", "ptnr1_auth_comp_id"));
                 label_auth_item_pairs.push_back(std::make_pair("ptnr1_label_asym_id", "ptnr1_auth_asym_id"));
                 label_auth_item_pairs.push_back(std::make_pair("ptnr2_label_comp_id", "ptnr2_auth_comp_id"));
                 label_auth_item_pairs.push_back(std::make_pair("ptnr2_label_asym_id", "ptnr2_auth_asym_id"));
                 for (std::vector<std::pair<std::string, std::string> >::const_iterator
                      pos = label_auth_item_pairs.begin(); pos != label_auth_item_pairs.end(); ++pos) {
                      if (t->IsColumnPresent(pos->first) && !t->IsColumnPresent(pos->second)) {
                           t->GetColumn(data, pos->first);
                           t->AddColumn(pos->second, data);
                      }
                 }

                 std::string type;
                 int rowNo = t->GetNumRows();
                 for (int i = 0; i < rowNo; ++i) {
                      get_value_clean_lower(type, t, i, "conn_type_id");
                      if (type == "covale" || type == "metalc" || type == "disulf" || type == "sltbrg" || type == "saltbr" || type == "hydrog")
                           t->UpdateCell(i, "conn_type_id", type);
                      else t->UpdateCell(i, "conn_type_id", "link");
                 }

                 block.WriteTable(t);
            }

            _read_struct_conn_category(t);
       }
}

void ValidateObj::Read_StructLink_and_StructConn()
{
       if (!_CifObj) return;

       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
       _read_pdbx_struct_link_and_struct_conn(_cifblock);
}

std::string ValidateObj::_get_leaving_flag(const RCSB::Atom* atom1, const RCSB::Atom* atom2)
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

void ValidateObj::_update_leaving_atom_flag()
{
       if (_links.empty()) return;

       for (std::list<_LINK>::iterator pos = _links.begin(); pos != _links.end(); ++pos) {
            if ((pos->type != "covale") /* || !pos->leaving_flag.empty() */) continue;
            pos->leaving_flag = _get_leaving_flag(pos->fstAtom, pos->sndAtom);
       }
}

void ValidateObj::_update_struct_conn(Block& block)
{
       if (_links.empty() && _ssbonds.empty() && _sltbrgs.empty() && _bspairs.empty()) {
            deleteTable(block, "struct_conn");
            deleteTable(block, "struct_conn_type");
            return;
       }

       try {
            _update_leaving_atom_flag();

            if ((!_links.empty() || !_ssbonds.empty()) && !_molecules.empty()) {
                 std::map<std::string, unsigned int> chain_order;
                 chain_order.clear();
                 RCSB::Chain* chain = _molecules[0]->GetFirstChain();
                 while (chain) {
                      if (chain_order.find(chain->PDB_ChainID()) == chain_order.end()) {
                           chain_order.insert(std::make_pair(chain->PDB_ChainID(), chain_order.size()));
                      }
                      chain = _molecules[0]->GetNextChain();
                 }

                 _re_order_links(chain_order, _ssbonds);
                 _re_order_links(chain_order, _links);
            }

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

void ValidateObj::_re_order_links(const std::map<std::string, unsigned int>& chain_order, std::list<_LINK>& links)
{
       if (links.empty()) return;

       int maximum_number = 0;
       for (std::list<_LINK>::const_iterator pos = links.begin(); pos != links.end(); ++pos) {
            if (atoi(pos->fstAtom->pdb_resnum().c_str()) > maximum_number) maximum_number = atoi(pos->fstAtom->pdb_resnum().c_str());
            if (atoi(pos->sndAtom->pdb_resnum().c_str()) > maximum_number) maximum_number = atoi(pos->sndAtom->pdb_resnum().c_str());
       }
       maximum_number += 1;

       std::multimap<int, _LINK> tmp_mapping;
       std::multimap<int, std::multimap<int, _LINK> > tmp_tmp_mapping;

       // first key: covale/metalc/disulf type
       std::map<std::string, std::multimap<int, std::multimap<int, _LINK> > > order_mapping;
       order_mapping.clear();
       for (std::list<_LINK>::iterator pos = links.begin(); pos != links.end(); ++pos) {
            int first_number = atoi(pos->fstAtom->pdb_resnum().c_str());
            std::map<std::string, unsigned int>::const_iterator mpos = chain_order.find(pos->fstAtom->pdb_chnid());
            if (mpos != chain_order.end()) first_number += mpos->second * maximum_number;

            int second_number = atoi(pos->sndAtom->pdb_resnum().c_str());
            mpos = chain_order.find(pos->sndAtom->pdb_chnid());
            if (mpos != chain_order.end()) second_number += mpos->second * maximum_number;

            if (first_number > second_number) {
                 int number = first_number;
                 first_number = second_number;
                 second_number = number;
                 RCSB::Atom* atom = pos->fstAtom;
                 std::string sym = pos->SymOP_1;
                 pos->fstAtom = pos->sndAtom;
                 pos->sndAtom = atom;
                 pos->SymOP_1 = pos->SymOP_2;
                 pos->SymOP_2 = sym;
            }

            std::map<std::string, std::multimap<int, std::multimap<int, _LINK> > >::iterator ompos = order_mapping.find(pos->type);
            if (ompos != order_mapping.end()) {
                 std::multimap<int, std::multimap<int, _LINK> >::iterator ttpos = ompos->second.find(first_number);
                 if (ttpos != ompos->second.end())
                      ttpos->second.insert(std::make_pair(second_number, *pos));
                 else {
                      tmp_mapping.clear();
                      tmp_mapping.insert(std::make_pair(second_number, *pos));
                      ompos->second.insert(std::make_pair(first_number, tmp_mapping));
                 }
            } else {
                 tmp_mapping.clear();
                 tmp_mapping.insert(std::make_pair(second_number, *pos));
                 tmp_tmp_mapping.clear();
                 tmp_tmp_mapping.insert(std::make_pair(first_number, tmp_mapping));
                 order_mapping.insert(std::make_pair(pos->type, tmp_tmp_mapping));
            }
       }

       links.clear();
       for (std::map<std::string, std::multimap<int, std::multimap<int, _LINK> > >::const_iterator
            ompos = order_mapping.begin(); ompos != order_mapping.end(); ++ompos) {
            for (std::multimap<int, std::multimap<int, _LINK> >::const_iterator tmpos = ompos->second.begin(); tmpos != ompos->second.end(); ++tmpos) {
                 for (std::multimap<int, _LINK>::const_iterator mpos = tmpos->second.begin(); mpos != tmpos->second.end(); ++mpos) {
                      links.push_back(mpos->second);
                 }
            }
       }
}

void ValidateObj::_update_struct_conn_table(ISTable* t, std::vector<std::string>& type_array, std::set<std::string>& type_set, const std::list<_LINK>& links)
{
       if (links.empty()) return;

       std::map<std::string, int> id_mapping;
       id_mapping.clear();

       int i = t->GetNumRows();
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
            t->UpdateCell(i, "pdbx_role", pos->pdbx_role);
            UpdateAtomInfo::UpdateTable_3(t, i, pos->fstAtom, "ptnr1_", NUM_ALL_ITEM);
            UpdateAtomInfo::UpdateTable_3(t, i, pos->sndAtom, "ptnr2_", NUM_ALL_ITEM);
            i++;
       }
}

