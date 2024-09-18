/*
FILE:     Ndb2Pdb_Remark_600s.C
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

#include "Element.h"
#include "GeometryUtil.h"
#include "Maxit.h"
#include "NdbToken.h"
#include "PdbWrite.h"
#include "Remark_define.h"
#include "Remark_extern.h"
#include "utillib.h"

typedef struct {
       const char *remark;
       int index;
} _REMARK_630_TEMPLATE;

#define NUM_REMARK_630 10

const static _REMARK_630_TEMPLATE _remark_630_template[NUM_REMARK_630] = {
       { "MOLECULE TYPE: ", 1 },
       { "MOLECULE NAME: ", 2 },
       { "(M=MODEL NUMBER; RES=RESIDUE NAME; C=CHAIN IDENTIFIER;", 0 },
       { " SSSEQ=SEQUENCE NUMBER; I=INSERTION CODE.)", 0 },
       { "", 0 },
       { "  M RES C SSSEQI", 0 },
       { "SOURCE: ", 3 },
       { "TAXONOMY: ", 4 },
       { "SUBCOMP:", 0 },
       { "DETAILS: ", 5 }
};

typedef struct {
       // vector[0]: ResName
       // vector[1]: Type
       // vector[2]: Chemical Name
       // vector[3]: Source Scientific Name
       // vector[4]: Taxonomy ID
       // vector[5]: Details
       std::vector<std::string> data;
       // Subcomponent list from _chem_comp.pdbx_subcomponent_list
       std::vector<std::string> subcomp;
       // vector[0]: Mol_ID
       // vector[1]: ResName
       // vector[2]: PDB Chain ID
       // vector[3]: ResNum
       // vector[4]: InsCode
       std::vector<std::vector<std::string> > instances;
} _REMARK_630;

static bool get_atom_pair(const _LINK& link, CrySymmetry& cell, RCSB::Atom& metal_atom, RCSB::Atom& polar_atom, std::string& SymOP);

void Maxit::_ndb_to_pdb_get_remark_610()
{
       if (_missing_non_polymer_residues.empty()) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_610, Remark_610, 0, 0, 0);

       for (std::map<int, std::list<std::vector<std::string> > >::const_iterator mpos =
            _missing_non_polymer_residues.begin(); mpos != _missing_non_polymer_residues.end(); ++mpos) {
            std::string mol_id = "";
            if (_missing_non_polymer_residues.size() > 1)
                 mol_id = String::IntToString(mpos->first);
            _ndb_to_pdb_get_remark_465_475(mol_id, 610, mpos->second);
       }
}

void Maxit::_ndb_to_pdb_get_remark_615()
{
       if (_zero_occ_non_polymer_residues.empty()) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_615, Remark_615, 0, 0, 0);

       for (std::map<int, std::list<std::vector<std::string> > >::const_iterator mpos =
            _zero_occ_non_polymer_residues.begin(); mpos != _zero_occ_non_polymer_residues.end(); ++mpos) {
            std::string mol_id = "";
            if (_zero_occ_non_polymer_residues.size() > 1)
                 mol_id = String::IntToString(mpos->first);
            _ndb_to_pdb_get_remark_465_475(mol_id, 615, mpos->second);
       }
}

void Maxit::_ndb_to_pdb_get_remark_620()
{
       if (_links.empty() || _molecules.empty()) return;

       // Atom lists: first atom is metal atom, the rest are polar atoms for each record
       std::vector<std::vector<RCSB::Atom> > MetalList;
       MetalList.clear();

       // index for MetalList: key=metal atom name index, val=position in MetalList
       std::map<std::string, unsigned int> MetalIndex;
       MetalIndex.clear();

       // Unique metal atom name index set
       std::set<std::string> atomIndex;
       atomIndex.clear();

       std::string SymOP;
       std::vector<RCSB::Atom*> atoms;
       std::vector<RCSB::Atom> tMetalList;
       RCSB::Atom metal_atom, polar_atom;
       for (std::list<_LINK>::const_iterator lpos = _links.begin(); lpos != _links.end(); ++lpos) {
            if (lpos->type != "metalc") continue;
            if (!get_atom_pair(*lpos, _cell, metal_atom, polar_atom, SymOP)) continue;

            std::string index = metal_atom.getAtomAllIndex();
            std::map<std::string, unsigned int>::const_iterator mpos = MetalIndex.find(index);
            if (mpos != MetalIndex.end())
                 // insert polar atom into existing metal atom cluster
                 MetalList[mpos->second].push_back(polar_atom);
            else {
                 // add new metal atom cluster
                 MetalIndex.insert(std::make_pair(index, MetalList.size()));
                 tMetalList.clear();
                 tMetalList.push_back(metal_atom);
                 tMetalList.push_back(polar_atom);
                 MetalList.push_back(tMetalList);
            }

            mpos = MetalIndex.find(index);
            if (mpos == MetalIndex.end()) continue;

            if (atomIndex.find(index) != atomIndex.end()) continue;

            atomIndex.insert(index);

            RCSB::Residue* res = _molecules[0]->find_pdb_residue(metal_atom.pdb_chnid(), metal_atom.pdb_resnam(), metal_atom.pdb_resnum(),
                                                                 metal_atom.ins_code()); 
            if (!res) continue;

            std::string alt_loc = polar_atom.alt_loc();
            try {
                 // find all linked atoms if metal atom is belongs to a metal-complex residue
                 const ConnectFormat& drug = _ccDic->find_drug(metal_atom.pdb_resnam());
                 std::map<std::string, std::vector<std::string> > LinkedAtoms = drug.getLinkedAtoms();
                 if (LinkedAtoms.empty()) continue;
                 std::map<std::string, std::vector<std::string> >::const_iterator mpos1 = LinkedAtoms.find(metal_atom.pdb_atmnam());
                 if (mpos1 == LinkedAtoms.end()) continue;
                 for (std::vector<std::string>::const_iterator rpos = mpos1->second.begin(); rpos != mpos1->second.end(); ++rpos) {
                      res->find_atom(*rpos, atoms);
                      if (atoms.empty()) continue;
                      for (std::vector<RCSB::Atom*>::const_iterator apos = atoms.begin(); apos != atoms.end(); ++apos) {
                           if ((*apos)->alt_loc() != alt_loc && (*apos)->alt_loc() != metal_atom.alt_loc() &&
                              !(*apos)->alt_loc().empty() && !alt_loc.empty() && !metal_atom.alt_loc().empty()) continue;
                           polar_atom = *(*apos);
                           _cell.symmetry_operation(polar_atom, SymOP);
                           MetalList[mpos->second].push_back(polar_atom);
                           break;
                      }
                 }
            } catch (const std::exception& exc) {}
       }

       bool found = false;
       for (std::vector<std::vector<RCSB::Atom> >::const_iterator pos = MetalList.begin(); pos != MetalList.end(); ++pos) {
            if (pos->size() > 2) {
                 found = true;
                 break;
            }
       }
       if (!found) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_620, Remark_620, 0, 0, 0);

       for (std::vector<std::vector<RCSB::Atom> >::const_iterator pos = MetalList.begin(); pos != MetalList.end(); ++pos) {
            if (pos->size() < 3) continue;

            _addNewRemark(620, "");
            _addNewRemark(620, "COORDINATION ANGLES FOR:  M RES CSSEQI METAL");

            std::string remark = "                            " + FormattedString((*pos)[0].pdb_resnam(), 3) + " "
                   + (*pos)[0].pdb_chnid_char() + FormattedString((*pos)[0].pdb_resnum(), 4, 0) + (*pos)[0].ins_code_char() + " "
                   + FormattedString(printAtomNameField(_ccDic, (*pos)[0].atom_type(), (*pos)[0].pdb_atmnam(), (*pos)[0].pdb_resnam()), 4);
            _addNewRemark(620, remark);
            _addNewRemark(620, "N RES CSSEQI ATOM");

            std::string remark_list = "";
            for (unsigned int i = 1; i < pos->size(); ++i) {
                 if (i == 10) break;
                 remark =  String::IntToString(i) + " " + FormattedString((*pos)[i].pdb_resnam(), 3) + " " + (*pos)[i].pdb_chnid_char()
                        + FormattedString((*pos)[i].pdb_resnum(), 4, 0) + (*pos)[i].ins_code_char() + " "
                        + FormattedString(printAtomNameField(_ccDic, (*pos)[i].atom_type(), (*pos)[i].pdb_atmnam(), (*pos)[i].pdb_resnam()), 4);

                 remark_list =  "N                 ";
                 for (unsigned int j = 1; j < i; ++j) {
                      double angle = cal_angle((*pos)[j], (*pos)[0], (*pos)[i]);
                      remark += FloatToString(angle, 6, 1);
                      remark_list += "   " + String::IntToString(j) + "  ";
                 }
                 _addNewRemark(620, remark);
            }
            _addNewRemark(620, remark_list);
       }
}

void Maxit::_ndb_to_pdb_get_remark_630()
{
       if (_molecules.empty()) return;

       std::vector<_REMARK_630> het_list;
       het_list.clear();

       std::map<std::string, unsigned int> het_index;
       het_index.clear();

       _REMARK_630 het;
       std::string name, type, source, tax_id, details;
       std::vector<std::string> data;
       std::vector<RCSB::Residue*> residues;

       for (std::vector<RCSB::Molecule*>::const_iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            std::string mol_ID = "   ";
            if (_molecules.size() > 1) mol_ID = String::IntToString((*mpos)->Mol_ID());
            RCSB::Chain* chain = (*mpos)->GetFirstChain();
            while (chain) {
                 if (chain->ResidueNumbers() > 1) {
                      chain = (*mpos)->GetNextChain();
                      continue;
                 }
                 chain->GetFirstResidueList(residues);
                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                      try {
                           const ConnectFormat& drug = _ccDic->find_drug((*rpos)->ResName());
                           std::string subcomponent = drug.getMetaData("pdbx_subcomponent_list");
                           std::map<std::string, std::vector<std::string> >::const_iterator
                               mpos1 = _non_polymer_prd_info.find((*rpos)->ResName());
                           if (subcomponent.empty() && mpos1 == _non_polymer_prd_info.end())
                                continue;

                           data.clear();
                           data.push_back(mol_ID);
                           data.push_back((*rpos)->ResName());
                           std::string chnid = " ";
                           if (!(*rpos)->pdb_chnid().empty()) chnid = (*rpos)->pdb_chnid(); 
                           data.push_back(chnid);
                           data.push_back((*rpos)->pdb_res_no());
                           std::string ins_code = " ";
                           if (!(*rpos)->ins_code().empty()) ins_code = (*rpos)->ins_code();
                           data.push_back(ins_code);

                           std::map<std::string, unsigned int>::const_iterator ipos = het_index.find((*rpos)->ResName());
                           if (ipos != het_index.end()) {
                                het_list[ipos->second].instances.push_back(data);
                           } else {
                                het_index.insert(std::make_pair((*rpos)->ResName(), het_list.size()));

                                type.clear();
                                source.clear();
                                tax_id.clear();
                                details.clear();

                                if (mpos1 != _non_polymer_prd_info.end()) {
                                     type = mpos1->second[1];
                                     details = mpos1->second[2];
                                }

                                std::map<int, Entity>::const_iterator mpos2 = _entities.find(chain->int_entity_id());
                                if (mpos2 != _entities.end()) {
                                     const std::vector<std::pair<std::string, std::map<std::string, std::string> > >& sources = mpos2->second.source();
                                     if (!sources.empty()) {
                                          std::map<std::string, std::string>::const_iterator
                                              mpos3 = sources[0].second.find("pdbx_gene_src_scientific_name");
                                          if (mpos3 == sources[0].second.end()) {
                                               mpos3 = sources[0].second.find("pdbx_organism_scientific");
                                               if (mpos3 == sources[0].second.end()) {
                                                    mpos3 = sources[0].second.find("organism_scientific");
                                                    if (mpos3 != sources[0].second.end()) source = mpos3->second;
                                               } else source = mpos3->second;
                                          } else source = mpos3->second;

                                          mpos3 = sources[0].second.find("pdbx_gene_src_ncbi_taxonomy_id");
                                          if (mpos3 == sources[0].second.end()) {
                                               mpos3 = sources[0].second.find("pdbx_ncbi_taxonomy_id");
                                               if (mpos3 == sources[0].second.end()) {
                                                    mpos3 = sources[0].second.find("ncbi_taxonomy_id");
                                                    if (mpos3 != sources[0].second.end())
                                                         tax_id = mpos3->second;
                                               } else tax_id = mpos3->second;
                                          } else tax_id = mpos3->second;
                                     }
                                }

                                het.data.clear();
                                het.subcomp.clear();
                                het.instances.clear();

                                String::UpperCase(drug.chemical_name(), name);
                                String::UpperCase(type);
                                String::UpperCase(source);
                                String::UpperCase(details);

                                het.data.push_back((*rpos)->ResName());
                                het.data.push_back(type);
                                het.data.push_back(name);
                                het.data.push_back(source);
                                het.data.push_back(tax_id);
                                het.data.push_back(details);
                                get_wordarray(het.subcomp, subcomponent, " ");
                                het.instances.push_back(data);
                                het_list.push_back(het);
                           }
                      } catch (const std::exception& exc) {}
                 }
                 chain = (*mpos)->GetNextChain();
            }
       }

       if (het_list.empty()) return;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
       int width = ndbformat.FieldList[2].FieldWidth - 1;

       for (std::vector<_REMARK_630>::const_iterator vpos = het_list.begin(); vpos != het_list.end(); ++vpos) {
            for (int i = 0; i < NUM_REMARK_630; ++i) {
                 std::string remark = _remark_630_template[i].remark;
                 if (remark == "  M RES C SSSEQI") {
                      _addNewRemark(630, remark);
                      for (std::vector<std::vector<std::string> >::const_iterator vvpos =
                           vpos->instances.begin(); vvpos != vpos->instances.end(); ++vvpos) {
                           remark = FormattedString((*vvpos)[0], 3) + " " + FormattedString((*vvpos)[1], 3) + " "
                                  + (*vvpos)[2] + " " + FormattedString((*vvpos)[3], 5) + (*vvpos)[4];
                           _addNewRemark(630, remark);
                      }
                 } else if (remark == "SUBCOMP:") {
                      if (!vpos->subcomp.empty()) {
                           remark += "    ";
                           int count = 0;
                           for (std::vector<std::string>::const_iterator vvpos = vpos->subcomp.begin(); vvpos != vpos->subcomp.end(); ++vvpos) {
                                if (count == 10) {
                                     _addNewRemark(630, remark);
                                     count = 0;
                                     remark = "          2 ";
                                }
                                remark += FormattedString(*vvpos, 3) + " ";
                                count++;
                           }
                      } else remark += " NULL";
                      _addNewRemark(630, remark);
                 } else if (_remark_630_template[i].index) {
                      if (!vpos->data[_remark_630_template[i].index].empty())
                           remark += vpos->data[_remark_630_template[i].index];
                      else remark += "NULL";
                      get_max_length_words(data, remark, width);
                      _ndb_to_pdb_add_remark(630, data);
                 } else _addNewRemark(630, remark);
            }
       }
}

static bool get_atom_pair(const _LINK& link, CrySymmetry& cell, RCSB::Atom& metal_atom, RCSB::Atom& polar_atom, std::string& SymOP)
{
       if (link.mol_index != 0) return false;

       if (Element::findMetalFlag(link.fstAtom->atom_type())) {
            SymOP = link.SymOP_1;
            metal_atom = *(link.fstAtom);
            cell.symmetry_operation(metal_atom, link.SymOP_1);
            polar_atom = *(link.sndAtom);
            cell.symmetry_operation(polar_atom, link.SymOP_2);
       } else if (Element::findMetalFlag(link.sndAtom->atom_type())) {
            SymOP = link.SymOP_2;
            metal_atom = *(link.sndAtom);
            cell.symmetry_operation(metal_atom, link.SymOP_2);
            polar_atom = *(link.fstAtom);
            cell.symmetry_operation(polar_atom, link.SymOP_1);
       } else return false;

       return true;
}
