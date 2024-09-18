/*
FILE:     FileObj_Util.C
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
#include <math.h>

#include "CompositeIndex.h"
#include "FileObj.h"
#include "utillib.h"

std::string FileObj::_insert_symm_matrix(std::vector<SymmMatrix>& matrices)
{
       std::vector<std::string> ids;
       ids.clear();
       std::set<std::string> id_set;
       id_set.clear();

       for (std::vector<SymmMatrix>::iterator pos = matrices.begin(); pos != matrices.end(); ++pos) {
            std::string id = _insert_symm_matrix(*pos);
            if (id_set.find(id) != id_set.end()) continue;
            ids.push_back(id);
            id_set.insert(id);
       }
       std::string id_list = "";
       for (std::vector<std::string>::const_iterator pos = ids.begin(); pos != ids.end(); ++pos) {
            if (!id_list.empty()) id_list += ",";
            id_list += *pos;
       }
       return id_list;
}

string FileObj::_insert_symm_matrix(SymmMatrix& symmtrix)
{
       for (unsigned int i = 0; i < _SymmMatrices.size(); ++i) {
            if (_SymmMatrices[i] == symmtrix)
                 return _SymmMatrices[i].getValue("id");
       }

       int i = _SymmMatrices.size() + 1;
       std::string id = String::IntToString(i);
       while (_SymmMatrix_Mapping.find(id) != _SymmMatrix_Mapping.end()) {
            i++;
            id = String::IntToString(i);
       }

       _update_symm_matrix(symmtrix);
       symmtrix.setValue("id", id);
       _SymmMatrix_Mapping.insert(std::make_pair(id, _SymmMatrices.size()));
       _SymmMatrices.push_back(symmtrix);
       return id;
}

void FileObj::_update_symm_matrix(SymmMatrix& symmtrix)
{
       double maxtrix[3][3], vect[3];
       if (_cell.is_artifical()) return;
       if (!symmtrix.getMatrix(maxtrix)) return;
       if (!symmtrix.getVector(vect)) return;

       std::string sym_oper, sym_as_xyz;
       _cell.ndb_get_symmetry_operation(maxtrix, vect, sym_oper);
       if (sym_oper.empty()) return;

       symmtrix.setValue("name", sym_oper);
       if (sym_oper == "1_555")
            symmtrix.setValue("type", "identity operation");
       else symmtrix.setValue("type", "crystal symmetry operation");
       _cell.ndb_get_symmetry_operation_name(sym_oper, sym_as_xyz);
       symmtrix.setValue("symmetry_operation", sym_as_xyz);

       _cell.ndb_get_symmetry_matrix(sym_oper, maxtrix, vect);
       for (int i = 0; i < 3; ++i) {
            std::string item = "vector[" + String::IntToString(i + 1) + "]"; 
            symmtrix.setValue(item, FloatToString(vect[i]));
            for (int j = 0; j < 3; ++j) {
                 item = "matrix[" + String::IntToString(i + 1)
                      + "][" + String::IntToString(j + 1) + "]";
                 symmtrix.setValue(item, FloatToString(maxtrix[i][j]));
            }
       }
}

void FileObj::_read_struct_assembly_categories(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_struct_assembly");
       if (t == NULL) return;

       _assemblies.clear();

       std::map<std::string, std::string> item_mapping;
       item_mapping.clear(); 
       item_mapping.insert(std::make_pair("id", "id"));
       item_mapping.insert(std::make_pair("details", "method"));
       item_mapping.insert(std::make_pair("method_details", "software"));
       item_mapping.insert(std::make_pair("oligomeric_details", "details"));
       item_mapping.insert(std::make_pair("oligomeric_count", "oligomeric"));

       Assembly assembly;
       std::string cs1, cs2, cs3, cs4;

       std::map<std::string, unsigned int> assembly_index;
       assembly_index.clear();

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            assembly.clear();
            for (std::map<std::string, std::string>::const_iterator mpos = item_mapping.begin(); mpos != item_mapping.end(); ++mpos) {
                 get_value_clean(cs1, t, i, mpos->first);
                 if (cs1.empty()) continue;
                 assembly.setValue(mpos->second, cs1);
            }
            cs1 = assembly.getValue("id");
            if (cs1.empty()) continue;

            assembly_index.insert(std::make_pair(cs1, _assemblies.size()));
            _assemblies.push_back(assembly);
       }

       t = getTablePtr(block, "pdbx_struct_assembly_prop");
       if (t) {
            std::string cs1, cs2, cs3;
            rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs1, t, i, "biol_id");
                 get_value_clean(cs2, t, i, "type");
                 get_value_clean(cs3, t, i, "value");
                 if (cs1.empty() || cs2.empty() || cs3.empty()) continue;
                 std::map<std::string, unsigned int>::const_iterator mpos = assembly_index.find(cs1);
                 if (mpos == assembly_index.end()) continue;

                 if (cs2 == "ABSA (A^2)")
                      _assemblies[mpos->second].setValue("absa", cs3);
                 else if (cs2 == "SSA (A^2)")
                      _assemblies[mpos->second].setValue("ssa", cs3);
                 else if (cs2 == "MORE")
                      _assemblies[mpos->second].setValue("more", cs3);
            }
       }

       t = getTablePtr(block, "pdbx_struct_assembly_gen");
       if (!t) {
            _assemblies.clear();
            return;
       }

       rowNo = t->GetNumRows();
       std::vector<std::string> pdb_chain_ids;
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs1, t, i, "assembly_id");
            get_value_clean(cs2, t, i, "oper_expression");
            get_value_clean(cs3, t, i, "auth_asym_id_list");
            get_value_clean(cs4, t, i, "asym_id_list");
            if (cs1.empty() || cs2.empty()) continue;
            std::map<std::string, unsigned int>::const_iterator mpos = assembly_index.find(cs1);
            if (mpos == assembly_index.end()) continue;

            if (!cs3.empty() && cs4.empty()) {
                 get_wordarray(pdb_chain_ids, cs3, ";, ");
                 _assemblies[mpos->second].InsertChains(cs2, pdb_chain_ids);
            } else if (!cs4.empty()) {
                 _get_chain_ids_from_asym_id_list(cs4, pdb_chain_ids);
                 if (!pdb_chain_ids.empty())
                      _assemblies[mpos->second].InsertChains(cs2, pdb_chain_ids);
            }
       }

       std::set<std::string> linear_polymeric_chain_id_set;
       linear_polymeric_chain_id_set.clear();
       if (!_molecules.empty()) {
            RCSB::Chain* chain = _molecules[0]->GetFirstChain();
            while (chain) {
                 if ((chain->chain_type() == "ATOMP") || (chain->chain_type() == "ATOMN")) linear_polymeric_chain_id_set.insert(chain->PDB_ChainID());
                 chain = _molecules[0]->GetNextChain();
            }
       }
       for (std::vector<Assembly>::iterator pos = _assemblies.begin(); pos != _assemblies.end(); ++pos) {
            pos->UpdateOligomericStatus(linear_polymeric_chain_id_set);
       }
}

void FileObj::_read_pdbx_struct_oper_list(Block& block)
{
       ISTable *t = getTablePtr(block, "pdbx_struct_oper_list");
       if (t == NULL) return;

       _SymmMatrices.clear();
       _SymmMatrix_Mapping.clear();

       SymmMatrix symmtrix;
       std::string cs;

       const std::vector<std::string>& itemNames = t->GetColumnNames();
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            symmtrix.clear();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 get_value_clean(cs, t, i, *pos);
                 symmtrix.setValue(*pos, cs);
            }
            std::string sym_oper = symmtrix.getValue("name");
            std::string type = symmtrix.getValue("type");
            std::string sym_as_xyz = symmtrix.getValue("symmetry_operation");
            if (!sym_oper.empty() && (type == "identity operation" || type == "crystal symmetry operation") && sym_as_xyz.empty()) {
                 _cell.ndb_get_symmetry_operation_name(sym_oper, sym_as_xyz);
                 symmtrix.setValue("symmetry_operation", sym_as_xyz);
            }
            if (sym_oper.empty()) symmtrix.setIdentityOperation();
            _SymmMatrix_Mapping.insert(std::make_pair(symmtrix.getValue("id"), _SymmMatrices.size()));
            _SymmMatrices.push_back(symmtrix);
       }
}

void FileObj::_write_struct_assembly_categories(Block& block)
{
       if (_assemblies.empty() || _molecules.empty()) {
            deleteTable(block, "pdbx_struct_assembly");
            deleteTable(block, "pdbx_struct_assembly_gen");
            deleteTable(block, "pdbx_struct_assembly_prop");
            return;
       }

       ISTable *a = _newTablePtr("pdbx_struct_assembly");
       ISTable *g = _newTablePtr("pdbx_struct_assembly_gen");
       ISTable *p = _newTablePtr("pdbx_struct_assembly_prop");

       std::map<std::string, std::string> item_mapping, type_mapping;
       item_mapping.clear();
       item_mapping.insert(std::make_pair("details", "method"));
       item_mapping.insert(std::make_pair("method_details", "software"));
       item_mapping.insert(std::make_pair("oligomeric_details", "details"));
       item_mapping.insert(std::make_pair("oligomeric_count", "oligomeric"));

       type_mapping.clear();
       type_mapping.insert(std::make_pair("ABSA (A^2)", "absa"));
       type_mapping.insert(std::make_pair("SSA (A^2)", "ssa"));
       type_mapping.insert(std::make_pair("MORE", "more"));

       int a_id = 0;
       int arow = 0;
       int grow = 0;
       int prow = 0;
       for (std::vector<Assembly>::const_iterator apos = _assemblies.begin(); apos != _assemblies.end(); ++apos) {
            a_id++;
            std::string assembly_id = String::IntToString(a_id);

            a->AddRow();
            a->UpdateCell(arow, "id", assembly_id);
            for (std::map<std::string, std::string>::const_iterator mpos = item_mapping.begin(); mpos != item_mapping.end(); ++mpos) {
                 a->UpdateCell(arow, mpos->first, apos->getValue(mpos->second));
            }
            arow++;

            for (std::map<std::string, std::string>::const_iterator mpos = type_mapping.begin(); mpos != type_mapping.end(); ++mpos) {
                 std::string value = apos->getValue(mpos->second);
                 if (value.empty()) continue;

                 p->AddRow();
                 p->UpdateCell(prow, "biol_id", assembly_id);
                 p->UpdateCell(prow, "type", mpos->first);
                 p->UpdateCell(prow, "value", value);
                 prow++;
            }

            const std::vector<std::pair<std::string, std::vector<std::string> > >& assembly_chains = apos->assembly_chains();
            for (std::vector<std::pair<std::string, std::vector<std::string> > >::const_iterator
                 cpos = assembly_chains.begin(); cpos != assembly_chains.end(); ++cpos) {
                 std::string asym_id_list = _get_asym_id_list_from_chain_ids(cpos->second);
                 if (asym_id_list.empty()) continue;
                 g->AddRow();
                 g->UpdateCell(grow, "assembly_id", assembly_id);
                 g->UpdateCell(grow, "oper_expression", cpos->first);
                 g->UpdateCell(grow, "asym_id_list", asym_id_list);
                 grow++;
            }
       }

       block.WriteTable(a);

       if (grow  == 0) {
            delete g;
            deleteTable(block, "pdbx_struct_assembly_gen");
       } else block.WriteTable(g);

       if (prow  == 0) {
            delete p;
            deleteTable(block, "pdbx_struct_assembly_prop");
       } else block.WriteTable(p);
}

void FileObj::_write_pdbx_struct_oper_list(Block& block)
{
       if (_SymmMatrices.empty()) {
            deleteTable(block, "pdbx_struct_oper_list");
            return;
       }

       ISTable *t = _newTablePtr("pdbx_struct_oper_list");
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       int row = 0;
       for (std::vector<SymmMatrix>::const_iterator pos = _SymmMatrices.begin(); pos != _SymmMatrices.end(); ++pos) {
            t->AddRow();
            for (std::vector<std::string>::const_iterator vpos = itemNames.begin(); vpos != itemNames.end(); ++vpos) {
                 t->UpdateCell(row, *vpos, pos->getValue(*vpos));
            }
            row++;
       }
       block.WriteTable(t);
}

std::string FileObj::_get_asym_id_list_from_chain_ids(const std::vector<std::string>& chain_id_list)
{
       if (chain_id_list.empty()) return "";

       if (_molecules.empty()) return "";

       std::multimap<int, std::string> order_asym_ids;
       order_asym_ids.clear();

       std::set<std::string> unique_set;
       unique_set.clear();

       std::vector<RCSB::Chain*> chain_list;
       for (std::vector<std::string>::const_iterator pos = chain_id_list.begin(); pos != chain_id_list.end(); ++pos) {
            _molecules[0]->GetPdbChains(*pos, chain_list);
            if (chain_list.empty()) continue;
            for (std::vector<RCSB::Chain*>::const_iterator cpos = chain_list.begin(); cpos != chain_list.end(); ++cpos) {
                 if ((*cpos)->chain_type() == "ATOMS") {
                      RCSB::Residue* res = (*cpos)->GetFirstResidue();
                      while (res) {
                           if (unique_set.find(res->chnid()) == unique_set.end()) {
                                order_asym_ids.insert(std::make_pair((*cpos)->order(), res->chnid()));
                                unique_set.insert(res->chnid());
                           }
                           res = (*cpos)->GetNextResidue();
                      }
                 } else if (unique_set.find((*cpos)->ChainID()) == unique_set.end()) {
                      unique_set.insert((*cpos)->ChainID());
                      order_asym_ids.insert(std::make_pair((*cpos)->order(), (*cpos)->ChainID()));
                 }
            }
       }

       std::string asym_id_list = "";
       for (std::multimap<int, std::string>::const_iterator mpos = order_asym_ids.begin(); mpos != order_asym_ids.end(); ++mpos) {
            if (!asym_id_list.empty()) asym_id_list += ",";
            asym_id_list += mpos->second;
       }

       return asym_id_list;
}

void FileObj::_get_chain_ids_from_asym_id_list(const std::string& asym_id_list, std::vector<std::string>& chain_id_list)
{
       chain_id_list.clear();

       if (asym_id_list.empty()) return;

       if (_molecules.empty()) return;

       std::set<std::string> chain_id_set;
       chain_id_set.clear();

       std::vector<std::string> data;
       get_wordarray(data, asym_id_list, ";, ");
       for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            RCSB::Chain* chain = _molecules[0]->GetAsymChain(*pos);
            if (!chain) continue;
            if (chain_id_set.find(chain->PDB_ChainID()) != chain_id_set.end()) continue;
            chain_id_set.insert(chain->PDB_ChainID());
            chain_id_list.push_back(chain->PDB_ChainID());
       }
}

void FileObj::_check_atom_number_count(std::vector<int>& atom_counts)
{
       atom_counts.clear();

       if (_molecules.empty()) return;

       int atom_no = _molecules[0]->count_atom_no("ATOMP");
       atom_counts.push_back(atom_no);

       int total = atom_no;

       atom_no = _molecules[0]->count_atom_no("ATOMN");
       atom_counts.push_back(atom_no);
       total += atom_no;

       atom_no = _molecules[0]->count_atom_no("HETAD") + _molecules[0]->count_atom_no("HETAC")
               + _molecules[0]->count_atom_no("HETAI") + _molecules[0]->count_atom_no("HETIC")
               + _molecules[0]->count_atom_no("HETAIN") + _molecules[0]->count_atom_no("ATOMS");
       atom_counts.push_back(atom_no);
       total += atom_no;

       atom_no = _molecules[0]->count_atom_no("HETAS");
       atom_counts.push_back(atom_no);
       total += atom_no;

       atom_counts.push_back(total);
}

bool FileObj::_readInstanceAssignFile(const std::string& assignfile, std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > >&
                                      atom_mapping, const bool& remove_mol_id_flag)
{
       // key: Instance ID
       // pair.first: Reference Comp_ID
       // pair.second: Instance vs. Reference atom mapping;
       atom_mapping.clear();

       CifFile *fobj = get_fobj(_logIo, "", assignfile);
       if (!fobj) return false;

       Block& block = fobj->GetBlock(fobj->GetFirstBlockName());

       std::map<std::string, std::string> selectedMatch, tmp_map;
       selectedMatch.clear();

       std::set<std::string> index_set;
       index_set.clear();

       std::string cs1, cs2;

       ISTable *t = getTablePtr(block, "pdbx_instance_assignment");
       if (!t) {
            _logIo->message("category 'pdbx_instance_assignment' not found.\n");
            delete fobj;
            return false;
       }

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs1, t, i, "inst_id");
            if (cs1.empty()) continue;
            get_value_clean(cs2, t, i, "het_id");
            if (cs2.empty()) continue;
            selectedMatch.insert(std::make_pair(cs1, cs2));
            index_set.insert(cs1 + "_" + cs2);
       }
       if (selectedMatch.empty()) {
            _logIo->message("No assignment found.\n");
            delete fobj;
            return false;
       }

       t = getTablePtr(block, "pdbx_match_list");
       if (t) {
            tmp_map.clear();
            rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(cs1, t, i, "inst_id");
                 if (cs1.empty()) continue;
                 if (tmp_map.find(cs1) != tmp_map.end()) continue;

                 get_value_clean(cs2, t, i, "reference_id");
                 if (cs2.empty()) continue;

                 tmp_map.insert(std::make_pair(cs1, cs2));
            }

            index_set.clear();
            for (std::map<std::string, std::string>::iterator mpos = selectedMatch.begin(); mpos != selectedMatch.end(); ++mpos) {
                 std::map<std::string, std::string>::const_iterator mpos1 = tmp_map.find(mpos->first);
                 if (mpos1 != tmp_map.end()) mpos->second = mpos1->second;
                 index_set.insert(mpos->first + "_" + mpos->second);
            }
       }

       t = getTablePtr(block, "pdbx_atom_mapping");
       if (!t) {
            _logIo->message("category 'pdbx_atom_mapping' not found.\n");
            delete fobj;
            return false;
       }

       std::string cs3, cs4;

       rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(cs1, t, i, "inst_id");
            if (cs1.empty()) continue;
            std::map<std::string, std::string>::const_iterator spos = selectedMatch.find(cs1);
            if (spos == selectedMatch.end()) continue;

            get_value_clean(cs2, t, i, "reference_id");
            if (cs2 != spos->second) continue;

            get_value_clean(cs3, t, i, "inst_atom_name");
            if (cs3.empty()) continue;

            get_value_clean(cs4, t, i, "reference_atom_name");
            if (cs4.empty()) continue;

            if (index_set.find(cs1 + "_" + cs2) == index_set.end()) continue;

            std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > >::iterator mpos = atom_mapping.find(cs1);
            if (mpos != atom_mapping.end())
                 mpos->second.second.insert(std::make_pair(cs3, cs4));
            else {
                 tmp_map.clear();
                 tmp_map.insert(std::make_pair(cs3, cs4));
                 atom_mapping.insert(std::make_pair(cs1, std::make_pair(cs2, tmp_map)));
            }
       }

       delete fobj;

       if (atom_mapping.empty()) {
            _logIo->message("No atom mapping found.\n");
            return false;
       }

       bool found_error = false;
       for (std::map<std::string, std::string>::const_iterator mpos = selectedMatch.begin(); mpos != selectedMatch.end(); ++mpos) {
            if (atom_mapping.find(mpos->first) == atom_mapping.end()) {
                 found_error = true;
                 std::string error = "No atom mapping found for instance " + mpos->first + ".\n";
                 _logIo->message(error.c_str());
            }
       }
       if (found_error) return false;

       if (remove_mol_id_flag) {
            std::vector<std::string> data /* , data1 */ ;
            std::map<std::string, std::string> tmp_map;
            // key: PDB Chain ID
            // value_map.key: Residue Name_Residue Number_Insertion Code
            // value_map.value: New Residue Name
            std::map<std::string, std::map<std::string, std::string> > residue_name_changing_mapping;
            residue_name_changing_mapping.clear();
            std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > > tmp_mapping;
            tmp_mapping.clear();
            for (std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > >::iterator
                 mpos = atom_mapping.begin(); mpos != atom_mapping.end(); ++mpos) {
/*
                 get_wordarray(data, mpos->first, "_");
                 data1.clear();
                 for (unsigned int i = 1; i < data.size(); ++i) data1.push_back(data[i]);
                 if (data1.size() == 3) data1.push_back("");
                 std::string index = CompositeIndex::getIndex(data1);
*/
                 get_wordarray_delimit_by_string(data, mpos->first, "_");
                 std::string index = CompositeIndex::getIndex(data[1], data[2], data[3], data[4]);
                 tmp_mapping.insert(std::make_pair(index, mpos->second));
                 if (data[2] != mpos->second.first) {
                      std::map<std::string, std::map<std::string, std::string> >::iterator mpos1 = residue_name_changing_mapping.find(data[1]);
                      if (mpos1 != residue_name_changing_mapping.end()) {
                           mpos1->second.insert(std::make_pair(CompositeIndex::getIndex(data[2], data[3], data[4]), mpos->second.first));
                      } else {
                           tmp_map.clear();
                           tmp_map.insert(std::make_pair(CompositeIndex::getIndex(data[2], data[3], data[4]), mpos->second.first));
                           residue_name_changing_mapping.insert(std::make_pair(data[1], tmp_map));
                      }
                 }
            }
            atom_mapping = tmp_mapping;

            for (std::map<std::string, std::map<std::string, std::string> >::const_iterator mpos = residue_name_changing_mapping.begin();
                 mpos != residue_name_changing_mapping.end(); ++mpos) {
                 std::map<std::string, std::vector<std::vector<std::string> > >::iterator mspos = _Seq_Scheme_Mapping.find(mpos->first); 
                 if (mspos == _Seq_Scheme_Mapping.end()) continue;

                 for (std::vector<std::vector<std::string> >::iterator vpos = mspos->second.begin(); vpos != mspos->second.end(); ++vpos) {
                      std::map<std::string, std::string>::const_iterator mmpos = mpos->second.find(CompositeIndex::getIndex((*vpos)[2], (*vpos)[3], (*vpos)[4]));
                      if (mmpos != mpos->second.end()) {
                           (*vpos)[1] = mmpos->second;
                           (*vpos)[2] = mmpos->second;
                      }
                 }
            }
       }

       return true;
}
