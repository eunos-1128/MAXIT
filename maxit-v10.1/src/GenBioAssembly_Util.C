/*
FILE:     GenBioAssembly_Util.C
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
#include <string.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <math.h>

#include <set>
#include <string>
#include <vector>

#include "Assembly_Parser_Util.h"
#include "CifCoordWrite.h"
#include "GenBioAssembly.h"
#include "PdbWrite.h"
#include "SeqCodeUtil.h"
#include "TypeDef.h"
#include "utillib.h"

GenBioAssembly::GenBioAssembly(): Maxit()
{
       _has_symmetry_related_operation = false;
       _remove_non_public_item_flag = false;
}

GenBioAssembly::~GenBioAssembly() {}

void GenBioAssembly::write_biological_pdb_files(const std::string& rootname, const std::string& inputfile,  const std::string& indexfile)
{
       if (_molecules.empty()) {
            _logIo->message("No coordinates\n");
            return;
       }

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);
       if (!getTablePtr(_cifblock, "pdbx_struct_assembly") || !getTablePtr(_cifblock, "pdbx_struct_assembly_gen")) {
            _logIo->message("No assembly information\n");
            return;
       }

       std::vector<std::string> removed_tokens, data, id_list, assembly_files;
       removed_tokens.clear();
       removed_tokens.push_back("CAVEAT");
       removed_tokens.push_back("COMPND");
       removed_tokens.push_back("SOURCE");
       removed_tokens.push_back("REVDAT");
       removed_tokens.push_back("SPRSDE");
       removed_tokens.push_back("DBREF");
       removed_tokens.push_back("SEQADV");
       removed_tokens.push_back("MODRES");
       removed_tokens.push_back("HET");
       removed_tokens.push_back("MASTER");

       for (std::vector<std::string>::const_iterator pos = removed_tokens.begin(); pos != removed_tokens.end(); ++pos) {
            _pdb_records.erase(*pos);
       }

       _updateRecordFront("HEADER", 3, "XXXX");

       std::string filename_root = rootname;
       if (filename_root.empty()) {
            // using pdb_id as filename root 
            _getRecordFront("PDBFIL", 1, filename_root);
            if (!filename_root.empty())
                 String::LowerCase(filename_root);
            else {
            // using inputfile file name root
                 get_wordarray(data, inputfile, "/");
                 filename_root = data[data.size() - 1];
                 get_wordarray(data, filename_root, ".");
                 filename_root = data[0];
            }
       }

       PdbWrite writer;
       writer.setLog(_logIo);
       writer.setCCDic(_ccDic);
       writer.setOutputFormat(NDB_FILE_FORMAT_PDB);

       double mat[4][4];
       std::vector<RCSB::Chain*> original_chains;
       std::set<std::string> chain_id_set, entity_id_set;

       NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_PDB);

       assembly_files.clear();
       int assembly_id = 0;
       for (std::vector<Assembly>::const_iterator apos = _assemblies.begin(); apos != _assemblies.end(); ++apos) {
            std::string method = apos->getValue("method");
            if (method != "representative helical assembly" && method != "complete icosahedral assembly" && method != "complete point assembly" &&
                method != "author_and_software_defined_assembly" && method != "author_defined_assembly" && method != "software_defined_assembly") continue;

            assembly_id++;
            std::string filename = filename_root + ".pdb" + String::IntToString(assembly_id);
            assembly_files.push_back(filename);
            FILE *fp = fopen(filename.c_str(), "w");
 
            std::map<std::string, std::list<std::vector<std::string> > > records = _pdb_records;
            writer.writeGeneralRecords(fp, &records);

            chain_id_set.clear();
            entity_id_set.clear();

            int mol_id = 1;
            for (std::vector<std::pair<std::string, std::vector<std::string> > >::const_iterator
                 cpos = apos->assembly_chains().begin(); cpos != apos->assembly_chains().end(); ++cpos) {
                 data.clear();
                 if (parseString(cpos->first, data) != 0) continue;

                 _get_original_chains(cpos->second, original_chains, chain_id_set, entity_id_set);
                 if (original_chains.empty()) continue;

                 for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
                      get_wordarray(id_list, *pos, " ");
                      if (!getFinalMatrix(mat, id_list, _SymmMatrices, _SymmMatrix_Mapping, _cell)) continue;
                 
                      RCSB::Molecule* mol = _get_bio_molecule(original_chains, mat);
                      if (!mol) continue;

                      fprintf(fp, "MODEL     %4d                          ", mol_id);
                      fprintf(fp, "                                        \n");
                      writer.ResetAtomCounting();
                      writer.WriteMolecule(fp, mol);
                      fprintf(fp, "ENDMDL                                  ");
                      fprintf(fp, "                                        \n");
                      mol_id++;

                      delete mol;
                 }
            }

            writer.writeMasterRecord(fp);

            fclose(fp);
       }

       NdbToken::resetAtomTokenField(NDB_FILE_FORMAT_NDB);

       if (assembly_files.empty() || indexfile.empty()) return;

       FILE *fp = fopen(indexfile.c_str(), "w");
       for (std::vector<std::string>::const_iterator pos = assembly_files.begin(); pos != assembly_files.end(); ++pos) {
            fprintf(fp, "%s\n", pos->c_str());
       }
       fclose(fp);
}

void GenBioAssembly::write_biological_cif_files(const std::string& rootname, const std::string& indexfile, const bool& public_flag,
                                                const bool& single_model_flag)
{
       if (_molecules.empty()) {
            _logIo->message("No coordinates\n");
            return;
       }

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);
       if (!getTablePtr(_cifblock, "pdbx_struct_assembly") || !getTablePtr(_cifblock, "pdbx_struct_assembly_gen")) {
            _logIo->message("No assembly information\n");
            return;
       }

       std::vector<std::string> copy_categories, data, id_list, order_list;

       std::string order_file = _rcsbroot  + "/data/ascii/order_list";

       if (public_flag && _dictUtil.Read()) _remove_non_public_item_flag = true;

       copy_categories.clear();
       copy_categories.push_back("audit_author");
       copy_categories.push_back("struct");
       copy_categories.push_back("struct_keywords");
       copy_categories.push_back("exptl");
       copy_categories.push_back("citation");
       copy_categories.push_back("citation_author");
       // copy_categories.push_back("cell");
       // copy_categories.push_back("symmetry");
       copy_categories.push_back("chem_comp");
       copy_categories.push_back("atom_type");
       copy_categories.push_back("atom_sites");
       if (!single_model_flag) {
            copy_categories.push_back("struct_conf");
            copy_categories.push_back("struct_conf_type");
            copy_categories.push_back("struct_conn");
            copy_categories.push_back("struct_conn_type");
            copy_categories.push_back("struct_sheet");
            copy_categories.push_back("struct_sheet_order");
            copy_categories.push_back("struct_sheet_range");
            copy_categories.push_back("pdbx_struct_sheet_hbond");
       }

       CifFile *fobj = create_fobj("", "XXXX");
       Block& block = fobj->GetBlock("XXXX");
       for (std::vector<std::string>::const_iterator pos = copy_categories.begin(); pos != copy_categories.end(); ++pos) {
            ISTable *t = getTableCopy(_cifblock, *pos);
            if (!t) continue;
            if (_remove_non_public_item_flag) {
                 std::vector<std::string> ColumnNames = t->GetColumnNames();
                 for (std::vector<std::string>::const_iterator cpos = ColumnNames.begin(); cpos != ColumnNames.end(); ++cpos) {
                      if (_dictUtil.isPublicItem(*pos, *cpos)) continue;
                      t->DeleteColumn(*cpos);
                 }
            }
            block.WriteTable(t);
       }
       reset_entry_id(block, "XXXX");

       FILE *fp = fopen(indexfile.c_str(), "w");

       double mat[4][4];
       std::vector<RCSB::Chain*> original_chains;
       std::set<std::string> chain_id_set, entity_id_set;
       std::vector<std::map<std::string, std::string> > values, chain_mapping_values;
       std::vector<std::set<int> > chain_index_list;
 
       // key: chain index
       // value[0]: asym ID
       // value[1]: chain ID
       // value[2]: symmetric operators
       std::map<int, std::vector<std::string> > original_asym_chain_id_mapping;

       int assembly_id = 0;
       for (std::vector<Assembly>::const_iterator apos = _assemblies.begin(); apos != _assemblies.end(); ++apos) {
            std::string method = apos->getValue("method");
            if (method != "representative helical assembly" && method != "complete icosahedral assembly" && method != "complete point assembly" &&
                method != "author_and_software_defined_assembly" && method != "author_defined_assembly" && method != "software_defined_assembly") continue;

            ISTable* t1 = _newTablePtr("atom_site");
            if (public_flag) {
                 t1->DeleteColumn("pdbx_auth_seq_id");
                 t1->DeleteColumn("pdbx_auth_comp_id");
                 t1->DeleteColumn("pdbx_auth_asym_id");
                 t1->DeleteColumn("pdbx_auth_atom_name");
            }
            ISTable* t2 = _newTablePtr("atom_site_anisotrop");
            if (public_flag) {
                 t2->DeleteColumn("pdbx_PDB_residue_no");
                 t2->DeleteColumn("pdbx_PDB_residue_name");
                 t2->DeleteColumn("pdbx_PDB_strand_id");
                 t2->DeleteColumn("pdbx_PDB_atom_name");
            }

            chain_id_set.clear();
            entity_id_set.clear();

            chain_mapping_values.clear();

            RCSB::Molecule* single_mol = NULL;
            if (single_model_flag) {
                 single_mol = new RCSB::Molecule;
                 single_mol->setCCDic(_ccDic);
                 single_mol->setMessage(&_error_messages);
                 single_mol->setLog(_logIo);
                 chain_index_list.clear();
                 original_asym_chain_id_mapping.clear();
            }

            int mol_id = 1;
            int atomSerialNo = 0;
            for (std::vector<std::pair<std::string, std::vector<std::string> > >::const_iterator
                 cpos = apos->assembly_chains().begin(); cpos != apos->assembly_chains().end(); ++cpos) {
                 data.clear();
                 if (parseString(cpos->first, data) != 0) continue;

                 _get_original_chains(cpos->second, original_chains, chain_id_set, entity_id_set);
                 if (original_chains.empty()) continue;

                 for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
                      get_wordarray(id_list, *pos, " ");
                      if (!getFinalMatrix(mat, id_list, _SymmMatrices, _SymmMatrix_Mapping, _cell)) continue;

                      if (single_model_flag) {
                           _insert_assembly_chains(true, original_chains, *pos, mat, single_mol, chain_index_list, original_asym_chain_id_mapping);
                      } else {
                           RCSB::Molecule* mol = _get_bio_molecule(original_chains, mat);
                           if (!mol) continue;

                           _update_atom_site_categories(mol_id, mol, atomSerialNo, t1, t2);
                           mol_id++;
                           delete mol;
                      }
                 }
            }

            if (single_model_flag && single_mol) {
                 single_mol->InternalOrder();
                 // single_mol->AssignAsymId();

                 _update_non_atom_site_categories(_cifblock, original_asym_chain_id_mapping, block, single_mol);
                 _update_struct_conf_categories(_cifblock, original_asym_chain_id_mapping, chain_id_set, chain_index_list, single_mol, block);
                 _update_struct_sheet_categories(_cifblock, original_asym_chain_id_mapping, chain_id_set, chain_index_list, single_mol, block);
                 // _update_struct_conn_category(_cifblock, original_asym_chain_id_mapping, chain_index_list, block, single_mol);
                 _update_atom_site_categories(mol_id, single_mol, atomSerialNo, t1, t2);

                 ISTable *t = _newTablePtr("struct_asym");
                 if (public_flag) {
                      t->DeleteColumn("pdbx_PDB_id");
                      t->DeleteColumn("pdbx_alt_id");
                      t->DeleteColumn("pdbx_type");
                      t->DeleteColumn("pdbx_order");
                 }
                 int row = 0;
                 RCSB::Chain* chain = single_mol->GetFirstChain();
                 while (chain) {
                      t->AddRow();
                      t->UpdateCell(row, "id", chain->ChainID());
                      if (!chain->PDB_ChainID_Flag().empty())
                           t->UpdateCell(row, "pdbx_blank_PDB_chainid_flag", chain->PDB_ChainID_Flag());
                      else t->UpdateCell(row, "pdbx_blank_PDB_chainid_flag", "N");
                      if (chain->isModified())
                           t->UpdateCell(row, "pdbx_modified", "Y");
                      else t->UpdateCell(row, "pdbx_modified", "N");
                      t->UpdateCell(row, "entity_id", chain->entity_id());
                      t->UpdateCell(row, "details", chain->details());
                      if (!public_flag) {
                           t->UpdateCell(row, "pdbx_PDB_id", chain->PDB_ChainID());
                           t->UpdateCell(row, "pdbx_alt_id", chain->PUB_ChainID());
                           t->UpdateCell(row, "pdbx_type", chain->chain_type());
                           t->UpdateCell(row, "pdbx_order", String::IntToString(row + 1));
                      }
                      row++;
                      chain = single_mol->GetNextChain();
                 }
                 block.WriteTable(t);

                 delete single_mol;

                 if (_has_symmetry_related_operation) {
                      t = getTableCopy(_cifblock, "pdbx_struct_oper_list");
                      if (t) block.WriteTable(t);
                 }
            } else {
                 _get_entity_poly_values(_cifblock, chain_id_set, data, values);
                 _update_values(block, "entity_poly", data, values);

                 _get_entity_values(_cifblock, "entity", "id", entity_id_set, data, values);
                 _update_values(block, "entity", data, values);

                 _get_entity_values(_cifblock, "entity_poly_seq", "entity_id", entity_id_set, data, values);
                 _update_values(block, "entity_poly_seq", data, values);
            }

            if (t1->GetNumRows() > 0)
                 block.WriteTable(t1);
            else {
                 delete t1;
                 deleteTable(block, "atom_site");
            }

            if (t2->GetNumRows() > 0)
                 block.WriteTable(t2);
            else {
                 delete t2;
                 deleteTable(block, "atom_site_anisotrop");
            }

            assembly_id++;
            std::string filename = rootname + "_assembly-model_P" + String::IntToString(assembly_id) + ".cif.V1";
            if (public_flag) filename = rootname + "-assembly" + String::IntToString(assembly_id) + ".cif";
            fprintf(fp, "%s\n", filename.c_str());

            get_order_list(order_file, block, order_list);

            fobj->SetQuoting(CifFile::eDOUBLE);
            if (!order_list.empty())
                 fobj->Write(filename, order_list);
            else fobj->Write(filename);
       }

       fclose (fp);
       delete fobj;
}

void GenBioAssembly::_get_original_chains(const std::vector<std::string>& pdb_chains, std::vector<RCSB::Chain*>& original_chains,
                                          std::set<std::string>& chain_id_set, std::set<std::string>& entity_id_set)
{
       original_chains.clear();

       if (_molecules.empty()) return;

       std::string asym_ids = _get_asym_id_list_from_chain_ids(pdb_chains);
       if (asym_ids.empty()) return;

       std::vector<std::string> asym_list;
       get_wordarray(asym_list, asym_ids, ",");

       for (std::vector<std::string>::const_iterator pos = asym_list.begin(); pos != asym_list.end(); ++pos) {
            RCSB::Chain* chain = _molecules[0]->GetAsymChain(*pos);
            if (chain) {
                 if (chain->chain_type() == "ATOMN" || chain->chain_type() == "ATOMP") chain_id_set.insert(chain->PDB_ChainID());
                 entity_id_set.insert(chain->entity_id());
                 original_chains.push_back(chain);
            }
       }
}

RCSB::Molecule* GenBioAssembly::_get_bio_molecule(const std::vector<RCSB::Chain*>& original_chains, double mat[4][4])
{
       std::vector<std::set<int> > chain_index_list;
       chain_index_list.clear();
       std::map<int, std::vector<std::string> > original_asym_chain_id_mapping;
       original_asym_chain_id_mapping.clear();

       RCSB::Molecule* mol = new RCSB::Molecule;
       mol->setCCDic(_ccDic);
       mol->setMessage(&_error_messages);
       mol->setLog(_logIo);
       _insert_assembly_chains(false, original_chains, "", mat, mol, chain_index_list, original_asym_chain_id_mapping);
       return mol;
}

void GenBioAssembly::_insert_assembly_chains(const bool& single_model_flag, const std::vector<RCSB::Chain*>& original_chains, const std::string& operators,
                                             double mat[4][4], RCSB::Molecule* mol, std::vector<std::set<int> >& chain_index_list,
                                             std::map<int, std::vector<std::string> >& original_asym_chain_id_mapping)
{
       bool is_identical_operator = true;
       std::string oper_suffix = "";
       if (single_model_flag && !operators.empty()) {
            for (int i = 0; i < 4; ++i) {
                 for (int j = 0; j < 4; ++j) {
                      if (i == j) {
                           if (fabs(mat[i][j] - 1.0) > 0.00001) is_identical_operator = false;
                      } else if (fabs(mat[i][j]) > 0.00001) is_identical_operator = false;
                 }
            }
            if (!operators.empty()) {
                 oper_suffix = operators;
                 replace_string(oper_suffix, " ", "-");
                 
            }
            if (!is_identical_operator) _has_symmetry_related_operation = true;
       }

       std::vector<RCSB::Residue*> residues;
       std::vector<std::string> data;

       std::map<std::string, std::string> chain_id_mapping;
       chain_id_mapping.clear();

       std::set<int> chain_index_set;
       chain_index_set.clear();

       for (std::vector<RCSB::Chain*>::const_iterator cpos = original_chains.begin(); cpos != original_chains.end(); ++cpos) {
            (*cpos)->GetFirstResidueList(residues);
            RCSB::Chain* chain = new RCSB::Chain;
            chain->setCCDic(_ccDic);
            chain->setMessage(&_error_messages);
            chain->setLog(_logIo);

            std::string asym_id = (*cpos)->ChainID();
            std::string chain_id = (*cpos)->PDB_ChainID();
            if (!oper_suffix.empty() && !is_identical_operator) {
                 asym_id += "-" + oper_suffix;
                 chain_id += "-" + oper_suffix;
            }

            chain->set_ChainID(asym_id);
            chain->set_PDB_ChainID(chain_id);

            chain->set_chain_type((*cpos)->chain_type());
            chain->set_entity_id((*cpos)->entity_id());
            chain->set_entity_key((*cpos)->entity_key());

            while (!residues.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                      RCSB::Residue* res = new RCSB::Residue;
                      res->setCCDic(_ccDic);
                      res->setLog(_logIo);
                      res->setMessage(&_error_messages);
                      res->set_token((*rpos)->token());
                      res->set_ResName((*rpos)->ResName());
                      res->set_chnid(asym_id);
                      res->set_res_no((*rpos)->res_no());
                      res->set_pdb_chnid(chain_id);
                      res->set_pdb_res_no((*rpos)->pdb_res_no());
                      res->set_ins_code((*rpos)->ins_code());
                      res->set_chain_type((*rpos)->chain_type());

                      RCSB::Atom* atm = (*rpos)->GetFirstAtom();
                      while (atm) {
                           RCSB::Atom* atom = new RCSB::Atom;
                           *atom = *atm;
                           atom->set_chnid(asym_id);
                           atom->set_pdb_chnid(chain_id);
                           res->insert_a_atom(atom, 1);
                           atm = (*rpos)->GetNextAtom();
                      }
                      res->transformation(mat);

                      chain->insert_a_residue(res);
                      mol->insert_a_residue(res);
                 }
                 (*cpos)->GetNextResidueList(residues);
            }
            mol->insert_a_chain(chain);
            chain_index_set.insert(chain->index());

            data.clear();
            data.push_back((*cpos)->ChainID());
            data.push_back((*cpos)->PDB_ChainID());
            data.push_back(oper_suffix);
            original_asym_chain_id_mapping.insert(std::make_pair(chain->index(), data));
       }
       chain_index_list.push_back(chain_index_set);
}

void GenBioAssembly::_update_atom_site_categories(const int& mol_id, RCSB::Molecule* mol, int& atomSerialNo, ISTable* atomTable, ISTable* anisotropTable)
{
       const std::vector<std::string>& atomItemNames = atomTable->GetColumnNames();
       const std::vector<std::string>& anisotropItemNames = anisotropTable->GetColumnNames();

       std::vector<RCSB::Residue*> residue_list;

       RCSB::Chain* chain = mol->GetFirstChain();
       while (chain) {
            chain->GetFirstResidueList(residue_list);
            while (!residue_list.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                      std::string PDB_atom_token = "HETATM";
                      if ((chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") && (SeqCodeUtil::is_a_standard_residue((*rpos)->ResName()) ||
                          (*rpos)->ResName() == "T")) PDB_atom_token = "ATOM";
                      RCSB::Atom* atom = (*rpos)->GetFirstAtom();
                      while (atom) {
                           atom->setValue(PDB_atom_token, 0);
                           atom->setValue(chain->entity_id(), 23);
                           atomSerialNo++;
                           CifCoordWrite::update_atom_site(atomItemNames, mol_id, atomSerialNo, atom, atomTable);
                           _update_atom_site_anisotrop(anisotropItemNames, atomSerialNo, atom, anisotropTable);
                           atom = (*rpos)->GetNextAtom();
                      }
                 }
                 chain->GetNextResidueList(residue_list);
            }
            chain = mol->GetNextChain();
       }
}

void GenBioAssembly::_get_entity_poly_values(Block& block, const std::set<std::string>& chain_id_set, std::vector<std::string>&
                                             items, std::vector<std::map<std::string, std::string> >& values)
{
       items.clear();
       values.clear();

       if (chain_id_set.empty()) return;

       ISTable *t = getTablePtr(block, "entity_poly");
       if (!t) return;

       std::vector<std::map<std::string, std::string> > all_values;
       get_values(t, all_values);
       if (all_values.empty()) return;

       items = t->GetColumnNames();

       std::vector<std::string> data;
       std::string pdbx_strand_id;
       for (std::vector<std::map<std::string, std::string> >::iterator pos = all_values.begin(); pos != all_values.end(); ++pos) {
            get_wordarray(data, (*pos)["pdbx_strand_id"], ",;");
            pdbx_strand_id.clear();
            for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                 if (chain_id_set.find(*vpos) == chain_id_set.end()) continue;
                 if (!pdbx_strand_id.empty()) pdbx_strand_id += ",";
                 pdbx_strand_id += *vpos;
            }
            if (pdbx_strand_id.empty()) continue;
            (*pos)["pdbx_strand_id"] = pdbx_strand_id;
            values.push_back(*pos);
       }
}

void GenBioAssembly::_get_entity_values(Block& block, const std::string& category, const std::string& key_item, const std::set<std::string>& entity_id_set,
                                        std::vector<std::string>& items, std::vector<std::map<std::string, std::string> >& values)
{
       items.clear();
       values.clear();

       ISTable *t = getTablePtr(block, category);
       if (!t) return;

       std::vector<std::map<std::string, std::string> > all_values;
       get_values(t, all_values);
       if (all_values.empty()) return;

       const std::vector<std::string>& item_names = t->GetColumnNames();
       for (std::vector<std::string>::const_iterator vpos = item_names.begin(); vpos != item_names.end(); ++vpos) {
            if (*vpos == "pdbx_number_of_molecules") continue;
            items.push_back(*vpos);
       }

       for (std::vector<std::map<std::string, std::string> >::iterator pos = all_values.begin(); pos != all_values.end(); ++pos) {
            if (entity_id_set.find((*pos)[key_item]) == entity_id_set.end()) continue;
            values.push_back(*pos);
       }
}

void GenBioAssembly::_update_values(Block& block, const std::string& category, const std::vector<std::string>& items,
                                    const std::vector<std::map<std::string, std::string> >& values)
{
       if (values.empty()) {
            deleteTable(block, category);
            return;
       }

       ISTable *t = add_new_table(category, items);
       int row = 0;
       for (std::vector<std::map<std::string, std::string> >::const_iterator pos = values.begin(); pos != values.end(); ++pos) {
            t->AddRow();
            for (std::map<std::string, std::string>::const_iterator mpos = pos->begin(); mpos != pos->end(); ++mpos) {
                 if (t->IsColumnPresent(mpos->first)) t->UpdateCell(row, mpos->first, mpos->second);
            }
            row++;
       }
       if (_remove_non_public_item_flag) {
            for (std::vector<std::string>::const_iterator cpos = items.begin(); cpos != items.end(); ++cpos) {
                 if (_dictUtil.isPublicItem(category, *cpos)) continue;
                 t->DeleteColumn(*cpos);
            }
       }
       block.WriteTable(t);
}
