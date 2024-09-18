/*
FILE:     mergeSiteUtil.C
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

#include "AnnotationObj.h"
#include "utillib.h"

static std::string generate_next_site_id(const std::set<std::string> &used_id);

void AnnotationObj::merge_struct_sites(const std::string& sitefile)
{
       if (_skip_task_set.find("site") != _skip_task_set.end()) return;

       CifFile *fobj = get_fobj(_logIo, "", sitefile);
       if (!fobj) return;

       std::list<_SITE> sites;

       Block &block = fobj->GetBlock(fobj->GetFirstBlockName());
       _read_struct_sites(block, sites);

       delete fobj;

       if (sites.empty()) return;

       // removed existing software generating site record
       std::set<std::string> used_ids;
       used_ids.clear();
       for (std::list<_SITE>::iterator lpos = _site.begin(); lpos != _site.end(); ++lpos) {
            if (lpos->code.empty() || String::IsEqual(lpos->code, "SOFTWARE", Char::eCASE_INSENSITIVE) ||
                String::IsEqual(lpos->code, "SOFTWARE SUPPLIED", Char::eCASE_INSENSITIVE)) {
                 lpos = _site.erase(lpos);
                 --lpos;
            } else used_ids.insert(lpos->SiteID);
       }

       for (std::list<_SITE>::iterator lpos = sites.begin(); lpos != sites.end(); ++lpos) {
            std::string site_id = generate_next_site_id(used_ids);
            lpos->SiteID = site_id;
            used_ids.insert(site_id);
            _site.push_back(*lpos);
       }
}

void AnnotationObj::_read_struct_sites(Block& block, std::list<_SITE>& sites)
{
       sites.clear();

       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "struct_site_gen");
       if (!t) return;

       std::map<int, std::string> site_order;
       site_order.clear();

       std::map<std::string, _SITE> site_mapping;
       site_mapping.clear();

       _SITE site;
       site.mol_index = 0;
       site.SiteID.clear();
       site.code.clear();
       site.details.clear();
       site.self_chain_indices.clear();
       site.self_res_indices.clear();
       site.associated_residues.clear();

       std::string pdb_chnid, res_name, res_num, ins_code, site_id, symmetry, cs;
       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            get_value_clean(pdb_chnid, t, i, "auth_asym_id");
            get_value_clean(res_name, t, i, "auth_comp_id");
            get_value_clean(res_num, t, i, "auth_seq_id");
            get_value_clean(ins_code, t, i, "pdbx_auth_ins_code");
            RCSB::Residue* res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
            if (!res) continue;

            get_value_clean(site_id, t, i, "site_id");
            get_value_clean(symmetry, t, i, "symmetry");
            if (symmetry == "1_555") symmetry.clear();
            std::map<std::string, _SITE>::iterator
                mpos = site_mapping.find(site_id);
            if (mpos != site_mapping.end())
                 mpos->second.associated_residues.push_back(std::make_pair(res->index(), symmetry));
            else {
                 site.SiteID = site_id;
                 site.associated_residues.clear();
                 site.associated_residues.push_back(std::make_pair(res->index(), symmetry));
                 site_order.insert(std::make_pair(site_mapping.size(), site_id));
                 site_mapping.insert(std::make_pair(site_id, site));
            }
       }

       if (site_mapping.empty()) return;

       std::map<std::string, std::vector<std::string> > index_mapping;
       index_mapping.clear();

       std::vector<std::string> data;

       t = getTablePtr(block, "pdbx_struct_site_contents");
       if (t && t->IsColumnPresent("site_id") && t->IsColumnPresent("auth_asym_id") && t->IsColumnPresent("auth_comp_id") &&
           t->IsColumnPresent("auth_seq_id") && t->IsColumnPresent("auth_ins_code")) {
            std::vector<std::string> items;
            items.clear();
            items.push_back("auth_asym_id");
            items.push_back("auth_comp_id");
            items.push_back("auth_seq_id");
            items.push_back("auth_ins_code");
            rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 data.clear();
                 for (std::vector<std::string>::const_iterator vpos = items.begin(); vpos != items.end(); ++vpos) {
                      get_value_clean(cs, t, i, *vpos);
                      data.push_back(cs);
                 }
                 get_value_clean(cs, t, i, "site_id");
                 index_mapping.insert(std::make_pair(cs, data));
            }
       }

       t = getTablePtr(block, "struct_site");
       if (t) {
            rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean(site_id, t, i, "id");
                 std::map<std::string, _SITE>::iterator mpos = site_mapping.find(site_id);
                 if (mpos == site_mapping.end()) continue;

                 get_value_clean_upper(cs, t, i, "pdbx_evidence_code");
                 mpos->second.code = cs;
                 get_value_clean(cs, t, i, "details");
                 mpos->second.details = cs;

                 pdb_chnid.clear();
                 res_name.clear();
                 res_num.clear();
                 ins_code.clear();
                 std::map<std::string, std::vector<std::string> >::const_iterator ipos = index_mapping.find(site_id);
                 if (ipos != index_mapping.end()) {
                      pdb_chnid = ipos->second[0];
                      res_name  = ipos->second[1];
                      res_num   = ipos->second[2];
                      ins_code  = ipos->second[3];
                 } else {
                      get_value_clean(pdb_chnid, t, i, "pdbx_auth_asym_id");
                      get_value_clean(res_name, t, i, "pdbx_auth_comp_id");
                      get_value_clean(res_num, t, i, "pdbx_auth_seq_id");
                      get_value_clean(ins_code, t, i, "pdbx_auth_ins_code");
                      if (pdb_chnid.empty() && res_name.empty() && res_num.empty() && !mpos->second.details.empty()) {
                           get_wordarray(data, mpos->second.details, " ");
                           if ((data.size() == 7) && String::IsEqual(data[0], "binding", Char::eCASE_INSENSITIVE) &&
                                String::IsEqual(data[1], "site", Char::eCASE_INSENSITIVE) && String::IsEqual(data[2], "for", Char::eCASE_INSENSITIVE) &&
                                String::IsEqual(data[3], "residue", Char::eCASE_INSENSITIVE)) {
                                pdb_chnid = data[5];
                                res_name = data[4];
                                if (isalpha(data[6][data[6].size() - 1])) {
                                     res_num = data[6].substr(0, data[6].size() - 1);
                                     ins_code = data[6].substr(data[6].size() - 1);
                                } else {
                                     res_num = data[6];
                                     ins_code.clear();
                                }
                           }
                      }
                 }
                 if (!pdb_chnid.empty() && !res_name.empty() && !res_num.empty()) {
                      RCSB::Residue* res = _molecules[0]->find_pdb_residue(pdb_chnid, res_name, res_num, ins_code);
                      if (res) mpos->second.self_res_indices.push_back(res->index());
                 } else if (!pdb_chnid.empty()) {
                      RCSB::Chain* chain = _molecules[0]->GetPolyChain(pdb_chnid);
                      if (chain) mpos->second.self_chain_indices.push_back(chain->index());
                 }
            }
       }

       for (std::map<int, std::string>::const_iterator pos = site_order.begin(); pos != site_order.end(); ++pos) {
            std::map<std::string, _SITE>::const_iterator mpos = site_mapping.find(pos->second);
            if (mpos == site_mapping.end()) continue;

            sites.push_back(mpos->second);
       }
}

void AnnotationObj::_write_struct_sites(Block& block)
{
       if (_site.empty() || _molecules.empty()) {
            deleteTable(block, "struct_site");
            deleteTable(block, "struct_site_gen");
            return;
       }

       std::vector<std::string> data;

       ISTable *t1 = _newTablePtr("struct_site");
       ISTable *t2 = _newTablePtr("struct_site_gen");
       int row1 = 0;
       int row2 = 0;
       for (std::list<_SITE>::const_iterator lpos = _site.begin(); lpos != _site.end(); ++lpos) {
            int num_residues = 0;
            for (std::vector<std::pair<int, std::string> >::const_iterator ppos =
                 lpos->associated_residues.begin(); ppos != lpos->associated_residues.end(); ++ppos) {
                 bool is_removed = false;
                 RCSB::Residue* res = _molecules[0]->find_residue(ppos->first, is_removed);
                 if (!res) continue;
                 num_residues++;
            }
            if (!num_residues) continue;

            if (!lpos->self_chain_indices.empty()) {
                 bool is_removed = false;
                 RCSB::Chain* chain = _molecules[0]->GetIndexChain(lpos->self_chain_indices[0], is_removed);
                 if (!chain) continue;
            }
            if (!lpos->self_res_indices.empty()) {
                 bool is_removed = false;
                 RCSB::Residue* res = _molecules[0]->find_residue(lpos->self_res_indices[0], is_removed);
                 if (!res) continue;
            }

            t1->AddRow();
            t1->UpdateCell(row1, "id", lpos->SiteID);
            if (lpos->code == "SOFTWARE")
                 t1->UpdateCell(row1, "pdbx_evidence_code", "Software");
            else if (lpos->code == "AUTHOR")
                 t1->UpdateCell(row1, "pdbx_evidence_code", "Author");
            else if (lpos->code == "UNKNOWN")
                 t1->UpdateCell(row1, "pdbx_evidence_code", "Unknown");
            else t1->UpdateCell(row1, "pdbx_evidence_code", lpos->code);
            t1->UpdateCell(row1, "pdbx_num_residues", String::IntToString(num_residues));
            t1->UpdateCell(row1, "details", lpos->details);
            if (!lpos->self_chain_indices.empty()) {
                 bool is_removed = false;
                 RCSB::Chain* chain = _molecules[0]->GetIndexChain(lpos->self_chain_indices[0], is_removed);
                 if (chain) t1->UpdateCell(row1, "pdbx_auth_asym_id", chain->PDB_ChainID());
            } else if (!lpos->self_res_indices.empty()) {
                 bool  is_removed = false;
                 RCSB::Residue* res = _molecules[0]->find_residue(lpos->self_res_indices[0], is_removed);
                 if (res) {
                      get_wordarray(data, lpos->details, " ");
                      if (data.size() > 4 && String::IsEqual(data[0], "binding", Char::eCASE_INSENSITIVE) &&
                          String::IsEqual(data[1], "site", Char::eCASE_INSENSITIVE) && String::IsEqual(data[2], "for", Char::eCASE_INSENSITIVE) &&
                          String::IsEqual(data[3], "residue", Char::eCASE_INSENSITIVE)) {
                           std::string cs = data[0] + " " + data[1] + " " + data[2] + " " + data[3] + " " + res->ResName() + " " + res->pdb_chnid() + " "
                                          + res->pdb_res_no() + res->ins_code();
                           t1->UpdateCell(row1, "details", cs);
                      }
                      t1->UpdateCell(row1, "pdbx_auth_asym_id", res->pdb_chnid());
                      t1->UpdateCell(row1, "pdbx_auth_comp_id", res->ResName());
                      t1->UpdateCell(row1, "pdbx_auth_seq_id", res->pdb_res_no());
                      t1->UpdateCell(row1, "pdbx_auth_ins_code", res->ins_code());
                 }
            }
            row1++;

            for (std::vector<std::pair<int, std::string> >::const_iterator ppos =
                 lpos->associated_residues.begin(); ppos != lpos->associated_residues.end(); ++ppos) {
                 bool is_removed = false;
                 RCSB::Residue* res = _molecules[0]->find_residue(ppos->first, is_removed);
                 if (!res) continue;

                 t2->AddRow();
                 t2->UpdateCell(row2, "id", String::IntToString(row2 + 1));
                 t2->UpdateCell(row2, "site_id", lpos->SiteID);
                 t2->UpdateCell(row2, "pdbx_num_res", String::IntToString(num_residues)); 
                 t2->UpdateCell(row2, "label_comp_id", res->ResName());
                 t2->UpdateCell(row2, "label_asym_id", res->chnid());
                 if (!res->res_no().empty())
                      t2->UpdateCell(row2, "label_seq_id", res->res_no());
                 else t2->UpdateCell(row2, "label_seq_id", ".");
                 t2->UpdateCell(row2, "pdbx_auth_ins_code", res->ins_code());
                 t2->UpdateCell(row2, "auth_comp_id", res->ResName());
                 t2->UpdateCell(row2, "auth_asym_id", res->pdb_chnid());
                 t2->UpdateCell(row2, "auth_seq_id", res->pdb_res_no());
                 t2->UpdateCell(row2, "label_atom_id", ".");
                 if (!ppos->second.empty())
                      t2->UpdateCell(row2, "symmetry", ppos->second);
                 else t2->UpdateCell(row2, "symmetry", "1_555");
                 row2++;
            }
       }

       if (row1 == 0 || row2 == 0) {
            delete t1;
            delete t2;
            deleteTable(block, "struct_site");
            deleteTable(block, "struct_site_gen");
       } else {
            block.WriteTable(t1);
            block.WriteTable(t2);
       }
}

void AnnotationObj::_remove_software_generated_site_record()
{
       if (_site.empty()) return;

       for (std::list<_SITE>::iterator lpos = _site.begin(); lpos != _site.end(); ++lpos) {
            if (String::IsEqual(lpos->code, "SOFTWARE", Char::eCASE_INSENSITIVE) || String::IsEqual(lpos->code, "SOFTWARE SUPPLIED", Char::eCASE_INSENSITIVE)) {
                 lpos = _site.erase(lpos);
                 --lpos;
            }
       }
}

static std::string generate_next_site_id(const std::set<std::string> &used_id)
{
       const std::string firstcolumn = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
       const std::string secondcolumn = "CDEFGHIJKLMNOPQRSTUVWXYZ";
       const std::string thirdcolumn = "123456789";

       std::string id;
       if (used_id.size() < 5616) {
            for (unsigned int i = 0; i < firstcolumn.size(); i++) {
                 for (unsigned int j = 0; j < secondcolumn.size(); j++) {
                      for (unsigned int k = 0; k < thirdcolumn.size(); k++) {
                           id.clear();
                           id += firstcolumn[i];
                           id += secondcolumn[j];
                           id += thirdcolumn[k];
                           if (used_id.find(id) == used_id.end()) return id;
                      }
                 }
            }
       } if (used_id.size() < 23192) {
            std::string sndcolumn = firstcolumn;
            std::string trdcolumn = firstcolumn;
            for (unsigned int i = 0; i < firstcolumn.size(); i++) {
                 for (unsigned int j = 0; j < sndcolumn.size(); j++) {
                      for (unsigned int k = 0; k < trdcolumn.size(); k++) {
                           id.clear();
                           id += firstcolumn[i];
                           id += sndcolumn[j];
                           id += trdcolumn[k];
                           if (used_id.find(id) == used_id.end()) return id;
                      }
                 }
            }
       } else {
            std::string sndcolumn = firstcolumn;
            std::string trdcolumn = firstcolumn;
            for (unsigned int i = 0; i < firstcolumn.size(); i++) {
                 for (unsigned int j = 0; j < sndcolumn.size(); j++) {
                      for (unsigned int k = 0; k < trdcolumn.size(); k++) {
                           for (unsigned int l = 0; l < thirdcolumn.size(); ++l) {
                                id.clear();
                                id += firstcolumn[i];
                                id += sndcolumn[j];
                                id += trdcolumn[k];
                                id += thirdcolumn[l];
                                if (used_id.find(id) == used_id.end()) return id;
                           }
                      }
                 }
            }
       }

       return "ZZZ0";
}
