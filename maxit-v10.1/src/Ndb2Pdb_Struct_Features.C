/*
FILE:     Ndb2Pdb_Struct_Features.C
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
#include "Maxit.h"
#include "NdbToken.h"
#include "SeqCodeUtil.h"
#include "utillib.h"

void Maxit::_ndb_to_pdb_update_HELIX()
{
       _pdb_records.erase("HELIX");

       if (_helix.empty() || _molecules.empty()) return;

       std::vector<std::string> FieldInfo;
       FieldInfo.clear();
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("HELIX");
       for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
       FieldInfo[0] = "HELIX";

       std::list<std::vector<std::string> > records;
       records.clear();

       for (std::list<_HELIX>::const_iterator lpos = _helix.begin(); lpos != _helix.end(); ++lpos) {
            bool is_removed = false;
            RCSB::Residue* init = _molecules[0]->find_residue(lpos->initRes, is_removed);
            if (!init) continue;
            RCSB::Residue* end = _molecules[0]->find_residue(lpos->endRes, is_removed);
            if (!end) continue;
            for (unsigned int i = 1; i < FieldInfo.size(); ++i) {
                 FieldInfo[i].clear();
            }
            FieldInfo[1] = String::IntToString(records.size() + 1);
            FieldInfo[2] = lpos->ID;
            FieldInfo[3] = init->ResName();
            FieldInfo[4] = init->pdb_chnid();
            FieldInfo[5] = init->pdb_res_no();
            FieldInfo[6] = init->ins_code();
            FieldInfo[7] = end->ResName();
            FieldInfo[8] = end->pdb_chnid();
            FieldInfo[9] = end->pdb_res_no();
            FieldInfo[10] = end->ins_code();
            FieldInfo[11] = String::IntToString(lpos->helixClass);
            FieldInfo[12] = lpos->comment;
            FieldInfo[13] = String::IntToString(end->position() - init->position() + 1);
            records.push_back(FieldInfo);
       }
       if (!records.empty()) _pdb_records.insert(std::make_pair("HELIX", records));
}

void Maxit::_ndb_to_pdb_update_SHEET()
{
       _pdb_records.erase("SHEET");

       if (_sheet.empty() || _molecules.empty()) return;

       for (std::list<_SHEET>::const_iterator lpos = _sheet.begin(); lpos != _sheet.end(); ++lpos) {
            if (lpos->complicateFlag) return;
       }

       std::vector<std::string> FieldInfo;
       FieldInfo.clear();
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("SHEET");
       for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
       FieldInfo[0] = "SHEET";

       std::list<std::vector<std::string> > all_records, sheet_records;
       all_records.clear();

       bool is_removed = false;
       for (std::list<_SHEET>::const_iterator lpos = _sheet.begin(); lpos != _sheet.end(); ++lpos) {
            sheet_records.clear();
            for (std::vector<_SHEET_STRAND>::const_iterator sspos = lpos->_strands.begin(); sspos != lpos->_strands.end(); ++sspos) {
                 RCSB::Residue* beginRes = _molecules[0]->find_residue(sspos->begin_res_index, is_removed);
                 if (!beginRes) {
                      sheet_records.clear();
                      break;
                 }
                 RCSB::Residue* endRes = _molecules[0]->find_residue(sspos->end_res_index, is_removed);
                 if (!endRes) {
                      sheet_records.clear();
                      break;
                 }
                 RCSB::Residue* currRes = NULL;
                 if (sspos->curr_hbond_res_index >= 0) {
                      currRes = _molecules[0]->find_residue(sspos->curr_hbond_res_index, is_removed);
                      if (!currRes) {
                           sheet_records.clear();
                           break;
                      }
                 }
                 RCSB::Residue* prevRes = NULL;
                 if (sspos->prev_hbond_res_index >= 0) {
                      prevRes = _molecules[0]->find_residue(sspos->prev_hbond_res_index, is_removed);
                      if (!prevRes) {
                           sheet_records.clear();
                           break;
                      }
                 }

                 for (unsigned int i = 1; i < FieldInfo.size(); ++i) {
                      FieldInfo[i].clear();
                 }
                 FieldInfo[1] = String::IntToString(sspos->strand_id);
                 FieldInfo[2] = lpos->sheetID;
                 FieldInfo[3] = String::IntToString(lpos->numStrands);
                 FieldInfo[4] = beginRes->ResName();
                 FieldInfo[5] = beginRes->pdb_chnid();
                 FieldInfo[6] = beginRes->pdb_res_no();
                 FieldInfo[7] = beginRes->ins_code();
                 FieldInfo[8] = endRes->ResName();
                 FieldInfo[9] = endRes->pdb_chnid();
                 FieldInfo[10] = endRes->pdb_res_no();
                 FieldInfo[11] = endRes->ins_code();
                 FieldInfo[12] = String::IntToString(sspos->sense);
                 if (currRes && prevRes) {
                      FieldInfo[13] = sspos->curr_hbond_atom_name;
                      FieldInfo[14] = currRes->ResName();
                      FieldInfo[15] = currRes->pdb_chnid();
                      FieldInfo[16] = currRes->pdb_res_no();
                      FieldInfo[17] = currRes->ins_code();
                      FieldInfo[18] = sspos->prev_hbond_atom_name;
                      FieldInfo[19] = prevRes->ResName();
                      FieldInfo[20] = prevRes->pdb_chnid();
                      FieldInfo[21] = prevRes->pdb_res_no();
                      FieldInfo[22] = prevRes->ins_code();
                 }
                 sheet_records.push_back(FieldInfo);
            }
            for (std::list<std::vector<std::string> >::const_iterator lrpos = sheet_records.begin(); lrpos != sheet_records.end(); ++lrpos) {
                 all_records.push_back(*lrpos);
            }
       }
       if (!all_records.empty()) _pdb_records.insert(std::make_pair("SHEET", all_records));
}

void Maxit::_ndb_to_pdb_update_TURN()
{
       _pdb_records.erase("TURN");

       if (_turn.empty() || _molecules.empty()) return;

       std::vector<std::string> FieldInfo;
       FieldInfo.clear();
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("TURN");
       for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
       FieldInfo[0] = "TURN";

       std::list<std::vector<std::string> > records;
       records.clear();

       for (std::list<_HELIX>::const_iterator lpos = _turn.begin(); lpos != _turn.end(); ++lpos) {
            bool is_removed = false;
            RCSB::Residue* init = _molecules[0]->find_residue(lpos->initRes, is_removed);
            if (!init) continue;
            RCSB::Residue* end = _molecules[0]->find_residue(lpos->endRes, is_removed);
            if (!end) continue;
            for (unsigned int i = 1; i < FieldInfo.size(); ++i) {
                 FieldInfo[i].clear();
            }
            FieldInfo[1] = String::IntToString(records.size() + 1);
            FieldInfo[2] = lpos->ID;
            FieldInfo[3] = init->ResName();
            FieldInfo[4] = init->pdb_chnid();
            FieldInfo[5] = init->pdb_res_no();
            FieldInfo[6] = init->ins_code();
            FieldInfo[7] = end->ResName();
            FieldInfo[8] = end->pdb_chnid();
            FieldInfo[9] = end->pdb_res_no();
            FieldInfo[10] = end->ins_code();
            FieldInfo[11] = lpos->comment;
            records.push_back(FieldInfo);
       }
       if (!records.empty()) _pdb_records.insert(std::make_pair("TURN", records));
}

void Maxit::_ndb_to_pdb_update_SSBOND()
{
       _pdb_records.erase("SSBOND");

       if (_ssbonds.empty()) return;

       std::vector<std::string> FieldInfo;
       FieldInfo.clear();
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("SSBOND");
       for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
       FieldInfo[0] = "SSBOND";

       std::list<std::vector<std::string> > records;
       records.clear();

       std::set<std::string> index_set;
       index_set.clear();
       for (std::list<_SSBOND>::const_iterator lpos = _ssbonds.begin(); lpos != _ssbonds.end(); ++lpos) {
            if (lpos->mol_index != 0) continue;

            std::string index = lpos->fstAtom->pdb_chnid() + "_" + lpos->fstAtom->pdb_resnam() + "_" + lpos->fstAtom->pdb_resnum() + "_"
                              + lpos->fstAtom->ins_code() + "_" + lpos->sndAtom->pdb_chnid() + "_" + lpos->sndAtom->pdb_resnam() + "_"
                              + lpos->sndAtom->pdb_resnum() + "_" + lpos->sndAtom->ins_code();
            if (index_set.find(index) != index_set.end()) continue;
            index_set.insert(index);
            index = lpos->sndAtom->pdb_chnid()  + "_" + lpos->sndAtom->pdb_resnam() + "_" + lpos->sndAtom->pdb_resnum() + "_" 
                  + lpos->sndAtom->ins_code() + "_" + lpos->fstAtom->pdb_chnid()  + "_" + lpos->fstAtom->pdb_resnam() + "_"
                  + lpos->fstAtom->pdb_resnum() + "_" + lpos->fstAtom->ins_code();
            index_set.insert(index);
      
            FieldInfo[1] = String::IntToString(records.size() + 1);
            FieldInfo[2] = lpos->fstAtom->pdb_resnam();
            FieldInfo[3] = lpos->fstAtom->pdb_chnid();
            FieldInfo[4] = lpos->fstAtom->pdb_resnum();
            FieldInfo[5] = lpos->fstAtom->ins_code();
            FieldInfo[6] = lpos->sndAtom->pdb_resnam();
            FieldInfo[7] = lpos->sndAtom->pdb_chnid();
            FieldInfo[8] = lpos->sndAtom->pdb_resnum();
            FieldInfo[9] = lpos->sndAtom->ins_code();
            FieldInfo[10] = getPDBSymmetry(lpos->SymOP_1);
            FieldInfo[11] = getPDBSymmetry(lpos->SymOP_2);
            FieldInfo[12] = FormattedFieldValue(lpos->dist, 2, 5, 2, false, true);
            records.push_back(FieldInfo);
       }
       if (!records.empty()) _pdb_records.insert(std::make_pair("SSBOND", records));
}

void Maxit::_ndb_to_pdb_update_LINK()
{
       _pdb_records.erase("LINK");

       if (_links.empty()) return;

       std::vector<std::string> FieldInfo;
       FieldInfo.clear();
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("LINK");
       for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
       FieldInfo[0] = "LINK";

       std::list<std::vector<std::string> > records;
       records.clear();

       for (std::list<_LINK>::const_iterator lpos = _links.begin(); lpos != _links.end(); ++lpos) {
            if (lpos->mol_index != 0) continue;

            FieldInfo[1] = lpos->fstAtom->pdb_atmnam();
            FieldInfo[2] = lpos->fstAtom->alt_loc();
            FieldInfo[3] = lpos->fstAtom->pdb_resnam();
            FieldInfo[4] = lpos->fstAtom->pdb_chnid();
            FieldInfo[5] = lpos->fstAtom->pdb_resnum();
            FieldInfo[6] = lpos->fstAtom->ins_code();
            FieldInfo[7] = lpos->sndAtom->pdb_atmnam();
            FieldInfo[8] = lpos->sndAtom->alt_loc();
            FieldInfo[9] = lpos->sndAtom->pdb_resnam();
            FieldInfo[10] = lpos->sndAtom->pdb_chnid();
            FieldInfo[11] = lpos->sndAtom->pdb_resnum();
            FieldInfo[12] = lpos->sndAtom->ins_code();
            FieldInfo[13] = getPDBSymmetry(lpos->SymOP_1);
            FieldInfo[14] = getPDBSymmetry(lpos->SymOP_2);
            FieldInfo[15] = FormattedFieldValue(lpos->dist, 2, 5, 2, false, true);
            records.push_back(FieldInfo);
       }
       if (!records.empty()) _pdb_records.insert(std::make_pair("LINK", records));
}

void Maxit::_ndb_to_pdb_update_SLTBRG()
{
       _pdb_records.erase("SLTBRG");

       if (_sltbrgs.empty()) return;

       std::vector<std::string> FieldInfo;
       FieldInfo.clear();
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("SLTBRG");
       for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
       FieldInfo[0] = "SLTBRG";

       std::list<std::vector<std::string> > records;
       records.clear();

       for (std::list<_SLTBRG>::const_iterator lpos = _sltbrgs.begin(); lpos != _sltbrgs.end(); ++lpos) {
            if (lpos->mol_index != 0) continue;

            FieldInfo[1] = lpos->fstAtom->pdb_atmnam();
            FieldInfo[2] = lpos->fstAtom->alt_loc();
            FieldInfo[3] = lpos->fstAtom->pdb_resnam();
            FieldInfo[4] = lpos->fstAtom->pdb_chnid();
            FieldInfo[5] = lpos->fstAtom->pdb_resnum();
            FieldInfo[6] = lpos->fstAtom->ins_code();
            FieldInfo[7] = lpos->sndAtom->pdb_atmnam();
            FieldInfo[8] = lpos->sndAtom->alt_loc();
            FieldInfo[9] = lpos->sndAtom->pdb_resnam();
            FieldInfo[10] = lpos->sndAtom->pdb_chnid();
            FieldInfo[11] = lpos->sndAtom->pdb_resnum();
            FieldInfo[12] = lpos->sndAtom->ins_code();
            FieldInfo[13] = getPDBSymmetry(lpos->SymOP_1);
            FieldInfo[14] = getPDBSymmetry(lpos->SymOP_2);
            records.push_back(FieldInfo);
       }
       if (!records.empty()) _pdb_records.insert(std::make_pair("SLTBRG", records));
}

void Maxit::_ndb_to_pdb_update_MODRES()
{
       _pdb_records.erase("MODRES");

       if (_modres.empty() || _molecules.empty()) return;

       std::string pdb_id;
       _getRecordFront("PDBFIL", 1, pdb_id);

       std::vector<std::string> FieldInfo;
       FieldInfo.clear();
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("MODRES");
       for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
       FieldInfo[0] = "MODRES";
       FieldInfo[1] = pdb_id;

       std::list<std::vector<std::string> > records;
       records.clear();

       for (std::list<_MODRES>::const_iterator lpos = _modres.begin(); lpos != _modres.end(); ++lpos) {
            bool is_removed = false;
            RCSB::Residue*res = _molecules[0]->find_residue(lpos->res_index, is_removed);
            if (!res) continue;

            FieldInfo[2] = res->ResName();
            FieldInfo[3] = res->pdb_chnid();
            FieldInfo[4] = res->pdb_res_no();
            FieldInfo[5] = res->ins_code();
            FieldInfo[6] = lpos->Standard_Name;
            FieldInfo[7] = lpos->details;
            records.push_back(FieldInfo);
       }
       if (!records.empty()) _pdb_records.insert(std::make_pair("MODRES", records));
}

void Maxit::_ndb_to_pdb_update_SITE()
{
       _pdb_records.erase("SITELS");
       _pdb_records.erase("SITE");

       if (_site.empty() || _molecules.empty()) return;

       std::set<int> index_set;
       for (std::list<_SITE>::const_iterator lpos = _site.begin(); lpos != _site.end(); ++lpos) {
            int num_residues = 0;
            bool is_removed = false;
            index_set.clear();
            for (std::vector<std::pair<int, std::string> >::const_iterator ppos =
                 lpos->associated_residues.begin(); ppos != lpos->associated_residues.end(); ++ppos) {
                 if (index_set.find(ppos->first) != index_set.end()) continue;
                 index_set.insert(ppos->first);
                 RCSB::Residue*res = _molecules[0]->find_residue(ppos->first, is_removed);
                 if (!res) continue;
                 num_residues++;
            }
            if (!num_residues) continue;

            _addNewRecord("SITELS");
            _updateRecordBack("SITELS", 2, lpos->SiteID);
            _updateRecordBack("SITELS", 3, lpos->details);
            _updateRecordBack("SITELS", 4, lpos->code);

            int serial_no = 0;
            int siteId = 0;
            index_set.clear();
            for (std::vector<std::pair<int, std::string> >::const_iterator ppos =
                 lpos->associated_residues.begin(); ppos != lpos->associated_residues.end(); ++ppos) {
                 if (index_set.find(ppos->first) != index_set.end()) continue;
                 index_set.insert(ppos->first);
                 RCSB::Residue*res = _molecules[0]->find_residue(ppos->first, is_removed);
                 if (!res) continue;
                 if (siteId == 0) {
                      serial_no++;
                      _addNewRecord("SITE");
                      _updateRecordBack("SITE", 1, String::IntToString(serial_no));
                      _updateRecordBack("SITE", 2, lpos->SiteID);
                      _updateRecordBack("SITE", 3, String::IntToString(num_residues));
                 }
                 _updateRecordBack("SITE", siteId * 4 + 4, res->ResName());
                 _updateRecordBack("SITE", siteId * 4 + 5, res->pdb_chnid());
                 _updateRecordBack("SITE", siteId * 4 + 6, res->pdb_res_no());
                 _updateRecordBack("SITE", siteId * 4 + 7, res->ins_code());
                 siteId++;
                 if (siteId == 4) siteId = 0;
            }
       }
}

void Maxit::_ndb_to_pdb_get_atom_connects(const bool& checking_linkage_flag, std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& atom_connects)
{
       atom_connects.clear();

       if (_molecules.empty()) return;

       for (std::list<_SSBOND>::const_iterator lpos = _ssbonds.begin(); lpos != _ssbonds.end(); ++lpos) {
            if (lpos->mol_index != 0 || (!lpos->SymOP_1.empty() && lpos->SymOP_1 != "1_555") ||
               (!lpos->SymOP_2.empty() && lpos->SymOP_2 != "1_555")) continue;
            atom_connects.push_back(std::make_pair(lpos->fstAtom, lpos->sndAtom));
       }

       for (std::list<_LINK>::const_iterator lpos = _links.begin(); lpos != _links.end(); ++lpos) {
            if (lpos->mol_index != 0 || (!lpos->SymOP_1.empty() && lpos->SymOP_1 != "1_555") ||
               (!lpos->SymOP_2.empty() && lpos->SymOP_2 != "1_555")) continue;
            atom_connects.push_back(std::make_pair(lpos->fstAtom, lpos->sndAtom));
       }

       std::vector<RCSB::Residue*> residues;
       std::vector<std::pair<RCSB::Residue*, std::string> > res_atom_pair_list;
       std::vector<std::vector<RCSB::Atom*> > all_pair_lists, pair_lists;
       std::set<std::string> allAtomNameSet, bondedAtomNameSet;

       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            chain->GetFirstResidueList(residues);
            while (!residues.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                      if ((*rpos)->ResName() == "HOH" || (*rpos)->ResName() == "DOD" || SeqCodeUtil::is_a_standard_residue((*rpos)->ResName())) continue;
                      all_pair_lists.clear();
                      try {
                           const ConnectFormat& drug = _ccDic->find_drug((*rpos)->ResName());
                           const std::vector<std::vector<std::string> >& bondlist = drug.bonds();
                           if (bondlist.empty()) continue;

                           bondedAtomNameSet.clear();
                           for (std::vector<std::vector<std::string> >::const_iterator bpos = bondlist.begin(); bpos != bondlist.end(); ++bpos) {
                                res_atom_pair_list.clear();
                                res_atom_pair_list.push_back(std::make_pair(*rpos, (*bpos)[0]));
                                res_atom_pair_list.push_back(std::make_pair(*rpos, (*bpos)[1]));
                                _get_atom_pair_list(res_atom_pair_list, pair_lists);
                                if (pair_lists.empty()) continue;

                                bondedAtomNameSet.insert((*bpos)[0]);
                                bondedAtomNameSet.insert((*bpos)[1]);
                                for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = pair_lists.begin(); apos != pair_lists.end(); ++apos) {
                                     if ((!checking_linkage_flag) || (checking_linkage_flag && BondUtil::is_a_bond_loose((*apos)[0], (*apos)[1]))) {
                                          all_pair_lists.push_back(*apos);
                                     }
                                }
                           }

                           if (checking_linkage_flag) {
                                allAtomNameSet.clear();
                                const std::vector<RCSB::Atom*>& atomList = (*rpos)->atoms();
                                for (std::vector<RCSB::Atom*>::const_iterator apos = atomList.begin(); apos != atomList.end(); ++apos) {
                                     allAtomNameSet.insert((*apos)->atmtype());
                                }
                                if (allAtomNameSet != bondedAtomNameSet) (*rpos)->GetAllBonds(all_pair_lists);
                           }
                      } catch (const std::exception& exc) {
                           if (checking_linkage_flag) (*rpos)->GetAllBonds(all_pair_lists);
                      }

                      for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = all_pair_lists.begin(); apos != all_pair_lists.end(); ++apos) {
                           atom_connects.push_back(std::make_pair((*apos)[0], (*apos)[1]));
                      }
                 }
                 chain->GetNextResidueList(residues);
            }
            chain = _molecules[0]->GetNextChain();
       }
}
