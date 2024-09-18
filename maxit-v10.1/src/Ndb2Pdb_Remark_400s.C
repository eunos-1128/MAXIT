/*
FILE:     Ndb2Pdb_Remark_400s.C
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

#include "Maxit.h"
#include "NdbToken.h"
#include "PdbWrite.h"
#include "Remark_define.h"
#include "Remark_extern.h"
#include "utillib.h"

static void get_set(const std::list<std::vector<std::string> >& rlist,
                    std::set<std::string>& index);

void Maxit::_get_missing_or_zero_occupancy_residues_or_atoms_for_remarks()
{
       _missing_polymer_residues.clear();
       _missing_polymer_atoms.clear();
       _zero_occ_polymer_residues.clear();
       _zero_occ_polymer_atoms.clear();
       _missing_non_polymer_residues.clear();
       _zero_occ_non_polymer_residues.clear();

       _get_missing_or_zero_occupancy_residues_or_atoms();

       _analysis_missing_or_zero_occupancy_residues();
       _analysis_missing_or_zero_occupancy_atoms();
}

void Maxit::_analysis_missing_or_zero_occupancy_residues()
{
       if (_unobs_or_zero_occ_residues.empty()) return;

       std::vector<std::string> data;
       std::list<std::vector<std::string> > tlist;

       for (std::list<std::vector<std::string> >::const_iterator
            lpos = _unobs_or_zero_occ_residues.begin();
            lpos != _unobs_or_zero_occ_residues.end(); ++lpos) {
            int mol_id = atoi((*lpos)[0].c_str());

            data.clear();
            data.push_back((*lpos)[3]);
            data.push_back((*lpos)[4]);
            data.push_back((*lpos)[5]);
            data.push_back((*lpos)[6]);
            if ((*lpos)[1] == "Y") {
                 if ((*lpos)[2] == "1") {
                      std::map<int, std::list<std::vector<std::string> > >::iterator
                          mpos = _missing_polymer_residues.find(mol_id);
                      if (mpos != _missing_polymer_residues.end())
                           mpos->second.push_back(data);
                      else {
                           tlist.clear();
                           tlist.push_back(data);
                           _missing_polymer_residues.insert(std::make_pair(mol_id, tlist));
                      }
                 } else if ((*lpos)[2] == "0") {
                      std::map<int, std::list<std::vector<std::string> > >::iterator
                          mpos = _zero_occ_polymer_residues.find(mol_id);
                      if (mpos != _zero_occ_polymer_residues.end())
                           mpos->second.push_back(data);
                      else {
                           tlist.clear();
                           tlist.push_back(data);
                           _zero_occ_polymer_residues.insert(std::make_pair(mol_id, tlist));
                      }
                 }
            } else if ((*lpos)[1] == "N") {
                 if ((*lpos)[2] == "1") {
                      std::map<int, std::list<std::vector<std::string> > >::iterator
                          mpos = _missing_non_polymer_residues.find(mol_id);
                      if (mpos != _missing_non_polymer_residues.end())
                           mpos->second.push_back(data);
                      else {
                           tlist.clear();
                           tlist.push_back(data);
                           _missing_non_polymer_residues.insert(std::make_pair(mol_id, tlist));
                      }
                 } else if ((*lpos)[2] == "0") {
                      std::map<int, std::list<std::vector<std::string> > >::iterator
                          mpos = _zero_occ_non_polymer_residues.find(mol_id);
                      if (mpos != _zero_occ_non_polymer_residues.end())
                           mpos->second.push_back(data);
                      else {
                           tlist.clear();
                           tlist.push_back(data);
                           _zero_occ_non_polymer_residues.insert(std::make_pair(mol_id, tlist));
                      }
                 }
            }
       }
}

void Maxit::_analysis_missing_or_zero_occupancy_atoms()
{
       if (_unobs_or_zero_occ_atoms.empty()) return;

       std::vector<std::string> data;
       std::list<std::vector<std::string> > tlist;

       for (std::list<std::vector<std::string> >::const_iterator
            lpos = _unobs_or_zero_occ_atoms.begin();
            lpos != _unobs_or_zero_occ_atoms.end(); ++lpos) {
            int mol_id = atoi((*lpos)[0].c_str());

            data.clear();
            data.push_back((*lpos)[3]);
            data.push_back((*lpos)[4]);
            data.push_back((*lpos)[5]);
            data.push_back((*lpos)[6]);
            if ((*lpos)[1] == "Y") {
                 if ((*lpos)[2] == "1") {
                      std::map<int, std::list<std::vector<std::string> > >::iterator
                          mpos = _missing_polymer_atoms.find(mol_id);
                      if (mpos != _missing_polymer_atoms.end()) {
                           std::vector<std::string>& field = mpos->second.back();
                           if (field[0] == data[0] && field[1] == data[1] &&
                               field[2] == data[2] && field[3] == data[3]) {
                               if ((*lpos)[7] != field[field.size() - 1])
                                    field.push_back((*lpos)[7]);
                           } else {
                                data.push_back((*lpos)[7]);
                                mpos->second.push_back(data);
                           }
                      } else {
                           data.push_back((*lpos)[7]);
                           tlist.clear();
                           tlist.push_back(data);
                           _missing_polymer_atoms.insert(std::make_pair(mol_id, tlist));
                      }
                 } else if ((*lpos)[2] == "0") {
                      std::map<int, std::list<std::vector<std::string> > >::iterator
                          mpos = _zero_occ_polymer_atoms.find(mol_id);
                      if (mpos != _zero_occ_polymer_atoms.end()) {
                           std::vector<std::string>& field = mpos->second.back();
                           if (field[0] == data[0] && field[1] == data[1] &&
                               field[2] == data[2] && field[3] == data[3]) {
                               if ((*lpos)[7] != field[field.size() - 1])
                                    field.push_back((*lpos)[7]);
                           } else {
                                data.push_back((*lpos)[7]);
                                mpos->second.push_back(data);
                           }
                      } else {
                           data.push_back((*lpos)[7]);
                           tlist.clear();
                           tlist.push_back(data);
                           _zero_occ_polymer_atoms.insert(std::make_pair(mol_id, tlist));
                      }
                 }
            } else if ((*lpos)[1] == "N") {
                 if ((*lpos)[2] == "1") {
                      std::map<int, std::list<std::vector<std::string> > >::iterator
                          mpos = _missing_non_polymer_residues.find(mol_id);
                      if (mpos != _missing_non_polymer_residues.end()) {
                           std::vector<std::string>& field = mpos->second.back();
                           if (field != data) mpos->second.push_back(data);
                      } else {
                           tlist.clear();
                           tlist.push_back(data);
                           _missing_non_polymer_residues.insert(std::make_pair(mol_id, tlist));
                      }
                 } else if ((*lpos)[2] == "0") {
                      std::map<int, std::list<std::vector<std::string> > >::iterator
                          mpos = _zero_occ_non_polymer_residues.find(mol_id);
                      if (mpos != _zero_occ_non_polymer_residues.end()) {
                           std::vector<std::string>& field = mpos->second.back();
                           if (field != data) mpos->second.push_back(data);
                      } else {
                           tlist.clear();
                           tlist.push_back(data);
                           _zero_occ_non_polymer_residues.insert(std::make_pair(mol_id, tlist));
                      }
                 }
            }
       }
}

void Maxit::_ndb_to_pdb_get_remark(const int& remark_no, const int& field_no,
                                   const std::string& keyword)
{
       std::string block_remark;
       _getRecordFront("ENTDTL", field_no, block_remark);
       if (block_remark.empty()) return;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
       int width = ndbformat.FieldList[2].FieldWidth;

       std::vector<std::string> remark_array;
       get_text_array_from_block(remark_array, block_remark, width);
       if (remark_array.empty()) return;

       _addNewRemark(remark_no, keyword);
       if (remark_no == 600) _addNewRemark(remark_no, "");
       _ndb_to_pdb_add_remark(remark_no, remark_array);
}

void Maxit::_ndb_to_pdb_get_remark_465()
{
       if (_missing_polymer_residues.empty()) return;

       bool need_new_nmr_template = _check_new_nmr_template(_missing_polymer_residues);

       if (need_new_nmr_template) {
            _ndb_to_pdb_get_general_remark(Num_Remark_465_NMR, Remark_465_NMR, 0, 0, 0);
            if (_molecules.size() > 1) {
                 std::string cs1 = String::IntToString(_molecules[0]->Mol_ID());
                 std::string cs2 = String::IntToString(_molecules[_molecules.size()-1]->Mol_ID());
                 _addNewRemark(465, "  MODELS " + cs1 + "-" + cs2);
            }
            _addNewRemark(465, "    RES C SSSEQI");
       } else _ndb_to_pdb_get_general_remark(Num_Remark_465, Remark_465, 0, 0, 0);

       for (std::map<int, std::list<std::vector<std::string> > >::const_iterator mpos =
            _missing_polymer_residues.begin(); mpos != _missing_polymer_residues.end(); ++mpos) {
            std::string mol_id = "";
            if (_missing_polymer_residues.size() > 1 && !need_new_nmr_template)
                 mol_id = String::IntToString(mpos->first);
            _ndb_to_pdb_get_remark_465_475(mol_id, 465, mpos->second);
            if (need_new_nmr_template) break;
       }
}

void Maxit::_ndb_to_pdb_get_remark_470()
{
       if (_missing_polymer_atoms.empty()) return;

       bool need_new_nmr_template = _check_new_nmr_template(_missing_polymer_atoms);

       if (need_new_nmr_template) {
            _ndb_to_pdb_get_general_remark(Num_Remark_470_NMR, Remark_470_NMR, 0, 0, 0);
            if (_molecules.size() > 1) {
                 std::string cs1 = String::IntToString(_molecules[0]->Mol_ID());
                 std::string cs2 = String::IntToString(_molecules[_molecules.size()-1]->Mol_ID());
                 _addNewRemark(470, "  MODELS " + cs1 + "-" + cs2);
            }
            _addNewRemark(470, "    RES CSSEQI  ATOMS");
       } else _ndb_to_pdb_get_general_remark(Num_Remark_470, Remark_470, 0, 0, 0);

       for (std::map<int, std::list<std::vector<std::string> > >::const_iterator mpos =
            _missing_polymer_atoms.begin(); mpos != _missing_polymer_atoms.end(); ++mpos) {
            std::string mol_id = "";
            if (_missing_polymer_atoms.size() > 1 && !need_new_nmr_template)
                 mol_id = String::IntToString(mpos->first);
            _ndb_to_pdb_get_remark_470_480(mol_id, 470, mpos->second);
            if (need_new_nmr_template) break;
       }
}

void Maxit::_ndb_to_pdb_get_remark_475()
{
       if (_zero_occ_polymer_residues.empty()) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_475, Remark_475, 0, 0, 0);

       for (std::map<int, std::list<std::vector<std::string> > >::const_iterator mpos =
            _zero_occ_polymer_residues.begin(); mpos != _zero_occ_polymer_residues.end(); ++mpos) {
            std::string mol_id = "";
            if (_zero_occ_polymer_residues.size() > 1)
                 mol_id = String::IntToString(mpos->first);
            _ndb_to_pdb_get_remark_465_475(mol_id, 475, mpos->second);
       }
}

void Maxit::_ndb_to_pdb_get_remark_480()
{
       if (_zero_occ_polymer_atoms.empty()) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_480, Remark_480, 0, 0, 0);

       for (std::map<int, std::list<std::vector<std::string> > >::const_iterator mpos =
            _zero_occ_polymer_atoms.begin(); mpos != _zero_occ_polymer_atoms.end(); ++mpos) {
            std::string mol_id = "";
            if (_zero_occ_polymer_atoms.size() > 1)
                 mol_id = String::IntToString(mpos->first);
            _ndb_to_pdb_get_remark_470_480(mol_id, 480, mpos->second);
       }
}

bool Maxit::_check_new_nmr_template(const std::map<int, std::list<std::vector<std::string> > >& r_maps)
{
       if (!(_experiment_type & EXPERIMENT_TYPE_NMR) &&
           !(_experiment_type & EXPERIMENT_TYPE_NMR_SOLID)) 
            return false;

       if (r_maps.size() < 2) return true;

       std::set<std::string> Index1, Index2;
       Index1.clear(); Index2.clear();

       std::map<int, std::list<std::vector<std::string> > >::const_iterator
           mpos = r_maps.begin();

       get_set(mpos->second, Index1);

       ++mpos;
       while (mpos != r_maps.end()) {
            get_set(mpos->second, Index2);
            if (Index2 != Index1) return false; 
            ++mpos;
       }

       return true;
}

void Maxit::_ndb_to_pdb_get_remark_465_475(const std::string& mol_id, const int& remark_no,
                                           const std::list<std::vector<std::string> >& l_list)
{
       int width = 5;
       if (remark_no == 610) width = 4;
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = l_list.begin(); lpos != l_list.end(); ++lpos) {
            std::string remark = "   ";
            if (!mol_id.empty())
                 remark = FormattedFieldValue(mol_id, 3, 3, 0, false, true);
            std::string chain_id = " ";
            if (!(*lpos)[0].empty()) chain_id = (*lpos)[0].substr(0, 1);
            std::string insert_code = " ";
            if (!(*lpos)[3].empty()) insert_code = (*lpos)[3].substr(0, 1);
            remark += " " + FormattedFieldValue((*lpos)[1], 3, 3, 0, false, true) + " "
                    + chain_id + " " + FormattedFieldValue((*lpos)[2], 3, width, 0,
                    false, true) + insert_code;
            _addNewRemark(remark_no, remark);
       }
}

void Maxit::_ndb_to_pdb_get_remark_470_480(const std::string& mol_id, const int& remark_no,
                                           const std::list<std::vector<std::string> >& l_list)
{
       std::string before = "";
       std::string after  = " ";
       if (remark_no == 480) {
            before = " ";
            after  = "";
       }
        
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = l_list.begin(); lpos != l_list.end(); ++lpos) {
            std::string head = "   ";
            if (!mol_id.empty()) head = FormattedFieldValue(mol_id, 3, 3, 0, false, true);
            std::string chain_id = " ";
            if (!(*lpos)[0].empty()) chain_id = (*lpos)[0].substr(0, 1);
            std::string insert_code = " ";
            if (!(*lpos)[3].empty()) insert_code = (*lpos)[3].substr(0, 1);
            head += " " + FormattedFieldValue((*lpos)[1], 3, 3, 0, false, true) + " "
                  + chain_id + before + FormattedFieldValue((*lpos)[2], 3, 4, 0, false, true)
                  + insert_code + after;
            std::string context = "";
            int num = 0;
            for (unsigned int i = 4; i < lpos->size(); ++i) {
                 if ((*lpos)[i].empty()) continue;
                 if (num == 7) {
                      if (!context.empty()) _addNewRemark(remark_no, head + context);
                      context.clear();
                      num = 0;
                 }
                 context += " " + printAtomNameField(_ccDic, "", (*lpos)[i], (*lpos)[1]);
                 num++;
            }
            if (!context.empty()) _addNewRemark(remark_no, head + context);
       }
}

void Maxit::_ndb_to_pdb_get_remark_900()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("RENTRY");
       if (ppos == _pdb_records.end()) return;

       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
       int width = ndbformat.FieldList[2].FieldWidth;

       _addNewRemark(900, "RELATED ENTRIES");

       std::vector<std::string> data;
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
            _addNewRemark(900, "RELATED ID: " + (*lpos)[3] + "   RELATED DB: " + (*lpos)[2]); 

            if ((*lpos)[4].empty()) continue;

            std::string buffer = (*lpos)[4];
            String::StripAndCompressWs(buffer);
            if (buffer.empty()) continue;

            String::UpperCase(buffer);
            get_max_length_words(data, buffer, width);
            _ndb_to_pdb_add_remark(900, data);
       }
}

static void get_set(const std::list<std::vector<std::string> >& rlist,
                    std::set<std::string>& index)
{
       index.clear();

       std::set<std::string> t_set;
       std::string cs;
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = rlist.begin(); lpos != rlist.end(); ++lpos) {
            t_set.clear();
            for (std::vector<std::string>::const_iterator
                 pos = lpos->begin(); pos != lpos->end(); ++pos) {
                 t_set.insert(*pos);
            }
            cs.clear();
            for (std::set<std::string>::const_iterator
                 spos = t_set.begin(); spos != t_set.end(); ++spos) {
                 if (!cs.empty()) cs += "_";
                 cs += *spos;
            }
            index.insert(cs);
       }
}
