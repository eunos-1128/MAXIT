/*
FILE:     Pdb2Ndb_Remark.C
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
#include "Remark_global.h"
#include "utillib.h"

#define REFMAC_TEMPLATE  1
#define PHENIX_TEMPLATE  2

#define NUM_REMARK  25

static const int remark_nums[NUM_REMARK] = { 0, 2, 3, 100, 200, 210, 220, 225, 230, 240,
       245, 250, 265, 280, 300, 350, 400, 450, 465, 600, 650, 700, 800, 900, 999
};

#define NUM_REFMAC_NCS_GROUP_TOKEN  6

static const char *_refmac_ncs_group_token[NUM_REFMAC_NCS_GROUP_TOKEN] = {
       "RFTPOS", "RFMPOS", "RFLPOS", "RFTTHR", "RFMTHR", "RFLTHR"
};

#define NUM_MERGE_EXP_TOKEN  12

static const char *_merge_exp_tokens[NUM_MERGE_EXP_TOKEN] = {
       "DIFPTL", "DTMEAS", "DTMETH", "DTTEMP", "PERCOM", "RADIAT", "REFLEC", "RESTOT", "SHELC", "WAVLEN", "XFILE1", "XFILE2"
};

#define NUM_MERGE_ID_TOKEN 4

static const char *_merge_diffrn_id_tokens[NUM_MERGE_ID_TOKEN] = {
       "PERCOM", "REFLEC", "RESTOT", "SHELC"
};

#define NUM_XFILE_TOKEN  2

static const char *_xfile_tokens[NUM_XFILE_TOKEN] = {
       "XFILE1", "XFILE2"
};

#define PHENIX_REPEAT_TOKEN  3

static const char *_phenix_repeat_tokens[PHENIX_REPEAT_TOKEN] = {
       "REFREM", "RMSTN1", "RMSTN2"
};

void Maxit::_pdb_to_ndb_process_remarks()
{
       _get_remark_map("REMARK", 2, _remarks);
       if (_remarks.empty()) return;

       int diffrn_id = 1;
       std::map<std::string, std::list<std::vector<std::string> > > merge_records;
       merge_records.clear();

       for (int i = 0; i < NUM_REMARK; ++i) {
            std::map<int, std::vector<std::string> >::iterator
                mpos = _remarks.find(remark_nums[i]);
            if (mpos == _remarks.end()) continue;

            if (remark_nums[i] == 0) {
                 if (!_DepUI_Flag) _pdb_to_ndb_process_block_remark(remark_nums[i], mpos->second);
            } else if (remark_nums[i] == 2)
                 _pdb_to_ndb_process_remark(remark_nums[i], mpos->second, 1, &Remark_2, 0);
            else if (remark_nums[i] == 3)
                 _pdb_to_ndb_process_remark_3(mpos->second);
            else if (remark_nums[i] == 100) {
                 bool bmrb_flag = false;
                 bool pdbe_flag = false;
                 bool rcsb_flag = false;
                 for (std::vector<std::string>::iterator
                      pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                      if ((*pos)[pos->size() - 1] == '.') pos->erase(pos->size() - 1);
                      if (pos->find("THE BMRB ID CODE IS") != std::string::npos)
                           bmrb_flag = true;
                      else if (pos->find("THE PDBE ID CODE IS") != std::string::npos)
                           pdbe_flag = true;
                      else if (pos->find("THE RCSB ID CODE IS") != std::string::npos)
                           rcsb_flag = true;
                 }
                 if (bmrb_flag)
                      _pdb_to_ndb_process_remark(remark_nums[i], mpos->second,
                                     Num_Remark_100_BMRB, Remark_100_BMRB, 0);
                 else if (pdbe_flag)
                      _pdb_to_ndb_process_remark(remark_nums[i], mpos->second,
                                     Num_Remark_100_EBI, Remark_100_EBI, 0);
                 else if (rcsb_flag)
                      _pdb_to_ndb_process_remark(remark_nums[i], mpos->second,
                                     Num_Remark_100_RCSB, Remark_100_RCSB, 0);
            } else if (remark_nums[i] == 200) {
                 _pdb_to_ndb_process_remark(remark_nums[i], mpos->second, Num_Remark_200,
                                            Remark_200, 0);
                 std::string exp_type = "X-RAY DIFFRACTION";
                 if (_experiment_type & EXPERIMENT_TYPE_FIBER)
                      exp_type = "FIBER DIFFRACTION";
                 diffrn_id = _pdb_to_ndb_process_experimental_TOKENs(diffrn_id, exp_type, "x-ray", merge_records);
                 if (_experiment_type == EXPERIMENT_TYPE_BASIC)
                      _experiment_type = EXPERIMENT_TYPE_XRAY;
            } else if (remark_nums[i] == 210) {
                 _pdb_to_ndb_process_remark(remark_nums[i], mpos->second, Num_Remark_210,
                                            Remark_210, 0);
                 if (_experiment_type == EXPERIMENT_TYPE_BASIC)
                      _experiment_type = EXPERIMENT_TYPE_NMR;
            } else if (remark_nums[i] == 215) {
                 if (_experiment_type == EXPERIMENT_TYPE_BASIC)
                      _experiment_type = EXPERIMENT_TYPE_NMR;
            } else if (remark_nums[i] == 217) {
                 if (_experiment_type == EXPERIMENT_TYPE_BASIC)
                      _experiment_type = EXPERIMENT_TYPE_NMR_SOLID;
            } else if (remark_nums[i] == 220) {
                 _pdb_to_ndb_process_remark(remark_nums[i], mpos->second, Num_Remark_220,
                                            Remark_220, 0);
                 if (_experiment_type == EXPERIMENT_TYPE_BASIC)
                      _experiment_type = EXPERIMENT_TYPE_MODEL;
            } else if (remark_nums[i] == 225) {
                 if (_experiment_type == EXPERIMENT_TYPE_BASIC)
                       _experiment_type = EXPERIMENT_TYPE_MODEL;
            } else if (remark_nums[i] == 230) {
                 _pdb_to_ndb_process_remark(remark_nums[i], mpos->second, Num_Remark_230,
                                            Remark_230, 0);
                 diffrn_id = _pdb_to_ndb_process_experimental_TOKENs(diffrn_id, "NEUTRON DIFFRACTION", "neutron", merge_records);
                 if (_experiment_type == EXPERIMENT_TYPE_BASIC)
                      _experiment_type = EXPERIMENT_TYPE_NEUTRON;
            } else if (remark_nums[i] == 240) {
                 _pdb_to_ndb_process_remark(remark_nums[i], mpos->second, Num_Remark_240,
                                            Remark_240, 0);
                 diffrn_id = _pdb_to_ndb_process_experimental_TOKENs(diffrn_id, "ELECTRON CRYSTALLOGRAPHY", "electron", merge_records);
                 if (_experiment_type == EXPERIMENT_TYPE_BASIC)
                      _experiment_type = EXPERIMENT_TYPE_ELECTRON;
            } else if (remark_nums[i] == 245)
                 _pdb_to_ndb_process_remark(remark_nums[i], mpos->second, Num_Remark_245,
                                            Remark_245, 0);
            else if (remark_nums[i] == 250)
                 _pdb_to_ndb_process_remark(remark_nums[i], mpos->second, Num_Remark_250,
                                            Remark_250, 0);
            else if (remark_nums[i] == 265)
                 _pdb_to_ndb_process_remark_265(mpos->second);
            else if (remark_nums[i] == 280)
                 _pdb_to_ndb_process_remark(remark_nums[i], mpos->second, Num_Remark_280,
                                            Remark_280, 0);
            else if (remark_nums[i] == 300) {
                 if ((_experiment_type & EXPERIMENT_TYPE_EM) || (_experiment_type & EXPERIMENT_TYPE_ET)) {
                      _pdb_to_ndb_process_repeat_remark(300, mpos->second, 1,
                                                        Remark_Cryo_Em_300);
                 }
            } else if (remark_nums[i] == 350)
                 _pdb_to_ndb_process_remark_350(mpos->second);
            else if (remark_nums[i] == 400)
                 _pdb_to_ndb_process_pdbx_entry_details_remark(mpos->second, 1, "COMPOUND");
            else if (remark_nums[i] == 450)
                 _pdb_to_ndb_process_pdbx_entry_details_remark(mpos->second, 2, "SOURCE");
            else if (remark_nums[i] == 465)
                 _pdb_to_ndb_process_remark_465(mpos->second);
            else if (remark_nums[i] == 600)
                 _pdb_to_ndb_process_pdbx_entry_details_remark(mpos->second, 3, "HETEROGEN");
            else if (remark_nums[i] == 650)
                 _pdb_to_ndb_process_block_remark(remark_nums[i], mpos->second);
            else if (remark_nums[i] == 700)
                 _pdb_to_ndb_process_block_remark(remark_nums[i], mpos->second);
            else if (remark_nums[i] == 800) {
                 std::vector<std::vector<std::string> > bRemarkCard;
                 std::set<std::string> template_tokens;
                 template_tokens.clear();
                 for (int j = 0; j < Num_Remark_800; ++j) {
                      for (int k = 0; k < Remark_800[j].NumField; ++k) {
                           if (strcmp(Remark_800[j].FieldList[k].Text, "") == 0) continue;
                           template_tokens.insert(Remark_800[j].FieldList[k].Text);
                      }
                 }
                 _pdb_to_ndb_get_block_remarks(mpos->second, template_tokens, bRemarkCard);
                 for (std::vector<std::vector<std::string> >::const_iterator
                      vpos = bRemarkCard.begin(); vpos != bRemarkCard.end(); ++vpos) {
                      _pdb_to_ndb_process_remark(remark_nums[i], *vpos, Num_Remark_800,
                                                 Remark_800, 1);
                 }
            } else if (remark_nums[i] == 900)
                 _pdb_to_ndb_process_900(mpos->second);
            else if (remark_nums[i] == 999)
                 _pdb_to_ndb_process_pdbx_entry_details_remark(mpos->second, 4, "SEQUENCE");
       }

       _pdb_records.erase("REMARK");

       if (!merge_records.empty()) {
            for (std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = merge_records.begin();
                 mpos != merge_records.end(); ++mpos) {
                 _pdb_records.erase(mpos->first);
                 _pdb_records.insert(std::make_pair(mpos->first, mpos->second));
            }
       }

       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           xpos = _pdb_records.find("XFILES");
       if (xpos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator
                 lxpos = xpos->second.begin(); lxpos != xpos->second.end(); ++lxpos) {
                 if ((*lxpos)[2].empty() && (*lxpos)[3].empty() && (*lxpos)[4].empty()) {
                      lxpos = xpos->second.erase(lxpos);
                      --lxpos;
                 }
            }
            if (xpos->second.empty()) _pdb_records.erase("XFILES");
       }

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("SITELS");
       if (ppos == _pdb_records.end() || _site.empty()) return;

       std::map<std::string, std::pair<std::string, std::string> > site_meta_data;
       site_meta_data.clear();

       for (std::list<std::vector<std::string> >::const_iterator
            pos = ppos->second.begin(); pos != ppos->second.end(); ++pos) {
            site_meta_data.insert(std::make_pair((*pos)[2],
                                  std::make_pair((*pos)[4], (*pos)[3])));
       }

       for (std::list<_SITE>::iterator
            lpos = _site.begin(); lpos != _site.end(); ++lpos) {
            std::map<std::string, std::pair<std::string, std::string> >::const_iterator
                mpos = site_meta_data.find(lpos->SiteID);
            if (mpos == site_meta_data.end()) continue;
            lpos->code = mpos->second.first;
            lpos->details = mpos->second.second;
       }
}

void Maxit::_pdb_to_ndb_process_remark_3(const std::vector<std::string>& remarkCards)
{
       if ((_experiment_type & EXPERIMENT_TYPE_EM) || (_experiment_type & EXPERIMENT_TYPE_ET))
            _pdb_to_ndb_process_remark(3, remarkCards, Num_Remark_3_Cryo_Em_In,
                                      Remark_3_Cryo_Em_In, 0);
       else if (_experiment_type & EXPERIMENT_TYPE_SOLN_SCT)
            _pdb_to_ndb_process_repeat_remark(3, remarkCards, 0, Remark_3_Solution);
       else if (_experiment_type & EXPERIMENT_TYPE_NMR ||
                _experiment_type & EXPERIMENT_TYPE_NMR_SOLID ||
                _experiment_type & EXPERIMENT_TYPE_MODEL)
            _pdb_to_ndb_process_remark(3, remarkCards, Num_Remark_3_NMR, Remark_3_NMR, 0);
       else {
            std::vector<std::vector<std::string> > bRemarkCard;
            std::set<std::string> template_tokens;
            template_tokens.clear();
            for (int n = 0; n < NUM_GROUP_REMARKS; ++n) {
                 for (int i = 0; i < group_remarks_in[n]->Remarks_No; ++i) {
                      for (int j = 0; j < group_remarks_in[n]->num_remarks[i]; ++j) {
                           for (int k = 0; k < group_remarks_in[n]->remarks[i][j].NumField; ++k) {
                                if (strcmp(group_remarks_in[n]->remarks[i][j].FieldList[k].Text, "") == 0) continue;
                                template_tokens.insert(group_remarks_in[n]->remarks[i][j].FieldList[k].Text);
                           }
                      }
                 }
            }
            for (int l = 0; l < NUM_METHOD; ++l) template_tokens.insert(_method_remark[l]);
            _pdb_to_ndb_get_block_remarks(remarkCards, template_tokens, bRemarkCard);
            if (bRemarkCard.empty()) return;

            for (std::vector<std::vector<std::string> >::const_iterator
                 bpos = bRemarkCard.begin(); bpos != bRemarkCard.end(); ++bpos) {
                 if (_pdb_to_ndb_process_remark(3, *bpos, Num_Remark_3_1, Remark_3_1, 0)) break;
                 if (_pdb_to_ndb_process_remark(3, *bpos, Num_Remark_3_2, Remark_3_2, 0)) break;
                 if (_pdb_to_ndb_process_remark(3, *bpos, Num_Remark_3_3, Remark_3_3, 0)) break;
            }
            int group_index = _pdb_to_ndb_get_refinement_index();

            if (group_index == 5)
                 _pdb_to_ndb_get_block_remarks(bRemarkCard, REFMAC_TEMPLATE);
            else if (group_index == 10)
                 _pdb_to_ndb_get_block_remarks(bRemarkCard, PHENIX_TEMPLATE);

            _pdb_to_ndb_process_group_remark_3(group_index, bRemarkCard);
       }
}

void Maxit::_pdb_to_ndb_process_repeat_remark(const int& remark_no, const std::vector<std::string>&
                                              remarkCards, const unsigned int& start_block,
                                              const GROUP_REMARKS& group_template)
{
       std::vector<std::vector<std::string> > bRemarkCard;
       std::set<std::string> template_tokens;
       template_tokens.clear();
       for (int i = 0; i < group_template.Remarks_No; ++i) {
            for (int j = 0; j < group_template.num_remarks[i]; ++j) {
                 for (int k = 0; k < group_template.remarks[i][j].NumField; ++k) {
                      if (strcmp(group_template.remarks[i][j].FieldList[k].Text, "") == 0)
                           continue;
                      template_tokens.insert(group_template.remarks[i][j].FieldList[k].Text);
                 }
            }
       }
       _pdb_to_ndb_get_block_remarks(remarkCards, template_tokens, bRemarkCard);
       if (bRemarkCard.empty()) return;

       int k = 0;
       for (unsigned int j = start_block; j < bRemarkCard.size(); ++j) {
            bool found = false;
            while (k < group_template.Remarks_No) {
                 found = _pdb_to_ndb_process_remark(remark_no, bRemarkCard[j],
                                     group_template.num_remarks[k],
                                     group_template.remarks[k],
                                     group_template.repeat[k]);
                 if (found) break;
                 k++;
            }
            if (found) {
                 if (k < group_template.Remarks_No &&
                     !group_template.repeat[k]) k++;
            } else k = 0;
       }
}

void Maxit::_pdb_to_ndb_get_block_remarks(const std::vector<std::string>& remarkCards,
                                          const std::set<std::string>& template_tokens,
                                          std::vector<std::vector<std::string> >& blocks)
{
       blocks.clear();

       std::vector<std::vector<std::string> > tmp_blocks;
       tmp_blocks.clear();

       std::vector<std::string> data;
       data.clear();

       for (std::vector<std::string>::const_iterator
            pos = remarkCards.begin(); pos != remarkCards.end(); ++pos) {
            if (pos->empty()) {
                 if (!data.empty()) tmp_blocks.push_back(data);
                 data.clear();
            } else data.push_back(*pos);
       }
       if (!data.empty()) tmp_blocks.push_back(data);

       if (tmp_blocks.empty()) return;

       blocks.push_back(tmp_blocks[0]);

       if (tmp_blocks.size() == 1) return;

       for (unsigned int i = 1; i < tmp_blocks.size(); ++i) {
            bool break_point = false;
            unsigned int count = 0;
            for (std::vector<std::string>::const_iterator
                 pos = tmp_blocks[i].begin(); pos != tmp_blocks[i].end(); ++pos) {
                 std::string cs = *pos;
                 if (cs.find(":") != std::string::npos) count++;
                 String::StripLeadingWs(cs);
                 for (std::set<std::string>::const_iterator spos =
                      template_tokens.begin(); spos != template_tokens.end(); ++spos) {
                      if (cs.substr(0, spos->size()) == *spos) {
                           break_point = true;
                           break;
                      }
                 }
                 if (break_point) break;
            }
            if (!break_point) {
                 if ((count * 2) >= tmp_blocks[i].size()) break_point = true;
            }
            if (break_point) blocks.push_back(tmp_blocks[i]);
            else {
                  blocks[blocks.size()-1].push_back("");
                  for (std::vector<std::string>::const_iterator
                       pos = tmp_blocks[i].begin(); pos != tmp_blocks[i].end(); ++pos) {
                       blocks[blocks.size()-1].push_back(*pos);
                  }
            }
       }
}

void Maxit::_pdb_to_ndb_get_block_remarks(std::vector<std::vector<std::string> >& blocks,
                                          const int& templete_type)
{
       std::vector<std::vector<std::string> > new_blocks;
       new_blocks.clear();

       std::vector<std::string> data;
       for (std::vector<std::vector<std::string> >::const_iterator
            ppos = blocks.begin(); ppos != blocks.end(); ++ppos) {
            data.clear();
            bool has_new_selection = false;
            bool has_new_torsion = false;
            bool has_limit_rmsd = false;
            bool has_restaint_rmsd = false;
            bool has_histogram = false;
            for (std::vector<std::string>::const_iterator
                 pos = ppos->begin(); pos != ppos->end(); ++pos) {
                 if ((templete_type == REFMAC_TEMPLATE &&
                     (pos->find("NCS OPERATOR :") != std::string::npos ||
                      pos->find("TWIN DOMAIN   :") != std::string::npos)) ||
                     (templete_type == PHENIX_TEMPLATE &&
                     (pos->find("NCS OPERATOR :") != std::string::npos ||
                      pos->find("NCS GROUP :") != std::string::npos ||
                      pos->find("TLS GROUP :") != std::string::npos))) {
                      if (templete_type == PHENIX_TEMPLATE && has_new_selection &&
                         (has_histogram || has_new_torsion || has_limit_rmsd ||
                          has_restaint_rmsd)) data.clear();
                      if (!data.empty()) new_blocks.push_back(data);
                      data.clear();
                 }
                 if (pos->find("SELECTION          :") != std::string::npos)
                      has_new_selection = true;
                 else if (pos->find("RESTRAINED TORSIONS:") != std::string::npos)
                      has_new_torsion = true;
                 else if (pos->find("BELOW LIMIT RMSD   :") != std::string::npos)
                      has_limit_rmsd = true;
                 else if (pos->find("ALL RESTRAINT RMSD :") != std::string::npos)
                      has_restaint_rmsd = true;
                 else if (pos->find("HISTOGRAM OF DIFFERENCES") != std::string::npos)
                      has_histogram = true;
                 data.push_back(*pos);
            }
            if (templete_type == PHENIX_TEMPLATE && has_new_selection && (has_histogram ||
                has_new_torsion || has_limit_rmsd || has_restaint_rmsd)) data.clear();
            if (!data.empty()) new_blocks.push_back(data);
       }
       blocks = new_blocks;
}

void Maxit::_pdb_to_ndb_process_group_remark_3(const int& group_index, const
                                std::vector<std::vector<std::string> >& bRemarkCard)
{
       _current_method.clear();

       GROUP_REMARKS *group_remark = group_remarks_in[group_index];

       std::vector<int> used_index;
       used_index.reserve(group_remark->Remarks_No);
       used_index.clear();
       for (int j = 0; j < group_remark->Remarks_No; ++j) used_index.push_back(0);

       int k = 0;
       for (unsigned int j = 1; j < bRemarkCard.size(); ++j) {
            if (bRemarkCard[j].size() == 1) {
                 bool found_method = false;
                 for (int l = 0; l < NUM_METHOD; ++l) {
                      if (bRemarkCard[j][0].find(_method_remark[l]) != std::string::npos) {
                           if (!_current_method.empty()) _updateRefineMethod(_current_method);
                           _current_method = _method_name[l];
                           _joint_methods.push_back(_current_method);
                           found_method = true;
                           break;
                      }
                 }
                 if (found_method) continue;
            }

            bool found = false;
            while (k < group_remark->Remarks_No) {
                 found = _pdb_to_ndb_process_remark(3, bRemarkCard[j],
                               group_remark->num_remarks[k], group_remark->remarks[k],
                               group_remark->repeat[k]);
                 if (found) {
                      used_index[k] = 1;
                      if (group_remark->remarks[k] == Remark_3_REFMAC_NCS_GROUP)
                           _pdb_to_ndb_update_REFMAC_NCS_GROUP();
                      break;
                 }
                 k++;
            }
            if (found) {
                 if (k >= group_remark->Remarks_No || (k < group_remark->Remarks_No &&
                     group_remark->repeat[k] == 0)) {
                      for (k = 0; k < group_remark->Remarks_No; ++k) {
                           if (!used_index[k]) break;
                      }
                 }
                 if (group_index == 10) k = 0;
            } else k = 0;
       }

       if (!_current_method.empty()) _updateRefineMethod(_current_method);
       if (_joint_methods.size() > 1) {
            _orderRecordWithFloatKey("SHELL", 1);
            const ndb_token_format& ndbformat = NdbToken::getTokenFormat("SHELL");
            _orderRecordWithOrderKey("SHELL", ndbformat.RefineField - 1, _joint_methods);
            if (group_index == 10) {
                 std::set<std::string> method_set;
                 std::list<std::vector<std::string> > vlist;
                 for (int i = 0; i < PHENIX_REPEAT_TOKEN; ++i) {
                      std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = _pdb_records.find(_phenix_repeat_tokens[i]);
                      if (mpos == _pdb_records.end()) continue;
                      method_set.clear();
                      const ndb_token_format& ndbformat1 = NdbToken::getTokenFormat(_phenix_repeat_tokens[i]);
                      for (std::list<std::vector<std::string> >::const_iterator lmpos = mpos->second.begin(); lmpos != mpos->second.end(); ++lmpos) {
                           if (!(*lmpos)[ndbformat1.RefineField - 1].empty()) method_set.insert((*lmpos)[ndbformat1.RefineField - 1]);
                      }
                      if (method_set.size() == _joint_methods.size()) continue;

                      vlist.clear();
                      for (std::vector<std::string>::const_iterator jpos = _joint_methods.begin(); jpos != _joint_methods.end(); ++jpos) {
                           for (std::list<std::vector<std::string> >::const_iterator lmpos = mpos->second.begin(); lmpos != mpos->second.end(); ++lmpos) {
                                std::vector<std::string> vec = *lmpos;
                                vec[ndbformat1.RefineField - 1] = *jpos;
                                vlist.push_back(vec);
                           }
                      }
                      _pdb_records.erase(_phenix_repeat_tokens[i]);
                      _pdb_records.insert(std::make_pair(_phenix_repeat_tokens[i], vlist));
                 }
            }
       }
}

int Maxit::_pdb_to_ndb_process_experimental_TOKENs(const int& diffrn_id, const std::string& type, const std::string& scattering_type,
                                                   std::map<std::string, std::list<std::vector<std::string> > >& merge_records)
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator ppos = _pdb_records.find("CRMETH");
       if (ppos != _pdb_records.end()) {
            for (std::list<std::vector<std::string> >::iterator lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 (*lpos)[5] = type;
            }
       }

       int return_id = _pdb_to_ndb_update_200_align(diffrn_id);

       std::string diffrn_ids = "";
       for (int i = diffrn_id; i <= return_id; ++i) {
            if (!diffrn_ids.empty()) diffrn_ids += ",";
            diffrn_ids += String::IntToString(i);
       }

       for (int i = 0; i < NUM_MERGE_ID_TOKEN; ++i) {
            ppos = _pdb_records.find(_merge_diffrn_id_tokens[i]);
            if (ppos == _pdb_records.end()) continue;

            const ndb_token_format& ndbformat = NdbToken::getTokenFormat(_merge_diffrn_id_tokens[i]);
            for (std::list<std::vector<std::string> >::iterator lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 (*lpos)[ndbformat.SeqField - 1] = diffrn_ids;
            }
       }
       
       for (int i = 0; i < NUM_XFILE_TOKEN; ++i) {
            ppos = _pdb_records.find(_xfile_tokens[i]);
            if (ppos == _pdb_records.end()) continue;

            const ndb_token_format& ndbformat = NdbToken::getTokenFormat(_xfile_tokens[i]);
            for (std::list<std::vector<std::string> >::iterator lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 (*lpos)[ndbformat.RefineField - 1] = type;
            }
       }

       bool found_merge = false;
       for (int i = 0; i < NUM_MERGE_EXP_TOKEN; ++i) {
            std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = _pdb_records.find(_merge_exp_tokens[i]);
            if (mpos == _pdb_records.end()) continue;

            found_merge = true;
            ppos = merge_records.find(_merge_exp_tokens[i]);
            if (ppos == merge_records.end())
                 merge_records.insert(std::make_pair(_merge_exp_tokens[i], mpos->second));
            else {
                 for (std::list<std::vector<std::string> >::const_iterator lpos = mpos->second.begin(); lpos != mpos->second.end(); ++lpos) {
                      ppos->second.push_back(*lpos);
                 }
            }
            _pdb_records.erase(_merge_exp_tokens[i]);
       }

       if (found_merge) {
            std::vector<std::string> FieldInfo;
            std::list<std::vector<std::string> > rlist;
            const ndb_token_format& ndbformat = NdbToken::getTokenFormat("SCTYPE");
            for (int id = diffrn_id; id <= return_id; ++id) {
                 FieldInfo.clear();
                 for (int i = 0; i < ndbformat.NumField; ++i) FieldInfo.push_back("");
                 FieldInfo[0] = "SCTYPE";
                 FieldInfo[ndbformat.SeqField - 1] = String::IntToString(id);
                 FieldInfo[2] = scattering_type;

                 ppos = merge_records.find("SCTYPE");
                 if (ppos == merge_records.end()) {
                      rlist.clear();
                      rlist.push_back(FieldInfo);
                      merge_records.insert(std::make_pair("SCTYPE", rlist));
                 } else ppos->second.push_back(FieldInfo);
            }
       }

       return (return_id + 1);
}

void Maxit::_pdb_to_ndb_process_remark_265(const std::vector<std::string>& remarkCards)
{
       if (!(_experiment_type & EXPERIMENT_TYPE_SOLN_SCT)) return;

       _pdb_to_ndb_process_repeat_remark(265, remarkCards, 0, Remark_265);

       int num = 0;
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("SOLEXP");
       if (ppos != _pdb_records.end()) {
            std::list<std::vector<std::string> >::iterator prev = ppos->second.end();
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 num++;
                 (*lpos)[1] = String::IntToString(num);
                 if ((*lpos)[2].empty() && prev != ppos->second.end() &&
                     !(*prev)[2].empty()) (*lpos)[2] = (*prev)[2];
                 prev = lpos; 
            }
       }

       ppos = _pdb_records.find("SOLMDL");
       if (ppos != _pdb_records.end()) {
            int j = 0;
            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 num++;
                 (*lpos)[9] = String::IntToString(num);
                 j++;
                 (*lpos)[10] = String::IntToString(j);

                 _addNewRecord("SOLEXP");
                 _updateRecordBack("SOLEXP", 1, (*lpos)[9]);
                 _updateRecordBack("SOLEXP", 2, "THEORETICAL MODELLING");
            }
       }
}

void Maxit::_pdb_to_ndb_update_REFMAC_NCS_GROUP()
{
       std::string value;
       _getRecordBack("NCSGRO", 3, value);
       if (value.empty()) return;

       String::StripAndCompressWs(value);
       if (value.empty()) return;

       std::string::size_type p = value.find_first_of(" ");
       if (p != std::string::npos) value.erase(p);

       for (int i = 0; i < NUM_REFMAC_NCS_GROUP_TOKEN; ++i) {
            std::map<std::string, std::list<std::vector<std::string> > >::iterator
                ppos = _pdb_records.find(_refmac_ncs_group_token[i]);
            if (ppos == _pdb_records.end()) continue;

            for (std::list<std::vector<std::string> >::iterator
                 lpos = ppos->second.begin(); lpos != ppos->second.end(); ++lpos) {
                 if ((*lpos)[2].empty()) (*lpos)[2] = value;
            }
       }
}
