/*
FILE:     NA_Classification.C
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
#include <math.h>

#include <set>

#include "BasePairInfo.h"
#include "SeqCodeUtil.h"

#define PROC  1
#define USED  2

typedef struct {
        int fst_chn_index;
        int fst_res_index;
        int snd_chn_index;
        int snd_res_index;
        std::string SymOP;
        int flag;
} _SIMPLE_BSPAIR;

typedef struct {
        std::string base;
        int  chn_idx;
        int  res_idx;
        std::string SymOP;
        std::vector<int> base_pair;
        int  is_break;
} BASE_PAIR_INFO;

typedef struct {
        char label;
        int  first;
        int  last;
        int  index;
} HELIX_INFO;

const static char stchpa[27] = "abcdefghijklmnopqrstuvwxyz";
const static char stchap[27] = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

typedef struct {
        int type;
        const char *type_string;
} _CLASSIFICATION;

#define NUM_CLASSIFICATION    14

const static _CLASSIFICATION classification_mapping[NUM_CLASSIFICATION] = {
       { DOUBLE_HELIX,    "double helix"         },
       { A_DOUBLE_HELIX,  "a-form double helix"  },
       { B_DOUBLE_HELIX,  "b-form double helix"  },
       { Z_DOUBLE_HELIX,  "z-form double helix"  },
       { PARALLEL_HELIX,  "parallel strands"     },
       { HAIRPIN_LOOP,    "hairpin loop"         },
       { TETRA_LOOP,      "tetraloop"            },
       { BULGE_LOOP,      "bulge loop"           },
       { MISMATCH,        "mismatched base pair" },
       { INTERNAL_LOOP,   "internal loop"        },
       { TRIPLE_HELIX,    "triple helix"         },
       { QUADRUPLE_HELIX, "quadruple helix"      },
       { THREE_STEM,      "three-way junction"   },
       { FOUR_STEM,       "four-way junction"    }
};

static void na_sec_cls(const std::vector<BASE_PAIR_INFO> &pair_info, int &classification, std::vector<HELIX_INFO> &helix_info);
static int whstrd(const std::vector<std::vector<int> > &strcd, const int stsht, const int finsht, int &nxstr, int &nxpl, const int lststd);
static int strdpl(const int stpl, const int finpl, const int thcode, const std::vector<std::vector<int> > &strand);
static void fdshus(const int stpt, const std::vector<std::vector<char> > &code, const std::vector<int> &shtcod, int &stsht, int &finsht, const int n_pair_info);
static int fdshps(const int stpt, const int shtnum, const std::vector<std::vector<char> > &code, std::vector<int> &shtcod,
                  const std::vector<std::vector<int> > &bridpt, const int n_pair_info);
static void setsht(const int refpl, const int shtnum, std::vector<int> &shtcod, const std::vector<std::vector<char> > &code, const int n_pair_info);
static int is_same_strand(const std::vector<char> &code_a, const std::vector<char> &code_b);
static int junction(int &count, const char label, const int fst, const std::vector<HELIX_INFO> &helix_info);
static bool is_tetraloop(const std::string &a1, const std::string &a2, const std::string &a3, const std::string &a4);
static bool non_watson_crick_pair(const std::string &a, const std::string &b);

void BasePairInfo::getClassification(std::vector<std::string>& classifications)
{
       classifications.clear();

       _find_additional_classification();

       for (int i = 0; i < NUM_CLASSIFICATION; ++i) {
            if (_classification & classification_mapping[i].type) {
                 classifications.push_back(classification_mapping[i].type_string);
            }
       }
}

void BasePairInfo::_find_additional_classification()
{
       if (_basepair_contacts.empty()) return;

       /* _BSPAIR defined base std::pair info. based on hydrogen bonding. There could be 
        * multiple hydrogen bondings for each pair. Remove the redundancy here.
        */
       std::vector<_SIMPLE_BSPAIR> base_info;
       base_info.clear();
       std::set<std::string> _BaseIndex;
       _BaseIndex.clear();
       _SIMPLE_BSPAIR bp;

       for (std::multimap<BPIndex, CONTACT>::const_iterator mpos = _basepair_contacts.begin(); mpos != _basepair_contacts.end(); ++mpos) {
            std::string SymOP = "";
            if (mpos->second.sym > 1 || mpos->second.lx || mpos->second.ly ||
                mpos->second.lz) {
                 SymOP = String::IntToString(mpos->second.sym) + "_" + String::IntToString(mpos->second.lx + 5)
                       + String::IntToString(mpos->second.ly + 5) + String::IntToString(mpos->second.lz + 5);
            }
            std::string cs = mpos->second.a_atm->pdb_chnid() + "_" + mpos->second.a_atm->pdb_resnam() + "_" + mpos->second.a_atm->pdb_resnum() + "_"
                           + mpos->second.a_atm->ins_code() + "_" + mpos->second.b_atm->pdb_chnid() + "_" + mpos->second.b_atm->pdb_resnam() + "_"
                           + mpos->second.b_atm->pdb_resnum() + "_" + mpos->second.b_atm->ins_code() + "_" + SymOP;
            if (_BaseIndex.find(cs) != _BaseIndex.end()) continue;

            _BaseIndex.insert(cs);

            RCSB::Residue* a_res = _mol->find_pdb_residue(mpos->second.a_atm->pdb_chnid(), mpos->second.a_atm->pdb_resnam(),
                                                          mpos->second.a_atm->pdb_resnum(), mpos->second.a_atm->ins_code());
            if (!a_res) continue;

            RCSB::Residue* b_res = _mol->find_pdb_residue(mpos->second.b_atm->pdb_chnid(), mpos->second.b_atm->pdb_resnam(),
                                                          mpos->second.b_atm->pdb_resnum(), mpos->second.b_atm->ins_code());
            if (!b_res) continue;

            bp.fst_chn_index = a_res->chain_index();
            bp.fst_res_index = a_res->position();
            bp.snd_chn_index = b_res->chain_index();
            bp.snd_res_index = b_res->position();
            bp.SymOP = SymOP;
            bp.flag = 0;
            base_info.push_back(bp);
       }

       if (base_info.empty()) return;

       std::set<IntegerStringOrder> chain_count_unique;
       std::multimap<IntegerStringOrder, int> chain_count;
       std::pair<std::multimap<IntegerStringOrder, int>::iterator,
            std::multimap<IntegerStringOrder, int>::iterator> range;
       std::multimap<IntegerStringOrder, int>::iterator mpos;
       std::set<IntegerStringOrder>::iterator spos; 
       IntegerStringOrder Index;
       std::set<int> second_chain;
       std::vector<BASE_PAIR_INFO> pair_info;
       BASE_PAIR_INFO bp_info;
       std::vector<RCSB::Residue*> residue_list;
       std::map<std::string, int> pair_info_index;
       std::vector<HELIX_INFO> helix_info;
       while (true) {
            chain_count_unique.clear();
            chain_count.clear();
            int fst_chn_index = -1;
            /* 
             * find first unprocessed chain and related chain(s) which
             * form base-pair with it
             */
            for (unsigned int i = 0; i < base_info.size(); ++i) {
                 if (base_info[i].flag == USED) continue;
                 if (base_info[i].fst_chn_index != fst_chn_index && fst_chn_index == -1) {
                      base_info[i].flag = PROC;
                      fst_chn_index = base_info[i].fst_chn_index;
                      Index.set_int_value(base_info[i].snd_chn_index);
                      Index.set_str_value(base_info[i].SymOP);
                      chain_count_unique.insert(Index);
                      chain_count.insert(std::make_pair(Index, i));
                 } else if (base_info[i].fst_chn_index == fst_chn_index) {
                      base_info[i].flag = PROC;
                      Index.set_int_value(base_info[i].snd_chn_index);
                      Index.set_str_value(base_info[i].SymOP);
                      chain_count_unique.insert(Index);
                      chain_count.insert(std::make_pair(Index, i));
                 }
            }

            if (chain_count_unique.empty()) break;

            second_chain.clear();
            /*
             * Remove isolated pairs. If the std::pair number of two chains is less than
             * one tenth of total std::pair number, these pairs will be removed.
             */
            for (spos = chain_count_unique.begin(); spos != chain_count_unique.end(); ++spos) {
                 if (chain_count.count(*spos) * 10 > chain_count.size()) {
                      second_chain.insert(spos->int_value());
                      continue;
                 }
                 range = chain_count.equal_range(*spos);
                 for (mpos = range.first; mpos != range.second; ++mpos) {
                      base_info[mpos->second].flag = USED;
                 }
            }

            /*
             * find all possible non symmetry-related pairs for second chains
             */
            for (unsigned int i = 0; i < base_info.size(); ++i) {
                 if (base_info[i].flag == USED || base_info[i].flag == PROC) continue;
                 if (second_chain.find(base_info[i].fst_chn_index) != second_chain.end() ||
                     (second_chain.find(base_info[i].snd_chn_index) != second_chain.end() && base_info[i].SymOP.empty())) base_info[i].flag = PROC;
            }

            /*
             * find all chains which are connected by base std::pair
             */
            chain_count_unique.clear();
            second_chain.clear();
            for (unsigned int i = 0; i < base_info.size(); ++i) {
                 if (base_info[i].flag != PROC) continue;
                 Index.set_int_value(base_info[i].fst_chn_index);
                 Index.set_str_value("");
                 chain_count_unique.insert(Index);
                 Index.set_int_value(base_info[i].snd_chn_index);
                 Index.set_str_value(base_info[i].SymOP);
                 chain_count_unique.insert(Index);
                 second_chain.insert(base_info[i].fst_res_index);
                 second_chain.insert(base_info[i].snd_res_index);
            }

            if (chain_count_unique.empty()) break;

            /*
             * put all found chains into a linear chain
             */
            pair_info.clear();
            pair_info_index.clear();
            for (spos = chain_count_unique.begin(); spos != chain_count_unique.end(); ++spos) {
                 bool is_removed = false;
                 RCSB::Chain* chn = _mol->GetIndexChain(spos->int_value(), is_removed);
                 if (!chn) continue;
                 for (unsigned int i = 0; i < chn->SeqLen(); ++i) {
                      if (chn->SeqRes(i)->ResIndex < 0) continue;
                      std::string base = chn->SeqRes(i)->Field[0];
                      bp_info.chn_idx = spos->int_value();
                      bp_info.res_idx = -1;
                      bp_info.SymOP = spos->str_value();
                      bp_info.base_pair.clear();
                      bp_info.is_break = 0;
                      chn->GetResidueListByIndex(chn->SeqRes(i)->ResIndex, residue_list);
                      for (std::vector<RCSB::Residue*>::iterator vpos = residue_list.begin(); vpos != residue_list.end(); ++vpos) { 
                           if (second_chain.find((*vpos)->index()) != second_chain.end()) {
                                bp_info.res_idx = (*vpos)->index();
                                base = (*vpos)->ResName();
                                break;
                           }
                      }
                      // get_parent_name(base);
                      bp_info.base = base;
                      std::string cs = String::IntToString(bp_info.chn_idx) + "_" + String::IntToString(bp_info.res_idx) + "_" + bp_info.SymOP;
                      pair_info_index.insert(std::make_pair(cs, pair_info.size()));
                      pair_info.push_back(bp_info);
                 }
                 if (!pair_info.empty()) pair_info[pair_info.size() - 1].is_break = 1;
            }

            /*
             * mapping std::pair info. into pair_info list
             */ 
            bool found = false;
            for (unsigned int i = 0; i < base_info.size(); ++i) {
                 if (base_info[i].flag != PROC) continue;
                 base_info[i].flag = USED;

                 std::string cs = String::IntToString(base_info[i].fst_chn_index) + "_" + String::IntToString(base_info[i].fst_res_index) + "_";
                 std::map<std::string, int>::iterator pos1 = pair_info_index.find(cs);
                 cs = String::IntToString(base_info[i].snd_chn_index) + "_" + String::IntToString(base_info[i].snd_res_index) + "_" + base_info[i].SymOP;
                 std::map<std::string, int>::iterator pos2 = pair_info_index.find(cs);
                 if (pos1 != pair_info_index.end() && pos2 != pair_info_index.end()) {
                      found = true;
                      pair_info[pos1->second].base_pair.push_back(pos2->second);
                      pair_info[pos2->second].base_pair.push_back(pos1->second);
                 }
            }

            if (!found) break;

            na_sec_cls(pair_info, _classification, helix_info);
            if (helix_info.empty()) continue;

            for (unsigned int i = 0; i < helix_info.size(); ++i) {
                 if (helix_info[i].label == ' ') continue;
                 for (unsigned int j = i + 1; j < helix_info.size(); ++j) {
                      if (helix_info[j].label != helix_info[i].label) continue;
                      helix_info[i].label = ' ';
                      helix_info[j].label = ' ';
                      for (int k = 0; k < helix_info[i].last-helix_info[i].first+1; ++k) {
                           int t = helix_info[i].first + k;
                           int t1 = pair_info[t].base_pair[0];
                           if (non_watson_crick_pair(pair_info[t].base, pair_info[t1].base)) _classification |= MISMATCH;
                      }
                      break;
                 }
            }
       }
}

static void na_sec_cls(const std::vector<BASE_PAIR_INFO> &pair_info, int &classification, std::vector<HELIX_INFO> &helix_info)
{
       helix_info.clear();

       std::vector<int> shtcod, data, data1;
       std::vector<std::vector<int> > strcd, bridpt;
       std::vector<char> data2;
       std::vector<std::vector<char> > code;
       shtcod.clear();
       strcd.clear();
       bridpt.clear();
       code.clear();

       int length = 4;
       for (unsigned int i = 0; i < pair_info.size(); ++i) {
            if ((int) pair_info[i].base_pair.size() > length)
                 length = pair_info[i].base_pair.size();
       }

       data.clear();
       data1.clear();
       data2.clear();
       for (int i = 0; i < length; ++i) {
            data.push_back(0);
            data1.push_back(-1);
       }
       for (int i = 0; i < 7; ++i) data2.push_back(' ');

       for (unsigned int i = 0; i < pair_info.size(); ++i) {
            shtcod.push_back(0);
            strcd.push_back(data);
            code.push_back(data2);
            bridpt.push_back(data1);
       }

       int nxstrn = 0;
       for (int i = 0; i < (int) pair_info.size() - 1; ++i) {
            for (int k = 0; k < (int) pair_info[i].base_pair.size(); ++k) {
                 int j = pair_info[i].base_pair[k];
                 if (j < i) continue;
                 if (pair_info[i+1].is_break) continue;

                 for (int l = 0; l < (int) pair_info[i+1].base_pair.size(); ++l) {
                      int j1 = pair_info[i+1].base_pair[l];
                      int dir = 0;
                      if (j1 == (j+1)) dir = 1;
                      else if (j1 == (j-1)) dir = -1;
                      if (dir) {
                           if (!strcd[i][k]) {
                                nxstrn++;
                                strcd[i][k] = nxstrn * dir;
                           }
                           strcd[i+1][l] = strcd[i][k];
                           for (unsigned int m = 0; m < pair_info[j].base_pair.size(); ++m) {
                                if (pair_info[j].base_pair[m] == i) {
                                     strcd[j][m] = strcd[i][k];
                                     break;
                                }
                           }
                           for (unsigned int m = 0; m < pair_info[j1].base_pair.size(); ++m) {
                                if (pair_info[j1].base_pair[m] == (i+1)) {
                                     strcd[j1][m] = strcd[i+1][l];
                                     break;
                                }
                           }
                           code[i][0] = 'H';
                           code[i+1][0] = 'H';
                           code[j][0] = 'H';
                           code[j1][0] = 'H';
                      }
                 }
            }
       }

       int i = 0;
       while (i < (int) pair_info.size()) {
            while (i < (int) pair_info.size() && code[i][0] != 'H') i++;
            if (i >= (int) pair_info.size()) break;
            int stsht = i;
            i++;
            while (i < (int) pair_info.size() && code[i][0] == 'H') i++;
            int finsht = i;

            int lststd = 0;
            while (1) {
                 int nxstr = 0;
                 int nxpl = 0;
                 int bstval = whstrd(strcd, stsht, finsht, nxstr, nxpl, lststd);
                 if (!bstval) break;
                 lststd = bstval;
                 int chcode = lststd % 26;
                 if (chcode == 0) chcode = 26;
                 char strdch = stchpa[chcode-1];
                 if (strcd[nxstr][nxpl] < 0)
                      strdch = stchap[chcode-1];
                 int stchpl = 1;
                 while (stchpl < 4 && code[nxstr][stchpl] != ' ') stchpl++;
                 if (stchpl == 1) {
                      int thstpl = nxstr;
                      while (1) {
                           int nxstpl = strdpl(thstpl, finsht, strcd[nxstr][nxpl], strcd);
                           if (nxstpl < 0) break;
                           while (stchpl < 4 && code[nxstpl][stchpl] != ' ') stchpl++;
                           thstpl = nxstpl;
                           if (stchpl > 1) break;
                      }
                 }
                 code[nxstr][stchpl] = strdch;
                 bridpt[nxstr][stchpl-1] = pair_info[nxstr].base_pair[nxpl];
                 int thstpl = nxstr;
                 while (1) {
                      int nxstpl = strdpl(thstpl, finsht, strcd[nxstr][nxpl], strcd);
                      if (nxstpl < 0) break;
                      code[nxstpl][stchpl] = strdch;
                      int strpl = 0;
                      while (strpl < 4 && strcd[nxstr][nxpl] != strcd[nxstpl][strpl]) strpl++;
                      bridpt[nxstpl][stchpl-1] = pair_info[nxstpl].base_pair[strpl];
                      thstpl = nxstpl;
                 }
            }
       }

       int shtnum = 0;
       i = 0;
       int stsht = -1, finsht = -1;
       while (i < (int) pair_info.size()) {
            fdshus(i, code, shtcod, stsht, finsht, (int) pair_info.size());
            if (stsht < 0) break;
            if (i == stsht && i == finsht) break;
            shtnum++;
            for (int j = stsht; j < finsht; ++j) {
                 shtcod[j] = shtnum;
                 for (int k = 0; k < 4; ++k) {
                      if (bridpt[j][k] >= 0 && !shtcod[bridpt[j][k]])
                           setsht(bridpt[j][k], shtnum, shtcod, code,
                               (int) pair_info.size());
                 }
            }
            while (fdshps(finsht, shtnum, code, shtcod, bridpt, (int) pair_info.size())) ;
            i = finsht;
       }

       if (!shtnum) return;

       std::vector<std::vector<char> > helix;
       helix.clear();
       for (i = 0; i < shtnum; ++i) {
            data2.clear();
            for (unsigned int j = 0; j < pair_info.size(); ++j) {
                 if (shtcod[j] != (i+1)) continue;
                 for (int k = 1; k < 5; ++k) {
                      if (code[j][k] == ' ') continue;
                      bool found = false;
                      for (unsigned l = 0; l < data2.size(); ++l) {
                           if (code[j][k] == data2[l]) {
                                found = true;
                                break;
                           }
                      }
                      if (!found) data2.push_back(code[j][k]);
                 }
            }
            helix.push_back(data2);
       }

       for (i = 0; i < shtnum; ++i) {
            if (helix[i].size() == 3 || helix[i].size() == 4)
                 classification |= QUADRUPLE_HELIX;
            else if (helix[i].size() == 2)
                 classification |= TRIPLE_HELIX;
            else if (helix[i].size() == 1 && helix[i][0] >= 'a' &&
                     helix[i][0] <= 'z')
                 classification |= PARALLEL_HELIX;
       }

       int first = -1;
       int last = -1;
       int last_h = -1;
       char label = ' ';
       int index = 1;

       HELIX_INFO h_info;
       for (i = 0; i < (int) pair_info.size(); ++i) {
            if (pair_info[i].is_break || code[i][1] == ' ') {
                 if (first >= 0 && last_h >= 0 && (last_h-first+1) >= 3) {
                      h_info.label = label;
                      h_info.first = first;
                      h_info.last  = last_h;
                      h_info.index = index;
                      helix_info.push_back(h_info);
                 }
                 first = -1; last_h = -1; label = ' ';
                 if (pair_info[i].is_break) {
                      index++;
                      last = -1;
                 }
                 if (code[i][1] == ' ') continue;
            }
            if (code[i][2] != ' ' || code[i][3] != ' ' || code[i][4] != ' ') continue;

            if (code[i][1] != label) {
                 if (first >= 0 && last_h >= 0 && (last_h-first+1) >= 3) {
                      h_info.label = label;
                      h_info.first = first;
                      h_info.last  = last_h;
                      h_info.index = index;
                      helix_info.push_back(h_info);
                 }
                 if (last >= 0) {
                      int last_j = pair_info[last].base_pair[0];
                      int ij = pair_info[i].base_pair[0];
                      if ((last_j >= last || last_j >= i) && (ij >= last || ij >= i) && pair_info[last_j].chn_idx == pair_info[ij].chn_idx &&
                          pair_info[last_j].SymOP == pair_info[ij].SymOP) {
                           if (last_j == i && ij == last) {
                                if ((i - last) == 3 && is_tetraloop(pair_info[i-3].base, pair_info[i-2].base, pair_info[i-1].base, pair_info[i].base))
                                     classification |= TETRA_LOOP;
                                else classification |= HAIRPIN_LOOP;
                           } else {
                                int diff_i = i - last;
                                int diff_j = abs(last_j - ij);
                                if (diff_i == 1 || diff_j == 1)
                                     classification |= BULGE_LOOP;
                                else if (diff_i == 2 && diff_j == 2)
                                     classification |= MISMATCH;
                                else {
                                     int start = last_j;
                                     int end = ij;
                                     if (start > end) {
                                          start = ij;
                                          end = last_j;
                                     }
                                     bool found = false;
                                     for (int j = start + 1; j < end; ++j) {
                                          if (code[j][0] == 'H') {
                                               found = true;
                                               break;
                                          }
                                     }
                                     if (!found) classification |= INTERNAL_LOOP;
                                }
                           }
                      }
                 }
                 label = code[i][1];
                 first = last_h = i;
                 last = -1;
            } else { last = i; last_h = i; }
       }
       if (first >= 0 && last_h >= 0 && (last_h-first+1) >= 3) {
            h_info.label = label;
            h_info.first = first;
            h_info.last  = last_h;
            h_info.index = index;
            helix_info.push_back(h_info);
       }

       for (i = 0; i < (int) helix_info.size() - 1; ++i) {
            int count = 0;
            if (junction(count, helix_info[i].label, i, helix_info)) {
                 if (count == 6) classification |= THREE_STEM;
                 else if (count == 8) classification |= FOUR_STEM;
            }
       }
}

static int whstrd(const std::vector<std::vector<int> > &strcd, const int stsht, const int finsht, int &nxstr, int &nxpl, const int lststd)
{
       int bstval = 0;
       for (int i = stsht; i < finsht; ++i) {
            for (int j = 0; j < 4; ++j) {
                 if (abs(strcd[i][j]) > lststd) {
                      if (bstval == 0) {
                           bstval = abs(strcd[i][j]);
                           nxstr = i;
                           nxpl = j;
                      } else if (abs(strcd[i][j]) < bstval) {
                           bstval = abs(strcd[i][j]);
                           nxstr = i;
                           nxpl = j;
                      }
                 }
            }
       }
       return bstval;
}

static int strdpl(const int stpl, const int finpl, const int thcode, const std::vector<std::vector<int> > &strand)
{
       int nxpl = -1;
       for (int i = stpl + 1; i < finpl; ++i) {
            bool found = false;
            for (int j = 0; j < 4; ++j) {
                 if (strand[i][j] == thcode) {
                      nxpl = i;
                      found = true;
                      break;
                 }
            }
            if (found) break;
       }
       return nxpl;
}

static void fdshus(const int stpt, const std::vector<std::vector<char> > &code, const std::vector<int> &shtcod, int &stsht, int &finsht, const int n_pair_info)
{
       int nxres = stpt;
       stsht = -1;
       while (nxres < n_pair_info && (code[nxres][0] == ' ' ||
               shtcod[nxres] != 0)) nxres++;
       if (nxres >= n_pair_info) return;
       stsht = nxres;
       while (nxres < n_pair_info && is_same_strand(code[stsht], code[nxres]))
            nxres++;
       finsht = nxres;
}

static int fdshps(const int stpt, const int shtnum, const std::vector<std::vector<char> > &code, std::vector<int> &shtcod,
                  const std::vector<std::vector<int> > &bridpt, const int n_pair_info)
{
       int nxres = stpt;
       while (nxres < n_pair_info && shtcod[nxres] != (-shtnum)) nxres++;
       if (nxres >= n_pair_info) return 0;
       int stsht = nxres;
       while (nxres < n_pair_info && shtcod[nxres] == (-shtnum)) nxres++;
       int finsht = nxres;
       for (int j = stsht; j < finsht; ++j) {
            shtcod[j] = shtnum;
            for (int k = 0; k < 4; ++k) {
                 if (bridpt[j][k] >= 0 && !shtcod[bridpt[j][k]]) setsht(bridpt[j][k], shtnum, shtcod, code, n_pair_info);
            }
       }
       return 1;
}

static void setsht(const int refpl, const int shtnum, std::vector<int> &shtcod, const std::vector<std::vector<char> > &code, const int n_pair_info)
{
       int nxres = refpl;
       // while (nxres < n_pair_info && code[nxres][0] != ' ') nxres++;
       while (nxres < n_pair_info && is_same_strand(code[refpl], code[nxres])) nxres++;
       int finsht = nxres;
       nxres = refpl;
       // while (nxres >= 0 && code[nxres][0] != ' ') nxres--;
       while (nxres >= 0 && is_same_strand(code[refpl], code[nxres])) nxres--;
       int stsht = nxres + 1;
       for (int i = stsht; i < finsht; ++i) shtcod[i] = -shtnum;
}

static int is_same_strand(const std::vector<char> &code_a, const std::vector<char> &code_b)
{
       int is_same = 0;
       for (int i = 1; i < 5; ++i) {
            if (code_a[i] != ' ' && code_a[i] == code_b[i]) {
                 is_same = 1;
                 break;
            }
       }
       return is_same;
}

static int junction(int &count, const char label, const int fst, const std::vector<HELIX_INFO> &helix_info)
{
       if ((fst + 1) >= (int) helix_info.size()) return 0;
       if (helix_info[fst].index != helix_info[fst+1].index) return 0;
       count += 2;
       if (helix_info[fst+1].label == label) return 1;

       bool found = false;
       int j = 0;
       for (int i = fst + 2; i < (int) helix_info.size(); ++i) {
            if (helix_info[fst+1].label == helix_info[i].label) {
                 found = true;
                 j = i;
                 break;
            }
       }
       if (!found) return 0;
       return (junction(count, label, j, helix_info));
}

static bool is_tetraloop(const std::string &a1, const std::string &a2, const std::string &a3, const std::string &a4)
{
       bool found = false;
       if (SeqCodeUtil::isSameResName(a1, "G") && SeqCodeUtil::isSameResName(a4, "A") &&
          (SeqCodeUtil::isSameResName(a3, "A") || SeqCodeUtil::isSameResName(a3, "G"))) found = true;
       if (SeqCodeUtil::isSameResName(a1, "U") && SeqCodeUtil::isSameResName(a4, "G") && SeqCodeUtil::isSameResName(a3, "C")) found = true;
       if (SeqCodeUtil::isSameResName(a1, "C") && SeqCodeUtil::isSameResName(a4, "G") && SeqCodeUtil::isSameResName(a2, "U")) found = true;
       return found;
}

static bool non_watson_crick_pair(const std::string &a, const std::string &b)
{
       if (SeqCodeUtil::isSameResName(a, "A") && (SeqCodeUtil::isSameResName(b, "T") || SeqCodeUtil::isSameResName(b, "U"))) return false;
       if (SeqCodeUtil::isSameResName(b, "A") && (SeqCodeUtil::isSameResName(a, "T") || SeqCodeUtil::isSameResName(a, "U"))) return false;
       if ((SeqCodeUtil::isSameResName(a, "G") && SeqCodeUtil::isSameResName(b, "C")) ||
           (SeqCodeUtil::isSameResName(a, "C") && SeqCodeUtil::isSameResName(b, "G")))
            return false;
       return true;
}
