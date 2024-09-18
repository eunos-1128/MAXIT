/*
FILE:     Promotif_Sheet.C
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

#include "Promotif.h"
#include "utillib.h"

void Promotif::_fineSheets(std::string& complicated_sheet_warning)
{
       std::list<_SHEET> _complicated_sheets;
       _complicated_sheets.clear();

       std::vector<int> pshtlt;
       pshtlt.clear();
       pshtlt.reserve(_residues.size() * 2 / 3);
       int maxstr = 0;
       int start = -1;
       for (unsigned int i = 0; i < _residues.size(); ++i) {
            if (_residues[i].ssbond[10] == 'E') {
                 if (start < 0) {
                      start = i;
                      maxstr++;
                 }
                 if (_residues[i].ssbond1) {
                      bool found = false;
                      for (unsigned int j = 0; j < pshtlt.size(); ++j) {
                           if (pshtlt[j] == _residues[i].ssbond1) {
                                found = true; break;
                           }
                      }
                      if (!found) pshtlt.push_back(_residues[i].ssbond1);
                 }
            } else start = -1;
       }

       if (!maxstr) return;

       std::vector<std::vector<int> > strks, shtks;
       strks.reserve(maxstr);
       shtks.reserve(maxstr);
       std::vector<int> data;
       data.clear();
       data.push_back(0);
       data.push_back(0);
       for (int i = 0; i < maxstr; ++i) {
            strks.push_back(data);
            shtks.push_back(data);
       }

       int nstr = 0;
       for (unsigned int i = 1; i < _residues.size(); ++i) {
            if (_residues[i].ssbond[10] == 'E') {
                 if (_residues[i-1].ssbond[10] != 'E') {
                      nstr++;
                      strks[nstr-1][0] = i;
                      if (_residues[i].ssbond1)
                           shtks[nstr-1][0] = _residues[i].ssbond1;
                 } else {
                      if (_residues[i].ssbond1) {
                           if (shtks[nstr-1][0] == 0)
                                shtks[nstr-1][0] = _residues[i].ssbond1;
                           else if (shtks[nstr-1][0] != _residues[i].ssbond1)
                                shtks[nstr-1][1] = _residues[i].ssbond1;
                      }
                 }
                 if (i == (_residues.size() - 1) || (i < (_residues.size() - 1) && _residues[i+1].ssbond[10] != 'E')) {
                      strks[nstr-1][1] = i;
                 }
            }
       }

       std::map<int, std::vector<std::string> > _c_sheet;
       _getUserDefinedSheet(_c_sheet);

       std::vector<STRAND> strand;
       strand.clear();
       strand.reserve(nstr);
       for (int i = 0; i < nstr; ++i) {
            STRAND str;
            strand.push_back(str);
       }

       _SHEET _a_sheet;
       _SHEET_STRAND _a_strand;
       _SHEET_TOPOLOGY _a_order;

       bool _complicated_sheet_flag = false;
       bool exist_incorrect_index = false;
       std::string message;
       message.clear();
       char cur_atom[2], prv_atom[2];
       memset(cur_atom, 0, 2);
       memset(prv_atom, 0, 2);
       for (unsigned int i = 0; i < pshtlt.size(); ++i) {
            for (int j = 0; j < nstr; ++j) {
                 strand[j].flag = -1;
                 strand[j].length = 0;
                 strand[j].init = NULL;
                 strand[j].term = NULL;
                 strand[j].num_of_partners = 0;
                 for (int k = 0; k < MAX_PARTNERS; ++k) {
                      strand[j].partners[k] = -1;
                      strand[j].type[k] = 0;
                      strand[j].cur_atom[k] = '\0';
                      strand[j].prv_atom[k] = '\0';
                      strand[j].cur_res[k] = NULL;
                      strand[j].prv_res[k] = NULL;
                 }
                 strand[j].cur_atom[MAX_PARTNERS] = '\0';
                 strand[j].prv_atom[MAX_PARTNERS] = '\0';
            }
            int numstr = 0;
            for (int j = 0; j < 2; ++j) {
                 for (int k = 0; k < nstr; ++k) {
                      if (pshtlt[i] == shtks[k][j]) {
                           strand[numstr].flag = k;
                           strand[numstr].length = strks[k][1] - strks[k][0] + 1;
                           strand[numstr].init = _residues[strks[k][0]].res;
                           strand[numstr].term = _residues[strks[k][1]].res;
                           numstr++;
                      }
                 }
            }
            if (numstr < 2) continue;

            for (int j = 0; j < numstr - 1; ++j) {
                 for (int k = j + 1; k < numstr; ++k) {
                      int jj = strand[j].flag;
                      int kk = strand[k].flag;
                      double energy = 0;
                      int type = 0;
                      RCSB::Residue* res_i = NULL;
                      RCSB::Residue* res_j = NULL;
                      for (int l = strks[jj][0]; l <= strks[jj][1]; ++l) {
                           for (int m = 0; m < 2; ++m) {
                                int t = _residues[l].bridpt[m];
                                if (t && t >= (strks[kk][0]+1) && t <= (strks[kk][1]+1)) {
                                     if (_residues[l].ssbond[4+m] >= 'A' && _residues[l].ssbond[4+m] <= 'Z') {
                                          type = -1;
                                          break;
                                     } else if (_residues[l].ssbond[4+m] >= 'a' && _residues[l].ssbond[4+m] <= 'z') {
                                          type = 1;
                                          break;
                                     }

                                     for (int n = l - 1; n >= strks[jj][0]; --n) {
                                          if (_residues[n].ssbond[4+m] >= 'A' && _residues[n].ssbond[4+m] <= 'Z') {
                                               type = -1;
                                               break;
                                          } else if (_residues[n].ssbond[4+m] >= 'a' && _residues[n].ssbond[4+m] <= 'z') {
                                               type = 1;
                                               break;
                                          }
                                     }
                                     if (type) break;

                                     for (int n = l + 1; n <= strks[jj][1]; ++n) {
                                          if (_residues[n].ssbond[4+m] >= 'A' && _residues[n].ssbond[4+m] <= 'Z') {
                                               type = -1;
                                               break;
                                          } else if (_residues[n].ssbond[4+m] >= 'a' && _residues[n].ssbond[4+m] <= 'z') {
                                               type = 1;
                                               break;
                                          }
                                     }
                                     if (type) break;
                                }
                           }
                           if (type) break;
                      }
                      if (!type) {
                           for (int l = strks[kk][0]; l <= strks[kk][1]; ++l) {
                                for (int m = 0; m < 2; ++m) {
                                     int t = _residues[l].bridpt[m];
                                     if (t && t >= (strks[jj][0]+1) && t <= (strks[jj][1]+1)) {
                                          if (_residues[l].ssbond[4+m] >= 'A' && _residues[l].ssbond[4+m] <= 'Z') {
                                               type = -1;
                                               break;
                                          } else if (_residues[l].ssbond[4+m] >= 'a' && _residues[l].ssbond[4+m] <= 'z') {
                                               type = 1;
                                               break;
                                          }

                                          for (int n = l - 1; n >= strks[kk][0]; --n) {
                                               if (_residues[n].ssbond[4+m] >= 'A' && _residues[n].ssbond[4+m] <= 'Z') {
                                                    type = -1;
                                                    break;
                                               } else if (_residues[n].ssbond[4+m] >= 'a' && _residues[n].ssbond[4+m] <= 'z') {
                                                    type = 1;
                                                    break;
                                               }
                                          }
                                          if (type) break;

                                          for (int n = l + 1; n <= strks[kk][1]; ++n) {
                                               if (_residues[n].ssbond[4+m] >= 'A' && _residues[n].ssbond[4+m] <= 'Z') {
                                                    type = -1;
                                                    break;
                                               } else if (_residues[n].ssbond[4+m] >= 'a' && _residues[n].ssbond[4+m] <= 'z') {
                                                    type = 1;
                                                    break;
                                               }
                                          }
                                          if (type) break;
                                     }
                                }
                                if (type) break;
                           }
                      }
                      if (type) {
                           for (int l = strks[jj][0]; l <= strks[jj][1]; ++l) {
                                for (int m = 0; m < 4; ++m) {
                                     int t = _residues[l].hbond[m];
                                     if (t && t >= (strks[kk][0]+1) && t <= (strks[kk][1]+1) && _residues[l].energy[m] < energy) {
                                          energy = _residues[l].energy[m];
                                          res_i = _residues[l].res;
                                          res_j = _residues[t-1].res;
                                          if (!(m % 2)) {
                                               cur_atom[0] = 'N';
                                               prv_atom[0] = 'O';
                                          } else {
                                               cur_atom[0] = 'O';
                                               prv_atom[0] = 'N';
                                          } 
                                     }
                                }
                           }
                           if (!res_i || !res_j) {
                                for (int l = strks[kk][0]; l <= strks[kk][1]; ++l) {
                                     for (int m = 0; m < 4; ++m) {
                                          int t = _residues[l].hbond[m];
                                          if (t && t >= (strks[jj][0]+1) && t <= (strks[jj][1]+1) && _residues[l].energy[m] < energy) {
                                               energy = _residues[l].energy[m];
                                               res_i = _residues[t-1].res;
                                               res_j = _residues[l].res;
                                               if (!(m % 2)) {
                                                    cur_atom[0] = 'O';
                                                    prv_atom[0] = 'N';
                                               } else {
                                                    cur_atom[0] = 'N';
                                                    prv_atom[0] = 'O';
                                               } 
                                          }
                                     }
                                }
                           }
                      }			   

                      if (type && res_i && res_j) {
                           if (strand[j].num_of_partners < MAX_PARTNERS) {
                                strand[j].partners[strand[j].num_of_partners] = k;
                                strand[j].type[strand[j].num_of_partners] = type;
                                strand[j].cur_res[strand[j].num_of_partners] = res_i;
                                strand[j].prv_res[strand[j].num_of_partners] = res_j;
     
                                strand[k].partners[strand[k].num_of_partners] = j;
                                strand[k].type[strand[k].num_of_partners] = type;
                                strand[k].cur_res[strand[k].num_of_partners] = res_j;
                                strand[k].prv_res[strand[k].num_of_partners] = res_i;
                                strand[j].cur_atom[strand[j].num_of_partners] = cur_atom[0];
                                strand[j].prv_atom[strand[j].num_of_partners] = prv_atom[0];
                                strand[k].cur_atom[strand[k].num_of_partners] = prv_atom[0];
                                strand[k].prv_atom[strand[k].num_of_partners] = cur_atom[0];
                                strand[j].num_of_partners++;
                                strand[k].num_of_partners++;
                           }
                      }
                 }
            }

            for (int j = 0; j < numstr; ++j) {
                 if (strand[j].num_of_partners == 0) {
                      for (int k = 0; k < numstr; ++k) {
                           for (int l = 0; l < strand[k].num_of_partners; ++l) {
                                if (strand[k].partners[l] > j) {
                                     strand[k].partners[l]--;
                                }
                           }
                      }
                      for (int k = j; k < numstr - 1; ++k) {
                           strand[k].num_of_partners = strand[k+1].num_of_partners;
                           strand[k].flag            = strand[k+1].flag;
                           strand[k].length          = strand[k+1].length;
                           strand[k].init            = strand[k+1].init;
                           strand[k].term            = strand[k+1].term;
                           for (int l = 0; l < strand[k].num_of_partners; ++l) {
                                strand[k].partners[l] = strand[k+1].partners[l];
                                strand[k].type[l]     = strand[k+1].type[l];
                                strand[k].cur_atom[l] = strand[k+1].cur_atom[l];
                                strand[k].prv_atom[l] = strand[k+1].prv_atom[l];
                                strand[k].cur_res[l]  = strand[k+1].cur_res[l];
                                strand[k].prv_res[l]  = strand[k+1].prv_res[l];
                           } 
                      }
                      j--;
                      numstr--;
                 }
            }
            if (numstr < 2) continue;
 
            std::vector<std::vector<int> > index, sheet;
            index.clear();
            data.clear();
            for (int j = 0; j < 5; ++j) index.push_back(data);

            for (int j = 0; j < numstr; ++j) {
                 strand[j].flag = 0;
                 if (strand[j].num_of_partners == 1)
                      index[0].push_back(j);
                 else if (strand[j].num_of_partners == 2)
                      index[1].push_back(j);
                 else if (strand[j].num_of_partners == 3)
                      index[2].push_back(j);
                 else if (strand[j].num_of_partners == 4)
                      index[3].push_back(j);
                 else break;
            }

            sheet.clear();
            std::map<int, std::vector<std::string> >::iterator pos = _c_sheet.find(i);
            if (pos != _c_sheet.end()) {
                 std::vector<std::string> string_data;
                 for (std::vector<std::string>::const_iterator vpos = pos->second.begin(); vpos != pos->second.end(); ++vpos) {
                      get_wordarray(string_data, *vpos, " ");
                      if (string_data.empty()) continue;
                      data.clear();
                      for (::vector<std::string>::const_iterator vpos1 = string_data.begin(); vpos1 != string_data.end(); ++vpos1) {
                           data.push_back(atoi(vpos1->c_str()) - 1);
                      }
                      sheet.push_back(data);
                 }
            }

            if (sheet.empty()) {
                 if (index[0].size() == 2 && (int) index[1].size() == (numstr - 2)) {
                      _find_strands(index[0][0], strand, data);
                      if (!data.empty()) sheet.push_back(data);
                 } else if ((int) index[1].size() == numstr) {
                      _find_strands(0, strand, data);
                      data.push_back(0);
                      sheet.push_back(data);
                 } else if (index[0].size() == 3 && index[2].size() == 1 && (int) index[1].size() == (numstr - 4)) {
                      for (int j = 0; j < 3; ++j) {
                           _find_strands(index[0][j], strand, index[j + 1]);
                           strand[index[2][0]].flag = 0;
                      }
                      sheet.push_back(index[1]);
                      sheet.push_back(index[1]);
                      for (int j = 0; j < 2; ++j) {
                           for (int k = (int) index[j + 2].size() - 2; k >= 0; --k) {
                                sheet[j].push_back(index[j + 2][k]);
                           }
                      }
                 } else if (index[0].size() == 4 && index[3].size() == 1 && (int) index[1].size() == (numstr - 5)) {
                      for (int j = 0; j < 4; ++j) {
                           _find_strands(index[0][j], strand, index[j + 1]);
                           strand[index[3][0]].flag = 0;
                      }
                      for (int j = 0; j < 2; ++j) {
                           for (int k = (int) index[j + 3].size() - 2; k >=0; --k) {
                                index[j + 1].push_back(index[j + 3][k]);
                           }
                           sheet.push_back(index[j + 1]);
                      }
                 } else if (index[0].size() == 4 && index[2].size() == 2 && (int) index[1].size() == (numstr - 6)) {
                      for (int j = 0; j < 2; ++j) {
                           _find_strands(index[0][j], strand, data);
                           sheet.push_back(data);
                           for (unsigned int k = 0; k < index[2].size(); ++k) {
                                strand[index[2][k]].flag = 0;
                           }
                           _find_strands(index[0][j+2], strand, index[j + 3]);
                           for (unsigned int k = 0; k < index[2].size(); ++k) {
                                strand[index[2][k]].flag = 0;
                           }
                      }
                      _find_strands(index[2][0], strand, index[0]);
                      for (unsigned int j = 0; j < sheet.size(); ++j) {
                           if (!sheet[j].empty() && !index[0].empty()) {
                                if (sheet[j][sheet[j].size() - 1] == index[0][0]) {
                                     for (unsigned int k = 1; k < index[0].size(); ++k) {
                                          sheet[j].push_back(index[0][k]);
                                     }
                                } else if (sheet[j][sheet[j].size() - 1] == index[0][index[0].size() - 1]) {
                                     for (int k = (int) index[0].size() - 2; k >= 0; --k) {
                                          sheet[j].push_back(index[0][k]);
                                     }
                                }
                           }
                           for (int k = (int) index[j + 3].size() - 2; k >= 0; --k) {
                                sheet[j].push_back(index[j + 3][k]);
                           }
                      }
                 } else if (index[0].size() == 1 && index[2].size() == 1 && (int) index[1].size() == (numstr - 2)) {
                      _find_strands(index[0][0], strand, data);
                      sheet.push_back(data);
                      strand[index[2][0]].flag = 0;
                      _find_strands(index[2][0], strand, index[3]);
                      for (int k = (int) index[3].size() - 1; k >= 0; --k) {
                           sheet[0].push_back(index[3][k]);
                      }
                 } else if (index[0].size() == 2 && index[2].size() == 2 && (int) index[1].size() == (numstr - 4)) {
                      _find_strands(index[0][0], strand, data);
                      sheet.push_back(data);
                      sheet.push_back(data);
                      for (unsigned int k = 0; k < index[2].size(); ++k) { 
                           strand[index[2][k]].flag = 0;
                      }
                      _find_strands(index[0][1], strand, index[1]);
                      for (unsigned int k = 0; k < index[2].size(); ++k) { 
                           strand[index[2][k]].flag = 0;
                      }
                      int k = sheet[0][sheet[0].size() - 1];
                      _find_strands(k, strand, index[3]);
                      k = index[3][index[3].size() - 1];
                      strand[k].flag = 0;
                      _find_strands(k, strand, index[4]);
                      for (k = 1; k < (int) index[3].size(); ++k) {
                           sheet[0].push_back(index[3][k]);
                      }
                      for (k = (int) index[4].size() - 1; k >= 0; --k) {
                           sheet[1].push_back(index[4][k]);
                      }
                      for (k = (int) index[1].size() - 2; k >= 0; --k) {
                           for (unsigned int j = 0; j < sheet.size(); ++j) {
                                sheet[j].push_back(index[1][k]);
                           }
                      }
                 } else if (index[3].size() == 1 &&
                           (int) index[1].size() == (numstr - 1)) {
                      for (int j = 0; j < 2; ++j) {
                           _find_strands(index[3][0], strand, index[j]);
                           strand[index[3][0]].flag = 0;
                      }
                      sheet.push_back(index[0]);
                      for (unsigned int k = 0; k < index[1].size(); ++k) {
                           sheet[0].push_back(index[1][k]);
                      }
                      sheet[0].push_back(index[3][0]);
                 } else if (index[3].size() == 1 && index[0].size() == 2 && (int) index[1].size() == (numstr - 3)) {
                      for (int j = 0; j < 2; ++j) {
                           _find_strands(index[0][j], strand, index[j + 1]);
                           strand[index[3][0]].flag = 0;
                      }
                      _find_strands(index[3][0], strand, index[4]);
                      sheet.push_back(index[0]);
                      for (unsigned int k = 1; k < index[4].size(); ++k) {
                           sheet[0].push_back(index[4][k]);
                      }
                      for (int k = (int) index[2].size() - 1; k >= 0; --k) {
                           sheet[0].push_back(index[2][k]);
                      }
                 } else if (index[2].size() == 2 && (int) index[1].size() == (numstr - 2)) {
                      for (int j = 0; j < 2; ++j) {
                           strand[index[2][0]].flag = 0;
                           _find_strands(index[2][0], strand, index[j]);
                      }
                      if (index[2].size() > 1 && !index[0].empty() &&
                          index[2][1] == index[0][index[0].size() - 1]) {
                           sheet.push_back(index[1]);
                           for (unsigned k = 0; k < index[0].size(); ++k) {
                                sheet[0].push_back(index[0][k]);
                           }
                      } else if (index[2].size() > 1 && !index[1].empty() && index[2][1] == index[1][index[1].size() - 1]) {
                           sheet.push_back(index[0]);
                           for (unsigned int k = 0; k < index[1].size(); ++k) {
                                sheet[0].push_back(index[1][k]);
                           }
                      }
                      if (!sheet.empty()) {
                           strand[index[2][1]].flag = 0;
                           _find_strands(index[2][1], strand, index[3]);
                           for (int k = (int) index[3].size() - 1; k >= 0; --k) {
                                sheet[0].push_back(index[3][k]);
                           }
                      }
                 }
            }

            index[4].clear();
            
            for (int j = 0; j < numstr; ++j) index[4].push_back(0);
            bool found_incorrect_index = false;
            for (unsigned int j = 0; j < sheet.size(); ++j) {
                 for (unsigned int k = 0; k < sheet[j].size(); ++k) {
                      if (sheet[j][k] >= numstr) {
                           found_incorrect_index = true;
                           exist_incorrect_index = true;
                           _logIo->messageError("Incorrect strand number " + String::IntToString(sheet[j][k] + 1));
                           continue;
                      }
                      index[4][sheet[j][k]] = 1;
                 }
            }
            if (found_incorrect_index) {
                 index.clear();
                 sheet.clear();
                 continue;
            }


            bool complicated_sheet = false;
            for (int j = 0; j < numstr; ++j) {
                 if (!index[4][j]) {
                      complicated_sheet = true;
                      break;
                 }
            }

            if (complicated_sheet || sheet.size() > 1) {
                 _a_sheet.mol_index  = _mol->index();
                 _a_sheet.sheetID    = String::IntToString(_complicated_sheets.size());
                 _a_sheet.numStrands = numstr;
                 _a_sheet.complicateFlag = true;
                 _a_sheet._strands.clear();
                 _a_sheet._strands.reserve(numstr);
                 _a_sheet._strand_orders.clear();
                 _a_sheet._strand_orders.reserve(numstr * (numstr + 1) / 2 + numstr);
                 for (int j = 0; j < numstr; ++j) {
                      _a_strand.strand_id = j + 1;
                      _a_strand.begin_res_index = strand[j].init->index();
                      _a_strand.end_res_index = strand[j].term->index();
                      _a_strand.sense = 0;
                      _a_strand.curr_hbond_res_index = -1;
                      _a_strand.curr_hbond_atom_name.clear();
                      _a_strand.prev_hbond_res_index = -1;
                      _a_strand.prev_hbond_atom_name.clear();
                      _a_sheet._strands.push_back(_a_strand);

                      for (int k = 0; k < strand[j].num_of_partners; ++k) {
                           if (strand[j].partners[k] < j) continue;

                           _a_order.sense_number = strand[j].type[k];
                           if (_a_order.sense_number == 1)
                                _a_order.sense_string = "parallel";
                           else if (_a_order.sense_number == (-1))
                                _a_order.sense_string = "anti-parallel";
                           _a_order.range_id_1 = String::IntToString(j + 1);
                           _a_order.range_id_1_res_index = strand[j].cur_res[k]->index();  // cur/prv are reversed
                           _a_order.range_id_1_atom_name = strand[j].cur_atom[k];
                           _a_order.range_id_2 = String::IntToString(strand[j].partners[k] + 1);
                           _a_order.range_id_2_res_index = strand[j].prv_res[k]->index();
                           _a_order.range_id_2_atom_name = strand[j].prv_atom[k];
                           _a_sheet._strand_orders.push_back(_a_order);
                      }
                 }
                 _complicated_sheets.push_back(_a_sheet);
            }

            if (complicated_sheet) {
                 message += "\nFor SHEET  " + String::IntToString(i + 1) + "\n";
                 message += "\tIndex_of_Strand Num_of_Partners Index_of_Partners\n";
                 for (int j = 0; j < numstr; ++j) {
                      message += "\t" + FloatToString(j + 1, 5, 0) + "               " + String::IntToString(strand[j].num_of_partners) + "       ";
                      for (int k = 0; k < strand[j].num_of_partners; ++k) {
                           message += " " + FloatToString(strand[j].partners[k]+1, 5, 0);
                      }
                      message += "\n";
                 }

                 _complicated_sheet_flag = true;
                 index.clear();
                 sheet.clear();
                 continue;
            }

            for (unsigned int j = 0; j < sheet.size(); ++j) {
                 _a_sheet.mol_index  = _mol->index();
                 _a_sheet.sheetID    = String::IntToString(_sheets.size());
                 _a_sheet.numStrands = sheet[j].size();
                 _a_sheet.complicateFlag = false;
                 _a_sheet._strands.clear();
                 _a_sheet._strands.reserve(sheet[j].size());
                 _a_sheet._strand_orders.clear();
                 _a_sheet._strand_orders.reserve(sheet[j].size());
                 for (unsigned int k = 0; k < sheet[j].size(); ++k) {
                      _a_strand.strand_id = k + 1;
                      _a_strand.begin_res_index = strand[sheet[j][k]].init->index();
                      _a_strand.end_res_index = strand[sheet[j][k]].term->index();
                      _a_strand.sense = 0;
                      _a_strand.curr_hbond_res_index = -1;
                      _a_strand.curr_hbond_atom_name.clear();
                      _a_strand.prev_hbond_res_index = -1;
                      _a_strand.prev_hbond_atom_name.clear();
                      if (k == 0) {
                           _a_sheet._strands.push_back(_a_strand);
                           continue;
                      }

                      _a_order.sense_number = 0;
                      _a_order.sense_string.clear();
                      for (int l = 0; l < strand[sheet[j][k]].num_of_partners; ++l) {
                           if (strand[sheet[j][k]].partners[l] != sheet[j][k-1]) continue;
                           _a_strand.curr_hbond_res_index = strand[sheet[j][k]].cur_res[l]->index();
                           _a_strand.curr_hbond_atom_name = strand[sheet[j][k]].cur_atom[l];
                           _a_strand.prev_hbond_res_index = strand[sheet[j][k]].prv_res[l]->index();
                           _a_strand.prev_hbond_atom_name = strand[sheet[j][k]].prv_atom[l];
                           _a_strand.sense = strand[sheet[j][k]].type[l];

                           _a_order.sense_number = _a_strand.sense;
                           if (_a_order.sense_number == 1)
                                _a_order.sense_string = "parallel";
                           else if (_a_order.sense_number == (-1))
                                _a_order.sense_string = "anti-parallel";
                           _a_order.range_id_1 = String::IntToString(k);
                           _a_order.range_id_1_res_index = strand[sheet[j][k]].prv_res[l]->index();
                           _a_order.range_id_1_atom_name = strand[sheet[j][k]].prv_atom[l];
                           _a_order.range_id_2 = String::IntToString(k + 1);
                           _a_order.range_id_2_res_index = strand[sheet[j][k]].cur_res[l]->index();
                           _a_order.range_id_2_atom_name = strand[sheet[j][k]].cur_atom[l];
                           break;
                      }
                      _a_sheet._strands.push_back(_a_strand);
                      if (!_a_order.sense_string.empty()) _a_sheet._strand_orders.push_back(_a_order);
                 }
                 _sheets.push_back(_a_sheet);
                 if (sheet.size() == 1) {
                     _a_sheet.sheetID = String::IntToString(_complicated_sheets.size());
                     _complicated_sheets.push_back(_a_sheet);
                 }
            }
            index.clear();
            sheet.clear();
       }

       if (_complicated_sheet_flag) {
            std::string message_head = "Warning: There is a complicated topology beta sheet. The SHEET records\n";
            message_head += "        that are automatically created by MAXIT may be incorrect.\n";
            complicated_sheet_warning = "Warning: Maxit failed to produced beta strands/sheets records and no sheet\n";
            complicated_sheet_warning += "        records were found in the entry.\n" + message;
            // _logIo->messageWarning(message_head + message);
            _sheets = _complicated_sheets;
       }

       if (exist_incorrect_index) _sheets.clear();

       strand.clear();
       strks.clear();
       shtks.clear();
       pshtlt.clear();
       _c_sheet.clear();

       if (_sheets.empty()) return;

       std::string character_1, character_2;
       if (_sheets.size() > 6084) {
            character_1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";
            character_2 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ123456789";
       } else {
            character_1 = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
            character_2 = "123456789";
       } 

       std::map<std::string, std::string> id_map;
       id_map.clear();

       std::set<std::string> id_set;
       id_set.clear();
       for (unsigned int i = 0; i < _sheets.size(); ++i) {
            std::string id = get_next_id(character_1, character_1, character_2, id_set);
            id_map.insert(std::make_pair(String::IntToString(i), id));
       }

       for (std::list<_SHEET>::iterator pos = _sheets.begin(); pos != _sheets.end(); ++pos) {
            std::map<std::string, std::string>::const_iterator mpos = id_map.find(pos->sheetID);
            if (mpos != id_map.end()) pos->sheetID = mpos->second; 
       }
}

void Promotif::_getUserDefinedSheet(std::map<int, std::vector<std::string> >& _c_sheet)
{
       _c_sheet.clear();
       if (_supportfile.empty()) return;

       FILE *fp = fopen(_supportfile.c_str(), "r");
       if (fp == NULL) return;

       std::string line;
       int sheet_id = 0;
       std::vector<std::string> string_data;
       string_data.clear();
       while (!feof(fp)) {
            get_line_from_file(fp, line);
            if (line.empty()) continue;
            if (line.find("SHEET") != std::string::npos) {
                 if (!string_data.empty()) _c_sheet.insert(std::make_pair(sheet_id, string_data));
                 get_wordarray(string_data, line, " ");
                 sheet_id = atoi(string_data[1].c_str()) - 1;
                 string_data.clear();
            } else string_data.push_back(line);
       }
       fclose (fp);
       if (!string_data.empty()) _c_sheet.insert(std::make_pair(sheet_id, string_data));
}

void Promotif::_find_strands(const int& start_index, std::vector<STRAND>& strand, std::vector<int>& sheet)
{
       sheet.clear();

       int index = start_index;
       while (index >= 0 && index < (int) strand.size()) {
            sheet.push_back(index);
            if (strand[index].num_of_partners > 2 && index != start_index) break;
            strand[index].flag = 1;
            index = _find_next_strand(strand, index);
       }
}

int Promotif::_find_next_strand(const std::vector<STRAND>& strand, const int& index)
{
       if (index >= 0 && index < (int) strand.size()) {
            for (int i = 0; i < strand[index].num_of_partners; ++i) {
                 int j = strand[index].partners[i];
                 if (!strand[j].flag) return j;
            }
       }
       return (-1);
}
