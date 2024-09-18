/*
FILE:     BasePair3DNA.C
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

#include "BasePairInfo.h"
#include "MatrixUtil.h"
#include "VectorUtil.h"

#define  elem(i, j) ((i) < (j)? (j)*((j)-1)/2+(i): (i)*((i)-1)/2+(j))

#define HTWIST0   0.05
#define XBIG      1.0e+18
#define HELIX_CHG 7.5
#define CNUM      8
#define O3P_UPPER 2.5

static int _leontis_westhof_index[2][3][3] = {
        { { 1,  3,  5 },  { 3,  7,  9 },  { 5,  9, 11 } },
        { { 2,  4,  6 },  { 4,  8, 10 },  { 6, 10, 12 } }
};

void BasePairInfo::_re_ordering()
{
       if (_best_basepairs.empty()) return;

       std::vector<std::vector<int> > bp_order, end_list;
       int num_ends = _get_bp_context(bp_order, end_list);

       std::vector<int> bp_idx, data;
       std::vector<std::vector<int> > helix_idx = _locate_helix(num_ends, bp_order, end_list, bp_idx); 
       if (helix_idx.empty()) return;

       _five2three(helix_idx, bp_idx);

       std::vector<std::vector<int> > index;
       for (std::vector<std::vector<int> >::const_iterator pos = helix_idx.begin(); pos != helix_idx.end(); ++pos) {
            index.clear();
            for (int j = (*pos)[0]; j < (*pos)[1]; ++j) {
                 data.clear();
                 data.push_back(_best_basepairs[bp_idx[j]].i);
                 data.push_back(_best_basepairs[bp_idx[j]].j);
                 data.push_back(_best_basepairs[bp_idx[j]].type);
                 data.push_back(_best_basepairs[bp_idx[j]].type_i);
                 data.push_back(_best_basepairs[bp_idx[j]].type_j);
                 data.push_back(_best_basepairs[bp_idx[j]].cis_or_trans);
                 index.push_back(data);
            }
            if (!index.empty()) _get_parameters(index);
       }
}

int BasePairInfo::_get_bp_context(std::vector<std::vector<int> >& bp_order, std::vector<std::vector<int> >& end_list)
{
       bp_order.clear();
       end_list.clear();

       int num_ends = 0;

       std::vector<int> data;
       for (int i = 0; i < 3; ++i) data.push_back(-1);
       for (unsigned int i = 0; i < _best_basepairs.size(); ++i) {
            end_list.push_back(data);
            bp_order.push_back(data);
       }

       COORD txyz, txyz2;
       std::vector<double> dist_matrix;
       dist_matrix.clear();
       for (unsigned int i = 0; i < elem(_best_basepairs.size(), 0); ++i) dist_matrix.push_back(XBIG);
       for (unsigned int i = 0; i < _best_basepairs.size() - 1; ++i) {
            for (unsigned int j = i + 1; j < _best_basepairs.size(); ++j) {
                 vector_difference(txyz, _best_basepairs[i].oave, _best_basepairs[j].oave);
                 dist_matrix[elem(i, j)] = vector_length(txyz);
            }
       }

       std::vector<int> ddidx;
       std::vector<double> ddmin;
       for (unsigned int i = 0; i < _best_basepairs.size(); ++i) {
            _get_idx_and_min(i, 0, _best_basepairs.size(), dist_matrix, ddidx, ddmin);

            if (ddidx[0] >= 0 && ddidx[1] >= 0) {
                 if (ddmin[0] > HELIX_CHG) {
                      end_list[num_ends++][0] = i;
                 } else {
                      vector_difference(txyz, _best_basepairs[i].oave, _best_basepairs[ddidx[0]].oave);
                      double d = vector_dot_product(_best_basepairs[i].zave, txyz);
                      if (ddidx[2] >= 0 && ddmin[1] <= HELIX_CHG && ddmin[2] <= HELIX_CHG) {
                           vector_difference(txyz, _best_basepairs[i].oave, _best_basepairs[ddidx[1]].oave);
                           vector_difference(txyz2, _best_basepairs[i].oave, _best_basepairs[ddidx[2]].oave);
                           double d2 = vector_dot_product(_best_basepairs[i].zave, txyz);
                           double d3 = vector_dot_product(_best_basepairs[i].zave, txyz2);
                           if (d * d2 < 0.0 && d * d3 < 0.0 && fabs(d2) > fabs(d3)) {
                                int k = ddidx[1]; ddidx[1] = ddidx[2]; ddidx[2] = k;
                                double temp1 = ddmin[1]; ddmin[1] = ddmin[2];
                                ddmin[2] = temp1;
                           }
                      }
                      int n = -1;
                      for (unsigned int j = 1; j < ddidx.size() && ddidx[j] >= 0; ++j) {
                           if (ddmin[j] > HELIX_CHG) break;
                           vector_difference(txyz2, _best_basepairs[i].oave, _best_basepairs[ddidx[j]].oave);
                           double d2 = vector_dot_product(_best_basepairs[i].zave, txyz2);
                           if (d * d2 < 0.0) {
                                n = j;
                                bp_order[i][0] = -2;
                                bp_order[i][1] = ddidx[0];
                                bp_order[i][2] = ddidx[j];
                                break;
                           }
                      }
                      if (n < 0) {
                           end_list[num_ends][0] = i;
                           end_list[num_ends][1] = ddidx[0];
                           bp_order[i][1] = ddidx[0];
                           vector_difference(txyz2, _best_basepairs[ddidx[1]].oave, _best_basepairs[ddidx[0]].oave);
                           double d2 = vector_dot_product(_best_basepairs[i].zave, txyz2);
                           if (d * d2 < 0.0 && vector_length(txyz2) <= HELIX_CHG) {
                                end_list[num_ends][2] = ddidx[1];
                                bp_order[i][2] = ddidx[1];
                           }
                           num_ends++;
                      }
                 }
            }
       }

       if (!num_ends && _is_circle_helix(bp_order)) {
            _get_idx_and_min(0, 0, _best_basepairs.size() - 1, dist_matrix, ddidx, ddmin);
            bp_order[0][0] = -1;
            bp_order[0][1] = ddidx[0];
            bp_order[0][2] = ddidx[1];
            end_list[0][0] = 0;
            end_list[0][1] = ddidx[0];
            end_list[0][2] = ddidx[1];
            _get_idx_and_min(_best_basepairs.size() - 1, 1, _best_basepairs.size(), dist_matrix, ddidx, ddmin);
            bp_order[_best_basepairs.size()-1][0] = -1;
            bp_order[_best_basepairs.size()-1][1] = ddidx[0];
            bp_order[_best_basepairs.size()-1][2] = ddidx[1];
            end_list[1][0] = _best_basepairs.size() - 1;
            end_list[1][1] = ddidx[0];
            end_list[1][2] = ddidx[1];
       }

       if (!num_ends) {
            end_list[num_ends][0] = 0;
            if (_best_basepairs.size() == 2) {
                 if (ddmin[0] <= HELIX_CHG) {
                      end_list[num_ends][1] = 1;
                      num_ends++;
                      end_list[num_ends][0] = 1;
                      end_list[num_ends][1] = 0;
                      num_ends++;
                 } else {
                      num_ends++;
                      end_list[num_ends][0] = 1;
                      num_ends++;
                 }
            } else num_ends++;
       }

       return num_ends;
}

std::vector<std::vector<int> > BasePairInfo::_locate_helix(const int& num_ends, const std::vector<std::vector<int> >& bp_order,
                                                      const std::vector<std::vector<int> >& end_list, std::vector<int>& bp_idx)
{
       bp_idx.clear();

       std::vector<std::vector<int> > helix_idx;
       helix_idx.clear();

       int num_pair = _best_basepairs.size();

       std::vector<int> matched_idx; // indicator for used bps
       matched_idx.clear();
       for (int i = 0; i < num_pair; ++i) {
            matched_idx.push_back(0);
            bp_idx.push_back(i);
       }

       std::vector<int> data;
       data.clear();
       data.push_back(0);
       data.push_back(0);

       int ip = 0;
       for (int i = 0; i < num_ends && ip < num_pair; ++i) {
            int k = 0;
            int k0 = 0;
            for (int j = 0; j < 3; ++j) {
                 if (end_list[i][j] >= 0) {
                      k += matched_idx[end_list[i][j]];
                      k0++;
                 }
            }
            if (k == k0) continue; // end point of a processed helix
            for (int j = 0; j < 3 && ip < num_pair; ++j) {
                 if (end_list[i][j] >= 0 && !matched_idx[end_list[i][j]]) {
                      bp_idx[ip++] = end_list[i][j];
                      matched_idx[end_list[i][j]] = 1;
                 }
            }
            for (int j = 0; j < num_pair && ip && ip < num_pair; ++j) {
                 int k = bp_idx[ip - 1];
                 int k2 = bp_order[k][1];
                 int k3 = bp_order[k][2];
                 if (bp_order[k][0] == -1) { // end-point
                      if (k2 >= 0 && !matched_idx[k2] && k3 == -1) {
                           bp_idx[ip++] = k2;
                           matched_idx[k2] = 1;
                      }
                      break; // normal case
                 }
                 int m = matched_idx[k2] + matched_idx[k3];
                 if (m == 2 || m == 0) break; // chain terminates
                 if (k2 == bp_idx[ip - 2]) {
                      bp_idx[ip++] = k3;
                      matched_idx[k3] = 1;
                 } else if (k3 == bp_idx[ip - 2]) {
                      bp_idx[ip++] = k2;
                      matched_idx[k2] = 1;
                 } else break; // no direct connection
            }
            data[1] = ip;
            helix_idx.push_back(data);
            if (ip && ip < num_pair) {
/*
                 if (((int) helix_idx.size()) >= num_ends) {
                      helix_idx.clear();
                      return helix_idx;
                 }
*/
                 data[0] = ip;
            }
       }
/*
       if (ip < num_pair) {
            data[1] = num_pair;
            helix_idx.push_back(data);
            for (int i = 0; i < num_pair; ++i) {
                 if (ip == num_pair - 1) break;
                 if (!matched_idx[i]) bp_idx[++ip] = i;
            }
       }
*/
       return helix_idx;
}

void BasePairInfo::_five2three(const std::vector<std::vector<int> >& helix_idx, std::vector<int>& bp_idx)
{
       for (std::vector<std::vector<int> >::const_iterator pos = helix_idx.begin(); pos != helix_idx.end(); ++pos) {
            for (int j = (*pos)[0]; j < (*pos)[1]; ++j) {
                 double di1_i2 = XBIG;
                 if (_distance_ab(_baseframes[_best_basepairs[bp_idx[j]].i], "P", _baseframes[_best_basepairs[bp_idx[j]].i], "O", di1_i2)) {
                      if (di1_i2 <= O3P_UPPER) return; /* Wrong: O3' connected to P, ignore 5'-->3' re-arrangement */
                 }
            }
       }

       std::vector<int> swapped;
       swapped.clear();
       for (unsigned int i = 0; i < _best_basepairs.size(); ++i) swapped.push_back(0);

       int i1, j1, i2, j2, direction[6];
       for (std::vector<std::vector<int> >::const_iterator pos = helix_idx.begin(); pos != helix_idx.end(); ++pos) {
            if ((*pos)[0] >= ((*pos)[1] - 1)) continue;
            _first_step(*pos, swapped, bp_idx);
            for (int j = (*pos)[0]; j < (*pos)[1] - 1; ++j) {
                 _get_ij(bp_idx[j], swapped[bp_idx[j]], i1, j1);
                 _get_ij(bp_idx[j + 1], swapped[bp_idx[j + 1]], i2, j2);

                 if (_wc_bporien(bp_idx[j], bp_idx[j + 1], i1, j1, i2, j2)) swapped[bp_idx[j + 1]] = !swapped[bp_idx[j + 1]];
                 else if (_check_o3dist(i1, j1, i2, j2)) swapped[bp_idx[j + 1]] = !swapped[bp_idx[j + 1]];
                 else if (!_is_linked(_baseframes[i1], _baseframes[i2]) && !_is_linked(_baseframes[j1], _baseframes[j2]) &&
                         (_is_linked(_baseframes[i1], _baseframes[j2]) || _is_linked(_baseframes[j1], _baseframes[i2])))
                      swapped[bp_idx[j + 1]] = !swapped[bp_idx[j + 1]];
                 else if (_check_others(i1, j1, i2, j2)) swapped[bp_idx[j + 1]] = !swapped[bp_idx[j + 1]];

                 _get_ij(bp_idx[j], swapped[bp_idx[j]], i1, j1);
                 _get_ij(bp_idx[j + 1], swapped[bp_idx[j + 1]], i2, j2);
                 if (_is_linked(_baseframes[i1], _baseframes[i2]) == -1) swapped[bp_idx[j + 1]] = !swapped[bp_idx[j + 1]];

                 _get_ij(bp_idx[j], swapped[bp_idx[j]], i1, j1);
                 _get_ij(bp_idx[j + 1], swapped[bp_idx[j + 1]], i2, j2);
                 if (_wc_bporien(bp_idx[j], bp_idx[j + 1], i1, j1, i2, j2)) swapped[bp_idx[j + 1]] = !swapped[bp_idx[j + 1]];
            }

            if (_check_direction(*pos, bp_idx, swapped, direction)) {
                 if (direction[0] + direction[1] + direction[3] + direction[4]) {
                      for (int j = (*pos)[0]; j < (*pos)[1] - 1; ++j) {
                           _get_ij(bp_idx[j], swapped[bp_idx[j]], i1, j1);
                           _get_ij(bp_idx[j + 1], swapped[bp_idx[j + 1]], i2, j2);
                           if (_wc_bporien(bp_idx[j], bp_idx[j + 1], i1, j1, i2, j2)) continue;
                           if (!_is_linked(_baseframes[i1], _baseframes[i2]) && !_is_linked(_baseframes[j1], _baseframes[j2]) &&
                              ((_is_linked(_baseframes[i1], _baseframes[j2]) == 1) || (_is_linked(_baseframes[i1], _baseframes[j2]) &&
                                _is_linked(_baseframes[j1], _baseframes[i2])))) swapped[bp_idx[j + 1]] = !swapped[bp_idx[j + 1]];
                      }
                 }
            } else {
                 int anti_parallel = direction[0] > direction[1] && direction[3] < direction[4];
                 int parallel      = direction[0] > direction[1] && direction[3] > direction[4];
                 for (int j = (*pos)[0]; j < (*pos)[1] - 1; ++j) {
                      _get_ij(bp_idx[j], swapped[bp_idx[j]], i1, j1);
                      _get_ij(bp_idx[j + 1], swapped[bp_idx[j + 1]], i2, j2);

                      int k = _is_linked(_baseframes[j1], _baseframes[j2]);
                      if (!_is_linked(_baseframes[i1], _baseframes[i2]) && ((anti_parallel && k == 1) || (parallel && (k == -1))))
                           swapped[bp_idx[j + 1]] = !swapped[bp_idx[j + 1]];

                      _get_ij(bp_idx[j + 1], swapped[bp_idx[j + 1]], i2, j2);
                      if (!_is_linked(_baseframes[i1], _baseframes[i2]) && !_is_linked(_baseframes[j1], _baseframes[j2])) {
                           if ((anti_parallel && _is_linked(_baseframes[j1], _baseframes[i2]) == 1) ||
                               (parallel && _is_linked(_baseframes[i1], _baseframes[j2]) == -1))
                                swapped[bp_idx[j]] = !swapped[bp_idx[j]];
                           else if ((parallel && _is_linked(_baseframes[j1], _baseframes[i2]) == -1) ||
                               (anti_parallel && _is_linked(_baseframes[i1], _baseframes[j2]) == 1))
                                swapped[bp_idx[j + 1]] = !swapped[bp_idx[j + 1]];
                      }
                 }
            }
            _check_direction(*pos, bp_idx, swapped, direction);

            for (int j = (*pos)[0]; j < (*pos)[1]; ++j) {
                 if (!swapped[bp_idx[j]]) continue;

                 i1 = _best_basepairs[bp_idx[j]].i;
                 _best_basepairs[bp_idx[j]].i = _best_basepairs[bp_idx[j]].j;
                 _best_basepairs[bp_idx[j]].j = i1;

                 i1 = _best_basepairs[bp_idx[j]].type_i;
                 _best_basepairs[bp_idx[j]].type_i = _best_basepairs[bp_idx[j]].type_j;
                 _best_basepairs[bp_idx[j]].type_j = i1;

                 _best_basepairs[bp_idx[j]].zave.x = -_best_basepairs[bp_idx[j]].zave.x;
                 _best_basepairs[bp_idx[j]].zave.y = -_best_basepairs[bp_idx[j]].zave.y;
                 _best_basepairs[bp_idx[j]].zave.z = -_best_basepairs[bp_idx[j]].zave.z;
            }
       }
}

void BasePairInfo::_get_idx_and_min(const int& i, const int& start, const int& end, const std::vector<double>& dist_matrix,
                                    std::vector<int>& ddidx, std::vector<double>& ddmin)
{
       ddidx.clear();
       ddmin.clear();

       std::multimap<double, int> dist_map;
       dist_map.clear();
       for (int j = start; j < end; ++j) {
            if (j == i) continue;
            dist_map.insert(std::make_pair(dist_matrix[elem(i, j)], j));
       }

       std::vector<std::pair<int, double> > dist_pair;
       dist_pair.clear();
       for (std::multimap<double, int>::const_iterator mpos = dist_map.begin(); mpos != dist_map.end(); ++mpos) {
            ddidx.push_back(mpos->second);
            ddmin.push_back(mpos->first);
            if (ddmin.size() >= CNUM) break;
       }
       for (unsigned int i = ddidx.size(); i < 3; ++i) {
            ddidx.push_back(-1);
            ddmin.push_back(XBIG);
       }
}

bool BasePairInfo::_is_circle_helix(std::vector<std::vector<int> >& bp_order)
{
       std::vector<int> circle;
       circle.clear();
       for (int i = 0; i < (int) bp_order.size(); ++i) {
            if (bp_order[i][0] != -2) return false;
            if (circle.empty()) {
                 circle.push_back(bp_order[i][1]);
                 circle.push_back(i);
                 circle.push_back(bp_order[i][2]);
            } else if (circle[0] == i) {
                 if (circle[1] == bp_order[i][1]) {
                      if (circle[circle.size() - 1] == bp_order[i][2]) break;
                      else circle.insert(circle.begin(), bp_order[i][2]);
                 } else if (circle[1] == bp_order[i][2]) {
                      if (circle[circle.size() - 1] == bp_order[i][1]) break;
                      else circle.insert(circle.begin(), bp_order[i][1]);
                 } else return false;
            } else if (circle[circle.size() - 1] == i) {
                 if (circle[circle.size() - 2] == bp_order[i][1]) {
                      if (circle[0] == bp_order[i][2]) break;
                      else circle.push_back(bp_order[i][2]);
                 } else if (circle[circle.size() - 2] == bp_order[i][2]) {
                      if (circle[0] == bp_order[i][1]) break;
                      else circle.push_back(bp_order[i][1]);
                 } else return false;
            } else return false; 
       }
       if (circle.size() != bp_order.size()) return false;
       return true;
}

bool BasePairInfo::_distance_ab(const _BASE& a, const std::string& ia, const _BASE& b, const std::string& ib, double& dist)
{
       std::map<std::string, COORD>::const_iterator
           mpos_a = a.coord_xyz.find(ia);
       if (mpos_a == a.coord_xyz.end()) return false;
       std::map<std::string, COORD>::const_iterator
           mpos_b = b.coord_xyz.find(ib);
       if (mpos_b == b.coord_xyz.end()) return false;

       COORD pc;
       vector_difference(pc, mpos_a->second, mpos_b->second);
       dist = vector_length(pc);
       return true;
}

void BasePairInfo::_first_step(const std::vector<int>& helix_idx, std::vector<int>& swapped, std::vector<int>& bp_idx)
/* make strand I of the first step in 5'--->3' direction */
{
       int i1, j1, i2, j2;
       _get_ij(bp_idx[helix_idx[0]], swapped[bp_idx[helix_idx[0]]], i1, j1);
       _get_ij(bp_idx[helix_idx[0] + 1], swapped[bp_idx[helix_idx[0] + 1]], i2, j2);
       int k = _is_linked(_baseframes[i1], _baseframes[i2]);
       if (k == -1) swapped[bp_idx[helix_idx[0]]] = !swapped[bp_idx[helix_idx[0]]];
       else if (!k) {
            _reverse(helix_idx[0], helix_idx[1] - helix_idx[0], bp_idx);
            _get_ij(bp_idx[helix_idx[0]], swapped[bp_idx[helix_idx[0]]], i1, j1);
            _get_ij(bp_idx[helix_idx[0] + 1], swapped[bp_idx[helix_idx[0] + 1]], i2, j2);
            k = _is_linked(_baseframes[i1], _baseframes[i2]);
           if (k == -1) swapped[bp_idx[helix_idx[0]]] = !swapped[bp_idx[helix_idx[0]]];
            else if (!k) _reverse(helix_idx[0], helix_idx[1] - helix_idx[0], bp_idx);
       }
}

void BasePairInfo::_get_ij(const int& m, const int& swapped, int &i, int &j)
{
       i = _best_basepairs[m].i;
       j = _best_basepairs[m].j;
       if (swapped) {
            i = _best_basepairs[m].j;
            j = _best_basepairs[m].i;
       }
}

bool BasePairInfo::_wc_bporien(const int& fst, const int& snd, const int& i1, const int& j1, const int& i2, const int& j2)
{
       if (_is_wc_geometry(fst) && _is_wc_geometry(snd)) {
            COORD txyz, txyz2;
            vector_sum(txyz, _baseframes[i1].orien[0], _baseframes[j1].orien[0]);
            vector_sum(txyz2, _baseframes[i2].orien[0], _baseframes[j2].orien[0]);
            if (magang(txyz, txyz2) > 125.0 || _is_linked(_baseframes[i1], _baseframes[i2]) || _is_linked(_baseframes[j1], _baseframes[j2])) return false;

            vector_difference(txyz, _baseframes[i1].orien[2], _baseframes[j1].orien[2]);
            vector_difference(txyz2, _baseframes[i2].orien[2], _baseframes[j2].orien[2]);
            vector_normalize(txyz);
            vector_normalize(txyz2);
            double di1_i2 = vector_dot_product(txyz, txyz2);
            vector_difference(txyz2, _baseframes[j2].orien[2], _baseframes[i2].orien[2]);
            vector_normalize(txyz2);
            double di1_j2 = vector_dot_product(txyz, txyz2);
            if (di1_i2 < 0.0 && di1_j2 > 0.0) return true;
       }
       return false;
}

bool BasePairInfo::_check_o3dist(const int& i1, const int& j1, const int& i2, const int& j2)
{
       double di1_i2 = XBIG;
       double di1_j2 = XBIG;
       double dj1_i2 = XBIG;
       double dj1_j2 = XBIG;
       if (_distance_ab(_baseframes[i1], "O", _baseframes[i2], "O", di1_i2) && _distance_ab(_baseframes[i1], "O", _baseframes[j2], "O", di1_j2) &&
           _distance_ab(_baseframes[j1], "O", _baseframes[i2], "O", dj1_i2) && _distance_ab(_baseframes[j1], "O", _baseframes[j2], "O", dj1_j2)) {
            if (di1_i2 > di1_j2 && dj1_j2 > dj1_i2) return true;
       }
       return false;
}

void BasePairInfo::_reverse(const int& st, const int& n, std::vector<int>& vec)
{
       std::vector<int> vtmp;
       vtmp.clear();
       for (int i = 0; i < n; ++i) {
            vtmp.push_back(vec[n-1+st-i]);
       }
       for (int i = 0; i < n; ++i) {
            vec[st + i] = vtmp[i];
       }
}

int BasePairInfo::_is_linked(const _BASE& base1, const _BASE& base2)
{
       double dist = XBIG;
       if (_distance_ab(base1, "O", base2, "P", dist)) {
            if (dist <= O3P_UPPER) return 1;
       }
       if (_distance_ab(base1, "P", base2, "O", dist)) {
            if (dist <= O3P_UPPER) return (-1);
       }
       return 0;
}

bool BasePairInfo::_is_wc_geometry(const int& bs_idx)
{
       return (_best_basepairs[bs_idx].dir.x > 0.0 && _best_basepairs[bs_idx].dir.y < 0.0 && _best_basepairs[bs_idx].dir.z < 0.0);
}

int BasePairInfo::_check_others(const int& i1, const int& j1, const int& i2, const int& j2)
{
       double a1[3], a2[3], r1[3], r2[3], d[4];

       if (_is_linked(_baseframes[i1], _baseframes[i2]) || _is_linked(_baseframes[j1], _baseframes[j2]) ||
           _is_linked(_baseframes[i1], _baseframes[j2]) || _is_linked(_baseframes[j1], _baseframes[i2]))
            return 0;

       a1[0] = vector_dot_product(_baseframes[i1].orien[0], _baseframes[i2].orien[0]);
       a1[1] = vector_dot_product(_baseframes[i1].orien[1], _baseframes[i2].orien[1]);
       a1[2] = vector_dot_product(_baseframes[i1].orien[2], _baseframes[i2].orien[2]);
       a2[0] = vector_dot_product(_baseframes[j1].orien[0], _baseframes[j2].orien[0]);
       a2[1] = vector_dot_product(_baseframes[j1].orien[1], _baseframes[j2].orien[1]);
       a2[2] = vector_dot_product(_baseframes[j1].orien[2], _baseframes[j2].orien[2]);
       int ii1 = a1[0] > 0 && a1[1] > 0 && a1[2] > 0;
       int ii2 = a2[0] > 0 && a2[1] > 0 && a2[2] > 0;
       if (ii1 && ii2) return 0;

       r1[0] = vector_dot_product(_baseframes[i1].orien[0], _baseframes[j2].orien[0]);
       r1[1] = vector_dot_product(_baseframes[i1].orien[1], _baseframes[j2].orien[1]);
       r1[2] = vector_dot_product(_baseframes[i1].orien[2], _baseframes[j2].orien[2]);
       r2[0] = vector_dot_product(_baseframes[j1].orien[0], _baseframes[i2].orien[0]);
       r2[1] = vector_dot_product(_baseframes[j1].orien[1], _baseframes[i2].orien[1]);
       r2[2] = vector_dot_product(_baseframes[j1].orien[2], _baseframes[i2].orien[2]);
       int jj1 = r1[0] > 0 && r1[1] > 0 && r1[2] > 0;
       int jj2 = r2[0] > 0 && r2[1] > 0 && r2[2] > 0;
       if (!ii1 && !ii2) {
            if (jj1 || jj2) return 1;
            else if (!jj1 && !jj2) return 0;
       }

       d[0] = dot2ang(a1[0]) + dot2ang(a1[1]) + dot2ang(a1[2]);
       d[1] = dot2ang(a2[0]) + dot2ang(a2[1]) + dot2ang(a2[2]);
       d[2] = dot2ang(r1[0]) + dot2ang(r1[1]) + dot2ang(r1[2]);
       d[3] = dot2ang(r2[0]) + dot2ang(r2[1]) + dot2ang(r2[2]);
       if (ii1 && jj1)
            return (d[0] > d[2]) ? 1 : 0;
       if (ii1 && jj2)
            return (d[0] > d[3]) ? 1 : 0;
       if (ii2 && jj1)
            return (d[1] > d[2]) ? 1 : 0;
       if (ii2 && jj2)
            return (d[1] > d[3]) ? 1 : 0;
       return 0;
}

bool BasePairInfo::_check_direction(const std::vector<int>& helix_idx, std::vector<int>& bp_idx, std::vector<int>& swapped, int *direction)
{
       int i1, j1, i2, j2;
       for (int j = 0; j < 6; ++j) direction[j] = 0;
       for (int j = helix_idx[0]; j < helix_idx[1] - 1; ++j) {
            int m = bp_idx[j];
            int n = bp_idx[j+1];
            _get_ij(m, swapped[m], i1, j1);
            _get_ij(n, swapped[n], i2, j2);
            int k = _is_linked(_baseframes[i1], _baseframes[i2]);
            (k == 1) ? ++direction[0] : (k == -1) ? ++direction[1] : ++direction[2];
            k = _is_linked(_baseframes[j1], _baseframes[j2]);
            (k == 1) ? ++direction[3] : (k == -1) ? ++direction[4] : ++direction[5];
       }
       if ((direction[0] && direction[1]) || (direction[3] && direction[4])) return false;

       if ((direction[0] + direction[1] + direction[3] + direction[4]) == 0) return true;

       int m = bp_idx[helix_idx[0]];
       int n = bp_idx[helix_idx[1] - 1];
       _get_ij(m, swapped[m], i1, j1);
       _get_ij(n, swapped[n], i2, j2);
       if (direction[0] && !direction[1]) {
            if (!direction[3] && direction[4]) {
                 if (i1 > j2) {
                      for (int j = helix_idx[0]; j < helix_idx[1]; ++j) swapped[bp_idx[j]] = !swapped[bp_idx[j]];
                      _reverse(helix_idx[0], helix_idx[1] - helix_idx[0], bp_idx);
                 }
            } else if (direction[3] && !direction[4]) {
                 if (i1 > j1) {
                      for (int j = helix_idx[0]; j < helix_idx[1]; ++j) swapped[bp_idx[j]] = !swapped[bp_idx[j]];
                 }
            }
       }
       return true;
}

void BasePairInfo::_get_parameters(const std::vector<std::vector<int> >& pair_num)
{
       int num_bp = pair_num.size();

       COORD mfoi_prev, mst_org;
       mfoi_prev.x = 0;
       mfoi_prev.y = 0;
       mfoi_prev.z = 0;

       for (int j = 0; j < num_bp; ++j) {
            if (vector_dot_product(_baseframes[pair_num[j][0]].orien[2], _baseframes[pair_num[j][1]].orien[2]) < 0.0) {
                 vector_reverse(_baseframes[pair_num[j][1]].orien[1]);
                 vector_reverse(_baseframes[pair_num[j][1]].orien[2]);
            }
       }

       int xdir = 0, ds = 2;
       mst_org.x = mst_org.y = mst_org.z = 0;
       for (int i = 0; i < ds; ++i) {
            for (int j = 0; j < num_bp; ++j) {
                 if (j < (num_bp - 1)) vector_difference(mst_org, _baseframes[pair_num[j + 1][i]].org, _baseframes[pair_num[j][i]].org);
                 if (vector_dot_product(mst_org, _baseframes[pair_num[j][i]].orien[2]) < 0.0) xdir++;
            }
       }
       if (xdir == ds * num_bp) {
            for (int i = 0; i < ds; ++i) {
                 for (int j = 0; j < num_bp; ++j) {
                      vector_reverse(_baseframes[pair_num[j][i]].orien[0]);
                      vector_reverse(_baseframes[pair_num[j][i]].orien[2]);
                 }
            }
       }

       COORD mfi[3], mfoi, mfi_prev[3];
       COORD mst_orien[3], mst_orienH[3], mst_orgH;
       double spi[6], hpi[6], dsum = 0;

       std::vector<double> twist;
       twist.clear();
       std::vector<COORD> aveS, aveH;
       aveS.clear();
       aveH.clear();
       vector_init(mfoi, XBIG, XBIG, XBIG);
       for (int i = 0; i < num_bp - 1; ++i) {
            aveS.push_back(mfoi);
            aveH.push_back(mfoi);    
       }

       _BASEBASE_PARAMS bb_params;
       _INTERBASE_PARAMS ib_params;
       for (int i = 0; i < num_bp; ++i) {
            _bpstep_par(_baseframes[pair_num[i][1]].orien, _baseframes[pair_num[i][1]].org, _baseframes[pair_num[i][0]].orien,
                        _baseframes[pair_num[i][0]].org, spi, mfi, mfoi);

            bb_params.I_index = _baseframes[pair_num[i][0]].res->index();
            bb_params.J_index = _baseframes[pair_num[i][1]].res->index();
            bb_params.Sym_I   = "1_555";
            bb_params.Sym_J   = "1_555";
            if (_baseframes[pair_num[i][0]].op || _baseframes[pair_num[i][0]].lx || _baseframes[pair_num[i][0]].ly || _baseframes[pair_num[i][0]].lz) {
                 bb_params.Sym_I = String::IntToString(_baseframes[pair_num[i][0]].op+1) + "_" + String::IntToString(_baseframes[pair_num[i][0]].lx+5)
                                 + String::IntToString(_baseframes[pair_num[i][0]].ly+5) + String::IntToString(_baseframes[pair_num[i][0]].lz+5);
            }
            if (_baseframes[pair_num[i][1]].op || _baseframes[pair_num[i][1]].lx || _baseframes[pair_num[i][1]].ly || _baseframes[pair_num[i][1]].lz) {
                 bb_params.Sym_J = String::IntToString(_baseframes[pair_num[i][1]].op+1) + "_" + String::IntToString(_baseframes[pair_num[i][1]].lx+5)
                                 + String::IntToString(_baseframes[pair_num[i][1]].ly+5) + String::IntToString(_baseframes[pair_num[i][1]].lz+5);
            }
            bb_params.Shear   = spi[0];
            bb_params.Stretch = spi[1];
            bb_params.Stagger = spi[2];
            bb_params.Buckle  = spi[3];
            bb_params.Propel  = spi[4];
            bb_params.Opening = spi[5];
            bb_params._saenger_classification.clear();
            if (pair_num[i][2] > 0 && pair_num[i][2] < 30)
                 bb_params._saenger_classification = String::IntToString(pair_num[i][2]);
            bb_params._leontis_westhof_classification.clear();
            if (pair_num[i][3] && pair_num[i][4] && pair_num[i][5]) {
                 bb_params._leontis_westhof_classification = String::IntToString(_leontis_westhof_index[pair_num[i][5]-1][pair_num[i][3]-1][pair_num[i][4]-1]);
            } else if (pair_num[i][2] == 19 || pair_num[i][2] == 20)
                 bb_params._leontis_westhof_classification = "1";
            _basebase_params.push_back(bb_params);

            if (i == 0) {
                 mfoi_prev = mfoi;
                 for (int j = 0; j < 3; ++j) mfi_prev[j] = mfi[j];
                 continue;
            }

            _bpstep_par(mfi_prev, mfoi_prev, mfi, mfoi, spi, mst_orien, mst_org);
            _helical_par(mfi_prev, mfoi_prev, mfi, mfoi, hpi, mst_orienH, mst_orgH);

            mfoi_prev = mfoi;
            for (int j = 0; j < 3; ++j) mfi_prev[j] = mfi[j];

            std::map<std::string, COORD>::const_iterator mpos0 = _baseframes[pair_num[i][0]].coord_xyz.find("P");
            std::map<std::string, COORD>::const_iterator mpos1 = _baseframes[pair_num[i-1][1]].coord_xyz.find("P");
            if (mpos0 != _baseframes[pair_num[i][0]].coord_xyz.end() && mpos1 != _baseframes[pair_num[i-1][1]].coord_xyz.end()) {
                _project_xyzP(mpos0->second, mpos1->second, mst_orien, mst_org, aveS[i-1]);
                _project_xyzP(mpos0->second, mpos1->second, mst_orienH, mst_orgH, aveH[i-1]);
            }

            dsum += spi[5];
            twist.push_back(spi[5]);

            ib_params.I1_index = _baseframes[pair_num[i-1][0]].res->index();
            ib_params.Sym_I1   = "1_555";
            if (_baseframes[pair_num[i-1][0]].op || _baseframes[pair_num[i-1][0]].lx || _baseframes[pair_num[i-1][0]].ly || _baseframes[pair_num[i-1][0]].lz) {
                 ib_params.Sym_I1 = String::IntToString(_baseframes[pair_num[i-1][0]].op+1) + "_" + String::IntToString(_baseframes[pair_num[i-1][0]].lx+5)
                                  + String::IntToString(_baseframes[pair_num[i-1][0]].ly+5) + String::IntToString(_baseframes[pair_num[i-1][0]].lz+5);
            }
            ib_params.J1_index = _baseframes[pair_num[i-1][1]].res->index();
            ib_params.Sym_J1   = "1_555";
            if (_baseframes[pair_num[i-1][1]].op || _baseframes[pair_num[i-1][1]].lx || _baseframes[pair_num[i-1][1]].ly || _baseframes[pair_num[i-1][1]].lz) {
                 ib_params.Sym_J1 = String::IntToString(_baseframes[pair_num[i-1][1]].op+1) + "_" + String::IntToString(_baseframes[pair_num[i-1][1]].lx+5)
                                  + String::IntToString(_baseframes[pair_num[i-1][1]].ly+5) + String::IntToString(_baseframes[pair_num[i-1][1]].lz+5);
            }
            ib_params.I2_index = _baseframes[pair_num[i][0]].res->index();
            ib_params.Sym_I2   = "1_555";
            if (_baseframes[pair_num[i][0]].op || _baseframes[pair_num[i][0]].lx || _baseframes[pair_num[i][0]].ly || _baseframes[pair_num[i][0]].lz) {
                 ib_params.Sym_I2 = String::IntToString(_baseframes[pair_num[i][0]].op+1) + "_" + String::IntToString(_baseframes[pair_num[i][0]].lx+5)
                                  + String::IntToString(_baseframes[pair_num[i][0]].ly+5) + String::IntToString(_baseframes[pair_num[i][0]].lz+5);
            }
            ib_params.J2_index = _baseframes[pair_num[i][1]].res->index();
            ib_params.Sym_J2   = "1_555";
            if (_baseframes[pair_num[i][1]].op || _baseframes[pair_num[i][1]].lx || _baseframes[pair_num[i][1]].ly || _baseframes[pair_num[i][1]].lz) {
                 ib_params.Sym_J2 = String::IntToString(_baseframes[pair_num[i][1]].op+1) + "_" + String::IntToString(_baseframes[pair_num[i][1]].lx+5)
                                  + String::IntToString(_baseframes[pair_num[i][1]].ly+5) + String::IntToString(_baseframes[pair_num[i][1]].lz+5);
            }
            ib_params.Shift    = spi[0];
            ib_params.Slide    = spi[1];
            ib_params.Rise     = spi[2];
            ib_params.Tilt     = spi[3];
            ib_params.Roll     = spi[4];
            ib_params.Twist    = spi[5];
            ib_params.X_disp   = hpi[0];
            ib_params.Y_disp   = hpi[1];
            ib_params.H_rise   = hpi[2];
            ib_params.Incl     = hpi[3];
            ib_params.Tip      = hpi[4];
            ib_params.H_twist  = hpi[5];
            _interbase_params.push_back(ib_params);
       }

       dsum /= (double) (num_bp - 1);
       if (dsum > 10.0 && dsum < 60.0) {
            int helix_value = DOUBLE_HELIX;
            for (int i = 0; i < num_bp - 1; ++i) {
                 if (twist[i] > 0.0 && aveS[i].x > -5.0 && aveS[i].x < -0.5 && aveS[i].y > 7.5 && aveS[i].y < 10.0 && aveS[i].z > -2.0 && aveS[i].z < 3.5 &&
                     aveH[i].x > -11.5 && aveH[i].x < 2.5 && aveH[i].y > 1.5 && aveH[i].y < 10.0 && aveH[i].z > -3.0 && aveH[i].z < 9.0) {
                      if (aveS[i].z >= 1.5) {
                           if (helix_value == A_DOUBLE_HELIX)
                                _classification |= A_DOUBLE_HELIX;
                           helix_value = A_DOUBLE_HELIX;
                      } else if (aveS[i].z <= 0.5) {
                           if (helix_value == B_DOUBLE_HELIX)
                                _classification |= B_DOUBLE_HELIX;
                           helix_value = B_DOUBLE_HELIX;
                      } else {
                           _classification |= DOUBLE_HELIX;
                           helix_value = DOUBLE_HELIX;
                      }
                 } else {
                      _classification |= DOUBLE_HELIX;
                      helix_value = DOUBLE_HELIX;
                 }
            }
       } else if (dsum < 0) {
            if (xdir == (ds * num_bp)) _classification |= Z_DOUBLE_HELIX;
            else _classification |= DOUBLE_HELIX;
       } else _classification |= DOUBLE_HELIX;
}

void BasePairInfo::_bpstep_par(COORD rot1[3], const COORD& org1, COORD rot2[3], const COORD& org2, double *pars, COORD mst_orien[3], COORD &mst_org)
{
       double phi, rolltilt;
       COORD t1, hinge, mstx, msty, para_bp1[3], para_bp2[3], temp[3];

       vector_cross_product(hinge, rot1[2], rot2[2]);
       rolltilt = magang(rot1[2], rot2[2]);

       arb_rotation(temp, hinge, -0.5 * rolltilt);
       multi_matrix(para_bp2, temp, rot2);
       arb_rotation(temp, hinge, 0.5 * rolltilt);
       multi_matrix(para_bp1, temp, rot1);

       pars[5] = vector_angle(para_bp1[1], para_bp2[1], para_bp2[2]);

       if (fabs(pars[5] - 180) < XEPS) {
            msty = para_bp2[0];
       } else if (fabs(pars[6] + 180) < XEPS) {
            msty = para_bp2[0];
            vector_reverse(msty);
       } else {
            vector_sum(msty, para_bp1[1], para_bp2[1]);
            vector_normalize(msty);
       }

       vector_cross_product(mstx, msty, para_bp2[2]);
       vector_middle(mst_org, org1, org2);
       vector_difference(t1, org2, org1);
       mst_orien[0] = mstx;
       mst_orien[1] = msty;
       mst_orien[2] = para_bp2[2];

       pars[0] = vector_dot_product(t1, mstx);
       pars[1] = vector_dot_product(t1, msty);
       pars[2] = vector_dot_product(t1, para_bp2[2]);

       phi = deg2rad(vector_angle(hinge, msty, para_bp2[2]));
       pars[4] = rolltilt * cos(phi);
       pars[3] = rolltilt * sin(phi);
}

void BasePairInfo::_helical_par(COORD rot1[3], const COORD& org1, COORD rot2[3], const COORD& org2, double *pars, COORD mst_orien[3], COORD &mst_org)
{
       double AD_mag, phi, TipInc1, TipInc2, vlen;
       COORD rot1_h[3], rot2_h[3], temp[3];
       COORD t1, t2, axis_h, hinge1, hinge2, AD_axis, org1_h, org2_h;

       vector_difference(t1, rot2[0], rot1[0]);
       vector_difference(t2, rot2[1], rot1[1]);
       vector_cross_product(axis_h, t1, t2);
       vlen = vector_length(axis_h);
       if (vlen < XEPS) {
            axis_h.x = 0.0; axis_h.y = 0.0; axis_h.z = 1.0;
       } else {
            axis_h.x /= vlen; axis_h.y /= vlen; axis_h.z /= vlen;
       }

       TipInc1 = magang(axis_h, rot1[2]);
       vector_cross_product(hinge1, axis_h, rot1[2]);
       arb_rotation(temp, hinge1, -TipInc1);
       multi_matrix(rot1_h, temp, rot1);

       TipInc2 = magang(axis_h, rot2[2]);
       vector_cross_product(hinge2, axis_h, rot2[2]);
       arb_rotation(temp, hinge2, -TipInc2);
       multi_matrix(rot2_h, temp, rot2);

       vector_sum(t1, rot1_h[0], rot2_h[0]);
       vector_sum(t2, rot1_h[1], rot2_h[1]);
       vector_normalize(t1);
       vector_normalize(t2);

       mst_orien[0] = t1;
       mst_orien[1] = t2;
       mst_orien[2] = axis_h;

       pars[5] = vector_angle(rot1_h[1], rot2_h[1], axis_h);

       vector_difference(t2, org2, org1);
       pars[2] = vector_dot_product(t2, axis_h);

       phi = deg2rad(vector_angle(hinge1, rot1_h[1], axis_h));
       pars[4] = TipInc1 * cos(phi);
       pars[3] = TipInc1 * sin(phi);

       vector_sum_scaled(t1, t2, -pars[2], axis_h);
       if (fabs(pars[5]) < HTWIST0)
            vector_sum_scaled(org1_h, org1, 0.5, t1);
       else {
            get_vector(t1, axis_h, 90 - pars[5] / 2.0, AD_axis);
            AD_mag = 0.5 * vector_length(t1) / sin(deg2rad(pars[5] / 2));
            vector_sum_scaled(org1_h, org1, AD_mag, AD_axis);
       }

       vector_sum_scaled(org2_h, org1_h, pars[2], axis_h);
       vector_middle(mst_org, org1_h, org2_h);
       vector_difference(t1, org1, org1_h);

       pars[0] = vector_dot_product(t1, rot1_h[0]);
       pars[1] = vector_dot_product(t1, rot1_h[1]);
}

void BasePairInfo::_project_xyzP(const COORD& P_1_i_plus_1, const COORD& P_2_i, COORD mst_orien[3], const COORD& mst_org, COORD& aveP)
{
       COORD temp, P_mst1, P_mst2;
       vector_difference(temp, P_1_i_plus_1, mst_org);
       P_mst1.x = vector_dot_product(temp, mst_orien[0]);
       P_mst1.y = vector_dot_product(temp, mst_orien[1]);
       P_mst1.z = vector_dot_product(temp, mst_orien[2]);
       vector_difference(temp, P_2_i, mst_org);
       P_mst2.x = vector_dot_product(temp, mst_orien[0]);
       P_mst2.y = -vector_dot_product(temp, mst_orien[1]);
       P_mst2.z = -vector_dot_product(temp, mst_orien[2]);
       vector_middle(aveP, P_mst1, P_mst2);
}
