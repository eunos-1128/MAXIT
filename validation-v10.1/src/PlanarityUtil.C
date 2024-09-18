/*
FILE:     PlanarityUtil.C
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
#include <math.h>

#include "GetPairList.h"
#include "GetPairList.C"
#include "PlanarityUtil.h"
#include "PlanarityUtil_global.h"

#define SIGN(a, b) ((b) < 0 ? -fabs(a) : fabs(a))

PlanarityUtil::PlanarityUtil()
{
       _clear();
       _init();
}

PlanarityUtil::~PlanarityUtil()
{
       _clear();
}

void PlanarityUtil::_clear()
{
       _planeIndex.clear();
}

void PlanarityUtil::_init()
{
       for (int i = 0; i < NUM_PLANE_RESIDUE; ++i) {
            _planeIndex.insert(std::make_pair(_plane_atoms[i].resname, i));
       }
}

bool PlanarityUtil::CheckPlane(RCSB::Residue* residue, double& rmsd, unsigned int& plane_count)
{
       rmsd = 0;

       std::map<std::string, int>::iterator mpos = _planeIndex.find(residue->ResName());
       if (mpos == _planeIndex.end()) return true;

       std::vector<RCSB::Atom*> atom_list;
       std::vector<std::vector<RCSB::Atom*> > atom_lists, pair_lists;

       atom_lists.clear();
       for (int i = 0; i < _plane_atoms[mpos->second].num_atom; ++i) {
            residue->find_atom(_plane_atoms[mpos->second].atoms[i], atom_list);
            if (atom_list.empty()) {
                 atom_lists.clear();
                 break;
            }
            atom_lists.push_back(atom_list);
       }
       if (atom_lists.empty()) return true;

       GetPairList(atom_lists, pair_lists);
       if (pair_lists.empty()) return true;

       plane_count++;

       bool return_flag = true;
       int count = 0;
       for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = pair_lists.begin(); apos != pair_lists.end(); ++apos) {
            double rmsd1 = _calculatePlaneRMSD(*apos, count);
            if (count) return_flag = false;
            if (rmsd1 > rmsd) rmsd = rmsd1;
       }

       if (rmsd > PLANECUTOFF) return_flag = false;
       return return_flag;
}

double PlanarityUtil::_calculatePlaneRMSD(const std::vector<RCSB::Atom*>& atoms, int& count)
{
       count = 0;
       double rmsd = 0;

       if (atoms.empty()) return rmsd;

       double center[3], d[3], e[3], mat[3][3];
       center[0] = center[1] = center[2] = 0;
       for (std::vector<RCSB::Atom*>::const_iterator pos = atoms.begin(); pos != atoms.end(); ++pos) {
            center[0] += (*pos)->orig().x;
            center[1] += (*pos)->orig().y;
            center[2] += (*pos)->orig().z;
       }
       for (int i = 0; i < 3; ++i) {
            center[i] /= (double) atoms.size(); 
       }
       mat[0][0] = mat[0][1] = mat[0][2] = mat[1][1] =
       mat[1][2] = mat[2][2] = 0;
       for (std::vector<RCSB::Atom*>::const_iterator pos = atoms.begin(); pos != atoms.end(); ++pos) {
            mat[0][0] += ((*pos)->orig().x  - center[0])
                       * ((*pos)->orig().x  - center[0]);
            mat[0][1] += ((*pos)->orig().x  - center[0])
                       * ((*pos)->orig().y  - center[1]);
            mat[0][2] += ((*pos)->orig().x  - center[0])
                       * ((*pos)->orig().z  - center[2]);
            mat[1][1] += ((*pos)->orig().y  - center[1])
                       * ((*pos)->orig().y  - center[1]);
            mat[1][2] += ((*pos)->orig().y  - center[1])
                       * ((*pos)->orig().z  - center[2]);
            mat[2][2] += ((*pos)->orig().z  - center[2])
                       * ((*pos)->orig().z  - center[2]);
       }
       mat[1][0] = mat[0][1];
       mat[2][0] = mat[0][2];
       mat[2][1] = mat[1][2];
       tred2(3, &mat[0][0], d, e);
       if (!tqli(3, &mat[0][0], d, e)) {
            eigsrt(3, &mat[0][0], d);
            double cal_rmsd = 0;
            for (std::vector<RCSB::Atom*>::const_iterator pos = atoms.begin(); pos != atoms.end(); ++pos) {
                 double dot = ((*pos)->orig().x  - center[0]) * mat[0][2]
                            + ((*pos)->orig().y  - center[1]) * mat[1][2]
                            + ((*pos)->orig().z  - center[2]) * mat[2][2];
                 cal_rmsd += dot * dot;
                 if (fabs(dot) > PLANECUTOFF) count++;
            }
            cal_rmsd = sqrt(cal_rmsd / (double) atoms.size());
            if (cal_rmsd > rmsd) rmsd = cal_rmsd;
       }
       return rmsd;
}

void PlanarityUtil::tred2(const int n, double *a, double *d, double *e)
{
       for (int i = n - 1; i > 0; i--) {
            int l = i - 1;
            double h = 0;
            double scale = 0;
            if (l > 0) {
                 for (int k = 0; k <= l; k++) scale += fabs(a[i * n + k]);
                 if (scale == 0.0)
                      e[i] = a[i * n + l];
                 else {
                      for (int k = 0; k <= l; k++) {
                           a[i * n + k] /= scale;
                           h += a[i * n + k] * a[i * n + k]; 
                      }
                      double f = a[i * n + l];
                      double g = f > 0 ? -sqrt(h) : sqrt(h);
                      e[i] = scale * g;
                      h -= f * g;
                      a[i * n + l] = f - g;
                      f = 0.0;
                      for (int j = 0; j <= l; ++j) {
                           a[j * n + i] = a[i * n + j] / h;
                           g = 0.0;
                           for (int k = 0; k <= j; k++) {
                                g += a[j * n + k] * a[i * n + k];
                           }
                           for (int k = j + 1; k <= l; k++) {
                                g += a[k * n + j] * a[i * n + k];
                           }
                           e[j] = g / h;
                           f += e[j] * a[i * n + j];
                      }
                      double hh = f / ( h + h);
                      for (int j = 0; j <= l; ++j) {
                           f = a[i * n + j];
                           e[j] = g = e[j] - hh * f;
                           for (int k = 0; k <= j; k++) {
                                a[j * n + k] -= (f * e[k] + g * a[i * n + k]);
                           }
                      }
                 }
            } else e[i] = a[i * n + l];
            d[i] = h;
       }

       d[0] = 0;
       e[0] = 0;
       for (int i = 0; i < n; ++i) {
            int l = i - 1;
            if (fabs(d[i]) > 0.0) {
                 for (int j = 0; j <= l; ++j) {
                      double g = 0;
                      for (int k = 0; k <= l; k++) {
                           g += a[i * n + k] * a[k * n + j];
                      }
                      for (int k = 0; k <= l; k++) {
                           a[k * n + j] -= g * a[k * n + i];
                      }
                 }
            }
            d[i] = a[i * n + i];
            a[i * n + i] = 1;
            for (int j = 0; j <= l; ++j) {
                 a[j * n + i] = a[i * n + j] = 0;
            }
       }
}

int PlanarityUtil::tqli(const int n, double *a, double *d, double *e)
{
       int m = 0;
       for (int i = 1; i < n; ++i) e[i - 1] = e[i];
       e[n - 1] = 0.0;
       for (int l = 0; l < n; ++l) {
            int iter = 0;
            do {
                for (m = l; m < n - 1; ++m) {
                     double dd = fabs(d[m]) + fabs(d[m + 1]);
                     if ((float) (fabs(e[m]) + dd) == (float) dd) break;
                }
                if (m != l) {
                     iter++;
                     if (iter == 30) return 1;
                     double g = (d[l + 1] - d[l]) / (2 * e[l]);
                     double r = sqrt(g * g + 1);
                     g = d[m] - d[l] + e[l] / (g + SIGN(r, g));
                     double s = 1;
                     double c = 1;
                     double p = 0;
                     for (int i = m - 1; i >= l; i--) {
                          double f = s * e[i];
                          double b = c * e[i];
                          if (fabs(f) >= fabs(g)) {
                               c = g / f;
                               r = sqrt(c * c + 1);
                               e[i + 1] = f * r;
                               c *= (s = 1 / r);
                          } else {
                               s = f / g;
                               r = sqrt(s * s + 1);
                               e[i + 1] = g * r;
                               s *= (c = 1 / r);
                          }
                          g = d[i + 1] - p;
                          r = (d[i] - g) * s + 2 * c * b;
                          p = s * r;
                          d[i + 1] = g + p;
                          g = c * r - b;
                          for (int k = 0; k < n; k++) {
                               f = a[k * n + i + 1];
                               a[k * n + i + 1] = s * a[k * n + i] + c * f;
                               a[k * n + i] = c * a[k * n + i] - s * f;
                          }
                     }
                     d[l] -= p;
                     e[l] = g;
                     e[m] = 0;
                }
            } while (m != l);
       }
       return 0;
}

void PlanarityUtil::eigsrt(const int n, double *a, double *d)
{
       int i, j, k;
       double p;

       for (i = 0; i < n - 1; ++i) {
            p = d[k = i];
            for (j = i + 1; j < n; ++j) {
                 if (d[j] >= p) p = d[k = j];
            }
            if (k != i) {
                 d[k] = d[i];
                 d[i] = p;
                 for (j = 0; j < n; ++j) {
                      p = a[j * n + i];
                      a[j * n + i] = a[j * n + k];
                      a[j * n + k] = p;
                 }
            }
       }
}
