/*
FILE:     CrySymmetry.C
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
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "CrySymmetry.h"
#include "MatrixUtil.h"
#include "utillib.h"

CrySymmetry::CrySymmetry()
{
       _init();
}

CrySymmetry::~CrySymmetry()
{
       _init();
}

void CrySymmetry::_init()
{
       __is_artifical = 0;
       __a.clear(); __b.clear(); __c.clear();
       __alpha.clear(); __beta.clear(); __gamma.clear();
       __space_group.clear();
       __symops.clear();
}

void CrySymmetry::clear()
{
       _init();
}

void CrySymmetry::setArtificalFlag(const int& flag)     { __is_artifical = flag; }
void CrySymmetry::setA(const std::string& str)          { __a = str; }
void CrySymmetry::setB(const std::string& str)          { __b = str; }
void CrySymmetry::setC(const std::string& str)          { __c = str; }
void CrySymmetry::setAlpha(const std::string& str)      { __alpha = str; }
void CrySymmetry::setBeta(const std::string& str)       { __beta = str; }
void CrySymmetry::setGamma(const std::string& str)      { __gamma = str; }
void CrySymmetry::setSpaceGroup(const std::string& str) { __space_group = str; }

void CrySymmetry::check_artifical_flag(const bool& is_crystal_entry)
{
       if (_is_empty_or_1(__a) && _is_empty_or_1(__b) && _is_empty_or_1(__c) && (_is_empty_or_1(__alpha) || _is_empty_or_90(__alpha)) &&
          (_is_empty_or_1(__beta) ||  _is_empty_or_90(__beta)) && (_is_empty_or_1(__gamma) || _is_empty_or_90(__gamma)) &&
          (_is_empty_or_P1(__space_group) || !is_crystal_entry)) __is_artifical = 1;
}

void CrySymmetry::SetCrySymmetry(std::vector<std::vector<int> >& sym_data)
{
       _set_sym_data(sym_data);
       _ndb_ortho_to_frac_matrix();
       _ndb_frac_to_ortho_matrix();
}

bool CrySymmetry::check_scale_matrix(double scale[3][3], double vec[3])
{
       bool ierror = false;
       bool wrong_scale = false;
       for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                 if (fabs(scale[i][j]) > 0.00001) {
                      if (fabs(scale[i][j]) > 0.8) wrong_scale = true; 
                      ierror = true; break;
                 }
            }
            if (ierror) break;
       }
       if (!ierror || wrong_scale) return false;

       ierror = false;
       for (int i = 0; i < 3; i++) o_to_f_v[i] = vec[i];
       for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                 if (fabs(scale[i][j] - o_to_f[i][j]) > 0.00005)
                      ierror = true;
            }
       }
       if (!ierror) return true;

       for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                 o_to_f[i][j] = scale[i][j];
            }
       }

       double v = scale[0][0] * scale[1][1] * scale[2][2] + scale[1][0] * scale[2][1] * scale[0][2] + scale[2][0] * scale[1][2] * scale[0][1] -
                 (scale[0][2] * scale[1][1] * scale[2][0] + scale[0][0] * scale[1][2] * scale[2][1] + scale[1][0] * scale[2][2] * scale[0][1]);

       if (fabs(v) > 0.0) {
            f_to_o[0][0] = (scale[1][1] * scale[2][2] - scale[1][2] * scale[2][1]) / v;
            f_to_o[0][1] = (scale[2][1] * scale[0][2] - scale[2][2] * scale[0][1]) / v;
            f_to_o[0][2] = (scale[0][1] * scale[1][2] - scale[1][1] * scale[0][2] ) / v;
            f_to_o[1][0] = (scale[1][2] * scale[2][0] - scale[2][2] * scale[1][0] ) / v;
            f_to_o[1][1] = (scale[0][0] * scale[2][2] - scale[0][2] * scale[2][0] ) / v;
            f_to_o[1][2] = (scale[0][2] * scale[1][0] - scale[1][2] * scale[0][0] ) / v;
            f_to_o[2][0] = (scale[1][0] * scale[2][1] - scale[2][0] * scale[1][1] ) / v;
            f_to_o[2][1] = (scale[2][0] * scale[0][1] - scale[2][1] * scale[0][0] ) / v;
            f_to_o[2][2] = (scale[0][0] * scale[1][1] - scale[1][0] * scale[0][1] ) / v;
       }
       return true;
}

void CrySymmetry::symmetry_operation(COORD &coord, const int& op, const int& lx, const int& ly, const int& lz)
{
       ndb_trans_coord_3by3(coord, NDB_TRANS_ORTHO_TO_FRAC);
       do_symmetry_operation(coord, coord, __symops[op]);
       coord.x += (double) lx;
       coord.y += (double) ly;
       coord.z += (double) lz;
       ndb_trans_coord_3by3(coord, NDB_TRANS_FRAC_TO_ORTHO);
}

void CrySymmetry::symmetry_operation(COORD &coord, const std::string &sym_op)
{
       int op, lx, ly, lz;
       get_symmetry_number(sym_op, op, lx, ly, lz);
       symmetry_operation(coord, op, lx, ly, lz);
}

void CrySymmetry::symmetry_operation(RCSB::Atom& atom, const std::string &sym_op)
{
       if (sym_op.empty() || sym_op == "1_555") return;

       COORD coord = atom.orig();
       symmetry_operation(coord, sym_op);
       atom.set_orig(coord);
}

void CrySymmetry::symmetry_operation(double &x, double &y, double &z, const int& op, const int& lx, const int& ly, const int& lz)
{
       COORD coord;
       coord.x = x;
       coord.y = y;
       coord.z = z;
       symmetry_operation(coord, op, lx, ly, lz);
       x = coord.x;
       y = coord.y;
       z = coord.z;
}

void CrySymmetry::do_symmetry_operation(COORD &new_coord, COORD old_coord, const NDBSYMMETRY& op)
{
       new_coord.x = old_coord.x * op.rot[0][0] + old_coord.y * op.rot[0][1] + old_coord.z * op.rot[0][2] + op.trans[0];
       new_coord.y = old_coord.x * op.rot[1][0] + old_coord.y * op.rot[1][1] + old_coord.z * op.rot[1][2] + op.trans[1];
       new_coord.z = old_coord.x * op.rot[2][0] + old_coord.y * op.rot[2][1] + old_coord.z * op.rot[2][2] + op.trans[2];
}

void CrySymmetry::ndb_trans_coord_3by3(COORD &coord, const int& iflag)
{
       double x0[3], t[3];

       x0[0] = coord.x;
       x0[1] = coord.y;
       x0[2] = coord.z;

       if (iflag == NDB_TRANS_ORTHO_TO_FRAC) {
            for (int i = 0; i < 3; i++) {
                t[i] = 0.0;
                for (int j = 0; j < 3; j++) {
                    t[i] = t[i] + o_to_f[i][j] * x0[j];
                }
            }
            coord.x = t[0] + o_to_f_v[0];
            coord.y = t[1] + o_to_f_v[1];
            coord.z = t[2] + o_to_f_v[2];
       } else if (iflag == NDB_TRANS_FRAC_TO_ORTHO) {
            for (int i = 0; i < 3; i++) {
                 x0[i] -= o_to_f_v[i];
            }
            for (int i = 0; i < 3; i++) {
                t[i] = 0.0;
                for (int j = 0; j < 3; j++) {
                    t[i] = t[i] + f_to_o[i][j] * x0[j];
                }
            }
            coord.x = t[0];
            coord.y = t[1];
            coord.z = t[2];
       }
}

void CrySymmetry::ndb_trans_coord_3by3(double &x, double &y, double &z, const int& iflag)
{
       COORD coord;
       coord.x = x;
       coord.y = y;
       coord.z = z;
       ndb_trans_coord_3by3(coord, iflag);
       x = coord.x;
       y = coord.y;
       z = coord.z;
}

void CrySymmetry::ndb_get_symmetry_operation_name(const int& op, const int& lx, const int& ly, const int& lz, std::string& sym_name)
{
       sym_name.clear();
       if (op >= (int) __symops.size()) return;

       int t[3], rot[3][3], tn[3], td[3];
       char buffer[80];

       std::string rotation, translation;
       t[0] = lx; t[1] = ly; t[2] = lz;
       for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                 rot[i][j] = (int) __symops[op].rot[i][j];
            }
            tn[i] = __symops[op].tn[i];
            td[i] = __symops[op].td[i];
            if (td[i] != 0) tn[i] += t[i] * td[i];
            else tn[i] = t[i];
            memset(buffer, 0, 80);
            _ndb_sym_rot_op(&rot[i][0], buffer);
            rotation = buffer;
            translation.clear();
            if (tn[i] != 0) {
                 if (td[i]) {
                      if (tn[i] / td[i] * td[i] == tn[i]) {
                           tn[i] = tn[i] / td[i];
                           translation = String::IntToString(tn[i]);
                      } else translation = String::IntToString(tn[i]) + "/" + String::IntToString(td[i]);
                 } else translation = String::IntToString(tn[i]);
                 if (translation[0] == '-')
                      rotation += translation;
                 else rotation += "+" + translation;
            }
            sym_name += rotation;
            if (i < 2) sym_name += ",";
       }
}

void CrySymmetry::ndb_get_symmetry_operation_name(const std::string &sym_op, std::string &sym_name)
{
       sym_name.clear();
       int op, lx, ly, lz;
       get_symmetry_number(sym_op, op, lx, ly, lz);
       ndb_get_symmetry_operation_name(op, lx, ly, lz, sym_name);
}

void CrySymmetry::ndb_get_symmetry_matrix(const int& op, const int& lx, const int& ly, const int& lz, double matrix[3][3], double dvector[3])
{
       if (op >= (int) __symops.size()) return;

       double matrix1[3][3], tt[3];

       tt[0] = (double) lx + __symops[op].trans[0];
       tt[1] = (double) ly + __symops[op].trans[1];
       tt[2] = (double) lz + __symops[op].trans[2];

       multi_matrix(matrix1, f_to_o, __symops[op].rot);
       multi_matrix(matrix, matrix1, o_to_f);

       for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                 if (fabs(matrix[i][j]) < 0.0001) matrix[i][j] = 0.0;
                 if (fabs(1.0 - fabs(matrix[i][j])) < 0.0001)
                      matrix[i][j] = matrix[i][j] / fabs(matrix[i][j]);
            }
       }

       multi_matrix_vector(dvector, f_to_o, tt);

       for (int i = 0; i < 3; i++) {
            if (fabs(dvector[i]) < 0.0001) dvector[i] = 0.0;
       }
}

void CrySymmetry::ndb_get_symmetry_matrix(const std::string &sym_op, double matrix[3][3], double dvector[3])
{
       int op, lx, ly, lz;
       get_symmetry_number(sym_op, op, lx, ly, lz);
       ndb_get_symmetry_matrix(op, lx, ly, lz, matrix, dvector);
}

void CrySymmetry::ndb_get_symmetry_operation(double matrix[3][3], double dvector[3], std::string &symmetry)
{
       symmetry.clear();

       double matrix1[3][3], rot[3][3], tt[3];
       multi_matrix(matrix1, o_to_f, matrix);
       multi_matrix(rot, matrix1, f_to_o);
       multi_matrix_vector(tt, o_to_f, dvector);
       ndb_get_symmetry_operator(rot, tt, symmetry);
}

void CrySymmetry::ndb_get_symmetry_operator(double matrix[3][3], double dvector[3], std::string &symmetry)
{
       symmetry.clear();
       int op = 0;
       for (unsigned int i = 0; i < __symops.size(); ++i) {
            bool found = true;
            for (int j = 0; j < 3; j++) {
                 for (int k = 0; k < 3; k++) {
                      int x = int(rint(matrix[j][k]));
                      int y = int(rint(__symops[i].rot[j][k]));
                      if (x != y) {
                           found = false;
                           break;
                      }
                 }
                 if (!found) break;
            }
            if (found) {
                 double r = dvector[0] - __symops[i].trans[0];
                 int x = int(rint(r));
                 if (fabs(r - (double) x) > 0.1) continue;
                 r = dvector[1] - __symops[i].trans[1];
                 x = int(rint(r));
                 if (fabs(r - (double) x) > 0.1) continue;
                 r = dvector[2] - __symops[i].trans[2];
                 x = int(rint(r));
                 if (fabs(r - (double) x) > 0.1) continue;
                 op = i + 1;
                 break;
            }
       }
       if (op == 0) return;

       int lx = int(rint(dvector[0] - __symops[op-1].trans[0]));
       int ly = int(rint(dvector[1] - __symops[op-1].trans[1]));
       int lz = int(rint(dvector[2] - __symops[op-1].trans[2]));
       symmetry = String::IntToString(op) + "_" + String::IntToString(lx + 5) + String::IntToString(ly + 5) + String::IntToString(lz + 5);
}

double CrySymmetry::cal_symmetry_distance(RCSB::Atom *atom1, RCSB::Atom *atom2, const std::string &SymOP1, const std::string &SymOP2)
{
       double dist = 0.0;
       if (__is_artifical || ((SymOP1.empty() || SymOP1 == "1_555") && (SymOP2.empty() || SymOP2 == "1_555")))
            dist = cal_distance(atom1, atom2);
       else if (SymOP2.empty() || SymOP2 == "1_555")
            dist = _cal_symm_distance(atom2->orig(), atom1->orig(), SymOP1);
       else dist = _cal_symm_distance(atom1->orig(), atom2->orig(), SymOP2);
       return dist;
}

void CrySymmetry::get_symmetry_number(const std::string &sym_operator, int &op, int &lx, int &ly, int &lz)
{
       op = lx = ly = lz = 0;
       std::string::size_type p = sym_operator.find("_");
       if (p != std::string::npos && sym_operator.substr(p).size() == 4) {
            std::string cs = sym_operator.substr(0, p);
            op = atoi(cs.c_str()) - 1;
            lx = sym_operator[p + 1] - '5';
            ly = sym_operator[p + 2] - '5';
            lz = sym_operator[p + 3] - '5';
       }
}

bool CrySymmetry::_is_empty_or_1(const std::string& length)
{
       if (length.empty()) return true;
       if (fabs(atof(length.c_str())) < 0.001) return true;
       if (fabs(atof(length.c_str()) - 1.0) < 0.001) return true;
       return false;
}

bool CrySymmetry::_is_empty_or_90(const std::string& angle)
{
       if (angle.empty()) return true;
       if (fabs(atof(angle.c_str())) < 0.001) return true;
       if (fabs(atof(angle.c_str()) - 90.0) < 0.001) return true;
       return false;
}

bool CrySymmetry::_is_empty_or_P1(const std::string& spg)
{
       if (spg.empty() || spg == "P 1") return true;
       return false;
}

void CrySymmetry::_set_sym_data(std::vector<std::vector<int> >& sym_data)
{
       if (sym_data.empty()) return;

       NDBSYMMETRY symop;
       for (std::vector<std::vector<int> >::iterator
            pos = sym_data.begin(); pos != sym_data.end(); ++pos) {
            for (int i = 0; i < 3; ++i) {
                 for (int j = 0; j < 3; ++j) {
                      symop.rot[i][j] = (double) (*pos)[i * 3 + j];
                 }
                 symop.tn[i] = (*pos)[9 + i * 2];
                 symop.td[i] = (*pos)[9 + i * 2 + 1];
                 symop.trans[i] = 0;
                 if (symop.td[i])
                      symop.trans[i] = (double) symop.tn[i] / (double) symop.td[i]; 
            }
            __symops.push_back(symop);
       }
}

void CrySymmetry::_ndb_ortho_to_frac_matrix()
{
       double aa[3], ang[3], alphr[3], cosa[3], sina[3];

       if (__is_artifical) {
            aa[0] = 1000.0;
            aa[1] = 1000.0;
            aa[2] = 1000.0;
            ang[0] = 90.0;
            ang[1] = 90.0;
            ang[2] = 90.0;
       } else {
            aa[0] = atof(__a.c_str());
            aa[1] = atof(__b.c_str());
            aa[2] = atof(__c.c_str());
            ang[0] = atof(__alpha.c_str());
            ang[1] = atof(__beta.c_str());
            ang[2] = atof(__gamma.c_str());
       }

       double rad = acos(-1.0)/180.0;
       for (int i = 0; i < 3; i++) {
            alphr[i] = ang[i] * rad;
            cosa[i]  = cos(alphr[i]);
            sina[i]  = sin(alphr[i]);
            if (fabs(ang[i] - 90.0) < 0.001) {
                 cosa[i] = 0.0;
                 sina[i] = 1.0;
            }
       }
       double va = sqrt(1.00 + 2.0 * cosa[0] * cosa[1] * cosa[2] - cosa[0] * cosa[0] - cosa[1] * cosa[1] - cosa[2] * cosa[2]);
     
       o_to_fc[0][0] = 1.0 / aa[0];
       o_to_fc[0][1] = -cosa[2] / (aa[0] * sina[2]);
       o_to_fc[0][2] = (cosa[2] * cosa[0]-cosa[1]) / (aa[0] * va * sina[2]);
       o_to_fc[1][0] = 0.0;
       o_to_fc[1][1] = 1.0 / (aa[1] * sina[2]);
       o_to_fc[1][2] = (cosa[1] * cosa[2]-cosa[0]) / (aa[1] * va * sina[2]);
       o_to_fc[2][0] = 0.0;
       o_to_fc[2][1] = 0.0;
       o_to_fc[2][2] = sina[2] / (aa[2] * va);

       o_to_f[0][0] = 1.0 / aa[0];
       o_to_f[0][1] = -cosa[2] / (aa[0] * sina[2]);
       o_to_f[0][2] = (cosa[2] * cosa[0]-cosa[1]) / (aa[0] * va * sina[2]);
       o_to_f[1][0] = 0.0;
       o_to_f[1][1] = 1.0 / (aa[1] * sina[2]);
       o_to_f[1][2] = (cosa[1] * cosa[2]-cosa[0]) / (aa[1] * va * sina[2]);
       o_to_f[2][0] = 0.0;
       o_to_f[2][1] = 0.0;
       o_to_f[2][2] = sina[2] / (aa[2] * va);

       o_to_f_v[0] = o_to_f_v[1] = o_to_f_v[2] = 0.0;
}

void CrySymmetry::_ndb_frac_to_ortho_matrix()
{
       double aa[3], ang[3], alphr[3], cosa[3], sina[3];

       if (__is_artifical) {
            aa[0] = 1000.0;
            aa[1] = 1000.0;
            aa[2] = 1000.0;
            ang[0] = 90.0;
            ang[1] = 90.0;
            ang[2] = 90.0;
       } else {
            aa[0] = atof(__a.c_str());
            aa[1] = atof(__b.c_str());
            aa[2] = atof(__c.c_str());
            ang[0] = atof(__alpha.c_str());
            ang[1] = atof(__beta.c_str());
            ang[2] = atof(__gamma.c_str());
       }

       double rad = acos(-1.0) / 180.0;
       for (int i = 0; i < 3; i++) {
            alphr[i] = ang[i] * rad;
            cosa[i]  = cos(alphr[i]);
            sina[i]  = sin(alphr[i]);
            if (fabs(ang[i] - 90.0) < 0.001) {
                 cosa[i]=0.0;
                 sina[i]=1.0;
            }
       }
       double va = sqrt(1.0 + 2.0 * cosa[0] * cosa[1] * cosa[2] - cosa[0] * cosa[0] - cosa[1] * cosa[1] - cosa[2] * cosa[2]);
       double v = aa[0] * aa[1] * aa[2] * va;

       f_to_o[0][0] = aa[0];
       f_to_o[0][1] = aa[1] * cosa[2];
       f_to_o[0][2] = aa[2] * cosa[1];
       f_to_o[1][0] = 0.0;
       f_to_o[1][1] = aa[1] * sina[2];
       f_to_o[1][2] = aa[2] * (cosa[0] - cosa[1] * cosa[2]) / sina[2];
       f_to_o[2][0] = 0.0;
       f_to_o[2][1] = 0.0;
       f_to_o[2][2] = v / (aa[0] * aa[1] * sina[2]);
}

void CrySymmetry::_ndb_sym_rot_op(const int *rot, char *sym_name)
{
       int i, len;
       for (i = 0; i < 3; i++) {
            len = strlen(sym_name);
            if (rot[i] > 0) {
                 if (i == 0) sprintf(&sym_name[len], "%c", 120+i);
                 else sprintf(&sym_name[len], "+%c", 120+i);
            } else if (rot[i] < 0) {
                 sprintf(&sym_name[len], "-%c", 120+i);
            }
       }
       if (sym_name[0] == '+') {
            len = strlen(sym_name);
            for (i = 0; i < len; i++) {
                 sym_name[i] = sym_name[i+1];
            }
       }
}

double CrySymmetry::_cal_symm_distance(const COORD &coord1, const COORD &coord2, const std::string &SymOP)
{
       COORD coord = coord2;

       int op, lx, ly, lz;
       get_symmetry_number(SymOP, op, lx, ly, lz);
       symmetry_operation(coord, op, lx, ly, lz);

       double dx = coord1.x - coord.x;
       double dy = coord1.y - coord.y;
       double dz = coord1.z - coord.z;

       return (sqrt(dx * dx + dy * dy + dz * dz));
}
