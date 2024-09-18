/*
FILE:     CrySymmetry.h
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
#ifndef _CRYSTALL_CELL_
#define _CRYSTALL_CELL_

#include <string>
#include <vector>

#define NDB_TRANS_ORTHO_TO_FRAC   1
#define NDB_TRANS_FRAC_TO_ORTHO   2

#define NDB_BASE_SYMMETRY         "1_555"

#include "Atom.h"

typedef struct {
        double rot[3][3];
        double trans[3];
        int    tn[3], td[3];
} NDBSYMMETRY;

class CrySymmetry
{
   private:
       int __is_artifical;
       std::string __a, __b, __c, __alpha, __beta, __gamma;
       std::string __space_group;

       // double o_to_s[3][3];
       // double o_to_s_v[3];
       double f_to_o[3][3];
       double o_to_fc[3][3];
       std::vector<NDBSYMMETRY> __symops;
       void _init();
       bool _is_empty_or_1(const std::string& length);
       bool _is_empty_or_90(const std::string& angle);
       bool _is_empty_or_P1(const std::string& spg);
       void _set_sym_data(std::vector<std::vector<int> >& sym_data);
       void _ndb_ortho_to_frac_matrix();
       void _ndb_frac_to_ortho_matrix();
       void _ndb_sym_rot_op(const int *rot, char *sym_name);
       double _cal_symm_distance(const COORD &coord1, const COORD &coord2,
                                const std::string &SymOP);
   public:
       CrySymmetry();
       ~CrySymmetry();
       void clear();
       const int& is_artifical() const { return __is_artifical; }
       const std::string& a() const { return __a; }
       const std::string& b() const { return __b; }
       const std::string& c() const { return __c; }
       const std::string& alpha() const { return __alpha; }
       const std::string& beta() const { return __beta; }
       const std::string& gamma() const { return __gamma; }
       const std::string& space_group() const { return __space_group; }
       const std::vector<NDBSYMMETRY>& symops() const { return __symops; }
       void setArtificalFlag(const int& flag);
       void setA(const std::string& str);
       void setB(const std::string& str);
       void setC(const std::string& str);
       void setAlpha(const std::string& str);
       void setBeta(const std::string& str);
       void setGamma(const std::string& str);
       void setSpaceGroup(const std::string& str);

       void check_artifical_flag(const bool& is_crystal_entry);
       void SetCrySymmetry(std::vector<std::vector<int> >& sym_data);
       bool check_scale_matrix(double scale[3][3], double vec[3]);
       void symmetry_operation(COORD &coord, const int& op, const int& lx, const int& ly, const int& lz);
       void symmetry_operation(COORD &coord, const std::string &sym_op);
       void symmetry_operation(RCSB::Atom& atom, const std::string &sym_op);
       void symmetry_operation(double &x, double &y, double &z, const int& op, const int& lx, const int& ly, const int& lz);
       void do_symmetry_operation(COORD &new_coord, COORD old_coord, const NDBSYMMETRY& op);
       void ndb_trans_coord_3by3(COORD &coord, const int& iflag);
       void ndb_trans_coord_3by3(double &x, double &y, double &z, const int& iflag);
       void ndb_get_symmetry_operation_name(const int& op, const int& lx, const int& ly, const int& lz, std::string& sym_name);
       void ndb_get_symmetry_operation_name(const std::string &sym_op, std::string& sym_name);
       void ndb_get_symmetry_matrix(const int& op, const int& lx, const int& ly, const int& lz, double matrix[3][3], double dvector[3]);
       void ndb_get_symmetry_matrix(const std::string &sym_op, double matrix[3][3], double dvector[3]);
       void ndb_get_symmetry_operation(double matrix[3][3], double dvector[3], std::string &symmetry);
       void ndb_get_symmetry_operator(double matrix[3][3], double dvector[3], std::string &symmetry);
       double cal_symmetry_distance(RCSB::Atom *atom1, RCSB::Atom *atom2, const std::string &SymOP1, const std::string &SymOP2);
       void get_symmetry_number(const std::string &sym_operator, int &op, int &lx, int &ly, int &lz);
       double o_to_f[3][3];
       double o_to_f_v[3];
};

#endif
