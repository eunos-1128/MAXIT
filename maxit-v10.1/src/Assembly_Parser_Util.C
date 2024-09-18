/*
FILE:     Assembly_Parser_Util.C
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

#include "Assembly_Parser_Util.h"
#include "MatrixUtil.h"
#include "utillib.h"

// static void separate_string_by_parenthesis(std::string &cifstring, std::vector<std::string> &list);
// static void separate_string_by_comma(std::string &cifstring, std::vector<std::string> &list);
static bool get_matrix_from_oper_list(const std::string& id_value, const std::vector<SymmMatrix>& SymmMatrices,
                                      const std::map<std::string, unsigned int>& SymmMatrix_Mapping, CrySymmetry& cell, double matrix[4][4]);

#if 0
move to other-util.C 
int parseString(const std::string &orginal_string, std::vector<std::string> &list)
{
       std::string op_string;
       String::RemoveWhiteSpace(orginal_string, op_string);
       int has_special_character = 0;
       for (unsigned int i = 0; i < op_string.size(); ++i) {
            if (op_string[i] == '(' || op_string[i] == ')' ||
                op_string[i] == ',' || op_string[i] == '-') {
                 has_special_character = 1;
                 break;
            }
       }
       if (!has_special_character) {
            list.push_back(op_string);
            return 0;
       }

       int numlp = 0;
       int numrp = 0;
       int numd = 0;
       int firstl = -1;
       int firstr = -1;
       int firstc = -1;
      
       for (unsigned int i = 0; i < op_string.size(); ++i) {
            if (op_string[i] == '(') {
                 numlp++;
                 if (firstl < 0) firstl = i;
            } else if (op_string[i] == ')') {
                 numrp++;
                 if (firstr < 0) firstr = i;
            } else if (op_string[i] == ',') {
                 if (firstc < 0) firstc = i;
            } else if (op_string[i] == '-') {
                 numd++;
            }
       } 

       if (numlp != numrp) return 2;
       if (numlp && numrp && firstr < firstl) return 2;

       int choice = 0;
       if (firstl >= 0) {
            if (firstc >= 0 && firstc < firstl)
                 choice = 2;
            else choice = 1;
       } else if (firstc >= 0) {
            if (firstl >= 0 && firstl < firstc)
                 choice = 1;
            else choice = 2;
       } else if (numd == 1) choice = 3;

       std::vector<std::string> list1, list2;
       if (choice == 1) {
            std::vector<std::vector<std::string> > lists;
            lists.clear();
            separate_string_by_parenthesis(op_string, list1); 
            if (list1.empty()) return 2;
            for (unsigned int i = 0; i < list1.size(); ++i) {
                 list2.clear();
                 int error = parseString(list1[i], list2);
                 if (error) return error;
                 lists.push_back(list2);
            }
            if (lists.empty()) return 2;

            list2 = lists[lists.size() - 1];
            for (int i = (int) lists.size() - 2; i >= 0; i--) {
                 list1.clear();
                 for (unsigned int j = 0; j < lists[i].size(); ++j) {
                      for (unsigned int k = 0; k < list2.size(); ++k) {
                           std::string cs = lists[i][j] + " " + list2[k];
                           list1.push_back(cs);
                      }
                 }
                 list2 = list1;
            }
            for (unsigned int i = 0; i < list2.size(); ++i) {
                 list.push_back(list2[i]);
            }
            return 0;
       } else if (choice == 2) {
            separate_string_by_comma(op_string, list1);
            if (list1.empty()) return 3;
            for (unsigned int i = 0; i < list1.size(); ++i) {
                 list2.clear();
                 int error = parseString(list1[i], list2);
                 if (error) return error;
                 for (unsigned int j = 0; j < list2.size(); ++j) {
                      list.push_back(list2[j]);
                 }
            }
            return 0;
       } else if (choice == 3) {
            std::string::size_type pos = op_string.find("-");
            if (pos == std::string::npos) return 3;
            int start = atoi(op_string.substr(0, pos).c_str());
            int end = atoi(op_string.substr(pos + 1).c_str());
            if (end <= start) return 3;
            if (end == 0) return 3;
            for (int i = start; i <= end; ++i) {
                 list.push_back(String::IntToString(i));
            }
            return 0;
       } else return 4;
}
#endif

bool getFinalMatrix(double mat[4][4], const std::vector<std::string> &id_list, const std::vector<SymmMatrix>& SymmMatrices,
                    const std::map<std::string, unsigned int>& SymmMatrix_Mapping, CrySymmetry& cell)
{
       double mat1[4][4], mat2[4][4];

       if (!id_list.empty()) {
            if (!get_matrix_from_oper_list(id_list[id_list.size()-1], SymmMatrices, SymmMatrix_Mapping, cell, mat))
                 return false;
            for (int i = (int) id_list.size() - 2; i >= 0; i--) {
                 if (!get_matrix_from_oper_list(id_list[i], SymmMatrices, SymmMatrix_Mapping, cell, mat1))
                      return false;
                 multi_matrix_4(mat2, mat1, mat);
                 for (int j = 0; j < 4; j++) {
                      for (int k = 0; k < 4; k++) {
                           mat[j][k] = mat2[j][k];
                      }
                 }
            }
            return true;
       }
       return false;
}

#if 0
static void separate_string_by_parenthesis(std::string &cifstring, std::vector<std::string> &list)
{
       list.clear();
       std::string cs;
       cs.clear();
       int count = 0;
       for (unsigned int i = 0; i < cifstring.size(); i++) {
            if (cifstring[i] == ' ' || cifstring[i] == '\t' ||
                cifstring[i] == '\n') continue;
            if (cifstring[i] == '(') {
                 if (count != 0) cs += cifstring[i];
                 count++;
            } else if (cifstring[i] == ')') {
                 count--;
                 if (count == 0) {
                      if (cs != "") list.push_back(cs);
                      cs.clear();
                 } else cs += cifstring[i];
            } else cs += cifstring[i];
       }
}

static void separate_string_by_comma(std::string &cifstring, std::vector<std::string> &list)
{
       list.clear();
       std::string cs;
       cs.clear();
       int count = 0;
       for (unsigned int i = 0; i < cifstring.size(); ++i) {
            if (cifstring[i] == ' ' || cifstring[i] == '\t' ||
                cifstring[i] == '\n') continue;
            if (cifstring[i] == ',' && count == 0) {
                 if (cs != "") list.push_back(cs);
                 cs.clear();
                 i++;
            }
            cs += cifstring[i];
            if (cifstring[i] == '(') count++;
            if (cifstring[i] == ')') count--;
       }
       if (cs != "") list.push_back(cs);
}

static void mult_matrixes(const double matrix1[4][4], const double matrix2[4][4], double result[4][4])
{
       int i, j, k;
       for (i = 0; i < 4; i++) {
            for (j = 0; j < 4; j++) {
                 result[i][j] = 0.0;
                 for (k = 0; k < 4; k++)
                      result[i][j] += (matrix1[i][k] * matrix2[k][j]);
            }
       }
}
#endif

static bool get_matrix_from_oper_list(const std::string& id_value, const std::vector<SymmMatrix>& SymmMatrices,
                                      const std::map<std::string, unsigned int>& SymmMatrix_Mapping, CrySymmetry& cell, double matrix[4][4])
{
       std::map<std::string, unsigned int>::const_iterator mpos = SymmMatrix_Mapping.find(id_value);
       if (mpos == SymmMatrix_Mapping.end()) return false;

       double matrix3[3][3], vect[3];
       bool return_flag = false;

       std::string sym_oper = SymmMatrices[mpos->second].getValue("name");
       std::string type = SymmMatrices[mpos->second].getValue("type");
       String::StripAndCompressWs(sym_oper);
       String::StripAndCompressWs(type);
       String::LowerCase(type);

       if (!cell.is_artifical() && type == "crystal symmetry operation" && !sym_oper.empty() && sym_oper != "1_555") {
            std::vector<std::string> str_list;
            get_wordarray(str_list, sym_oper, "_");
            if ((str_list.size() == 2) && (str_list[1].size() == 3)) {
                 cell.ndb_get_symmetry_matrix(sym_oper, matrix3, vect);
                 return_flag = true;
            }
       }
       if (!return_flag && SymmMatrices[mpos->second].getMatrix(matrix3) && SymmMatrices[mpos->second].getVector(vect)) return_flag = true;

       if (!return_flag) return return_flag;

       for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) matrix[i][j] = matrix3[i][j];
            matrix[i][3] = vect[i];
       }
       matrix[3][0] = matrix[3][1] = matrix[3][2] = 0;
       matrix[3][3] = 1;

       return return_flag;
}
