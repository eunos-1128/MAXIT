/*
FILE:     SeqCodeUtil.C
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
#include <unistd.h>
#include <sys/stat.h>

#include "SeqCodeUtil.h"
#include "SeqCodeUtil_global.h"
#include "utillib.h"

void SeqCodeUtil::initialize()
{
       _a1_a3_mapping.clear();
       _a3_a1_mapping.clear();
       _standard_residue_list_set1.clear();
       _standard_residue_list_set2.clear();
       _standard_aa_residue_list_set.clear();
       _standard_Laa_residue_list_set.clear();
       _standard_Daa_residue_list_set.clear();
       _standard_na_residue_list_set.clear();
       _standard_dna_residue_list_set.clear();
       _standard_rna_residue_list_set.clear();
       _water_residue_set.clear();
       _d1_d2_mapping.clear();
       _similar_residue_set.clear();
       _special_na_residue_list_set.clear();
       _special_na_residue_mapping.clear();

       for (int i = 0; i < AANUMBER; ++i) {
            std::string cs = aacodes[i].a1let; cs += "_"; cs += aacodes[i].type;
            _a1_a3_mapping.insert(std::make_pair(cs, aacodes[i].a3let));
            _a3_a1_mapping.insert(std::make_pair(aacodes[i].a3let, aacodes[i].a1let));
       }
  
       for (int i = 0; i < NUM_DNA_RESIDUE; ++i) {
            _standard_residue_list_set1.insert(_dna_residues[i]);
            _standard_residue_list_set2.insert(_dna_residues[i]);
            _standard_na_residue_list_set.insert(_dna_residues[i]);
            _standard_dna_residue_list_set.insert(_dna_residues[i]);
            _d1_d2_mapping.insert(std::make_pair(_dna_residues[i], _s_dna_residues[i]));
            _d1_d2_mapping.insert(std::make_pair(_s_dna_residues[i], _dna_residues[i]));
       }
  
       for (int i = 0; i < NUM_RNA_RESIDUE; ++i) {
            _standard_residue_list_set1.insert(_rna_residues[i]);
            _standard_residue_list_set2.insert(_rna_residues[i]);
            _standard_na_residue_list_set.insert(_rna_residues[i]);
            _standard_rna_residue_list_set.insert(_rna_residues[i]);
       }

       for (int i = 0; i < NUM_AA_RESIDUE; ++i) {
            _standard_residue_list_set1.insert(_aa_residues[i]);
            _standard_residue_list_set2.insert(_aa_residues[i]);
            _standard_aa_residue_list_set.insert(_aa_residues[i]);
            _standard_Laa_residue_list_set.insert(_aa_residues[i]);
       }

       for (int i = 0; i < NUM_DAA_RESIDUE; ++i) {
            _standard_aa_residue_list_set.insert(_daa_residues[i]);
            _standard_Daa_residue_list_set.insert(_daa_residues[i]);
       }

       for (int i = 0; i < NUM_MIXED_RESIDUE; ++i) {
            _standard_residue_list_set1.insert(_mixed_residues[i]);
            _standard_residue_list_set2.insert(_mixed_residues[i]);
            _standard_Laa_residue_list_set.insert(_mixed_residues[i]);
       }

       for (int i = 0; i < NUM_UNKNOWN_RESIDUE; ++i) {
            _standard_residue_list_set1.insert(_unknown_residues[i]);
            _standard_residue_list_set2.insert(_unknown_residues[i]);
       }

       for (int i = 0; i < NUM_WATER_RESIDUE; ++i) {
            _standard_residue_list_set2.insert(_water_residues[i]);
            _water_residue_set.insert(_water_residues[i]);
       }

       for (int i = 0; i < NUM_NN_RESIDUE; ++i) {
            _standard_residue_list_set2.insert(_nn_residues[i]);
       }

       for (int i = 0; i < NUM_SIMILAR_RESIDUE; ++i) {
            std::string cs = _similar_residues[i][0]; cs += "_"; cs += _similar_residues[i][1];
            _similar_residue_set.insert(cs);
            cs = _similar_residues[i][1]; cs += "_"; cs += _similar_residues[i][0];
            _similar_residue_set.insert(cs);
       }

       for (int i = 0; i < NUM_SPECIAL_NN_RESIDUE; ++i) {
            _special_na_residue_list_set.insert(_special_nn_residues[i]);
       }

       for (int i = 0; i < NUM_MAPPING_NN_RESIDUE; ++i) {
            _special_na_residue_mapping.insert(std::make_pair(_special_nn_residues[i], _mapping_nn_residues[i]));
       }
}

std::string SeqCodeUtil::GetOneLetterCode(const std::string& code)
{
       std::string cs;
       String::UpperCase(code, cs);
       std::map<std::string, std::string>::const_iterator pos = _a3_a1_mapping.find(cs);
       if (pos != _a3_a1_mapping.end()) return pos->second;
       return cs; 
}

std::string SeqCodeUtil::GetOneLetterCodeWithX(const std::string& code)
{
       if (code.size() <= 1) return code;
       else {
            std::string code1 = GetOneLetterCode(code);
            if (code1.size() == 1) return code1;
       }
       return "X";
}

std::string SeqCodeUtil::GetOneLetterCodeWithMSE(const std::string& code)
{
       if (code.size() <= 1) return code;
       else if (String::IsEqual(code, "MSE", Char::eCASE_INSENSITIVE))
            return "M";
       else {
            std::string code1 = GetOneLetterCode(code);
            if (code1.size() == 1) return code1;
       }
       return "X";
}

std::string SeqCodeUtil::GetThreeLetterCode(const char code, const std::string& type, bool& is_standard_flag)
{
       is_standard_flag = false;
       std::string cs, typeU;
       cs.clear(); cs += code;
       String::UpperCase(cs);
       std::string cs1 = cs + "_" + type;
       String::UpperCase(cs1);
       std::map<std::string, std::string>::const_iterator pos = _a1_a3_mapping.find(cs1);
       if (pos != _a1_a3_mapping.end()) {
            is_standard_flag = true;
            return pos->second;
       }
       String::UpperCase(type, typeU);
       if (typeU == "POLYPEPTIDE(D)") {
            cs1 = cs + "_polypeptide(L)";
            String::UpperCase(cs1);
            pos = _a1_a3_mapping.find(cs1);
            if (pos != _a1_a3_mapping.end()) {
                 is_standard_flag = true;
                 return pos->second;
            }
       }
       return cs;
}

std::string SeqCodeUtil::GetOneLetterCodeSeq(const std::vector<std::string>& seqs, const int length)
{
       std::string a1seq;
       a1seq.clear();
       int count = 0;
       for (std::vector<std::string>::const_iterator
            pos = seqs.begin(); pos != seqs.end(); ++pos) {
            std::string cs = GetOneLetterCode(*pos);
            if (cs.size() > 1) {
                 a1seq += "(" + cs + ")";
                 count += cs.size() + 2;
            } else {
                 a1seq += cs;
                 count += cs.size();
            }
            if (length > 0 && count >= length) {
                 count = 0;
                 a1seq += "\n";
            }
       }
       return a1seq;
}

std::string SeqCodeUtil::GetStandardNaCode(const std::string& special_code)
{
       std::string cs;
       String::UpperCase(special_code, cs);
       std::map<std::string, std::string>::const_iterator
           pos = _special_na_residue_mapping.find(cs);
       if (pos != _special_na_residue_mapping.end()) return pos->second;
       return ""; 
}

std::string SeqCodeUtil::GetOneLetterCodeSeqArray(const std::string& a1seq, std::vector<std::string>& one_seqs)
{
       one_seqs.clear();

       std::list<std::string> seq_list;
       std::string err = GetOneLetterCodeSeqList(a1seq, seq_list);

       one_seqs.reserve(seq_list.size());
       for (std::list<std::string>::const_iterator spos = seq_list.begin(); spos != seq_list.end(); ++spos) {
            one_seqs.push_back(*spos);
       }

       return err;
}

std::string SeqCodeUtil::GetOneLetterCodeSeqList(const std::string& a1seq, std::list<std::string>& one_seqs)
{
       one_seqs.clear();

       std::set<std::string> missing_set;
       std::vector<std::string> missing_vec;
       std::string cs;
       missing_set.clear();
       missing_vec.clear();

       for (unsigned int j = 0; j < a1seq.size(); ++j) {
            if ((a1seq[j] == ' ') || (a1seq[j] == '\t') || (a1seq[j] == '\n') || (a1seq[j] == '\r')) continue;
            if (a1seq[j] != '(') {
                 if (a1seq[j] == ')') {
                      if (missing_set.find("(") == missing_set.end()) {
                           missing_set.insert("(");
                           missing_vec.push_back("(");
                      }
                      continue;
                 }
                 cs.clear();
                 cs += a1seq[j];
                 one_seqs.push_back(cs);
            } else {
                 if (j == a1seq.size() - 1) {
                      if (missing_set.find(")") == missing_set.end()) {
                           missing_set.insert(")");
                           missing_vec.push_back(")");
                      }
                      continue;
                 }

                 ++j;
                 cs.clear();
                 while (a1seq[j] != ')') {
                      if (j == a1seq.size() - 1) {
                           if (missing_set.find(")") == missing_set.end()) {
                                missing_set.insert(")");
                                missing_vec.push_back(")");
                           }
                           break;
                      } else if (a1seq[j] == '(') {
                           if (missing_set.find(")") == missing_set.end()) {
                                missing_set.insert(")");
                                missing_vec.push_back(")");
                           }
                      }
                      if (a1seq[j] != ' ' && a1seq[j] != '\t' && a1seq[j] != '\n' && a1seq[j] != '(') cs += a1seq[j];
                      ++j;
                 }
                 one_seqs.push_back(cs);
            }
       }

       return _getParenthesesError(missing_vec);
}

std::string SeqCodeUtil::GetThreeLetterCodeSeq(const std::string& a1seq, const std::string& type, std::vector<std::string>& seqs)
{
       seqs.clear();
       std::vector<std::string> one_seqs;
       std::string error = GetOneLetterCodeSeqArray(a1seq, one_seqs);

       std::string poly_type = type;
       if (poly_type.empty()) {
            for (std::vector<std::string>::const_iterator
                 pos = one_seqs.begin(); pos != one_seqs.end(); ++pos) {
                 if (pos->size() == 1 && *pos != "A" && *pos != "C" && *pos != "G" && *pos != "I" && *pos != "T" && *pos != "U") {
                      poly_type = "polypeptide(L)";
                      break;
                 }
            }
            if (poly_type.empty()) poly_type = "polyribonucleotide";
       }
 
       bool is_standard_flag = false;
       seqs.reserve(one_seqs.size());
       for (std::vector<std::string>::const_iterator pos = one_seqs.begin(); pos != one_seqs.end(); ++pos) {
            if (pos->size() == 1) {
                 std::string cs = GetThreeLetterCode((*pos)[0], poly_type, is_standard_flag);
                 seqs.push_back(cs);
            } else seqs.push_back(*pos);
       }
       return error;
}

std::string SeqCodeUtil::GetOneAndThreeLetterSeq(const std::string& a1seq, const std::string& type, std::vector<std::vector<std::string> >& seqs)
{
       seqs.clear();
       std::vector<std::string> one_seqs, tmp_seq;
       std::string error = GetOneLetterCodeSeqArray(a1seq, one_seqs);

       std::string poly_type = type;
       if (poly_type.empty()) {
            for (std::vector<std::string>::const_iterator
                 pos = one_seqs.begin(); pos != one_seqs.end(); ++pos) {
                 if (pos->size() == 1 && *pos != "A" && *pos != "C" && *pos != "G" && *pos != "I" && *pos != "T" && *pos != "U") {
                      poly_type = "polypeptide(L)";
                      break;
                 }
            }
            if (poly_type.empty()) poly_type = "polyribonucleotide";
       }
 
       bool is_standard_flag = false;
       seqs.reserve(one_seqs.size());
       for (std::vector<std::string>::const_iterator pos = one_seqs.begin(); pos != one_seqs.end(); ++pos) {
            tmp_seq.clear();
            if (pos->size() == 1) {
                 std::string cs = GetThreeLetterCode((*pos)[0], poly_type, is_standard_flag);
                 tmp_seq.push_back(*pos);
                 tmp_seq.push_back(cs);
            } else {
                 std::string code1 = GetOneLetterCodeWithX(*pos);
                 tmp_seq.push_back(code1);
                 tmp_seq.push_back(*pos);
            }
            seqs.push_back(tmp_seq);
       }
       return error;
}

bool SeqCodeUtil::is_a_standard_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_standard_residue_list_set1.find(cs) != _standard_residue_list_set1.end())
            return true;
       else return false;
}

bool SeqCodeUtil::is_a_standard_residue_plus_MSE(const std::string& resname)
{
       if (String::IsEqual(resname, "MSE", Char::eCASE_INSENSITIVE)) return true;
       return is_a_standard_residue(resname);
}

bool SeqCodeUtil::is_a_nonstandard_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_standard_residue_list_set2.find(cs) != _standard_residue_list_set2.end())
            return false;
       else return true;
}

bool SeqCodeUtil::is_standard_aa_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_standard_Laa_residue_list_set.find(cs) != _standard_Laa_residue_list_set.end())
            return true;
       else return false;
}

bool SeqCodeUtil::is_standard_aa_residue_plus_MSE(const std::string& resname)
{
       if (String::IsEqual(resname, "MSE", Char::eCASE_INSENSITIVE)) return true;
       return is_standard_aa_residue(resname);
}

bool SeqCodeUtil::is_standard_Daa_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_standard_Daa_residue_list_set.find(cs) != _standard_Daa_residue_list_set.end())
            return true;
       else return false;
}

bool SeqCodeUtil::is_aa_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_standard_aa_residue_list_set.find(cs) != _standard_aa_residue_list_set.end())
            return true;
       else return false;
}

bool SeqCodeUtil::is_aa_residue_plus_MSE(const std::string& resname)
{
       if (String::IsEqual(resname, "MSE", Char::eCASE_INSENSITIVE)) return true;
       return is_aa_residue(resname);
}

bool SeqCodeUtil::is_dna_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_standard_dna_residue_list_set.find(cs) != _standard_dna_residue_list_set.end())
            return true;
       else return false;
}

bool SeqCodeUtil::is_rna_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_standard_rna_residue_list_set.find(cs) != _standard_rna_residue_list_set.end())
            return true;
       else return false;
}

bool SeqCodeUtil::is_na_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_standard_na_residue_list_set.find(cs) != _standard_na_residue_list_set.end())
            return true;
       else return false;
}

bool SeqCodeUtil::is_special_na_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_special_na_residue_list_set.find(cs) != _special_na_residue_list_set.end())
            return true;
       else return false;
}

bool SeqCodeUtil::is_water_residue(const std::string& resname)
{
       std::string cs;
       String::UpperCase(resname, cs);
       if (_water_residue_set.find(cs) != _water_residue_set.end())
            return true;
       else return false;
}

bool SeqCodeUtil::isSameResName(const std::string& resname1, const std::string& resname2)
{
       std::string cs1, cs2;
       String::UpperCase(resname1, cs1);
       String::UpperCase(resname2, cs2);
       if (cs1 == cs2) return true;

       std::map<std::string, std::string>::const_iterator pos = _d1_d2_mapping.find(cs1);
       if (pos != _d1_d2_mapping.end() && pos->second == cs2) return true;
       return false;
}

bool SeqCodeUtil::isSimilarResidue(const std::string& resname1, const std::string& resname2)
{
       std::string cs = resname1 + "_" + resname2;
       if (_similar_residue_set.find(cs) != _similar_residue_set.end())
            return true;
       else return false;
}

std::string SeqCodeUtil::_getParenthesesError(const std::vector<std::string>& vec)
{
       if (vec.empty()) return "";

       std::string error = "";
       for (std::vector<std::string>::const_iterator pos = vec.begin(); pos != vec.end(); ++pos) {
            if (!error.empty()) error += " and ";
            if ((*pos) == "(")
                 error += "open parenthesis '('";
            else if ((*pos) == ")")
                 error += "close parenthesis ')'";
       }
       if (error.empty()) return "";
       return "missing " + error;
}
