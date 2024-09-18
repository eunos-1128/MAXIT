/*
FILE:     Pdb2Ndb_Remark_350.C
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
#include "utillib.h"

void Maxit::_pdb_to_ndb_process_remark_350(const std::vector<std::string>& remark_array)
{
       if (remark_array.empty()) return;

       std::string method = "";

       std::map<int, std::vector<std::string> >::const_iterator mpos = _remarks.find(300);
       if (mpos != _remarks.end()) {
            for (std::vector<std::string>::const_iterator pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
                 if (pos->find("ICOSAHEDRAL POINT SYMMETRY") != std::string::npos)
                      method = "complete icosahedral assembly";
                 else if (pos->find("HELICAL SYMMETRY") != std::string::npos)
                      method = "representative helical assembly";
                 else if (pos->find("POINT SYMMETRY") != std::string::npos)
                      method = "complete point assembly";
            }
       }

       bool is_author_defined = false;
       bool is_software_defined = false;
       bool is_chain_id_line = false;

       _assemblies.clear();
       _SymmMatrices.clear();

       std::vector<SymmMatrix> matrices;
       matrices.clear();

       Assembly assembly;
       SymmMatrix symmmatrix;

       std::string cs;

       std::vector<std::string> data, chain_ids, matrix;
       chain_ids.clear();
       matrix.clear();

       std::string::size_type p;
       for (std::vector<std::string>::const_iterator pos = remark_array.begin(); pos != remark_array.end(); ++pos) {
            if (pos->empty()) continue;
            if ((p = pos->find("BIOMOLECULE:")) != std::string::npos) {
                 if (!chain_ids.empty() && !matrices.empty()) {
                      cs = _insert_symm_matrix(matrices);
                      assembly.InsertChains(cs, chain_ids);
                 }
                 matrices.clear();
                 chain_ids.clear();

                 if (!assembly.empty()) {
                      cs.clear();
                      if (is_author_defined && is_software_defined)
                           cs = "author_and_software_defined_assembly";
                      else if (is_author_defined)
                           cs = "author_defined_assembly";
                      else if (is_software_defined)
                           cs = "software_defined_assembly";
                      if (!cs.empty()) assembly.setValue("method", cs);
                      else if (!method.empty()) assembly.setValue("method", method);
                      assembly.Merge();
                      _assemblies.push_back(assembly);
                 }
                 is_author_defined = false;
                 is_software_defined = false;
                 assembly.clear();
                 get_wordarray(data, pos->substr(p + 12), " ");
                 if (!data.empty()) assembly.setValue("id", data[0]);
            } else if ((p = pos->find("AUTHOR DETERMINED BIOLOGICAL UNIT:")) != std::string::npos) {
                 is_author_defined = true;
                 cs = pos->substr(p + 34);
                 String::StripAndCompressWs(cs);
                 assembly.setValue("details", cs);
            } else if ((p = pos->find("QUATERNARY STRUCTURE FOR THIS ENTRY:")) != std::string::npos) {
                 is_author_defined = true;
                 cs = pos->substr(p + 36);
                 String::StripAndCompressWs(cs);
                 assembly.setValue("details", cs);
            } else if ((p = pos->find("SOFTWARE DETERMINED QUATERNARY STRUCTURE:")) != std::string::npos) {
                 is_software_defined = true;
                 cs = pos->substr(p + 41);
                 String::StripAndCompressWs(cs);
                 assembly.setValue("details", cs);
            } else if ((p = pos->find("SOFTWARE USED:")) != std::string::npos) {
                 is_software_defined = true;
                 cs = pos->substr(p + 14);
                 String::StripAndCompressWs(cs);
                 assembly.setValue("software", cs);
            } else if ((p = pos->find("AVERAGE BURIED SURFACE AREA:")) != std::string::npos) {
                 get_wordarray(data, pos->substr(p + 28), " ");
                 if (!data.empty()) assembly.setValue("absa", data[0]);
            } else if ((p = pos->find("BURIED SURFACE AREA FOR THE COMPLEX:"))
                                 != std::string::npos) {
                 get_wordarray(data, pos->substr(p + 36), " ");
                 if (!data.empty()) assembly.setValue("absa", data[0]);
            } else if ((p = pos->find("TOTAL BURIED SURFACE AREA:")) != std::string::npos) {
                 get_wordarray(data, pos->substr(p + 26), " ");
                 if (!data.empty()) assembly.setValue("absa", data[0]);
            } else if ((p = pos->find("TOTAL BURIED SURFACE AREA FOR THE COMPLEX:")) != std::string::npos) {
                 get_wordarray(data, pos->substr(p + 42), " ");
                 if (!data.empty()) assembly.setValue("absa", data[0]);
            } else if ((p = pos->find("TOTAL SURFACE AREA FOR THE COMPLEX:"))
                                 != std::string::npos) {
                 get_wordarray(data, pos->substr(p + 35), " ");
                 if (!data.empty()) assembly.setValue("ssa", data[0]);
            } else if ((p = pos->find("SURFACE AREA OF THE COMPLEX:")) != std::string::npos) {
                 get_wordarray(data, pos->substr(p + 28), " ");
                 if (!data.empty()) assembly.setValue("ssa", data[0]);
            } else if ((p = pos->find("SURFACE AREA FOR THE COMPLEX:")) != std::string::npos) {
                 get_wordarray(data, pos->substr(p + 29), " ");
                 if (!data.empty()) assembly.setValue("ssa", data[0]);
            } else if ((p = pos->find("GAIN IN SOLVENT FREE ENERGY:")) != std::string::npos) {
                 get_wordarray(data, pos->substr(p + 28), " ");
                 if (!data.empty()) assembly.setValue("more", data[0]);
            } else if ((p = pos->find("CHANGE IN SOLVENT FREE ENERGY:")) != std::string::npos) {
                 get_wordarray(data, pos->substr(p + 30), " ");
                 if (!data.empty()) assembly.setValue("more", data[0]);
            } else if ((p = pos->find("APPLY THE FOLLOWING TO CHAINS:")) != std::string::npos) {
                 if (!chain_ids.empty() && !matrices.empty()) {
                      cs = _insert_symm_matrix(matrices);
                      assembly.InsertChains(cs, chain_ids);
                 }

                 matrices.clear();
                 chain_ids.clear();

                 get_wordarray(chain_ids, pos->substr(p + 30), ";, ");
                 is_chain_id_line = true;
            } else if (pos->find("BIOMT1") != std::string::npos) {
                 is_chain_id_line = false;
                 get_wordarray(data, *pos, " ");
                 if (data.size() == 6) {
                      for (int i = 0; i < 4; ++i) matrix.push_back(data[i + 2]);
                 }
            } else if (pos->find("BIOMT2") != std::string::npos) {
                 is_chain_id_line = false;
                 get_wordarray(data, *pos, " ");
                 if (data.size() == 6) {
                      for (int i = 0; i < 4; ++i) matrix.push_back(data[i + 2]);
                 }
            } else if (pos->find("BIOMT3") != std::string::npos) {
                 is_chain_id_line = false;
                 get_wordarray(data, *pos, " ");
                 if (data.size() == 6) {
                      for (int i = 0; i < 4; ++i) matrix.push_back(data[i + 2]);
                 }
                 if (matrix.size() == 12) {
                      symmmatrix.clear();
                      for (unsigned int i = 0; i < matrix.size(); ++i) {
                           symmmatrix.setValue(matrix_items[i + 4], matrix[i]);
                      }
                      matrices.push_back(symmmatrix);
                 }
                 matrix.clear();
            } else if (is_chain_id_line && pos->find("NULL") == std::string::npos) {
                 cs = *pos;
                 p = pos->find("AND CHAINS:");
                 if (p != std::string::npos)
                      cs = pos->substr(p + 11);
                 get_wordarray(data, cs, ";, ");
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      chain_ids.push_back(*vpos);
                 }
            }
       }

       if (!chain_ids.empty() && !matrices.empty()) {
            cs = _insert_symm_matrix(matrices);
            assembly.InsertChains(cs, chain_ids);
       }

       if (!assembly.empty()) {
            cs.clear();
            if (is_author_defined && is_software_defined)
                 cs = "author_and_software_defined_assembly";
            else if (is_author_defined)
                 cs = "author_defined_assembly";
            else if (is_software_defined)
                 cs = "software_defined_assembly";
            if (!cs.empty()) assembly.setValue("method", cs);
            else if (!method.empty()) assembly.setValue("method", method);
            assembly.Merge();
            _assemblies.push_back(assembly);
       }
}

void Maxit::_pdb_to_ndb_process_pdbx_entry_details_remark(const std::vector<std::string>& remarkCards, const int& field_no, const std::string& keyword)
{
       std::string block_remark;
       _pdb_to_ndb_get_block_remark(remarkCards, keyword, block_remark);
       if (block_remark.empty()) return;

       _updateRecordFront("ENTDTL", field_no, block_remark);
}

void Maxit::_pdb_to_ndb_process_block_remark(const int& remark_no, const std::vector<std::string>& remarkCards)
{
       std::string block_remark;
       _pdb_to_ndb_get_block_remark(remarkCards, "", block_remark);
       if (block_remark.empty()) return;

       _addNewRecord("PDBRMK");
       _updateRecordBack("PDBRMK", 1, String::IntToString(remark_no));
       _updateRecordBack("PDBRMK", 3, block_remark);
}

void Maxit::_pdb_to_ndb_get_block_remark(const std::vector<std::string>& remarkCards, const std::string& keyword, std::string& block_remark)
{
       block_remark.clear();

       std::string cs;
       bool first = true;
       for (std::vector<std::string>::const_iterator pos = remarkCards.begin(); pos != remarkCards.end(); ++pos) {
            if (first && pos->empty()) continue;
            if (!keyword.empty()) {
                 String::UpperCase(*pos, cs);
                 String::StripAndCompressWs(cs);
                 if (keyword == cs) continue;
            }

            first = false;

            if (keyword == "COMPOUND" && pos->substr(0, 7) == " GROUP:") break;

            if (!block_remark.empty()) block_remark += "\n";
            block_remark += *pos;
       }

       String::StripTrailingWs(block_remark);
}

void Maxit::_pdb_to_ndb_process_remark_465(const std::vector<std::string>& remarkCards)
{
       std::map<std::string, std::vector<std::vector<std::string> > > missing_residue_mapping;
       missing_residue_mapping.clear();

       std::vector<std::vector<std::string> > data_vector;
       std::vector<std::string> data;

       bool start = false;
       for (std::vector<std::string>::const_iterator pos = remarkCards.begin(); pos != remarkCards.end(); ++pos) {
            if (pos->empty()) continue;
            if (pos->find("M RES C SSSEQI") != std::string::npos) {
                 start = true;
                 continue;
            }
            if (!start) continue;
            std::string name = pos->substr(4, 3);
            String::StripAndCompressWs(name);
            std::string chnid = pos->substr(8, 1);
            String::StripAndCompressWs(chnid);
            std::string num = pos->substr(10, 5);
            String::StripAndCompressWs(num);
            std::string ins_code = pos->substr(15, 1);
            String::StripAndCompressWs(ins_code);

            data.clear();
            data.push_back(name);
            data.push_back(num);
            data.push_back(ins_code);

            std::map<std::string, std::vector<std::vector<std::string> > >::iterator mpos = missing_residue_mapping.find(chnid);
            if (mpos != missing_residue_mapping.end())
                 mpos->second.push_back(data);
            else {
                 data_vector.clear();
                 data_vector.push_back(data);
                 missing_residue_mapping.insert(std::make_pair(chnid, data_vector));
            }
       }

       if (missing_residue_mapping.empty()) return;

       for (std::vector<RCSB::Molecule*>::iterator pos = _molecules.begin(); pos != _molecules.end(); ++pos) {
            RCSB::Chain* chain = (*pos)->GetFirstChain();
            while (chain) {
                 if (chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") {
                      std::map<std::string, std::vector<std::vector<std::string> > >::const_iterator mpos = missing_residue_mapping.find(chain->PDB_ChainID());
                      if (mpos != missing_residue_mapping.end())
                           chain->update_missing_residues(mpos->second);
                 }
                 chain = (*pos)->GetNextChain();
            }
       }
}

void Maxit::_pdb_to_ndb_process_900(const std::vector<std::string>& remarkCards)
{
       std::string db_id = "";
       std::string db_name = "";
       std::string details = "";
       for (std::vector<std::string>::const_iterator pos = remarkCards.begin(); pos != remarkCards.end(); ++pos) {
            std::string::size_type p = pos->find("RELATED ID:");
            std::string::size_type q = pos->find("RELATED DB:");
            if (p != std::string::npos && q != std::string::npos) {
                 if (!db_id.empty() && !db_name.empty()) {
                      String::StripAndCompressWs(details);
                      _addNewRecord("RENTRY");
                      _updateRecordBack("RENTRY", 2, db_name);
                      _updateRecordBack("RENTRY", 3, db_id);
                      _updateRecordBack("RENTRY", 4, details);
                 }
                 db_id = pos->substr(p + 11, q - p - 11);
                 String::StripAndCompressWs(db_id);
                 db_name = pos->substr(q + 11);
                 String::StripAndCompressWs(db_name); 
                 details.clear();
            } else {
                 if (!details.empty()) details += " ";
                 details += *pos;
            }
       }

       if (!db_id.empty() && !db_name.empty()) {
            String::StripAndCompressWs(details);
            _addNewRecord("RENTRY");
            _updateRecordBack("RENTRY", 2, db_name);
            _updateRecordBack("RENTRY", 3, db_id);
            _updateRecordBack("RENTRY", 4, details);
       }
}
