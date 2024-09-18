/*
FILE:     Ndb2Pdb_HET_Util.C
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
#include <stdexcept>

#include "GenString.h"
#include "Maxit.h"
#include "SeqCodeUtil.h"

static void insert_formul_card(std::vector<std::pair<std::string, int> >& het_compNum_array, std::map<std::string, int>& het_counting_mapping,
                               const std::string &hetID, const int& serialNo, const int& num);

void Maxit::_ndb_to_pdb_get_HET_info()
{
       _pdb_records.erase("HET");
       _pdb_records.erase("HETNAM");
       _pdb_records.erase("HETSYN");
       _pdb_records.erase("FORMUL");

       if (_molecules.empty()) return;

       std::map<std::string, int> het_counting_mapping;
       het_counting_mapping.clear();
       std::vector<std::pair<std::string, int> > het_compNum_array;
       het_compNum_array.clear();

       int serialNo = 0;
       std::vector<RCSB::Residue*> residues;
       RCSB::Chain* chain = _molecules[0]->GetFirstChain();
       while (chain) {
            serialNo++;
            if (chain->chain_type() == "HETAS" && (chain->SeqRes(0)->Field[0] == "HOH" || chain->SeqRes(0)->Field[0] == "DOD")) {
                 if (chain->SeqRes(0)->Field[0] == "HOH")
                      insert_formul_card(het_compNum_array, het_counting_mapping, "HOH", serialNo, chain->ResidueNumbers());
                 else insert_formul_card(het_compNum_array, het_counting_mapping, "DOD", serialNo, chain->ResidueNumbers());
                 chain = _molecules[0]->GetNextChain();
                 continue;
            }

            for (unsigned int i = 0; i < chain->SeqLen(); ++i) {
                 if (chain->SeqRes(i)->ResIndex < 0) continue;
                 chain->GetResidueListByIndex(chain->SeqRes(i)->ResIndex, residues);
                 for (std::vector<RCSB::Residue*>::const_iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                      if (!SeqCodeUtil::is_a_standard_residue((*rpos)->ResName()) || (SeqCodeUtil::is_a_standard_residue((*rpos)->ResName()) &&
                         ((chain->chain_type() != "ATOMP" && chain->chain_type() != "ATOMN") || chain->SeqLen() == 1))) {
                           insert_formul_card(het_compNum_array, het_counting_mapping, (*rpos)->ResName(), serialNo, 1);
                           if ((*rpos)->ResName() != "DUM" && (*rpos)->ResName() != "UNK" && (*rpos)->ResName() != "N") {
                                 _addNewRecord("HET");
                                 _updateRecordBack("HET", 1, (*rpos)->ResName());
                                 _updateRecordBack("HET", 2, (*rpos)->pdb_chnid());
                                 _updateRecordBack("HET", 3, (*rpos)->pdb_res_no());
                                 _updateRecordBack("HET", 4, (*rpos)->ins_code());
                                 _updateRecordBack("HET", 5, String::IntToString((*rpos)->AtomNumbers()));
                           }
                      }
                 }
            }
            chain = _molecules[0]->GetNextChain();
       }

       if (het_compNum_array.empty()) return;

       for (std::vector<std::pair<std::string, int> >::const_iterator pos = het_compNum_array.begin(); pos != het_compNum_array.end(); ++pos) {
            if (pos->first == "UNL") {
                 _addNewRecord("HETNAM");
                 _updateRecordBack("HETNAM", 2, pos->first);
                 _updateRecordBack("HETNAM", 3, "UNKNOWN LIGAND");
                 continue;
            }

            _addNewRecord("FORMUL");
            _updateRecordBack("FORMUL", 1, String::IntToString(pos->second));
            _updateRecordBack("FORMUL", 2, pos->first);
            std::map<std::string, int>::const_iterator
                mpos = het_counting_mapping.find(pos->first);

            if (pos->first == "HOH" || pos->first == "DOD") {
                 _updateRecordBack("FORMUL", 4, "*");
                 std::string formula = "(H2 O)";
                 if (pos->first == "DOD") formula = "(D2 O)";
                 if (mpos != het_counting_mapping.end() && mpos->second > 1)
                      _updateRecordBack("FORMUL", 5, String::IntToString(mpos->second) + formula); 
                 else _updateRecordBack("FORMUL", 5, formula);
                 continue;
            }

            try {
                 const ConnectFormat& drug = _ccDic->find_drug(pos->first);
                 if (!drug.chemical_name().empty()) {
                      _addNewRecord("HETNAM");
                      _updateRecordBack("HETNAM", 2, pos->first);
                      _updateRecordBack("HETNAM", 3, drug.chemical_name());
                 }
                 if (!drug.synonym().empty()) {
                      _addNewRecord("HETSYN");
                      _updateRecordBack("HETSYN", 2, pos->first);
                      _updateRecordBack("HETSYN", 3, drug.synonym());
                 }

                 std::string formula = drug.formula();
                 std::string charge = drug.getMetaData("pdbx_formal_charge");
                 if (!charge.empty() && charge != "0") {
                      int chg = abs(atoi(charge.c_str()));
                      formula += " " + String::IntToString(chg);
                      if (charge.find("-") != std::string::npos)
                           formula += "-";
                      else formula += "+";
                 }
                 if (mpos != het_counting_mapping.end() && mpos->second > 1)
                      _updateRecordBack("FORMUL", 5, String::IntToString(mpos->second) + "(" + formula + ")"); 
                 else _updateRecordBack("FORMUL", 5, formula);
            } catch (const std::exception& exc) {}
       }
}

static void insert_formul_card(std::vector<std::pair<std::string, int> >& het_compNum_array, std::map<std::string, int>& het_counting_mapping,
                               const std::string &hetID, const int& serialNo, const int& num)
{
       if (hetID == "DUM" || hetID == "UNK" || hetID == "N") return;

       std::map<std::string, int>::iterator mpos = het_counting_mapping.find(hetID);
       if (mpos != het_counting_mapping.end()) {
            mpos->second += num;
            return;
       }

       het_counting_mapping.insert(std::make_pair(hetID, num));

       het_compNum_array.push_back(std::make_pair(hetID, serialNo));
}
