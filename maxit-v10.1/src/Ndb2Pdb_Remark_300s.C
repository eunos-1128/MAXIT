/*
FILE:     Ndb2Pdb_Remark_300s.C
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

#include "Assembly_Parser_Util.h"
#include "Maxit.h"
#include "NdbToken.h"
#include "Remark_define.h"
#include "Remark_extern.h"
#include "utillib.h"

void Maxit::_ndb_to_pdb_get_remark_300()
{
       std::string author_details = _ndb_to_pdb_get_struct_biol_details();
 
       if (_ndb_to_pdb_get_remark_300_split()) {
            if (!author_details.empty()) _get_author_defined_remark_300_details(author_details);
            return;
       }

       std::string biol_ids = _ndb_to_pdb_get_biol_ids();

       if (biol_ids.empty()) {
            if (!author_details.empty()) {
                 _ndb_to_pdb_get_general_remark(Num_Remark_300, Remark_300, 0, 0, 0);
                 _get_author_defined_remark_300_details(author_details);
            }
            return;
       }

       _updateRecordFront("CHAINS", 2, biol_ids);
       _ndb_to_pdb_get_general_remark(Num_Remark_300, Remark_300, 0, 0, 0);
       if (!author_details.empty()) _get_author_defined_remark_300_details(author_details);
       _pdb_records.erase("CHAINS");

       std::vector<std::string> remark_array;

       _ndb_to_pdb_get_remark_300_point_symmetry(remark_array);
       _ndb_to_pdb_add_remark(300, remark_array);

       _ndb_to_pdb_get_remark_300_helical_symmetry(remark_array);
       _ndb_to_pdb_add_remark(300, remark_array);
}

std::string Maxit::_ndb_to_pdb_get_struct_biol_details()
{
       if (!_CifObj) return "";

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);
       ISTable *t = getTablePtr(_cifblock, "struct_biol");
       if (!t) return "";
       
       std::string author_details, cs;
       author_details.clear();

       get_value_upper(cs, t, 0, "details");
       if (!cs.empty()) author_details = "REMARK: " + cs;

       return author_details;
}

bool Maxit::_ndb_to_pdb_get_remark_300_split()
{
       std::string value;
       _getRecordFront("SPLIT", 2, value);
       if (value.empty()) return false;

       _updateRecordFront("CHAINS", 2, "1");
       _ndb_to_pdb_get_general_remark(Num_Remark_300, Remark_300, 0, 0, 0);
       _pdb_records.erase("CHAINS");

       return true;
}

void Maxit::_get_author_defined_remark_300_details(const std::string& details)
{
       std::vector<std::string> remark_array;
       const ndb_token_format& ndbformat = NdbToken::getTokenFormat("REMARK");
       int width = ndbformat.FieldList[2].FieldWidth;

       get_text_array_from_block(remark_array, details, width);

       _ndb_to_pdb_add_remark(300, remark_array);
}

std::string Maxit::_ndb_to_pdb_get_biol_ids()
{
       if (_assemblies.empty()) return "";

       std::string biol_ids;
       biol_ids.clear();

       int biol_unit = 0;
       for (std::vector<Assembly>::const_iterator apos = _assemblies.begin(); apos != _assemblies.end(); ++apos) {
            std::string method = apos->getValue("method");
            if (method != "representative helical assembly" && method != "complete icosahedral assembly" && method != "complete point assembly" &&
                method != "author_and_software_defined_assembly" && method != "author_defined_assembly" && method != "software_defined_assembly") continue;

            biol_unit++;
            if (!biol_ids.empty()) biol_ids += ", ";
            biol_ids += String::IntToString(biol_unit);
       }
       return biol_ids;
}

void Maxit::_ndb_to_pdb_get_remark_300_point_symmetry(std::vector<std::string>& remark_array)
{
       remark_array.clear();

       if (!_CifObj) return;

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);
       ISTable *t = getTablePtr(_cifblock, "pdbx_point_symmetry");
       if (!t) return;

       std::string cs1, cs2;
       get_value(cs1, t, 0, "Schoenflies_symbol");
       get_value(cs2, t, 0, "circular_symmetry");
       if (cs1 == "I" || cs1 == "O" || cs1 == "T" || cs1 == "D" || cs1 == "C") {
            std::string remark = "THE ASSEMBLY REPRESENTED IN THIS ENTRY HAS REGULAR";
            remark_array.push_back(remark);

            if (cs1 == "I")
                 remark = "ICOSAHEDRAL POINT SYMMETRY (SCHOENFLIES SYMBOL = ";
            else if (cs1 == "O")
                 remark = "OCTAHEDRAL POINT SYMMETRY (SCHOENFLIES SYMBOL = ";
            else if (cs1 == "T")
                 remark = "TETRAHEDRAL POINT SYMMETRY (SCHOENFLIES SYMBOL = ";
            else if (cs1 == "D")
                 remark = "DIHEDRAL POINT SYMMETRY (SCHOENFLIES SYMBOL = ";
            else if (cs1 == "C")
                 remark = "CYCLIC POINT SYMMETRY (SCHOENFLIES SYMBOL = ";
            remark += cs1;
            if (cs1 == "D" || cs1 == "C") remark += cs2;
            remark += ").";
            remark_array.push_back(remark);
       }
}

void Maxit::_ndb_to_pdb_get_remark_300_helical_symmetry(std::vector<std::string>& remark_array)
{
       remark_array.clear();

       if (!_CifObj) return;

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);
       ISTable *t = getTablePtr(_cifblock, "pdbx_helical_symmetry");
       if (!t) return;

       std::string cs1, cs2, cs3, cs4, cs5, twist, height;
       get_value(cs1, t, 0, "rotation_per_n_subunits");
       get_value(cs2, t, 0, "rise_per_n_subunits");
       get_value(cs3, t, 0, "n_subunits_divisor");
       get_value_lower(cs4, t, 0, "dyad_axis");
       get_value(cs5, t, 0, "circular_symmetry");

       twist.clear();
       height.clear();
       if (atoi(cs3.c_str()) > 1) {
            twist = FloatToString(atof(cs1.c_str()), 0, 2) + "/" + cs3;
            height = FloatToString(atof(cs2.c_str()), 0, 2) + "/" + cs3; 
       } else {
            twist = FloatToString(atof(cs1.c_str()), 0, 2);
            height = FloatToString(atof(cs2.c_str()), 0, 2);
       }

       remark_array.push_back("THE ASSEMBLY REPRESENTED IN THIS ENTRY HAS REGULAR");
       remark_array.push_back("HELICAL SYMMETRY WITH THE FOLLOWING PARAMETERS:");

       std::string remark = "ROTATION PER SUBUNIT (TWIST) = " + twist + " DEGREES";
       remark_array.push_back(remark);
       remark = "RISE PER SUBUNIT (HEIGHT) = " + height + " ANGSTROMS";
       remark_array.push_back(remark);

       if (atoi(cs5.c_str()) > 1) {
            remark = "IN ADDITION, THERE IS " + cs5 + "-FOLD CIRCULAR";
            remark_array.push_back(remark);
            remark_array.push_back("SYMMETRY AROUND THE HELIX AXIS");
       }

       if (cs4 == "yes") {
            remark_array.push_back("IN ADDITION, THERE IS DYAD SYMMETRY");
            remark_array.push_back("PERPENDICULAR TO THE HELIX AXIS");
       }
}

void Maxit::_ndb_to_pdb_get_remark_350()
{
       if (_ndb_to_pdb_get_remark_350_split()) return;

       if (_assemblies.empty()) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_350 - 1, Remark_350, 0, 0, 0);

       std::vector<std::string> remark_array, id_list;

       std::set<std::string> linear_polymeric_chain_id_set;
       linear_polymeric_chain_id_set.clear();
       if (!_molecules.empty()) {
            RCSB::Chain* chain = _molecules[0]->GetFirstChain();
            while (chain) {
                 if ((chain->chain_type() == "ATOMP") || (chain->chain_type() == "ATOMN")) linear_polymeric_chain_id_set.insert(chain->PDB_ChainID());
                 chain = _molecules[0]->GetNextChain();
            }
       }

       double mat[4][4];
       int assembly_id = 0;
       for (std::vector<Assembly>::iterator apos = _assemblies.begin(); apos != _assemblies.end(); ++apos) {
            apos->UpdateOligomericStatus(linear_polymeric_chain_id_set);
            std::string method = apos->getValue("method");
            if (method != "representative helical assembly" && method != "complete icosahedral assembly" && method != "complete point assembly" &&
                method != "author_and_software_defined_assembly" && method != "author_defined_assembly" && method != "software_defined_assembly") continue;

            assembly_id++;
 
            remark_array.clear();
            remark_array.push_back("");
            remark_array.push_back("BIOMOLECULE: " + String::IntToString(assembly_id));

            std::string biol_unit = apos->getValue("details");
            String::UpperCase(biol_unit);
            if (biol_unit.empty()) {
                 biol_unit = apos->getValue("oligomeric");
                 int oligomeric_count = atoi(biol_unit.c_str());
                 if (oligomeric_count < 0) oligomeric_count = 0;
                 if (oligomeric_count <= NUM_MER)
                      biol_unit = _num_mers[oligomeric_count];
                 else biol_unit = String::IntToString(oligomeric_count) + "-MERIC";
            }
            std::string software = apos->getValue("software");
            if (software.empty()) software = "PISA";
            String::UpperCase(software);
            std::string absa = apos->getValue("absa");
            std::string ssa = apos->getValue("ssa");
            std::string more = apos->getValue("more");

            if (method == "author_and_software_defined_assembly" || method == "author_defined_assembly")
                 remark_array.push_back("AUTHOR DETERMINED BIOLOGICAL UNIT: " + biol_unit);

            if (method == "author_and_software_defined_assembly" || method == "software_defined_assembly") {
                 remark_array.push_back("SOFTWARE DETERMINED QUATERNARY STRUCTURE: " +  biol_unit);
                 remark_array.push_back("SOFTWARE USED: " + software);
                 if (!absa.empty())
                      remark_array.push_back("TOTAL BURIED SURFACE AREA: " + String::IntToString(atoi(absa.c_str())) + " ANGSTROM**2");
                 if (!ssa.empty())
                      remark_array.push_back("SURFACE AREA OF THE COMPLEX: " + String::IntToString(atoi(ssa.c_str())) + " ANGSTROM**2");
                 if (!more.empty())
                      remark_array.push_back("CHANGE IN SOLVENT FREE ENERGY: " + FloatToString(atof(more.c_str()), 0, 1) + " KCAL/MOL");
            }

            _ndb_to_pdb_add_remark(350, remark_array);

            int serialNo = 0;
            const std::vector<std::pair<std::string, std::vector<std::string> > >& assembly_chains = apos->assembly_chains();
            for (std::vector<std::pair<std::string, std::vector<std::string> > >::const_iterator
                 cpos = assembly_chains.begin(); cpos != assembly_chains.end(); ++cpos) {
                 _ndb_to_pdb_get_remark_350_APPLY_CHAINS(cpos->second);
                 remark_array.clear();
                 if (parseString(cpos->first, remark_array) == 0) {
                      for (std::vector<std::string>::const_iterator pos = remark_array.begin(); pos != remark_array.end(); ++pos) {
                           get_wordarray(id_list, *pos, " "); 
                           if (getFinalMatrix(mat, id_list, _SymmMatrices, _SymmMatrix_Mapping, _cell)) {
                                serialNo++;
                                _ndb_to_pdb_get_remark_350_Remark_SYMMA(serialNo, mat);
                           }
                      }
                 }
            }
       }
}

bool Maxit::_ndb_to_pdb_get_remark_350_split()
{
       std::string value;
       _getRecordFront("SPLIT", 2, value);
       if (value.empty()) return false;

       if (!_assemblies.empty()) return false;

       std::vector<std::string> remark_array;
       remark_array.clear();
       remark_array.push_back("BIOMOLECULE: 1");

       std::string cs;

       std::vector<std::string> id_list;
       id_list.clear();
       if (_CifObj) {
            Block &_cifblock = _CifObj->GetBlock(_firstBlockName);
            ISTable *t = getTablePtr(_cifblock, "struct_biol_gen");
            if (t) {
                 std::set<std::string> id_set;
                 id_set.clear();
     
                 int rowNo = t->GetNumRows();
                 for (int i = 0; i < rowNo; ++i) {
                      get_value(cs, t, i, "pdbx_new_pdb_asym_id");
                      if (cs.empty()) continue;
                      if (id_set.find(cs) != id_set.end()) continue;
                      id_set.insert(cs);
                      id_list.push_back(cs);
                 }
            }
       }

       if (id_list.size() <= NUM_MER) {
            cs = "QUATERNARY STRUCTURE FOR THIS ENTRY: ";
            cs += _num_mers[id_list.size()];
       } else cs = "QUATERNARY STRUCTURE FOR THIS ENTRY: "
                 + String::IntToString(id_list.size()) + "-MERIC";
       remark_array.push_back(cs);

       _ndb_to_pdb_add_remark(350, remark_array);

       _ndb_to_pdb_get_remark_350_APPLY_CHAINS(id_list);

       double mat[4][4];
       mat[0][0] = mat[1][1] = mat[2][2] = 1.0;
       mat[0][1] = mat[0][2] = mat[0][3] = 0.0;
       mat[1][0] = mat[1][2] = mat[1][3] = 0.0;
       mat[2][0] = mat[2][1] = mat[2][3] = 0.0;
       _ndb_to_pdb_get_remark_350_Remark_SYMMA(1, mat);

       return true;
}

void Maxit::_ndb_to_pdb_get_remark_350_APPLY_CHAINS(const std::vector<std::string> &id_list)
{
       std::string cs;
       cs.clear();
       for (unsigned int i = 0; i < id_list.size(); ++i) {
            if (!cs.empty()) cs += ", ";
            cs += id_list[i];
       }

       std::vector<std::string> data, remark_array;
       get_max_length_words(data, cs, 29);

       remark_array.clear();
       cs = "APPLY THE FOLLOWING TO CHAINS: ";
       for (std::vector<std::string>::iterator pos = data.begin(); pos != data.end(); ++pos) {
            remark_array.push_back(cs + *pos);
            cs = "                   AND CHAINS: ";
       }
       _ndb_to_pdb_add_remark(350, remark_array);
}

void Maxit::_ndb_to_pdb_get_remark_350_Remark_SYMMA(const int& serialNo, double mat[4][4])
{
       std::vector<std::string> remark_array;
       remark_array.clear();

       std::string cs;
       for (int k = 0; k < 3; ++k) {
            cs.clear();
            for (int j = 0; j < Remark_SYMMA.NumField; ++j) {
                 cs += get_space_between_remark_field(Remark_SYMMA, j);

                 bool left_adjust = false;
                 if (Remark_SYMMA.FieldList[j].FieldJustification == PDB_REMARK_LEFT_JUSTIFIED)
                      left_adjust = true;

                 int field_no = Remark_SYMMA.FieldList[j].FieldId-1;
                 if (!strcmp("BIOMT", Remark_SYMMA.FieldList[j].Text))
                      cs += Remark_SYMMA.FieldList[j].Text + String::IntToString(k + 1);
                 else if (field_no == 1)
                      cs += FloatToString((double) serialNo, Remark_SYMMA.FieldList[j].FieldWidth, Remark_SYMMA.FieldList[j].FieldPrec, left_adjust);
                 else if (j != Remark_SYMMA.NumField - 1)
                      cs += FloatToString(mat[k][j - 2], Remark_SYMMA.FieldList[j].FieldWidth, Remark_SYMMA.FieldList[j].FieldPrec, left_adjust);
                 else cs += FloatToString(mat[k][3], Remark_SYMMA.FieldList[j].FieldWidth, Remark_SYMMA.FieldList[j].FieldPrec, left_adjust);
            }
            remark_array.push_back(cs);
       }
       _ndb_to_pdb_add_remark(Remark_SYMMA.Remark_No, remark_array);
}

/*
void Maxit::_ndb_to_pdb_get_remark_375()
{
       _ndb_to_pdb_get_remark_375(_special_position_atoms);
}
*/

void Maxit::_ndb_to_pdb_get_remark_375(const std::list<Value>& special_position_atoms)
{
       if (special_position_atoms.empty()) return;

       _ndb_to_pdb_get_general_remark(Num_Remark_375, Remark_375, 0, 0, 0);
       for (std::list<Value>::const_iterator lpos = special_position_atoms.begin(); lpos != special_position_atoms.end(); ++lpos) {
            const std::vector<RCSB::Atom*>& atoms = lpos->atoms();
            std::string remark = "    ";
            if (atoms[0]->pdb_resnam() != "HOH")
                 remark = FormattedFieldValue(atoms[0]->pdb_atmnam(), 3, 4, 0, true, true);
            remark += " " + FormattedFieldValue(atoms[0]->pdb_resnam(), 3, 3, 0, false, true) + " ";
            remark += atoms[0]->pdb_chnid_char();
            remark += FormattedFieldValue(atoms[0]->pdb_resnum(), 3, 4, 0, false, true);
            remark += atoms[0]->ins_code_char();
            remark += " LIES ON A SPECIAL POSITION.";

            _addNewRemark(375, remark);
       }
}

void Maxit::_ndb_to_pdb_get_remark_375()
{
       std::list<Value> special_position_atoms = _get_special_position_atoms();
       _ndb_to_pdb_get_remark_375(special_position_atoms);
}
