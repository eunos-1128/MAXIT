/*
FILE:     Maxit_Util.C
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

#include "Maxit.h"
#include "NdbToken.h"
#include "utillib.h"

Maxit::Maxit(): AnnotationObj()
{
       Maxit::clear();
}

Maxit::~Maxit()
{
       Maxit::clear();
}

void Maxit::clear()
{
       _remark_decimal_precision = FORCE_DECIMAL_PRECISION;
       _pdb_header_only_flag = false;
       _joint_methods.clear();
       _current_method.clear();
       _remarks.clear();
       _original_entry_ids.clear();
       _missing_polymer_residues.clear();
       _missing_polymer_atoms.clear();
       _zero_occ_polymer_residues.clear();
       _zero_occ_polymer_atoms.clear();
       _missing_non_polymer_residues.clear();
       _zero_occ_non_polymer_residues.clear();
       _non_polymer_prd_info.clear();
}

void Maxit::clear_pdb_records()
{
       _pdb_records.clear();
}

void Maxit::set_pdb_header_only_flag()
{
       _pdb_header_only_flag = true;
}

std::string Maxit::_get_refinement_program(bool& is_refmac5, bool& is_phenix, bool& is_cns_xplor, bool& is_buster)
{
       is_refmac5 = false;
       is_phenix = false;
       is_cns_xplor = false;
       is_buster = false;

       std::vector<std::string> refinement_programs;
       refinement_programs.clear();

       std::string cs;
       _getRecordFront("REFMET", 3, cs);
       if (!cs.empty()) {
            String::UpperCase(cs);
            refinement_programs.push_back(cs);
       }

       if (_CifObj) {
            Block& block = _CifObj->GetBlock(_firstBlockName);
            ISTable *t = getTablePtr(block, "computing");
            if (t) {
                 get_value_clean_upper(cs, t, 0, "structure_refinement");
                 if (!cs.empty()) refinement_programs.push_back(cs);
            }
            t = getTablePtr(block, "software");
            if (t) {
                 std::string cs1, cs2;
                 int rowNo = t->GetNumRows();
                 for (int i = 0; i < rowNo; ++i) {
                      get_value_clean_lower(cs1, t, i, "classification");
                      if (cs1 == "refinement") {
                           get_value_clean_upper(cs1, t, i, "name");
                           get_value_clean_upper(cs2, t, i, "version");
                           cs = cs1 + " " + cs2;
                           String::StripAndCompressWs(cs);
                           if (!cs.empty()) refinement_programs.push_back(cs);
                      }
                 }
            }
       }

       cs = _select_correct_refinement_program(refinement_programs);

       std::string cs1 = _get_known_refinement_name(cs);

       if (cs1.substr(0, 6) == "X-PLOR" || cs1.substr(0, 3) == "CNS" ||
           cs1.substr(0, 3) == "CNX") is_cns_xplor = true;
       if (cs1.substr(0, 8) == "REFMAC 5") is_refmac5 = true;
       if (cs1.substr(0, 6) == "PHENIX") is_phenix = true;
       if (cs1.substr(0, 6) == "BUSTER") is_buster = true;

       return cs;
}

int Maxit::_pdb_to_ndb_get_refinement_index(const std::string& program_name)
{
       std::string value1 = program_name;
       if (value1.empty()) _getRecordFront("REFMET", 3, value1);
       if (value1.empty()) return 0;

       String::UpperCase(value1);

       std::string value = _get_known_refinement_name(value1);

       if (value.substr(0, 3) == "CNS")
            return 1;
       else if (value.substr(0, 6) == "NUCLSQ")
            return 2;
       else if (value.substr(0, 6) == "PROLSQ" ||
                value.substr(0, 6) == "PROFFT" ||
                value.substr(0, 4) == "CCP4" ||
                value.substr(0, 7) == "PROTEIN" ||
                value.substr(0, 8) == "RESTRAIN" ||
                value.substr(0, 6) == "GPRLSA")
            return 3;
       else if (value.substr(0, 8) == "REFMAC 5")
            return 5;
       else if (value.substr(0, 6) == "REFMAC")
            return 4;
       else if (value.substr(0, 6) == "BUSTER")
            return 7;
       else if (value.substr(0, 3) == "TNT")
            return 6;
       else if (value.substr(0, 5) == "SHELX")
            return 8;
       else if (value.substr(0, 3) == "CNX")
            return 9;
       else if (value.substr(0, 6) == "PHENIX")
            return 10;
       else if (value.substr(0, 6) == "PRIMEX")
            return 11;

       return 0;
}

std::string Maxit::_get_known_refinement_name(const std::string& prog_name)
{
       std::vector<std::string> data;

       get_wordarray(data, prog_name, ",");
       for (std::vector<std::string>::const_iterator pos = data.begin(); pos != data.end(); ++pos) {
            if (pos->substr(0, 6) == "X-PLOR"   || pos->substr(0, 3) == "CNS"      ||
                pos->substr(0, 6) == "NUCLSQ"   || pos->substr(0, 6) == "PROLSQ"   ||
                pos->substr(0, 6) == "PROFFT"   || pos->substr(0, 4) == "CCP4"     ||
                pos->substr(0, 7) == "PROTEIN"  || pos->substr(0, 8) == "RESTRAIN" ||
                pos->substr(0, 6) == "GPRLSA"   || pos->substr(0, 8) == "REFMAC 5" ||
                pos->substr(0, 6) == "REFMAC"   || pos->substr(0, 6) == "BUSTER"   ||
                pos->substr(0, 3) == "TNT"      || pos->substr(0, 5) == "SHELX"    ||
                pos->substr(0, 3) == "CNX"      || pos->substr(0, 6) == "PHENIX"   ||
                pos->substr(0, 6) == "PRIMEX") return *pos;
       }
       return "";
}

std::string Maxit::_select_correct_refinement_program(const std::vector<std::string>& programs)
{
       if (programs.empty()) return ""; 

       for (std::vector<std::string>::const_iterator pos = programs.begin(); pos != programs.end(); ++pos) {
            if (pos->substr(0, 8) == "REFMAC 5" || pos->substr(0, 6) == "PHENIX") return *pos;
       }

       for (std::vector<std::string>::const_iterator pos = programs.begin(); pos != programs.end(); ++pos) {
            if (pos->substr(0, 6) == "X-PLOR"   || pos->substr(0, 3) == "CNS"      ||
                pos->substr(0, 6) == "NUCLSQ"   || pos->substr(0, 6) == "PROLSQ"   ||
                pos->substr(0, 6) == "PROFFT"   || pos->substr(0, 4) == "CCP4"     ||
                pos->substr(0, 7) == "PROTEIN"  || pos->substr(0, 8) == "RESTRAIN" ||
                pos->substr(0, 6) == "GPRLSA"   || pos->substr(0, 6) == "REFMAC"   ||
                pos->substr(0, 6) == "BUSTER"   || pos->substr(0, 3) == "TNT"      ||
                pos->substr(0, 5) == "SHELX"    || pos->substr(0, 3) == "CNX"      ||
                pos->substr(0, 6) == "PRIMEX") return *pos;
       }
 
       return programs[0];
}
