/*
FILE:     Element.C
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
#include <ctype.h>

#include "Element.h"
#include "Element_global.h"
#include "GenString.h"

void Element::initialize()
{
       _element_mapping.clear();
       for (int i = 0; i < NUM_PERIODIC_TABLE_ATOM; ++i) {
            std::string cs = periodic_table[i].Atom_Symbol;
            String::UpperCase(cs);
            _element_mapping.insert(std::make_pair(cs, i));
       }
}

void Element::getAllowElement(std::set<std::string>& e_set)
{
       e_set.clear();
       e_set.insert("D");

       for (std::map<std::string, int>::iterator pos = _element_mapping.begin(); pos != _element_mapping.end(); ++pos) {
            if (pos->first == "??") continue;
            e_set.insert(pos->first);
       }
}

std::string Element::getAtomSymbol(const std::string& atomName)
{
       if (atomName.empty()) return "";

       std::string atomType, tmpAtom, tmpAtom2;

       if (!isalnum(atomName[0]))
            tmpAtom = atomName.substr(1);
       else tmpAtom = atomName;

       String::StripAndCompressWs(tmpAtom);
       String::UpperCase(tmpAtom);
       if (tmpAtom.size() > 1 && isdigit(tmpAtom[0]) && tmpAtom[1] == 'H') {
            atomType = "H";
       }
       for (unsigned int i = 0; i < tmpAtom.size(); ++i) {
            if (!isdigit(tmpAtom[i]) && tmpAtom[i] != 'X' && tmpAtom[i] != ' ') {
                 tmpAtom2 = tmpAtom.substr(i, 2);
                 break;
            }
       }
    
       if (tmpAtom2.size() > 1) {
            if (isdigit(tmpAtom2[1])) tmpAtom2.erase(1);
            else if (_element_mapping.find(tmpAtom2) == _element_mapping.end() && tmpAtom2 != "Q" && tmpAtom2 != "D") {
                 tmpAtom = tmpAtom2.substr(1, 1); 
                 tmpAtom2.erase(1);
                 if (_element_mapping.find(tmpAtom2) == _element_mapping.end() && tmpAtom2 != "Q" && tmpAtom2 != "D") tmpAtom2 = tmpAtom;
            }
       }
       atomType = tmpAtom2;
       return atomType;
}

void Element::getAtomSymbols(std::map<std::string, std::string>& Name_Type)
{
       if (Name_Type.empty()) return;

       std::map<std::string, std::string> name_type;
       name_type.clear();
       std::string type;
       for (std::map<std::string, std::string>::iterator pos = Name_Type.begin(); pos != Name_Type.end(); ++pos) {
            std::string type = getAtomSymbol(pos->first);
            if ((type == "CA" || type == "CD" || type == "CE") && Name_Type.size() > 2)
                 type = "C";
            else if (type == "CM")
                 type = "C";
            else if (type == "NA" && Name_Type.size() > 2)
                 type = "N";
            else if (type == "NB" || type == "ND" || type == "NE")
                 type = "N";
            else if (type == "PA" || type == "PB" || type == "PD" || type == "PO")
                 type = "P";
            else if (type == "SG")
                 type = "S";
            pos->second = type;
       }

       _refineAtomType(Name_Type);
}

int Element::getAtomNumber(const std::string& atomtype)
{
       std::string cs;
       String::UpperCase(atomtype, cs);
       std::map<std::string, int>::iterator pos = _element_mapping.find(cs);
       if (pos != _element_mapping.end()) return (pos->second + 1);
       return (-1);
}

double Element::getIonicRadii(const int atom_number)
{
       if (atom_number >= 1 && atom_number <= NUM_PERIODIC_TABLE_ATOM)
            return periodic_table[atom_number-1].Ionic_Radii;
       else return 0;
}

int Element::findMetalFlag(const std::string& atomType)
{
       std::string cs;
       String::UpperCase(atomType, cs);
       std::map<std::string, int>::iterator pos = _element_mapping.find(cs);
       if (pos != _element_mapping.end()) return (periodic_table[pos->second].is_metal_element);
       return 0;
}

double Element::getWeight(const std::string& atomType)
{
       std::string cs;
       String::UpperCase(atomType, cs);
       std::map<std::string, int>::iterator pos = _element_mapping.find(cs);
       if (pos != _element_mapping.end()) return (periodic_table[pos->second].Atom_Weight);
       return 0;
}

int Element::getAtomAlignPosition(const std::string& atomName, const std::string& atomType)
{
       int align = 1;
       std::string::size_type p = atomName.find(atomType);
       if (p == 1 || atomType.size() == 2 || (atomType.size() == 1 && atomName.size() > 1 && atomName[1] == atomType[0] &&
           atomType[0] != 'D' && atomType[0] != 'H' && atomType[0] != 'C' && atomType[0] != 'N'))
            align = 0;
       return align;
}

void Element::_refineAtomType(std::map<std::string, std::string>& name_type)
{
       if (name_type.empty()) return;

       int n_n = 0, n_c = 0, n_ac = 0, n_op = 0;
       for (std::map<std::string, std::string>::iterator mpos = name_type.begin(); mpos != name_type.end(); ++mpos) {
            if (mpos->second == "CA" || mpos->second == "CD" || mpos->second == "CE") n_c++;
            if (mpos->second == "AC") n_ac++;
            if (mpos->second == "NP" || mpos->second == "NO") n_op++;
            if (mpos->second == "N") {
                 std::string atom_name = mpos->first;
                 std::string::size_type p = atom_name.find(mpos->second);
                 if (p != std::string::npos && p < atom_name.size() - 1) {
                      std::string type = atom_name.substr(p + 1);
                      if (type.size() > 1) type.erase(1);
                      if (isalpha(type[0]) && getAtomNumber(type) > 0) n_n++;
                 }
            }
       }
       if (n_c < 2 && n_n < 4 && n_ac < 4 && n_op < 4) return;

       if (n_c >= 2) {
            for (std::map<std::string, std::string>::iterator mpos = name_type.begin(); mpos != name_type.end(); ++mpos) {
                 if (mpos->second == "CA" || mpos->second == "CD" || mpos->second == "CE") mpos->second = "C";
            }
       }
       if (n_n >= 4) {
            for (std::map<std::string, std::string>::iterator mpos = name_type.begin(); mpos != name_type.end(); ++mpos) {
                 if (mpos->second != "N") continue;
                 std::string atom_name = mpos->first;
                 std::string::size_type p = atom_name.find(mpos->second);
                 if (p != std::string::npos && p < atom_name.size() - 1) {
                      std::string type = atom_name.substr(p + 1);
                      if (type.size() > 1) type.erase(1);
                      if (isalpha(type[0]) && getAtomNumber(type) > 0) {
                           mpos->second = type;
                      }
                 }
            }
       }
       if (n_ac >= 4) {
            for (std::map<std::string, std::string>::iterator mpos = name_type.begin(); mpos != name_type.end(); ++mpos) {
                 if (mpos->second == "AC") mpos->second = "C";
            }
       }
       if (n_op >= 4) {
            for (std::map<std::string, std::string>::iterator mpos = name_type.begin(); mpos != name_type.end(); ++mpos) {
                 if (mpos->second == "NP" || mpos->second == "NO") {
                      mpos->second = mpos->second.substr(1);
                 }
            }
       }
}
