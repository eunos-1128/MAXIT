/*
FILE:     BondUtil.C
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
#include <math.h>

#include "BondUtil.h"
#include "BondUtil_global.h"
#include "Element.h"
#include "GenString.h"

void BondUtil::initialize()
{
       _bond_mapping.clear();
       for (int i = 0; i < NUM_BOND; i++) {
            _bond_mapping.insert(std::make_pair(_allowed_bonds[i].label, i));
       }
}

std::string BondUtil::_get_label(const std::string& _type1, const std::string& _type2)
{
       std::string type1 = _type1;
       if (_type1 == "D") type1 = "H";
       std::string type2 = _type2;
       if (_type2 == "D") type2 = "H";

       if (type1 < type2)
            return (type1 + "_" + type2);
       else return (type2 + "_" + type1);
}

_ALLOWED_BOND* BondUtil::_find_bond(const std::string& label)
{
       std::string cs;
       String::UpperCase(label, cs);
       std::map<std::string, int>::iterator pos = _bond_mapping.find(cs);
       if (pos != _bond_mapping.end())
            return &_allowed_bonds[pos->second];
       else return NULL;
}

int BondUtil::is_a_bond(RCSB::Atom* atom1, RCSB::Atom* atom2)
{
       std::string type1 = atom1->atom_type();
       std::string type2 = atom2->atom_type();
       if (type1.empty()) type1 = Element::getAtomSymbol(atom1->atmtype());
       if (type2.empty()) type2 = Element::getAtomSymbol(atom2->atmtype());

       return (is_a_bond(atom1, type1, atom2, type2));
}

int BondUtil::is_a_bond_with_range(RCSB::Atom* atom1, RCSB::Atom* atom2, double& dist, double& lower_limit, double& upper_limit)
{
       std::string type1 = atom1->atom_type();
       std::string type2 = atom2->atom_type();
       if (type1.empty()) type1 = Element::getAtomSymbol(atom1->atmtype());
       if (type2.empty()) type2 = Element::getAtomSymbol(atom2->atmtype());

       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            dist = cal_distance(atom1, atom2);
            if (dist >= exist->lower_limit && dist <= exist->upper_limit) {
                 lower_limit = exist->lower_limit;
                 upper_limit = exist->upper_limit;
                 return exist->type;
            }
       }

       return 0;
}

int BondUtil::get_bond_type(const std::string& type1, const std::string& type2)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) return exist->type;
       return 0;
}

void BondUtil::get_bond_range(const std::string& type1, const std::string& type2, double& lower_limit, double& upper_limit)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            lower_limit = exist->lower_limit;
            upper_limit = exist->upper_limit;
       }
}

int BondUtil::is_a_bond(RCSB::Atom* atom1, const std::string& type1, RCSB::Atom* atom2, const std::string& type2)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            double dist = cal_distance(atom1, atom2);
            if (dist >= exist->lower_limit && dist <= exist->upper_limit) {
                 return exist->type;
            }
       }

       return 0;
}

int BondUtil::is_a_bond_with_extension(RCSB::Atom* atom1, const std::string& type1, RCSB::Atom* atom2, const std::string& type2)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            double dist = cal_distance(atom1, atom2);
            if (dist >= exist->lower_limit && dist <= exist->upper_limit) {
                 return exist->type;
            }
            if (exist->type == 1 && dist < 2.0 && (fabs(dist - exist->lower_limit) < 0.2 || fabs(dist - exist->upper_limit) < 0.2)) {
                 return exist->type;
            }
       }

       return 0;
}

int BondUtil::is_a_bond_with_extension(const std::string& type1, const std::string& type2, const double& dist)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            if (dist >= exist->lower_limit && dist <= exist->upper_limit) {
                 return exist->type;
            }
            if (exist->type == 1 && dist < 2.0 && (fabs(dist - exist->lower_limit) < 0.2 || fabs(dist - exist->upper_limit) < 0.2)) {
                 return exist->type;
            }
       }

       return 0;
}

int BondUtil::is_a_link(const std::string& type1, const std::string& type2, const double& dist)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            if (dist >= exist->lower_limit && dist <= exist->upper_limit) {
                 return exist->type;
            }
       }

       return 0;
}

int BondUtil::is_a_link_loose(const std::string& type1, const std::string& type2, const double& dist)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            if (dist >= (exist->lower_limit * 0.75 ) && dist <= (exist->upper_limit * 1.5)) {
                 return exist->type;
            }
       }

       return 0;
}

int BondUtil::is_a_bond_loose(RCSB::Atom* atom1, RCSB::Atom* atom2)
{
       std::string label = _get_label(atom1->atom_type(), atom2->atom_type());
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            double dist = cal_distance(atom1, atom2);
            return is_a_bond_loose(atom1->atom_type(), atom2->atom_type(), dist);
       }

       return 0;
}

int BondUtil::is_a_bond_loose(const std::string& type1, const std::string& type2, const double& dist)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            double lower_limit = exist->lower_limit * 0.80;
            double upper_limit = exist->upper_limit * 1.25;
            if (exist->type == 2) {
                 lower_limit = exist->lower_limit * 0.5;
                 upper_limit = exist->upper_limit * 1.5;
            }
            if ((dist >= lower_limit) && (dist <= upper_limit)) return exist->type;
       }

       return 0;
}

int BondUtil::is_a_hydrogen_bonding(const std::string& type1, const std::string& type2, const double& dist)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            if (dist <= (exist->upper_limit * 1.5)) {
                 return exist->type;
            }
       }

       return 0;
}

int BondUtil::is_a_link_flexible(const std::string& type1, const std::string& type2, const double& dist)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            if (exist->type == 1 && dist >= (exist->lower_limit * 0.75 ) && dist <= (exist->upper_limit * 1.5)) {
                 return exist->type;
            }
       }

       return 0;
}

int BondUtil::is_a_link_with_extension(const std::string& type1, const std::string& type2, const double& dist)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (!exist) return 0;

       double extension = 0.15;
       if (exist->large_flag) extension = 0.25;
       if (dist >= (exist->lower_limit - extension) && dist <= (exist->upper_limit + extension)) return exist->type;
       if ((type1 == "H" || type1 == "D" || type2 == "H" || type2 == "D") && ( dist <= (exist->upper_limit + extension))) return exist->type;

       return 0;
}

int BondUtil::is_a_link_with_upper_limit(const std::string& type1, const std::string& type2, const double& dist)
{
       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (!exist) return 0;

       double extension = 0.15;
       if (exist->large_flag) extension = 0.25;
       if (dist <= (exist->upper_limit + extension)) return exist->type;

       return 0;
}

int BondUtil::is_glycosidic_linkage(const std::string& name1, const std::string& type1, const std::string& name2, const std::string& type2)
{
       if ((type1 == "H") || (type2 == "H") || (type1 == "D") || (type2 == "D")) return 0;
       if (((type1 == "C") && (type2 == "C")) || ((type1 != "C") && (type2 != "C"))) return 0;

       std::set<std::string> c_atom_set, o_atom_set;
       c_atom_set.clear();
       o_atom_set.clear();
       for (int i = 1; i < 10; ++i) {
            c_atom_set.insert("C" + String::IntToString(i));
            o_atom_set.insert("O" + String::IntToString(i));
       }   

       if ((type1 == "C") && (c_atom_set.find(name1) == c_atom_set.end())) return 0;
       if ((type2 == "C") && (c_atom_set.find(name2) == c_atom_set.end())) return 0;
       if ((type1 == "O") && (o_atom_set.find(name1) == o_atom_set.end())) return 0;
       if ((type2 == "O") && (o_atom_set.find(name2) == o_atom_set.end())) return 0;

       return 1;
}

int BondUtil::is_glycosidic_bond(RCSB::Atom* atom1, RCSB::Atom* atom2)
{
       std::string type1 = atom1->atom_type();
       std::string type2 = atom2->atom_type();
       if (type1.empty()) type1 = Element::getAtomSymbol(atom1->atmtype());
       if (type2.empty()) type2 = Element::getAtomSymbol(atom2->atmtype());

       int ret = is_glycosidic_linkage(atom1->atmtype(), type1, atom2->atmtype(), type2);
       if (ret == 0) return 0;

       return (is_a_bond(atom1, type1, atom2, type2));
}

int BondUtil::is_a_missing_linkage(RCSB::Atom* atom1, RCSB::Atom* atom2, const double& dist)
{
       std::string type1 = atom1->atom_type();
       std::string type2 = atom2->atom_type();
       if (type1.empty()) type1 = Element::getAtomSymbol(atom1->atmtype());
       if (type2.empty()) type2 = Element::getAtomSymbol(atom2->atmtype());

       std::string label = _get_label(type1, type2);
       _ALLOWED_BOND *exist = _find_bond(label);
       if (exist) {
            if (dist >= exist->lower_limit && dist <= exist->upper_limit) {
                 return 1; // good
            }
            if (dist > exist->upper_limit) {
                 double upper_limit = exist->upper_limit * 1.1;
                 if (dist <= upper_limit) return 1; // good
                 upper_limit = exist->upper_limit * 1.2;
                 if (dist <= upper_limit) return 2; // ok
            } else if (dist < exist->lower_limit) {
                 double lower_limit = exist->lower_limit * 0.9;
                 if (dist >= lower_limit) return 1; // good
                 lower_limit = exist->lower_limit * 0.8;
                 if (dist >= lower_limit) return 2; // ok
            }
            return (-1); // dist is too far away from allowed range
       }

       return 0; // not found defined link type
}
