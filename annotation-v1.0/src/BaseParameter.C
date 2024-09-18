/*
FILE:     BaseParameter.C
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
#include <math.h>

#include "BaseParameter.h"
#include "BaseParameter_global.h"
#include "SuperposeUtil.h"

BaseParameter::BaseParameter()
{
       _clear();
       _init();
}

BaseParameter::~BaseParameter()
{
       _clear();
}

void BaseParameter::_clear()
{
       _cell = NULL;
       _ccDic = NULL;
       _ref_na_type_index.clear();
       _res_mapping.clear();
       _code_mapping.clear();
       _ring_atoms.clear();
       _res_linked_atom_mapping.clear();
}

void BaseParameter::_init()
{
       for (int i = 0; i < NUM_REF_NA; ++i) {
            _ref_na_type_index.insert(std::make_pair(ref_na_type[i].na_type, i));
       }

       for (int i = 0; i < NUM_RES_MAPPING; ++i) {
            _res_mapping.insert(std::make_pair(_res_code_mapping[i][0], _res_code_mapping[i][1]));
       }

       for (int i = 0; i < NUM_CODE_MAPPING; ++i) {
            _code_mapping.insert(std::make_pair(_code_code_mapping[i][0], _code_code_mapping[i][1]));
       }

       for (int i = STAND_R.num_atom - 1; i >= 0; --i) {
            _ring_atoms.insert(std::make_pair(STAND_R.atomname[i], _ring_atoms.size()));
       }
}

void BaseParameter::setCell(CrySymmetry* cell)
{
       _cell = cell;
}

void BaseParameter::setCCDic(ConnectDic* ccdic)
{
       _ccDic = ccdic;
}

void BaseParameter::getBaseFrame(_BASE& base, RCSB::Residue* res, const int& op, const int& lx, const int& ly, const int& lz)
{
       base.res = res;
       base.code.clear();
       base.has_base = false;
       base.op = op;
       base.lx = lx;
       base.ly = ly;
       base.lz = lz;
       base.coord_xyz.clear();
       base.polygons.clear();

       std::string RY_code;
       double eRing_xyz[27], sRing_xyz[27];
       int nmatch = _find_match(res, eRing_xyz, sRing_xyz, base.code, RY_code, base.polygons);

       RCSB::Atom* atom = res->find_atom("P");
       if (atom) base.coord_xyz.insert(std::make_pair("P", atom->orig()));

       atom = res->find_atom("O3'");
       if (atom) base.coord_xyz.insert(std::make_pair("O", atom->orig()));

       atom = res->find_atom("C1'");
       if (atom) base.coord_xyz.insert(std::make_pair("C", atom->orig()));

       if (RY_code == "R")
            atom = res->find_atom("N9");
       else atom = res->find_atom("N1");
       if (atom) base.coord_xyz.insert(std::make_pair("N", atom->orig()));

       if (op || lx || ly || lz) {
            for (std::map<std::string, COORD>::iterator mpos = base.coord_xyz.begin(); mpos != base.coord_xyz.end(); ++mpos) {
                 _cell->symmetry_operation(mpos->second, op, lx, ly, lz);
            }
            for (std::vector<COORD>::iterator pos = base.polygons.begin(); pos != base.polygons.end(); ++pos) {
                 _cell->symmetry_operation(*pos, op, lx, ly, lz);
            }
       }

       if (nmatch < 3) return; 

       base.has_base = true;

       if (op || lx || ly || lz) {
            for (int k = 0; k < nmatch; ++k) {
                 _cell->symmetry_operation(eRing_xyz[3 * k], eRing_xyz[3 * k + 1], eRing_xyz[3 * k + 2], op, lx, ly, lz);
            }
       }

       double rt[12], ave_exyz[4], ave_sxyz[4];
       getSuperposeMatrix(nmatch, eRing_xyz, sRing_xyz, rt, ave_exyz, ave_sxyz);

       for (int k = 0; k < 3; ++k) {
            for (int l = 0; l < 3; ++l)
                 ave_exyz[k] -= rt[k * 3 + l] * ave_sxyz[l];
       }

       base.org.x = ave_exyz[0];
       base.org.y = ave_exyz[1];
       base.org.z = ave_exyz[2];
       for (int k = 0; k < 3; ++k) {
            base.orien[k].x = rt[k];
            base.orien[k].y = rt[3 + k];
            base.orien[k].z = rt[6 + k];
       }
}

int BaseParameter::_find_match(RCSB::Residue* residue, double* eRing_xyz, double* sRing_xyz, std::string& code, std::string& RY_code,
                               std::vector<COORD>& polygons)
{
       RY_code.clear();
       code = _getNACode(residue->ResName());
       if (code.empty()) return 0;

       std::map<std::string, std::string>::const_iterator pos = _code_mapping.find(code);
       if (pos != _code_mapping.end()) RY_code = pos->second;

       if (RY_code.empty()) return 0;

       REF_NA *block = NULL;
       std::map<std::string, int>::const_iterator mpos = _ref_na_type_index.find(RY_code);
       if (mpos != _ref_na_type_index.end()) block = ref_na_type[mpos->second].ref_type;
       if (block == NULL) return 0;

       REF_NA *ref = NULL;
       mpos = _ref_na_type_index.find(code);
       if (mpos != _ref_na_type_index.end()) ref = ref_na_type[mpos->second].ref_type;
       if (ref == NULL) return 0;

       std::map<std::string, int> index;
       index.clear();
       for (int i = 0; i < ref->num_atom; ++i) {
            index.insert(std::make_pair(ref->atomname[i], i));
       }

       std::map<unsigned int, COORD> polygon_coords;
       polygon_coords.clear();

       int nmatch = 0;
       for (int k = 0; k < block->num_atom; ++k) {
            mpos = index.find(block->atomname[k]);
            if (mpos == index.end()) continue;

            RCSB::Atom* atom = residue->find_atom(block->atomname[k]);
            if (atom == NULL) continue;

            std::map<std::string, unsigned int>::const_iterator mmpos = _ring_atoms.find(block->atomname[k]);
            if (mmpos != _ring_atoms.end()) {
                 std::string atom_name = _getPolygonAtomName(residue->ResName(), block->atomname[k]);
                 if (atom_name == mmpos->first) {
                      polygon_coords.insert(std::make_pair(mmpos->second, atom->orig())); 
                 } else {
                      RCSB::Atom* atom1 = residue->find_atom(atom_name);
                      if (atom1) polygon_coords.insert(std::make_pair(mmpos->second, atom1->orig()));
                      else polygon_coords.insert(std::make_pair(mmpos->second, atom->orig()));
                 }
            }

            eRing_xyz[nmatch * 3] = atom->orig().x;
            eRing_xyz[nmatch*3+1] = atom->orig().y;
            eRing_xyz[nmatch*3+2] = atom->orig().z;
            for (int l = 0; l < 3; ++l) {
                 sRing_xyz[nmatch * 3 + l] = ref->ref[mpos->second][l];
            }
            nmatch++;
       }

       for (std::map<unsigned int, COORD>::const_iterator mpos = polygon_coords.begin(); mpos != polygon_coords.end(); ++mpos) {
            polygons.push_back(mpos->second);
       }

       return nmatch;
}

std::string BaseParameter::_getNACode(const std::string& resname)
{
       std::map<std::string, std::string>::const_iterator mpos = _res_mapping.find(resname);
       if (mpos != _res_mapping.end()) return mpos->second;

       try {
            const ConnectFormat& drug = _ccDic->find_drug(resname);

            std::string code = drug.getMetaData("one_letter_code");
            std::string parent = drug.getMetaData("mon_nstd_parent_comp_id");
            if (!code.empty()) {
                 if (_ref_na_type_index.find(code) != _ref_na_type_index.end())
                      return code;
            } else if (!parent.empty()) {
                 const ConnectFormat& drug1 = _ccDic->find_drug(parent);
                 code = drug1.getMetaData("one_letter_code");
                 if (!code.empty() && _ref_na_type_index.find(code) !=
                     _ref_na_type_index.end()) 
                      return code;
            }
       } catch (const std::exception& exc) {}
       return "";
}

std::string BaseParameter::_getPolygonAtomName(const std::string& resname, const std::string& atom_name)
{
       std::map<std::string, std::map<std::string, std::vector<std::string> > >::const_iterator rpos = _res_linked_atom_mapping.find(resname);
       if (rpos == _res_linked_atom_mapping.end()) {
            std::map<std::string, std::vector<std::string> > linked_atom_mapping;
            linked_atom_mapping.clear();
            try {
                 const ConnectFormat& drug = _ccDic->find_drug(resname);
                 linked_atom_mapping = drug.getLinkedAtoms();
            } catch (const std::exception& exc) {}
            _res_linked_atom_mapping.insert(std::make_pair(resname, linked_atom_mapping));
            rpos = _res_linked_atom_mapping.find(resname);
       }

       std::map<std::string, std::vector<std::string> >::const_iterator mpos = rpos->second.find(atom_name);
       if (mpos == rpos->second.end()) return atom_name;

       for (std::vector<std::string>::const_iterator pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
            if ((*pos)[0] == 'H') continue;
            if (_ring_atoms.find(*pos) != _ring_atoms.end()) continue;
            return *pos;
       }
       return atom_name;
}
