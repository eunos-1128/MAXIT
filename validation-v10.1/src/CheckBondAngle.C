/*
FILE:     CheckBondAngle.C
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
#include <stdlib.h>
#include <math.h>

#include "CompositeIndex.h"
#include "GeometryUtil.h"
#include "GetPairList.h"
#include "GetPairList.C"
#include "ValidateObj.h"

#define AA_BOND_CUTOFF 1.65
#define NA_BOND_CUTOFF 1.85

const std::vector<_MEAN_STD>& ValidateObj::_getInterResidueBondStandard(RCSB::Residue* curr, const std::string& polymer_type, double& cutoff_value)
{
       std::string bond = "";

       cutoff_value = AA_BOND_CUTOFF;

       if (polymer_type == "ATOMP") {
            bond = "AA_LINK";
            if (curr->ResName() == "PRO" ||
                curr->ResName() == "DPR") {
                 bond = "PRO_LINK";
            }
       } else if (polymer_type == "ATOMN") {
            cutoff_value = NA_BOND_CUTOFF;
            bond = "NN_LINK";
       }

       return StandardUtil::GetBond(bond);
}

void ValidateObj::_checkInterResidueBond(const int& mol_id, RCSB::Residue* prev, RCSB::Residue* curr, const std::string& polymer_type, std::set<std::string>&
                                         found_allowed_bonds, const bool& check_linkage_only)
{
       found_allowed_bonds.clear();

       double cutoff_value;
       const std::vector<_MEAN_STD>& bondValues = _getInterResidueBondStandard(curr, polymer_type, cutoff_value);
       if (bondValues.empty()) return;

       std::vector<std::pair<RCSB::Residue*, std::string> > res_atom_pair_list;
       std::vector<std::vector<RCSB::Atom*> > pair_lists;
       for (std::vector<_MEAN_STD>::const_iterator pos = bondValues.begin(); pos != bondValues.end(); ++pos) {
            res_atom_pair_list.clear();
            res_atom_pair_list.push_back(std::make_pair(prev, pos->ats[0]));
            res_atom_pair_list.push_back(std::make_pair(curr, pos->ats[1]));
            _get_atom_pair_list(res_atom_pair_list, pair_lists);
            if (pair_lists.empty()) continue;

            for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = pair_lists.begin(); apos != pair_lists.end(); ++apos) {
                 if (_checkBondOutlier(mol_id, *pos, *apos, "", "Y", cutoff_value, check_linkage_only)) {
                      std::string index = CompositeIndex::getIndex((*apos)[0]->atmtype(), (*apos)[0]->alt_loc(), (*apos)[1]->atmtype(), (*apos)[1]->alt_loc());
                      found_allowed_bonds.insert(index);
                 }
            }
       }
}

const std::vector<_MEAN_STD>& ValidateObj::_getInterResidueAngleStandard(const int& type, RCSB::Residue* curr, const std::string& polymer_type)
{
       std::string angle = "";

       if (type == ANGLE_TYPE_1) {
            if (polymer_type == "ATOMP") {
                 angle = "AA_ANG1";
                 if (curr->ResName() == "GLY")
                      angle = "GLY_ANG1";
                 else if (curr->ResName() == "PRO" ||
                          curr->ResName() == "DPR")
                      angle = "PRO_ANG1";
            } else if (polymer_type == "ATOMN")
                 angle = "NN_ANG1";
       } else if (type == ANGLE_TYPE_2) {
            if (polymer_type == "ATOMP") {
                 angle = "AA_ANG2";
                 if (curr->ResName() == "GLY")
                      angle = "GLY_ANG2";
                 else if (curr->ResName() == "PRO" ||
                          curr->ResName() == "DPR")
                      angle = "PRO_ANG2";
            } else if (polymer_type == "ATOMN")
                 angle = "NN_ANG2";
       }
       return StandardUtil::GetAngle(angle);
}

void ValidateObj::_checkInterResidueAngle(const int& type, const int& mol_id, RCSB::Residue* prev, RCSB::Residue* curr, const std::string& polymer_type,
                                          const std::set<std::string>& found_allowed_bonds, const std::map<std::string, std::string>&
                                          conformer_cis_trans_mapping)
{

       const std::vector<_MEAN_STD>& angleValues = _getInterResidueAngleStandard(type, curr, polymer_type);
       if (angleValues.empty()) return;

       std::vector<std::pair<RCSB::Residue*, std::string> > res_atom_pair_list;
       std::vector<std::vector<RCSB::Atom*> > pair_lists;
       for (std::vector<_MEAN_STD>::const_iterator pos = angleValues.begin(); pos != angleValues.end(); ++pos) {
            res_atom_pair_list.clear();
            res_atom_pair_list.push_back(std::make_pair(prev, pos->ats[0]));
            if (type == ANGLE_TYPE_1)
                 res_atom_pair_list.push_back(std::make_pair(prev, pos->ats[1]));
            else if (type == ANGLE_TYPE_2)
                 res_atom_pair_list.push_back(std::make_pair(curr, pos->ats[1]));
            else continue;
            res_atom_pair_list.push_back(std::make_pair(curr, pos->ats[2]));
            _get_atom_pair_list(res_atom_pair_list, pair_lists);
            if (pair_lists.empty()) continue;

            for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = pair_lists.begin(); apos != pair_lists.end(); ++apos) {
                 std::string index = "";
                 if (type == ANGLE_TYPE_1) {
                      index = CompositeIndex::getIndex((*apos)[1]->atmtype(), (*apos)[1]->alt_loc(), (*apos)[2]->atmtype(), (*apos)[2]->alt_loc());
                 } else if (type == ANGLE_TYPE_2) {
                      index = CompositeIndex::getIndex((*apos)[0]->atmtype(), (*apos)[0]->alt_loc(), (*apos)[1]->atmtype(), (*apos)[1]->alt_loc());
                 }
                 if (found_allowed_bonds.find(index) == found_allowed_bonds.end()) continue;

                 std::string cis_trans_type = "";
                 std::string alt_loc = GetAltLoc(*apos);
                 std::map<std::string, std::string>::const_iterator mpos = conformer_cis_trans_mapping.find(alt_loc);
                 if (mpos != conformer_cis_trans_mapping.end()) cis_trans_type = mpos->second;

                 _checkAngleOutlier(mol_id, *pos, *apos, cis_trans_type, "Y");
            }
       }
}

void ValidateObj::_checkIntraResidueValues(const int& mol_id, RCSB::Residue* residue, const std::map<std::string, std::string>& conformer_cis_trans_mapping)
{
       const std::vector<_MEAN_STD>& bondValues = StandardUtil::GetBond(residue->ResName());
       const std::vector<_MEAN_STD>& angleValues = StandardUtil::GetAngle(residue->ResName());
       if (bondValues.empty() && angleValues.empty()) return;

       _checkIntraResidueValues(BOND_TYPE, mol_id, bondValues, residue, conformer_cis_trans_mapping);
       _checkIntraResidueValues(ANGLE_TYPE, mol_id, angleValues, residue, conformer_cis_trans_mapping);
}

void ValidateObj::_checkIntraResidueValues(const int& type, const int& mol_id, const std::vector<_MEAN_STD>& standardValues, RCSB::Residue* residue,
                                           const std::map<std::string, std::string>& conformer_cis_trans_mapping)
{
       if (standardValues.empty()) return;

       std::vector<std::pair<RCSB::Residue*, std::string> > res_atom_pair_list;
       std::vector<std::vector<RCSB::Atom*> > pair_lists;
       for (std::vector<_MEAN_STD>::const_iterator pos = standardValues.begin(); pos != standardValues.end(); ++pos) {
            res_atom_pair_list.clear();
            for (std::vector<std::string>::const_iterator apos = pos->ats.begin(); apos != pos->ats.end(); ++apos) {
                 res_atom_pair_list.push_back(std::make_pair(residue, *apos));
            }
            _get_atom_pair_list(res_atom_pair_list, pair_lists);
            if (pair_lists.empty()) continue;

            for (std::vector<std::vector<RCSB::Atom*> >::const_iterator apos = pair_lists.begin(); apos != pair_lists.end(); ++apos) {
                 std::string cis_trans_type = "";
                 std::string alt_loc = GetAltLoc(*apos);
                 std::map<std::string, std::string>::const_iterator mpos = conformer_cis_trans_mapping.find(alt_loc);
                 if (mpos != conformer_cis_trans_mapping.end()) cis_trans_type = mpos->second;

                 if (type == BOND_TYPE)
                      _checkBondOutlier(mol_id, *pos, *apos, cis_trans_type, "N");
                 else if (type == ANGLE_TYPE)
                      _checkAngleOutlier(mol_id, *pos, *apos, cis_trans_type, "N");
            }
       }
}

bool ValidateObj::_checkBondOutlier(const int& mol_id, const _MEAN_STD& standard, const std::vector<RCSB::Atom*>& atoms, const std::string& cis_trans_type,
                                    const std::string& linker_flag, const double& cutoff_distance, const bool& check_linkage_only)
{
       double val = cal_distance(atoms[0], atoms[1]);
       if (cutoff_distance > 0 && val > cutoff_distance) return false;
       if (check_linkage_only) return true;

       double ept = 0;
       _checkOutlier(mol_id, standard, cis_trans_type, linker_flag, val, ept, atoms, _violated_bond);

       _bond_count++;
       double d = val - ept;
       _rmsbond += d * d;

       return true;
}

void ValidateObj::_checkAngleOutlier(const int& mol_id, const _MEAN_STD& standard, const std::vector<RCSB::Atom*>& atoms, const std::string&
                                     cis_trans_type, const std::string& linker_flag)
{
       double val = cal_angle(atoms[0], atoms[1], atoms[2]);

       double ept = 0;
       _checkOutlier(mol_id, standard, cis_trans_type, linker_flag, val, ept, atoms, _violated_angle);

       _angle_count++;
       double d = val - ept;
       _rmsangle += d * d;
}

void ValidateObj::_checkOutlier(const int& mol_id, const _MEAN_STD& standard, const std::string& cis_trans_type, const std::string& linker_flag,
                                const double& val, double& ept, const std::vector<RCSB::Atom*>& atoms, std::list<Value>& _outlier_list)
{
       double std = 0;
       if (StandardUtil::isOutlier(val, cis_trans_type, ept, std, standard.mean_std)) {
            Value value;
            value.set_Mol_ID(mol_id);
            value.set_val(val);
            value.set_ept(ept);
            value.set_std(std);
            value.set_sval(linker_flag);
            for (std::vector<RCSB::Atom*>::const_iterator pos = atoms.begin(); pos != atoms.end(); ++pos) {
                 value.insert_atom(*pos);
            }
            _outlier_list.push_back(value);
       }
}
