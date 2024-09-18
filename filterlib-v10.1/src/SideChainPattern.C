/*
FILE:     SideChainPattern.C
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

#include <vector>

#include "SideChainPattern.h"

#define NUM_ALA_LIKE_ATOMS    6
static const char *__ala_atom_names[NUM_ALA_LIKE_ATOMS]   = { "N", "CA", "C", "O", "CB", "OXT" };

#define NUM_ALA_LIKE_RESIDUES 24
static const char *__ala_residue_names_1[NUM_ALA_LIKE_RESIDUES] = {
       "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
       "MET", "MSE", "PHE", "PRO", "PYL", "SEC", "SER", "THR", "TRP", "TYR", "UNK", "VAL"
};

#define NUM_C_GAMMA_LIKE_ATOMS_1 7
static const char *__c_gamma_atom_names_1[NUM_C_GAMMA_LIKE_ATOMS_1] = { "CG", "N", "CA", "C", "O", "CB", "OXT" };

#define NUM_C_GAMMA_LIKE_RESIDUES_1 14
static const char *__c_gamma_residue_names_1[NUM_C_GAMMA_LIKE_RESIDUES_1] = {
       "ARG", "ASN", "ASP", "GLN", "GLU", "HIS", "LEU", "LYS", "MET", "MSE", "PHE", "PYL", "TRP", "TYR"
};

#define NUM_C_GAMMA_LIKE_ATOMS_2 7
static const char *__c_gamma_atom_names_2[NUM_C_GAMMA_LIKE_ATOMS_2] = { "CG1", "N", "CA", "C", "O", "CB", "OXT" };

#define NUM_C_GAMMA_LIKE_RESIDUES_2 2
static const char *__c_gamma_residue_names_2[NUM_C_GAMMA_LIKE_RESIDUES_2] = { "ILE", "VAL" };

#define NUM_C_GAMMA_LIKE_ATOMS_3 7
static const char *__c_gamma_atom_names_3[NUM_C_GAMMA_LIKE_ATOMS_3] = { "CG2", "N", "CA", "C", "O", "CB", "OXT" };

#define NUM_C_GAMMA_LIKE_RESIDUES_3 3
static const char *__c_gamma_residue_names_3[NUM_C_GAMMA_LIKE_RESIDUES_3] = { "ILE", "THR", "VAL" };

#define NUM_C_DELTA_LIKE_ATOMS_1 8
static const char *__c_delta_atom_names_1[NUM_C_DELTA_LIKE_ATOMS_1] = { "CD", "N", "CA", "C", "O", "CB", "CG", "OXT" };

#define NUM_C_DELTA_LIKE_RESIDUES_1 5
static const char *__c_delta_residue_names_1[NUM_C_DELTA_LIKE_RESIDUES_1] = { "ARG", "GLN", "GLU", "LYS", "PYL" };

#define NUM_C_DELTA_LIKE_ATOMS_2 9
static const char *__c_delta_atom_names_2[NUM_C_DELTA_LIKE_ATOMS_2] = { "CD1", "N", "CA", "C", "O", "CB", "CG", "CG1", "OXT" };

#define NUM_C_DELTA_LIKE_RESIDUES_2 4
static const char *__c_delta_residue_names_2[NUM_C_DELTA_LIKE_RESIDUES_2] = { "ILE", "PHE", "TRP", "TYR" }; 

#define NUM_C_DELTA_LIKE_ATOMS_3 8
static const char *__c_delta_atom_names_3[NUM_C_DELTA_LIKE_ATOMS_3] = { "CD2", "N", "CA", "C", "O", "CB", "CG", "OXT" };

#define NUM_C_DELTA_LIKE_RESIDUES_3 4
static const char *__c_delta_residue_names_3[NUM_C_DELTA_LIKE_RESIDUES_3] = { "HIS", "PHE", "TRP", "TYR" }; 

#define NUM_ILE_LIKE_ATOMS  8
static const char *__ile_atom_names[NUM_ILE_LIKE_ATOMS] = { "N", "CA", "C", "O", "CB", "CG1", "CG2", "OXT" };

#define NUM_ILE_LIKE_RESIDUES 2
static const char *__ile_residue_names[NUM_ILE_LIKE_RESIDUES] = { "ILE", "VAL" };

#define NUM_SER_LIKE_ATOMS  8
static const char *__ser_atom_names[NUM_SER_LIKE_ATOMS] = { "N", "CA", "C", "O", "CB", "OG", "OG1", "OXT" };

#define NUM_SER_LIKE_RESIDUES 2
static const char *__ser_residue_names[NUM_SER_LIKE_RESIDUES] = { "SER", "THR" };

#define NUM_PHE_LIKE_ATOMS 12
static const char *__phe_atom_names[NUM_PHE_LIKE_ATOMS] = { "N", "CA", "C", "O", "CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OXT" };

#define NUM_PHE_LIKE_RESIDUES 2
static const char *__phe_residue_names[NUM_PHE_LIKE_RESIDUES] = { "PHE", "TYR" };

#define NUM_LYS_LIKE_ATOMS 10
static const char *__lys_atom_names[NUM_LYS_LIKE_ATOMS] = { "N", "CA", "C", "O", "CB", "CG", "CD", "CE", "NZ", "OXT" };

#define NUM_LYS_LIKE_RESIDUES 2
static const char *__lys_residue_names[NUM_LYS_LIKE_RESIDUES] = { "LYS", "PYL" };

// first: pattern name
// second.first: allowed residue set
// second.second: allowed atom set
static std::vector<std::pair<std::string, std::vector<std::pair<std::set<std::string>, std::set<std::string> > > > > __side_chain_pattern_list;

// key: pattern name
// value: allowed residue set
static std::map<std::string, std::set<std::string> > __pattern_residue_set_mapping;

// key: pattern name
// value: first - key atom name
// value: second - allowed residue set
static std::map<std::string, std::vector<std::pair<std::string, std::set<std::string> > > > __pattern_key_atom_mapping;

void SideChainPattern::initialize()
{
       __side_chain_pattern_list.clear();
       __pattern_residue_set_mapping.clear();
       __pattern_key_atom_mapping.clear();

       std::set<std::string> atom_set, residue_set;
       std::vector<std::pair<std::set<std::string>, std::set<std::string> > > vec;
       std::vector<std::pair<std::string, std::set<std::string> > > vec1;

       // ala-gly-like
       vec.clear();

       atom_set.clear();
       for (int i = 0; i < NUM_ALA_LIKE_ATOMS; ++i) atom_set.insert(__ala_atom_names[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_ALA_LIKE_RESIDUES; ++i) residue_set.insert(__ala_residue_names_1[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));

       __side_chain_pattern_list.push_back(std::make_pair("ala-gly-like", vec));

       // c-gamma-like
       vec.clear();
       vec1.clear();

       atom_set.clear();
       for (int i = 0; i < NUM_C_GAMMA_LIKE_ATOMS_1; ++i) atom_set.insert(__c_gamma_atom_names_1[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_C_GAMMA_LIKE_RESIDUES_1; ++i) residue_set.insert(__c_gamma_residue_names_1[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));
       vec1.push_back(std::make_pair(__c_gamma_atom_names_1[0], residue_set));

       atom_set.clear();
       for (int i = 0; i < NUM_C_GAMMA_LIKE_ATOMS_2; ++i) atom_set.insert(__c_gamma_atom_names_2[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_C_GAMMA_LIKE_RESIDUES_2; ++i) residue_set.insert(__c_gamma_residue_names_2[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));
       vec1.push_back(std::make_pair(__c_gamma_atom_names_2[0], residue_set));

       atom_set.clear();
       for (int i = 0; i < NUM_C_GAMMA_LIKE_ATOMS_3; ++i) atom_set.insert(__c_gamma_atom_names_3[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_C_GAMMA_LIKE_RESIDUES_3; ++i) residue_set.insert(__c_gamma_residue_names_3[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));
       vec1.push_back(std::make_pair(__c_gamma_atom_names_3[0], residue_set));

       __side_chain_pattern_list.push_back(std::make_pair("c-gamma-like", vec));
       __pattern_key_atom_mapping.insert(std::make_pair("c-gamma-like", vec1));

       // c-delta-like
       vec.clear();
       vec1.clear();

       atom_set.clear();
       for (int i = 0; i < NUM_C_DELTA_LIKE_ATOMS_1; ++i) atom_set.insert(__c_delta_atom_names_1[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_C_DELTA_LIKE_RESIDUES_1; ++i) residue_set.insert(__c_delta_residue_names_1[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));
       vec1.push_back(std::make_pair(__c_delta_atom_names_1[0], residue_set));

       atom_set.clear();
       for (int i = 0; i < NUM_C_DELTA_LIKE_ATOMS_2; ++i) atom_set.insert(__c_delta_atom_names_2[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_C_DELTA_LIKE_RESIDUES_2; ++i) residue_set.insert(__c_delta_residue_names_2[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));
       vec1.push_back(std::make_pair(__c_delta_atom_names_2[0], residue_set));

       atom_set.clear();
       for (int i = 0; i < NUM_C_DELTA_LIKE_ATOMS_3; ++i) atom_set.insert(__c_delta_atom_names_3[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_C_DELTA_LIKE_RESIDUES_3; ++i) residue_set.insert(__c_delta_residue_names_3[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));
       vec1.push_back(std::make_pair(__c_delta_atom_names_3[0], residue_set));

       __side_chain_pattern_list.push_back(std::make_pair("c-delta-like", vec));
       __pattern_key_atom_mapping.insert(std::make_pair("c-delta-like", vec1));

       // ile-val-like
       vec.clear();

       atom_set.clear();
       for (int i = 0; i < NUM_ILE_LIKE_ATOMS; ++i) atom_set.insert(__ile_atom_names[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_ILE_LIKE_RESIDUES; ++i) residue_set.insert(__ile_residue_names[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));

       __side_chain_pattern_list.push_back(std::make_pair("ile-val-like", vec));

       // ser-thr-like
       vec.clear();

       atom_set.clear();
       for (int i = 0; i < NUM_SER_LIKE_ATOMS; ++i) atom_set.insert(__ser_atom_names[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_SER_LIKE_RESIDUES; ++i) residue_set.insert(__ser_residue_names[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));

       __side_chain_pattern_list.push_back(std::make_pair("ser-thr-like", vec));

       // phe-tyr-like
       vec.clear();

       atom_set.clear();
       for (int i = 0; i < NUM_PHE_LIKE_ATOMS; ++i) atom_set.insert(__phe_atom_names[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_PHE_LIKE_RESIDUES; ++i) residue_set.insert(__phe_residue_names[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));

       __side_chain_pattern_list.push_back(std::make_pair("phe-tyr-like", vec));

       // lys-pyl-like
       vec.clear();

       atom_set.clear();
       for (int i = 0; i < NUM_LYS_LIKE_ATOMS; ++i) atom_set.insert(__lys_atom_names[i]);
       residue_set.clear();
       for (int i = 0; i < NUM_LYS_LIKE_RESIDUES; ++i) residue_set.insert(__lys_residue_names[i]);
       vec.push_back(std::make_pair(residue_set, atom_set));

       __side_chain_pattern_list.push_back(std::make_pair("lys-pyl-like", vec));

       __pattern_residue_set_mapping.clear();
       for (std::vector<std::pair<std::string, std::vector<std::pair<std::set<std::string>, std::set<std::string> > > > >::const_iterator
            pos = __side_chain_pattern_list.begin(); pos != __side_chain_pattern_list.end(); ++pos) {
            residue_set = pos->second[0].first;
            for (unsigned int i = 1; i < pos->second.size(); ++i) {
                 for (std::set<std::string>::const_iterator spos = pos->second[i].first.begin(); spos != pos->second[i].first.end(); ++spos) {
                      residue_set.insert(*spos);
                 }
            }
            __pattern_residue_set_mapping.insert(make_pair(pos->first, residue_set));
       }
}

std::string SideChainPattern::get_pattern(const std::string& resname, const std::set<std::string>& atom_name_set)
{
       for (std::vector<std::pair<std::string, std::vector<std::pair<std::set<std::string>, std::set<std::string> > > > >::const_iterator
            pos = __side_chain_pattern_list.begin(); pos != __side_chain_pattern_list.end(); ++pos) {
            for (std::vector<std::pair<std::set<std::string>, std::set<std::string> > >::const_iterator
                 vpos = pos->second.begin(); vpos != pos->second.end(); ++vpos) {
                 if (vpos->first.find(resname) == vpos->first.end()) continue;

                 bool found_extra_atom = false;
                 for (std::set<std::string>::const_iterator spos = atom_name_set.begin(); spos != atom_name_set.end(); ++spos) {
                      if (vpos->second.find(*spos) == vpos->second.end()) {
                           found_extra_atom = true;
                           break;
                      }
                 }
                 if (found_extra_atom) continue;

                 return pos->first;
            }
       }
       return ""; 
}

bool SideChainPattern::is_a_valid_rename(const std::string& pattern, const std::string& old_resname, const std::string& new_resname)
{
       std::map<std::string, std::set<std::string> >::const_iterator mpos = __pattern_residue_set_mapping.find(pattern);
       if (mpos == __pattern_residue_set_mapping.end()) return false;
       if ((mpos->second.find(old_resname) == mpos->second.end()) || (mpos->second.find(new_resname) == mpos->second.end())) return false;
       return true;
}

void SideChainPattern::get_key_atoms(const std::string& pattern, std::set<std::string>& atom_set)
{
       atom_set.clear();

       std::map<std::string, std::vector<std::pair<std::string, std::set<std::string> > > >::const_iterator mpos = __pattern_key_atom_mapping.find(pattern);
       if (mpos == __pattern_key_atom_mapping.end()) return;

       for (std::vector<std::pair<std::string, std::set<std::string> > >::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
            atom_set.insert(vpos->first);
       }
}

bool SideChainPattern::get_change_info(const std::string& key_atom, const std::string& pattern, const std::string& old_resname, const std::string& new_resname,
                                       std::set<std::string>& delete_atom_set, std::map<std::string, std::pair<std::string, std::string> >& rename_atom_map)
{
       delete_atom_set.clear();
       rename_atom_map.clear();

       if (new_resname == "UNK") return true;

       if (!pattern.empty()) {
            if (pattern == "ala-gly-like") {
                 if (new_resname == "GLY") {
                      delete_atom_set.insert("CB");
                      delete_atom_set.insert("HB1");
                      delete_atom_set.insert("HB2");
                      delete_atom_set.insert("HB3");
                      delete_atom_set.insert("DB1");
                      delete_atom_set.insert("DB2");
                      delete_atom_set.insert("DB3");
                 }
                 return true;
            } else if (pattern == "ser-thr-like") {
                 if (old_resname == "SER") rename_atom_map.insert(std::make_pair("OG", std::make_pair("OG1", "")));
                 else if (new_resname == "SER") rename_atom_map.insert(std::make_pair("OG1", std::make_pair("OG", "")));
                 return true;
            }
                 
            std::map<std::string, std::vector<std::pair<std::string, std::set<std::string> > > >::const_iterator
                mpos = __pattern_key_atom_mapping.find(pattern);
            if (mpos == __pattern_key_atom_mapping.end()) return true;

            int old_index = -1, new_index = -1;
            for (unsigned int i = 0; i < mpos->second.size(); ++i) {
                 if ((key_atom == mpos->second[i].first) && (mpos->second[i].second.find(old_resname) != mpos->second[i].second.end())) old_index = i;
                 if (mpos->second[i].second.find(new_resname) != mpos->second[i].second.end()) new_index = i;
            }
            if ((old_index < 0) || (new_index < 0)) return false;

            if (key_atom != mpos->second[new_index].first) rename_atom_map.insert(std::make_pair(key_atom, std::make_pair(mpos->second[new_index].first, "")));
            if (pattern == "c-delta-like") {
                 if (old_resname == "ILE") rename_atom_map.insert(std::make_pair("CG1", std::make_pair("CG", "")));
                 else if (new_resname == "ILE") rename_atom_map.insert(std::make_pair("CG", std::make_pair("CG1", "")));
            }
       }

       if (old_resname == "ASP" && new_resname == "ASN")      rename_atom_map.insert(std::make_pair("OD2", std::make_pair("ND2", "N")));
       else if (old_resname == "ASN" && new_resname == "ASP") rename_atom_map.insert(std::make_pair("ND2", std::make_pair("OD2", "O")));
       else if (old_resname == "GLU" && new_resname == "GLN") rename_atom_map.insert(std::make_pair("OE2", std::make_pair("NE2", "N")));
       else if (old_resname == "GLN" && new_resname == "GLU") rename_atom_map.insert(std::make_pair("NE2", std::make_pair("OE2", "O")));

       return true;
}
