/*
FILE:     ValidateObj.h
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
#ifndef _H_VALIDATEOBJ_H_
#define _H_VALIDATEOBJ_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "FileObj.h"
#include "Link.h"
#include "StandardUtil.h"
#include "PlanarityUtil.h"
#include "Value.h"

#define BOND_TYPE    1
#define ANGLE_TYPE   2
#define ANGLE_TYPE_1 3
#define ANGLE_TYPE_2 4 

#define WATER_CUTOFF 5.805  // 5.8 + 0.005 (add 0.005 to round up to .xx)

class ValidateObj: public FileObj {
   public:

       /**
       **  Constructs ValidateObj
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       ValidateObj();

       /**
       **  Destructs ValidateObj
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       ~ValidateObj();

       void calculate_all_contact();
       void CheckGeometry(const std::string& sigma = "");
       void CheckChiralityAndPlanarity();
       void WriteOutput();

       /**
       **  Convert PDB format LINK records into _LINK and insert to _links.
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void InputPDBLINKRecord();

       /**
       **  Convert PDB format SSBOND records into _SSBOND and insert to _ssbonds.
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void InputPDBSSBONDRecord();

       /**
       **  Read struct_conn & pdbx_struct_link categories and input
       **  values into _links, _ssbonds, _sltbrgs & _bspairs 
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void Read_StructLink_and_StructConn();

       /**
       **  Check out of range waters
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void CheckDistantWaters();

       /**
       **  Output out of range waters to cif category
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void OutputDistantWater();

       /**
       **  Calculate covalent bonds
       **
       **  \param[in]: link_radii - 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void CalculateCovalentBonds(const std::string& link_radii);

       void checkChainContact();

       void SwitchLabeling();

       void CheckLabeling();

   private:
       // key: pdb_chnid1_pdb_resnam1_pdb_resnum1_ins_code1_pdb_atmnam1_alt_loc1_pdb_chnid2_pdb_resnam2_pdb_resnum2_ins_code2_pdb_atmnam2_alt_loc2
       std::map<std::string, std::set<double> > _defined_links;

   protected:
       bool _exclude_extra_hydrogen_flag;
       double _rmsbond;
       unsigned int _bond_count;
       double _rmsangle;
       unsigned int _angle_count;
       unsigned int _cis_trans_count;
       unsigned int _main_planarity_count;
       unsigned int _phi_psi_count;
       unsigned int _side_planarity_count;
       bool _cal_a_contact_flag;
       bool _cal_s_contact_flag;
       bool _cal_bond_flag;
       bool _cal_angle_flag;
       bool _cal_cis_trans_flag;
       bool _cal_main_planarity_flag;
       bool _cal_phi_psi_flag;
       bool _cal_polymer_linkage_flag;
       bool _cal_side_planarity_flag;
       bool _cal_chiral_center_flag;
       bool _cal_out_range_water_flag;
       bool _cal_special_position_flag;
       bool _distance_checking_upper_limit_only_flag;
       std::list<Value> _a_contact;
       std::list<Value> _s_contact;
       std::list<Value> _violated_bond;
       std::list<Value> _violated_angle;
       std::list<Value> _violated_cis_trans;
       std::list<Value> _violated_main_planarity;
       std::list<Value> _violated_phi_psi;
       std::list<Value> _violated_polymer_linkage;
       std::list<Value> _violated_side_planarity;
       std::list<Value> _violated_chiral_center;
       // vector[0]: PDB_model_num
       // vector[1]: Distance
       // vector[2]: Neighbor atom type: polymer/non-polymer
       std::list<std::pair<RCSB::Atom*, std::vector<std::string> > > _out_range_waters;
       std::list<Value> _special_position_atoms;
       // vector[0]: PDB_model_num
       // vector[1]: polymer_flag
       // vector[2]: occupancy_flag
       // vector[3]: auth_asym_id
       // vector[4]: auth_comp_id
       // vector[5]: auth_seq_id
       // vector[6]: PDB_ins_code
       // vector[7]: label_asym_id
       // vector[8]: label_comp_id
       // vector[9]: label_seq_id
       std::list<std::vector<std::string> > _unobs_or_zero_occ_residues;
       // vector[0]:  PDB_model_num
       // vector[1]:  polymer_flag
       // vector[2]:  occupancy_flag
       // vector[3]:  auth_asym_id
       // vector[4]:  auth_comp_id
       // vector[5]:  auth_seq_id
       // vector[6]:  PDB_ins_code
       // vector[7]:  auth_atom_id
       // vector[8]:  label_alt_id
       // vector[9]:  label_asym_id
       // vector[10]: label_comp_id
       // vector[11]: label_seq_id
       // vector[12]: label_atom_id
       std::list<std::vector<std::string> > _unobs_or_zero_occ_atoms;
       // first[0]:  PDB_model_num
       // first[1]:  auth_asym_id
       // first[2]:  auth_comp_id
       // first[3]:  auth_seq_id
       // first[4]:  PDB_ins_code
       // second:    atom name
       std::list<std::pair<std::vector<std::string>, std::vector<std::string> > > _extra_atoms;

       std::list<_LINK>  _links;
       std::list<_SSBOND> _ssbonds;
       std::list<_SLTBRG> _sltbrgs;
       std::list<_BSPAIR> _bspairs;

       // key: pdb_chnid1_pdb_resnam1_pdb_resnum1_ins_code1_pdb_chnid2_pdb_resnam2_pdb_resnum2_ins_code2
       // value: first - pdb_atmnam1, second - pdb_atmnam2
       std::map<std::string, std::vector<std::pair<std::string, std::string> > > _covale_link_mapping;

       std::set<std::string> _amino_acid_backbone_atom_set;

       void clear();
       void _clear_contact();
       void _clear_geometry();
       void _clear_chirality_and_planarity();

       const std::vector<_MEAN_STD>& _getInterResidueBondStandard(RCSB::Residue* curr, const std::string& polymer_type, double& cutoff_value);
       void _checkInterResidueBond(const int& mol_id, RCSB::Residue* prev, RCSB::Residue* curr, const std::string& polymer_type,
                                   std::set<std::string>& found_allowed_bonds, const bool& check_linkage_only = false);
       const std::vector<_MEAN_STD>& _getInterResidueAngleStandard(const int& type, RCSB::Residue* curr, const std::string& polymer_type);
       void _checkInterResidueAngle(const int& type, const int& mol_id, RCSB::Residue* prev, RCSB::Residue* curr, const std::string& polymer_type,
                                    const std::set<std::string>& found_allowed_bonds, const std::map<std::string, std::string>&
                                    conformer_cis_trans_mapping);
       void _checkIntraResidueValues(const int& mol_id, RCSB::Residue* residue, const std::map<std::string, std::string>& conformer_cis_trans_mapping);
       void _checkIntraResidueValues(const int& type, const int& mol_id, const std::vector<_MEAN_STD>& standardValues, RCSB::Residue* residue,
                                     const std::map<std::string, std::string>& conformer_cis_trans_mapping);
       bool _checkBondOutlier(const int& mol_id, const _MEAN_STD& standard, const std::vector<RCSB::Atom*>& atoms, const std::string& cis_trans_type,
                              const std::string& linker_flag, const double& cutoff_distance = -1, const bool& check_linkage_only = false);
       void _checkAngleOutlier(const int& mol_id, const _MEAN_STD& standard, const std::vector<RCSB::Atom*>& atoms, const std::string&
                               cis_trans_type, const std::string& linker_flag);
       void _checkOutlier(const int& mol_id, const _MEAN_STD& standard, const std::string& cis_trans_type, const std::string& linker_flag,
                          const double& val, double& ept, const std::vector<RCSB::Atom*>& atoms, std::list<Value>& _outlier_list);
       void _checkInterResidueValues(const int& mol_id, const std::string& polymer_type, RCSB::Residue* prev, RCSB::Residue* curr);

       void _getPHITorsionAtomList(RCSB::Residue* prev, RCSB::Residue* curr, std::vector<std::vector<RCSB::Atom*> >& pair_lists);
       void _getPHITorsion(const int& mol_id, RCSB::Residue* prev, RCSB::Residue* curr,
                           const std::set<std::string>& found_allowed_bonds,
                           std::vector<Value>& phi);
       void _getPSITorsionAtomList(RCSB::Residue* prev, RCSB::Residue* curr, std::vector<std::vector<RCSB::Atom*> >& pair_lists);
       void _getPSITorsion(const int& mol_id, RCSB::Residue* prev, RCSB::Residue* curr,
                           const std::set<std::string>& found_allowed_bonds,
                           std::vector<Value>& psi);
       void _getOMEGATorsionAtomList(RCSB::Residue* prev, RCSB::Residue* curr, std::vector<std::vector<RCSB::Atom*> >& pair_lists);
       void _checkOMEGATorsion(const int& mol_id, RCSB::Residue* prev, RCSB::Residue* curr,
                               const std::set<std::string>& found_allowed_bonds,
                               std::map<std::string, std::string>&
                               pro_cis_trans_mapping);
       void _getCISOMEGATorsion(const int& mol_id, RCSB::Residue* prev, RCSB::Residue* curr,
                                const std::set<std::string>& found_allowed_bonds,
                                std::list<Value>& cis_peptide_list);
       void _checkPLANARTorsion(const int& mol_id, RCSB::Residue* prev, RCSB::Residue* curr,
                                const std::set<std::string>& found_allowed_bonds);
       void _checkPHIPSITorsion(const std::vector<Value>& phi_list,
                                const std::vector<Value>& psi_list);
       Value _getTorsionValue(const int& mol_id, const std::vector<RCSB::Atom*>& atoms);

       void _checkChirality(RCSB::Residue* residue, RCSB::Molecule* mol);
       void _checkPlanarity(const int& mol_id, RCSB::Residue* residue, PlanarityUtil& putil);

       std::string _getPDBId();
       void _updateValidationCategories(Block& block);
       void _update_pdbx_validate_close_contact(Block& block);
       void _update_pdbx_validate_symm_contact(Block& block);
       void _update_pdbx_validate_rmsd_bond(Block& block);
       void _update_pdbx_validate_rmsd_angle(Block& block);
       void _update_pdbx_validate_torsion(Block& block);
       void _update_pdbx_validate_peptide_omega(Block& block);
       void _update_pdbx_validate_main_chain_plane(Block& block);
       void _update_pdbx_validate_polymer_linkage(Block& block);
       void _update_pdbx_validate_planarity(Block& block);
       void _update_pdbx_validate_chirality(Block& block);
       void _get_software_ordinal(std::set<std::string>& software_ordinal_set);
       std::string _get_struct_descriptor(Block& block);

       /**
       **  Update pdbx_distant_solvent_atoms category
       **
       **  \param[out]: block - reference to data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_pdbx_distant_solvent_atoms(Block& block);

       void _update_pdbx_struct_special_symmetry(Block& block);

       void _update_pdbx_unobs_or_zero_occ_residues(Block& block);

       void _update_pdbx_unobs_or_zero_occ_atoms(Block& block);

       /**
       **  Reformat symmetry operator from 1555 to 1_555
       **
       **  \param[in]: old_symmetry - symmetry operator
       **
       **  \return reformated symmetry operator
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _reformat_symmetry(const std::string& old_symmetry);

       /**
       **  Insert _LINK into _links.
       **
       **  \param[in]: link - a _LINK value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_a_link(const _LINK& link);

       /**
       **  Remove _LINK from _links.
       **
       **  \param[in]: link - a _LINK value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _remove_a_link(const _LINK& link);

       /**
       **  Insert _LINK into _link_residue_set & _covale_link_mapping.
       **
       **  \param[in]: fstAtom - first linked atom
       **  \param[in]: sndAtom - second linked atom
       **  \param[in]: type    - link type
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_link_mapping(RCSB::Atom* fstAtom, RCSB::Atom* sndAtom, const std::string& type);

       /**
       **  Insert _SSBOND into _ssbonds.
       **
       **  \param[in]: ssbond - a _SSBOND value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_a_ssbond(const _SSBOND& ssbond);

       /**
       **  Insert _SLTBRG into _sltbrgs.
       **
       **  \param[in]: sltbrg - a _SLTBRG value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_a_sltbrg(const _SLTBRG& sltbrg);

       /**
       **  Insert _BSPAIR into _bspairs.
       **
       **  \param[in]: bspair - a _BSPAIR value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_a_bspair(const _BSPAIR& bspair);

       /**
       **  Read struct_conn category and insert _LINK/_SSBOND/_SLTBRG/_BSPAIR records
       **
       **  \param[in]: t - ISTable contains struct_conn values
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _read_struct_conn_category(ISTable* t);

       /**
       **  Read pdbx_struct_link category and insert _LINK records
       **
       **  \param[in]: t - ISTable contains pdbx_struct_link values
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _read_pdbx_struct_link_category(ISTable* t);

       /**
       **  Read struct_conn & pdbx_struct_link categories and input
       **  values into _links, _ssbonds, _sltbrgs & _bspairs 
       **
       **  \param[in]: block - reference to a coordinate data block 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _read_pdbx_struct_link_and_struct_conn(Block& block);

       /**
       **  Read _links & _ssbonds and create _defined_links mapping
       **
       **  \param Not applicable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _create_link_mapping();

       /**
       **  Insert link records into _defined_links
       **
       **  \param[in]: links - list of link records
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_defined_links(const std::list<_LINK>& links);

       /**
       **  Insert two linked atoms into _defined_links
       **
       **  \param[in]: atom1 - first atom
       **  \param[in]: atom2 - second atom
       **  \param[in]: dist  - distance between two atoms
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_links(RCSB::Atom* atom1, RCSB::Atom* atom2, const std::string& dist);

       /**
       **  Check if there is a link between two atoms
       **
       **  \param[in]: atom1 - first atom
       **  \param[in]: atom2 - second atom
       **  \param[in]: dist  - distance between two atoms
       **
       **  \return true -  if link exists between two atoms
       **  \return false - if link does not exists between two atoms
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _is_a_link(RCSB::Atom* atom1, RCSB::Atom* atom2, const double& dist);

       /**
       **  Check out of range waters
       **
       **  \param[in]: mol - point to a Molecule
       **  \param[in]: water_lists - list of waters found in that Molecule
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _check_distant_waters(RCSB::Molecule* mol, const std::list<RCSB::Residue*>& water_lists);

       /**
       **  Check out of range ligands
       **
       **  \param[in]: mol - point to a Molecule
       **  \param[in]: Mol_ID - Molecule ID
       **  \param[in]: residue_lists - list of ligands found in that Molecule
       **  \param[out]: distant_residue_lists - list of out of range ligands
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _check_distant_nonpolymer_residues(RCSB::Molecule* mol, const std::string& Mol_ID, const std::list<RCSB::Residue*>& residue_lists,
                                               std::list<std::pair<std::string, RCSB::Residue*> >& distant_residue_lists);

       void _get_missing_or_zero_occupancy_residues_or_atoms();

       void _finding_missing_and_extra_atoms(const bool& exclude_hydrogen = true);

       int _cal_z_value();

       bool _cal_matthew_and_solvent(std::vector<double>& return_values);

       /**
       **  Checking polymer linkage between adjacent residues 
       **
       **  \param[in]: chain  - pointer to polymer chain
       **  \param[in]: Mol_ID - Model serial number
       **  \param[in]: check_gap_distance - flag for checking gap distance
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _check_polymer_linkage(RCSB::Chain* chain, const int& Mol_ID, const bool& check_gap_distance = false);

       /**
       **  Checking polymer linkage between adjacent residues 
       **
       **  \param[in]: chain  - pointer to polymer chain
       **  \param[in]: Mol_ID - Model serial number
       **  \param[in]: check_gap_distance - flag for checking gap distance
       **  \param[out]: warning_messages - warning messages
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _run_polymer_linkage_checking(RCSB::Chain* chain, const int& Mol_ID, const bool& check_gap_distance, std::vector<std::string>& warning_messages);

       /**
       **  Checking polymer linkage between adjacent residue pairs
       **
       **  \param[in]: Mol_ID     - Model serial number
       **  \param[in]: chain_type - polymer chain type
       **  \param[in]: prev       - first residue
       **  \param[in]: curr       - second residue
       **  \param[in]: gap_flag   - flag to indicate if there is gap between two residues
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _run_adjacent_linkage_checking(const int& Mol_ID, const std::string& chain_type, RCSB::Residue* prev, RCSB::Residue* curr, const bool& gap_flag);

       /**
       **  Calculating possible terminal atom pairs distance
       **
       **  \param[in]: prev       - first residue
       **  \param[in]: curr       - second residue
       **  \param[in]: tailAtoms  - tail terminal atom set
       **  \param[in]: headAtoms  - head terminal atom set
       **  \param[out]: global_atom_pair_map - result atom distance pair map
       **  \param[out]: checked_pair_set - checked atom name pair set
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cal_adjacent_linkage_distances(RCSB::Residue* prev, RCSB::Residue* curr, const std::set<std::string>& tailAtoms, const std::set<std::string>& headAtoms,
                                            std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > >& global_atom_pair_map,
                                            std::set<std::string>& checked_pair_set);

       /**
       **  Calculating linked atom pairs distance from LINK records
       **
       **  \param[in]: prev       - first residue
       **  \param[in]: curr       - second residue
       **  \param[out]: global_atom_pair_map - result atom distance pair map
       **  \param[out]: checked_pair_set - checked atom name pair set
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cal_adjacent_link_record_distances(RCSB::Residue* prev, RCSB::Residue* curr, std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > >&
                                                global_atom_pair_map, std::set<std::string>& checked_pair_set);

       /**
       **  Calculating a terminal atom pair distance
       **
       **  \param[in]: prev       - first residue
       **  \param[in]: curr       - second residue
       **  \param[in]: tailAtom   - tail terminal atom
       **  \param[in]: headAtom   - head terminal atom
       **  \param[out]: atom_pair_map - result atom distance pair map
       **  \param[out]: checked_pair_set - checked atom name pair set
       **
       **  \return shortest distance
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       double _cal_adjacent_linkage_distance(RCSB::Residue* prev, RCSB::Residue* curr, const std::string& tailAtom, const std::string& headAtom,
                                            std::multimap<double, std::pair<bool, std::pair<RCSB::Atom*, RCSB::Atom*> > >& atom_pair_map,
                                            std::set<std::string>& checked_pair_set);

       /**
       **  Checking gap distance between residue pairs
       **
       **  \param[in]: Mol_ID - Model serial number
       **  \param[in]: prev   - first residue list
       **  \param[in]: curr   - second residue list
       **  \param[in]: gap_count - Gap number
       **
       **  \return problematical residue name pair
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _run_gap_distance_checking(const int& Mol_ID, const std::vector<RCSB::Residue*>& prev, const std::vector<RCSB::Residue*>& curr,
                                              const int& gap_count);

       /**
       **  Get regular and canonical one letter codes
       **
       **  \param[out]: code     - regular one letter code
       **  \param[out]: code_can - canonical one letter code
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_oneletter_code(std::string& code, std::string& code_can);

       /**
       **  Write ATOM information
       **
       **  \param[in]: atom_list - atom list
       **
       **  \return ATOM content
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _write_atom_list(const std::list<RCSB::Atom*>& atom_list);

       /**
       **  Check cell information for non crystal entry
       **
       **  \param[in]: block - reference to data block
       **
       **  \return error message
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _check_cell_info(Block& block, const bool& comment_flag);

       std::string _write_link_list(std::list<std::vector<std::string> >& links);

       std::string _get_zero_occupancy_atom_list_summary(std::list<RCSB::Atom*>& zero_occupancy_atom_list);

       std::list<Value> _get_special_position_atoms(const bool& pdb_format_flag = true);
       double _cal_distance_between_residues(RCSB::Residue* res1, const std::string& F_atom, RCSB::Residue* res2, const std::string& S_atom);
       void _get_atom_pair_list(const std::vector<std::pair<RCSB::Residue*, std::string> >& residue_atom_name_pair,
                                std::vector<std::vector<RCSB::Atom*> >& pair_lists);

       void _insert_covalent_bonding_info(const std::string& resIdx, const std::string& fstAtomIdx, const std::string& sndAtomIdx, const std::string&
                                          dist, std::map<std::string, std::vector<std::vector<std::string> > >& covalent_bonding);

       /**
       **  Checking dbref
       **
       **  \param[out]: block - reference to data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _checking_dbref(Block& block);

       /**
        **  Checking assembly categories
        **
        **  \param[out]: block - reference to data block
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
       bool _checking_assembly(Block& block);

       bool _check_struct_title(const std::string& title);

       /**
        **  Get pdbx_leaving_atom_flag for covalent LINK record
        **
        **  \param[in]: atom1 - first atom in LINK record
        **  \param[in]: aton2 - second atom in LINK record
        **
        **  \return both - Both atoms have displaced leaving atoms
        **  \return one  - One of them has displaced leaving atom
        **  \return none - None of them has displaced leaving atom
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
       std::string _get_leaving_flag(const RCSB::Atom* atom1, const RCSB::Atom* atom2);

       /**
        **  Update pdbx_leaving_atom_flag for existing covalent LINK records
        **
        **  \param: Not applicable
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
       void _update_leaving_atom_flag();

       /**
        **  Update struct_conn category using values from
        **  _links, _ssbonds, _sltbrgs & _bspairs 
        **
        **  \param[in]: block - reference to a coordinate data block 
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
       void _update_struct_conn(Block& block);

       /**
        **  Reorder _links & _ssbonds based on chain IDs & residue numbers
        **
        **  \param[out]: chain_order - chain ID order
        **  \param[out]: links       - list of link values
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
       void _re_order_links(const std::map<std::string, unsigned int>& chain_order, std::list<_LINK>& links);

       /**
        **  Update struct_conn category using values from _links, _ssbonds, _sltbrgs & _bspairs 
        **
        **  \param[out]: t          - point to struct_conn table
        **  \param[out]: type_array - unique struct_conn.type values stored in vector
        **  \param[out]: type_set   - unique struct_conn.type values stored in set
        **  \param[in]: links       - list of link values
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
       void _update_struct_conn_table(ISTable* t, std::vector<std::string>& type_array, std::set<std::string>& type_set, const std::list<_LINK>& links);

       /**
        **  Get missingbackbone atoms
        **
        **  \param[in]: res      - pointer to a residue
        **  \param[out]: atoms   - missing/extra atom name list
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
        void _get_missing_backbone_atoms(RCSB::Residue* res, std::vector<std::string>& atoms);
       
       /**
        **  Checking missing/extra atoms
        **
        **  \param[in]: mol_ID   - molecule ID
        **  \param[in]: res      - pointer to a residue
        **  \param[in]: atoms    - missing/extra atom name list
        **  \param[out]: message - output missing/extra atom text
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
        void _check_missing_extra_atoms(const std::string& mol_ID, RCSB::Residue* res, const std::vector<std::string>& atoms, std::string& message);

       /**
        **  Checking pdbx_database_related category
        **
        **  \param[out]: block - reference to data block
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
        void _checking_pdbx_database_related(Block& block, const bool& release_flag = false);

       /**
        **  Get struct_ncs_dom_lim category's item names
        **
        **  \param[in]: t - ISTable pointer to "struct_ncs_dom_lim" table
        **  \param[out]: all_beg_items - all beg_* item names
        **  \param[out]: beg_items - beg_* item names
        **  \param[out]: all_end_items - all end_* item names
        **  \param[out]: end_items - end_* item names
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
        void _get_struct_ncs_dom_lim_item_names(ISTable* t, std::vector<std::string>& all_beg_items, std::vector<std::vector<std::string> >& beg_items,
                                                std::vector<std::string>& all_end_items, std::vector<std::vector<std::string> >& end_items);


       /**
        **  Find the Begin & End residue(s) from struct_ncs_dom_lim category
        **
        **  \param[in]: t - ISTable pointer to "struct_ncs_dom_lim" table
        **  \param[in]: row - table row number
        **  \param[in]: items - item name lists
        **  \param[in]: all_items - all item name lists
        **  \param[out]: msg - message for missing residue
        **
        **  \return Not applicable
        **
        **  \pre None
        **
        **  \post None
        **
        **  \exception: None
        */
        RCSB::Residue* _find_struct_ncs_dom_lim_residue(ISTable* t, const unsigned int& row, const std::vector<std::vector<std::string> >& items,
                                                        const std::vector<std::string>& all_items, std::string& msg);
};

#endif
