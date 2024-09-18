/*
FILE:     AnnotationObj.h
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
#ifndef _H_ANNOTATIONOBJ_H_
#define _H_ANNOTATIONOBJ_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "BasePairParamDef.h"
#include "Contact.h"
#include "MovingAtom.h"
#include "SecondStruct.h"
#include "ValidateObj.h"
#include "xtal.h"

class AnnotationObj: public ValidateObj {
   public:

       /**
       **  Constructs AnnotationObj
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
       AnnotationObj();

       /**
       **  Destructs AnnotationObj
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
       ~AnnotationObj();

       /**
       **  Set _get_skip_option_flag to true 
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
       void set_get_skip_option_flag();

       /**
       **  Set _update_xray_v4_to_v5_flag to true 
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
       void set_update_xray_v4_to_v5_flag();

       /**
       **  Set _update_em_v4_to_v5_flag to true 
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
       void set_update_em_v4_to_v5_flag();

       /**
       **  Set _update_nmr_v4_to_v5_flag to true 
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
       void set_update_nmr_v4_to_v5_flag();

       /**
       **  Set _carbohydrate_gmml_track_mode_flag to true 
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
       void set_carbohydrate_gmml_track_mode_flag();

       /**
       **  Set _include_date_original_flag to true 
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
       void set_include_date_original_flag();

       /**
       **  Set _intra_molecular_connectivity_flag to true 
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
       void set_intra_molecular_connectivity_flag();

       /**
       **  Set PDB2Glycan program path
       **
       **  \param[in]: path - the full path name for PDB2Glycan program
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void set_pdb2glycan_program(const std::string& program);

       /**
       **  Update various cif categories
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
       void cif_update();

       /**
       **  Merge EBI status info. into pdbx_database_status category
       **
       **  \param[in]: status_file - status information cif file
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void merge_ebi_status(const std::string& status_file);

       /**
       **  Update initial deposit date and release date in pdbx_database_status category from database_PDB_rev category 
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
       void update_deposit_and_release_dates();

       /**
       **  Add dep ID and reset entry ID to dep ID for legacy entries
       **
       **  \param[in]: depid - dep ID
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void add_depid(const std::string& depid);

       /**
       **  Update various EM-related cif categories
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
       void cif_update_EM_V4_V5();

       /**
       **  Convert old database_PDB_rev & pdbx_version categories into pdbx_audit_revision* categories
       **
       **  \param[in]: depid   - deposition ID
       **  \param[in]: version - internal version number
       **
       **  \return true if it is successful
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool old_to_new_audit_revision_conversion(const std::string& depid = "", const std::string& version = "");

       /**
       **  Read all structural features from cif categories into memory.
       **
       **  \param[in]: skip_flag - Skip coordinate-related structural features
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void ReadStructuralFeatures(const bool& skip_flag = false);

       /**
       **  Write out all structural features to cif categories
       **
       **  \param[in]: update_z_value_flag - update Z value flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void WriteStructuralFeatures(const bool& update_z_value_flag = true, const bool& skip_flag = false, const bool& update_struct_flag = false);

       /**
       **  Update pdbx_nmr_sample_details & pdbx_nmr_exptl_sample categories
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
       void Update_NMR_sample_details();

       /**
       **  Assign seqeuntial entity ID to each chain
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
       void AssignEntityId();

       /**
       **  Assign asym ID to each chain
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
       void AssignAsymId();

       /**
       **  Update polymer sequence information based on Sequence Module output
       **
       **  \param[in]: assignfile - assignment file name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void UpdateSeqModuleInfo(const std::string& assignfile);

       /**
       **  Calculate all LINK & SSBOND recrods from coordinates
       **
       **  \param[in]: type - 1 (FIND_LINK) for LINK record
       **  \param[in]: type - 2 (FIND_SSBOND) for SSBOND record
       **  \param[in]: type - 3 (FIND_LINK | FIND_SSBOND) for both LINK & SSBOND records
       **  \param[in]: force_flag - force recalculation flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void Calculate_Link_and_SSBond(const int& type, const bool& force_flag);

       /**
       **  Convert selected close contact(s) into link record(s).
       **
       **  \param[in]: datafile - selected close contact list file
       **
       **  \return true if successful
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool Converting_Close_Contact_to_Link(const std::string& datafile);

       /**
       **  Remove selected link(s) from link record(s).
       **
       **  \param[in]: datafile - selected link list file
       **
       **  \return true if successful
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool Removing_Covalent_Link(const std::string& datafile);

       /**
       **  Update struct_conn & pdbx_struct_link categories using
       **  values from _links, _ssbonds, _sltbrgs & _bspairs 
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
       void Update_StructLink_and_StructConn();

       /**
       **  Call BasePairInfo class to calculate all base pair
       **  information for nucleic acid chain(s).
       **
       **  \param[in]: force_flag - force recalculation flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void CalculateBasePairInfo(const bool& force_flag);

       /**
       **  Update ndb_struct_na_base_pair & ndb_struct_na_base_pair_step categories
       **  using values from _basebase_params & _interbase_params
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
       void Update_Base_Pair_and_Step();

       /**
       **  Call Promotif class to calculate secondary structures for peptide chain(s).
       **
       **  \param[in]: supportfile - input file name which contains user defined sheet
       **                            topology information
       **  \param[in]: force_flag  - force recalculation flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void CalculateSecondStruct(const std::string& supportfile, const bool& force_flag);

       /**
       **  Update struct_conf, struct_sheet, struct_sheet_order, struct_sheet_range
       **  & pdbx_struct_sheet_hbond categories using values from _helix & _sheet
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
       void Update_SecondStruct();

       /**
       **  Moving waters closer to polymers based crystallographic symmetry
       **
       **  \param[in]: skip_flag - skip assignment flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void Moving_Waters(const bool& skip_flag = false);

       /**
       **  Assign non-polymer and water PDB chain IDs based on relative distance to polymer
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
       // void AssignPDBChainIDandNumbering();

       /**
       **  Update struct_mon_prot_cis category
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
       void UpdateStructMonProtCis();

       /**
       **  Update coordinates categories
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
       void UpdateCoordinateCategory();

       /**
       **  Merge Struct SITE records
       **
       **  \param[in]: sitefile - site record cif file name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void merge_struct_sites(const std::string& sitefile);

       /**
       **  Merge ligand assignments
       **
       **  \param[in]: assignfile - ligand assignment file name
       **
       **  \return true  - if successful
       **          false - else
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool updateInstance(const std::string& assignfile);

       /**
       **  Update CAVEAT
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
       void updateCAVEAT();

       /**
       **  Correct standard residue's atom name 
       **
       **  \param[in]: UpperCase_Flag - 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void CorrectAtomName(const bool& UpperCase_Flag);

       /**
       **  Delete/rename terminal atoms in middle of sequence
       **
       **  \param[in]: rename_flag - 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void correct_terminal_atoms(const bool& rename_flag);

       /**
       **  Move best representive model to the first model
       **
       **  \param: Not applicable
       **
       **  \return true  - if successful
       **          false - else
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool reorder_model();

       /**
       **  Move best representive model to the first model
       **
       **  \param[in]: mol_id - best representive model ID number
       **
       **  \return true  - if successful
       **          false - else
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool reorder_model(const std::string& mol_id);

       /**
       **  Check if it is a big entry
       **
       **  \param: Not applicable
       **
       **  \return true/false
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool is_large_entry();

       /**
       **  Extra clean up steps for V4 -> V5 conversion
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
       void V4_V5_Cleanup();

       /**
       **  Convert obsolete & withdrawal entries to V5 format
       **
       **  \param[in]: depid - deposition ID
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void convert_obsolete_and_withdrawal_entry(const std::string& depid);

       /**
       **  Re-assign atom name from _chem_comp_atom.atom_id to _chem_comp_atom.pdbx_stnd_atom_id
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
       void re_assign_residue_to_standard_name();

       /**
       **  Annotating carbohydrate entities
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
       void Carbohydrate_annotation();

       /**
       **  Output public cif categories and items only
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
       void Output_public_cif_only();

       /**
       **  Update oligosaccharide related PRD information.
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
       void UpdateOligosaccharidePrdInfo();

       /**
       **  Add default as is assembly information.
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
       void AddDefaultAssembly();

   protected:
       typedef struct {
             int major_revision;
             int minor_revision;
             std::string data_content_type;
             std::string revision_date;
             std::string internal_deposition_id;
             std::string internal_version;
             std::set<std::string> revision_group;
             // vector[0]: provider
             // vector[1]: type
             // vector[2]: description
             std::vector<std::vector<std::string> > revision_details;
             std::set<std::string> revision_category;
             std::set<std::string> revision_item;
       } _Audit_Revision;

       bool _merge_seq_module_flag;
       bool _coordinate_merge_flag;
       bool _get_skip_option_flag;
       bool _update_xray_v4_to_v5_flag;
       bool _update_em_v4_to_v5_flag;
       bool _update_nmr_v4_to_v5_flag;
       bool _carbohydrate_gmml_track_mode_flag;
       bool _numbering_based_order_flag;
       bool _branch_annotated_flag;
       bool _include_date_original_flag;
       bool _intra_molecular_connectivity_flag;
       std::string _pdb2glycanProgram;
       std::vector<std::string> _classifications;
       std::list<_BASEBASE_PARAMS> _basebase_params;
       std::list<_INTERBASE_PARAMS> _interbase_params;
       std::list<_MODRES> _modres;
       std::list<_HELIX> _helix;
       std::list<_HELIX> _turn;
       std::list<_SHEET> _sheet;
       std::list<_SITE> _site;
       std::list<_DBREF> _dbrefs;
       std::list<_SEQADV> _seqadvs;
       std::list<_MOVING_ATOM> _moving_water;

       // key: sugar_residue_id chnid_resname_resnum_ins_code
       // value[0]: glyco_site_residue chnid
       // value[1]: glyco_site_residue resname
       // value[2]: glyco_site_residue resnum
       // value[3]: glyco_site_residue ins_code
       std::map<std::string, std::vector<std::string> > _glyco_site_linkage_info;

       // key: residue_id chnid_resname_resnum_ins_code
       // value: linked polymer chain ID
       std::map<std::string, std::string> _nonpolymer_group_mapping;

       std::map<std::string, std::string> _new_chain_entity_id_mapping;

       // PRD information
       // map key: PRD ID
       // map value: key/value pair where keys are:
       //                name, type, class, details
       std::map<std::string, std::map<std::string, std::string> > _prd_features;

       // PRD instance information
       // pair.first: PRD ID
       // pair.second.first: polymer chain index
       // pair.second.second: non-polymer residue index
       std::vector<std::pair<std::string, std::pair<std::set<int>, std::set<int> > > > _prd_instances;

       std::set<std::string> _mandatory_updated_categories;
       std::string _default_exp_type;
       int _experimental_method_count;
       std::string _default_diffrn_id;
       std::map<std::string, std::string> _scattering_type_mapping;

       typedef struct {
              // int residue_index;
              bool is_instance;
              std::string feature_type;
              std::string details;
              std::vector<RCSB::Atom*> atoms;
       } _PDBX_ENTITY_INSTANCE_FEATURE;

       std::list<_PDBX_ENTITY_INSTANCE_FEATURE> _pdbx_entity_instance_feature;

       typedef struct {
              RCSB::Atom* atom;
              // vector[0]: pre_group_PDB
              // vector[1]: pre_auth_atom_id
              // vector[2]: pre_auth_comp_id
              // vector[3]: pre_auth_asym_id
              // vector[4]: pre_auth_seq_id
              // vector[5]: pre_auth_alt_id
              // vector[6]: pre_PDB_ins_code
              // vector[7]: pre_occupancy
              std::vector<std::string> pre_data;
       } _PDBX_REMEDIATION_ATOM_SITE_MAPPING;

       std::list<_PDBX_REMEDIATION_ATOM_SITE_MAPPING> _pdbx_remediation_atom_site_mapping;

       std::set<std::string> _branch_residue_set, _glyco_link_set;

       // pair.first: meta data map
       // pair.second: chain_index
       std::list<std::pair<std::map<std::string, std::string>, int> > _refine_ls_restr_ncs_list;

       // pair.first: meta data map
       // pair.second.first: begin_res_index
       // pair.second.second: end_res_index
       std::list<std::pair<std::map<std::string, std::string>, std::pair<int, int> > > _struct_ncs_dom_lim_list;

       /**
       **  Clear all data storage
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
       void clear();

       /**
       **  Clear data storage
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
       void _clear_annotation_obj_allocated_memories();

       /**
       **  Find LINK and/or SSBOND recrods from contact list
       **
       **  \param[in]: type - 1 (FIND_LINK) for LINK record
       **  \param[in]: type - 2 (FIND_SSBOND) for SSBOND record
       **  \param[in]: type - 3 (FIND_LINK | FIND_SSBOND) for both LINK & SSBOND records
       **  \param[in]: contact - contact information calculated from coordinates
       **  \param[in]: short_protein_chain_id_set - short protein chain's PDB ID set
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _find_Link_and_SSBond(const int& type, const std::vector<CONTACT>& contact, const std::set<std::string>& short_protein_chain_id_set);

       /**
       **  Find LINK and/or SSBOND recrods from contact list
       **
       **  \param[in]: bond_type - 1 for covalent bond
       **  \param[in]: bond_type - 2 for metal coordination
       **  \param[in]: contact     - contact information
       **  \param[in]: short_protein_chain_id_set - short protein chain's PDB ID set
       **  \param[out]: glycosylation_type - Glycosylation type
       **
       **  \return 1 (FIND_LINK)   - for LINK record
       **  \return 2 (FIND_SSBOND) - for SSBOND record
       **  \return 0               - none of the above
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _assign_contact_type(const int& bond_type, const CONTACT& contact, const std::set<std::string>& short_protein_chain_id_set,
                                std::string& glycosylation_type);

       /**
       **  Check if two linked atoms is a glyco like linkage
       **
       **  \param[in]: atom1 - first linked atom
       **  \param[in]: atom2 - second linked atom
       **
       **  \return ture - if two atoms are defined in allowed linakge list
       **  \return false - otherwise
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _find_glyco_like_linkage(RCSB::Atom* atom1, RCSB::Atom* atom2);

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
       // std::string _get_leaving_flag(const RCSB::Atom* atom1, const RCSB::Atom* atom2);

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
       // void _update_leaving_atom_flag();

       /**
       **  Update pdbx_struct_link category using value(s) from _links
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
       void _update_pdbx_struct_link(Block& block);

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
       // void _update_struct_conn(Block& block);

       /**
       **  Update struct_conn category using values from
       **  _links, _ssbonds, _sltbrgs & _bspairs 
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
       // void _update_struct_conn_table(ISTable* t, std::vector<std::string>& type_array, std::set<std::string>& type_set, const std::list<_LINK>& links);

       /**
       **  Read _basebase_params from ndb_struct_na_base_pair category
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
       void _read_ndb_struct_na_base_pair(Block& block);

       /**
       **  Update ndb_struct_na_base_pair category using values from _basebase_params
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
       void _update_ndb_struct_na_base_pair(Block& block);

       /**
       **  Read _interbase_params from ndb_struct_na_base_pair_step category
       **  _interbase_params
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
       void _read_ndb_struct_na_base_pair_step(Block& block);

       /**
       **  Update ndb_struct_na_base_pair_step category using values from
       **  _interbase_params
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
       void _update_ndb_struct_na_base_pair_step(Block& block);

       /**
       **  Update ndb_struct_conf_na category using values from _classifications
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
       void _update_ndb_struct_conf_na(Block& block);

       /**
       **  Read _helix & _turn from struct_conf category
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
       void _read_struct_conf(Block& block);

       /**
       **  Update struct_conf category using values from _helix & _turn
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
       void _update_struct_conf(Block& block);

       /**
       **  Update struct_conf category using values from _helix & _turn
       **
       **  \param[out]: t - point to struct_conf table
       **  \param[in]: type - HELX_P/TURN_P
       **  \param[in]: helices - list of helix values
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_struct_conf_table(ISTable* t, const std::string& type, const std::list<_HELIX>& helices);

       /**
       **  Read _sheet records from struct_sheet, struct_sheet_order, struct_sheet_range
       **      & pdbx_struct_sheet_hbond categories
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
       void _read_struct_sheet(Block& block);

       /**
       **  Update struct_sheet, struct_sheet_order, struct_sheet_range &
       **  pdbx_struct_sheet_hbond categories using values from _sheet
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
       void _update_sheet_categories(Block& block);

       /**
       **  Update struct_sheet category
       **
       **  \param[out]: t - point to struct_sheet table
       **  \param[in]: sheetID - sheet ID
       **  \param[in]: num_strands - number of strands in the sheet
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_sheet(ISTable *t, const std::string& sheetID, const std::string& num_strands);

       /**
       **  Get data for struct_sheet_range category
       **
       **  \param[in]: sheetID - sheet ID
       **  \param[in]: strands - strand list
       **  \param[out]: range_data - data for struct_sheet_range category
       **
       **  \return true if found data for struct_sheet_range category and no errors
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _get_range_data(const std::string& sheetID, const std::vector<_SHEET_STRAND>& strands, std::vector<std::map<std::string, std::string> >& range_data);

       /**
       **  Get data for struct_sheet_order & pdbx_struct_sheet_hbond categories
       **
       **  \param[in]: sheetID - sheet ID
       **  \param[in]: strand_orders - strand order list
       **  \param[out]: order_data - data for struct_sheet_order category
       **  \param[out]: hbond_data - data for pdbx_struct_sheet_hbond category
       **
       **  \return true if no errors found
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _get_order_and_hbond_data(const std::string& sheetID, const std::vector<_SHEET_TOPOLOGY>& strand_orders, std::vector<std::map<std::string,
                                      std::string> >& order_data, std::vector<std::map<std::string, std::string> >& hbond_data);
       /**
       **  Update struct_sheet_range category
       **
       **  \param[out]: t - point to struct_sheet_range table
       **  \param[in]: id - serial ID
       **  \param[in]: sheet - single strand record
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       // void _insert_sheet_range(ISTable *t, const int& id, const _SHEET& sheet);

       /**
       **  Update struct_sheet_order category
       **
       **  \param[out]: t - point to struct_sheet_order table
       **  \param[in]: id - serial ID
       **  \param[in]: sheetID - sheet ID
       **  \param[in]: sense - sense of strand with respect to previous strand
       **                      in the sheet. 1 if parallel, -1 if anti-parallel
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       // void _insert_sheet_order(ISTable *t, const int& id, const std::string& sheetID, const int& sense);

       /**
       **  Update pdbx_struct_sheet_hbond category
       **
       **  \param[out]: t - point to pdbx_struct_sheet_hbond table
       **  \param[in]: id - serial ID
       **  \param[in]: sheet - single strand record
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       // void _insert_sheet_hbond(ISTable *t, const int& id, const _SHEET& sheet);

       /**
       **  Get non polymers to polymer chain ID mapping from _link records
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
       // void _get_nonpolymer_group_mapping();

       /**
       **  Moving waters to a Molecule
       **
       **  \param[out]: mol - point to a Molecule
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       // void _reposition_waters(RCSB::Molecule* mol);

       /**
       **  Moving waters to a Molecule with defined water residues list
       **
       **  \param[out]: mol - point to a Molecule
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
       // void _reposition_waters(RCSB::Molecule* mol, const std::list<RCSB::Residue*>& water_lists);

       /**
       **  Update encap coordinate table
       **
       **  \param[out]: t - point to encap coordinate table
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateEncapCoordinateas(ISTable *t);

       /**
       **  Update atom_site & atom_site_anisotrop tables
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
       void _updateAtomSite(Block& block);

       /**
       **  Update pdbx_solvent_atom_site_mapping table
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
       void _updateSolventAtomSiteMapping(Block& block);

       /**
       **  Update atom_site_anisotrop table
       **
       **  \param[in]: items - item names for atom_site_anisotrop category
       **  \param[in]: atom_id - atom serial number
       **  \param[in]: atom - pointer to atom record
       **  \param[out]: t - pointer to atom_site_anisotrop table
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_atom_site_anisotrop(const std::vector<std::string>& items, const int& atom_id, RCSB::Atom* atom, ISTable* t);

       /**
       **  Update _cell.Z_PDB value
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
       void _update_Z_value(Block& block);

       /**
       **  Map sequence alignments outputed from Sequence Module with
       **    coordinates
       **
       **  \param[out]: seq_coord_mapping - auth sequence vs. coordinate residue
       **                                   alignment mapping
       **  \param[out]: chainid_chain_mapping - PDB chain ID vs. polymer chain mapping
       **  \param[out]: chainid_order_mapping - PDB chain ID vs. order mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateMoleculeWithSeqCoordMapping(std::map<std::string, std::vector<std::vector<std::string> > >& seq_coord_mapping,
               std::map<std::string, RCSB::Chain*>& chainid_chain_mapping, std::map<int, std::string>& chainid_order_mapping);

       /**
       **  Update source information
       **
       **  \param[in]: sources - source information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateSourceInfo(const std::map<std::string, std::vector<std::map<std::string, std::string> > >& sources);

       /**
       **  Update name information
       **
       **  \param[in]: names - name information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateNameInfo(const std::map<std::string, std::vector<std::string> >& names);

       /**
       **  Add self reference to struct_ref categories
       **
       **  \param[in]: pdb_id - PDB entry ID
       **  \param[in]: align_id - strut_ref_seq.align_id
       **  \param[out]: struct_ref -
       **  \param[out]: struct_ref_seq -
       **  \param[out]: self_ref_chains - self reference chain IDs
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateSelfStructRef(const std::string& pdb_id, int& align_id, std::map<std::string, std::map<std::string, std::string> >& struct_ref,
                 std::map<std::string, std::vector<std::map<std::string, std::string> > >& struct_ref_seq, std::set<std::string>& self_ref_chains);

       /**
       **  Update struct_ref_seq category
       **
       **  \param[out]: block - reference to data block
       **  \param[in]: struct_ref_seq -
       **  \param[in]: chainid_order_mapping - ordered PDB chain IDs
       **  \param[out]: new_old_ref_id_mapping - new vs old ref_id mapping
       **  \param[out]: new_old_align_id_mapping - new vs old align_id mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _outputStructRefSeq(Block &block, const std::map<std::string, std::vector<std::map<std::string, std::string> > >& struct_ref_seq,
                                const std::map<int, std::string>& chainid_order_mapping, std::map<int, std::string>& new_old_ref_id_mapping,
                                std::map<int, std::string>& new_old_align_id_mapping);

       /**
       **  Update struct_ref category
       **
       **  \param[out]: block - reference to data block
       **  \param[in]: struct_ref -
       **  \param[in]: new_old_ref_id_mapping - new vs old ref_id mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _outputStructRef(Block &block, const std::map<std::string, std::map<std::string, std::string> >& struct_ref,
                             const std::map<int, std::string>& new_old_ref_id_mapping);

       /**
       **  Update struct_ref_seq_dif category
       **
       **  \param[out]: block - reference to data block
       **  \param[in]: struct_ref_seq_dif -
       **  \param[in]: new_old_align_id_mapping - new vs old align_id mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _outputStructRefSeqDif(Block &block, const std::map<std::string, std::vector<std::map<std::string, std::string> > >& struct_ref_seq_dif,
                                   const std::map<int, std::string>& new_old_align_id_mapping);

       /**
       **  Update _modres record
       **
       **  \param[in]: modres_list - modified residue information
       **                            vector[0]: PDB_chain ID
       **                            vector[1]: Residue Name
       **                            vector[2]: Residue Number
       **                            vector[3]: Insertion Code
       **                            vector[4]: Parent Residue Name
       **                            vector[5]: Comment
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateModRes(const std::vector<std::vector<std::string> >& modres_list);

       /**
       **  Update entity IDs for all related categories
       **
       **  \param[in]: old_new_entity_mapping - old vs. new Entity IDs mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateEntityIDs(const std::map<std::string, std::string>& old_new_entity_mapping);

       /**
       **  Update struct_ref & struct_ref_seq categories
       **
       **  \param[in]: chain_entity_id_mapping - PDB chain ID vs. new Entity IDs mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateStructRef(const std::map<std::string, std::string>& chain_entity_id_mapping);

       /**
       **  Read MODRES from pdbx_struct_mod_residue category
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
       void _read_modres(Block& block);

       /**
       **  Write MODRES to pdbx_struct_mod_residue category
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
       void _write_modres(Block& block);

       /**
       **  Read PRD features from pdbx_molecule_features category
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
       void _read_prd_features(Block& block);

       /**
       **  Write PRD features to pdbx_molecule_features category
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
       void _write_prd_features(Block& block);

       /**
       **  Read PRD instances from pdbx_molecule category
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
       void _read_prd_instances(Block& block);

       /**
       **  Write PRD instances to pdbx_molecule category
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
       void _write_prd_instances(Block& block);

       /**
       **  Read user input features from pdbx_entity_instance_feature category
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
       void _read_pdbx_entity_instance_feature(Block& block);

       /**
       **  Write user input features to pdbx_entity_instance_feature category
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
       void _write_pdbx_entity_instance_feature(Block& block);

       /**
       **  Read remediation atom mapping changes from pdbx_remediation_atom_site_mapping category
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
       void _read_pdbx_remediation_atom_site_mapping(Block& block);

       /**
       **  Write remediation atom mapping changes to pdbx_remediation_atom_site_mapping category
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
       void _write_pdbx_remediation_atom_site_mapping(Block& block);

       /**
       **  Read DBREF instances from struct_ref_seq category
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
       void _read_struct_ref_seq(Block& block);

       /**
       **  Write DBREF instances to struct_ref_seq category
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
       void _write_struct_ref_seq(Block& block);

       /**
       **  Read SEQADV instances from struct_ref_seq_dif category
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
       void _read_struct_ref_seq_dif(Block& block);

       /**
       **  Write SEQADV instances to struct_ref_seq_dif category
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
       void _write_struct_ref_seq_dif(Block& block);

       /**
       **  Read re-positioned waters from pdbx_solvent_atom_site_mapping category
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
       void _read_pdbx_solvent_atom_site_mapping(Block& block);

       /**
       **  Read _refine_ls_restr_ncs_list from refine_ls_restr_ncs category
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
       void _read_refine_ls_restr_ncs_list(Block& block);

       /**
       **  Write _refine_ls_restr_ncs_list to refine_ls_restr_ncs category
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
       void _write_refine_ls_restr_ncs_list(Block& block);

       /**
       **  Read _struct_ncs_dom_lim_list from struct_ncs_dom_lim category
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
       void _read_struct_ncs_dom_lim_list(Block& block);

       /**
       **  Write _struct_ncs_dom_lim_list to struct_ncs_dom_lim category
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
       void _write_struct_ncs_dom_lim_list(Block& block);

       void _get_entity_category_definition(std::vector<std::pair<std::string, std::pair<std::string, std::map<std::string, std::string> > > >&
                                            entity_category);

       void _get_entity_source_category_definition(std::vector<std::pair<std::string, std::pair<std::string, std::map<std::string, std::string> > > >&
                                                   entity_category);

       /**
       **  Read entity information from entity, entity_keywords, entity_name_com,
       **      entity_name_sys, entity_poly, entity_src_gen, entity_src_nat and
       **      pdbx_entity_src_syn categories
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
       void _read_entity_info(Block& block);

       /**
       **  Read information from categories
       **
       **  \param[out]: block   - reference to data block
       **  \param[in]: category - category name
       **  \param[in]: item_id  - entity ID item name
       **  \param[in]: item_mapping - cif_item/key mapping
       **  \param[out]: info_mapping - entity ID vs. entity information mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _read_info_mapping(Block& block, const std::string& category, const std::string& item_id, const std::map<std::string, std::string>& item_mapping,
                               std::map<int, std::vector<std::map<std::string, std::string> > >& info_mapping);

       /**
       **  Read information from categories
       **
       **  \param[out]: block   - reference to data block
       **  \param[in]: category - category name
       **  \param[in]: item_id  - entity ID item name
       **  \param[in]: item_mapping - cif_item/key mapping
       **  \param[out]: info_mapping - entity ID vs. entity information mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _read_info_mapping_1(Block& block, const std::string& category, const std::string& item_id, const std::map<std::string, std::string>& item_mapping,
                                 std::map<int, std::vector<std::map<std::string, std::string> > >& info_mapping, const bool& source_flag = true);

       /**
       **  Write entity information to entity, entity_keywords, entity_name_com,
       **      entity_name_sys, entity_src_gen, entity_src_nat and pdbx_entity_src_syn
       **      categories
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
       void _write_entity_info(Block& block);

       /**
       **  Write entity information to entity, entity_keywords, entity_name_com and
       **      entity_name_sys categories
       **
       **  \param[out]: block       - reference to data block
       **  \param[in]: category     - category name
       **  \param[in]: id           - entity ID item name
       **  \param[in]: item_mapping - cif_item/key mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _write_entity_category(Block& block, const std::string& category, const std::string& id, const std::map<std::string, std::string>& item_mapping);

       /**
       **  Write pdbx_entity_nonpoly & pdbx_entity_branch categories
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
       void _write_entity_nonpoly_and_branch_categores(Block& block);

       /**
       **  Write source information to entity_src_nat, entity_src_gen & pdbx_entity_src_syn
       **      categories
       **
       **  \param[out]: block   - reference to data block
       **  \param[in]: category - category name
       **  \param[in]: type     - source type
       **  \param[in]: item_mapping - cif_item/key mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _write_entity_source_category(Block& block, const std::string& category, const std::string& type, const std::map<std::string, std::string>&
                                          item_mapping);

       /**
       **  Write pdbx_entity_branch_descriptor category
       **
       **  \param[out]: block   - reference to data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _write_pdbx_entity_branch_descriptor(Block& block);

       /**
       **  Write pdbx_entity_branch_link category
       **
       **  \param[out]: block   - reference to data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _write_pdbx_entity_branch_link(Block& block);

       void _assign_new_entity_ids(std::map<std::string, std::string>& new_chain_entity_id_mapping, std::map<std::string, std::string>&
                                   old_new_entity_mapping);

       /**
       **  Write entity_poly_seq, pdbx_poly_seq_scheme & pdbx_nonpoly_scheme categories
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
       void _write_entity_scheme_info(Block& block);

       /**
       **  Write atom_type & chem_comp categories
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
       void _write_atom_type_and_chem_comp(Block& block);

       /**
       **  Read feature key/value pairs from CCD definition
       **
       **  \param[in]: comp_id - comp ID
       **  \param[in]: checking_mapping - Checking key/value(s) pairs
       **  \param[in]: t - pointer to feature category table object
       **  \param[out]: branch_list - list of feature key/value pairs
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_chem_comp_feature_list(const std::string& comp_id, const std::map<std::string, std::set<std::string> >& checking_mapping,
                                           ISTable* t, std::vector<std::map<std::string, std::string> >& branch_list);


       /**
       **  Write chem_comp features to model file
       **
       **  \param[out]: block - reference to data block
       **  \param[in]: categoryName - category name
       **  \param[in]: branch_list - list of feature key/value pairs
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _write_chem_comp_feature_categories(Block& block, const std::string& categoryName, const std::vector<std::map<std::string, std::string> >& branch_list);

       /**
       **  Read information from struct_asym category
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
       void _read_struct_asym(Block& block);

       /**
       **  Write struct_asym category
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
       void _write_struct_asym(Block& block);

       /**
       **  Write entity_poly category
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
       void _write_entity_poly(Block& block);

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
       // void _get_oneletter_code(std::string& code, std::string& code_can);

       /**
       **  Update record row in pdbx_poly_seq_scheme and pdbx_nonpoly_scheme tables
       **
       **  \param[out]: t - pointer to scheme table
       **  \param[out]: key_val_mapping - key/value pair
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_scheme_table(ISTable *t, const std::map<std::string, std::string>& key_val_mapping);

       /**
       **  Update atom_sites category
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
       void _update_atom_sites(Block& block);

       /**
       **  Update pdbx_database_status.pdb_format_compatible flag
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
       void _update_pdb_format_compatible(Block& block);

       /**
       **  Update SG related catetories
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
       void _process_sg_entry(Block& block);

       /**
       **  Update KEYWRDS
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
       void _update_struct_keywords(Block& block);

       /**
       **  Update struct category
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
       void _update_struct_category(Block& block);

       /**
       **  Update pdbx_struct_conn_angle category
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
       void _update_pdbx_struct_conn_angle(Block& block);

       /**
       **  Update entity information for spliting polymer
       **
       **  \param[in]: entity_id - old entity id
       **  \param[in]: entityids - new entity ids
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _split_entity_info(const int& entity_id, const std::map<int, std::pair<std::vector<std::string>, std::vector<std::string> > >& entityids);

       /**
       **  Update entity information for merging polymer
       **
       **  \param[in]: entity_id - new entity id
       **  \param[in]: pdb_chainid - PDB chain id
       **  \param[in]: seqs - sequence record
       **  \param[in]: entityids - old entity ids
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _merge_entity_info(const int& entity_id, const std::string& pdb_chainid, const std::vector<std::string>& seqs,
                               const std::vector<std::vector<int> >& entityids);

       /**
       **  Update various cif categories in compliance with dictionary
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_dictionary_compliance(Block& block);

       /**
       **  Check if a table is empty without any dependent child categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: catName - category name
       **  \param[out]: t - pointer to ISTable
       **
       **  \return true if it is empty without any dependent child categories
       **          false otherwise
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _is_empty_table_with_dictionary_check(Block& block, const std::string& catName, ISTable *t);

       /**
       **  Update related parent categories item values
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: catName - child category name
       **  \param[out]: child - pointer to child category ISTable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_parent_categories(Block& block, const std::string& catName, ISTable* child);

       /**
       **  Update mandatory item[s]
       **
       **  \param[out]: t - pointer to ISTable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_mandatory_item(ISTable *t);

       /**
       **  Update enumeration values
       **
       **  \param[out]: t - pointer to ISTable
       **  \param[in]: itemName - table column name
       **  \param[in]: enums - enumeration values
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_enumeration_value(ISTable *t, const std::string& itemName, const std::map<std::string, std::string>& enums);

       /**
       **  Get default mandatory item value.
       **
       **  \param[in]: itemName - table column name
       **  \param[in]: keyFlag - kye item flag
       **  \param[in]: row - table row index
       **
       **  \return default mandatory item value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _cif_update_get_default_mandatory_item_value(const std::string& itemName, const bool& keyFlag, const int& row);

       /**
       **  Update non-valid values to "."
       **
       **  \param[out]: t - pointer to ISTable
       **  \param[in]: itemName - table column name
       **  \param[in]: regular_expression - regular expression string
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_check_value(ISTable *t, const std::string& itemName, const std::string& regular_expression);

       /**
       **  Update selected experimental data categories
       **
       **  \param[out]: t - pointer to ISTable
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_category(ISTable *t);

       /**
       **  Remove un-related empty experimental data categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_remove_exp_spcific_categories(Block& block);

       /**
       **  Update software category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_software(Block& block);

       /**
       **  Update specific categories with duplicate rows
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_repeat_row_categories(Block& block);

       /**
       **  Update diffrn_source & diffrn_radiation_wavelength categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_diffrn_source(Block& block);

       /**
       **  Update diffrn_radiation & diffrn_radiation_wavelength categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_diffrn_radiation(Block& block);

       /**
       **  Update diffrn_radiation & diffrn_radiation_wavelength categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_diffrn_categories(Block& block);

       /**
       **  Update target_item column in target_t table based on the values
       **     in source_item column in source_t table
       **
       **  \param[in]: source_t - pointer to source table 
       **  \param[in]: source_item - source column name
       **  \param[out]: target_t - pointer to target table 
       **  \param[in]: target_item - target column name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_table_index(ISTable *source_t, const std::string& source_item, ISTable *target_t, const std::string& target_item);

       /**
       **  Update column in table using values defined data_set
       **
       **  \param[in]: data_set - data value set
       **  \param[out]: t - pointer to table 
       **  \param[in]: item - column name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_table_index(const std::set<std::string>& data_set, ISTable *t, const std::string& item);

       /**
       **  Remove 'K' from temperature value
       **
       **  \param[out]: t - pointer to table 
       **  \param[in]: item - column name
       **  \param[in]: row - row index
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _remove_k(ISTable *t, const std::string &item, const int& row);

       /**
       **  Remove 'INFINI' from value
       **
       **  \param[out]: t - pointer to table 
       **  \param[in]: item - column name
       **  \param[in]: row - row index
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _remove_infinite(ISTable *t, const std::string &item, const int& row);

       /**
       **  Remove comma "," from value
       **
       **  \param[out]: t - pointer to table 
       **  \param[in]: item - column name
       **  \param[in]: row - row index
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _remove_comma(ISTable *t, const std::string &item, const int& row);

       /**
       **  Update pdbx_nmr_ensemble category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_nmr_ensemble(Block& block);

       /**
       **  Update pdbx_nmr_spectrometer category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_nmr_spectrometer(Block& block);

       /**
       **  Update pdbx_nmr_exptl_sample_conditions category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_nmr_exptl_sample_conditions(Block& block);

       /**
       **  Update pdbx_nmr_exptl category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_nmr_exptl(Block& block);

       /**
       **  Read old_new_version_mapping.cif
       **
       **  \param: Not applicable
       **
       **  \return old vs. new revisio type mapping
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::map<std::string, std::string> _read_old_new_version_mapping();

       /**
       **  Write pdbx_audit_revision* categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: audit_revision_records - new audit revision records
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _write_audit_revision_records(Block& block, const std::map<int, _Audit_Revision>& audit_revision_records);

       /**
       **  Convert X-RAY related categories from V4 dictionary to V5 dictionary
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_xray_V4_V5_conversion(Block& block);

       /**
       **  Convert EM categories from V4 dictionary to V5 dictionary
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_em_V4_V5_conversion(Block& block);

       /**
       **  Convert NMR categories from V4 dictionary to V5 dictionary
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_nmr_V4_V5_conversion(Block& block);

       /**
       **  Update all EM related categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_em_categories(Block& block);

       /**
       **  Update refine_ls_restr_ncs category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_refine_ls_restr_ncs(Block& block);

       /**
       **  Update struct_ncs_dom category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_struct_ncs_dom(Block& block);

       /**
       **  Update struct_ncs_ens category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_struct_ncs_ens(Block& block);

       /**
       **  Update refine category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_refine(Block& block);

       /**
       **  Update various refine-related categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_refine_related_categories(Block& block);

       /**
       **  Update refine_B_iso category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_refine_B_iso(Block& block);

       /**
       **  Update refine_occupancy category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_refine_occupancy(Block& block);

       /**
       **  Update refine_ls_shell category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_refine_ls_shell(Block& block);

       /**
       **  Update refine, reflns, reflns_shell and diffrn categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_refine_diffrn_categories(Block& block);

       /**
       **  Update exptl_crystal_grow category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_exptl_crystal_grow(Block& block);

       /**
       **  Update exptl_crystal_grow_comp category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_exptl_crystal_grow_comp(Block& block);

       /**
       **  Update exptl_crystal category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_exptl_crystal(Block& block);

       /**
       **  Update exptl_crystal category with Matthew coefficient and solvent content
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_matthew_and_solvent(Block& block);

       /**
       **  Update phasing_MAD_set, phasing_MAD_clust, phasing_MAD_clust,
       **     phasing_MAD_expt and phasing_set categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_phasing_MAD_set(Block& block);

       /**
       **  Update pdbx_phasing_MAD_set_shell and pdbx_phasing_MAD_shell categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_phasing_MAD_shell(Block& block);

       /**
       **  Update phasing_MIR_der category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_phasing_MIR_der(Block& block);

       /**
       **  Update phasing_set category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_phasing_set(Block& block);

       /**
       **  Update pdbx_coordinate_model category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_coordinate_model_type(Block& block);

       /**
       **  Update OSAKA -> PDBJ & EBI -> PDBE for process and deposit sites
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_deposit_process_sites(Block& block);

       /**
       **  Update pdbx_prerelease_seq category
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_prerelease_seq(Block& block);

       /**
       **  Update miscellaneous categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_miscellaneous_categories(Block& block);

       /**
       **  Remove empty row(s) & table(s), update ordinal & pdbx_ordinal values 
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_update_remove_empty_row_and_tables(Block& block);

       /**
       **  Read SITE records
       **
       **  \param[in]:  block - reference to a coordinate data block
       **  \param[out]: sites - SITE records
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _read_struct_sites(Block& block, std::list<_SITE>& sites);

       /**
       **  Write SITE records
       **
       **  \param[out]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _write_struct_sites(Block& block);

       /**
       **  Remove software generated SITE records
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
       void _remove_software_generated_site_record();

/*
       void _assignPDBChainIDToNonPolymer(RCSB::Molecule* mol, std::vector<RCSB::Chain*>& polymer_chains,
                                          std::vector<RCSB::Chain*>& nonpolymer_chains);
       void _assignSinglePDBChainID(const std::string& pdb_chain_id, RCSB::Molecule* mol);
       void _assignMultiplePDBChainIDs(xtal& mycrys, std::vector<RCSB::Chain*>& nonpolymer_chains);
       bool _is_numbering_unique(RCSB::Molecule* mol);
       void _assignPDBNumering(RCSB::Molecule* mol);
       void _checking_polymer_numbering(RCSB::Molecule* mol);
       void _checking_nonpolymer_numbering(RCSB::Molecule* mol);
*/

       void _update_linkage_information(const std::set<long>& atom_set, std::list<_LINK>& links);

       void _cif_update_pdbx_data_processing_status(Block& block);

       void _cif_update_add_value_to_pdbx_data_processing_status(Block& block, const std::set<std::string>& value_set);

       void _update_struct_ncs_oper(Block& block);

       void _get_branch_entity_links(RCSB::Chain* chain, const std::list<_LINK>& input_links, const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& bad_links);

       void _insert_link_maps(const std::string& bondtype, const int& first_index, RCSB::Atom* first_atom, const int& second_index, RCSB::Atom* second_atom,
                  std::multimap<int, std::multimap<int, std::map<std::string, std::string> > >& link_maps, std::set<std::string>& link_indices);

       std::string _get_leaving_atom(const std::string& compId, const std::string& atomId);

       std::string _get_stereo_config(const std::string& compId, const std::string& atomId);

       void _annotate_branch_polymers();

       void _add_default_assembly_info();

       void _get_glyco_site_linkage_info();

       RCSB::Residue* _get_glyco_site_residue(RCSB::Molecule* mol, RCSB::Chain* chain);

       void _get_carb_annotation(RCSB::Residue* glycoRes, RCSB::Chain* chain, std::list<_LINK>& model_related_links, int& sequential_count);
 
       void _run_assign_carbohydrate_chain_id();

       void _get_branch_entity_instance_categories(RCSB::Chain* chain, std::vector<ISTable*>& tableList);

       void _get_pdb2glycan_annotation(const std::string& fullPath, const std::string& uid, const std::vector<ISTable*>& tableList,
                                       std::vector<std::map<std::string, std::string> >& linear_descriptors);

       void _write_pdb_format_file(const std::string& fullPath, const std::string& uid, const std::list<_LINK>& model_related_links, const
                                   std::map<std::string, std::string>& pdb_chnid_map, const std::vector<std::pair<std::vector<RCSB::Residue*>,
                                   std::pair<std::string, std::string> > > sugar_residue_list);

       void _get_pdb_care_annotation(const std::string& fullPath, const std::string& uid, const std::string& removeName,
                                     std::vector<std::map<std::string, std::string> >& linear_descriptors);

       void _fill_in_atom_list(RCSB::Residue* res, const std::string& type, const std::string& pdb_chnid,
                               std::list<std::pair<std::string, RCSB::Atom*> >& atom_list,
                               std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& connect_list);

};

#endif
