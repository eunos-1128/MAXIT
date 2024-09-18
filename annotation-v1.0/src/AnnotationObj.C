/*
FILE:     AnnotationObj.C
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
#include <math.h> 

#include "AnnotationObj.h"
#include "CategoryMapping.h"
#include "CifCoordWrite.h"
#include "Element.h"
#include "GeometryUtil.h"
#include "GetPairList.h"
#include "GetPairList.C"
#include "PdbWrite.h"
#include "SeqCodeUtil.h"
#include "UpdateAtomInfo.h"
#include "utillib.h"

static bool get_atom_pair(const _LINK& link, CrySymmetry& cell, RCSB::Atom& metal_atom, RCSB::Atom& polar_atom,
                          std::string& SymOP_metal, std::string& SymOP_polar);

AnnotationObj::AnnotationObj(): ValidateObj()
{
       AnnotationObj::clear();
}

AnnotationObj::~AnnotationObj()
{
       AnnotationObj::clear();
}

void AnnotationObj::clear()
{
       _clear_annotation_obj_allocated_memories();
       _branch_annotated_flag = false;
       _pdb2glycanProgram.clear();
}

void AnnotationObj::_clear_annotation_obj_allocated_memories()
{
       _merge_seq_module_flag = false;
       _coordinate_merge_flag = false;
       _get_skip_option_flag = false;
       _update_xray_v4_to_v5_flag = false;
       _update_em_v4_to_v5_flag = false;
       _update_nmr_v4_to_v5_flag = false;
       _carbohydrate_gmml_track_mode_flag = false;
       _numbering_based_order_flag = false;
       _include_date_original_flag = false;
       _intra_molecular_connectivity_flag = false;
       _classifications.clear();
       _basebase_params.clear();
       _interbase_params.clear();
       _modres.clear();
       _helix.clear();
       _turn.clear();
       _sheet.clear();
       _site.clear();
       _dbrefs.clear();
       _seqadvs.clear();
       _moving_water.clear();
       _glyco_site_linkage_info.clear();
       _nonpolymer_group_mapping.clear();
       _new_chain_entity_id_mapping.clear();
       _prd_features.clear();
       _prd_instances.clear();
       _mandatory_updated_categories.clear();
       _default_exp_type.clear();
       _experimental_method_count = 0;
       _default_diffrn_id.clear();
       _scattering_type_mapping.clear();
       _pdbx_entity_instance_feature.clear();
       _pdbx_remediation_atom_site_mapping.clear();
       _branch_residue_set.clear();
       _glyco_link_set.clear();
       _refine_ls_restr_ncs_list.clear();
       _struct_ncs_dom_lim_list.clear();
}

void AnnotationObj::set_get_skip_option_flag()
{
       _get_skip_option_flag = true;
}

void AnnotationObj::set_update_xray_v4_to_v5_flag()
{
       _update_xray_v4_to_v5_flag = true;
}

void AnnotationObj::set_update_em_v4_to_v5_flag()
{
       _update_em_v4_to_v5_flag = true;
}

void AnnotationObj::set_update_nmr_v4_to_v5_flag()
{
       _update_nmr_v4_to_v5_flag = true;
}

void AnnotationObj::set_carbohydrate_gmml_track_mode_flag()
{
       _carbohydrate_gmml_track_mode_flag = true;
}

void AnnotationObj::set_include_date_original_flag()
{
       _include_date_original_flag = true;
}

void AnnotationObj::set_intra_molecular_connectivity_flag()
{
       _intra_molecular_connectivity_flag = true;
}

void AnnotationObj::set_pdb2glycan_program(const std::string& program)
{
       _pdb2glycanProgram = program;
}

void AnnotationObj::cif_update()
{
       _default_exp_type = _get_diffraction_experiment_method();
       if (_default_exp_type.empty()) _default_exp_type = ".";

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       // check carefully before turning the following functionality on
       _cif_update_remove_exp_spcific_categories(_cifblock);
       _cif_update_software(_cifblock);
       _cif_update_repeat_row_categories(_cifblock);

       if ((_experiment_type & EXPERIMENT_TYPE_NMR) || (_experiment_type & EXPERIMENT_TYPE_NMR_SOLID)) {
            if (_update_nmr_v4_to_v5_flag) _cif_update_nmr_V4_V5_conversion(_cifblock);
            _cif_update_nmr_ensemble(_cifblock);
            // _cif_update_nmr_spectrometer(_cifblock);
            // _cif_update_nmr_exptl_sample_conditions(_cifblock);
            // _cif_update_nmr_exptl(_cifblock);
       }
       if ((_experiment_type & EXPERIMENT_TYPE_EM) || (_experiment_type & EXPERIMENT_TYPE_EC) || (_experiment_type & EXPERIMENT_TYPE_ET)) {
            if (_update_em_v4_to_v5_flag) _cif_update_em_V4_V5_conversion(_cifblock);
            // _cif_update_em_imaging(_cifblock);
            // _cif_update_em_sample_preparation(_cifblock);
            // _cif_update_em_3d_fitting(_cifblock);
            _cif_update_em_categories(_cifblock);
       }
       if ((_experiment_type & EXPERIMENT_TYPE_XRAY) || (_experiment_type & EXPERIMENT_TYPE_NEUTRON) ||
          (_experiment_type & EXPERIMENT_TYPE_FIBER) || (_experiment_type & EXPERIMENT_TYPE_ELECTRON) ||
          (_experiment_type & EXPERIMENT_TYPE_EC)) {
            if (_update_xray_v4_to_v5_flag) _cif_update_xray_V4_V5_conversion(_cifblock);
            _cif_update_diffrn_source(_cifblock);
            _cif_update_diffrn_radiation(_cifblock);
            // _cif_update_diffrn_categories(_cifblock);
            if (_refine_ls_restr_ncs_list.empty()) _cif_update_refine_ls_restr_ncs(_cifblock);
            // _cif_update_struct_ncs_dom(_cifblock);
            // _cif_update_struct_ncs_ens(_cifblock);
            _cif_update_refine(_cifblock);
            // _cif_update_refine_related_categories(_cifblock);
            // _cif_update_refine_B_iso(_cifblock);
            // _cif_update_refine_occupancy(_cifblock);
            // _cif_update_refine_ls_shell(_cifblock);
            // _cif_update_refine_diffrn_categories(_cifblock);
            // _cif_update_exptl_crystal_grow(_cifblock);
            // _cif_update_exptl_crystal_grow_comp(_cifblock);
            // _cif_update_exptl_crystal(_cifblock);
            // _cif_update_phasing_MAD_set(_cifblock);
            // _cif_update_phasing_MAD_shell(_cifblock);
            // _cif_update_phasing_MIR_der(_cifblock);
            // _cif_update_phasing_set(_cifblock);
       }
       _cif_update_coordinate_model_type(_cifblock);
       _cif_update_deposit_process_sites(_cifblock);
       _cif_update_prerelease_seq(_cifblock);
       _cif_update_miscellaneous_categories(_cifblock);
       _cif_update_dictionary_compliance(_cifblock);
       _cif_update_remove_empty_row_and_tables(_cifblock);

       if (_get_skip_option_flag) _cif_update_pdbx_data_processing_status(_cifblock); 
}

void AnnotationObj::merge_ebi_status(const std::string& status_file)
{
       CifFile *fobj = get_fobj(_logIo, "", status_file);
       if (!fobj) return;

       Block &block = fobj->GetBlock(fobj->GetFirstBlockName());
       ISTable *old_t = getTablePtr(block, "ndb_database_status");
       if (!old_t) {
            delete fobj;
            return;
       }

       std::string cs;

       bool coordinates_hold_for_1_year = false;
       get_value_clean_upper(cs, old_t, 0, "coordinates_hold_for_1_year");
       if (cs == "Y") coordinates_hold_for_1_year = true;

       bool coordinates_hold_for_publication = false;
       get_value_clean_upper(cs, old_t, 0, "coordinates_hold_for_publication");
       if (cs == "Y") coordinates_hold_for_publication = true;

       bool exp_data_hold_for_1_year = false;
       get_value_clean_upper(cs, old_t, 0, "exp_data_hold_for_1_year");
       if (cs == "Y") exp_data_hold_for_1_year = true;

       bool exp_data_hold_for_4_years = false;
       get_value_clean_upper(cs, old_t, 0, "exp_data_hold_for_4_years");
       if (cs == "Y") exp_data_hold_for_4_years = true;

       bool exp_data_hold_for_publication = false;
       get_value_clean_upper(cs, old_t, 0, "exp_data_hold_for_publication");
       if (cs == "Y") exp_data_hold_for_publication = true;

       bool is_nmr = false;
       bool is_xray = false;
       ISTable *t = getTablePtr(block, "exptl");
       if (t) {
            get_value_clean_upper(cs, t, 0, "method");
            if (cs.find("NMR") != string::npos) is_nmr = true;
            else if (cs == "X-RAY DIFFRACTION" || cs == "X-RAY POWDER DIFFRACTION" || cs == "NEUTRON DIFFRACTION" ||
                     cs == "FIBER DIFFRACTION" || cs == "FIBRE DIFFRACTION") is_xray = true;
       }

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       ISTable *new_t = _getTablePtr(_cifblock, "pdbx_database_status");
       if (!new_t) {
            new_t = _newTablePtr("pdbx_database_status");
            new_t->AddRow();
       }

       cs.clear();
       if (coordinates_hold_for_1_year)
            cs = "HOLD FOR 1 YEAR";
       else if (coordinates_hold_for_publication)
            cs = "HOLD FOR PUBLICATION";
       if (!cs.empty()) {
            new_t->UpdateCell(0, "dep_release_code_coordinates", cs);
       }

       cs.clear();
       if (exp_data_hold_for_1_year)
            cs = "HOLD FOR 1 YEAR";
       else if (exp_data_hold_for_4_years)
            cs = "HOLD FOR 4 YEARS";
       else if (exp_data_hold_for_publication)
            cs = "HOLD FOR PUBLICATION";
       if (!cs.empty()) {
            if (is_xray) new_t->UpdateCell(0, "dep_release_code_struct_fact", cs);
            else if (is_nmr) new_t->UpdateCell(0, "dep_release_code_nmr_constraints", cs);
       }

       const std::vector<std::string>& columnNames = new_t->GetColumnNames();
       for (std::vector<std::string>::const_iterator vpos = columnNames.begin(); vpos != columnNames.end(); ++vpos) {
            get_value_clean(cs, new_t, 0, *vpos);
            if (!cs.empty()) continue;
            get_value_clean(cs, old_t, 0, *vpos);
            if (cs.empty()) continue;
            if (cs == "EBI")
                 new_t->UpdateCell(0, *vpos, "PDBE");
            else new_t->UpdateCell(0, *vpos, cs);
       }

       get_value_clean(cs, new_t, 0, "pdb_format_compatible");
       if (cs.empty()) new_t->UpdateCell(0, "pdb_format_compatible", "Y");

       _cifblock.WriteTable(new_t);

       t = getTableCopy(block, "audit_contact_author");
       if (t) _cifblock.WriteTable(t);

       delete fobj;
}

void AnnotationObj::update_deposit_and_release_dates()
{
       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       ISTable *t = _getTablePtr(_cifblock, "pdbx_database_status");
       ISTable *t1 = getTablePtr(_cifblock, "database_PDB_rev");
       if (!t || !t1) return;

       std::string deposit_date, release_date;
       get_value_clean(deposit_date, t, 0, "recvd_initial_deposition_date");
       get_value_clean(release_date, t, 0, "date_of_NDB_release");
       if (!deposit_date.empty() && !release_date.empty()) return;

       std::string found_deposit_date, found_release_date;
       get_value_clean(found_deposit_date, t1, 0, "date_original");
       get_value_clean(found_release_date, t1, 0, "date");

       if (deposit_date.empty() && !found_deposit_date.empty()) t->UpdateCell(0, "recvd_initial_deposition_date", found_deposit_date);
       if (release_date.empty() && !found_release_date.empty()) t->UpdateCell(0, "date_of_NDB_release", found_release_date);
       _cifblock.WriteTable(t);
}

void AnnotationObj::add_depid(const std::string& depid)
{
       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       ISTable *t = getTablePtr(_cifblock, "database_2");
       int row = t->GetNumRows();
       t->AddRow();
       t->UpdateCell(row, "database_id", "WWPDB");
       t->UpdateCell(row, "database_code", depid);
       _cifblock.WriteTable(t);

       t = getTablePtr(_cifblock, "pdbx_database_status");
       if (t) {
            std::string cs;
            get_value_clean_upper(cs, t, 0, "deposit_site");
            if (cs == "PDB") {
                 t->UpdateCell(0, "deposit_site", "BNL");
                 _cifblock.WriteTable(t);
            }
            get_value_clean_upper(cs, t, 0, "process_site");
            if (cs == "PDB") {
                 t->UpdateCell(0, "process_site", "BNL");
                 _cifblock.WriteTable(t);
            }
       }

       _resetFirstBlockName = depid;
}

void AnnotationObj::ReadStructuralFeatures(const bool& skip_flag)
{
       if (!_CifObj) return;
       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);

       // Read entity, entity_keywords, entity_name_com, entity_name_sys, entity_src_nat, entity_src_gen & pdbx_entity_src_syn
       _read_entity_info(_cifblock);

       // Read _struct_asym.pdbx_blank_PDB_chainid_flag & _struct_asym.details
       _read_struct_asym(_cifblock);

       // Read MODRES
       _read_modres(_cifblock);

       // Read HELIX & TURN
       if (!skip_flag) _read_struct_conf(_cifblock);

       // Read SHEET
       if (!skip_flag) _read_struct_sheet(_cifblock);

       // Read LINK, SSBOND, SLTBRG & BSPAIR
       if (!skip_flag) _read_pdbx_struct_link_and_struct_conn(_cifblock);
       _update_leaving_atom_flag();

       // Read SITE
       if (!skip_flag) _read_struct_sites(_cifblock, _site);

       // Read _basebase_params
       if (!skip_flag) _read_ndb_struct_na_base_pair(_cifblock);

       // Read _interbase_params
       if (!skip_flag) _read_ndb_struct_na_base_pair_step(_cifblock);

       // Read assembly info
       _read_struct_assembly_categories(_cifblock);
       _read_pdbx_struct_oper_list(_cifblock);
 
       // Read pdbx_molecule_features
       if (!skip_flag) _read_prd_features(_cifblock);

       // Read pdbx_molecule
       if (!skip_flag) _read_prd_instances(_cifblock);

       // Read pdbx_entity_instance_feature
       if (!skip_flag) _read_pdbx_entity_instance_feature(_cifblock);

       // Read pdbx_remediation_atom_site_mapping
       _read_pdbx_remediation_atom_site_mapping(_cifblock);

       // Read struct_ref_seq
       _read_struct_ref_seq(_cifblock);

       // Read struct_ref_seq_dif
       _read_struct_ref_seq_dif(_cifblock);

       // Read pdbx_solvent_atom_site_mapping
       if (!skip_flag) _read_pdbx_solvent_atom_site_mapping(_cifblock);

       _read_refine_ls_restr_ncs_list(_cifblock);

       _read_struct_ncs_dom_lim_list(_cifblock);
}

void AnnotationObj::WriteStructuralFeatures(const bool& update_z_value_flag, const bool& skip_flag, const bool& update_struct_flag)
{
       if (!_CifObj) return;

       for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            (*mpos)->CheckUniquePDBNumbering();
       }

       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);

       _update_pdb_format_compatible(_cifblock);

       _process_sg_entry(_cifblock);

       _update_struct_keywords(_cifblock);

       if (!_DepUI_Flag) _cif_update_matthew_and_solvent(_cifblock);

       if (update_z_value_flag) _update_Z_value(_cifblock);

       _update_atom_sites(_cifblock);

       if (update_struct_flag) _update_struct_category(_cifblock);

       // Write entity, entity_keywords, entity_name_com, entity_name_sys, entity_src_nat, entity_src_gen & pdbx_entity_src_syn
       _write_entity_info(_cifblock);

       // Write entity_poly, entity_poly_seq, pdbx_poly_seq_scheme & pdbx_nonpoly_scheme
       // Must run before _write_atom_type_and_chem_comp() function in order to create pdbx_chem_comp_identifier & pdbx_chem_comp_synonyms
       // categories for sugar residues
       _write_entity_scheme_info(_cifblock);

       // Write atom_type & chem_comp
       _write_atom_type_and_chem_comp(_cifblock);

       _write_struct_asym(_cifblock);

       // Write pdbx_struct_mod_residue
       _write_modres(_cifblock);

       // Write struct_conf_type & struct_conf
       _update_struct_conf(_cifblock);

       // Write struct_sheet, struct_sheet_order, struct_sheet_range & pdbx_struct_sheet_hbond
       _update_sheet_categories(_cifblock);

       // Delete pdbx_struct_link
       deleteTable(_cifblock, "pdbx_struct_link");
 
       // Write struct_conn & struct_conn_type
       _update_struct_conn(_cifblock);

       // Write struct_site & struct_site_gen
       _write_struct_sites(_cifblock);

       // Write ndb_struct_na_base_pair
       _update_ndb_struct_na_base_pair(_cifblock);

       // Write ndb_struct_na_base_pair_step
       _update_ndb_struct_na_base_pair_step(_cifblock);

       // Write assembly info
       _write_struct_assembly_categories(_cifblock);
       _write_pdbx_struct_oper_list(_cifblock);
 
       // Write pdbx_molecule_features
       _write_prd_features(_cifblock);

       // Write pdbx_molecule
       _write_prd_instances(_cifblock);

       // Write pdbx_entity_instance_feature
       _write_pdbx_entity_instance_feature(_cifblock);

       // Write struct_ref_seq
       if (!skip_flag) _write_struct_ref_seq(_cifblock);

       if (!skip_flag) _updateStructRef(_new_chain_entity_id_mapping);

       // Write struct_ref_seq_dif
       if (!skip_flag) _write_struct_ref_seq_dif(_cifblock);

       _get_missing_or_zero_occupancy_residues_or_atoms();

       // Write validation related categories
       _updateValidationCategories(_cifblock);

       // write pdbx_struct_conn_angle for REMARK 620
       _update_pdbx_struct_conn_angle(_cifblock);

       // _cif_update_prerelease_seq(_cifblock);
       // _cif_update_miscellaneous_categories(_cifblock);
       // _cif_update_remove_empty_row_and_tables(_cifblock);
       cif_update();

       // Update _struct_ncs_oper.id to integer number
       _update_struct_ncs_oper(_cifblock);

       _write_refine_ls_restr_ncs_list(_cifblock);

       _write_struct_ncs_dom_lim_list(_cifblock);
}

void AnnotationObj::Update_NMR_sample_details()
{
       if (!_CifObj) return;
       Block& block = _CifObj->GetBlock(_firstBlockName);

       ISTable *t1 = _getTablePtr(block, "pdbx_nmr_sample_details");
       if (!t1) return;

       ISTable *t2 = _getTablePtr(block, "pdbx_nmr_exptl_sample");
       if (!t2) return;

       int rowNo1 = t1->GetNumRows();
       int rowNo2 = t2->GetNumRows();

       std::string cs, cs1, contents, component, concentration, concentration_units, isotopic_labeling;
       for (int i = 0; i < rowNo1; i++) {
            contents.clear();
            get_value_clean(cs, t1, i, "solution_id");
            for (int j = 0; j < rowNo2; j++) {
                 get_value_clean(cs1, t2, j, "solution_id");
                 if (cs1 == cs) {
                      get_value_clean(component, t2, j, "component");
                      get_value_clean(concentration, t2, j, "concentration");
                      if (concentration.empty()) get_value_clean(concentration, t2, j, "concentration_range");
                      get_value_clean(concentration_units, t2, j, "concentration_units");
                      get_value_clean(isotopic_labeling, t2, j, "isotopic_labeling");
                      if (isotopic_labeling == "none" || isotopic_labeling == "natural abundance") {
                           isotopic_labeling = "";
                           // t2->UpdateCell(j, "isotopic_labeling", isotopic_labeling);
                      }
                      std::string itemName = concentration + " " + concentration_units;
                      if (!isotopic_labeling.empty()) itemName += " " + isotopic_labeling;
                      itemName += " " + component;
                      if (!contents.empty()) contents += ", ";
                      contents += itemName;
                 }
            }
            get_value_clean(cs, t1, i, "solvent_system");
            if (cs != "") {
                 if (contents != "") contents += ", ";
                 contents += cs;
            }
            t1->UpdateCell(i, "contents", contents);
       }
       block.WriteTable(t1);
       block.WriteTable(t2);
}

void AnnotationObj::CorrectAtomName(const bool& UpperCase_Flag)
{
       if (_molecules.empty()) return;

       std::vector<std::map<std::string, std::vector<std::string> > > pro_atom_name_mapping_list, other_atom_name_mapping_list;
       get_proline_n_terminal_hydrogen_mapping_list(pro_atom_name_mapping_list);
       get_other_n_terminal_hydrogen_mapping_list(other_atom_name_mapping_list);

       std::vector<RCSB::Residue*> residues;
       for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            RCSB::Chain* chain = (*mpos)->GetFirstChain();
            while (chain) {
                 bool first_residue = true;
                 for (unsigned int i = 0; i < chain->SeqLen(); ++i) {
                      chain->GetResidueList(i, residues);
                      if (residues.empty()) continue;
                      for (std::vector<RCSB::Residue*>::iterator pos = residues.begin(); pos != residues.end(); ++pos) {
                           if (UpperCase_Flag) (*pos)->AtomNameUpperCase();
                           if (SeqCodeUtil::is_a_standard_residue((*pos)->ResName()) || SeqCodeUtil::is_standard_Daa_residue((*pos)->ResName())) {
                                (*pos)->correction_atom_name((*mpos)->Mol_ID());
                           }
                           if (first_residue) {
                                if (chain->chain_type() == "ATOMP") {
                                     if ((*pos)->ResName() == "PRO")
                                          (*pos)->fix_n_terminal_hydrogen((i == 0), pro_atom_name_mapping_list);
                                     else (*pos)->fix_n_terminal_hydrogen((i == 0), other_atom_name_mapping_list);
                                } else if (chain->chain_type() == "ATOMN") (*pos)->fix_5_terminal_hydrogen();
                           }
                      }
                      first_residue = false;
                 }
                 chain = (*mpos)->GetNextChain();
            }
       }
}

void AnnotationObj::AssignAsymId()
{
       if (_molecules.empty()) return;

       for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            (*mpos)->AssignAsymId();
       }
}

void AnnotationObj::UpdateStructMonProtCis()
{
       if (_molecules.empty()) return;
       if (!_CifObj) return;

       StandardUtil::initialize();
       if (!StandardUtil::Read(*_logIo, _rcsbroot)) return;

       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);

       std::vector<RCSB::Residue*> prev, curr;
       std::vector<std::vector<RCSB::Residue*> > pair_list, residue_lists;
       // value: atomname1_altloc1_atomname2_altloc2
       std::set<std::string> allowed_inter_residue_bonds;

       std::list<Value> cis_peptide_list;
       cis_peptide_list.clear();

       for (std::vector<RCSB::Molecule*>::const_iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            RCSB::Chain* chain = (*mpos)->GetFirstChain();
            while (chain) {
                 if (chain->chain_type() != "ATOMP") {
                      chain = (*mpos)->GetNextChain();
                      continue;
                 }

                 prev.clear();
                 chain->GetFirstResidueList(curr);
                 while (!curr.empty()) {
                      if (!prev.empty()) {
                           residue_lists.clear();
                           residue_lists.push_back(prev);
                           residue_lists.push_back(curr);
                           GetPairList(residue_lists, pair_list);
                           for (std::vector<std::vector<RCSB::Residue*> >::const_iterator pos = pair_list.begin(); pos != pair_list.end(); ++pos) {
                                _checkInterResidueBond((*mpos)->Mol_ID(), (*pos)[0], (*pos)[1], chain->chain_type(), allowed_inter_residue_bonds, true);
                                if (allowed_inter_residue_bonds.empty()) continue;
                                _getCISOMEGATorsion((*mpos)->Mol_ID(), (*pos)[0], (*pos)[1], allowed_inter_residue_bonds, cis_peptide_list);
                           }
                      }

                      prev = curr;
                      chain->GetNextResidueList(curr);
                 }
                 chain = (*mpos)->GetNextChain();
            }
       }

       if (cis_peptide_list.empty()) {
            deleteTable(_cifblock, "struct_mon_prot_cis");
            return;
       }

       ISTable *t = _newTablePtr("struct_mon_prot_cis");

       int i = 0;
       for (std::list<Value>::const_iterator pos = cis_peptide_list.begin(); pos != cis_peptide_list.end(); ++pos) {
            const std::vector<RCSB::Atom*>& atoms = pos->atoms();
            t->AddRow();
            t->UpdateCell(i, "pdbx_id", String::IntToString(i + 1));
            t->UpdateCell(i, "label_comp_id", atoms[0]->restype());
            t->UpdateCell(i, "label_seq_id", atoms[0]->resnum());
            t->UpdateCell(i, "label_asym_id", atoms[0]->chnid());
            t->UpdateCell(i, "label_alt_id", ".");
            t->UpdateCell(i, "pdbx_PDB_ins_code", atoms[0]->ins_code());
            t->UpdateCell(i, "auth_comp_id", atoms[0]->pdb_resnam());
            t->UpdateCell(i, "auth_seq_id", atoms[0]->pdb_resnum());
            t->UpdateCell(i, "auth_asym_id", atoms[0]->pdb_chnid());
            t->UpdateCell(i, "pdbx_label_comp_id_2", atoms[3]->restype());
            t->UpdateCell(i, "pdbx_label_seq_id_2", atoms[3]->resnum());
            t->UpdateCell(i, "pdbx_label_asym_id_2", atoms[3]->chnid());
            t->UpdateCell(i, "pdbx_PDB_ins_code_2", atoms[3]->ins_code());
            t->UpdateCell(i, "pdbx_auth_comp_id_2", atoms[3]->pdb_resnam());
            t->UpdateCell(i, "pdbx_auth_seq_id_2", atoms[3]->pdb_resnum());
            t->UpdateCell(i, "pdbx_auth_asym_id_2", atoms[3]->pdb_chnid());
            t->UpdateCell(i, "pdbx_PDB_model_num", String::IntToString(pos->Mol_ID()));
            t->UpdateCell(i, "pdbx_omega_angle", FloatToString(pos->val(), 0, 2));
            i++;
       }

       _cifblock.WriteTable(t);
}

void AnnotationObj::UpdateCoordinateCategory()
{
       if (_molecules.empty()) return;
       if (!_CifObj) return;

       Block& _cifblock = _CifObj->GetBlock(_firstBlockName);
/*
       ISTable *t = getTablePtr(_cifblock, "pdbx_original_pdb_coordinates");
       if (!t) t = getTablePtr(_cifblock, "ndb_original_pdb_coordinates");
       if (!t) t = getTablePtr(_cifblock, "ndb_original_ndb_coordinates");
       if (!t) t = getTablePtr(_cifblock, "pdb_original_pdb_coordinates");

       if (t) {
            _updateEncapCoordinateas(t);
            _cifblock.WriteTable(t);
       } else */ _updateAtomSite(_cifblock);

       _updateSolventAtomSiteMapping(_cifblock);

       // Write pdbx_remediation_atom_site_mapping
       _write_pdbx_remediation_atom_site_mapping(_cifblock);
}

void AnnotationObj::AddDefaultAssembly()
{
       if (!_assemblies.empty()) return;

       Block &_cifblock = _CifObj->GetBlock(_firstBlockName);

       std::vector<std::string> linear_polymer_chains, all_chains;
       linear_polymer_chains.clear();
       all_chains.clear();

       std::string type, chain_id;
       ISTable *t = getTablePtr(_cifblock, "struct_asym");
       if (t) {
            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean(chain_id, t, i, "pdbx_PDB_id");
                 if (chain_id.empty()) continue;
                 get_value_clean_upper(type, t, i, "pdbx_type");
                 if ((type == "ATOMN") || (type == "ATOMP")) linear_polymer_chains.push_back(chain_id);
                 if ((type == "ATOMN") || (type == "ATOMP") || (type == "ATOMS")) all_chains.push_back(chain_id);
            }
       }

       if (all_chains.empty()) {
            std::vector<std::pair<std::string, std::string> > cat_item_list;
            cat_item_list.clear();
            cat_item_list.push_back(std::make_pair("pdbx_poly_seq_scheme", "pdb_strand_id"));
            cat_item_list.push_back(std::make_pair("pdbx_branch_scheme", "pdb_asym_id"));
 
            std::set<std::string> unique_set;
            unique_set.clear();

            for (std::vector<std::pair<std::string, std::string> >::const_iterator vpos = cat_item_list.begin(); vpos != cat_item_list.end(); ++vpos) {
                 t = getTablePtr(_cifblock, vpos->first);
                 if (!t) continue;

                 for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                      get_value_clean(chain_id, t, i, vpos->second);
                      if (chain_id.empty()) continue;
                      if (unique_set.find(chain_id) != unique_set.end()) continue;
                      unique_set.insert(chain_id);
                      all_chains.push_back(chain_id);
                      if (vpos->first == "pdbx_poly_seq_scheme") linear_polymer_chains.push_back(chain_id);
                 }
            }
       }
       if (all_chains.empty()) return;

       Assembly assembly;
       assembly.clear();
       assembly.setValue("id", "1");
       assembly.setValue("method", "author_defined_assembly");
       if (!linear_polymer_chains.empty())
            assembly.setValue("oligomeric", String::IntToString(linear_polymer_chains.size()));
       else assembly.setValue("oligomeric", String::IntToString(all_chains.size()));
       assembly.InsertChains("1", all_chains);
       assembly.UpdateOligomericDetails();
       _assemblies.push_back(assembly);

       SymmMatrix symmtrix;
       symmtrix.clear();
       symmtrix.setValue("id", "1");
       symmtrix.setValue("type", "identity operation");
       symmtrix.setValue("name", "1_555");
       symmtrix.setValue("symmetry_operation", "x,y,z");
       symmtrix.setValue("matrix[1][1]", "1.0000000000");
       symmtrix.setValue("matrix[1][2]", "0.0000000000");
       symmtrix.setValue("matrix[1][3]", "0.0000000000");
       symmtrix.setValue("vector[1]", "0.0000000000");
       symmtrix.setValue("matrix[2][1]", "0.0000000000");
       symmtrix.setValue("matrix[2][2]", "1.0000000000");
       symmtrix.setValue("matrix[2][3]", "0.0000000000");
       symmtrix.setValue("vector[2]", "0.0000000000");
       symmtrix.setValue("matrix[3][1]", "0.0000000000");
       symmtrix.setValue("matrix[3][2]", "0.0000000000");
       symmtrix.setValue("matrix[3][3]", "1.0000000000");
       symmtrix.setValue("vector[3]", "0.0000000000");

       _SymmMatrix_Mapping.clear();
       _SymmMatrices.clear();
       _SymmMatrix_Mapping.insert(std::make_pair("1", _SymmMatrices.size()));
       _SymmMatrices.push_back(symmtrix);
}

void AnnotationObj::_updateEncapCoordinateas(ISTable *t)
{
       std::string text = "";
       PdbWrite writer;
       writer.setLog(_logIo);
       writer.setCCDic(_ccDic);
       writer.setMolecule(&_molecules);
       writer.WriteCoordinates(text);
       t->UpdateCell(0, "coord_section", text);
}

void AnnotationObj::_updateAtomSite(Block& block)
{
       if (_molecules.empty()) return;

       deleteTable(block, "pdbx_original_pdb_coordinates");
       deleteTable(block, "ndb_original_pdb_coordinates");
       deleteTable(block, "ndb_original_ndb_coordinates");
       deleteTable(block, "pdb_original_pdb_coordinates");

       ISTable *t1 = _newTablePtr("atom_site");
       if (!_extra_atom_site_item_list.empty()) check_missing_item(t1, _extra_atom_site_item_list);
       const std::vector<std::string>& atomitemNames = t1->GetColumnNames();
       ISTable *t2 = _newTablePtr("atom_site_anisotrop");
       const std::vector<std::string>& anisotropitemNames = t2->GetColumnNames();

       std::vector<RCSB::Residue*> residue_list;
       int atomSerialNo = 0;
       bool has_sigatm = false;
       bool has_siguij = false;
       for (std::vector<RCSB::Molecule*>::const_iterator pos = _molecules.begin(); pos != _molecules.end(); ++pos) {
            RCSB::Chain* chain = (*pos)->GetFirstChain();
            while (chain) {
                 chain->GetFirstResidueList(residue_list);
                 while (!residue_list.empty()) {
                      for (std::vector<RCSB::Residue*>::const_iterator rpos = residue_list.begin(); rpos != residue_list.end(); ++rpos) {
                           std::string PDB_atom_token = "HETATM";
                           if ((chain->chain_type() == "ATOMP" || chain->chain_type() == "ATOMN") &&
                               (SeqCodeUtil::is_a_standard_residue((*rpos)->ResName()) || (*rpos)->ResName() == "T"))
                                PDB_atom_token = "ATOM";
                           RCSB::Atom* atom = (*rpos)->GetFirstAtom();
                           while (atom) {
                                atom->setValue(PDB_atom_token, 0);
                                if (!(*rpos)->entity_id().empty())
                                     atom->setValue((*rpos)->entity_id(), 23);
                                else atom->setValue(chain->entity_id(), 23);
                                atomSerialNo++;
                                if (atom->has_sigatm()) has_sigatm = true;
                                if (atom->has_siguij()) has_siguij = true;
                                CifCoordWrite::update_atom_site(atomitemNames, (*pos)->Mol_ID(), atomSerialNo, atom, t1);
                                if (atom->has_anisou()) _update_atom_site_anisotrop(anisotropitemNames, atomSerialNo, atom, t2);
                                atom = (*rpos)->GetNextAtom();
                           }
                      }
                      chain->GetNextResidueList(residue_list);
                 }
                 chain = (*pos)->GetNextChain();
            }
       }

       if (t1->GetNumRows() > 0) {
            if (!has_sigatm) {
                 const char *remove_sigatm[5] = { "Cartn_x_esd", "Cartn_y_esd", "Cartn_z_esd", "occupancy_esd", "B_iso_or_equiv_esd" };
                 for (int i = 0; i < 5; ++i) {
                      if (t1->IsColumnPresent(remove_sigatm[i])) t1->DeleteColumn(remove_sigatm[i]);
                 }
            }
            block.WriteTable(t1);
       } else { delete t1; if (_coordinate_merge_flag) deleteTable(block, "atom_site"); }

       if (t2->GetNumRows() > 0) {
            if (!has_siguij) {
                 const char *remove_siguij[6] = { "U[1][1]_esd", "U[2][2]_esd", "U[3][3]_esd", "U[1][2]_esd", "U[1][3]_esd", "U[2][3]_esd" };
                 for (int i = 0; i < 6; ++i) {
                      if (t2->IsColumnPresent(remove_siguij[i])) t2->DeleteColumn(remove_siguij[i]);
                 }
            }
            block.WriteTable(t2);
       } else { delete t2; if (_coordinate_merge_flag) deleteTable(block, "atom_site_anisotrop"); }
}

void AnnotationObj::_updateSolventAtomSiteMapping(Block& block)
{
       if (_moving_water.empty()) {
            deleteTable(block, "pdbx_solvent_atom_site_mapping");
            return;
       }

       ISTable *t = _newTablePtr("pdbx_solvent_atom_site_mapping");
       const std::vector<std::string>& itemNames = t->GetColumnNames();

       std::map<std::string, std::string> item_mapping;
       item_mapping.clear();
       item_mapping.insert(std::make_pair("auth_alt_id", "label_alt_id"));
       item_mapping.insert(std::make_pair("PDB_ins_code", "pdbx_PDB_ins_code"));
       item_mapping.insert(std::make_pair("pre_auth_asym_id", "pdbx_auth_asym_id"));
       item_mapping.insert(std::make_pair("pre_auth_atom_id", "pdbx_auth_atom_name"));
       item_mapping.insert(std::make_pair("pre_auth_comp_id", "pdbx_auth_comp_id"));
       item_mapping.insert(std::make_pair("pre_auth_seq_id", "pdbx_auth_seq_id"));
       item_mapping.insert(std::make_pair("pre_auth_alt_id", "label_alt_id"));
       item_mapping.insert(std::make_pair("pre_PDB_ins_code", "pdbx_PDB_ins_code"));

       std::string cs;
       int row = 0;
       for (std::list<_MOVING_ATOM>::const_iterator lpos = _moving_water.begin(); lpos != _moving_water.end(); ++lpos) {
            t->AddRow();
            for (std::vector<std::string>::const_iterator pos = itemNames.begin(); pos != itemNames.end(); ++pos) {
                 std::map<std::string, std::string>::const_iterator mpos = item_mapping.find(*pos);
                 if (*pos == "id")
                      t->UpdateCell(row, *pos, String::IntToString(row + 1));
                 else if (*pos == "symmetry") {
                      if (!lpos->symmetry.empty())
                           t->UpdateCell(row, *pos, lpos->symmetry);
                      else t->UpdateCell(row, *pos, String::IntToString(lpos->sym + 1) + "_" + String::IntToString(lpos->lx + 5) +
                                                     String::IntToString(lpos->ly + 5) + String::IntToString(lpos->lz + 5));
                 } else if (*pos == "symmetry_as_xyz") {
                      if (!lpos->symmetry_as_xyz.empty())
                           cs = lpos->symmetry_as_xyz;
                      else {
                           _cell.ndb_get_symmetry_operation_name(lpos->sym, lpos->lx, lpos->ly, lpos->lz, cs);
                           String::UpperCase(cs);
                      }
                      t->UpdateCell(row, *pos, cs);
                 } else if (*pos == "pre_Cartn_x")
                      t->UpdateCell(row, *pos, FloatToString(lpos->old_coord.x, 0, 3));
                 else if (*pos == "pre_Cartn_y")
                      t->UpdateCell(row, *pos, FloatToString(lpos->old_coord.y, 0, 3));
                 else if (*pos == "pre_Cartn_z")
                      t->UpdateCell(row, *pos, FloatToString(lpos->old_coord.z, 0, 3));
                 else {
                      std::string item = *pos;
                      if (mpos != item_mapping.end()) item = mpos->second;
                      cs = lpos->atom->getValue(item);
                      if (cs.empty() && (item == "label_alt_id" || item == "label_seq_id"))
                           cs =  ".";
                      t->UpdateCell(row, *pos, cs);
                 }
            }
            row++;
       }
       block.WriteTable(t);
}

void AnnotationObj::_update_atom_site_anisotrop(const std::vector<std::string>& items, const int& atom_id, RCSB::Atom* atom, ISTable* t)
{
       std::vector<std::string> data;
       atom->getAuxiliaryValue("ANISOU", data);
       if (data.empty()) return;

       int row = t->GetNumRows();
       t->AddRow();
       for (std::vector<std::string>::const_iterator pos = items.begin(); pos != items.end(); ++pos) {
            if (*pos == "id")
                 t->UpdateCell(row, *pos, String::IntToString(atom_id));
            else if (*pos == "pdbx_label_alt_id" || *pos == "pdbx_label_seq_id") {
                 std::string cs = atom->getAnisouValue(*pos);
                 if (cs.empty() || cs == "?") cs =  ".";
                 t->UpdateCell(row, *pos, cs);
            } else t->UpdateCell(row, *pos, atom->getAnisouValue(*pos));
       }
}

void AnnotationObj::_update_Z_value(Block& block)
{
       // always automatically update Z value based on John B.'s suggestion:
       // Z value = (highest number of polymer chains) X (number of symmetry operators) X (number of generate matrices or 1 if none are given)
       //
       if (_cell.is_artifical() || _molecules.empty()) return;

       ISTable *t = _getTablePtr(block, "cell");
       if (!t) return;

       std::string z;
       get_value_clean(z, t, 0, "Z_PDB");
       // if (!z.empty()) return;

       int num_ncs_oper = 0;
       ISTable *ncs = getTablePtr(block, "struct_ncs_oper");
       if (ncs) {
            int rowNo = ncs->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_lower(z, ncs, i, "code");
                 if (z == "generate") num_ncs_oper++;
            }
       }
       num_ncs_oper += 1;

       int z_value = _cell.symops().size() * _molecules[0]->MaxNumEntity() * num_ncs_oper;
       if (z_value) {
            t->UpdateCell(0, "Z_PDB", String::IntToString(z_value));
            block.WriteTable(t);
       }
}

void AnnotationObj::_read_struct_ref_seq(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "struct_ref_seq");
       if (!t) return;

       std::map<std::string, std::string> mapping;
       std::string cs, pdb_num, ins_code;

       _DBREF dbref;

       int rowNo = t->GetNumRows();
       const std::vector<std::string>& items = t->GetColumnNames();
       for (int i = 0; i < rowNo; ++i) {
            mapping.clear();
            for (std::vector<std::string>::const_iterator vpos = items.begin(); vpos != items.end(); ++vpos) {
                 get_value_clean(cs, t, i, *vpos);
                 if (cs.empty()) continue;
                 mapping.insert(std::make_pair(*vpos, cs));
            }
            if (mapping.empty()) continue;

            dbref.mol_index = -1;
            dbref.chain_index = -1;
            dbref.fst_pos_index = -1;
            dbref.snd_pos_index = -1;
            dbref.data_mapping = mapping;

            get_value_clean(cs, t, i, "pdbx_strand_id");
            if (!cs.empty()) {
                 RCSB::Chain* chain = _molecules[0]->GetPolyChain(cs);
                 if (chain) {
                      dbref.mol_index = 0;
                      dbref.chain_index = chain->index();
                      get_value_clean(pdb_num, t, i, "pdbx_auth_seq_align_beg");
                      get_value_clean(ins_code, t, i, "pdbx_seq_align_beg_ins_code");
                      dbref.fst_pos_index = chain->find_field_index(pdb_num, ins_code);
                      get_value_clean(pdb_num, t, i, "pdbx_auth_seq_align_end");
                      get_value_clean(ins_code, t, i, "pdbx_seq_align_end_ins_code");
                      dbref.snd_pos_index = chain->find_field_index(pdb_num, ins_code);
                 }
            }

            _dbrefs.push_back(dbref);
       }

       // if (!_dbrefs.empty()) block.DeleteTable("struct_ref_seq");
}

void AnnotationObj::_write_struct_ref_seq(Block& block)
{
       if (_merge_seq_module_flag || _dbrefs.empty()) return;

       ISTable *t = _newTablePtr("struct_ref_seq");
       int row = 0;
       for (std::list<_DBREF>::const_iterator lpos = _dbrefs.begin(); lpos != _dbrefs.end(); ++lpos) {
            t->AddRow();
            std::string pdbx_PDB_id_code = "";
            std::string pdbx_db_accession = "";
            for (std::map<std::string, std::string>::const_iterator mpos = lpos->data_mapping.begin(); mpos != lpos->data_mapping.end(); ++mpos) {
                 if (mpos->first == "pdbx_PDB_id_code") pdbx_PDB_id_code = mpos->second;
                 if (mpos->first == "pdbx_db_accession") pdbx_db_accession = mpos->second;
                 t->UpdateCell(row, mpos->first, mpos->second);
            }
            bool is_self_reference = false;
            if (!pdbx_PDB_id_code.empty() && (pdbx_PDB_id_code == pdbx_db_accession)) is_self_reference = true;
            if (lpos->mol_index >= 0 && lpos->chain_index >= 0) {
                 bool is_removed = false;
                 RCSB::Chain* chain = _molecules[lpos->mol_index]->GetIndexChain(lpos->chain_index, is_removed);
                 if (chain) {
                      t->UpdateCell(row, "pdbx_strand_id", chain->PDB_ChainID());
                      _FIELD* seqres = chain->SeqRes(lpos->fst_pos_index);
                      if (seqres) {
                           t->UpdateCell(row, "seq_align_beg", seqres->Field[2]);
                           t->UpdateCell(row, "pdbx_auth_seq_align_beg", seqres->Field[4]);
                           if (is_self_reference) t->UpdateCell(row, "db_align_beg", seqres->Field[4]);
                           t->UpdateCell(row, "pdbx_seq_align_beg_ins_code", seqres->InsCode);
                      }
                      seqres = chain->SeqRes(lpos->snd_pos_index);
                      if (seqres) {
                           t->UpdateCell(row, "seq_align_end", seqres->Field[2]);
                           t->UpdateCell(row, "pdbx_auth_seq_align_end", seqres->Field[4]);
                           if (is_self_reference) t->UpdateCell(row, "db_align_end", seqres->Field[4]);
                           t->UpdateCell(row, "pdbx_seq_align_end_ins_code", seqres->InsCode);
                      }
                 }
            }
            row++;
       }
       block.WriteTable(t);
}

void AnnotationObj::_read_struct_ref_seq_dif(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "struct_ref_seq_dif");
       if (!t) return;

       std::map<std::string, std::string> mapping;
       std::string cs, pdb_num, ins_code;

       _SEQADV seqadv;

       int rowNo = t->GetNumRows();
       const std::vector<std::string>& items = t->GetColumnNames();
       for (int i = 0; i < rowNo; ++i) {
            mapping.clear();
            for (std::vector<std::string>::const_iterator vpos = items.begin(); vpos != items.end(); ++vpos) {
                 get_value_clean(cs, t, i, *vpos);
                 if (cs.empty()) continue;
                 mapping.insert(std::make_pair(*vpos, cs));
            }
            if (mapping.empty()) continue;

            seqadv.mol_index = -1;
            seqadv.chain_index = -1;
            seqadv.pos_index = -1;
            seqadv.data_mapping = mapping;

            get_value_clean(cs, t, i, "pdbx_pdb_strand_id");
            if (!cs.empty()) {
                 RCSB::Chain* chain = _molecules[0]->GetPolyChain(cs);
                 if (chain) {
                      seqadv.mol_index = 0;
                      seqadv.chain_index = chain->index();
                      get_value_clean(pdb_num, t, i, "pdbx_auth_seq_num");
                      get_value_clean(ins_code, t, i, "pdbx_pdb_ins_code");
                      seqadv.pos_index = chain->find_field_index(pdb_num, ins_code);
                 }
            }

            _seqadvs.push_back(seqadv);
       }

       // if (!_seqadvs.empty()) block.DeleteTable("struct_ref_seq_dif");
}

void AnnotationObj::_write_struct_ref_seq_dif(Block& block)
{
       if (_merge_seq_module_flag || _seqadvs.empty()) return;

       ISTable *t = _newTablePtr("struct_ref_seq_dif");
       int row = 0;
       for (std::list<_SEQADV>::const_iterator lpos = _seqadvs.begin(); lpos != _seqadvs.end(); ++lpos) {
            t->AddRow();
            for (std::map<std::string, std::string>::const_iterator mpos = lpos->data_mapping.begin(); mpos != lpos->data_mapping.end(); ++mpos) {
                 t->UpdateCell(row, mpos->first, mpos->second);
            }
            if (lpos->mol_index >= 0 && lpos->chain_index >= 0) {
                 bool is_removed = false;
                 RCSB::Chain* chain = _molecules[lpos->mol_index]->GetIndexChain(lpos->chain_index, is_removed);
                 if (chain) {
                      t->UpdateCell(row, "pdbx_pdb_strand_id", chain->PDB_ChainID());
                      _FIELD* seqres = chain->SeqRes(lpos->pos_index);
                      if (seqres) {
                           t->UpdateCell(row, "seq_num", seqres->Field[2]);
                           t->UpdateCell(row, "pdbx_auth_seq_num", seqres->Field[4]);
                           t->UpdateCell(row, "pdbx_pdb_ins_code", seqres->InsCode);
                      }
                 }
            }
            t->UpdateCell(row, "pdbx_ordinal", String::IntToString(row + 1));
            row++;
       }
       block.WriteTable(t);
}

void AnnotationObj::_read_pdbx_solvent_atom_site_mapping(Block& block)
{
       if (_molecules.empty()) return;

       ISTable *t = getTablePtr(block, "pdbx_solvent_atom_site_mapping");
       if (!t) return;

       std::string cs;
       std::vector<std::string> items, data;

       items.clear();
       items.push_back("auth_asym_id");
       items.push_back("auth_comp_id");
       items.push_back("auth_seq_id");
       items.push_back("PDB_ins_code");
       items.push_back("auth_atom_id");
       items.push_back("auth_alt_id");
       items.push_back("symmetry");
       items.push_back("symmetry_as_xyz");
       items.push_back("pre_Cartn_x");
       items.push_back("pre_Cartn_y");
       items.push_back("pre_Cartn_z");

       _MOVING_ATOM change_atom;
       change_atom.sym = 0;
       change_atom.lx  = 0;
       change_atom.ly  = 0;
       change_atom.lz  = 0;

       int rowNo = t->GetNumRows();
       for (int i = 0; i < rowNo; ++i) {
            data.clear();
            for (std::vector<std::string>::const_iterator vpos = items.begin(); vpos != items.end(); ++vpos) {
                 get_value_clean(cs, t, i, *vpos);
                 data.push_back(cs);
            }

            RCSB::Residue* res = _molecules[0]->find_pdb_residue(data[0], data[1], data[2], data[3]);
            if (!res) continue;

            RCSB::Atom* atom = res->find_atom(data[4], data[5], true);
            if (!atom) continue;

            change_atom.atom = atom;
            change_atom.symmetry = data[6];
            change_atom.symmetry_as_xyz = data[7];
            change_atom.old_coord.x = atof(data[8].c_str());
            change_atom.old_coord.y = atof(data[9].c_str());
            change_atom.old_coord.z = atof(data[10].c_str());
            _moving_water.push_back(change_atom);
       }
}

void AnnotationObj::correct_terminal_atoms(const bool& rename_flag)
{
       std::set<long> delete_atom_set;
       delete_atom_set.clear();

       std::vector<RCSB::Residue*> residues;
       std::vector<RCSB::Atom*> atoms, oxt_atoms;
       std::map<unsigned int, RCSB::Residue*> residue_mapping;
       for (std::vector<RCSB::Molecule*>::iterator mpos = _molecules.begin(); mpos != _molecules.end(); ++mpos) {
            RCSB::Chain* chain = (*mpos)->GetFirstChain();
            while (chain) {
                 if (chain->chain_type() == "ATOMP") {
                      residue_mapping.clear();
                      for (unsigned int j = 0; j < chain->SeqLen(); ++j) {
                           chain->GetResidueList(j, residues);
                           if (residues.empty()) continue;
                           int end = NUM_AA_TERMINAL_ATOMS;
                           const char **terminal_atoms = _AA_Terminal_Atoms;
                           if (j == 0) {
                                terminal_atoms = _AA_C_Terminal_Atoms;
                                end = NUM_AA_C_TERMINAL_ATOMS;
                           } else if (j == (chain->SeqLen() - 1)) {
                                terminal_atoms = _AA_N_Terminal_Atoms;
                                end = NUM_AA_N_TERMINAL_ATOMS;
                           }
                           oxt_atoms.clear();
                           for (std::vector<RCSB::Residue*>::iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                                for (int k = 0; k < end; ++k) {
                                     if (!_ccDic->is_terminal_atom((*rpos)->ResName(), terminal_atoms[k])) continue;
     
                                     if (!strcmp(terminal_atoms[k], "H1")) {
                                          (*rpos)->find_atom(terminal_atoms[k], atoms);
                                          if (atoms.empty()) continue;
                                          for (std::vector<RCSB::Atom*>::iterator apos = atoms.begin(); apos != atoms.end(); ++apos) {
                                               (*apos)->set_atmtype("H");
                                               (*apos)->set_pdb_atmnam("H");
                                          }
                                          (*rpos)->UpdateIndices();
                                     } else if (/* rename_flag && */ !strcmp(terminal_atoms[k], "OXT")) {
                                          (*rpos)->remove_a_atom(terminal_atoms[k], atoms);
                                          for (std::vector<RCSB::Atom*>::iterator apos = atoms.begin(); apos != atoms.end(); ++apos) {
                                               oxt_atoms.push_back(*apos);
                                          }
                                     } else (*rpos)->delete_a_atom(terminal_atoms[k], delete_atom_set);
                                }
                           }

                           if (oxt_atoms.empty()) continue;

                           if (((j + 1) < chain->SeqLen()) && (chain->SeqRes(j+1)->ResIndex < 0) && rename_flag) {
                                _FIELD* seqres = chain->SeqRes(j+1);
                                RCSB::Residue* residue = new RCSB::Residue;
                                residue->setCCDic(_ccDic);
                                residue->setLog(_logIo);
                                residue->setMessage(&_error_messages);
                                residue->set_token(chain->chain_type());
                                residue->set_ResName(seqres->Field[0]);
                                residue->set_chnid(chain->ChainID());
                                residue->set_res_no(seqres->Field[2]);
                                residue->set_pdb_chnid(chain->PDB_ChainID());
                                residue->set_pdb_res_no(seqres->Field[4]);
                                residue->set_ins_code(seqres->InsCode);
                                for (std::vector<RCSB::Atom*>::iterator apos = oxt_atoms.begin(); apos != oxt_atoms.end(); ++apos) {
                                     (*apos)->set_atmtype("N");
                                     (*apos)->set_pdb_atmnam("N");
                                     (*apos)->set_atom_type("N");
                                     (*apos)->set_restype(seqres->Field[0]);
                                     (*apos)->set_pdb_resnam(seqres->Field[0]);
                                     (*apos)->set_resnum(seqres->Field[2]);
                                     (*apos)->set_pdb_resnum(seqres->Field[4]);
                                     (*apos)->set_ins_code(seqres->InsCode);
                                     residue->insert_a_atom(*apos, (*mpos)->Mol_ID());
                                }
                                (*mpos)->insert_a_residue(residue);
                                residue_mapping.insert(std::make_pair(j + 1, residue));
                           } else {
                                for (std::vector<RCSB::Atom*>::iterator apos = oxt_atoms.begin(); apos != oxt_atoms.end(); ++apos) {
                                     std::string msg = "Deleting atom '" + (*apos)->pdb_chnid() + " " + (*apos)->restype() + " " + (*apos)->pdb_resnum()
                                                     + (*apos)->ins_code() +  " " + (*apos)->pdb_atmnam() + "'.\n";
                                     _logIo->message(msg.c_str());
     
                                     delete_atom_set.insert((long) (*apos));
                                     delete *apos;
                                }
                           }
                      }
                      if (!residue_mapping.empty()) chain->insert_OXT2N_residue(residue_mapping);
                 } else if (chain->chain_type() == "ATOMN") {
                      for (unsigned int j = 0; j < chain->SeqLen(); ++j) {
                           chain->GetResidueList(j, residues);
                           if (residues.empty()) continue;
                           int end = NUM_NA_TERMINAL_ATOMS;
                           const char **terminal_atoms = _NA_Terminal_Atoms;
                           if (j == 0) {
                                terminal_atoms = _NA_3_Terminal_Atoms;
                                end = NUM_NA_3_TERMINAL_ATOMS;
                           } else if (j == (chain->SeqLen() - 1)) {
                                terminal_atoms = _NA_5_Terminal_Atoms;
                                end = NUM_NA_5_TERMINAL_ATOMS;
                           }
                           for (std::vector<RCSB::Residue*>::iterator rpos = residues.begin(); rpos != residues.end(); ++rpos) {
                                for (int k = 0; k < end; ++k) {
                                     if (!_ccDic->is_terminal_atom((*rpos)->ResName(), terminal_atoms[k])) continue;
     
                                     (*rpos)->delete_a_atom(terminal_atoms[k], delete_atom_set);
                                }
                           }
                      }
                 }
                 chain = (*mpos)->GetNextChain();
            }
       }

       _update_linkage_information(delete_atom_set, _links);
       _update_linkage_information(delete_atom_set, _ssbonds);
       _update_linkage_information(delete_atom_set, _sltbrgs);
       // _update_linkage_information(delete_atom_set, _bspairs);

       if (delete_atom_set.empty()) return;

       for (std::list<_PDBX_REMEDIATION_ATOM_SITE_MAPPING>::iterator lpos = _pdbx_remediation_atom_site_mapping.begin();
            lpos != _pdbx_remediation_atom_site_mapping.end(); ++lpos) {
            if (delete_atom_set.find((long) lpos->atom) != delete_atom_set.end()) {
                 lpos = _pdbx_remediation_atom_site_mapping.erase(lpos);
                 --lpos;
            }
       }

       for (std::list<_PDBX_ENTITY_INSTANCE_FEATURE>::iterator lpos = _pdbx_entity_instance_feature.begin();
            lpos != _pdbx_entity_instance_feature.end(); ++lpos) {
            atoms.clear();
            for (std::vector<RCSB::Atom*>::const_iterator apos = lpos->atoms.begin(); apos != lpos->atoms.end(); ++apos) {
                 if (delete_atom_set.find((long) (*apos)) != delete_atom_set.end()) continue;
                 atoms.push_back(*apos);
            }
            if (atoms.size() < lpos->atoms.size()) lpos->atoms = atoms;
       }
}

void AnnotationObj::_update_linkage_information(const std::set<long>& atom_set, std::list<_LINK>& links)
{
       if (atom_set.empty() || links.empty()) return;

       for (std::list<_LINK>::iterator lpos = links.begin(); lpos != links.end(); ++lpos) {
            if (atom_set.find((long) lpos->fstAtom) != atom_set.end() ||
                atom_set.find((long) lpos->sndAtom) != atom_set.end()) {
                 lpos = links.erase(lpos);
                 --lpos;
            }
       }
} 

void AnnotationObj::_update_struct_category(Block& block)
{
       ISTable *t = getTablePtr(block, "struct");
       if (!t) return;

       std::string descriptor;
       get_value_clean(descriptor, t, 0, "pdbx_descriptor");
       if (!descriptor.empty()) return;

       descriptor = _get_struct_descriptor(block);

       if (descriptor.empty()) return;

       std::vector<std::string> missing_item;
       missing_item.clear();
       missing_item.push_back("pdbx_descriptor");
       check_missing_item(t, missing_item);

       t->UpdateCell(0, "pdbx_descriptor", descriptor);
       block.WriteTable(t);
}

void AnnotationObj::_update_pdbx_struct_conn_angle(Block& block)
{
       if (_links.empty() || _molecules.empty()) {
            deleteTable(block, "pdbx_struct_conn_angle");
            return;
       }

       // Atom lists: first atom is metal atom, the rest are polar atoms for each record
       std::vector<std::vector<std::pair<RCSB::Atom, std::string> > > MetalList;
       MetalList.clear();

       // index for MetalList: key=metal atom name index, val=position in MetalList
       std::map<std::string, unsigned int> MetalIndex;
       MetalIndex.clear();

       // Unique metal atom name index set
       std::set<std::string> atomIndex;
       atomIndex.clear();

       RCSB::Atom metal_atom, polar_atom;
       std::string SymOP_metal, SymOP_polar;
       std::vector<RCSB::Atom*> atoms;
       std::vector<std::pair<RCSB::Atom, std::string> > tMetalList;
       for (std::list<_LINK>::const_iterator lpos = _links.begin(); lpos != _links.end(); ++lpos) {
            if (lpos->type != "metalc") continue;
            if (!get_atom_pair(*lpos, _cell, metal_atom, polar_atom, SymOP_metal, SymOP_polar)) continue;

            std::string index = metal_atom.getAtomAllIndex();
            std::map<std::string, unsigned int>::const_iterator
                mpos = MetalIndex.find(index);
            if (mpos != MetalIndex.end())
                 // insert polar atom into existing metal atom cluster
                 MetalList[mpos->second].push_back(std::make_pair(polar_atom, SymOP_polar));
            else {
                 // add new metal atom cluster
                 MetalIndex.insert(std::make_pair(index, MetalList.size()));
                 tMetalList.clear();
                 tMetalList.push_back(std::make_pair(metal_atom, SymOP_metal));
                 tMetalList.push_back(std::make_pair(polar_atom, SymOP_polar));
                 MetalList.push_back(tMetalList);
            }

            mpos = MetalIndex.find(index);
            if (mpos == MetalIndex.end()) continue;

            if (atomIndex.find(index) != atomIndex.end()) continue;

            atomIndex.insert(index);

            RCSB::Residue* res = _molecules[0]->find_pdb_residue(metal_atom.pdb_chnid(), metal_atom.pdb_resnam(), metal_atom.pdb_resnum(),
                                                                 metal_atom.ins_code());
            if (!res) continue;

            std::string alt_loc = polar_atom.alt_loc();
            try {
                 // find all linked atoms if metal atom is belongs to a metal-complex residue
                 const ConnectFormat& drug = _ccDic->find_drug(metal_atom.pdb_resnam());
                 std::map<std::string, std::vector<std::string> > LinkedAtoms = drug.getLinkedAtoms();
                 if (LinkedAtoms.empty()) continue;
                 std::map<std::string, std::vector<std::string> >::const_iterator mpos1 = LinkedAtoms.find(metal_atom.pdb_atmnam());
                 if (mpos1 == LinkedAtoms.end()) continue;
                 for (std::vector<std::string>::const_iterator rpos = mpos1->second.begin(); rpos != mpos1->second.end(); ++rpos) {
                      res->find_atom(*rpos, atoms);
                      if (atoms.empty()) continue;
                      for (std::vector<RCSB::Atom*>::const_iterator apos = atoms.begin(); apos != atoms.end(); ++apos) {
                           if ((*apos)->alt_loc() != alt_loc && (*apos)->alt_loc() != metal_atom.alt_loc() &&
                              !(*apos)->alt_loc().empty() && !alt_loc.empty() && !metal_atom.alt_loc().empty()) continue;
                           polar_atom = *(*apos);
                           _cell.symmetry_operation(polar_atom, SymOP_metal);
                           MetalList[mpos->second].push_back(std::make_pair(polar_atom, SymOP_metal));
                           break;
                      }
                 }
            } catch (const std::exception& exc) {}
       }

       int row = 0;
       ISTable *t = _newTablePtr("pdbx_struct_conn_angle");

       for (std::vector<std::vector<std::pair<RCSB::Atom, std::string> > >::iterator pos = MetalList.begin(); pos != MetalList.end(); ++pos) {
            if (pos->size() < 3) continue;

            for (unsigned int i = 1; i < pos->size(); ++i) {
                 for (unsigned int j = 1; j < i; ++j) {
                      double angle = cal_angle((*pos)[j].first, (*pos)[0].first, (*pos)[i].first);
                      t->AddRow();
                      t->UpdateCell(row, "id", String::IntToString(row + 1));
                      UpdateAtomInfo::UpdateTable(t, row, (*pos)[j].first, "ptnr1_", NUM_ALL_ITEM + 1);
                      if (!(*pos)[j].second.empty())
                           t->UpdateCell(row, "ptnr1_symmetry", (*pos)[j].second);
                      else t->UpdateCell(row, "ptnr1_symmetry", "1_555");
                      UpdateAtomInfo::UpdateTable(t, row, (*pos)[0].first, "ptnr2_", NUM_ALL_ITEM + 1);
                      if (!(*pos)[0].second.empty())
                           t->UpdateCell(row, "ptnr2_symmetry", (*pos)[0].second);
                      else t->UpdateCell(row, "ptnr2_symmetry", "1_555");
                      UpdateAtomInfo::UpdateTable(t, row, (*pos)[i].first, "ptnr3_", NUM_ALL_ITEM + 1);
                      if (!(*pos)[i].second.empty())
                           t->UpdateCell(row, "ptnr3_symmetry", (*pos)[i].second);
                      else t->UpdateCell(row, "ptnr3_symmetry", "1_555");
                      t->UpdateCell(row, "value", FloatToString(angle, 0, 1));
                      row++;
                 }
            }
       }

       if (row == 0) {
            delete t;
            deleteTable(block, "pdbx_struct_conn_angle");
       } else block.WriteTable(t);
}

void AnnotationObj::_cif_update_pdbx_data_processing_status(Block& block)
{
       std::set<std::string> allow_set;
       allow_set.clear();
       allow_set.insert("helix");
       allow_set.insert("link");
       allow_set.insert("sheet");
       allow_set.insert("site");
       allow_set.insert("solvent position");
       allow_set.insert("ssbond");

       _cif_update_add_value_to_pdbx_data_processing_status(block, allow_set);
}

void AnnotationObj::_cif_update_add_value_to_pdbx_data_processing_status(Block& block, const std::set<std::string>& value_set)
{
       std::map<std::string, std::set<std::string> > skip_task_names;
       skip_task_names.clear();

       std::set<std::string> t_set;

       ISTable *t = getTablePtr(block, "pdbx_data_processing_status");
       if (t) {
            std::string status, task_name;
            int rowNo = t->GetNumRows();
            for (int i = 0; i < rowNo; ++i) {
                 get_value_clean_lower(status, t, i, "status");
                 if (status.empty()) continue;
                 get_value_clean_lower(task_name, t, i, "task_name");
                 if (task_name.empty()) continue;
                 std::map<std::string, std::set<std::string> >::iterator mpos = skip_task_names.find(status);
                 if (mpos != skip_task_names.end()) mpos->second.insert(task_name);
                 else {
                      t_set.clear();
                      t_set.insert(task_name);
                      skip_task_names.insert(std::make_pair(status, t_set));
                 }
            }
       }

       std::map<std::string, std::set<std::string> >::iterator mpos = skip_task_names.find("skip");
       if (mpos != skip_task_names.end()) {
            for (std::set<std::string>::const_iterator spos = value_set.begin(); spos != value_set.end(); ++spos) {
                 mpos->second.insert(*spos);
            }
       } else skip_task_names.insert(std::make_pair("skip", value_set));

       t = new ISTable("pdbx_data_processing_status");
       t->AddColumn("status");
       t->AddColumn("task_name");
       int row = 0;
       for (std::map<std::string, std::set<std::string> >::const_iterator mpos = skip_task_names.begin(); mpos != skip_task_names.end(); ++mpos) {
            for (std::set<std::string>::const_iterator spos = mpos->second.begin(); spos != mpos->second.end(); ++spos) {
                 t->AddRow();
                 t->UpdateCell(row, "status", mpos->first);
                 t->UpdateCell(row, "task_name", *spos);
                 row++;
            }
       }
       block.WriteTable(t);
}

void AnnotationObj::_update_struct_ncs_oper(Block& block)
{
       ISTable *t = getTablePtr(block, "struct_ncs_oper");
       if (!t) return;

       int id_value;
       std::string cs;
       bool found_non_integer_value = false;
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(cs, t, i, "id");
            if (cs.empty() || !IsInteger(cs, id_value)) {
                 found_non_integer_value = true;
                 break;
            }
       }
       if (!found_non_integer_value) return;

       std::map<std::string, std::string> id_mapping;
       id_mapping.clear();
       for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
            get_value_clean(cs, t, i, "id");
            std::string id = String::IntToString(i + 1);
            id_mapping.insert(std::make_pair(cs, id));
            t->UpdateCell(i, "id", id);
       }
       block.WriteTable(t);

       std::vector<std::string> child_categories;
       child_categories.clear();
       child_categories.push_back("pdbx_struct_ncs_virus_gen");
       child_categories.push_back("struct_ncs_ens_gen");
       for (std::vector<std::string>::const_iterator vpos = child_categories.begin(); vpos != child_categories.end(); ++vpos) {
            t = getTablePtr(block, *vpos);
            if (!t) continue;

            for (unsigned int i = 0; i < t->GetNumRows(); ++i) {
                 get_value_clean(cs, t, i, "oper_id");
                 std::map<std::string, std::string>::const_iterator mpos = id_mapping.find(cs);
                 if (mpos == id_mapping.end()) continue;
                 t->UpdateCell(i, "oper_id", mpos->second);
            }
            block.WriteTable(t);
       }
}

static bool get_atom_pair(const _LINK& link, CrySymmetry& cell, RCSB::Atom& metal_atom, RCSB::Atom& polar_atom,
                          std::string& SymOP_metal, std::string& SymOP_polar)
{
       if (link.mol_index != 0) return false;

       if (Element::findMetalFlag(link.fstAtom->atom_type())) {
            SymOP_metal = link.SymOP_1;
            metal_atom = *(link.fstAtom);
            cell.symmetry_operation(metal_atom, link.SymOP_1);
            SymOP_polar = link.SymOP_2;
            polar_atom = *(link.sndAtom);
            cell.symmetry_operation(polar_atom, link.SymOP_2);
       } else if (Element::findMetalFlag(link.sndAtom->atom_type())) {
            SymOP_metal = link.SymOP_2;
            metal_atom = *(link.sndAtom);
            cell.symmetry_operation(metal_atom, link.SymOP_2);
            SymOP_polar = link.SymOP_1;
            polar_atom = *(link.fstAtom);
            cell.symmetry_operation(polar_atom, link.SymOP_1);
       } else return false;

       return true;
}
