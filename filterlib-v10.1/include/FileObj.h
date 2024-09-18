/*
FILE:     FileObj.h
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
#ifndef _FILEOBJ_H_
#define _FILEOBJ_H_

#include <stdio.h>
#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Assembly.h"
#include "CategoryMappingDef.h"
#include "CifFile.h"
#include "ConnectDic.h"
#include "CrySymmetry.h"
#include "DictSdbUtil.h"
#include "Entity.h"
#include "LogUtil.h"
#include "MessageUtil.h"
#include "Molecule.h"
#include "SymmMatrix.h"

/* Changing the experinemtal type values also need to change the experimental_type* values in merge_supplemental_category_list.cif accordingly. */

#define EXPERIMENT_TYPE_BASIC           1
#define EXPERIMENT_TYPE_XRAY            2
#define EXPERIMENT_TYPE_NEUTRON         4
#define EXPERIMENT_TYPE_FIBER           8
#define EXPERIMENT_TYPE_ELECTRON       16
#define EXPERIMENT_TYPE_EC             32
#define EXPERIMENT_TYPE_EM             64
#define EXPERIMENT_TYPE_ET            128
#define EXPERIMENT_TYPE_NMR           256
#define EXPERIMENT_TYPE_NMR_SOLID     512
#define EXPERIMENT_TYPE_SOLN_SCT     1024
#define EXPERIMENT_TYPE_MODEL        2048

class FileObj
{
   private:
       int _read_atom_site_record_number;
       int _write_atom_site_record_number;
       bool _check_atom_site_record_number_flag;

       // key: AAResName_AAAtomName_SugarResName_SugarAtomName
       // value: "N-Glycosylation", "O-Glycosylation", "C-Mannosylation"
       std::map<std::string, std::string> _glyco_site_map;

   protected:
       LogUtil *_logIo;
       ConnectDic *_ccDic;
       MessageUtil _error_messages;
       std::string _rcsbroot;
       std::string _input_filename;
       std::string _output_filename;
       std::string _StructureId;
       std::string _UppercasePdbId;
       std::string _LowercasePdbId;
       std::string _firstBlockName;
       std::string _resetFirstBlockName;
       int _input_format;
       int _output_format;
       int _experiment_type;
       bool _DepUI_Flag;
       std::map<std::string, std::list<std::vector<std::string> > > _pdb_records;
       std::vector<RCSB::Molecule*> _molecules;
       bool _caseSensitive;
       bool _rename_residue_flag;
       bool _find_sugar_chain_flag;
       bool _skip_split_flag;
       CifFile* _CifObj;
       DictUtil _dictUtil;

       // map key: PDB chain ID
       std::map<std::string, SEQ> _Seqs;

       // mak key: PDB chain ID
       // value[0]: pdb_strand_id
       // value[1]: mon_id
       // value[2]: pdb_mon_id
       // value[3]: pdb_seq_num
       // value[4]: pdb_ins_code
       // value[5]: entity_id
       std::map<std::string, std::vector<std::vector<std::string> > > _Seq_Scheme_Mapping;

       std::map<std::string, BRANCH_INFO> _Branch_Seq_Scheme_Mapping;

       std::map<std::string, int> _ExpTypeMapping;
       std::map<int, std::string> _TypeExpMapping;
       CrySymmetry _cell;
       double _scale[3][3], _vec[3];
       bool _exist_scale_matrix;
       std::vector<SymmMatrix> _SymmMatrices;
       std::map<std::string, unsigned int> _SymmMatrix_Mapping;
       std::vector<Assembly> _assemblies;
       bool _update_symmetry_flag;

       int _max_entity_id;
       // key: integer entity ID
       // value: entity information
       std::map<int, Entity> _entities;
       // key: entity_key
       // value: integer entity ID
       std::map<std::string, int> _entity_ids;

       std::set<std::string> _skip_task_set, _glyco_site_res_pair_set;

       // pair.first: new residue name
       // pair.second: old residue list
       //              vector[0]: PDB chain ID
       //              vector[1]: Residue Name
       //              vector[2]: Residue Number
       //              vector[3]: Insertion Code
       std::vector<std::pair<std::string, std::vector<std::vector<std::string> > > > _merge_residue_list;

       // value: pdb_chnid1_pdb_resnam1_pdb_resnum1_ins_code1_pdb_chnid2_pdb_resnam2_pdb_resnum2_ins_code2
       std::set<std::string> _link_residue_set;

       std::vector<std::string> _extra_atom_site_item_list;

       void _initialize_reference_variables();
       void _initialize_other_variables();
       void _initExpTypeMapping();
       std::string _get_diffraction_experiment_method();
       void _clear_allocated_memories();
       void _clearSEQ(SEQ& seq);

       void _fill_in_glyco_site_info();
       bool _is_glycosylation_site(const std::string& idx);
       bool _is_glycosylation_site(RCSB::Atom* atom1, RCSB::Atom* atom2);
       std::string _find_glycosylation_type(const std::string& idx);
       std::string _find_glycosylation_type(RCSB::Atom* atom1, RCSB::Atom* atom2);

       void _read_mmcif_metadata(const bool& coord_and_seq_only = false, const bool& update_database_PDB_rev_flag = true);
       void _legacy_vs_v5_conversion(Block& _cifblock, const int& status);
       void _read_extra_dataBlock();
       bool _verify_chem_comp_info(Block& block, const std::string& comp_id);
       std::string _check_chem_comp_error(const std::string& comp_id, const std::vector<std::string>& items,
                                          const std::vector<std::vector<std::string> >& values);
       void _read_mmcif_metadata_process(Block& block, const bool& coord_and_seq_only = false, const bool& update_database_PDB_rev_flag = true);
       void _getExperimentType(const std::string& exp_type);
       void _updateCrySymmetry();
       void _getSequenceType();
       std::string _getSequenceType(const std::vector<std::string>& res);

       void _getExperimentTypefromEXPDTA();
       void _getSequencefromSEQRES();
       void _getCrySymmetryfromCRYST1();
       void _getAuthorDefinedScaleMatrixfromSCALE();

       void _getExperimentTypefrom_exptl(Block& block);
       void _getSequenceInformation(Block& block);
       void _getSequenceInformationfrom_pdbx_poly_seq_scheme(Block& block);
       void _getSequenceInformationfrom_entity_poly(Block& block);
       void _getCrySymmetryfrom_symmetry_and_cell(Block& block);
       void _getAuthorDefinedScaleMatrixfrom_atom_sites(Block& block);
       void _update_database_PDB_rev(Block& block);
       void _reformat_one_letter_code_sequence(Block& block);
       void _switch_high_low_resolution(Block& block);
       void _remove_extra_space_in_name(Block& block);
       void _reorder_refine_ls_shell_reflns_shell(Block& block, const std::string& category);
       void _add_start_date(Block& block);
       void _getResiduesfromNames(const int& mol_id, const std::vector<std::vector<std::string> >&
                                  data, int& mol_index, std::vector<RCSB::Residue*>& residues, 
                                  std::vector<std::string>& errors);
       void _getAtomsfromNames(const int& mol_id, const std::vector<std::vector<std::string> >& data, int& mol_index, std::vector<RCSB::Atom*>& atoms,
                               std::vector<std::string>& errors, const bool& origin_flag=false);

       /**
       **  Retrieves a pointer to the table.
       **
       **  \param[in]: block - reference to a coordinate data block 
       **  \param[in]: catName - category name
       **
       **  \return Pointer to the table, if table was found
       **  \return NULL, if table was not found
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       ISTable* _getTablePtr(Block& block, const std::string& catName);

       /**
       **  Creates a new table.
       **
       **  \param[in]: cifcategory - category definition
       **
       **  \return Pointer to the table
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       ISTable* _newTablePtr(const CIF_CATEGORY& cifcategory);

       /**
       **  Creates a new table.
       **
       **  \param[in]: catName - category name
       **
       **  \return Pointer to the table, if catName is defined in mapping file ndb_cif.cif
       **  \return NULL, if table was not
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       ISTable* _newTablePtr(const std::string& catName);

       std::string _insert_symm_matrix(std::vector<SymmMatrix>& matrices);
       std::string _insert_symm_matrix(SymmMatrix& symmtrix);
       void _update_symm_matrix(SymmMatrix& symmtrix);
       void _read_struct_assembly_categories(Block& block);
       void _read_pdbx_struct_oper_list(Block& block);
       void _write_struct_assembly_categories(Block& block);
       void _write_pdbx_struct_oper_list(Block& block);
       std::string _get_asym_id_list_from_chain_ids(const std::vector<std::string>& chain_id_list);
       void _get_chain_ids_from_asym_id_list(const std::string& asym_id_list,
                                             std::vector<std::string>& chain_id_list);
       void _check_atom_number_count(std::vector<int>& atom_counts);
       bool _readInstanceAssignFile(const std::string& assignfile, std::map<std::string, std::pair<std::string,
                         std::map<std::string, std::string> > >& atom_mapping, const bool& remove_mol_id_flag);

       void _read_branch_polymer_info();
       void _read_branch_scheme_list_category(ISTable* t, const std::vector<std::string>& items, std::map<std::string,
                                              std::vector<std::vector<std::string> > >& mapping);
       void _read_branch_link_descriptor_category(ISTable*t, const std::vector<std::string>& items, std::map<std::string,
                                                  std::vector<std::map<std::string, std::string> > >& mapping);
       void _check_branch_list_with_link(const std::map<std::string, std::vector<std::map<std::string, std::string> > >& link_mapping,
                                         std::map<std::string, std::vector<std::vector<std::string> > >& list_mapping);
   public:
       FileObj();
       ~FileObj();
       const int& experiment_type() const { return _experiment_type; }
       void setLog(LogUtil* logPt);
       void setCCDic(ConnectDic* ccdic);
       void setRCSBROOT(const std::string&);
       void setDictSdbPath(const std::string&);
       void set_input_filename(const std::string&);
       void set_output_filename(const std::string&);
       void set_input_format(const int);
       void set_output_format(const int);
       void set_depUI_flag(const bool& flag);
       void set_CaseSensitive();
       void setRenameResidueFlag(const bool& flag);
       void setFindSugarChainFlag(const bool& flag);
       void setSkipSplitFlag(const bool& flag);
       void set_merge_residue_info(const std::vector<std::string>&);
       void read_pdb_file();
       bool parse_mmcif_file();
       void clear_branch_polymer_info();
       bool read_mmcif_file(const bool& coord_and_seq_only = false, const bool& check_format_flag = false);
       void read_mmcif_file_without_parsing(const bool& coord_and_seq_only = false);
       bool read_mmcif_metadata(const bool& coord_and_seq_only = false);
       void read_mmcif_coordinate(const bool& check_format_flag = false);
       bool is_pdb_format_compatible();
       bool build_molecule(const bool& skip_update_flag = false, const bool& split_flag = false, const std::string& assignfile = "");
       void build_molecule_only(const bool& using_coordinate_only = false);

       void assign_prev_entity_id();
       void write_mmcif_file(const bool& internal_flag = false, const bool& double_quoting_flag = true,
                             const bool& start_date_flag = true);
       void print_data();
};

#endif
