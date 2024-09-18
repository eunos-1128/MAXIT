/*
FILE:     Maxit.h
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
#ifndef _H_MAXIT_H_
#define _H_MAXIT_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "AnnotationObj.h"
#include "Remark.h"

#define FORCE_DECIMAL_PRECISION       0
#define UNFORCE_DECIMAL_PRECISION     1

#define NUM_PDB_SOURCE               36
#define NUM_MER                      20

class Maxit: public AnnotationObj {
   public:
       /**
       **  Constructs Maxit
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
       Maxit();

       /**
       **  Destructs Maxit
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
       ~Maxit();

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
       void clear();

       /**
       **  Clear PDB record data storage
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
       void clear_pdb_records();

       /**
       **  Set PDB Header file only flag
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
       void set_pdb_header_only_flag();

       /**
       **  Map NDB token records into cif categories
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
       void ndb_to_cif();

       /**
       **  Map cif categories into NDB token records
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
       void cif_to_ndb();

       /**
       **  Map PDB records into NDB token records
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
       void pdb_to_ndb();

       /**
       **  Map NDB records into PDB token records
       **
       **  \param[in]: biological_flag - flag to indicate to generate biological unit pdb file
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void ndb_to_pdb(const bool& biological_flag = false);

       /**
       **  Write PDB format file
       **
       **  \param[in]: segid_flag - Using SEGID flag
       **  \param[in]: checking_linkage_flag - checking bonded atoms' distance flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void write_pdb_file(const bool& segid_flag, const bool& checking_linkage_flag);

       /**
       **  Read PDB structural feature records
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
       void Read_PDB_StructuralFeatures();

   protected:
       int _remark_decimal_precision;
       bool _pdb_header_only_flag;

       // Joint refinement methods
       std::vector<std::string> _joint_methods;

       // refinement method
       std::string _current_method;

       // key: remark number
       // value: remark context
       std::map<int, std::vector<std::string> > _remarks;

       // Original PDB entry ID(s) for re-refinement entry
       std::vector<std::string> _original_entry_ids;

       std::map<int, std::list<std::vector<std::string> > > _missing_polymer_residues;
       std::map<int, std::list<std::vector<std::string> > > _missing_polymer_atoms;
       std::map<int, std::list<std::vector<std::string> > > _zero_occ_polymer_residues;
       std::map<int, std::list<std::vector<std::string> > > _zero_occ_polymer_atoms;
       std::map<int, std::list<std::vector<std::string> > > _missing_non_polymer_residues;
       std::map<int, std::list<std::vector<std::string> > > _zero_occ_non_polymer_residues;

       // key: HET ID
       // value: vector[0]: PRD ID
       //        vector[1]: Type + Class
       //        vector[2]: Details
       std::map<std::string, std::vector<std::string> > _non_polymer_prd_info;

       /**
       **  Get refinement program information
       **
       **  \param[out]: is_refmac5 - REFMAC 5+ program flag
       **  \param[out]: is_phenix  - PHENIX program flag
       **  \param[out]: is_cns_xplor - CNS/X-PLOR program flag
       **  \param[out]: is_buster  - BUSTER program flag
       **
       **  \return selected refinement program
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _get_refinement_program(bool& is_refmac5, bool& is_phenix, bool& is_cns_xplor, bool& is_buster);

       /**
       **  Update Ref_ID & Align_ID
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
       void _ndb_to_cif_update_Ref_ID_and_Align_ID();

       /**
       **  Update TMLPTL records from RFTPOS, RFMPOS, RFLPOS, RFTTHR, RFMTHR
       **             & RFLTHR records for REFMAC 5+ program
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
       void _ndb_to_cif_add_TMLPTL_token_for_REFMAC_5();

       /**
       **  Update NCSGRO & TMLPTL records from PHENCS records for PHENIX program
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
       void _ndb_to_cif_update_TMLPTL_and_NCSGRO_tokens_for_PHENIX();

       /**
       **  Update PDBRMK for pdbx_database_remark 0, 650 & 700
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
       void _ndb_to_cif_capture_pdbx_database_remark();

       /**
       **  Update Refine_ID for BVALUE record
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
       void _ndb_to_cif_update_BVALUE_Refine_ID();

       /**
       **  Mapping PDB/NDB records to cif categories based on mapping file ndb_cif.cif
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: refmac_phenix_flag - REFMAC 5+ or PHENIX refinement program
       **                                   flag
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_cif_mapping(Block& block, const bool& refmac_phenix_flag);

       /**
       **  Update refine_ls_restr category
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
       void _ndb_to_cif_process_refine_ls_restr(Block& block);

       /**
       **  Update database_PDB_rev & database_PDB_rev_record categories from
       **           REVDAT record
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
       void _ndb_to_cif_process_database_pdb_rev(Block& block);

       /**
       **  Update pdbx_database_related category from RENTRY & REVDAT records
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
       void _ndb_to_cif_process_pdbx_database_related(Block& block);

       /**
       **  Update matrix-related categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: Category - cif category definition
       **  \param[in]: tokenid - correspondence PDB record token name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_cif_process_matrix(Block& block, const CIF_CATEGORY &Category, const std::string& tokenid);

       /**
       **  Update name-related categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: Category - cif category definition
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_cif_process_name(Block& block, const CIF_CATEGORY &Category);

       /**
       **  Update refine_ls_restr_ncs category for REFMAC & PHENIX program
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: Category - cif category definition
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_cif_process_refmac_refine_ls_restr_ncs(Block& block,
                                         const CIF_CATEGORY &Category);

       /**
       **  Update general categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: Category - cif category definition
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_cif_process_general(Block& block, const CIF_CATEGORY &Category);

       /**
       **  Update database_PDB_remark category
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: Category - cif category definition
       **  \param[in]: tokenid - correspondence PDB record token name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_cif_process_remark(Block& block, const CIF_CATEGORY &Category);
                                       // const std::string& tokenid);

       /**
       **  Update general categories
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: Category - cif category definition
       **  \param[in]: tokenid - correspondence PDB record token name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_cif_process_general_category(Block& block, const CIF_CATEGORY& Category, const std::string& tokenid);

       /**
       **  Update _exptl_crystal_grow.pdbx_details to include _exptl_crystal_grow.pH, 
       **           _exptl_crystal_grow.method & _exptl_crystal_grow.temp info.
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
       void _cif_to_ndb_update_exptl_crystal_grow(Block& block);

       /**
       **  Update PDB numbering in struct_ref_seq & struct_ref_seq_dif categories
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
       void _cif_to_ndb_update_struct_ref_seq_and_diff(Block& block);

       /**
       **  Mapping cif categories to PDB/NDB records based on mapping file ndb_cif.cif
       **
       **  \param[out]: block - reference to a coordinate data block
       **  \param[in]: refmac_phenix_flag - REFMAC 5+ or PHENIX refinement program
       **                                   flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_mapping(Block& block, const bool& refmac_phenix_flag);

       /**
       **  Update various PDB/NDB records
       **
       **  \param[in]: is_cns_xplor - CNS/X-PLOR refinement program flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_post_processing(const bool& is_cns_xplor);

       /**
       **  Map database_PDB_remark to PREMRK & pdbx_database_remark to PDBRMK
       **
       **  \param[in]: t - pointer to cif category table
       **  \param[in]: tokenid - PDB record token name
       **  \param[in]: field_no - record field number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_remark(ISTable *t, const std::string& tokenid, const int& field_no);

       /**
       **  Map database_PDB_rev & database_PDB_rev_record gategories to REVDAT records
       **
       **  \param[in]: block - reference to a coordinate data block
       **  \param[in]: t - pointer to "database_PDB_rev" table
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_database_pdb_rev(Block& block, ISTable *t);

       /**
       **  Map pdbx_database_related category to RENTRY/SPLIT records
       **
       **  \param[in]: t - pointer to "pdbx_database_related" table
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_rcsb_database_related(ISTable *t);

       /**
       **  Map cif category to single PDB/NDB record
       **
       **  \param[in]: t - pointer to table
       **  \param[in]: Category - cif category definition
       **  \param[in]: clean_string_flag - flag to indicate whether need to compress
       **                                  while space
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_general_category(ISTable* t, const CIF_CATEGORY& Category, const bool& clean_string_flag);

       /**
       **  Map cif category to single PDB/NDB record
       **
       **  \param[in]: t - pointer to table
       **  \param[in]: field_item_mapping - field number vs. item name mapping
       **  \param[in]: clean_string_flag - flag to indicate whether need to compress
       **                                  while space
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_single_token(ISTable* t, const std::string& tokenid, const std::map<int, std::string>& field_item_mapping,
                                             const bool& clean_string_flag);

       /**
       **  Map refine_ls_restr category to various PDB/NDB records
       **
       **  \param[in]: t - pointer to "refine_ls_restr" table
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_refine_ls_restr(ISTable * t);

       /**
       **  Map refine_ls_restr_ncs category to TMLPTL records
       **
       **  \param[in]: t - pointer to "refine_ls_restr_ncs" table
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_refmac_refine_ls_restr_ncs(ISTable * t);

       /**
       **  Map cif category to PDB/NDB records
       **
       **  \param[in]: t - pointer to table
       **  \param[in]: Category - cif category definition
       **  \param[in]: is_special - special flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_category(ISTable* t, const CIF_CATEGORY& Category, const bool& is_special);

       /**
       **  Map cif category to PDB/NDB records
       **
       **  \param[in]: t         - pointer to table
       **  \param[in]: Category  - cif category definition
       **  \param[in]: irow      - table row index
       **  \param[in]: tokenid   - PDB/NDB token name
       **  \param[out]: idvalue  - key ids
       **  \param[out]: SerialNo - key values
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_find_serial_no(ISTable* t, const CIF_CATEGORY& Category, const int& irow, const std::string& tokenid, std::vector<std::string>&
                                       idvalue, std::vector<std::string>& SerialNo);

       /**
       **  Insert value to PDB/NDB records
       **
       **  \param[in]: tokenid       - PDB/NDB token name
       **  \param[in]: fieldno       - field number
       **  \param[in]: SerialNoField - key field numbers
       **  \param[in]: SerialNo      - key values
       **  \param[in]: value         - value need to be imported to PDB record
       **  \param[in]: IsJrnl        - JRNL token flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_general(const std::string& tokenid, const int& fieldno, const std::vector<int>& SerialNoField,
                                        const std::vector<std::string>& SerialNo, const std::string& value);

       /**
       **  Insert value to PDB/NDB records
       **
       **  \param[in]: tokenid       - PDB/NDB token name
       **  \param[in]: fieldno       - field number
       **  \param[in]: SerialNoField - key field numbers
       **  \param[in]: value         - value need to be imported to PDB record
       **  \param[in]: SeqField      - sequential field number
       **  \param[out]: Field        - reference to record Field
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_value(const std::string& tokenid, const int& fieldno, const std::vector<int>& SerialNoField,
                                      const std::string& value, const int& SeqField, std::vector<std::string>& Field);

       /**
       **  Insert value from struct_ref category to DBREF records
       **
       **  \param[in]: tokenid       - PDB/NDB token name
       **  \param[in]: fieldno       - field number
       **  \param[in]: SerialNoField - key field numbers
       **  \param[in]: SerialNo      - key values
       **  \param[in]: value         - value need to be imported to PDB record
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_special(const std::string& tokenid, const int& fieldno, const std::vector<int>& SerialNoField,
                                        const std::vector<std::string>& SerialNo, const std::string& value);

       /**
       **  Map matrix categories to records
       **
       **  \param[in]: t - pointer to table
       **  \param[in]: Category - cif category definition
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cif_to_ndb_process_matrix(ISTable *t, const CIF_CATEGORY& Category);

       /**
       **  Add new record to _pdb_records
       **
       **  \param[in]: token - PDB token name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _addNewRecord(const std::string& token);

       /**
       **  Update value at first record
       **
       **  \param[in]: token    - PDB token name
       **  \param[in]: field_no - field number
       **  \param[in]: value    - new value
       **  \param[in]: force_flag - ture: force to create record if it does not exist &
       **                                 overwrites existing value
       **                           false: copies value only if record exists and
       **                                  value does not exist
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateRecordFront(const std::string& token, const int& field_no, const std::string& value, const bool& force_flag = true);

       /**
       **  Insert value at first record
       **
       **  \param[in]: token    - PDB token name
       **  \param[in]: field_no - field number
       **  \param[in]: value    - new value
       **  \param[in]: delimiter - delimiter between concantenated strings
       **  \param[in]: force_flag - ture: force to create record if it does not exist &
       **                                 overwrites existing value
       **                           false: copies value only if record exists and
       **                                  value does not exist
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateRecordFront(const std::string& token, const int& field_no, const std::string& value, const std::string& delimiter,
                               const bool& force_flag = true);

       /**
       **  Update value at last record
       **
       **  \param[in]: token    - PDB token name
       **  \param[in]: field_no - field number
       **  \param[in]: value    - new value
       **  \param[in]: force_flag - ture: force to create record if it does not exist &
       **                                 overwrites existing value
       **                           false: copies value only if record exists and
       **                                  value does not exist
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateRecordBack(const std::string& token, const int& field_no, const std::string& value, const bool& force_flag = true);

       /**
       **  Insert value at last record
       **
       **  \param[in]: token    - PDB token name
       **  \param[in]: field_no - field number
       **  \param[in]: value    - new value
       **  \param[in]: delimiter - delimiter between concantenated strings
       **  \param[in]: force_flag - ture: force to create record if it does not exist &
       **                                 overwrites existing value
       **                           false: copies value only if record exists and
       **                                  value does not exist
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateRecordBack(const std::string& token, const int& field_no, const std::string& value, const std::string& delimiter,
                              const bool& force_flag = true);

       /**
       **  Get value at first record
       **
       **  \param[in]: token    - PDB token name
       **  \param[in]: field_no - field number
       **  \param[out]: value    - retrieved value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _getRecordFront(const std::string& token, const int& field_no, std::string& value);

       /**
       **  Get value at last record
       **
       **  \param[in]: token    - PDB token name
       **  \param[in]: field_no - field number
       **  \param[out]: value    - retrieved value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _getRecordBack(const std::string& token, const int& field_no, std::string& value);

       /**
       **  Copy value from source token to target token
       **
       **  \param[in]: source_token - source PDB token name
       **  \param[in]: source_field - source field number
       **  \param[in]: target_token - target PDB token name
       **  \param[in]: target_field - target field number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _copyValue(const std::string& source_token, const int& source_field, const std::string& target_token, const int& target_field,
                       const bool& force_flag = true);

       /**
       **  Copy value from source token to target token
       **
       **  \param[in]: source_token - source PDB token name
       **  \param[in]: target_token - target PDB token name
       **  \param[in]: target_field - target field number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _copyValue(const std::string& source_token, const std::string& target_token, const std::map<int, int>& field_mapping,
                       const bool& force_flag = true);

       /**
       **  Update values
       **
       **  \param[in]: token            - target PDB token name
       **  \param[in]: key_field        -  target key field number
       **  \param[in]: key_value        - key value
       **  \param[in]: field_value_pair - field number/value pair
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateValue(const std::string& token, const int& key_field,  const std::string& key_value, const std::map<int, std::string>& field_value_pair);

       /**
       **  Get value
       **
       **  \param[in]: token      - target PDB token name
       **  \param[in]: key_field  -  target key field number
       **  \param[in]: key_value  - key value
       **  \param[in]: field_no   - field number
       **  \param[out]: value     - field value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _getValue(const std::string& token, const int& key_field,  const std::string& key_value, const int& fieldno, std::string& value);

       /**
       **  Reorder PDB records 
       **
       **  \param[in]: token   - PDB token name
       **  \param[in]: fieldno - field number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _orderRecordWithIntegerKey(const std::string& token, const int& fieldno);

       /**
       **  Reorder PDB records 
       **
       **  \param[in]: token   - PDB token name
       **  \param[in]: fieldno - field number
       **  \param[in]: reverse_order - reverse order flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _orderRecordWithFloatKey(const std::string& token, const int& fieldno, const bool& reverse_order = false);

       /**
       **  Reorder PDB records 
       **
       **  \param[in]: token   - PDB token name
       **  \param[in]: fieldno - field number
       **  \param[in]: order   - order value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _orderRecordWithOrderKey(const std::string& token, const int& fieldno, const std::vector<std::string>& order);

       /**
       **  Strip and compress white space in string
       **
       **  \param[in]: token   - PDB token name
       **  \param[in]: fieldno - field number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _cleanString(const std::string& token, const int& fieldno);

       /**
       **  Remove white space in string
       **
       **  \param[in]: token   - PDB token name
       **  \param[in]: fieldno - field number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _removeWhiteSpace(const std::string& token, const int& fieldno);

       /**
       **  Remove white space in name string
       **
       **  \param[in]: token   - PDB token name
       **  \param[in]: fieldno - field number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _removeWhiteSpaceBetweenName(const std::string& token, const int& fieldno);

       /**
       **  Set String to MixCase
       **
       **  \param[in]: token   - PDB token name
       **  \param[in]: fieldno - field number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _mixCase(const std::string& token, const int& fieldno);

       /**
       **  Copy R_Value_Obs to R_Value_Work and vs versa for CNS & X-PLOR program
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
       void _updateRValueForCNS_XPLOR();

       /**
       **  Update all NDB tokens with RefineField
       **
       **  \param[in]: method - refinement method
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateRefineMethod(const std::string& method);

       /**
       **  Remove empty record(s)
       **
       **  \param[in]: token      - PDB token name
       **  \param[in]: skip_field - skip field number set
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       // void _removeEmptyRecord(const std::string& token, const std::set<unsigned int>& skip_field);

       /**
       **  Remove empty records
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
       void _removeEmptyRecords();

       /**
       **  Update records
       **
       **  \param[in]: format - NDB File Format
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateRecords(const int& format);

       /**
       **  Collapse records
       **
       **  \param[in]: token         - PDB/NDB token name
       **  \param[in]: field_no_set  - field number set
       **  \param[in]: delimiter     - delimiter between concantenated strings
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _collapseRecord(const std::string& token, const std::set<int>& field_no_set, const std::string& delimiter);

       void _checkRecord(const std::string& token);

       /**
       **  Convert REMARK/PDBRMK records into remark number (key) / context (value) pairs
       **
       **  \param[in]: token    - REMARK/PDBRMK tokens
       **  \param[in]: field_no - field number
       **  \param[out]: Remarks - number (key) / context (value) pairs
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_remark_map(const std::string& token, const int& field_no, std::map<int, std::vector<std::string> >& Remarks);

       /**
       **  Add new REMARK to _pdb_records
       **
       **  \param[in]: remark_no - REMARK number
       **  \param[in]: value     - REMARK context
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _addNewRemark(const int& remark_no, const std::string& value);

       /**
       **  Process PDB HEADER record
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
       void _pdb_to_ndb_process_header();

       /**
       **  Process PDB COMPND and SOURCE records
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
       void _pdb_to_ndb_process_compnd_and_source();

       /**
       **  Process PDB COMPND records
       **
       **  \param[in]: rlist        - COMPND records
       **  \param[out]: mol_compnds - COMPND information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_compnd(const std::list<std::vector<std::string> >& rlist, std::map<int, std::map<std::string, std::string> >& mol_compnds);

       /**
       **  Process PDB SOURCE records
       **
       **  \param[in]: rlist        - SOURCE records
       **  \param[out]: mol_sources - SOURCE information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_source(const std::list<std::vector<std::string> >& rlist, std::map<int, std::map<std::string, std::string> >& mol_sources);

       /**
       **  Process PDB record value
       **
       **  \param[out]: value - PDB record value
       **  \param[in]: remove_flag - flag to indicate if ';' need to be removed
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_record_value(std::string& value, const bool& remove_flag);

       /**
       **  Update _entities with COMPND & SOURCE information
       **
       **  \param[in]: mol_compnds - COMPND information
       **  \param[in]: mol_sources - SOURCE information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_update_entities(const std::map<int, std::map<std::string, std::string> >& mol_compnds,
                                        const std::map<int, std::map<std::string, std::string> >& mol_sources); 
       /**
       **  Update _entities with metadata & source information
       **
       **  \param[in]: entity_set - entity ID set
       **  \param[in]: metadata   - meta data information for Entity
       **  \param[in]: source     - source information for Entity
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_update_entities(const std::set<int>& entity_set, const std::map<std::string, std::string>& metadata,
                                        const std::map<std::string, std::string>& source);

       /**
       **  Process PDB HELIX records
       **
       **  \param[in]: records - HELIX record list
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_read_HELIX(const std::list<std::vector<std::string> >& records);

       /**
       **  Process PDB SHEET records
       **
       **  \param[in]: records - SHEET record list
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_read_SHEET(const std::list<std::vector<std::string> >& records);

       /**
       **  Process PDB TURN records
       **
       **  \param[in]: records - TURN record list
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_read_TURN(const std::list<std::vector<std::string> >& records);

       /**
       **  Process PDB SLTBRG records
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
       void _pdb_to_ndb_read_SLTBRG();

       /**
       **  Process PDB MODRES records
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
       void _pdb_to_ndb_read_MODRES();

       /**
       **  Process PDB SITE records
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
       void _pdb_to_ndb_read_SITE();

       /**
       **  Process PDB REMARK records
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
       void _pdb_to_ndb_process_remarks();

       /**
       **  Process REMARK 3
       **
       **  \param[in]: remarkCards - remark context array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_remark_3(const std::vector<std::string>& remarkCards);

       /**
       **  Get REMARK 3 refinement program template index
       **
       **  \param: Not applicable
       **
       **  \return index number
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _pdb_to_ndb_get_refinement_index(const std::string& prog_name = "");

       std::string _get_known_refinement_name(const std::string& prog_name);

       std::string _select_correct_refinement_program(const std::vector<std::string>& programs);

       /**
       **  Process repreat remarks
       **
       **  \param[in]: remark_no      - remark number
       **  \param[in]: remarkCards    - remark context array
       **  \param[in]: start_block    - start block number
       **  \param[in]: group_template - remark template
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_repeat_remark(const int& remark_no, const std::vector<std::string>& remarkCards, const unsigned int& start_block,
                                              const GROUP_REMARKS& group_template);

       /**
       **  Divide remarks into blocks
       **
       **  \param[in]: remarkCards - remark context array
       **  \param[in]: template_tokens - template token set
       **  \param[out]: blocks     - remark context blocks
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_get_block_remarks(const std::vector<std::string>& remarkCards, const std::set<std::string>& template_tokens,
                                          std::vector<std::vector<std::string> >& blocks);

       /**
       **  Divide remarks into blocks
       **
       **  \param[out]: blocks       - remark context blocks
       **  \param[in]: templete_type - REFMAC or PHENIX template
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_get_block_remarks(std::vector<std::vector<std::string> >& blocks, const int& templete_type);

       /**
       **  Process REMARK 3 group template
       **
       **  \param[in]: group_index - refinement program template index
       **  \param[in]: blocks      - remark context blocks
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_group_remark_3(const int& group_index, const std::vector<std::vector<std::string> >& bRemarkCard);

       /**
       **  Update CRMETH.Experimental_Technique field and other experimental token fields
       **
       **  \param[in]: diffrn_id - input diffrn ID
       **  \param[in]: type - experimental type 
       **  \param[in]: scattering_type - _diffrn_radiation.pdbx_scattering_type value
       **  \param[out]: merge_records - concatenated experimental records
       **
       **  \return output diffrn ID
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _pdb_to_ndb_process_experimental_TOKENs(const int& diffrn_id, const std::string& type, const std::string& scattering_type,
                                                   std::map<std::string, std::list<std::vector<std::string> > >& merge_records);

       /**
       **  Process REMARK 265
       **
       **  \param[in]: remarkCards - remark context array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_remark_265(const std::vector<std::string>& remarkCards);

       /**
       **  Process REMARK 350
       **
       **  \param[in]: remarkCards - remark context array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_remark_350(const std::vector<std::string>& remarkCards);

       /**
       **  Process REMARK 400, 450, 600 & 999
       **
       **  \param[in]: remarkCards - remark context array
       **  \param[in]: field_no    - field number for ENTDTL token
       **  \param[in]: keyword     - keyword for free text remark
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_pdbx_entry_details_remark(const std::vector<std::string>& remarkCards, const int& field_no, const std::string& keyword);

       /**
       **  Process REMARK 650 & 700
       **
       **  \param[in]: remark_no    - remark number
       **  \param[in]: remarkCards - remark context array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_block_remark(const int& remark_no, const std::vector<std::string>& remarkCards);

       /**
       **  Get block remark context
       **
       **  \param[in]: remarkCards   - remark context array
       **  \param[in]: keyword       - keyword for free text remark
       **  \param[out]: block_remark - block remark context
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_get_block_remark(const std::vector<std::string>& remarkCards, const std::string& keyword, std::string& block_remark);

       /**
       **  Process REMARK 465
       **
       **  \param[in]: remarkCards - remark context array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_remark_465(const std::vector<std::string>& remarkCards);

       /**
       **  Process REMARK 900
       **
       **  \param[in]: remarkCards - remark context array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_900(const std::vector<std::string>& remarkCards);

       /**
       **  Update "RFTPOS", "RFMPOS", "RFLPOS", "RFTTHR", "RFMTHR" and "RFLTHR" tokens
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
       void _pdb_to_ndb_update_REFMAC_NCS_GROUP();

       /**
       **  Main procedure for converting PDB REMARK records into NDB tokens records
       **
       **  \param[in]: remark_no    - remark number
       **  \param[in]: remarkRecord - remark context
       **  \param[in]: remark_num   - number of remark template items
       **  \param[in]: remark       - PDB REMARK template
       **  \param[in]: repeat       - repeat flag
       **
       **  \return true  - if remark context matches with template
       **  \return false - otherwise
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _pdb_to_ndb_process_remark(const int& remark_no, const std::vector<std::string>& remarkRecord, const int& remark_num, REMARKS *remark,
                                       const int& repeat);

       /**
       **  Find match index between remark context and template
       **
       **  \param[in]: remark_no    - remark number
       **  \param[in]: remarkRecord - remark context
       **  \param[in]: remark_num   - number of remark template items
       **  \param[in]: remark       - PDB REMARK template
       **  \param[in]: repeat       - repeat flag
       **  \param[out]: index       - match index array
       **
       **  \return count of matched number
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _pdb_to_ndb_get_match_index(const int& remark_no, const std::vector<std::string>& remarkRecord, const int& remark_num, REMARKS *remark,
                                       const int& repeat, std::vector<int>& index);

       /**
       **  Check if context matches template
       **
       **  \param[in]: context - remark context
       **  \param[in]: remark  - PDB REMARK template
       **
       **  \return true  - if remark context matches with template
       **  \return false - otherwise
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _pdb_to_ndb_find_a_match(const std::string& context, const REMARKS& remark);

       /**
       **  Process Phenix shell bin
       **
       **  \param[in]: remarkRecord - remark context
       **  \param[in]: index        - match index array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_process_Phenix_Shell_BIN(const std::vector<std::string>& remarkRecord, const std::vector<int>& index);

       /**
       **  Get PDB/NDB token name/context key/value pairs
       **
       **  \param[in]: remarkRecord        - remark context
       **  \param[in]: remark              - PDB REMARK template
       **  \param[in]: is_REFMAC_NCS_GROUP - REFMAC_NCS_GROUP template flag
       **  \param[in]: repeat              - repeat flag
       **  \param[out]: token_Index        - repeat NDB token set
       **  \param[out]: value_pairs        - PDB/NDB token name/context key/value pairs
       **  \param[out]: token              - correspondence PDB/NDB record token name
       **  \param[out]: field_no           - record field number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_get_value_pairs(const std::string& remarkRecord, const REMARKS& remark, const bool& is_REFMAC_NCS_GROUP, const int& repeat,
                                        const std::string& block_remark, std::set<std::string>& token_Index, std::map<std::string,
                                        std::map<int, std::string> >& value_pairs);

       /**
       **  Get context value for correspondence PDB/NDB record token
       **
       **  \param[out]: answer     - output value
       **  \param[in]: remarkCard  - remark context
       **  \param[in]: remark      - PDB REMARK template
       **  \param[in]: field_no    - record field number
       **  \param[in]: moving_flag - column moving flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_get_remark_answer(std::string& answer, const std::string& remarkCard, const REMARKS& remark, const int& field_no,
                                          const bool& moving_flag);

       /**
       **  Insert value into  PDB/NDB token name/context key/value pairs
       **
       **  \param[out]: value_pairs - PDB/NDB token name/context key/value pairs
       **  \param[in]: token        - correspondence PDB/NDB record token name
       **  \param[in]: field_no     - record field number
       **  \param[in]: value        - context value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_insert_value_pairs(std::map<std::string, std::map<int, std::string> >& value_pairs, const std::string& token,
                                           const int& field_no, const std::string& value);

       /**
       **  Insert values to front of _pdb_records
       **
       **  \param[in]: token       - correspondence PDB/NDB record token name
       **  \param[in]: value_pairs - field_no/context key/value pairs
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_insert_front(const std::string& token, const std::map<int, std::string>& value_pairs);

       /**
       **  Insert values to back of _pdb_records
       **
       **  \param[in]: token       - correspondence PDB/NDB record token name
       **  \param[in]: value_pairs - field_no/context key/value pairs
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_insert_back(const std::string& token, const std::map<int, std::string>& value_pairs);

       /**
       **  Insert value to NDB token with RefineID
       **
       **  \param[in]: Multiple_Treatment - Multiple treatment field no
       **  \param[in]: RefineFieldNo - RefineID field number
       **  \param[in]: token         - correspondence PDB/NDB record token name
       **  \param[in]: value_pairs   - field_no/context key/value pairs
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_insert_with_Multiple_Treatment(const int& Multiple_Treatment, const int& RefineFieldNo, const std::string& token,
                                                       const std::map<int, std::string>& value_pairs);

       /**
       **  Insert value to NDB token with RefineID
       **
       **  \param[in]: RefineFieldNo - RefineID field number
       **  \param[in]: token         - correspondence PDB/NDB record token name
       **  \param[in]: value_pairs   - field_no/context key/value pairs
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _pdb_to_ndb_insert_with_RefineField(const int& RefineFieldNo, const std::string& token, const std::map<int, std::string>& value_pairs);

       /**
       **  Update SFTWAR token
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
       void _pdb_to_ndb_update_SFTWAR();

       /**
       **  Update NMRSDT & NMREXP tokens
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
       void _pdb_to_ndb_update_NMRSDT_and_NMREXP();

       /**
       **  Update NMRSPM token
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
       void _pdb_to_ndb_update_NMRSPM();

       /**
       **  Update NMRSMP token
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
       void _pdb_to_ndb_update_NMRSMP();

       /**
       **  Update EMSFTW token
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
       void _pdb_to_ndb_update_EMSFTW();

       /**
       **  Update EMPIXL token
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
       void _pdb_to_ndb_update_EMPIXL();

       /**
       **  Post process various PDB/NDB tokens  
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
       void _pdb_to_ndb_postprocessing();

       /**
       **  Update EXPDTA token
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
       void _pdb_to_ndb_update_EXPDTA();

       /**
       **  Update NCSTLS, TMLPTL, CCPNCS tokens for REFMAC local NCS statistics
       **
       **  \param: Not applicable
       **
       **  \return true if it has local NCS statistics
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _pdb_to_ndb_update_local_NCS();

       /**
       **  Update NCSGRO token
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
       void _pdb_to_ndb_update_NCSGRO();

       /**
       **  Update CCPNCS token
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
       void _pdb_to_ndb_update_CCPNCS();

       /**
       **  Update PHENCS token
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
       void _pdb_to_ndb_update_PHENCS();

       /**
       **  Update TLSGRO & TLSRNG tokens
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
       void _pdb_to_ndb_update_TLSGRO_and_TLSRNG();

       /**
       **  Update RFACTR token
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
       void _pdb_to_ndb_update_RFACTR();

       /**
       **  Update RFACTR & FREERF tokens for SHELX program
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
       void _pdb_to_ndb_update_RFACTR_and_FREERF_for_SHELX();

       /**
       **  Update RFACTR token for TNT, NUCLSQ and PROLSQ programs
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
       void _pdb_to_ndb_update_RFACTR_for_TNT_NUCLSQ_PROLSQ();

       /**
       **  Align multiple NDB tokens for REMARK 200
       **
       **  \param[in]: start_id - input diffrn ID
       **
       **  \return output diffrn ID
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _pdb_to_ndb_update_200_align(const int& start_id);

       /**
       **  Update WAVLEN token
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
       void _pdb_to_ndb_update_WAVLEN();

       /**
       **  Update DTMETH, DTMEAS & RADIAT tokens
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
       void _pdb_to_ndb_update_DTMETH_DTMEAS_RADIAT();

       /**
       **  Update EMDATA token
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
       void _pdb_to_ndb_update_EMDATA();

       /**
       **  Update TLSGRO token
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
       void _pdb_to_ndb_update_TLSGRO();

       /**
       **  Update TWIN token
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
       void _pdb_to_ndb_update_TWIN();

       /**
       **  Update SOLMOD token
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
       void _pdb_to_ndb_update_SOLMOD();

       /**
       **  Update TLSGRO token
       **
       **  \param: Not applicable
       **
       **  \return number of TLSGRO records
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _ndb_to_pdb_update_TLSGRO();

       /**
       **  Update REFMAC 5 program related NDB tokens
       **
       **  \param: Not applicable
       **
       **  \return number of NCSGRO records
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _ndb_to_pdb_processing_refmac5();

       /**
       **  Update PHENIX/BUSTER programs related NDB tokens
       **
       **  \param[in]: is_buster - BUSTER program flag
       **
       **  \return number of NCSGRO records
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _ndb_to_pdb_processing_phenix(const bool& is_buster);

       /**
       **  Update SCTYPE token
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
       void _ndb_to_pdb_update_SCTYPE();

       /**
       **  Update CAVEAT token
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
       void _ndb_to_pdb_update_CAVEAT();

       /**
       **  Update EXPDTA token
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
       void _ndb_to_pdb_update_EXPDTA();

       /**
       **  Update HEADER token
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
       void _ndb_to_pdb_proc_header();

       /**
       **  Update AUTH & AUTHOR tokens
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
       void _ndb_to_pdb_processing_authors();

       /**
       **  Update RNOCUT token
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
       void _ndb_to_pdb_processing_shelx_and_tnt();

       /**
       **  Update diffration source related NDB tokens
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
       void _ndb_to_pdb_processing_synchrotron();

       /**
       **  Update DTMETH, DTMEAS & RADIAT tokens
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
       void _ndb_to_pdb_update_DTMETH_DTMEAS_RADIAT();

       /**
       **  Update WAVLEN token
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
       void _ndb_to_pdb_update_WAVLEN();

       /**
       **  Get synchrotron site & beamline
       **
       **  \param[in]: key_value - key value
       **  \param[out]: site     - synchrotron site
       **  \param[out]: beamline - beamline ID
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_beamline(const std::string& key_value, std::string& site, std::string& beamline);

       /**
       **  Get diffrn ID set
       **
       **  \param[out]: diffrn_ids - diffrn ID set
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_diffrn_ids(std::set<int>& diffrn_ids);

       /**
       **  Get diffrn ID vs. value mapping
       **
       **  \param[in]: token - correspondence PDB/NDB record token name
       **  \param[in]: key_field - key item field number
       **  \param[in]: value_field - value item field number
       **  \param[out]: diffrn_value - ID vs. value mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_diffrn_value(const std::string& token, const int& key_field, const int& value_field, std::map<int, std::string>& diffrn_value);

       /**
       **  Get diffraction type vs. diffrn IDs mapping
       **
       **  \param[in]: diffrn_ids         - diffrn ID set
       **  \param[in]: diffrn_sctype      - ID vs. type mapping
       **  \param[out]: sctype_diffrn_ids - type vs. ID set mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_sctype_diffrn_ids(const std::set<int>& diffrn_ids, const std::map<int, std::string>& diffrn_sctype,
                                              std::map<std::string, std::set<int> >& sctype_diffrn_ids);

       /**
       **  Get diffrn ID vs. crystal ID mapping
       **
       **  \param[out]: diffrn_crystal_ids - diffrn ID vs. crystal ID mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_diffrn_crystal_ids(std::map<int, int>& diffrn_crystal_ids);

       /**
       **  Align multiple NDB tokens for REMARK 200
       **
       **  \param[in]: sctype_diffrn_ids - type vs. ID set mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_update_200_align(const std::map<std::string, std::set<int> >& sctype_diffrn_ids);

       /**
       **  Update COMPND & SOURCE tokens
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
       void _ndb_to_pdb_processing_compnd_source();

       /**
       **  Update COMPND token
       **
       **  \param[in]: Mol_ID - Molecule Identifier
       **  \param[in]: type   - source type
       **  \param[in]: entity - entity information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_processing_compnd(const int& Mol_ID, const std::string& type, const Entity& entity);

       /**
       **  Update SOURCE token
       **
       **  \param[in]: Mol_ID - Molecule Identifier
       **  \param[in]: type   - source type
       **  \param[in]: entity - entity information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_processing_source(const int& Mol_ID, const std::string& type, const Entity& entity);

       /**
       **  Insert COMPND/SOURCE records
       **
       **  \param[in]: token - COMPND/SOURCE
       **  \param[in]: data  - value list
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_add_compnd_or_source(const std::string& token, const std::vector<std::string>& data);

       /**
       **  Update REMARK
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
       void _ndb_to_pdb_get_remarks();

       /**
       **  Generate PDB REMARK
       **
       **  \param[out]: RemarkNo - skip REMARK numbers
       **
       **  \return unprocessed theoretical model flag
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _ndb_to_pdb_get_skip_remark_number(std::set<int>& RemarkNo);

       /**
       **  Update REFMET & NMRSFT
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
       void _ndb_to_pdb_processing_nmr_software();

       /**
       **  Update NMRECD, NMRSPM, NMRSMP & NMRSDT
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
       void _ndb_to_pdb_processing_nmr_experiment();

       /**
       **  Update CMRFMT
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
       void _ndb_to_pdb_processing_em_software();

       /**
       **  Update CMRNST
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
       void _ndb_to_pdb_processing_em_pixel();

       /**
       **  Get EBI ID & version V.3.15 flag
       **
       **  \param[out]: ebi_id       - EBI ID
       **  \param[out]: version_flag - PDB format version V.3.15 flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_auxiliary_data(std::string& ebi_id, bool& version_flag);

       /**
       **  Get user defined REMARK
       **
       **  \param[in]: remark_no - REMARK number
       **  \param[in]: remarks   - REMARK context
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_user_defined_remark(const int& remark_no, const std::vector<std::string>& remarks);

       /**
       **  Update RFACTR
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
       void _ndb_to_pdb_processing_RFACTR();

       /**
       **  Generate general REMARKs based on REMARK template
       **
       **  \param[in]: num_remarks - number of Remarks
       **  \param[in]: Remarks     - pointer to REMARK template array
       **  \param[in]: repeat      - repeat flag
       **  \param[in]: index       - key index value
       **  \param[in]: empty_flag  - empty REMARK flag 
       **  \param[in]: sctype      - diffraction scattering type
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_general_remark(const int& num_remarks, REMARKS *Remarks, const int& repeat, const int& index, const int& empty_flag,
                                           const std::string& sctype = "");

       /**
       **  Update REMARK template field width
       **
       **  \param[in]: num_remarks - number of Remarks
       **  \param[in]: Remarks     - pointer to REMARK template array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_update_template_FieldWidth(const int& num_remarks, REMARKS *Remarks);

       /**
       **  Generate REFMAC TWIN REMARK
       **
       **  \param[in]: num_remarks - number of Remarks
       **  \param[in]: Remarks     - pointer to REMARK template array
       **  \param[in]: sctype      - diffraction scattering type
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_REFMAC_TWIN_remark(const int& num_remarks, REMARKS *Remarks, const std::string& sctype);

       /**
       **  Generate general repeat REMARK
       **
       **  \param[in]: card_id     - PDB token name
       **  \param[in]: num_remarks - number of Remarks
       **  \param[in]: Remarks     - pointer to REMARK template array
       **  \param[in]: empty_flag  - empty REMARK flag 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_repeat_remark(const std::string& card_id, const int& num_remarks, REMARKS *Remarks, const int& empty_flag);

       /**
       **  Generate general repeat REMARK
       **
       **  \param[in]: num_remarks - number of Remarks
       **  \param[in]: Remarks     - pointer to REMARK template array
       **  \param[in]: index       - key index value
       **  \param[in]: sctype      - diffraction scattering type
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_non_repeat_remark(const int& num_remarks, REMARKS *Remarks, const int& index, const std::string& sctype);

       /**
       **  Generate Refmac local NCS REMARK
       **
       **  \param[in]: Remark_No - Remark number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_REFMAC_NCS_LOCAL(const int& Remark_No);

       /**
       **  Generate Phenix BIN REMARK
       **
       **  \param[in]: Remark_No - Remark number
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_Phenix_BIN_remark(const int& Remark_No);

       /**
       **  Generate Phenix NCS REMARK
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
       void _ndb_to_pdb_get_PHENIX_NCS_remark();

       /**
       **  Generate general REMARKs based on REMARK template
       **
       **  \param[in]: Remark - REMARK template
       **  \param[in]: sctype - diffraction scattering type
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_general_remark(const REMARKS& Remark, const std::string& sctype);

       /**
       **  Generate general REMARKs based on REMARK template
       **
       **  \param[in]: Remark     - REMARK template
       **  \param[in]: sctype     - diffraction scattering type
       **  \param[in]: diffrn_id  - diffraction ids
       **  \param[in]: pCardField - PDB/NDB record
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_general_remark(const REMARKS& Remark, const std::string& sctype, const std::string& diffrn_id,
                                           const std::vector<std::string>& pCardField);

       /**
       **  Get Remark field value from related NDB token
       **
       **  \param[in]: FieldList  - REMARK Field template
       **  \param[in]: Remark_No  - Remark number
       **  \param[in]: sctype     - diffraction scattering type
       **  \param[in]: diffrn_id  - diffraction ids
       **  \param[in]: last_field - last field flag
       **
       **  \return field context
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _ndb_to_pdb_get_field_value(const REMARK_FIELD& FieldList, const int& Remark_No, const std::string& sctype, 
                                               const std::string& diffrn_id, const bool& last_field);

       /**
       **  Re-format numeric value based on REMARK template
       **
       **  \param[out]: value    - field context
       **  \param[in]: width     - Field width
       **  \param[in]: FieldList - REMARK Field template
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_reformat_numeric_value(std::string& value, const int& width, const REMARK_FIELD& FieldList);

       /**
       **  Format REMARK context
       **
       **  \param[in]: Remark_No  - Remark number
       **  \param[in]: FieldList0 - first REMARK Field template
       **  \param[in]: FieldList  - cirrent REMARK Field template
       **  \param[in]: value      - field context
       **  \param[out]: remark_string - final REMARK context
       **  \param[out]: remark_array  - final REMARK context array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_format_field_value(const int& Remark_No, const REMARK_FIELD& FieldList0, const REMARK_FIELD& FieldList, std::string& value,
                                           std::string& remark_string, std::vector<std::string>& remark_array);

       /**
       **  Generate REAMRK 3
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
       void _ndb_to_pdb_get_group_remark_3(const std::string& refinement_program);
       void _ndb_to_pdb_get_remark_3_top();
       void _ndb_to_pdb_get_remark_3_update_Xplor();
       void _ndb_to_pdb_update_joint_methods();
       void _ndb_to_pdb_get_remark_3(GROUP_REMARKS *group_remark);


       /**
       **  Generate REAMRK 100
       **
       **  \param[in]: is_unprocessed_model - unprocessed model flag
       **  \param[in]: EBI_ID - EBI ID
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_100(const bool& is_unprocessed_model, const std::string& EBI_ID);

       /**
       **  Generate REAMRK 265
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
       void _ndb_to_pdb_get_remark_265();

       /**
       **  Generate REAMRK 280
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
       void _ndb_to_pdb_get_remark_280();

       /**
       **  Generate REAMRK 290
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
       void _ndb_to_pdb_get_remark_290();

       /**
       **  Generate REAMRK 300
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
       void _ndb_to_pdb_get_remark_300();

       /**
       **  Get author defined REMARk 300 details
       **
       **  \param: Not applicable
       **
       **  \return detail context
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _ndb_to_pdb_get_struct_biol_details();

       /**
       **  Generate SPLIT entry's REAMRK 300
       **
       **  \param: Not applicable
       **
       **  \return true - if it's SPLIT entry
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _ndb_to_pdb_get_remark_300_split();

       /**
       **  Generate author's detail descriptions REMARK 300
       **
       **  \param[in]: details - detail descriptions
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_author_defined_remark_300_details(const std::string& details);

       /**
       **  Generate REMARK from context array
       **
       **  \param[in]: Remark_No    - REMARK number
       **  \param[in]: remark_array - remark_array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_add_remark(const int& Remark_No, const std::vector<std::string>& remark_array);

       /**
       **  Get Biological Molecule IDs
       **
       **  \param: Not applicable
       **
       **  \return Biol_IDs context
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _ndb_to_pdb_get_biol_ids();

       /**
       **  Get point symmetry related remark context
       **
       **  \param[out]: remark_array - remark context array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_300_point_symmetry(std::vector<std::string>& remark_array);

       /**
       **  Get helical symmetry related remark context
       **
       **  \param[out]: remark_array - remark context array
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_300_helical_symmetry(std::vector<std::string>& remark_array);

       /**
       **  Generate REAMRK 350
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
       void _ndb_to_pdb_get_remark_350();

       /**
       **  Generate SPLIT entry's REAMRK 350
       **
       **  \param: Not applicable
       **
       **  \return true - if it's SPLIT entry
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _ndb_to_pdb_get_remark_350_split();

       /**
       **  Generate REAMRK 350's PDB chain ID list
       **
       **  \param[in]: id_list - PDB chain ID list
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_350_APPLY_CHAINS(const std::vector<std::string> &id_list);

       /**
       **  Generate REAMRK 350's symmetry matrix
       **
       **  \param[in]: serialNo - serial number
       **  \param[in]: mat      - matrix value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_350_Remark_SYMMA(const int& serialNo, double mat[4][4]);

       /**
       **  Generate REAMRK 375
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
       void _ndb_to_pdb_get_remark_375();

       /**
       **  Generate REAMRK 375
       **
       **  \param: special_position_atoms - special position atom list
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_375(const std::list<Value>& special_position_atoms);

       /**
       **  Generate REAMRKs 400, 450, 600 & 999
       **
       **  \param[in]: remark_no - REMARK number
       **  \param[in]: field_no  - field number for ENTDTL token
       **  \param[in]: keyword   - keyword for free text remark
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark(const int& remark_no, const int& field_no, const std::string& keyword);

       /**
       **  Generate PRD group REAMRK 400
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
       void _ndb_to_pdb_get_remark_400_group();

       void _ndb_to_pdb_get_polymer_chain_info(std::map<std::string, int>& chain_order, std::map<std::string, RCSB::Chain*>& asym_id_chain);

       void _ndb_to_pdb_get_non_polymer_prd_info(ISTable *instances, ISTable *features, const std::map<std::string, RCSB::Chain*>& asym_id_chain,
                                                 const std::set<std::string>& polymer_prds);

       /**
       **  Generate REAMRK 465
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
       void _ndb_to_pdb_get_remark_465();

       /**
       **  Generate REAMRK 470
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
       void _ndb_to_pdb_get_remark_470();

       /**
       **  Generate REAMRK 475
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
       void _ndb_to_pdb_get_remark_475();

       /**
       **  Generate REAMRK 480
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
       void _ndb_to_pdb_get_remark_480();

       bool _check_new_nmr_template(const std::map<int, std::list<std::vector<std::string> > >& r_maps);

       void _ndb_to_pdb_get_remark_465_475(const std::string& mol_id, const int& remark_no, const std::list<std::vector<std::string> >& l_list);

       void _ndb_to_pdb_get_remark_470_480(const std::string& mol_id, const int& remark_no, const std::list<std::vector<std::string> >& l_list);

       /**
       **  Generate REAMRK 500
       **
       **  \param: Not applicable
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500();

       /**
       **  Generate REAMRK 500 close contact
       **
       **  \param: Not applicable
       **  \param[in]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500_asym_contact();
       void _ndb_to_pdb_get_remark_500_asym_contact(Block& block);

       /**
       **  Generate REAMRK 500 symmetry contact
       **
       **  \param: Not applicable
       **  \param[in]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500_symm_contact();
       void _ndb_to_pdb_get_remark_500_symm_contact(Block& block);

       /**
       **  Generate REAMRK 500 contacts
       **
       **  \param[in]: contacts - asym/sym contacts
       **  \param[in]: sym_flag - asym/sym flag
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500_close_contact(const std::list<Value>& contacts, const bool& sym_flag);
       void _ndb_to_pdb_get_remark_500_close_contact(ISTable *t, const bool& sym_flag);
       void _remark_500_close_contact(const std::vector<RCSB::Atom*>& atoms, const double& val, const std::string& symmetry);

       /**
       **  Generate REAMRK 500 bond rms deviation
       **
       **  \param: Not applicable
       **  \param[in]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500_bond_deviation();
       void _ndb_to_pdb_get_remark_500_bond_deviation(Block& block);
       void _remark_500_bond_deviation(const std::vector<RCSB::Atom*>& atoms, const double& val, const int& mol_id);

       /**
       **  Generate REAMRK 500 angle rms deviation
       **
       **  \param: Not applicable
       **  \param[in]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500_angle_deviation();
       void _ndb_to_pdb_get_remark_500_angle_deviation(Block& block);
       void _remark_500_angle_deviation(const std::vector<RCSB::Atom*>& atoms, const double& val, const int& mol_id);

       /**
       **  Generate REAMRK 500 Ramachandran outliers
       **
       **  \param: Not applicable
       **  \param[in]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500_Ramachandran_outliers();
       void _ndb_to_pdb_get_remark_500_Ramachandran_outliers(Block& block);
       void _remark_500_Ramachandran_outliers(RCSB::Atom* atom, const double& psi, const double& phi, const int& mol_id);

       /**
       **  Generate REAMRK 500 non cis and non trans torsion
       **
       **  \param: Not applicable
       **  \param[in]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500_non_cis_trans_torsions();
       void _ndb_to_pdb_get_remark_500_non_cis_trans_torsions(Block& block);
       void _remark_500_non_cis_trans_torsions(const std::vector<RCSB::Atom*>& atoms, const double& val, const int& mol_id);

       /**
       **  Generate REAMRK 500 side chain plane deviation
       **
       **  \param: Not applicable
       **  \param[in]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500_side_chain_plane();
       void _ndb_to_pdb_get_remark_500_side_chain_plane(Block& block);
       void _remark_500_side_chain_plane(RCSB::Atom* atom, const double& val, const int& mol_id, const std::string& details);

       /**
       **  Generate REAMRK 500 main chain planarity deviation
       **
       **  \param: Not applicable
       **  \param[in]: block - reference to a coordinate data block
       **
       **  \return Not applicable
       **
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_get_remark_500_main_chain_planarity();
       void _ndb_to_pdb_get_remark_500_main_chain_planarity(Block& block);
       void _remark_500_main_chain_planarity(RCSB::Atom* atom, const double& val, const int& mol_id);

       /**
       **  Generate REAMRK 525
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
       void _ndb_to_pdb_get_remark_525();
       void _remark_525(RCSB::Atom* atom, const double& val);

       /**
       **  Generate REAMRK 610
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
       void _ndb_to_pdb_get_remark_610();

       /**
       **  Generate REAMRK 615
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
       void _ndb_to_pdb_get_remark_615();

       /**
       **  Generate REAMRK 620
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
       void _ndb_to_pdb_get_remark_620();

       /**
       **  Generate REAMRK 630
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
       void _ndb_to_pdb_get_remark_630();

       /**
       **  Generate REAMRK 900
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
       void _ndb_to_pdb_get_remark_900();

       void _get_missing_or_zero_occupancy_residues_or_atoms_for_remarks();
       void _analysis_missing_or_zero_occupancy_residues();
       void _analysis_missing_or_zero_occupancy_atoms();
       /**
       **  Update MDLTYP record
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
       void _ndb_to_pdb_processing_Ca_and_P_atom_only();

       /**
       **  Update REF token
       **
       **  \param[out]: id_abbrev_mapping - Jrnl ID vs. abbrev. mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_update_REF(std::map<std::string, std::string>& id_abbrev_mapping);

       /**
       **  Update REFN token
       **
       **  \param[in]: id_abbrev_mapping - Jrnl ID vs. abbrev. mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _ndb_to_pdb_update_REFN(const std::map<std::string, std::string>& id_abbrev_mapping);

       /**
       **  Update PMID token
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
       void _ndb_to_pdb_update_PMID();

       /**
       **  Update DBREF token
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
       void _ndb_to_pdb_update_DBREF();

       /**
       **  Update SEQRES token
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
       void _ndb_to_pdb_update_SEQRES();

       /**
       **  Update HELIX token
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
       void _ndb_to_pdb_update_HELIX();

       /**
       **  Update SHEET token
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
       void _ndb_to_pdb_update_SHEET();

       /**
       **  Update TURN token
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
       void _ndb_to_pdb_update_TURN();

       /**
       **  Update SSBOND token
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
       void _ndb_to_pdb_update_SSBOND();

       /**
       **  Update LINK token
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
       void _ndb_to_pdb_update_LINK();

       /**
       **  Update SLTBRG token
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
       void _ndb_to_pdb_update_SLTBRG();

       /**
       **  Update CISPEP token
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
       void _ndb_to_pdb_update_CISPEP();

       /**
       **  Update HET, HETNAM, HETSYN & FORMUL tokens
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
       void _ndb_to_pdb_get_HET_info();

       /**
       **  Update MODRES token
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
       void _ndb_to_pdb_update_MODRES();

       /**
       **  Update SITE & SITELS tokens
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
       void _ndb_to_pdb_update_SITE();

       /**
       **  Add CRYST1 record if not exist
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
       void _ndb_to_pdb_update_CRYST1();

       /**
       **  Update ORIGX token
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
       void _ndb_to_pdb_update_ORIGX();

       void _ndb_to_pdb_update_OBSLTE();

       void _ndb_to_pdb_get_remark_5();

       void _ndb_to_pdb_get_atom_connects(const bool& checking_linkage_flag, std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& atom_connects);
};

#endif
