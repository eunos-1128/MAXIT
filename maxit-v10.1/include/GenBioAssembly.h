/*
FILE:     GenBioAssembly.h
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
#ifndef _H_GEN_BIOL_ASSEMBLY_H_
#define _H_GEN_BIOL_ASSEMBLY_H_

#include <string>
#include <vector>

#include "Maxit.h"

class GenBioAssembly: public Maxit {
   public:

       /**
       **  Constructs GenBioAssembly
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
       GenBioAssembly();

       /**
       **  Destructs GenBioAssembly
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
       ~GenBioAssembly();

       /**
       **  Write biological pdb files
       **
       **  \param[in]: rootname  - output file root name
       **  \param[in]: inputfile - input file name
       **  \param[in]: indexfile - output index file name
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void write_biological_pdb_files(const std::string& rootname, const std::string& inputfile, const std::string& indexfile);

       /**
       **  Write biological cif files
       **
       **  \param[in]: rootname  - output file root name
       **  \param[in]: indexfile - output index file name
       **  \param[in]: public_flag - flag to indicate if it is for public FTP site
       **  \param[in]: single_model_flag - flag to indicate if the symmetry-related chains merged into single model
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void write_biological_cif_files(const std::string& rootname, const std::string& indexfile, const bool& public_flag, const bool& single_model_flag);

   private:
       bool _has_symmetry_related_operation;
       bool _remove_non_public_item_flag;

       /**
       **  Get original chain list
       **
       **  \param[in]:  pdb_chains      - input file name
       **  \param[out]: original_chains - output file root name
       **  \param[out]: chain_id_set - polymer PDB chain ID set
       **  \param[out]: entity_id_set - entity ID set
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_original_chains(const std::vector<std::string>& pdb_chains, std::vector<RCSB::Chain*>& original_chains,
                                 std::set<std::string>& chain_id_set, std::set<std::string>& entity_id_set);

       /**
       **  Get biological unit molecule
       **
       **  \param[in]: original_chains - original chain list
       **  \param[in]: mat  - 4x4 transformation matrix
       **
       **  \return Moulecule pointer
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       RCSB::Molecule* _get_bio_molecule(const std::vector<RCSB::Chain*>& original_chains, double mat[4][4]);

       /**
       **  Add chains to biological unit molecule
       **
       **  \param[in]: single_model_flag - single model flag
       **  \param[in]: original_chains - original chain list
       **  \param[in]: operators - symmetric operators
       **  \param[in]: mat  - 4x4 transformation matrix
       **  \param[out]: mol - moulecule object pointer
       **  \param[out]: chain_index_list - chain index set list
       **  \param[out]: original_asym_chain_id_mapping - original asym ID & PDB chain ID mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_assembly_chains(const bool& single_model_flag, const std::vector<RCSB::Chain*>& original_chains, const std::string& operators,
                                    double mat[4][4], RCSB::Molecule* mol, std::vector<std::set<int> >& chain_index_list,
                                    std::map<int, std::vector<std::string> >& original_asym_chain_id_mapping);

       /**
       **  Update atom_site & atom_site_anisotrop tables
       **
       **  \param[in]: mol_id - moulecule ID number
       **  \param[in]: mol - moulecule object pointer
       **  \param[out]: atomSerialNo - atom serial number
       **  \param[out]: atomTable - atom_site table
       **  \param[out]: anisotropTable - atom_site_anisotrop table
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_atom_site_categories(const int& mol_id, RCSB::Molecule* mol, int& atomSerialNo, ISTable* atomTable, ISTable* anisotropTable);

       /**
       **  Get selected entity_poly category's value(s)
       **
       **  \param[in]:  block - Reference to data block
       **  \param[in]:  chain_id_set - selected chain ID set
       **  \param[out]: items - item names for the category
       **  \param[out]: values - selected values for the category
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_entity_poly_values(Block& block, const std::set<std::string>& chain_id_set, std::vector<std::string>& items,
                                    std::vector<std::map<std::string, std::string> >& values);

       /**
       **  Get selected category's value(s)
       **
       **  \param[in]:  block - Reference to data block
       **  \param[in]:  category - category name
       **  \param[in]:  key_item - key item name
       **  \param[in]:  entity_id_set - selected entity ID set
       **  \param[out]: items - item names for the category
       **  \param[out]: values - selected values for the category
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_entity_values(Block& block, const std::string& category, const std::string& key_item, const std::set<std::string>& entity_id_set,
                               std::vector<std::string>& items, std::vector<std::map<std::string, std::string> >& values);

       /**
       **  Update category values
       **
       **  \param[out]: block - Reference to data block
       **  \param[in]:  category - category name
       **  \param[in]:  items - item names for the category
       **  \param[in]:  values - values for the category
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_values(Block& block, const std::string& category, const std::vector<std::string>& items,
                           const std::vector<std::map<std::string, std::string> >& values);

       /**
       **  Add coordinate related meta categories
       **
       **  \param[in]: modelBlock - Reference to data block of model file
       **  \param[in]: original_asym_chain_id_mapping - original asym ID & PDB chain ID mapping
       **  \param[out]: biolBlock - Reference to data block of assembly file
       **  \param[out]: mol - moulecule object pointer
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_non_atom_site_categories(Block& modelBlock, const std::map<int, std::vector<std::string> >& original_asym_chain_id_mapping,
                                             Block& biolBlock, RCSB::Molecule* mol);

       /**
       **  Add struct_conf_type & struct_conf categories
       **
       **  \param[in]: modelBlock - Reference to data block of model file
       **  \param[in]: original_asym_chain_id_mapping - original asym ID & PDB chain ID mapping
       **  \param[in]: chain_id_set - polymer chain ID set
       **  \param[in]: chain_index_list - list of chain index set for each symmetric operator
       **  \param[in]: mol - moulecule object pointer
       **  \param[out]: biolBlock - Reference to data block of assembly file
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_struct_conf_categories(Block& modelBlock, const std::map<int, std::vector<std::string> >& original_asym_chain_id_mapping,
                                           const std::set<std::string>& chain_id_set, const std::vector<std::set<int> >& chain_index_list,
                                           RCSB::Molecule* mol, Block& biolBlock);

       /**
       **  Add struct_sheet, struct_sheet_order, struct_sheet_range, & pdbx_struct_sheet_hbond categories
       **
       **  \param[in]: modelBlock - Reference to data block of model file
       **  \param[in]: original_asym_chain_id_mapping - original asym ID & PDB chain ID mapping
       **  \param[in]: chain_id_set - polymer chain ID set
       **  \param[in]: chain_index_list - list of chain index set for each symmetric operator
       **  \param[in]: mol - moulecule object pointer
       **  \param[out]: biolBlock - Reference to data block of assembly file
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_struct_sheet_categories(Block& modelBlock, const std::map<int, std::vector<std::string> >& original_asym_chain_id_mapping,
                                            const std::set<std::string>& chain_id_set, const std::vector<std::set<int> >& chain_index_list,
                                            RCSB::Molecule* mol, Block& biolBlock);

       /**
       **  Update pdbx_entity_mapping category
       **
       **  \param[out]: block - Reference to data block
       **  \param[in]:  entity_id_mapping - new vs. old entity ID mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_entity_id_mapping_category(Block& block, const std::map<int, std::pair<std::string, int> >& entity_id_mapping);

       /**
       **  Update entity category
       **
       **  \param[in]: modelBlock - Reference to data block of model file
       **  \param[out]: biolBlock - Reference to data block of assembly file
       **  \param[in]:  entity_id_mapping - new vs. old entity ID mapping
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_entity_category(Block& modelBlock, Block& biolBlock, const std::map<int, std::pair<std::string, int> >& entity_id_mapping);

       /**
       **  Update polymer entity related categories ( entity_poly, entity_poly_seq etc.)
       **
       **  \param[in]: modelBlock - Reference to data block of model file
       **  \param[out]: biolBlock - Reference to data block of assembly file
       **  \param[in]:  polymer_entity_id_mapping - new vs. old entity ID mapping for polymeric entities
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_polymer_entity_related_categories(Block& modelBlock, Block& biolBlock, const std::map<int, std::pair<std::string,
                                                      std::vector<std::string> > >& polymer_entity_id_mapping);

       /**
       **  Update pdbx_chain_mapping category
       **
       **  \param[out]: block - Reference to data block
       **  \param[in]:  values - values for asym & PDB chain IDs mapping category
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_chain_id_mapping_category(Block& block, const std::vector<std::map<std::string, std::string> >& values);

       /**
       **  Update single key value categories
       **
       **  \param[in]: modelBlock - Reference to data block of model file
       **  \param[out]: biolBlock - Reference to data block of assembly file
       **  \param[in]: category - category name
       **  \param[in]: key_item - key item name
       **  \param[in]: updated_key_value_mapping -
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_single_key_value_categories(Block& modelBlock, Block& biolBlock, const std::string& category, const std::string& key_item, const std::map<int,
                                                std::pair<std::string, std::map<std::string, std::string> > >& updated_key_value_mapping);

       /**
       **  Update multiple key value categories
       **
       **  \param[in]: modelBlock - Reference to data block of model file
       **  \param[out]: biolBlock - Reference to data block of assembly file
       **  \param[in]: category - category name
       **  \param[in]: key_item - key item name
       **  \param[in]: ordinal_item - ordinal item name
       **  \param[in]: updated_key_value_mapping -
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_multiple_key_value_categories(Block& modelBlock, Block& biolBlock, const std::string& category, const std::string& key_item,
                                                  const std::string& ordinal_item, const std::map<int, std::pair<std::string,
                                                  std::map<std::string, std::string> > >& updated_key_value_mapping);

       /**
       **  Get single_key_value_map 
       **
       **  \param[in]: block - Reference to data block
       **  \param[in]: category - category name
       **  \param[in]: key_item - key item name
       **  \param[out]:  item_names - item name list
       **  \param[out]:  single_key_value_map - 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_single_key_value_map(Block& block, const std::string& category, const std::string& key_item, std::vector<std::string>& item_names,
                                      std::map<std::string, std::map<std::string, std::string> >& single_key_value_map);

       /**
       **  Get multiple_key_value_map 
       **
       **  \param[in]: block - Reference to data block
       **  \param[in]: category - category name
       **  \param[in]: key_item - key item name
       **  \param[out]:  item_names - item name list
       **  \param[out]:  multiple_key_value_map - 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_multiple_key_value_map(Block& block, const std::string& category, const std::string& key_item, std::vector<std::string>& item_names,
                                        std::map<std::string, std::vector<std::map<std::string, std::string> > >& multiple_key_value_map);

       /**
       **  Get helix record data from struct_conf category
       **
       **  \param[in]: block - Reference to data block
       **  \param[in]: chain_id_set - polymer chain ID set
       **  \param[in]: key_item - key item name
       **  \param[out]: item_names - item name list
       **  \param[out]: original_values - 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_original_struct_conf_values(Block& block, const std::set<std::string>& chain_id_set, std::vector<std::string>& item_names,
                                             std::map<std::string, std::vector<std::map<std::string, std::string> > >& original_values);

       /**
       **  Get sheet records data from struct_sheet* categories
       **
       **  \param[in]: block - Reference to data block
       **  \param[in]: category - category name
       **  \param[in]: sheet_id_item - sheet ID item name
       **  \param[in]: check_item_names - checking item name list
       **  \param[in]: check_value_set - checking value set
       **  \param[out]: item_names - item name list
       **  \param[out]: original_values - 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_sheet_values(Block& block, const std::string& category, const std::string& sheet_id_item, const std::vector<std::string>& check_item_names,
                              const std::set<std::string>& check_value_set, std::vector<std::string>& category_item_names,
                              std::vector<std::pair<std::string, std::vector<std::map<std::string, std::string> > > >& original_values);

       /**
       **  Check if selected values are valid
       **
       **  \param[in]: check_item_names - checking item name list
       **  \param[in]: check_value_set - checking value set
       **  \param[in]: selected_values - 
       **
       **  \return true if all checked values are defined in check_value_set
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _is_valid_selected_values(const std::vector<std::string>& check_item_names, const std::set<std::string>& check_value_set,
                                      const std::vector<std::map<std::string, std::string> >& selected_values);

       /**
       **  Get updated sheet records from original sheet records
       **
       **  \param[in]: input_values - original sheet records
       **  \param[in]: sheet_id_map - original vs new sheet ID map
       **  \param[in]: chain_id_item_map - auth vs label chain ID item name map
       **  \param[in]: chain_id_map - original vs new chain/asym IDs map
       **  \param[out]: output_values - updated sheet records
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _get_updated_sheet_value(const std::vector<std::pair<std::string, std::vector<std::map<std::string, std::string> > > >& input_values,
                                     const std::map<std::string, std::string>& sheet_id_map, const std::map<std::string, std::string>& chain_id_item_map,
                                     const std::map<std::string, std::pair<std::string, std::string> >& chain_id_map,
                                     std::vector<std::map<std::string, std::string> >& output_values);
};

#endif
