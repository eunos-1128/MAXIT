/*
FILE:     ConnectFormat.h
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

/*!
** \file ConnectFormat.h
**
** \brief Header file for ConnectFormat class.
*/

#ifndef _H_CONNECTFORMAT_H_
#define _H_CONNECTFORMAT_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#include "AtomFormat.h"
#include "TableFile.h"

/**
**  \class ConnectFormat
** 
**  \brief Public class that respresents a chemical component definition.
**   bond information.
** 
**  The class provides methods for construction, destruction, clearance,
**  assignment operator and retrieval.
*/

class ConnectFormat {
   private:
       std::string _empty_string;
       std::vector<std::string> _empty_vector;
       std::vector<std::map<std::string, std::string> > _empty_tree_list;

       // chiral center flag
       bool _has_chiral_center;
       bool _has_backbone_terminal_flag;
       // number of heavy atoms
       int  _nheavyatoms;
       // Meta data from chem_comp category
       std::map<std::string, std::string> _chem_comp_data;
       // atom records from chem_comp_atom category
       std::vector<AtomFormat> _atoms, _order_atoms;
       // atom records from chem_comp_atom category in vector
       std::vector<std::vector<std::string> > _atom_vectors;
       // atom pair list: first - atom_type, second - atom_name
       std::vector<std::pair<std::string, std::string> > _atom_pair_list;
       // leaving atoms
       std::set<std::string> _leaving_atom_set;
       // bond records from chem_comp_bond category
       std::vector<std::vector<std::string> > _bonds;
       // bond pair list: first - first atom, second - second atom
       std::vector<std::pair<std::string, std::string> > _bond_pair_list;
       // linked atom mapping
       std::map<std::string, std::vector<std::string> > _linked_atom_mapping;
       // linked atom mapping using _atoms indices
       std::map<unsigned int, std::vector<unsigned int> > _linked_atom_index_mapping;
       // atom name index for _atoms
       std::map<std::string, unsigned int> _atomOrder;
       std::map<std::string, unsigned int> _stndAtomOrder;

       std::set<std::string> _terminal_connected_atoms;

       // map key: atom type
       // value: terminal connected atom list
       std::map<std::string, std::vector<std::string> > _terminal_connected_atom_type_mapping;

       // map key: category name
       std::map<std::string, ISTable> _category_mapping;

       // map key: _chem_comp_atom.atom_id
       // value: _chem_comp_atom.pdbx_stnd_atom_id
       std::map<std::string, std::string> _atom_id_pdbx_stnd_atom_id_mapping;

       // key: "_pdbx_chem_comp_identifier.type"_"_pdbx_chem_comp_identifier.program" - all uppercase
       // value.first: _pdbx_chem_comp_identifier.program_version
       // value.second: _pdbx_chem_comp_identifier.identifier
       std::map<std::string, std::vector<std::pair<std::string, std::string> > > _chem_comp_identifier;

       // key: "_pdbx_chem_comp_feature.type"_"_pdbx_chem_comp_feature.source" - all uppercasemap 
       // value: _pdbx_chem_comp_feature.value
       std::map<std::string, std::string> _chem_comp_feature;

       // sugar anomeric carbon atom name
       std::string _anomeric_carbon;
       // sugar anomeric oxygen atom name, for modified residue, it may not be the oxygen
       std::string _anomeric_oxygen;
       // sugar ring atom set;
       std::set<std::string> _ring_atoms;

       // key: atom name
       std::map<std::string, int> _sugar_locant_mapping;

       // key: atom name
       // map-key: atom name
       // map-val: atom type
       //
       // vector[0]: first generation
       // vector[1]: second generation
       // ....
       // vector[N-1]: N-th generation
       std::map<std::string, std::vector<std::map<std::string, std::string> > > _sugar_side_atom_tree_mapping;

       /**
       **  Insert bonded atoms into _linked_atom_mapping
       **
       **  \param[in]: first_atom - bonded first atom
       **  \param[in]: second_atom - bonded second atom
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_linked_atom_mapping(const std::string& first_atom, const std::string& second_atom);

       /**
       **  Insert bonded atom indices into _linked_atom_index_mapping
       **
       **  \param[in]: index1 - index for bonded first atom
       **  \param[in]: index2 - index for bonded second atom
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_linked_atom_index_mapping(const unsigned int& index1, const unsigned int& index2);

       /**
       **  Get reordered atoms based on sub component list
       **
       **  \param: None
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _getOrderedAtoms();

       /**
       **  Find break between sub component and numbering each sub component
       **
       **  \param: Not applicable
       **
       **  \return true  - if no sub component or sub components are already numbered
       **          false - if sub components are not numbered 
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _findSubComponentBreak();

       /**
       **  Find break between sub components with same residue names
       **
       **  \param[in]:   begin - range begin
       **  \param[in]:   end   - range end
       **  \param[out]:  residue_breaks - break point set
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _getSameSubComponentBreak(const unsigned int& begin, const unsigned int& end, std::set<unsigned int>& residue_breaks);

       /**
       **  Get standardized bond type
       **
       **  \param: inBondType - input bond type
       **
       **  \return standardized bond type
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _get_bond_type(const std::string& inBondType);

       /**
       **  Retrieves connected terminal atoms with leaving attached to them
       **
       **  \param: Not applicable
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _findTerminalConnectedAtoms();

       /**
       **  Retrieves _chem_comp_identifier vaule by map key
       **
       **  \param[in]: key - "_pdbx_chem_comp_identifier.type"_"_pdbx_chem_comp_identifier.program" value
       **
       **  \return _pdbx_chem_comp_identifier.identifier value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& _get_chem_comp_identifier(const std::string& key) const;

       /**
       **  Retrieves _chem_comp_feature vaule by map key
       **
       **  \param[in]: key - "_pdbx_chem_comp_feature.type"_"_pdbx_chem_comp_feature.source" value
       **
       **  \return _pdbx_chem_comp_feature.value value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& _get_chem_comp_feature(const std::string& key) const;

       /**
       **  Update mono saccharide related information.
       **
       **  \param: Not applicable
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _updateSaccharideInformation();

       /**
       **  Get number part from atom name.
       **
       **  \param[in]: atom_name - atom name
       **
       **  \return the number if exists.
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _get_number_part_from_atom_name(const std::string& atom_name);

       /**
       **  Update mono saccharide side atom related information.
       **
       **  \param[in]: include_set - include atom set
       **  \param[in/out]: exclude_set - exclude atom set
       **  \param[out]: side_atom_tree_list - attached side atom tree list
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _update_side_atom_info(const std::map<std::string, int>& include_map, std::set<std::string>& exclude_set,
                                   std::vector<std::map<std::string, std::string> >& side_atom_tree_list);

   public:
       /**
       **  Constructs ConnectFormat
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
       ConnectFormat() { clear(); }

       /**
       **  Destructs ConnectFormat
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
       ~ConnectFormat() { clear(); }

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
       **  Copy a ConnectFormat to another ConnectFormat (assignment operator).
       **
       **  \param[in]: atom - reference to the source ConnectFormat
       **
       **  \return Reference to the destination ConnectFormat
       **
       **  \pre None
       **
       **  \post Constructed ConnectFormat is a clone as the ConnectFormat
       **        referenced by \e atom
       **
       **  \exception: None
       */
       ConnectFormat& operator=(const ConnectFormat& drug);

       /**
       **  Retrieves chiral center flag
       **
       **  \param: None
       **
       **  \return Constant reference to a bool that indicates if the ConnectFormat
       **             has chiral center(s)
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const bool& has_chiral_center() const;

       /**
       **  Retrieves backbone and/or terminal atom flag
       **
       **  \param: None
       **
       **  \return Constant reference to a bool that indicates if the ConnectFormat
       **             has backbone and/or terminal atom(s)
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const bool& has_backbone_terminal_flag() const;

       /**
       **  Retrieves heavy atom number
       **
       **  \param: None
       **
       **  \return Constant reference to a integer that contains heavy atom number
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const int& nheavyatoms() const;

       /**
       **  Retrieves total atom number
       **
       **  \param: None
       **
       **  \return total atom number
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const int natoms() const;

       /**
       **  Retrieves total bond number
       **
       **  \param: None
       **
       **  \return total bond number
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const int nbonds() const;

       /**
       **  Retrieves chemical component ID
       **
       **  \param: None
       **
       **  \return chemical component ID
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& drugname() const;

       /**
       **  Retrieves chemical formula
       **
       **  \param: None
       **
       **  \return chemical formula
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& formula() const;

       /**
       **  Retrieves chemical name
       **
       **  \param: None
       **
       **  \return chemical name
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& chemical_name() const;

       /**
       **  Retrieves synonym
       **
       **  \param: None
       **
       **  \return synonym
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& synonym() const;

       /**
       **  Retrieves "_pdbx_chem_comp_identifier.identifier" value with "_pdbx_chem_comp_identifier.type='IUPAC CARBOHYDRATE SYMBOL'"
       **       and "_pdbx_chem_comp_identifier.program='PDB-CARE'".
       **
       **  \param: None
       **
       **  \return "_pdbx_chem_comp_identifier.identifier" value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& sugar_iupac_symbol() const;

       /**
       **  Retrieves "_pdbx_chem_comp_identifier.identifier" value with "_pdbx_chem_comp_identifier.type='CONDENSED IUPAC CARBOHYDRATE SYMBOL'"
       **       and "_pdbx_chem_comp_identifier.program='GMML'".
       **
       **  \param: None
       **
       **  \return "_pdbx_chem_comp_identifier.identifier" value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& sugar_condensed_iupac_symbol() const;

       /**
       **  Retrieves "_pdbx_chem_comp_identifier.identifier" value with "_pdbx_chem_comp_identifier.type='SNFG CARBOHYDRATE SYMBOL'"
       **       and "_pdbx_chem_comp_identifier.program='GMML'".
       **
       **  \param: None
       **
       **  \return "_pdbx_chem_comp_identifier.identifier" value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& sugar_snfg_symbol() const;

       /**
       **  Retrieves "_pdbx_chem_comp_feature.value" value with "_pdbx_chem_comp_feature.type='CARBOHYDRATE ISOMER'"
       **       and "_pdbx_chem_comp_feature.source='PDB'".
       **
       **  \param: None
       **
       **  \return "_pdbx_chem_comp_feature.value" value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& sugar_isomer() const;

       /**
       **  Retrieves anomeric carbon atom name
       **
       **  \param: None
       **
       **  \return _anomeric_carbon
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& anomeric_carbon() const;

       /**
       **  Retrieves anomeric oxygen atom name
       **
       **  \param: None
       **
       **  \return _anomeric_oxygen
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& anomeric_oxygen() const;

       /**
       **  Retrieves "_pdbx_chem_comp_feature.value" value with "_pdbx_chem_comp_feature.type='CARBOHYDRATE RING'"
       **       and "_pdbx_chem_comp_feature.source='PDB'".
       **
       **  \param: None
       **
       **  \return "_pdbx_chem_comp_feature.value" value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& sugar_ring() const;

       /**
       **  Retrieves "_pdbx_chem_comp_feature.value" value with "_pdbx_chem_comp_feature.type='CARBOHYDRATE ANOMER'"
       **       and "_pdbx_chem_comp_feature.source='PDB'".
       **
       **  \param: None
       **
       **  \return "_pdbx_chem_comp_feature.value" value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& sugar_anomer() const;

       /**
       **  Retrieves "_pdbx_chem_comp_feature.value" value with "_pdbx_chem_comp_feature.type='CARBOHYDRATE PRIMARY CARBONYL GROUP'"
       **       and "_pdbx_chem_comp_feature.source='PDB'".
       **
       **  \param: None
       **
       **  \return "_pdbx_chem_comp_feature.value" value
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& sugar_primary_carbonyl_group() const;

       /**
       **  Retrieves the sugar atom locant number
       **
       **  \param[in]: atom_name atom name
       **
       **  \return the locant number if it exists.
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const int get_sugar_atom_locant_number(const std::string& atom_name) const;

       /**
       **  Retrieves the sugar side atom tree list
       **
       **  \param[in]: atom_name atom name
       **
       **  \return if it exists.
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::map<std::string, std::string> >& get_sugar_side_atom_tree_list(const std::string& atom_name) const;

       /**
       **  Retrieves meta data from chem_comp category
       **
       **  \param[in]: item - item name from chem_comp category
       **
       **  \return value for selected item
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& getMetaData(const std::string& item) const;

       /**
       **  Merge meta data from input information
       **
       **  \param[in]: meta_data - input meta data 
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void mergeMetaData(const std::map<std::string, std::string>& meta_data);

       /**
       **  Retrieves atom based on atom name
       **
       **  \param[in]: name - atom name
       **
       **  \return found atom
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: NotFoundException - if atom with name does not exist
       */
       const AtomFormat& find_atom(const std::string& name) const;

       /**
       **  Retrieves atom based on index
       **
       **  \param[in]: idx - index to atoms array
       **
       **  \return found atom
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: out_of_range - if idx is out of atoms array range
       */
       const AtomFormat& find_atom(const int& idx) const;

       /**
       **  Retrieves atom record index based on atom name
       **
       **  \param[in]: name - atom name
       **
       **  \return: idx to atom array - if found
       **           -1                - if failed
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception None
       */
       const int find_atom_pos(const std::string& name) const;

       /**
       **  Retrieves atoms
       **
       **  \param: None
       **
       **  \return Constant reference to a vector that contains atoms
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<AtomFormat>& atoms() const;

       /**
       **  Retrieves ordered atoms
       **
       **  \param: None
       **
       **  \return Constant reference to a vector that contains ordered atoms
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<AtomFormat>& order_atoms() const;

       /**
       **  Retrieves bonds
       **
       **  \param: None
       **
       **  \return Constant reference to a vector that contains bonds
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::vector<std::string> >& bonds() const;

       /**
       **  Retrieves atoms stored as std::vector<std::vector<std::string> >
       **
       **  \param: None
       **
       **  \return vector that contains bonds
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::vector<std::string> >& atomArray() const;

       /**
       **  Retrieves atom type/name pair list
       **
       **  \param: None
       **
       **  \return vector that contains atom type/name pair list
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::pair<std::string, std::string> >& getAtomList() const;

       /**
       **  Retrieves bond list represented as atom name pair list
       **
       **  \param: None
       **
       **  \return vector that contains bond list represented 
       **          as atom name pair list
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::pair<std::string, std::string> >& getBondList() const;

       /**
       **  Retrieves linked atoms represented as atom name/vector<atom name> mapping
       **
       **  \param: None
       **
       **  \return a map that contains name/vector<atom name> mapping
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::map<std::string, std::vector<std::string> >& getLinkedAtoms() const;

       /**
       **  Retrieves leaving atoms
       **
       **  \param: None
       **
       **  \return a set that contains leaving atoms' name
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::set<std::string>& getLeavingdAtoms() const;

       /**
       **  Retrieves the leaving atom for a linked atom
       **
       **  \param[in]: linked_atom - the linked atom name
       **
       **  \return a leaving atom name
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string getLeavingdAtom(const std::string& linked_atom) const;

       /**
       **  Retrieves chiral center related information
       **
       **  \param[in]: center_atom_name - chiral center atom name
       **  \param[out]: atomnames - atom names
       **  \param[out]: coords    - coordinates
       **
       **  \return None
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void getChiralCenterAtoms(const std::string& center_atom_name, std::vector<std::string>& atomnames,
                                 std::vector<std::vector<double> >& coords) const;

       /**
       **  Retrieves terminal atoms with leaving attached to them
       **
       **  \param: None
       **
       **  \return terminal atom name set
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::set<std::string>& getTerminalConnectedAtoms() const;

       /**
       **  Retrieves first terminal atom with leaving attached to it & required atom type
       **
       **  \param[in]: type - atom type
       **
       **  \return first terminal atom name
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& getFirstTerminalConnectedAtomByType(const std::string& type) const;

       /**
       **  Retrieves terminal atom(s) with leaving attached to it & required atom type
       **
       **  \param[in]: type - atom type
       **
       **  \return terminal atom name(s)
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::string>& getTerminalConnectedAtomByType(const std::string& type) const;

       /**
       **  Retrieves _chem_comp_atom.atom_id to _chem_comp_atom.pdbx_stnd_atom_id atom name mapping
       **
       **  \param: None
       **
       **  \return atom name mapping
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::map<std::string, std::string>& getAtomIdToStandardAtomIdMapping() const;

       /**
       **  Check if the given atom is a terminal atom
       **
       **  \param[in]: atomName - atom name
       **
       **  \return true/false
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool isTerminalAtom(const std::string& atomName) const;

       /**
       **  Check if it has valid phosphoryl group
       **
       **  \param: None
       **
       **  \return true/false
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool hasValidPhosphorylGroup() const;

       /**
       **  Retrieves bonds between sub component residues
       **
       **  \param[out]: links - atom type
       **          format: vector[0]: residue name 1
       **                  vector[1]: residue numbering 1
       **                  vector[2]: atom name 1
       **                  vector[3]: residue name 2
       **                  vector[4]: residue numbering 2
       **                  vector[5]: atom name 2
       **                  vector[6]: value_order
       **
       **  \return terminal atom name
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void getInterSubComponentBonds(std::vector<std::vector<std::string> >& links) const;

       /**
       **  Get category data table by name
       **
       **  \param[in]: categoryName - category name
       **
       **  \return ISTable pointer if category exists
       **          NULL if category does not exist
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       ISTable* getCategory(const std::string& categoryName);

       /**
       **  Read chemical component information from reference data block
       **
       **  \param[in]: block - reference to chemical component data block
       **
       **  \return true  - if successful
       **          false - if failed
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool Read(Block &block, const bool& with_chem_comp = true, const bool& with_preferred_name = false);

       /**
       **  Write chemical component information to reference data block
       **
       **  \param[in]: compId - default chemical component ID
       **  \param[out]: block - reference to chemical component data block
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void Write(const std::string& compId, Block &block);
};

#endif
