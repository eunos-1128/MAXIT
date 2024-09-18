/*
FILE:     ConnectFormat.C
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
** \file ConnectFormat.C
**
** \brief Implementation file for ConnectFormat class.
*/

#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "ConnectFormat.h"
#include "ConnectFormat_global.h"
#include "Exceptions.h"
#include "utillib.h"

void ConnectFormat::clear()
{
       _empty_string.clear();
       _empty_vector.clear();
       _empty_tree_list.clear();
       _has_chiral_center = false;
       _has_backbone_terminal_flag = false;
       _nheavyatoms = 0;
       _chem_comp_data.clear();
       _atoms.clear();
       _order_atoms.clear();
       _atom_vectors.clear();
       _atom_pair_list.clear();
       _leaving_atom_set.clear();
       _bonds.clear();
       _bond_pair_list.clear();
       _linked_atom_mapping.clear();
       _linked_atom_index_mapping.clear();
       _atomOrder.clear();
       _terminal_connected_atoms.clear();
       _terminal_connected_atom_type_mapping.clear();
       _category_mapping.clear();
       _atom_id_pdbx_stnd_atom_id_mapping.clear();
       _chem_comp_identifier.clear();
       _chem_comp_feature.clear();
       _anomeric_carbon.clear();
       _anomeric_oxygen.clear();
       _ring_atoms.clear();
       _sugar_locant_mapping.clear();
       _sugar_side_atom_tree_mapping.clear();
}

ConnectFormat& ConnectFormat::operator=(const ConnectFormat& drug)
{
       if (this != &drug) {
            clear();
            _has_chiral_center = drug._has_chiral_center;
            _has_backbone_terminal_flag = drug._has_backbone_terminal_flag;
            _nheavyatoms = drug._nheavyatoms;
            _chem_comp_data = drug._chem_comp_data;
            _atoms = drug._atoms;
            _order_atoms = drug._order_atoms;
            _atom_vectors = drug._atom_vectors;
            _atom_pair_list = drug._atom_pair_list;
            _leaving_atom_set = drug._leaving_atom_set;
            _bonds = drug._bonds;
            _bond_pair_list = drug._bond_pair_list;
            _linked_atom_mapping = drug._linked_atom_mapping;
            _linked_atom_index_mapping = drug._linked_atom_index_mapping;
            _atomOrder = drug._atomOrder;
            _terminal_connected_atoms = drug._terminal_connected_atoms;
            _terminal_connected_atom_type_mapping = drug._terminal_connected_atom_type_mapping;
            _category_mapping = drug._category_mapping;
            _atom_id_pdbx_stnd_atom_id_mapping = drug._atom_id_pdbx_stnd_atom_id_mapping;
            _chem_comp_identifier = drug._chem_comp_identifier;
            _chem_comp_feature = drug._chem_comp_feature;
            _anomeric_carbon = drug._anomeric_carbon;
            _anomeric_oxygen = drug._anomeric_oxygen;
            _ring_atoms = drug._ring_atoms;
            _sugar_locant_mapping = drug._sugar_locant_mapping;
            _sugar_side_atom_tree_mapping = drug._sugar_side_atom_tree_mapping;
       }
       return (*this);
}

const bool& ConnectFormat::has_chiral_center() const
{
       return _has_chiral_center;
}

const bool& ConnectFormat::has_backbone_terminal_flag() const
{
       return _has_backbone_terminal_flag;
}

const int& ConnectFormat::nheavyatoms() const
{
       return _nheavyatoms;
}

const int ConnectFormat::natoms() const
{
       return ((int) _atoms.size());
}

const int ConnectFormat::nbonds() const
{
       return ((int) _bonds.size());
}

const std::string& ConnectFormat::drugname() const
{
       return getMetaData("id");
}

const std::string& ConnectFormat::formula() const
{
       return getMetaData("formula");
}

const std::string& ConnectFormat::chemical_name() const
{
       return getMetaData("name");
}

const std::string& ConnectFormat::synonym() const
{
       return getMetaData("pdbx_synonyms");
}

const std::string& ConnectFormat::sugar_iupac_symbol() const
{
       return _get_chem_comp_identifier("IUPAC CARBOHYDRATE SYMBOL_PDB-CARE");
}

const std::string& ConnectFormat::sugar_condensed_iupac_symbol() const
{
       return _get_chem_comp_identifier("CONDENSED IUPAC CARBOHYDRATE SYMBOL_GMML");
}

const std::string& ConnectFormat::sugar_snfg_symbol() const
{
       return _get_chem_comp_identifier("SNFG CARBOHYDRATE SYMBOL_GMML");
}

const std::string& ConnectFormat::sugar_isomer() const
{
       return _get_chem_comp_feature("CARBOHYDRATE ISOMER_PDB");
}

const std::string& ConnectFormat::anomeric_carbon() const
{
       return _anomeric_carbon;
}

const std::string& ConnectFormat::anomeric_oxygen() const
{
       return _anomeric_oxygen;
}

const std::string& ConnectFormat::sugar_ring() const
{
       return _get_chem_comp_feature("CARBOHYDRATE RING_PDB");
}

const std::string& ConnectFormat::sugar_anomer() const
{
       return _get_chem_comp_feature("CARBOHYDRATE ANOMER_PDB");
}

const std::string& ConnectFormat::sugar_primary_carbonyl_group() const
{
       return _get_chem_comp_feature("CARBOHYDRATE PRIMARY CARBONYL GROUP_PDB");
}

const int ConnectFormat::get_sugar_atom_locant_number(const std::string& atom_name) const
{
       std::map<std::string, int>::const_iterator mpos = _sugar_locant_mapping.find(atom_name);
       if (mpos != _sugar_locant_mapping.end()) return mpos->second;

       return 0;
}

const std::vector<std::map<std::string, std::string> >& ConnectFormat::get_sugar_side_atom_tree_list(const std::string& atom_name) const
{
       std::map<std::string, std::vector<std::map<std::string, std::string> > >::const_iterator mpos = _sugar_side_atom_tree_mapping.find(atom_name);
       if (mpos != _sugar_side_atom_tree_mapping.end()) return mpos->second;

       return _empty_tree_list;
}

const std::string& ConnectFormat::getMetaData(const std::string& item) const
{
       std::map<std::string, std::string>::const_iterator mpos = _chem_comp_data.find(item);
       if (mpos != _chem_comp_data.end()) return mpos->second;
       return _empty_string;
}

void ConnectFormat::mergeMetaData(const std::map<std::string, std::string>& meta_data)
{
       for (std::map<std::string, std::string>::const_iterator mpos = meta_data.begin(); mpos != meta_data.end(); ++mpos) {
            std::map<std::string, std::string>::iterator cpos = _chem_comp_data.find(mpos->first);
            if (cpos != _chem_comp_data.end()) cpos->second = mpos->second;
            else _chem_comp_data.insert(std::make_pair(mpos->first, mpos->second));
       }
}

const AtomFormat& ConnectFormat::find_atom(const std::string& name) const
{
       std::map<std::string, unsigned int>::const_iterator pos = _atomOrder.find(name);
       if (pos != _atomOrder.end())
            return _atoms[pos->second];
       else throw NotFoundException();
}

const AtomFormat& ConnectFormat::find_atom(const int& pos) const
{
       if (pos >= 0 && pos < (int) _atoms.size())
            return _atoms[pos];
       else throw out_of_range("");
}

const int ConnectFormat::find_atom_pos(const std::string& name) const
{
       std::map<std::string, unsigned int>::const_iterator pos = _atomOrder.find(name);
       if (pos != _atomOrder.end())
            return (int) pos->second;
       else return -1;
}

const std::vector<AtomFormat>& ConnectFormat::atoms() const
{
       return _atoms;
}

const std::vector<AtomFormat>& ConnectFormat::order_atoms() const
{
       if (!_order_atoms.empty())
            return _order_atoms;
       else return _atoms;
}

const std::vector<std::vector<std::string> >& ConnectFormat::bonds() const
{
       return _bonds;
}

const std::vector<std::vector<std::string> >& ConnectFormat::atomArray() const
{
       return _atom_vectors;
}

const std::vector<std::pair<std::string, std::string> >& ConnectFormat::getAtomList() const
{
       return _atom_pair_list;
}

const std::vector<std::pair<std::string, std::string> >& ConnectFormat::getBondList() const
{
       return _bond_pair_list;
}

const std::map<std::string, std::vector<std::string> >& ConnectFormat::getLinkedAtoms() const
{
       return _linked_atom_mapping;
}

const std::set<std::string>& ConnectFormat::getLeavingdAtoms() const
{
       return _leaving_atom_set;
}

const std::string ConnectFormat::getLeavingdAtom(const std::string& linked_atom) const
{
       std::map<std::string, unsigned int>::const_iterator apos = _atomOrder.find(linked_atom);
       if (apos != _atomOrder.end()) {
            std::map<std::string, std::vector<std::string> >::const_iterator mpos = _linked_atom_mapping.find(linked_atom);
            if (mpos != _linked_atom_mapping.end()) {
                 std::vector<std::string> hydrogen_atoms;
                 hydrogen_atoms.clear();
                 for (std::vector<std::string>::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                      if (_leaving_atom_set.find(*vpos) != _leaving_atom_set.end()) return *vpos;

                      apos = _atomOrder.find(*vpos);
                      if ((apos != _atomOrder.end()) && _atoms[apos->second].ishydrogen()) hydrogen_atoms.push_back(*vpos);
                 }
                 if (!hydrogen_atoms.empty()) return hydrogen_atoms[0];
            }
       }

       return "";
}

void ConnectFormat::getChiralCenterAtoms(const std::string& center_atom_name, std::vector<std::string>& atomnames,
                                         std::vector<std::vector<double> >& coords) const
{
       atomnames.clear();
       coords.clear();

       std::map<std::string, unsigned int>::const_iterator pos = _atomOrder.find(center_atom_name);
       if (pos == _atomOrder.end()) return;

       std::vector<AtomFormat> Atoms;
       Atoms.clear();
       Atoms.push_back(_atoms[pos->second]);

       std::vector<std::vector<std::string> > bonds = Atoms[0].bonds();
       for (std::vector<std::vector<std::string> >::const_iterator vpos = bonds.begin(); vpos != bonds.end(); ++vpos) {
            pos = _atomOrder.find((*vpos)[0]);
            if (pos == _atomOrder.end()) continue;
            if (_atoms[pos->second].atomtype() == "H" || _atoms[pos->second].atomtype() == "D") continue;
            Atoms.push_back(_atoms[pos->second]);
       }
       if (Atoms.size() < 4) return;

       bool has_ideal_coor = true;
       bool has_model_coor = true;
       for (std::vector<AtomFormat>::const_iterator vpos = Atoms.begin(); vpos != Atoms.end(); ++vpos) {
            if (vpos->get_value(8).empty()) has_model_coor = false;
            if (vpos->get_value(9).empty()) has_model_coor = false;
            if (vpos->get_value(10).empty()) has_model_coor = false;
            if (vpos->get_value(11).empty()) has_ideal_coor = false;
            if (vpos->get_value(12).empty()) has_ideal_coor = false;
            if (vpos->get_value(13).empty()) has_ideal_coor = false;
       }
       if (!has_ideal_coor && !has_model_coor) return;

       std::vector<double> coord;
       for (std::vector<AtomFormat>::const_iterator vpos = Atoms.begin(); vpos != Atoms.end(); ++vpos) {
            atomnames.push_back(vpos->atomname());
            coord.clear();
/*
            if (has_ideal_coor) {
                 coord.push_back(atof(vpos->get_value(11).c_str()));
                 coord.push_back(atof(vpos->get_value(12).c_str()));
                 coord.push_back(atof(vpos->get_value(13).c_str()));
            } else {
                 coord.push_back(atof(vpos->get_value(8).c_str()));
                 coord.push_back(atof(vpos->get_value(9).c_str()));
                 coord.push_back(atof(vpos->get_value(10).c_str()));
            }
*/
            if (has_model_coor) {
                 coord.push_back(atof(vpos->get_value(8).c_str()));
                 coord.push_back(atof(vpos->get_value(9).c_str()));
                 coord.push_back(atof(vpos->get_value(10).c_str()));
            } else {
                 coord.push_back(atof(vpos->get_value(11).c_str()));
                 coord.push_back(atof(vpos->get_value(12).c_str()));
                 coord.push_back(atof(vpos->get_value(13).c_str()));
            }
            coords.push_back(coord);
       }
}

const std::set<std::string>& ConnectFormat::getTerminalConnectedAtoms() const
{
       return _terminal_connected_atoms;
}

const std::string& ConnectFormat::getFirstTerminalConnectedAtomByType(const std::string& type) const
{
       std::map<std::string, std::vector<std::string> >::const_iterator mpos = _terminal_connected_atom_type_mapping.find(type);
       if (mpos != _terminal_connected_atom_type_mapping.end()) return mpos->second[0];
       return _empty_string;
}

const std::vector<std::string>& ConnectFormat::getTerminalConnectedAtomByType(const std::string& type) const
{
       std::map<std::string, std::vector<std::string> >::const_iterator mpos = _terminal_connected_atom_type_mapping.find(type);
       if (mpos != _terminal_connected_atom_type_mapping.end()) return mpos->second;
       return _empty_vector;
}

const std::map<std::string, std::string>& ConnectFormat::getAtomIdToStandardAtomIdMapping() const
{
       return _atom_id_pdbx_stnd_atom_id_mapping;
}

bool ConnectFormat::isTerminalAtom(const std::string& atomName) const
{
       std::map<std::string, unsigned int>::const_iterator pos = _atomOrder.find(atomName);
       if (pos != _atomOrder.end()) {
            if (String::IsEqual(_atoms[pos->second].leaving_atom_flag(), "Y", Char::eCASE_INSENSITIVE)) return true;
       } else return true;
       return false;
}

bool ConnectFormat::hasValidPhosphorylGroup() const
{
       std::set<std::string> phosphoryl_group_atoms, found_oxygen_atoms;
       phosphoryl_group_atoms.clear();
       phosphoryl_group_atoms.insert("P");
       phosphoryl_group_atoms.insert("OP1");
       phosphoryl_group_atoms.insert("OP2");
       phosphoryl_group_atoms.insert("OP3");
       for (std::set<std::string>::const_iterator spos = phosphoryl_group_atoms.begin(); spos != phosphoryl_group_atoms.end(); ++spos) {
            std::map<std::string, unsigned int>::const_iterator mpos = _atomOrder.find(*spos);
            if (mpos == _atomOrder.end()) return false;
            if (((*spos) == "O3P") && !String::IsEqual(_atoms[mpos->second].leaving_atom_flag(), "Y", Char::eCASE_INSENSITIVE)) return false;
       }
       found_oxygen_atoms.clear();
       std::map<std::string, std::vector<std::string> >::const_iterator mpos = _linked_atom_mapping.find("P");
       if (mpos != _linked_atom_mapping.end()) {
            for (std::vector<std::string>::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 if (phosphoryl_group_atoms.find(*vpos) != phosphoryl_group_atoms.end()) found_oxygen_atoms.insert(*vpos);
            }
       }
       if (found_oxygen_atoms.size() == 3) return true;
       return false;
}

void ConnectFormat::getInterSubComponentBonds(std::vector<std::vector<std::string> >& links) const
{
       links.clear();

       std::vector<std::string> data;

       for (std::vector<std::vector<std::string> >::const_iterator bpos = _bonds.begin(); bpos != _bonds.end(); ++bpos) {
            std::map<std::string, unsigned int>::const_iterator mpos1 = _atomOrder.find((*bpos)[0]);
            if (mpos1 == _atomOrder.end()) continue;
            std::map<std::string, unsigned int>::const_iterator mpos2 = _atomOrder.find((*bpos)[1]);
            if (mpos2 == _atomOrder.end()) continue;
            if (_atoms[mpos1->second].get_value(14) == _atoms[mpos2->second].get_value(14) &&
                _atoms[mpos1->second].get_value(15) == _atoms[mpos2->second].get_value(15)) continue;
            // vector[0]: comp_id_1
            // vector[1]: numbering_1
            // vector[2]: atom_id_1
            // vector[3]: component_id_1
            // vector[4]: comp_id_2
            // vector[5]: numbering_2
            // vector[6]: atom_id_2
            // vector[7]: component_id_2
            // vector[8]: bond_type
            // vector[9]: polyper_type_1
            // vector[10]: polyper_type_2
            data.clear();
            data.push_back(_atoms[mpos1->second].get_value(14));
            data.push_back(_atoms[mpos1->second].get_value(15));
            data.push_back(_atoms[mpos1->second].get_value(16));
            data.push_back(_atoms[mpos1->second].get_value(19));
            data.push_back(_atoms[mpos2->second].get_value(14));
            data.push_back(_atoms[mpos2->second].get_value(15));
            data.push_back(_atoms[mpos2->second].get_value(16));
            data.push_back(_atoms[mpos2->second].get_value(19));
            data.push_back((*bpos)[2]);
            data.push_back(_atoms[mpos1->second].get_value(17));
            data.push_back(_atoms[mpos2->second].get_value(17));
            links.push_back(data);
       }
}

ISTable* ConnectFormat::getCategory(const std::string& categoryName)
{
       std::map<std::string, ISTable>::iterator mpos = _category_mapping.find(categoryName);
       if (mpos != _category_mapping.end()) return &(mpos->second);
       return NULL;
}

bool ConnectFormat::Read(Block &block, const bool& with_chem_comp, const bool& with_preferred_name)
{
       ISTable *Table = getTablePtr(block, "chem_comp");
       if (!Table && with_chem_comp) return false;

       std::string cs;
       if (Table) {
            for (int i = 0; i < NUM_CHEM_COMP; ++i) {
                 get_value_clean(cs, Table, 0, _chem_comp[i]);
                 if (cs.empty()) continue;
                 _chem_comp_data.insert(std::make_pair(_chem_comp[i], cs));
            }
            if (_chem_comp_data.find("pdbx_type") == _chem_comp_data.end()) {
                 _chem_comp_data.insert(std::make_pair("pdbx_type", "HETAIN"));
            }
       }

       Table = getTablePtr(block, "chem_comp_atom");
       if (!Table) return false;

       std::vector<std::string> data;
       AtomFormat atom;
       int rowNo = Table->GetNumRows();
       for (int i = 0; i < rowNo; i++) {
            data.clear();
            for (int j = 0; j < NUM_CHEM_COMP_ATOM; ++j) {
                 get_value_clean(cs, Table, i, _chem_comp_atom[j]);
                 if (cs.empty() && ((j == 8) || (j == 9) || (j == 10))) get_value_clean(cs, Table, i, _chem_comp_atom_varient[j]);
                 data.push_back(cs);
            }
            if (data[0].empty()) continue;

            if (!data[2].empty() && data[2] != "H" && data[2] != "D") _nheavyatoms++;
            if (!data[7].empty() && data[7] != "N") _has_chiral_center = true;
            if ((data[20] == "Y") || (data[21] == "Y") || (data[22] == "Y")) _has_backbone_terminal_flag = true;

            atom.clear();
            atom.set_atom(data);
            _atomOrder.insert(std::make_pair(data[0], _atoms.size()));
            _atoms.push_back(atom);
            _atom_vectors.push_back(data);
            _atom_pair_list.push_back(std::make_pair(data[2], data[0]));
            if (data[6] == "Y") _leaving_atom_set.insert(data[0]);
       }

       std::set<std::string> read_tables;
       read_tables.clear();
       read_tables.insert("chem_comp");
       read_tables.insert("chem_comp_atom");

       Table = getTablePtr(block, "chem_comp_bond");
       if (Table) {
            read_tables.insert("chem_comp_bond");

            std::vector<std::string> bond;
            rowNo = Table->GetNumRows();
            for (int i = 0; i < rowNo; i++) {
                 data.clear();
                 for (int j = 0; j < NUM_CHEM_COMP_BOND + 2; ++j) {
                      get_value_clean(cs, Table, i, _chem_comp_bond[j]);
                      if (j == 2) {
                           if (cs.empty()) get_value_clean(cs, Table, i, _chem_comp_bond_varient[j]);
                           String::UpperCase(cs);
                           cs = _get_bond_type(cs);
                      }
                      data.push_back(cs);
                 }
                 _bonds.push_back(data);
                 _bond_pair_list.push_back(std::make_pair(data[0], data[1]));

                 bond.clear();
                 for (unsigned int j = 1; j < data.size(); ++j) {
                      bond.push_back(data[j]);
                 }
                 std::map<std::string, unsigned int>::const_iterator mpos1 = _atomOrder.find(data[0]);
                 if (mpos1 != _atomOrder.end()) _atoms[mpos1->second].add_bond(bond);

                 std::map<std::string, unsigned int>::const_iterator mpos2 = _atomOrder.find(data[1]);
                 if (mpos2 != _atomOrder.end()) {
                      bond[0] = data[0];
                      _atoms[mpos2->second].add_bond(bond);
                 }
                 if ((mpos1 != _atomOrder.end()) && (mpos2 != _atomOrder.end())) {
                      _insert_linked_atom_mapping(data[0], data[1]);
                      _insert_linked_atom_mapping(data[1], data[0]);
                      _insert_linked_atom_index_mapping(mpos1->second, mpos2->second);
                      _insert_linked_atom_index_mapping(mpos2->second, mpos1->second);
                 }
            }
       }

       Table = getTablePtr(block, "pdbx_chem_comp_identifier");
       if (Table) {
            std::string type, program, version, identifier;
            std::vector<std::pair<std::string, std::string> > version_identifier_pair;
            for (unsigned int i = 0; i < Table->GetNumRows(); ++i) {
                 get_value_clean_upper(type, Table, i, "type");
                 get_value_clean_upper(program, Table, i, "program");
                 get_value_clean(version, Table, i, "program_version");
                 get_value_clean(identifier, Table, i, "identifier");
                 if (type.empty() || program.empty() || identifier.empty()) continue;

                 std::string key = type + "_" + program;
                 std::map<std::string, std::vector<std::pair<std::string, std::string> > >::iterator mpos = _chem_comp_identifier.find(key);
                 if (mpos != _chem_comp_identifier.end())
                      mpos->second.push_back(std::make_pair(version, identifier));
                 else {
                      version_identifier_pair.clear();
                      version_identifier_pair.push_back(std::make_pair(version, identifier));
                      _chem_comp_identifier.insert(std::make_pair(key, version_identifier_pair));
                 }
            }
       }

       Table = getTablePtr(block, "pdbx_chem_comp_feature");
       if (Table) {
            std::string type, source, value;
            for (unsigned int i = 0; i < Table->GetNumRows(); ++i) {
                 get_value_clean_upper(type, Table, i, "type");
                 get_value_clean_upper(source, Table, i, "source");
                 get_value_clean(value, Table, i, "value");
                 if (type.empty() || source.empty() || value.empty()) continue;

                 std::string key = type + "_" + source;
                 _chem_comp_feature.insert(std::make_pair(type + "_" + source, value));
            }
       }

       Table = getTablePtr(block, "chem_comp_atom_mapping");
       if (Table) {
            std::string atom_id;
            std::string stnd_atom_id;
            for (unsigned int i = 0; i < Table->GetNumRows(); ++i) {
                 get_value_clean(atom_id, Table, i, "atom_id");
                 get_value_clean(stnd_atom_id, Table, i, "pdbx_stnd_atom_id");
                 if (!atom_id.empty() && !stnd_atom_id.empty()) _atom_id_pdbx_stnd_atom_id_mapping.insert(std::make_pair(atom_id, stnd_atom_id));
            }
       }

       if (with_preferred_name) {
            Table = getTablePtr(block, "pdbx_chem_comp_synonyms");
            if (Table) {
                 std::string type, source, name;
                 for (unsigned int i = 0; i < Table->GetNumRows(); ++i) {
                      get_value_clean_lower(type, Table, i, "type");
                      get_value_clean_upper(source, Table, i, "provenance");
                      if ((type == "preferred") && (source == "PDB")) {
                           get_value_clean(name, Table, i, "name");
                           if (!name.empty()) {
                                std::map<std::string, std::string>::iterator mpos = _chem_comp_data.find("name");
                                if (mpos != _chem_comp_data.end()) mpos->second = name;
                                else _chem_comp_data.insert(std::make_pair("name", name));
                           }
                           break;
                      }
                 }
            }
       }

       std::vector<std::string> tableNames;
       block.GetTableNames(tableNames);
       for (std::vector<std::string>::const_iterator pos = tableNames.begin(); pos != tableNames.end(); ++pos) {
            if (read_tables.find(*pos) != read_tables.end()) continue;

            ISTable *t = getTablePtr(block, *pos);
            if (t) _category_mapping.insert(std::make_pair(*pos, *t));
       }

       _getOrderedAtoms();
       _findTerminalConnectedAtoms();
       _updateSaccharideInformation();

       return true;
}

void ConnectFormat::Write(const std::string& compId, Block &block)
{
       std::string id = compId;

       if (!_chem_comp_data.empty()) {
            std::map<std::string, std::string>::const_iterator mpos = _chem_comp_data.find("id");
            if ((mpos != _chem_comp_data.end()) && !mpos->second.empty()) id = mpos->second;

            ISTable* t = new ISTable("chem_comp");
            for (int i = 0; i < NUM_CHEM_COMP; ++i) t->AddColumn(_chem_comp[i]);
            t->AddRow();
            for (mpos = _chem_comp_data.begin(); mpos != _chem_comp_data.end(); ++mpos) {
                 t->UpdateCell(0, mpos->first, mpos->second);
            }
            block.WriteTable(t);
       }

       if (!_atoms.empty()) {
            ISTable* t = new ISTable("chem_comp_atom");
            t->AddColumn("comp_id");
            for (int i = 0; i < NUM_CHEM_COMP_ATOM; ++i) {
                 t->AddColumn(_chem_comp_atom[i]);
            }

            int row = 0;
            for (std::vector<AtomFormat>::const_iterator apos = _atoms.begin(); apos != _atoms.end(); ++apos) {
                 t->AddRow();
                 t->UpdateCell(row, "comp_id", id);
                 const std::vector<std::string>& data = apos->atom();
                 for (int i = 0; i < NUM_CHEM_COMP_ATOM; ++i) {
                      t->UpdateCell(row, _chem_comp_atom[i], data[i]);
                 }
                 row++;
            }
            block.WriteTable(t);
       }

       if (!_bonds.empty()) {
            ISTable* t = new ISTable("chem_comp_bond");
            t->AddColumn("comp_id");
            for (int i = 0; i < NUM_CHEM_COMP_BOND; ++i) {
                 t->AddColumn(_chem_comp_bond[i]);
            }

            int row = 0;
            for (std::vector<std::vector<std::string> >::const_iterator vpos = _bonds.begin(); vpos != _bonds.end(); ++vpos) {
                 t->AddRow();
                 t->UpdateCell(row, "comp_id", id);
                 for (int i = 0; i < NUM_CHEM_COMP_BOND; ++i) {
                      t->UpdateCell(row, _chem_comp_bond[i], (*vpos)[i]);
                 }
                 row++;
            }
            block.WriteTable(t);
       }
}

void ConnectFormat::_insert_linked_atom_mapping(const std::string& first_atom, const std::string& second_atom)
{
       std::map<std::string, std::vector<std::string> >::iterator mpos = _linked_atom_mapping.find(first_atom);
       if (mpos != _linked_atom_mapping.end()) mpos->second.push_back(second_atom);
       else {
            std::vector<std::string> data;
            data.clear();
            data.push_back(second_atom);
            _linked_atom_mapping.insert(std::make_pair(first_atom, data));
       }
}

void ConnectFormat::_insert_linked_atom_index_mapping(const unsigned int& index1, const unsigned int& index2)
{
       std::map<unsigned int, std::vector<unsigned int> >::iterator mpos = _linked_atom_index_mapping.find(index1);
       if (mpos != _linked_atom_index_mapping.end()) mpos->second.push_back(index2);
       else {
            std::vector<unsigned int> data;
            data.clear();
            data.push_back(index2);
            _linked_atom_index_mapping.insert(std::make_pair(index1, data));
       }
}

void ConnectFormat::_getOrderedAtoms()
{
       if (!_findSubComponentBreak()) return;

       std::set<std::string> component_id_set;
       component_id_set.clear();
       for (unsigned int i = 0; i < _atoms.size(); ++i) {
            component_id_set.insert(_atoms[i].get_value(19));
       }
 
       if (component_id_set.size() > 0) {
            std::vector<unsigned int> t_vec;
            std::map<int, std::vector<unsigned int> > t_map;
            std::map<int, std::map<int, std::vector<unsigned int> > > component_id_residue_number_mapping;
            component_id_residue_number_mapping.clear();
            for (unsigned int i = 0; i < _atoms.size(); ++i) {
                 int component_id = atoi(_atoms[i].get_value(19).c_str());
                 int residue_number = atoi(_atoms[i].get_value(15).c_str());
                 std::map<int, std::map<int, std::vector<unsigned int> > >::iterator cmpos = component_id_residue_number_mapping.find(component_id);
                 if (cmpos != component_id_residue_number_mapping.end()) {
                      std::map<int, std::vector<unsigned int> >::iterator rmpos = cmpos->second.find(residue_number);
                      if (rmpos != cmpos->second.end()) rmpos->second.push_back(i);
                      else {
                           t_vec.clear();
                           t_vec.push_back(i);
                           cmpos->second.insert(std::make_pair(residue_number, t_vec));
                      }
                 } else {
                      t_vec.clear();
                      t_vec.push_back(i);
                      t_map.clear();
                      t_map.insert(std::make_pair(residue_number, t_vec));
                      component_id_residue_number_mapping.insert(std::make_pair(component_id, t_map));
                 }
            }

            _order_atoms.clear();
            for (std::map<int, std::map<int, std::vector<unsigned int> > >::const_iterator cmpos = component_id_residue_number_mapping.begin();
                 cmpos != component_id_residue_number_mapping.end(); ++cmpos) {
                 for (std::map<int, std::vector<unsigned int> >::const_iterator rmpos = cmpos->second.begin(); rmpos != cmpos->second.end(); ++rmpos) {
                      for (std::vector<unsigned int>::const_iterator vpos = rmpos->second.begin(); vpos != rmpos->second.end(); ++vpos) {
                           _order_atoms.push_back(_atoms[*vpos]);
                      }
                 }
            }
       } else {
            std::map<int, int> index;
            index.clear();
            bool has_polymer = false;
            bool has_non_polymer = false;
            for (unsigned int i = 0; i < _atoms.size(); ++i) {
                 int j = atoi(_atoms[i].get_value(15).c_str()) * _atoms.size() + i;
                 index.insert(std::make_pair(j, i));
                 if (_atoms[i].get_value(17) == "polymer") has_polymer = true;
                 if (_atoms[i].get_value(17) == "non-polymer") has_non_polymer = true;
            }

            _order_atoms.clear();
            for (std::map<int, int>::iterator pos = index.begin(); pos != index.end(); ++pos) {
                 _order_atoms.push_back(_atoms[pos->second]);
            }
            if (!has_polymer || !has_non_polymer) return;

            std::vector<AtomFormat> atoms;
            atoms.clear();
            for (std::vector<AtomFormat>::const_iterator pos = _order_atoms.begin(); pos != _order_atoms.end(); ++pos) {
                 if (pos->get_value(17) == "non-polymer") continue;
                 atoms.push_back(*pos);
            }
            for (std::vector<AtomFormat>::const_iterator pos = _order_atoms.begin(); pos != _order_atoms.end(); ++pos) {
                 if (pos->get_value(17) == "polymer") continue;
                 atoms.push_back(*pos);
            }
            _order_atoms = atoms;
       }
}

bool ConnectFormat::_findSubComponentBreak()
{
       std::string numbering = "";
       int count = 0;
       bool is_hydrogen_break = false;
       std::set<unsigned int> residue_breaks;
       for (unsigned int i = 0; i < _atoms.size(); ++i) {
            if (_atoms[i].get_value(15) != numbering) {
                 numbering = _atoms[i].get_value(15);
                 count++;
            }
            // insert sub component break point if two sub components
            // have different residue names
            if (i && _atoms[i].get_value(14) != _atoms[i-1].get_value(14))
                 residue_breaks.insert(i);
            // insert sub component break point between heavy atom & hydrogen
            if (i && (_atoms[i].get_value(2) == "H" || _atoms[i].get_value(2) == "D") &&
                _atoms[i-1].get_value(2) != "H" && _atoms[i-1].get_value(2) != "D") {
                 residue_breaks.insert(i);
                 is_hydrogen_break = true;
            }
       }
       if (count > 0) return true;

       if (is_hydrogen_break && residue_breaks.size() == 1) residue_breaks.clear();

       // convert residue point break set into sub component range pair
       std::vector<std::pair<unsigned int, unsigned int> > pair_array;
       pair_array.clear();
       if (residue_breaks.empty()) // whole ligand has single range pair
            pair_array.push_back(std::make_pair(0, _atoms.size() - 1));
       else {
            std::vector<unsigned int> break_array;
            break_array.clear();
            for (std::set<unsigned int>::iterator
                 pos = residue_breaks.begin(); pos != residue_breaks.end(); ++pos) {
                 break_array.push_back(*pos);
            }

            // insert first sub component range pair
            pair_array.push_back(std::make_pair(0, break_array[0] - 1));
            // insert middle sub component(s) range pair
            for (unsigned int i = 0; i < break_array.size() - 1; ++i) {
                 pair_array.push_back(std::make_pair(break_array[i],
                                                     break_array[i+1]-1));
            }
            // insert last sub component range pair
            unsigned int i = break_array.size() - 1;
            pair_array.push_back(std::make_pair(break_array[i], _atoms.size() - 1));
       }

       for (std::vector<std::pair<unsigned int, unsigned int> >::iterator
            pos = pair_array.begin(); pos != pair_array.end(); ++pos) {
            _getSameSubComponentBreak(pos->first, pos->second, residue_breaks);
       }

       std::string type = "non-polymer";
       if (!residue_breaks.empty() && getMetaData("pdbx_release_status") == "REF_ONLY") type = "polymer";
       count = 1;
       bool first_hydrogen = true;
       for (unsigned int i = 0; i < _atoms.size(); ++i) {
            if (residue_breaks.find(i) != residue_breaks.end()) {
                 count++;
                 if (first_hydrogen && (_atoms[i].get_value(2) == "H" || _atoms[i].get_value(2) == "D")) {
                      count = 1;
                      first_hydrogen = false;
                 }
            }
            if (_atoms[i].get_value(15).empty())
                 _atoms[i].set_value(15, String::IntToString(count));
            if (_atoms[i].get_value(17).empty()) _atoms[i].set_value(17, type);
       }

       if (!residue_breaks.empty())
            return true;
       else return false;
}

void ConnectFormat::_getSameSubComponentBreak(const unsigned int& begin, const unsigned int& end, std::set<unsigned int>& residue_breaks)
{
       for (unsigned int i = begin + 1; i <= end; ++i) {
            if (_atoms[i].get_value(16) == _atoms[begin].get_value(16)) residue_breaks.insert(i);
       }
}

void ConnectFormat::_findTerminalConnectedAtoms()
{
       std::string id = "";
       std::map<std::string, std::string>::const_iterator mmpos = _chem_comp_data.find("id");
       if (mmpos != _chem_comp_data.end()) id = mmpos->second;

       std::vector<std::string> data;
       if (id == "ACE") {
            _terminal_connected_atoms.insert("C");
            data.clear();
            data.push_back("C");
            _terminal_connected_atom_type_mapping.insert(std::make_pair("C", data));
            return;
       } else if (id == "NH2") {
            _terminal_connected_atoms.insert("N");
            data.clear();
            data.push_back("N");
            _terminal_connected_atom_type_mapping.insert(std::make_pair("N", data));
            return;
       }

       if (_atoms.empty()) return;

       if (_atoms.size() == 1) {
            _terminal_connected_atoms.insert(_atoms[0].atomname());
            data.clear();
            data.push_back(_atoms[0].atomname());
            _terminal_connected_atom_type_mapping.insert(std::make_pair(_atoms[0].atomtype(), data));
            return;
       }

       for (std::vector<AtomFormat>::const_iterator pos = _atoms.begin(); pos != _atoms.end(); ++pos) {
            if (!String::IsEqual(pos->leaving_atom_flag(), "Y", Char::eCASE_INSENSITIVE)) continue;

            const std::vector<std::vector<std::string> >& bonds = pos->bonds();
            for (std::vector<std::vector<std::string> >::const_iterator bpos = bonds.begin(); bpos != bonds.end(); ++bpos) {
                 std::map<std::string, unsigned int>::const_iterator mpos = _atomOrder.find((*bpos)[0]);
                 if (mpos == _atomOrder.end()) continue;
                 if (_atoms[mpos->second].atomtype() == "H" || _atoms[mpos->second].atomtype() == "D") continue;

                 _terminal_connected_atoms.insert(_atoms[mpos->second].atomname());
                 if (String::IsEqual(_atoms[mpos->second].leaving_atom_flag(), "Y", Char::eCASE_INSENSITIVE)) continue;

                 std::map<std::string, std::vector<std::string> >::iterator tpos = _terminal_connected_atom_type_mapping.find(_atoms[mpos->second].atomtype());
                 if (tpos != _terminal_connected_atom_type_mapping.end())
                      tpos->second.push_back(_atoms[mpos->second].atomname());
                 else {
                      data.clear();
                      data.push_back(_atoms[mpos->second].atomname());
                      _terminal_connected_atom_type_mapping.insert(std::make_pair(_atoms[mpos->second].atomtype(), data));
                 }
            }
       }
}

std::string ConnectFormat::_get_bond_type(const std::string& inBondType)
{
       std::string outBondType = inBondType;

       if (inBondType.find("SING") == 0)      outBondType = "SING";
       else if (inBondType.find("DOUB") == 0) outBondType = "DOUB";
       else if (inBondType.find("TRIP") == 0) outBondType = "TRIP";
       else if (inBondType.find("QUAD") == 0) outBondType = "QUAD";
       else if (inBondType.find("AROM") == 0) outBondType = "AROM";
       else if (inBondType.find("POLY") == 0) outBondType = "POLY";
       else if (inBondType.find("DELO") == 0) outBondType = "DELO";
       else if (inBondType.find("PI") == 0)   outBondType = "PI";

       return outBondType;
}

const std::string& ConnectFormat::_get_chem_comp_identifier(const std::string& key) const
{
       std::map<std::string, std::vector<std::pair<std::string, std::string> > >::const_iterator mpos = _chem_comp_identifier.find(key);
       if (mpos != _chem_comp_identifier.end()) return mpos->second[0].second;
       return _empty_string;
}

const std::string& ConnectFormat::_get_chem_comp_feature(const std::string& key) const
{
       std::map<std::string, std::string>::const_iterator mpos = _chem_comp_feature.find(key);
       if (mpos != _chem_comp_feature.end()) return mpos->second;
       return _empty_string;
}

void ConnectFormat::_updateSaccharideInformation()
{
       std::string type = getMetaData("type");
       String::UpperCase(type);
       if (type.find("SACCHARIDE") == std::string::npos) return;

       std::string statusCode = getMetaData("pdbx_release_status");
       String::UpperCase(statusCode);
       if (statusCode == "REF_ONLY") return;

       std::vector<std::vector<unsigned int> > ringLists;
       find_rings(_atoms.size(), _linked_atom_index_mapping, ringLists);

       if (ringLists.empty()) return;

       if (sugar_primary_carbonyl_group() == "aldose") {
            std::map<std::string, unsigned int>::const_iterator mpos = _atomOrder.find("C1");
            if (mpos != _atomOrder.end()) {
                 _anomeric_carbon = "C1";
                 _anomeric_oxygen = "O1";
            }
       } else if (sugar_primary_carbonyl_group() == "ketose") {
            std::map<std::string, unsigned int>::const_iterator mpos = _atomOrder.find("C2");
            if (mpos != _atomOrder.end()) {
                 _anomeric_carbon = "C2";
                 _anomeric_oxygen = "O2";
            }
       }

       for (std::vector<std::vector<unsigned int> >::const_iterator rpos = ringLists.begin(); rpos != ringLists.end(); ++rpos) {
            bool found = false;
            for (std::vector<unsigned int>::const_iterator vpos = rpos->begin(); vpos != rpos->end(); ++vpos) {
                 _ring_atoms.insert(_atoms[*vpos].atomname());
                 if (_atoms[*vpos].atomname() == _anomeric_carbon) found = true;
            }
            if (found || (ringLists.size() == 1)) break;
            _ring_atoms.clear();
       }

       if (_ring_atoms.empty()) return;

       if (!_anomeric_carbon.empty()) {
            std::map<std::string, std::vector<std::string> >::const_iterator mpos = _linked_atom_mapping.find(_anomeric_carbon);
            if (mpos != _linked_atom_mapping.end()) {
                 for (std::vector<std::string>::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                      if (((*vpos)[0] == 'H') || (_ring_atoms.find(*vpos) != _ring_atoms.end())) continue;
                      if ((vpos->size() == 2) && ((*vpos)[1] == _anomeric_carbon[1])) {
                           _anomeric_oxygen = *vpos;
                           break;
                      }
                 }
            }
       }

       // map-key: atom name
       // map-val: atom type
       std::vector<std::map<std::string, std::string> > side_atom_tree_list;
       std::map<std::string, int> side_atom_map;
       std::set<std::string> exclude_atom_set;
       std::set<std::string> exclude_set = _ring_atoms;
       
       for (std::set<std::string>::const_iterator spos = _ring_atoms.begin(); spos != _ring_atoms.end(); ++spos) {
            int locant_number = _get_number_part_from_atom_name(*spos);
            _sugar_locant_mapping.insert(std::make_pair(*spos, locant_number));

            std::map<std::string, std::vector<std::string> >::const_iterator mpos = _linked_atom_mapping.find(*spos);
            if (mpos == _linked_atom_mapping.end()) continue;

            side_atom_map.clear();
            for (std::vector<std::string>::const_iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 if (exclude_set.find(*vpos) != exclude_set.end()) continue;
                 side_atom_map.insert(std::make_pair(*vpos, locant_number));
            }
            if (side_atom_map.empty()) continue;

            exclude_atom_set.clear();
            exclude_atom_set.insert(*spos);
            side_atom_tree_list.clear();
            _update_side_atom_info(side_atom_map, exclude_atom_set, side_atom_tree_list); 

            if (!side_atom_tree_list.empty()) _sugar_side_atom_tree_mapping.insert(std::make_pair(*spos, side_atom_tree_list));
       }
}

int ConnectFormat::_get_number_part_from_atom_name(const std::string& atom_name)
{
       if (isdigit(atom_name[atom_name.size() - 1]) != 0) {
            for (unsigned int i = atom_name.size() - 1; i >= 0; --i) {
                 if (isdigit(atom_name[i]) == 0) {
                      return atoi(atom_name.substr(i+1).c_str());
                 }
                 if (i == 0) break;
            }
       }
       return 0;
}

void ConnectFormat::_update_side_atom_info(const std::map<std::string, int>& include_map, std::set<std::string>& exclude_set,
                                           std::vector<std::map<std::string, std::string> >& side_atom_tree_list)
{
       std::map<std::string, std::string> side_atom_map;
       side_atom_map.clear();

       std::map<std::string, int> child_atom_map;
       child_atom_map.clear();

       for (std::map<std::string, int>::const_iterator impos = include_map.begin(); impos != include_map.end(); ++impos) {
            exclude_set.insert(impos->first);
            int locant_number = _get_number_part_from_atom_name(impos->first);
            if (locant_number == 0) locant_number = impos->second;
            _sugar_locant_mapping.insert(std::make_pair(impos->first, locant_number));

            std::map<std::string, unsigned int>::const_iterator mpos = _atomOrder.find(impos->first);
            if (mpos != _atomOrder.end()) side_atom_map.insert(std::make_pair(_atoms[mpos->second].atomname(), _atoms[mpos->second].atomtype()));

            std::map<std::string, std::vector<std::string> >::const_iterator lmpos = _linked_atom_mapping.find(impos->first);
            if (lmpos == _linked_atom_mapping.end()) continue;

            for (std::vector<std::string>::const_iterator vpos = lmpos->second.begin(); vpos != lmpos->second.end(); ++vpos) {
                 if (exclude_set.find(*vpos) != exclude_set.end()) continue;
                 child_atom_map.insert(std::make_pair(*vpos, locant_number));
            }
       }
       if (!side_atom_map.empty()) side_atom_tree_list.push_back(side_atom_map);

       if (child_atom_map.empty()) return;

       _update_side_atom_info(child_atom_map, exclude_set, side_atom_tree_list);
}
