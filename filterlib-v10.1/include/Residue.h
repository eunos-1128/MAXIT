/*
FILE:     Residue.h
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
#ifndef _H_RESIDUE_H_
#define _H_RESIDUE_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Atom.h"
#include "ConnectDic.h"
#include "CrySymmetry.h"
#include "LogUtil.h"
#include "MessageUtil.h"

#define NUM_AA_TERMINAL_ATOMS       5
#define NUM_AA_N_TERMINAL_ATOMS     3
#define NUM_AA_C_TERMINAL_ATOMS     2
#define NUM_AA_LEAVING_ATOMS        3
#define NUM_AA_N_LEAVING_ATOMS      3
#define NUM_NA_TERMINAL_ATOMS       4
#define NUM_NA_5_TERMINAL_ATOMS     3
#define NUM_NA_3_TERMINAL_ATOMS     1
#define NUM_NA_LEAVING_ATOMS        4
#define NUM_NA_5_LEAVING_ATOMS      9
#define NUM_NA_3_LEAVING_ATOMS      3

#define CAAtom_ONLY                 1
#define PAtom_ONLY                  2

#define ZERO_OCCUPANCY           0.001
#define ZERO_OCCUPANCY_LIGAND    0.101

#define LESS_ZERO_B_FACTOR       0.000
#define ZERO_B_FACTOR            0.001
#define LOW_B_FACTOR              10.0
#define HIGH_B_FACTOR            150.0

#define NUM_OVERLAP_CUTOFF          5

class IntegerStringOrder
{
   public:
      IntegerStringOrder();
      const int int_value() const;
      const std::string& str_value() const;
      void set_int_value(const int&);
      void set_str_value(const std::string&);
      IntegerStringOrder& operator=(const IntegerStringOrder&);
      bool operator() (const IntegerStringOrder&, const IntegerStringOrder&) const;
      bool operator<(const IntegerStringOrder&) const;
      bool operator==(const IntegerStringOrder&) const;
   private:
      int  _int_value;
      std::string _str_value;
};

namespace RCSB {
 
class Residue
{
   private:
       ConnectDic* _ccDic;
       LogUtil* _logIo;
       MessageUtil* _messageIo;

       std::string _token;
       std::string _ResName;
       std::string _OrigResName;
       std::string _chnid;
       std::string _res_no;
       std::string _pdb_chnid;
       std::string _pdb_res_no;
       std::string _ins_code;
       std::string _alt_loc;
       std::vector<std::string> _alt_loc_list;
       int _ter_flag;
       int _type;
       std::string _chain_type;
       std::string _entity_id;
       int _chain_break;
       bool _wrong_connectivity;
       bool _disorder_flag;
       bool _full_disorder_flag;
       int _index;
       int _tls_group;
       int _chain_index;
       int _position;
       std::vector<RCSB::Atom*> _atoms;
       std::set<std::string> _uniqueAtomNames;
       std::multimap<std::string, int> _atomNames;
       std::multimap<std::string, int> _atomNameAlts;
       std::multimap<std::string, int> _atomNameAlts_Orininal;
       std::multimap<IntegerStringOrder, int> _atomOrder;
       std::vector<std::string> _missing;
       std::vector<std::string> _extras;
       std::multimap<IntegerStringOrder, int>::iterator _atom_pos;
       void _reset();
       bool _UpdateIndices();
       void _writeAtomRecord(RCSB::Atom* atom, std::string& record);
       void _writeErrorMessage(RCSB::Atom* atom1, RCSB::Atom* atom2, const std::string& type, const std::string& message, const int& Mol_ID); 
       void _writeOccupancyMessage(const std::string &name, const float& occupancy, const int& Mol_ID);
       void _writeOccupancyError(const std::string &name, const float& occupancy, const int& Mol_ID);
       void _find_missing_or_extra_atoms(const int&, const int&, std::set<std::string>&, std::set<std::string>&,
                                         std::set<std::string>&, const ConnectFormat&, const bool&);
       bool _is_OP123_permutation(const std::string&, const std::map<std::string, std::string>&, std::map<std::string, std::string>&);
       bool _connected_to_same_atom(const std::string& first_atom, const std::string& second_atom);
       RCSB::Atom* _find_connected_atom(RCSB::Atom* atom);
       std::vector<RCSB::Atom*> _get_terminal_hydrogens(const std::string& atom_name);
       void _delete_atoms(const std::set<int>& index_set, std::set<long>& atom_set);
       int  _correction_atom_name(const std::string& resname, const int& Mol_ID, const bool& special_na_flag=false, const bool& has_checking_atoms=false);
       void _update_atom_nomenclature(const std::string& old_name, const std::string& new_name, const std::string& atom_type);
       void _update_atom_nomenclature(const std::pair<std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator>& range,
                                      const std::string& new_name, const std::string& atom_type);
       void _createScratchGraph(const int& Mol_ID, std::vector<std::pair<std::string, std::string> >& atomlist, std::vector<std::pair<std::string,
                                std::string> >& bondlist, std::map<std::string, std::pair<std::string, std::string> >& hydrogenMap,
                                const bool& heavy_atom_only_flag = false);
       void _change_hydrogen_names();
       void _getAtomNameTypeMapping(std::map<std::string, std::string>& name_type);
       void _getAtoms(const std::map<std::string, std::string>& name_type, std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string,
                      RCSB::Atom*> > > >& heavyAtoms, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >& hydrogenAtoms);
       void _getAllConformerLists(const std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& heavyAtoms,
                                  std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& allConformerLists);
       void _getPairLists(const unsigned int& index, const std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& a_list,
                          std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& p_list);
       std::set<std::string> _get_alt_loc_set(const std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >& a_list);
       std::set<std::string> _get_set_common_element(const std::set<std::string>& set1, const std::set<std::string>& set2);
       void _getHeavyAtomAndBondList(const int& Mol_ID, const std::vector<std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >&
                                     allConformerLists, std::vector<std::pair<std::string, std::string> >& atomlist,
                                     std::vector<std::pair<std::string, std::string> >& bondlist);
       void _calDistMap(const std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >& atomList, std::vector<std::pair<std::string,
                        std::string> >& atomlist, std::vector<std::pair<std::string, std::string> >& bondlist, std::map<std::string,
                        std::pair<std::pair<std::string, std::string>, double> >& dist_map, std::vector<std::set<std::string> >& atom_sets);
       void _findLinkBetweenSubGraph(const std::map<std::string, std::pair<std::pair<std::string, std::string>, double> >& dist_map, const 
                                     std::vector<std::pair<std::string, std::string> >& atomlist, std::vector<std::set<std::string> >& 
                                     atom_sets, std::vector<std::pair<std::string, std::string> >& bondlist);
       void _getHydrogenBondingInfo(const std::map<std::string, std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > > >& heavyAtoms,
                                    const std::vector<std::pair<std::string, std::pair<std::string, RCSB::Atom*> > >& hydrogenAtoms, std::map<std::string,
                                    std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > > >& hydrogenBondingInfo);
       void _getHydrogenBondList(const std::vector<std::pair<std::string, std::vector<std::pair<std::string, std::string> > > >& hydrogens,
                                 std::vector<std::pair<std::string, std::string> >& linkedHydrogens,
                                 std::map<std::string, std::pair<std::string, std::string> >& hydrogenMap);
       std::string _check_existing_linkage(std::map<std::string, std::vector<std::string> >& bond_info_map);
       std::string _check_missing_linkage(const std::map<std::string, std::vector<std::string> >& bond_info_map,
                                          const std::map<std::string, std::string>& atom_name_type_map);
       void _get_heavy_atom_only_list(std::vector<std::pair<std::string, std::string> >& atomlist, std::vector<std::pair<std::string, std::string> >& bondlist,
                                      const std::set<std::string>& terminal_atom_set);
   public:
       Residue();
       ~Residue() { _reset(); }
       void clear();
       void setCCDic(ConnectDic *ccdic);
       void setLog(LogUtil *logPt);
       void setMessage(MessageUtil *message);
       
       void UpdateIndices();
       const std::string& token() const;
       const std::string& ResName() const;
       const std::string& OrigResName() const;
       const std::string& chnid() const;
       const std::string& res_no() const;
       const std::string& pdb_chnid() const;
       const std::string& pdb_res_no() const;
       const std::string& ins_code() const;
       const std::string& alt_loc() const;
       const std::vector<RCSB::Atom*>& atoms() const;
       const bool& has_alt_loc() const;
       const std::vector<std::string>& alt_loc_list() const;
       const int& ter_flag() const;
       const unsigned int AtomNumbers() const;
       const unsigned int NumAtoms() const;
       const int NonHydrogenAtomNumbers();
       const int& type() const;
       const std::string& chain_type() const;
       const std::string& entity_id() const;
       const int& chain_break() const;
       const bool& wrong_connectivity() const;
       const bool& disorder_flag() const;
       const int& index() const;
       const int& tls_group() const;
       const int& chain_index() const;
       const int& position() const;
       const std::vector<std::string>& missing() const;
       const std::vector<std::string>& extras() const;
       const bool is_single_atom() const;
       const int ca_or_p_atom_only() const;
       void  set_token(const std::string&);
       void  set_ResName(const std::string&);
       void  set_chnid(const std::string&);
       void  set_res_no(const std::string&);
       void  set_pdb_chnid(const std::string&);
       void  set_pdb_res_no(const std::string&);
       void  set_ins_code(const std::string& a);
       void  set_type(const int&);
       void  set_chain_type(const std::string&);
       void  set_entity_id(const std::string&);
       void  set_chain_break(const int&);
       void  set_index(const int&);
       void  set_tls_group(const int&);
       void  set_chain_index(const int&);
       void  set_position(const int&);
       void  set_alt_loc();
       std::string get_side_chain_pattern();
       RCSB::Atom* atoms(const int&);
       RCSB::Atom* GetFirstAtom();
       RCSB::Atom* GetNextAtom();
       void  GetAtomNameList(std::list<std::string>&);
       void  GetAllBonds(std::vector<std::vector<RCSB::Atom*> >& bond_list);
       void  insert_a_atom(RCSB::Atom*, const int&);
       void  delete_a_atom(const std::string&, std::set<long>&);
       void  remove_a_atom(const std::string& name, std::vector<RCSB::Atom*>& atom_list);
       void  find_atom(const std::string&, std::vector<RCSB::Atom*>&);
       void  find_atom_by_type(const std::string&, std::vector<RCSB::Atom*>&);
       RCSB::Atom* find_atom(const std::string&);
       RCSB::Atom* find_atom(const std::string& name, const std::string& alt_loc, const bool& origin_flag=false);
       void  get_atom_type();
       void  reorder_atoms();
       void  find_extra_H_atom();
       void  find_missing_or_extra_atoms(const int&, const int&, std::set<std::string>&, std::set<std::string>&, std::set<std::string>&, const bool&);
       void  find_missing_or_extra_atoms(const ConnectFormat& drug);
       void  refine_missing_atoms(const std::set<std::string>& bonded_atoms);
       void  find_atom_type_mismatch(const int&);
       void  check_hydrogen_bond_distance(const int&);
       void  find_wrong_alt_id(const int&, std::string&);
       void  check_same_coordinates(const int&);
       bool  is_zero_occupancy_residue(std::vector<std::vector<std::string> >&);
       bool  is_zero_occupancy_ligand();
       bool  changed_zero_occupancy();
       float get_partial_occupancy();
       float get_q_score();
       bool  is_pdb_format_compatible();
       void  check_occupancy(const int& Mol_ID, std::list<RCSB::Atom*>& zero_occupancy_atom_list);
       bool  check_occupancy_and_altloc(const int&);
       std::string check_b_factor(int&, int&, int&, int&, std::list<RCSB::Atom*>&);
       void  check_occupancy_and_b_factor(std::vector<std::string>&, std::vector<std::string>&);
       void  update_occupancy(const std::string& new_occ, const bool& check_full_occ_flag = true);
       void  CorresInfoChecking(bool&, std::list<RCSB::Atom*>&, std::list<RCSB::Atom*>&, std::list<RCSB::Atom*>&);
       bool  change_water_nomenclature();
       void  correction_special_na_name(const bool& first_flag, const int& Mol_ID);
       void  AtomNameUpperCase();
       void  correction_atom_name(const int&);
       bool  is_matched_with_CCD(std::vector<std::string>& extra_atom_list, const bool& check_extra_atom_flag = false);
       bool  has_modification(const std::string& chain_type);
       void  check_linkage(const ConnectFormat& drug, std::list<std::string>& error_list, const bool& reset_flag = true);
       void  update_nomenclature(const std::string&, const std::string&);
       void  update_nomenclature(const std::string&, const std::string&, const std::string&);
       void  update_nomenclature(const std::string&, const std::string&, const std::string&, const std::string&, const std::string&, const std::string&);
       void  update_nomenclature(const std::map<std::string, std::string>& mapping);
       void  update_pdb_chnid(const std::string &pdb_chainid);
       void  correction_name(const std::string&);
       void  convert_between_orthogonal_and_fractional(CrySymmetry& crySymm, const int&);
       void  symmetry_operations(CrySymmetry& crySymm, const int&, const int&, const int&, const NDBSYMMETRY&, const bool&);
       void  transformation(double m[4][4]);
       bool  rename(const std::string&, const std::string&);
       void  fix_n_terminal_hydrogen(const bool&, const std::vector<std::map<std::string, std::vector<std::string> > >&);
       void  fix_n_terminal_hydrogen_based_on_mapping();
       void  fix_5_terminal_hydrogen();
       void  fix_5_terminal_hydrogen_based_on_mapping();
       void  check_connectivity(const int&);
       bool  update_atom_name(const std::map<std::string, std::string>& atom_mapping, const int& Mol_ID, const bool& index_flag = true);
       void  update_partial_atom_name(const std::map<std::string, std::string>& atom_mapping, const int& Mol_ID, const std::string& New_ResName);
       bool  update_atom_mapping(const int& Mol_ID, const std::string& New_ResName, const std::map<std::string, std::string>& atom_mapping);
       bool  is_connect(RCSB::Residue*);
       bool  is_overlap(RCSB::Residue*);
};

}

extern const char *_AA_Terminal_Atoms[NUM_AA_TERMINAL_ATOMS];
extern const char *_AA_N_Terminal_Atoms[NUM_AA_N_TERMINAL_ATOMS];
extern const char *_AA_C_Terminal_Atoms[NUM_AA_C_TERMINAL_ATOMS];
extern const char *_AA_Leaving_Atoms[NUM_AA_LEAVING_ATOMS];
extern const char *_AA_N_Leaving_Atoms[NUM_AA_N_LEAVING_ATOMS];
extern const char *_NA_Terminal_Atoms[NUM_NA_TERMINAL_ATOMS];
extern const char *_NA_5_Terminal_Atoms[NUM_NA_5_TERMINAL_ATOMS];
extern const char *_NA_3_Terminal_Atoms[NUM_NA_3_TERMINAL_ATOMS];
extern const char *_NA_Leaving_Atoms[NUM_NA_LEAVING_ATOMS];
extern const char *_NA_5_Leaving_Atoms[NUM_NA_5_LEAVING_ATOMS];
extern const char *_NA_3_Leaving_Atoms[NUM_NA_3_LEAVING_ATOMS];

extern bool is_connect(const std::string &atm_name_a, RCSB::Residue* res_a, const std::string &atm_name_b, RCSB::Residue* res_b);
extern bool is_connect(const std::string &atm_name_a, RCSB::Residue* res_a, const std::string &atm_name_b, RCSB::Residue* res_b, const double cut_off);
extern bool is_connect(RCSB::Residue* res_a, const std::string &atm_type_a, RCSB::Residue* res_b, const std::string &atm_type_b);
extern bool is_glycosidic_link(RCSB::Residue* res_a, RCSB::Residue* res_b);
extern bool is_connect(RCSB::Residue* res_a, const std::string &atm_type_a, RCSB::Residue* res_b, const std::string &atm_type_b,
                       std::vector<std::pair<RCSB::Atom*, RCSB::Atom*> >& pair_list);
extern bool is_polymer_connect(const std::string& chain_type, RCSB::Residue* res1, RCSB::Residue* res2);
extern void get_proline_n_terminal_hydrogen_mapping_list(std::vector<std::map<std::string, std::vector<std::string> > >& atom_name_mapping_list);
extern void get_other_n_terminal_hydrogen_mapping_list(std::vector<std::map<std::string, std::vector<std::string> > >& atom_name_mapping_list);

#endif

// error type:
// alt_code
// alt_conf
// atom_type
// covalent_h_bond
// dissociated_residue
// occupancy
// residue_match
// same_coor
