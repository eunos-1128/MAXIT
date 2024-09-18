/*
FILE:     Molecule.h
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
#ifndef _HMolecule_H_
#define _HMolecule_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "Chain.h"
#include "ConnectDic.h"
#include "LogUtil.h"
#include "MessageUtil.h"

typedef struct {
        std::string PDB_ChainId;
        std::vector<std::vector<std::string> > seqs;
        std::vector<std::map<std::string, std::string> > linkages;
        std::vector<std::map<std::string, std::string> > descriptors;
} BRANCH_INFO;

namespace RCSB {

class Molecule
{
   private:
       ConnectDic* _ccDic;
       LogUtil* _logIo;
       MessageUtil* _messageIo;
       bool _format_checking;
       bool _residue_reorder;
       bool _is_pdb_input_format;
       bool _found_sugar_flag;
       bool _rename_residue_flag;
       bool _find_sugar_chain_flag;
       int  _Mol_ID;
       int  _index;
       int  _chain_index;
       int  _type;
       std::vector<RCSB::Chain*> _chains;
       std::vector<RCSB::Residue*> _residues;
       std::list<RCSB::Atom*> _atoms;

       // key: residue index
       // value: _residues position index
       std::map<int, int> _resIndex;

       // key: chain order
       // value: _chains position index
       std::multimap<int, int> _chainOrder;

       // key: chain index
       // value: _chains position index
       std::multimap<int, int> _chainIndex;
       std::multimap<int, int>::iterator _order_pos;

       // key: pdbchnid_resnme_pdbresnum_inscode 
       // value: _residues position index
       std::multimap<std::string, int> _pdbIndex;

       // key: pdbchnid_resnme_pdbresnum_inscode 
       // value: _residues position index
       std::multimap<std::string, int> _origPdbIndex;

       // key: pdbchnid_pdbresnum_inscode 
       // value: _residues position index
       std::multimap<std::string, int> _origPdbNumIndex;

       // key: asymid_resnme_labelresnum_
       // value: _residues position index
       // std::multimap<std::string, int> _cifIndex;

       // key: asymid_resnme_labelresnum_inscode
       // value: _residues position index
       std::multimap<std::string, int> _origCifIndex;

       // key: PDB_ChainID
       // value: _chains position index
       std::multimap<std::string, int> _chainIdIndex;

       // key: ChainID (Asym ID)
       // value: _chains position index
       std::multimap<std::string, int> _asymIdIndex;

       // key: asymId
       // value: original asymId
       std::map<std::string, std::string> _sugarMergedAsmyIDMap;

       // key: new PDB chain ID
       // value: old PDB chain IDs
       std::map<std::string, std::set<std::string> > _PDBchainID_changed_map;

       std::map<std::string, BRANCH_INFO> _Branch_Seq_Scheme_Mapping;

       std::set<int> _removed_residue_index_set, _removed_chain_index_set;

       // key: pdbchnid_resnme_pdbresnum_inscode
       // value: new ResName
       std::map<std::string, std::string> _residue_re_name_mapping;

       void Reset();
       void _find_sugar_chains(std::map<unsigned int, std::list<std::list<unsigned int> > >& sugar_chains, const std::vector<unsigned int>&
                               sugar_residue_set, const std::set<std::string>& link_residue_set, const bool& find_sugar_link_flag);
       // void find_sugar_links();
       void _remove_empty_chains(const std::set<int>& chain_set);
       void _update_special_na_residues(const std::vector<std::vector<unsigned int> >& idx_lists);
       void _insert_a_atom(RCSB::Atom* atom, std::map<std::string, RCSB::Residue*>& residue_mapping, std::map<std::string, int>& last_number_mapping,
                           std::list<RCSB::Residue*>& residue_list, const bool& split_flag, const std::string& polymer_type);
       bool _is_split_residue(RCSB::Residue*, std::map<std::string, int>&);
       void _get_PRD_involved_chains(const std::vector<std::vector<std::string> >& mapping, std::set<int>& chain_set, bool& _successful);
       void _get_PRD_involved_residue_and_atoms(const std::set<int>& chain_set, std::multimap<std::string, int>& residue_mapping,
                                                std::multimap<std::string, RCSB::Atom*>& found_atoms, bool& _successful);
       void _get_involved_atoms(RCSB::Residue* res, std::multimap<std::string, RCSB::Atom*>& found_atoms);
       bool _get_new_residue_list(const std::vector<std::vector<std::string> >& mapping, const std::multimap<std::string, RCSB::Atom*>& found_atoms, const
                                  std::vector<std::vector<std::string> >& in_links, std::list<RCSB::Residue*>& residue_list, std::map<std::string, std::string>&
                                  polymer_type_mapping, std::vector<std::pair<std::string, std::vector<std::vector<RCSB::Atom*> > > >& out_links);
       void _remove_old_residues_and_chains(const std::multimap<std::string, int>& residue_mapping, const std::set<int>& chain_set);
       void _add_new_residues_and_chains(std::list<RCSB::Residue*>& residue_list, const std::map<std::string, SEQ>& seqs, const std::string& assigned_chain_type);
       bool _update_chain_nomenclature(const std::map<int, std::list<std::pair<int, std::string> > >& mapping);
       void _analysis_chain_group(const std::string PDB_Chain_ID, const std::vector<std::pair<unsigned int, unsigned int> >& pair_array,
                                  const bool& has_ter_card, const int& polymer_status, const std::string& chain_type_from_seq,
                                  const std::vector<std::vector<std::string> >& mapping, const std::set<std::string>& branch_chain_id_set,
                                  std::set<unsigned int>& chain_breaks, std::set<std::string>& polymer_index,
                                  std::map<std::string, std::string>& chain_type_mapping);
       std::string _analysis_chain_type(const unsigned int& start, const unsigned int& end, bool& has_ATOM_token, bool& is_water_chain, bool& is_100_percent);
       bool _check_linkage_between_residues(const std::string& chain_type, const unsigned int& start, const unsigned int& end);
       void _insert_overlap_residue_index_map(const unsigned int& r_index_1, const unsigned int& r_index_2,
                                     std::map<unsigned int, std::set<unsigned int> >& overlap_index_map);
       void _insert_chains(std::map<std::string, std::vector<RCSB::Chain*> >& mapping_chains, RCSB::Chain* chain);
       void _get_last_number(int& last_number, const std::vector<RCSB::Chain*>& chains);
   public:
       Molecule();
       ~Molecule() { Reset(); }
       void clear();
       const int& Mol_ID() const;
       const int& index() const;
       const int Num_Chain() const;
       const int& type() const;
       const bool has_sugar_entity();
       const bool need_sugar_entity_update();
       const std::vector<RCSB::Residue*>& Residues() const;
       const std::list<RCSB::Atom*>& atoms() const;
       void setCCDic(ConnectDic *ccdic);
       void setLog(LogUtil *logPt);
       void setMessage(MessageUtil *message);
       void setFormatChecking();
       void unsetResidueReroder();
       void setPDBInputFormat();
       void set_Carbohydrate_Annotation_Info(const std::map<std::string, BRANCH_INFO>& mapping);
       void setRenameResidueFlag(const bool& flag);
       void setFindSugarChainFlag(const bool& flag);
       void set_Mol_ID(const int& mol_id);
       void set_index(const int& index_id);
       void set_chain_terminal_card();
       const int MaxNumEntity();
       void update_residue_indices();
       void insert_a_atom(RCSB::Atom* atom);
       bool insert_a_atom(const std::vector<std::string>& vals, const bool&, const int&);
       void insert_a_residue(RCSB::Residue* residue);
       void changed_zero_occupancy();
       bool find_residues(const bool& UNX_flag);
       bool find_residues(const bool& split_flag, const std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > >& atom_mapping,
                          std::set<std::string>& instanceid_set, const bool& UNX_flag);
       void find_chains(const std::set<std::string>& link_residue_set, const bool& find_sugar_link_flag);
       void find_chains(const std::map<std::string, std::vector<std::vector<std::string> > >& scheme_mapping, std::map<std::string, SEQ>& seqs,
                        const std::set<std::string>& input_link_residue_set, const bool& keep_solvent_position_flag, const bool& find_sugar_link_flag);
       void InternalOrder(const bool& numbering_based_order_flag = false);
       void GenerateInternalOrderIndex();  
       void AssignAsymId();
       void CheckUniquePDBNumbering();
       int  ClearNonAsymAndWaterChains();
       RCSB::Chain* GetFirstChain();
       RCSB::Chain* GetNextChain();
       RCSB::Chain* GetIndexChain(const int&, bool&);
       RCSB::Chain* GetPolyChain(const std::string&);
       RCSB::Chain* GetAsymChain(const std::string&);
       RCSB::Chain* GetOrigAsymSugarChain(const std::string&);
       void GetPdbChains(const std::string& pdb_chnid, std::vector<RCSB::Chain*>& chain_list);
       RCSB::Residue* find_pdb_residue(const std::string&, const std::string&, const int&, const std::string&);
       RCSB::Residue* find_pdb_residue(const std::string&, const std::string&, const std::string&, const std::string&);
       RCSB::Residue* find_pdb_residue(const std::string&);
       RCSB::Residue* find_orig_pdb_residue(const std::string& chnid, const std::string& resname, const std::string& resnum, const std::string& ins_code);
       RCSB::Residue* find_orig_pdb_residue(const std::string& chnid, const std::string& resnum, const std::string& ins_code);
       // RCSB::Residue* find_cif_residue(const std::string&, const std::string&, const int);
       // RCSB::Residue* find_cif_residue(const std::string&, const std::string&, const std::string&);
       RCSB::Residue* find_orig_cif_residue(const std::string& chnid, const std::string& resname, const std::string& resnum);
       RCSB::Residue* find_residue(const int&, bool&);
       RCSB::Atom* find_atom(RCSB::Atom*);
       // void GenResidueIndex();
       void insert_a_chain(RCSB::Chain*, const bool& set_index_flag = true);
       void insert_a_chain(const std::list<RCSB::Residue*>& residue_list);
       RCSB::Chain* insert_a_chain(RCSB::Residue* res, const bool& set_prev_entity_id_flag = false);
       void build_residue_index();
       int  count_atom_no(const std::string&);
       int  count_atom_no(const std::set<std::string>& token_set);
       int  count_total_atom_number();
       void update_cif_nomenclature();
       void GetChainsAndWaters(std::vector<RCSB::Chain*>& polymer_chains, std::vector<RCSB::Chain*>& carbohydrate_chains,
                               std::vector<RCSB::Chain*>& nonpolymer_chains, std::list<RCSB::Residue*>& water_lists);
       void GetChains(const std::string& pdb_chainid, std::vector<RCSB::Chain*>& polymer_chains, std::vector<RCSB::Chain*>& nonpolymer_chains,
                      std::vector<RCSB::Chain*>& waters);
       void GetChains(std::set<std::string>& pdb_chainids, std::map<std::string, std::vector<RCSB::Chain*> >& polymer_chains,
                      std::map<std::string, std::vector<RCSB::Chain*> >& other_chains);
       void GetChains(std::vector<std::string>& pdb_chain_ids, std::map<std::string, std::vector<RCSB::Chain*> >& polymer_chains,
                      std::map<std::string, std::vector<RCSB::Chain*> >& non_polymer_chains, std::map<std::string, std::vector<RCSB::Chain*> >& water_chains);
       void GetMolInfo(std::vector<std::string>& pdb_chain_ids, std::map<std::string, RCSB::Chain*>& polymers, std::map<std::string, std::list<RCSB::Residue*> >&
                       nonpolymers, std::map<std::string, std::list<RCSB::Residue*> >& waters);
       bool PolymerUpdate(const std::vector<std::vector<std::string> >& mapping, const std::vector<std::vector<std::string> >& in_links,
                          std::vector<std::pair<std::string, std::vector<std::vector<RCSB::Atom*> > > >& out_links);
       bool PrdUpdate(const std::vector<std::vector<std::string> >& mapping, const std::map<std::string, SEQ>& seqs, const
                      std::vector<std::vector<std::string> >& in_links, std::vector<std::pair<std::string, std::vector<std::vector<RCSB::Atom*> > > >& out_links);
       bool PrdUpdateAtom(const std::map<std::string, std::map<std::string, std::string> >& atom_mapping);
       bool PrdUpdateNumbering(const std::map<std::string, std::string>& num_mapping);
       bool MergePolymer(const int& new_entity_id, const std::string& new_asym_id, const std::string& new_chain_id,
                         const std::vector<std::vector<std::string> >& group, std::vector<std::string>& seqs,
                         std::vector<std::vector<int> >& merged_entityids);
       bool SplitPolymer(const std::vector<std::vector<std::string> >& splits, const std::vector<std::vector<std::string> >& chain_info,
                         std::map<int, std::pair<std::vector<std::string>, std::vector<std::string> > >& splited_entityids);
       bool EditPolymer(const std::vector<std::vector<std::string> >& deletes);
       void MergeResidue(const std::vector<std::pair<std::string, std::vector<std::vector<std::string> > > >& merge_residue_list);
       bool UpdateInstance(const std::map<std::string, std::pair<std::string, std::map<std::string, std::string> > >& atom_mapping);
       bool update_residue_nomenclature(const std::map<std::string, std::map<std::string, std::string> >& mapping);
       void update_numbering(const std::map<std::string, std::map<std::string, std::string> >& mapping);
       void update_numbering(const std::map<std::string, std::vector<std::string> >& mapping);
       void update_chainid(const std::map<std::string, std::string>& mapping);
       void assign_carbohydrate_chain_id(std::map<std::string, std::string>& sugar_chain_id_mapping);
       void merged_linked_sugar_chains(const std::set<std::string>& residue_set, const std::set<std::string>& link_set);
};

}
/*
extern void get_residue_indices(const std::vector<Molecule*>&, const int&, int&, std::vector<int>&,
                                const std::vector<std::vector<std::string> >&, std::vector<std::string>&);
extern void get_atom_indices(const std::vector<Molecule*>&, const int&, int&, std::vector<RCSB::Atom*>&,
                             const std::vector<std::vector<std::string> >&, std::vector<std::string>&);
*/
#endif
