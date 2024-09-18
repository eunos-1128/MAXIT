/*
FILE:     Chain.h
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
#ifndef _H_CHAIN_H_
#define _H_CHAIN_H_

#include <list>
#include <map>
#include <set>
#include <string>

#include "ConnectDic.h"
#include "Residue.h"

#define NUM_FIELD                 7

#define ATOMN_TYPE_DNA_ONLY       1
#define ATOMN_TYPE_RNA_ONLY       2
#define ATOMN_TYPE_DNA_RNA        3
#define ATOMN_TYPE_TRNA           6

typedef struct {
        std::string ChainId;
        std::string PDB_ChainId;
        std::string PolyType;
        std::string EvidenceCode;
        std::string chain_type;
        std::string entity_id;
        bool never_used;
        int order;
        std::vector<std::string> res;
} SEQ;

typedef struct {
        std::string Field[NUM_FIELD];
        std::string InsCode;
        int  ResIndex;
} _FIELD;

typedef struct {
        std::vector<std::string> _seqa;
        std::vector<std::string> _seqb;
} _DataContainer;


namespace RCSB {

class Chain
{
   private:
        ConnectDic* _ccDic;
        LogUtil* _logIo;
        MessageUtil* _messageIo;

        std::string _ChainID;
        std::string _PDB_ChainID;
        std::string _PUB_ChainID;
        std::string _PDB_ChainID_Flag;
        std::string _details;
        std::string _PolyType;
        std::string _EvidenceCode;
        std::string _chain_type;
        std::string _prev_entity_id;
        std::string _descriptor;
        std::string _entity_id;
        std::string _entity_key;
        bool _has_sequence;
        bool _empty_chain;
        bool _missing_sequence;
        int  _ca_or_p_atom_only;
        bool _isModified;
        int _index;
        int _order;
        int _na_type;
        int _ResidueNumbers;
        std::vector<std::vector<RCSB::Residue*> > _residues;
        std::vector<_FIELD*> _SeqRes;
        std::map<std::string, int> _numIndex;
        std::map<std::string, std::pair<int, int> > _pdbIndex;
        // std::map<std::string, std::pair<int, int> > _cifIndex;
        std::vector<std::map<std::string, std::string> > _linear_descriptors;
        std::vector<std::map<std::string, std::string> > _branch_links;
        std::list<int> _firstresidueIndex;
        std::list<int>::iterator _first_pos;
        void clear();
        void Reset();
        void _get_entity_key();
        void __set_chain_index();
        void check_alignment(const _DataContainer&, const std::vector<int>&, std::vector<std::vector<int> >&);
        bool find_match_list(const _DataContainer&, const std::vector<int>&, const std::vector<std::vector<int> >&, const std::vector<int>&,
                             const std::vector<int>&, const std::vector<int>&, const std::vector<int>&, std::vector<int>&);
        void check_unique_number();
        bool _IsConnect(RCSB::Residue *res1, RCSB::Residue *res2);
        void _update_seqres();
        bool _check_need_renumber();
   public:
        Chain()  { clear(); }
        ~Chain() { Reset(); }
        const bool empty() const;
        const std::string& ChainID() const;
        const std::string& PDB_ChainID() const;
        const std::string& PUB_ChainID() const;
        const std::string& PDB_ChainID_Flag() const;
        const std::string& details() const;
        const std::string& PolyType() const;
        const std::string& EvidenceCode() const;
        const std::string& chain_type() const;
        const std::string& prev_entity_id() const;
        const std::string& entity_id() const;
        const int int_prev_entity_id() const;
        const int int_entity_id() const;
        const std::string& entity_key();
        const int& ResidueNumbers() const;
        const bool& has_sequence() const;
        const bool& empty_chain() const;
        const bool& missing_sequence() const;
        const int& ca_or_p_atom_only() const;
        const bool& isModified() const;
        const int& index() const;
        const int& order() const;
        const int& na_type() const;
        const std::vector<std::map<std::string, std::string> >& get_linear_descriptors() const;
        const std::vector<std::map<std::string, std::string> >& get_branch_links() const;
        const unsigned int SeqLen() const;
        const unsigned int NumAtoms();
        const bool IsLastResidue() const;
        bool need_carbohydrate_annotation();
        void setCCDic(ConnectDic *ccdic);
        void setLog(LogUtil *logPt);
        void setMessage(MessageUtil *message);
        void set_ChainID(const std::string& a);
        void set_PDB_ChainID(const std::string& a);
        void set_PUB_ChainID(const std::string& a);
        void set_PDB_ChainID_Flag(const std::string& a);
        void set_details(const std::string& a);
        void set_PolyType(const std::string& a);
        void set_EvidenceCode(const std::string& a);
        void set_chain_type(const std::string& a);
        void set_missing_sequence();
        void set_index(const int& a);
        void set_order(const int& a);
        void set_na_type(const int& a);
        void set_isModified(const bool& a);
        void set_entity_id(const int& a);
        void set_entity_id(const std::string& a);
        void set_prev_entity_id(const std::string& a);
        void set_descriptor(const std::string& a);
        void set_entity_key(const std::string& a);
        void set_linear_descriptors(const std::vector<std::map<std::string, std::string> >& linear_descriptors);
        void set_branch_links(const std::vector<std::map<std::string, std::string> >& branch_links);
        void set_has_sequence();
        void setReserve(const int& number);
        RCSB::Residue* GetFirstResidue();
        RCSB::Residue* GetNextResidue();
        void GetFirstResidueList(std::vector<RCSB::Residue*>&);
        void GetNextResidueList(std::vector<RCSB::Residue*>&);
        void GetLastResidueList(std::vector<RCSB::Residue*>&);
        void GetResidueList(const int&, std::vector<RCSB::Residue*>&);
        void GetResidueListByIndex(const int&, std::vector<RCSB::Residue*>&);
        _FIELD* SeqRes(const int&);
        void merge_chain(Chain* chn);
        void merge_residue(RCSB::Residue* residue);
        void insert_a_residue(RCSB::Residue* residue, const bool& skip_alt_flag = false);
        void insert_hetero_residues(const std::vector<RCSB::Residue*>&);
        void insert_OXT2N_residue(const std::map<unsigned int, RCSB::Residue*>& residue_mapping);
        void remove_residue(const std::set<int>& residue_set);
        int find_field_index(const int&, const std::string&);
        int find_field_index(const std::string&, const std::string&);
        // RCSB::Residue* find_cif_residue(const std::string&, const int);
        // RCSB::Residue* find_cif_residue(const std::string&, const std::string&);
        RCSB::Residue* find_pdb_residue(const std::string&, const int, const std::string&);
        RCSB::Residue* find_pdb_residue(const std::string&, const std::string&, const std::string&);
        std::pair<int, RCSB::Residue*> find_pdb_residue_with_sort_index(const std::string&, const std::string&, const std::string&, const std::string&);
        RCSB::Residue* find_prev_residue(const std::string&, const int, const std::string&);
        RCSB::Residue* find_prev_residue(const std::string&, const std::string&, const std::string&);
        RCSB::Residue* find_next_residue(const std::string&, const std::string&, const std::string&);
        void insert_sequence(const std::vector<std::vector<std::string> >&);
        void insert_sequence(const SEQ&);
        bool seq_alignment(const std::vector<std::vector<std::string> >& mapping, const bool& rename_residue_flag = false);
        bool Mapping_SeqTool_Alignment(std::vector<std::vector<std::string> >&);
        void seq_alignment(const SEQ&, const std::map<std::string, std::pair<std::string, std::string> >&);
        bool check_seq_alignment(const std::vector<std::string>&, std::vector<std::string>&);
        std::string checking_missing_residue_in_coordinates();
        unsigned int count_missing_residue_in_coordinates();
        void update_missing_residues(const std::vector<std::vector<std::string> >&);
        void update_indices();
        void update_number();
        bool has_sequence_mismatch();
        void extract_alignment(std::vector<std::vector<std::string> >&);
        void check_polymer_sequence(const int&);
        std::string print_alignment(const int&);
        void update_asymId();
        void update_residues_nomenclature(const bool&, const bool&, int&, const bool& insertion_flag = false);
        void update_residues_nomenclature(const std::vector<std::string>&, const bool& insertion_flag = false);
        void set_PDB_ChainID_to_residues(const std::string& chain_id);
        void renumbering_polymer_chain(const bool& check_uniqueness_flag = false);
        void renumbering_polymer_chain(Chain *template_chain);
        void renumbering_polymer_chain_starting_with_one();
        void correction_residues_name();
        float weight();
        void get_seq(std::vector<std::string>& seqs);
        void get_seq(std::list<std::string>& seqs);
        void get_seq(std::list<std::vector<std::string> >& seqs);
        void get_seqres(const int&, const int&, std::list<std::vector<std::string> >&);
        std::string get_poly_type();
        void set_chain_type_to_residues();
        void check_ca_or_p_atom_only();
        void update_na_type();
        void get_descriptor(std::string&);
        void check_inscode(const int&);
        // void rename_dna_residues(const int&);
        void get_prev_entity_id();
        void get_split_chains(const std::vector<std::vector<std::string> >& splits, const std::vector<std::vector<std::string> >& chain_info,
                              int& chain_index, std::vector<Chain*>& split_chains);
        void get_nonpolymer_chains(const std::vector<std::vector<std::string> >& deletes, int& chain_index, std::vector<Chain*>& nonpolymer_chains);
        bool Update_SeqRes(const std::string& res_name, const std::string& res_num, const std::string& ins_code, const std::string& new_name);
        bool Update_SeqRes(const std::string& pdb_chnid, const std::string& res_name, const std::string& res_num, const std::string& ins_code,
                           const std::string& new_name);
        bool Update_SeqRes(const std::string& res_name, const std::string& res_num, const std::string&
                           ins_code, const std::map<std::string, std::string>& mapping);
        void update_nomenclature(const std::list<std::pair<int, std::string> >& idx_list, std::list<std::list<RCSB::Residue*> >& res_list);
        void merge_waters(const std::list<RCSB::Residue*>& res_list);
        bool Update(const std::string& resname, const std::string& resnum, const std::string& ins_code, const std::vector<RCSB::Residue*>& res_list);
        void reorder_residues(const std::map<int, std::string>& index_mapping);
        bool found_missing_single_residues();
        void check_missing_misplace_residues(const std::map<std::string, std::vector<RCSB::Chain*> >& res_name_chain_mapping,
                                             std::set<int>& removed_chain_index_set);
};

}

#endif

// error type:
// ins_code
