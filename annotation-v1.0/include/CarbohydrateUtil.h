/*
FILE:     CarbohydrateUtil.h
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
#ifndef _H_CARBOHYDRATE_UTIL_H_
#define _H_CARBOHYDRATE_UTIL_H_

#include <list>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "ConnectDic.h"
#include "Link.h"
#include "LogUtil.h"
#include "Residue.h"

class CarbohydrateUtil {
   public:

       /**
       **  Constructs CarbohydrateUtil
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
       CarbohydrateUtil();

       /**
       **  Destructs CarbohydrateUtil
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
       ~CarbohydrateUtil();

       /**
       **  Set LogUtil object
       **
       **  \param[in]: logPt - LogUtil object
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setLog(LogUtil* logPt);

       /**
       **  Set chemical component dictionary object
       **
       **  \param[in]: ccdic - Chemical component dictionary object
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setCCDic(ConnectDic* ccdic);

       void setPDBID(const std::string& id);

       /**
       **  Add a residue to oligosaccharide list
       **
       **  \param[in]: residue - input residue pointor
       **  \param[in]: type - residue type (ATOMP/ATOMS)
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void addResidue(const std::vector<RCSB::Residue*>& residues, const std::string& type = "ATOMS");

       /**
       **  Extract the linkage information from entry's linkage records and add to the sugar instance
       **
       **  \param[in]: links - linkage records for the entry
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void addLinks(const std::list<_LINK>& links);

       /**
       **  Build oligosaccharide structure.
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
       void buildOligoSaccharide();

       /**
       **  Get process status
       **
       **  \param: Not applicable
       **
       **  \return !_error_flag
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool getStatus();

       /**
       **  Get entity chemical descriptor.
       **
       **  \param: Not applicable
       **
       **  \return entity descriptor
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string getEntityDescriptor();

       /**
       **  Get condensed linear descriptor.
       **
       **  \param: Not applicable
       **
       **  \return the descriptor if exists
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string getCondensedDescriptor();

       /**
       **  Get new residue numbering index map.
       **
       **  \param: Not applicable
       **
       **  \return _new_index_map
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::map<int, std::string>& getReorderIndexMapping() const;

       /**
       **  Get additonal link from oligo saccharide.
       **
       **  \param: Not applicable
       **
       **  \return _additional_links
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::list<_LINK>& getAdditionalLinks() const;

       /**
       **  Get detected bad link list
       **
       **  \param: Not applicable
       **
       **  \return _bad_links
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::list<std::pair<RCSB::Atom*, RCSB::Atom*> >& getBadLinks() const;

       /**
       **  Get entity key
       **
       **  \param: Not applicable
       **
       **  \return entity key
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string getEntityKey();

   private:
       struct SugarLinkage {
            std::string _current_atom;
            unsigned int _current_number;
            unsigned int _linked_residue_index;
            std::string _linked_atom;
            unsigned int _linked_number;
            std::string _linkage_type;
            std::string _inverse_linkage_type;
       };

       struct SugarResidue {
            std::vector<RCSB::Residue*> _residues;
            RCSB::Residue* _attached_protein_residue;
            std::string _id;
            std::string _name;
            std::string _chemical_name;
            std::string _condensed_symbol;
            std::string _snfg_symbol;
            std::string _primary_carbonyl_group;
            std::string _anomeric_carbon;
            std::string _anomeric_oxygen;
            bool _is_visited;
            bool _is_root;
            int _new_index;

            std::vector<SugarLinkage> _reducing_linkage_list;
            std::vector<SugarLinkage> _non_reducing_linkage_list;
            // for anomeric to anomeric or other no standard reducing end linkages
            std::vector<SugarLinkage> _anomeric_anomeric_linkage_list;
       };

       struct OxygenTransfer {
            RCSB::Residue* src_residue;
            std::string src_carbon;
            std::string src_oxygen;
            RCSB::Residue* tgt_residue;
            std::string tgt_carbon;
            std::string tgt_oxygen;
       };

       bool _error_flag;
       bool _is_cyclic;
       std::string _pdb_id;
       LogUtil *_logIo;
       ConnectDic *_ccDic;
       RCSB::Residue* _protein_residue;
       std::vector<SugarResidue> _sugar_residue_list;
       std::map<std::string, unsigned int> _sugar_residue_index_map;

       std::set<std::string> _missing_info_ccd_id_set;

       // pair.first.first: first residue index
       // pair.first.second: first linked atom name
       // pair.second.first: second residue index
       // pair.second.second: second linked atom name
       std::vector<std::pair<std::pair<unsigned int, std::string>, std::pair<unsigned int, std::string> > > _link_list;

       std::list<_LINK> _additional_links;

       std::list<std::pair<RCSB::Atom*, RCSB::Atom*> > _bad_links;

       std::string _cyclic_descriptor;
       std::string _descriptor_in_full_form;
       std::string _descriptor_in_condensed_form;
       int _max_branch_length;
       std::set<std::string> _descriptor_id_set;
       int _global_index;

       std::map<int, std::string> _new_index_map;

       std::map<unsigned int, std::vector<unsigned int> > _linked_residue_index_mapping;

       std::list<OxygenTransfer> _oxygen_transfer_list;

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
       void _clear();

       /**
       **  Clear SugarResidue data storage
       **
       **  \param[out]: sugar_residue - reference to SugarResidue
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _clear_sugar_residue(SugarResidue& sugar_residue);

       /**
       **  Clear SugarLinkage data storage
       **
       **  \param[out]: sugar_linkage - reference to SugarLinkage
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _clear_sugar_linkage(SugarLinkage& sugar_linkage);

       /**
       **  Check if O1/O2 anomeric oxygen atom need to be moved to reducing end residue.
       **
       **  \param[in]: src_residue - source residue
       **  \param[in]: src_carbon - anomeric carbon name in source residue
       **  \param[in]: src_oxygen - oxygen atom name attached to anomeric carbon in source residue
       **  \param[in]: tgt_residue - target residue
       **  \param[in]: tgt_carbon - target carbon name for oxygen to attached.
       **  \param[in]: tgt_oxygen - target oxygen name
       **
       **  \return true if the moving operation is needed.
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _transfer_oxygen_atom(RCSB::Residue* src_residue, const std::string& src_carbon, const std::string& src_oxygen,
                                  RCSB::Residue* tgt_residue, const std::string& tgt_carbon, const std::string& tgt_oxygen);

       /**
       **  Move all possible O1/O2 anomeric oxygen atoms to the reducing end residue.
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
       void _transfer_oxygen_atom_list();

       /**
       **  Find the additional covalent bond(s) which are not defined struct_conn category.
       **
       **  \param[in]: i - first SugarResidue index
       **  \param[in]: j - second SugarResidue index
       **  \param[in]: type_i - first SugarResidue atom type
       **  \param[in]: type_j - second SugarResidue atom type
       **
       **  \return true if found the C-O bond
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _find_additional_links(unsigned int& i, unsigned int& j, const std::string& type_i, const std::string& type_j);

       /**
       **  Update SugarResidue._linked_residue_list based on _link_list
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
       void _update_linked_residue_list();

       /**
       **  Get number part from atom name.
       **
       **  \param[in]: res_name - residue name
       **  \param[in]: atom_name - atom name
       **
       **  \return number if it exists
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       int _get_number_part_from_atom_name(const std::string& res_name, const std::string& atom_name);

       /**
       **  Insert residue index into _linked_residue_index_mapping
       **
       **  \param[in]: first_index - first residue index
       **  \param[in]: second_index - second residue index
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _insert_linked_residue_index_mapping(const unsigned int& first_index, const unsigned int& second_index);

       /**
       **  Update cyclic saccharide information.
       **
       **  \param[out]: root - root index
       **
       **  \return true if successful.
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       bool _update_cyclic_saccharide_info(unsigned int& root);

       /**
       **  Traverse oligosaccharide graph.
       **
       **  \param[in]: idx - current SugarResidue index
       **  \param[out]: full_form_descriptor - full form pattern descriptor
       **  \param[out]: short_form_descriptor - condensed form pattern descriptor
       **  \param[out]: max_length - maximum sub branch length
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _traverse_graph(const unsigned int& idx, std::string& full_form_descriptor, std::string& short_form_descriptor, int& max_length);

       /**
       **  Re-index the residue numbering order based on the linear descriptor
       **
       **  \param[in]: descriptor - linear descriptor
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _reindex_residues(const std::string& descriptor);

       /**
       **  Generate _new_index_map.
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
       void _generate_new_index_map();

       /**
       **  Re-assign anomeric_anomeric_linkage_list to _non_reducing_linkage_list/_reducing_linkage_list.
       **
       **  \param[in]: idx - current SugarResidue index
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _reassign_anomeric_anomeric_linkage_list(const unsigned int& idx);

       /**
       **  Re-assign _reducing_linkage_list to _non_reducing_linkage_list/_reducing_linkage_list.
       **
       **  \param[in]: idx - current SugarResidue index
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _reassign_reducing_linkage_list(const unsigned int& idx);

       /**
       **  Re-order _non_reducing_linkage_list based on atom's order number
       **
       **  \param[in]: idx - current SugarResidue index
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void _reorder_non_reducing_linkage_list(const unsigned int& idx);

       /**
       **  Get reducing end terminal
       **
       **  \param[in]: ccd_id - chemical component ID
       **  \param[in]: anomeric_carbon - anomeric carbon atom name
       **
       **  \return the terminal name
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       std::string _get_reducing_terminal(const std::string& ccd_id, const std::string& anomeric_carbon);

};

#endif
