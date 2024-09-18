/*
FILE:     Entity.h
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
#ifndef _H_ENTITY_H_
#define _H_ENTITY_H_

#include <map>
#include <set>
#include <string>
#include <vector>

#define POLYMER_TEXT              "polymer"
#define BRANCHED_TEXT             "branched"
#define NON_POLYMER_TEXT          "non-polymer"
#define WATER_TEXT                "water"
#define DNA_TEXT                  "DNA"
#define RNA_TEXT                  "RNA"
#define PR_TEXT                   "PROTEIN"
#define T_RNA_TEXT                "T-RNA"
#define HYBRID_TEXT               "DNA/RNA"

class Entity
{
   private:
       // entity ID
       std::string _entity_id;

       // entity chain type
       std::string _chain_type;

       // number of copies
       int _number;

       // CompositeIndex of polymer seqs
       std::string _entity_key;

       // polymer sequence
       std::vector<std::string> _seqs;

       // Array of unique PDB chain IDs
       std::vector<std::string> _PDB_chainIDs;

       // Set of unique PDB chain IDs
       std::set<std::string> _PDB_chainID_set;

       // key/value pair to store data from entity, entity_keywords, entity_name_com,
       //     entity_name_sys and entity_poly categories where key is the cif item name
       //     except: sys_name for _entity_name_sys.name 
       //             com_name for _entity_name_com.name
       //             poly_type for _entity_poly.type
       std::map<std::string, std::string> _meta_data;

       // pair.first: source type: nat, man, syn
       // pair.second: key/value pair to store data from entity_src_nat, entity_src_gen
       //              and pdbx_entity_src_syn categories where key is the cif item name
       std::vector<std::pair<std::string, std::map<std::string, std::string> > > _source;

       // key: items in pdbx_entity_descriptor category
       std::vector<std::map<std::string, std::string> > _linear_descriptors;

       // key: items in pdbx_entity_branch_link category
       std::vector<std::map<std::string, std::string> > _entity_links;

       /**
       **  Update multiple source information
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
       void _updateSource();
   public:
       /**
       **  Constructs Entity
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
       Entity();

       /**
       **  Destructs Entity
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
       ~Entity();

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
       **  Check if entity has real data
       **
       **  \param: Not applicable
       **
       **  \return true  - if has _meta_data or _source data
       **          false - if has no _meta_data and _source data
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const bool empty() const;

       /**
       **  Copies a entity to another entity (assignment operator).
       **
       **  \param[in]: entity - reference to the source entity
       **
       **  \return Reference to the destination entity
       **
       **  \pre None
       **
       **  \post Constructed entity is a clone as the entity referenced by \e entity
       **
       **  \exception: None
       */
       Entity& operator=(const Entity& entity);

       /**
       **  Retrieves entity ID
       **
       **  \param: None
       **
       **  \return Constant reference to a string that contains entity ID.
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& entity_id() const;

       /**
       **  Retrieves chain type
       **
       **  \param: None
       **
       **  \return Constant reference to a string that contains chain type.
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& chain_type() const;

       /**
       **  Retrieves entity Key
       **
       **  \param: None
       **
       **  \return Constant reference to a string that contains entity Key.
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string& entity_key() const;

       /**
       **  Retrieves polymer sequence
       **
       **  \param: None
       **
       **  \return Constant reference to a vector that contains polymer sequence
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::string>& Seqs() const;

       /**
       **  Retrieves PDB chain IDs
       **
       **  \param: None
       **
       **  \return Constant reference to a vector that contains PDB chain IDs
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::string>& PDB_chainID() const;

       /**
       **  Retrieves source information
       **
       **  \param: None
       **
       **  \return Constant reference to a vector that contains source information
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::pair<std::string, std::map<std::string, std::string> > >& source() const;

       /**
       **  Retrieves linear descriptor information
       **
       **  \param: None
       **
       **  \return Constant reference to a vector that contains linear descriptor information
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::map<std::string, std::string> >& linear_descriptors() const;

       /**
       **  Retrieves intra-entity linkage information
       **
       **  \param: None
       **
       **  \return Constant reference to a vector that contains intra-entity linkage information
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::vector<std::map<std::string, std::string> >& entity_links() const;

       /**
       **  Get value based on key item name
       **
       **  \param[in]: key - key item name
       **
       **  \return the value corresponds to key item name
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string getValue(const std::string& key) const;

       /**
       **  Get first residue name from _seqs
       **
       **  \param: None
       **
       **  \return first residue name
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       const std::string getFirstResName() const;

       /**
       **  Set entity ID
       **
       **  \param[in]: id - entity ID
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setEntityID(const std::string& id);

       /**
       **  Set chain type
       **
       **  \param[in]: type - chain type
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setChainType(const std::string& type);

       /**
       **  Add number of entity copies
       **
       **  \param[in]: num - number of copies
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void addNumber(const int& num);

       /**
       **  Set sequence
       **
       **  \param[in]: seqs - sequence
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setSeqs(const std::vector<std::string>& seqs);

       /**
       **  Insert PDB chain ID
       **
       **  \param[in]: chnid - PDB chain ID
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void insertPDBChainID(const std::string& chnid);

       /**
       **  Insert source information
       **
       **  \param[in]: type - source type
       **  \param[in]: data - source data
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void insertSource(const std::string& type, const std::map<std::string, std::string>& data);

       /**
       **  Update source information
       **
       **  \param[in]: type - source type
       **  \param[in]: key  - item name
       **  \param[in]: val  - source value
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void updateSource(const std::string& type, const std::string& key, const std::string& val);

       /**
       **  Insert key/value pair to _meta_data
       **
       **  \param[in]: key - item key
       **  \param[in]: val - item value
       **  \param[in]: force_flag - force to update pdbx_description flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void insertValue(const std::string& key, const std::string& val, const bool& force_flag = false);

       /**
       **  Update key/value pair to _meta_data
       **
       **  \param[in]: key - item key
       **  \param[in]: val - item value
       **  \param[in]: attach_flag - attach value to existing value flag
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void updateValue(const std::string& key, const std::string& val, const bool& attach_flag = false);

       /**
       **  Delete key/value pair to _meta_data
       **
       **  \param[in]: key - item key
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void deleteValue(const std::string& key);

       /**
       **  Copy meta data, sources, linear descriptors and entity linkages information from another ENTITY
       **
       **  \param[in]: origin - original ENTITY
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void CopyAllMetaData(const Entity& origin);

       /**
       **  Copy meta data information from another ENTITY
       **
       **  \param[in]: origin - original ENTITY
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void CopyMetaData(const Entity& origin);

       /**
       **  merge source information from another ENTITY
       **
       **  \param[in]: origin - original ENTITY
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void MergeSource(const Entity& origin);

       /**
       **  merge source information from another ENTITY
       **
       **  \param[in]: origin - original ENTITY
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void MergeSource_1(const Entity& origin);

       /**
       **  Merge meta data information from another ENTITY
       **
       **  \param[in]: origin - original ENTITY
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void MergeMetaData(const Entity& origin, const bool& concat_flag = true);

       /**
       **  Update _entity.type, _entity.src_method, _entity.pdbx_number_of_molecules
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
       void update();

       /**
       **  Replace source information
       **
       **  \param[in]: new_sources - new source information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setSource(const std::vector<std::pair<std::string, std::map<std::string, std::string> > >& new_sources);

       /**
       **  Replace linear descriptor information
       **
       **  \param[in]: descriptors - new linear descriptor information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setLinearDescriptors(const std::vector<std::map<std::string, std::string> >& descriptors);

       /**
       **  Replace intra-entity linkage information
       **
       **  \param[in]: descriptors - new intra-entity linkage information
       **
       **  \return Not applicable
       **
       **  \pre None
       **
       **  \post None
       **
       **  \exception: None
       */
       void setEntityLinks(const std::vector<std::map<std::string, std::string> >& entity_links);
};

#endif
