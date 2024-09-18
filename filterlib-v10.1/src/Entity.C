/*
FILE:     Entity.C
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
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "CompositeIndex.h"
#include "Entity.h"
#include "utillib.h"

Entity::Entity()
{
       clear();
}

Entity::~Entity()
{
       clear();
}

void Entity::clear()
{
       _entity_id.clear();
       _chain_type.clear();
       _number = 0;
       _entity_key.clear();
       _seqs.clear();
       _PDB_chainIDs.clear();
       _PDB_chainID_set.clear();
       _meta_data.clear();
       _source.clear();
       _linear_descriptors.clear();
       _entity_links.clear();
}

const bool Entity::empty() const
{
       return (_meta_data.empty() && _source.empty());
}

void Entity::_updateSource()
{
       if (_source.size() < 2) return;

       // set _entity.src_method to "man"
       std::map<std::string, std::string>::iterator
           mpos = _meta_data.find("src_method");
       if (mpos != _meta_data.end())
            mpos->second = "man";
       else _meta_data.insert(std::make_pair("src_method", "man"));
}

Entity& Entity::operator=(const Entity& entity)
{
       if (this != &entity) {
            clear();
            _entity_id       = entity._entity_id;
            _entity_key      = entity._entity_key;
            _seqs            = entity._seqs;
            _PDB_chainIDs    = entity._PDB_chainIDs;
            _PDB_chainID_set = entity._PDB_chainID_set;
            _meta_data       = entity._meta_data;
            _source          = entity._source;
            _linear_descriptors = entity._linear_descriptors;
            _entity_links    = entity._entity_links;
       }
       return (*this);
}

const std::string& Entity::entity_id() const
{
       return _entity_id;
}

const std::string& Entity::chain_type() const
{
       return _chain_type;
}

const std::string& Entity::entity_key() const
{
       return _entity_key;
}

const std::vector<std::string>& Entity::Seqs() const
{
       return _seqs;
}

const std::vector<std::string>& Entity::PDB_chainID() const
{
       return _PDB_chainIDs;
}

const std::vector<std::pair<std::string, std::map<std::string, std::string> > >& Entity::source() const
{
       return _source;
}

const std::vector<std::map<std::string, std::string> >& Entity::linear_descriptors() const
{
       return _linear_descriptors;
}

const std::vector<std::map<std::string, std::string> >& Entity::entity_links() const
{
       return _entity_links;
}

const std::string Entity::getValue(const std::string& key) const
{
       std::map<std::string, std::string>::const_iterator mpos = _meta_data.find(key);
       if (mpos != _meta_data.end()) return mpos->second;
       return "";
}

const std::string Entity::getFirstResName() const
{
       if (!_seqs.empty()) return _seqs[0];
       return "";
}

void Entity::setEntityID(const std::string& id)
{
       _entity_id = id;
}

void Entity::setChainType(const std::string& type)
{
       _chain_type = type;
}

void Entity::addNumber(const int& num)
{
       _number += num;
}

void Entity::setSeqs(const std::vector<std::string>& seqs)
{
       if (seqs.empty()) return;

       _seqs = seqs;
       if (_seqs.size() == 1)
            _entity_key = _seqs[0];
       else _entity_key = CompositeIndex::getIndex(_seqs);
}

void Entity::insertPDBChainID(const std::string& chnid)
{
       if (_PDB_chainID_set.find(chnid) != _PDB_chainID_set.end()) return;

       _PDB_chainID_set.insert(chnid);
       _PDB_chainIDs.push_back(chnid);
}

void Entity::insertSource(const std::string& type, const std::map<std::string, std::string>& data)
{
       if (type.empty() || data.empty()) return;

       _source.push_back(std::make_pair(type, data));
       // _updateSource();
}

void Entity::updateSource(const std::string& type, const std::string& key, const std::string& val)
{
       for (std::vector<std::pair<std::string, std::map<std::string, std::string> > >::iterator pos = _source.begin(); pos != _source.end(); ++pos) {
            if (pos->first == type) {
                 std::map<std::string, std::string>::const_iterator mpos = pos->second.find(key);
                 // only insert new value, do not change existing value
                 if (mpos == pos->second.end()) pos->second.insert(std::make_pair(key, val));
                 return;
            }
       }

       std::map<std::string, std::string> tmp_map;
       tmp_map.clear();
       tmp_map.insert(std::make_pair(key, val));
       _source.push_back(std::make_pair(type, tmp_map));
}

void Entity::insertValue(const std::string& key, const std::string& val, const bool& force_flag)
{
       if (key.empty() || val.empty()) return;

       std::map<std::string, std::string>::iterator mpos = _meta_data.find(key);
       if (mpos != _meta_data.end() && key == "pdbx_description" && !force_flag) return;

       if (mpos != _meta_data.end()) mpos->second = val;
       else _meta_data.insert(std::make_pair(key, val));
}

void Entity::updateValue(const std::string& key, const std::string& val, const bool& attach_flag)
{
       if (key.empty() || val.empty()) return;

       std::map<std::string, std::string>::iterator mpos = _meta_data.find(key);
       if (mpos == _meta_data.end())
            _meta_data.insert(std::make_pair(key, val));
       else {
            if (attach_flag)
                 mpos->second += ", " + val;
            else mpos->second = val;
       }
}

void Entity::deleteValue(const std::string& key)
{
       if (key.empty() || _meta_data.empty()) return;

       if (_meta_data.find(key) == _meta_data.end()) return;

       _meta_data.erase(key);
}

void Entity::CopyAllMetaData(const Entity& origin)
{
       if (!origin._meta_data.empty()) _meta_data = origin._meta_data;
       if (!origin._source.empty()) _source = origin._source;
       if (!origin._linear_descriptors.empty()) _linear_descriptors = origin._linear_descriptors;
       if (!origin._entity_links.empty()) _entity_links = origin._entity_links;
}

void Entity::CopyMetaData(const Entity& origin)
{
       if (origin._meta_data.empty()) return;

       _meta_data = origin._meta_data;
}

void Entity::MergeSource(const Entity& origin)
{
       if (origin._source.empty()) return;

       if (_source.empty()) {
            _source = origin._source;
            return;
       }

       for (std::vector<std::pair<std::string, std::map<std::string, std::string> > >::
            const_iterator vpos = origin._source.begin(); vpos != origin._source.end(); ++vpos) {
            _source.push_back(*vpos);
       }
       // _updateSource();
       for (unsigned int i = 0; i < _source.size(); ++i) {
            _source[i].second["pdbx_src_id"] = String::IntToString(i + 1);
       }
}

void Entity::MergeSource_1(const Entity& origin)
{
       if (origin._source.empty()) return;

       if (_source.empty()) {
            _source = origin._source;
            return;
       }

       // only works for simply case
       if ((_source.size() != 1) || (origin._source.size() != 1)) return;

       for (std::map<std::string, std::string>::const_iterator mpos = origin._source[0].second.begin(); mpos != origin._source[0].second.end(); ++mpos) {
            std::map<std::string, std::string>::iterator mpos1 = _source[0].second.find(mpos->first);
            if (mpos1 != _source[0].second.end()) continue;
            _source[0].second.insert(std::make_pair(mpos->first, mpos->second));
       }
}

void Entity::MergeMetaData(const Entity& origin, const bool& concat_flag)
{
       if (origin._meta_data.empty()) return;
    
       if (_meta_data.empty()) {
            _meta_data = origin._meta_data;
       }

       std::string cs1, cs2;
       for (std::map<std::string, std::string>::const_iterator mpos = origin._meta_data.begin(); mpos != origin._meta_data.end(); ++mpos) {
            std::map<std::string, std::string>::iterator mpos1 = _meta_data.find(mpos->first);
            if (mpos1 == _meta_data.end())
                 _meta_data.insert(std::make_pair(mpos->first, mpos->second));
            else if (concat_flag) {
                 if (mpos->first == "formula_weight") continue;

                 String::UpperCase(mpos->second, cs1);
                 String::UpperCase(mpos1->second, cs2);
                 if (cs1 !=  cs2) {
                      mpos1->second += ", " + mpos->second;
                 }
            }
       }
}

void Entity::update()
{
       std::string type = NON_POLYMER_TEXT;
       std::string method = "man";

       if (_chain_type == "ATOMP" || _chain_type == "ATOMN") {
            type = POLYMER_TEXT;
            if (!_source.empty()) method = _source[0].first;
       } else if (_chain_type == "ATOMS" /* && _seqs.size() > 1 */ )
            type = BRANCHED_TEXT;
       else if (_seqs.size() == 1 && (_seqs[0] == "HOH" || _seqs[0] == "DOD")) {
            type = WATER_TEXT;
            method = "nat";
       } else if (_seqs.size() == 1 && _chain_type != "ATOMS") {
            method = "syn";
       }

       std::map<std::string, std::string>::iterator mpos = _meta_data.find("type");
       if (mpos != _meta_data.end())
            mpos->second = type;
       else _meta_data.insert(std::make_pair("type", type));

       mpos = _meta_data.find("pdbx_number_of_molecules");
       if (mpos != _meta_data.end())
            mpos->second = String::IntToString(_number);
       else _meta_data.insert(std::make_pair("pdbx_number_of_molecules", String::IntToString(_number)));

       mpos = _meta_data.find("src_method");
       if (mpos != _meta_data.end()) return;

       _meta_data.insert(std::make_pair("src_method", method));
}

void Entity::setSource(const std::vector<std::pair<std::string, std::map<std::string, std::string> > >& new_sources)
{
       _source = new_sources;
}

void Entity::setLinearDescriptors(const std::vector<std::map<std::string, std::string> >& descriptors)
{
       _linear_descriptors = descriptors;
}
void Entity::setEntityLinks(const std::vector<std::map<std::string, std::string> >& entity_links)
{
       _entity_links = entity_links;
}
