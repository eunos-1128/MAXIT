/*
FILE:     Pdb2Ndb_Header.C
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Maxit.h"
#include "Source_global.h"
#include "utillib.h"

void Maxit::pdb_to_ndb()
{
       _pdb_to_ndb_process_header();
       _pdb_to_ndb_process_compnd_and_source();
       _pdb_to_ndb_process_remarks();
       if (_experiment_type & EXPERIMENT_TYPE_NMR ||
           _experiment_type & EXPERIMENT_TYPE_NMR_SOLID) {
            _pdb_to_ndb_update_SFTWAR();
            _pdb_to_ndb_update_NMRSDT_and_NMREXP();
            _pdb_to_ndb_update_NMRSPM();
            _pdb_to_ndb_update_NMRSMP();
       }
       _pdb_to_ndb_update_EMSFTW();
       _pdb_to_ndb_update_EMPIXL();
       _pdb_to_ndb_postprocessing();
}

void Maxit::_pdb_to_ndb_process_header()
{
       std::map<std::string, std::list<std::vector<std::string> > >::iterator
           ppos = _pdb_records.find("HEADER");
       if (ppos == _pdb_records.end()) return;

       std::vector<std::string>& Field = ppos->second.front();
       String::StripAndCompressWs(Field[3]);
       if (Field[3].empty()) return;

       if (Field[3].find(" ") != std::string::npos) {
            Field[3].clear();
            return;
       }
 
       String::UpperCase(Field[3]);

       _updateRecordFront("PDBFIL", 1, Field[3], true);
       _updateRecordFront("NDBFIL", 1, Field[3], true);

       _StructureId = Field[3];
  
       if (!Field[2].empty()) _updateRecordFront("STATUS", 6, Field[2], true);
}

void Maxit::_pdb_to_ndb_process_compnd_and_source()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator
           ppos = _pdb_records.find("COMPND");
       if (ppos == _pdb_records.end()) return;

       // first key: mol_ID as integer
       // second key: pdbx_description (molecule), chains, pdbx_fragment,
       //             com_name (synonym), pdbx_ec, engineered, pdbx_mutation,
       //             biol_unit & details
       std::map<int, std::map<std::string, std::string> > mol_compnds;

       _pdb_to_ndb_process_compnd(ppos->second, mol_compnds);

       // remove COMPND records
       _pdb_records.erase("COMPND");

       if (mol_compnds.empty()) return;

       ppos = _pdb_records.find("SOURCE");
       if (ppos == _pdb_records.end()) return;

       // first key: mol_ID as integer
       std::map<int, std::map<std::string, std::string> > mol_sources;
       _pdb_to_ndb_process_source(ppos->second, mol_sources);

       // remove SOURCE records
       _pdb_records.erase("SOURCE");

       _pdb_to_ndb_update_entities(mol_compnds, mol_sources);
}

void Maxit::_pdb_to_ndb_process_compnd(const std::list<std::vector<std::string> >& rlist,
                        std::map<int, std::map<std::string, std::string> >& mol_compnds)
{
       mol_compnds.clear();

       std::map<std::string, std::pair<std::string, int> > key_mapping;
       key_mapping.clear();
       key_mapping.insert(std::make_pair("MOLECULE:", std::make_pair("pdbx_description", 9)));
       key_mapping.insert(std::make_pair("CHAIN:", std::make_pair("chains", 6)));
       key_mapping.insert(std::make_pair("FRAGMENT:", std::make_pair("pdbx_fragment", 9)));
       key_mapping.insert(std::make_pair("SYNONYM:", std::make_pair("com_name", 8)));
       key_mapping.insert(std::make_pair("EC:", std::make_pair("pdbx_ec", 3)));
       key_mapping.insert(std::make_pair("ENGINEERED:", std::make_pair("engineered", 11)));
       key_mapping.insert(std::make_pair("MUTATION:", std::make_pair("pdbx_mutation", 9)));
       key_mapping.insert(std::make_pair("BIOLOGICAL_UNIT:", std::make_pair("biol_unit", 16)));
       key_mapping.insert(std::make_pair("OTHER_DETAILS:", std::make_pair("details", 14)));

       std::list<std::vector<std::string> > tmp_list;
       std::map<std::string, std::string> tmp_map;

       int mol_id = -1;
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = rlist.begin(); lpos != rlist.end(); ++lpos) {

            // copy COMPND records to PCOMPN records
            std::map<std::string, std::list<std::vector<std::string> > >::iterator
                qpos = _pdb_records.find("PCOMPN");
            if (qpos == _pdb_records.end()) {
                 tmp_list.clear();
                 tmp_list.push_back(*lpos);
                 _pdb_records.insert(std::make_pair("PCOMPN", tmp_list));
            } else qpos->second.push_back(*lpos);
            
            if ((*lpos)[2].substr(0, 7) == "MOL_ID:") {
                 mol_id = atoi((*lpos)[2].substr(7).c_str());
            } else {
                 if (mol_id < 0) continue;
                 for (std::map<std::string, std::pair<std::string, int> >::const_iterator
                      kpos = key_mapping.begin(); kpos != key_mapping.end(); ++kpos) {
                      if ((*lpos)[2].substr(0, kpos->second.second) == kpos->first) {
                           std::string cs = (*lpos)[2].substr(kpos->second.second);
                           _pdb_to_ndb_process_record_value(cs, true);
                           if (cs.empty()) break;
                           if (kpos->first == "MOLECULE:") {
                                if (cs.size() > 9 && cs.substr(0, 9) == "PROTEIN (" &&
                                    cs[cs.size() - 1] == ')') {
                                     cs.erase(cs.size() - 1);
                                     cs.erase(0, 9);
                                     String::StripAndCompressWs(cs);
                                }
                           }
                           std::map<int, std::map<std::string, std::string> >::iterator
                               mpos = mol_compnds.find(mol_id);
                           if  (mpos == mol_compnds.end()) {
                                 tmp_map.clear();
                                 tmp_map.insert(std::make_pair(kpos->second.first, cs));
                                 mol_compnds.insert(std::make_pair(mol_id, tmp_map));
                           } else mpos->second.insert(std::make_pair(kpos->second.first, cs));
                           break;
                      }
                 }
            }
       }
}

void Maxit::_pdb_to_ndb_process_source(const std::list<std::vector<std::string> >& rlist,
                        std::map<int, std::map<std::string, std::string> >& mol_sources)
{
       mol_sources.clear();

       std::list<std::vector<std::string> > tmp_list;
       std::map<std::string, std::string> tmp_map;

       int mol_id = -1;
       for (std::list<std::vector<std::string> >::const_iterator
            lpos = rlist.begin(); lpos != rlist.end(); ++lpos) {

            // copy SOURCE records to PSOURC records
            std::map<std::string, std::list<std::vector<std::string> > >::iterator
                qpos = _pdb_records.find("PSOURC");
            if (qpos == _pdb_records.end()) {
                 tmp_list.clear();
                 tmp_list.push_back(*lpos);
                 _pdb_records.insert(std::make_pair("PSOURC", tmp_list));
            } else qpos->second.push_back(*lpos);

            if ((*lpos)[2].substr(0, 7) == "MOL_ID:") {
                 mol_id = atoi((*lpos)[2].substr(7).c_str());
            } else {
                 if (mol_id < 0) continue;
                 for (int i = 0; i < NUM_PDB_SOURCE; ++i) {
                      int len = strlen(_pdb_source_tokens[i][0]);
                      if ((*lpos)[2].substr(0, len) != _pdb_source_tokens[i][0]) continue;
                      std::string cs = (*lpos)[2].substr(len);
                      _pdb_to_ndb_process_record_value(cs, true);
                      if (cs.empty()) break;
                      std::map<int, std::map<std::string, std::string> >::iterator
                          mpos = mol_sources.find(mol_id);
                      if  (mpos == mol_sources.end()) {
                            tmp_map.clear();
                            tmp_map.insert(std::make_pair(_pdb_source_tokens[i][1], cs));
                            mol_sources.insert(std::make_pair(mol_id, tmp_map));
                      } else mpos->second.insert(std::make_pair(_pdb_source_tokens[i][1], cs));
                      break;
                 }
            }
       }
}

void Maxit::_pdb_to_ndb_process_record_value(std::string& value, const bool& remove_flag)
{
       if (value.empty()) return;

       if (remove_flag && value[value.size() - 1] == ';') value.erase(value.size() - 1);
       String::StripAndCompressWs(value);
       if (value == "NULL" || value == "N/A" || value == "NONE") value.clear();
}

void Maxit::_pdb_to_ndb_update_entities(
                 const std::map<int, std::map<std::string, std::string> >& mol_compnds,
                 const std::map<int, std::map<std::string, std::string> >& mol_sources)
{
       std::map<std::string, int> chain_id_entity_id_mapping;
       chain_id_entity_id_mapping.clear();
       for (std::map<int, Entity>::const_iterator
            ipos = _entities.begin(); ipos != _entities.end(); ++ipos) {
            const std::vector<std::string>& PDB_chainIDs = ipos->second.PDB_chainID();
            for (std::vector<std::string>::const_iterator
                 vpos = PDB_chainIDs.begin(); vpos != PDB_chainIDs.end(); ++vpos) {
                 chain_id_entity_id_mapping.insert(std::make_pair(*vpos, ipos->first));
            }
       }

       std::set<int> updated_set, new_entity_set;
       updated_set.clear();

       // entity item tokens
       std::vector<std::string> items;
       items.clear();
       items.push_back("pdbx_description");
       items.push_back("pdbx_fragment");
       items.push_back("com_name");
       items.push_back("pdbx_ec");
       items.push_back("pdbx_mutation");
       items.push_back("details");

       std::vector<std::string> chain_ids;
       std::map<std::string, std::string> metadata, source;

       for (std::map<int, std::map<std::string, std::string> >::const_iterator
            cpos = mol_compnds.begin(); cpos != mol_compnds.end(); ++cpos) {
            std::map<std::string, std::string>::const_iterator
                mpos = cpos->second.find("chains");
            if (mpos == cpos->second.end()) continue;

            get_wordarray(chain_ids, mpos->second, ", ");
            if (chain_ids.empty()) continue;
     
            new_entity_set.clear();
            for (std::vector<std::string>::const_iterator
                 vpos = chain_ids.begin(); vpos != chain_ids.end(); ++vpos) {
                 std::map<std::string, int>::const_iterator
                     mpos1 = chain_id_entity_id_mapping.find(*vpos);
                 if (mpos1 == chain_id_entity_id_mapping.end()) continue;
                 if (updated_set.find(mpos1->second) != updated_set.end()) continue;
                 updated_set.insert(mpos1->second);
                 new_entity_set.insert(mpos1->second);
            }
            if (new_entity_set.empty()) continue;

            metadata.clear();
            source.clear();
            
            std::string type = "nat";
            mpos = cpos->second.find("engineered");
            if (mpos != cpos->second.end() && mpos->second == "YES")
                 type = "man";
            std::map<int, std::map<std::string, std::string> >::const_iterator
                spos = mol_sources.find(cpos->first);
            if (spos != mol_sources.end()) {
                 mpos = spos->second.find("source_type");
                 if (mpos != spos->second.end() && mpos->second == "YES") type = "syn";

                 source = spos->second;
                 source.erase("source_type");
            }

            metadata.insert(std::make_pair("src_method", type));
            for (std::vector<std::string>::const_iterator
                 vpos = items.begin(); vpos != items.end(); ++vpos) {
                 mpos = cpos->second.find(*vpos);
                 if (mpos == cpos->second.end()) continue;
                 metadata.insert(std::make_pair(*vpos, mpos->second));
            }

            _pdb_to_ndb_update_entities(new_entity_set, metadata, source);
       }
}

void Maxit::_pdb_to_ndb_update_entities(const std::set<int>& entity_set,
                               const std::map<std::string, std::string>& metadata,
                               const std::map<std::string, std::string>& source)
{
       std::string type = "";
       std::map<std::string, std::string>::const_iterator
           mpos = metadata.find("src_method");
       if (mpos != metadata.end()) type = mpos->second;

       for (std::set<int>::const_iterator
            spos = entity_set.begin(); spos != entity_set.end(); ++spos) {
            std::map<int, Entity>::iterator ipos = _entities.find(*spos);
            if (ipos == _entities.end()) continue;
            if (!source.empty() && !type.empty()) ipos->second.insertSource(type, source);
            for (mpos = metadata.begin(); mpos != metadata.end(); ++mpos) {
                 ipos->second.insertValue(mpos->first, mpos->second, true);
            }
       }
}
