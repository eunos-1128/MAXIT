/*
FILE:     Assembly.C
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
#include <stdlib.h>

#include "Assembly.h"
#include "utillib.h"

const char *_num_mers[NUM_MER+1] = {
       "UNKNOWN",
       "MONOMERIC",
       "DIMERIC",
       "TRIMERIC",
       "TETRAMERIC",
       "PENTAMERIC",
       "HEXAMERIC",
       "HEPTAMERIC",
       "OCTAMERIC",
       "NONAMERIC",
       "DECAMERIC",
       "UNDECAMERIC",
       "DODECAMERIC",
       "TRIDECAMERIC",
       "TETRADECAMERIC",
       "PENTADECAMERIC",
       "HEXADECAMERIC",
       "HEPTADECAMERIC",
       "OCTADECAMERIC",
       "NONADECAMERIC",
       "EICOSAMERIC"
};

void Assembly::clear()
{
       _meta_data.clear();
       _assembly_chains.clear();
       _assembly_index.clear();
}

Assembly& Assembly::operator=(const Assembly& asmbly)
{
       if (this != &asmbly) {
            clear();
            _meta_data = asmbly._meta_data;
            _assembly_chains = asmbly._assembly_chains;
            _assembly_index = asmbly._assembly_index;
       }
       return (*this);
}

const bool Assembly::empty() const
{
       return (_meta_data.empty() && _assembly_chains.empty());
}

const std::string Assembly::getValue(const std::string& key) const
{
       std::map<std::string, std::string>::const_iterator pos = _meta_data.find(key);
       if (pos != _meta_data.end()) return pos->second;
       return "";
}

void Assembly::setValue(const std::string& key, const std::string& val)
{
       std::map<std::string, std::string>::iterator mpos = _meta_data.find(key);
       if (mpos != _meta_data.end())
            mpos->second = val;
       else _meta_data.insert(std::make_pair(key, val));
}

const std::vector<std::pair<std::string, std::vector<std::string> > >& Assembly::assembly_chains() const
{
       return _assembly_chains;
}

void Assembly::InsertChain(const std::string& oper_id, const std::string& chain_id)
{
       std::map<std::string, std::pair<unsigned int, std::set<std::string> > >::iterator mpos = _assembly_index.find(oper_id);
       if (mpos != _assembly_index.end()) {
            if (mpos->second.second.find(chain_id) == mpos->second.second.end()) {
                 mpos->second.second.insert(chain_id);
                 _assembly_chains[mpos->second.first].second.push_back(chain_id);
            }
       } else {
            std::vector<std::string> data_vector;
            data_vector.clear();
            data_vector.push_back(chain_id);
            std::set<std::string> data_set;
            data_set.clear();
            data_set.insert(chain_id);
            _assembly_index.insert(std::make_pair(oper_id, std::make_pair(_assembly_chains.size(), data_set)));
            _assembly_chains.push_back(std::make_pair(oper_id, data_vector));
       }
}

void Assembly::InsertChains(const std::string& oper_id, const std::vector<std::string>& chain_ids)
{
       for (std::vector<std::string>::const_iterator pos = chain_ids.begin(); pos != chain_ids.end(); ++pos) {
            InsertChain(oper_id, *pos);
       }
}

void Assembly::UpdateOligomericStatus(const std::set<std::string>& linear_polymeric_chain_id_set)
{
       if (_meta_data.find("oligomeric") != _meta_data.end() && _meta_data.find("details") != _meta_data.end()) return;

       _update_oligomeric_count(linear_polymeric_chain_id_set);

       if (_meta_data.find("details") == _meta_data.end()) UpdateOligomericDetails();
}

void Assembly::RenameReoderChainID(const std::map<std::string, int>& chain_order, const std::map<std::string, std::string>& chain_id_mapping)
{
       _assembly_index.clear();
       unsigned int count = 0;
       std::multimap<int, std::string> ordered_chains;
       std::vector<std::string> data_vector;
       std::set<std::string> data_set;
       for (std::vector<std::pair<std::string, std::vector<std::string> > >::iterator pos = _assembly_chains.begin(); pos != _assembly_chains.end(); ++pos) {
            ordered_chains.clear();
            int order = 0;
            for (std::vector<std::string>::const_iterator vpos = pos->second.begin(); vpos != pos->second.end(); ++vpos) {
                 std::string chain_id = *vpos;
                 // rename PDB chain ID if the original PDB chain ID has been changed
                 if (!chain_id_mapping.empty()) {
                      std::map<std::string, std::string>::const_iterator mpos = chain_id_mapping.find(chain_id);
                      if (mpos != chain_id_mapping.end()) chain_id = mpos->second;
                 }
                 if (!chain_order.empty()) {
                      std::map<std::string, int>::const_iterator sipos = chain_order.find(chain_id);
                      if (sipos == chain_order.end()) continue; // skip non existing PDB chain ID
                      ordered_chains.insert(std::make_pair(sipos->second, sipos->first));
                 } else ordered_chains.insert(std::make_pair(order, chain_id));
                 order++;
            }

            // rewrite PDB chain IDs
            data_vector.clear();
            data_set.clear();
            for (std::multimap<int, std::string>::const_iterator ispos = ordered_chains.begin(); ispos != ordered_chains.end(); ++ispos) {
                 data_vector.push_back(ispos->second);
                 data_set.insert(ispos->second);
            }
            pos->second = data_vector;
            _assembly_index.insert(std::make_pair(pos->first, std::make_pair(count, data_set)));
            count++;
       }
}

void Assembly::Merge()
{
       std::vector<std::vector<unsigned int> > index;
       index.clear();

       std::vector<unsigned int> tmp_vector;
       unsigned int count = 0;
       for (std::vector<std::pair<std::string, std::vector<std::string> > >::const_iterator
            pos = _assembly_chains.begin(); pos != _assembly_chains.end(); ++pos) {
            bool found = false;
            for (std::vector<std::vector<unsigned int> >::iterator ipos = index.begin(); ipos != index.end(); ++ipos) {
                 if (pos->second == _assembly_chains[(*ipos)[0]].second) {
                      found = true;
                      ipos->push_back(count);
                      break;
                 }
            }
            if (!found) {
                 tmp_vector.clear();
                 tmp_vector.push_back(count);
                 index.push_back(tmp_vector);
            }
            count++;
       }

       if (index.size() == _assembly_chains.size()) return;

       std::set<std::string> data_set;

       std::vector<std::pair<std::string, std::vector<std::string> > > tmp_assembly_chains;
       tmp_assembly_chains.clear();
       _assembly_index.clear();

       for (std::vector<std::vector<unsigned int> >::const_iterator ipos = index.begin(); ipos != index.end(); ++ipos) {
            std::string oper_id = "";
            for (std::vector<unsigned int>::const_iterator iipos = ipos->begin(); iipos != ipos->end(); ++iipos) {
                 if (!oper_id.empty()) oper_id += ",";
                 oper_id += _assembly_chains[*iipos].first;
            }
            data_set.clear();
            for (std::vector<std::string>::const_iterator vpos = _assembly_chains[(*ipos)[0]].second.begin();
                 vpos != _assembly_chains[(*ipos)[0]].second.end(); ++vpos) {
                 data_set.insert(*vpos);
            }
            _assembly_index.insert(std::make_pair(oper_id, std::make_pair(tmp_assembly_chains.size(), data_set))); 
            tmp_assembly_chains.push_back(std::make_pair(oper_id, _assembly_chains[(*ipos)[0]].second));
       }
       _assembly_chains = tmp_assembly_chains;
}

void Assembly::AddSugarChains(const std::map<std::string, std::vector<std::string> >& sugar_old_new_chain_id_mapping)
{
       if (sugar_old_new_chain_id_mapping.empty()) return;

       std::set<std::string> existing_chain_id_set;
       for (std::vector<std::pair<std::string, std::vector<std::string> > >::iterator pos = _assembly_chains.begin(); pos != _assembly_chains.end(); ++pos) {
            if (pos->second.empty()) continue;

            existing_chain_id_set.clear();
            for (std::vector<std::string>::const_iterator ppos = pos->second.begin(); ppos != pos->second.end(); ++ppos) {
                 existing_chain_id_set.insert(*ppos);
            }
            unsigned int size = pos->second.size();
            for (unsigned int i = 0; i < size; ++i) {
                 std::map<std::string, std::vector<std::string> >::const_iterator mpos = sugar_old_new_chain_id_mapping.find(pos->second[i]);
                 if (mpos == sugar_old_new_chain_id_mapping.end()) continue;
                 for (std::vector<std::string>::const_iterator ppos = mpos->second.begin(); ppos != mpos->second.end(); ++ppos) {
                      if (existing_chain_id_set.find(*ppos) != existing_chain_id_set.end()) continue;
                      pos->second.push_back(*ppos);
                 }
            }
       }
}

void Assembly::UpdateOligomericDetails()
{
       std::map<std::string, std::string>::iterator mpos = _meta_data.find("oligomeric");
       if (mpos == _meta_data.end()) return;

       int count = atoi(mpos->second.c_str());
       if (count == 0) return;

       std::string cs;
       if (count <= NUM_MER) {
            cs = _num_mers[count];
            String::LowerCase(cs);
       } else cs = String::IntToString(count) + "-meric";

       mpos = _meta_data.find("details");
       if (mpos != _meta_data.end())
            mpos->second = cs;
       else _meta_data.insert(std::make_pair("details", cs));
}

void Assembly::_update_oligomeric_count(const std::set<std::string>& linear_polymeric_chain_id_set)
{
       if (_meta_data.find("oligomeric") != _meta_data.end()) return;

       unsigned int count = 0;
       std::vector<std::string> oper_list;
       for (std::vector<std::pair<std::string, std::vector<std::string> > >::const_iterator
            pos = _assembly_chains.begin(); pos != _assembly_chains.end(); ++pos) {
            oper_list.clear();
            if (parseString(pos->first, oper_list) == 0) {
                 // count += oper_list.size() * pos->second.size();
                 unsigned int linear_polymeric_chain_count = 0;
                 for (std::vector<std::string>::const_iterator ppos = pos->second.begin(); ppos != pos->second.end(); ++ppos) {
                      if (linear_polymeric_chain_id_set.find(*ppos) != linear_polymeric_chain_id_set.end()) linear_polymeric_chain_count++;
                 }
                 count += oper_list.size() * linear_polymeric_chain_count;
            }
       }
       if (count) _meta_data.insert(std::make_pair("oligomeric", String::IntToString(count)));

       UpdateOligomericDetails();
}
