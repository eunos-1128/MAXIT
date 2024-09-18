/*
FILE:     CNA_Assign_Numbering.C
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

#include "ChainIDNumberAssignment.h"

void ChainIDNumberAssignment::_assign_numering(RCSB::Molecule* mol)
{
       _minumum_number = -9999;
       _non_polymer_numbering.clear();
       _unique_index.clear();
       _last_numbers.clear();

       std::vector<std::string> pdb_chain_ids;
       std::map<std::string, std::vector<RCSB::Chain*> > polymer_chains, non_polymer_chains, water_chains;
       mol->GetChains(pdb_chain_ids, polymer_chains, non_polymer_chains, water_chains);

       _checking_polymer_numbering(polymer_chains);
       // if (!_checking_unique_non_polymer_numbering(non_polymer_chains, water_chains))
            _checking_non_polymer_numbering(non_polymer_chains);
       _checking_non_polymer_numbering(water_chains);
}

void ChainIDNumberAssignment::_checking_polymer_numbering(std::map<std::string, std::vector<RCSB::Chain*> >& polymers)
{
       if (polymers.empty()) return;

       std::set<std::string> index;
       for (std::map<std::string, std::vector<RCSB::Chain*> >::iterator mpos = polymers.begin(); mpos != polymers.end(); ++mpos) {
            _minumum_number = -9999;
            for (std::vector<RCSB::Chain*>::iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 index.clear();
                 if (_found_redundant_number(index, *vpos)) (*vpos)->renumbering_polymer_chain();
                 _generate_index(*vpos);
            }
            _last_numbers.insert(std::make_pair(mpos->first, _minumum_number));
       }
}

bool ChainIDNumberAssignment::_checking_unique_non_polymer_numbering(std::map<std::string, std::vector<RCSB::Chain*> >& non_polymer_chains,
                                     std::map<std::string, std::vector<RCSB::Chain*> >& water_chains)
{
       if (non_polymer_chains.empty()) return true;

       std::map<std::string, int> copied_last_numbers = _last_numbers;

       std::set<std::string> index;
       for (std::map<std::string, std::vector<RCSB::Chain*> >::iterator mpos = non_polymer_chains.begin(); mpos != non_polymer_chains.end(); ++mpos) {
            _minumum_number = -9999;
            _non_polymer_numbering.clear();
            index.clear();
            for (std::vector<RCSB::Chain*>::iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 if (_found_redundant_number(index, *vpos)) return false;
                 _generate_numbering(*vpos);
            }
            std::map<std::string, int>::iterator pos = copied_last_numbers.find(mpos->first);
            if (pos != copied_last_numbers.end() && _minumum_number >  pos->second) pos->second = _minumum_number;
       }

       if (water_chains.empty()) return true;

       for (std::map<std::string, std::vector<RCSB::Chain*> >::iterator mpos = water_chains.begin(); mpos != water_chains.end(); ++mpos) {
            _minumum_number = -9999;
            _non_polymer_numbering.clear();
            for (std::vector<RCSB::Chain*>::iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 _generate_numbering(*vpos);
            }

            std::map<std::string, int>::iterator pos = copied_last_numbers.find(mpos->first);
            if (pos == copied_last_numbers.end()) continue;

            int num = (pos->second / 100 + 1) * 100;
            int start_num = num;
            if ((((_non_polymer_numbering[0] - 1) / 100) * 100 == (_non_polymer_numbering[0] - 1)) &&
                _non_polymer_numbering[0] >= num && _non_polymer_numbering[0] < 9999) start_num = _non_polymer_numbering[0] - 1;

            for (std::vector<RCSB::Chain*>::iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 start_num += (*vpos)->SeqLen();
            }
            if (start_num > 9999) return false;
       }

       return true;
}
                                                                    

void ChainIDNumberAssignment::_checking_non_polymer_numbering(std::map<std::string, std::vector<RCSB::Chain*> >& chains)
{
       if (chains.empty()) return;

       std::set<std::string> index;
       for (std::map<std::string, std::vector<RCSB::Chain*> >::iterator mpos = chains.begin(); mpos != chains.end(); ++mpos) {
            bool found = false;
            _minumum_number = -9999;
            _non_polymer_numbering.clear();
            index.clear();
            for (std::vector<RCSB::Chain*>::iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 if (_found_redundant_number(index, *vpos)) found = true;
                 _generate_numbering(*vpos);
            }
            if (_non_polymer_numbering.empty()) continue;

            std::map<std::string, int>::iterator pos = _last_numbers.find(mpos->first);
            if (pos == _last_numbers.end()) continue;
            int num = (pos->second / 100 + 1) * 100;

            if (((_non_polymer_numbering[0] / 100) * 100 != _non_polymer_numbering[0]) || _non_polymer_numbering[0] < num) found = true;
            if (!found) {
                 for (unsigned int i = 1; i < _non_polymer_numbering.size(); ++i) {
                      if (_non_polymer_numbering[i] != (_non_polymer_numbering[i - 1] + 1)) {
                           found = true;
                           break;
                      }
                 }
            }

            if (_minumum_number > pos->second) pos->second = _minumum_number;

            if (!found) continue;

            int start_num = num + 1;
            if ((((_non_polymer_numbering[0] - 1) / 100) * 100 == (_non_polymer_numbering[0] - 1)) &&
                _non_polymer_numbering[0] >= num && _non_polymer_numbering[0] < 9999) start_num = _non_polymer_numbering[0];

            for (std::vector<RCSB::Chain*>::iterator vpos = mpos->second.begin(); vpos != mpos->second.end(); ++vpos) {
                 (*vpos)->update_residues_nomenclature(false, true, start_num, true);
                 _generate_index(*vpos);
            }
            pos->second = start_num;
       }
}

void ChainIDNumberAssignment::_generate_index(RCSB::Chain* chain)
{
       std::map<std::string, std::set<std::string> >::iterator mpos = _unique_index.find(chain->PDB_ChainID());
       if (mpos == _unique_index.end()) {
            std::set<std::string> index_set;
            index_set.clear();
            _unique_index.insert(std::make_pair(chain->PDB_ChainID(), index_set));
            mpos = _unique_index.find(chain->PDB_ChainID());
       }

       for (unsigned int i = 0; i < chain->SeqLen(); ++i) {
            _FIELD *seq = chain->SeqRes(i);
            if (seq == NULL) continue;
            std::string key = seq->Field[4] + "_" + seq->InsCode;
            mpos->second.insert(key);
            int num = atoi(seq->Field[4].c_str());
            if (num > _minumum_number) _minumum_number = num;
       }
}

void ChainIDNumberAssignment::_generate_numbering(RCSB::Chain* chain)
{
       for (unsigned int i = 0; i < chain->SeqLen(); ++i) {
            _FIELD *seq = chain->SeqRes(i);
            if (seq == NULL) continue;
            int num = atoi(seq->Field[4].c_str());
            if (num > _minumum_number) _minumum_number = num;
            _non_polymer_numbering.push_back(num);
       }
}
