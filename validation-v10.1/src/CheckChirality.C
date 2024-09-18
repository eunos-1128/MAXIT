/*
FILE:     CheckChirality.C
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
#include <math.h>

#include "BondUtil.h"
#include "ChiralVolumeUtil.h"
#include "CompositeIndex.h"
#include "GetPairList.h"
#include "GetPairList.C"
#include "ValidateObj.h"

static bool is_linked(const std::vector<RCSB::Atom*>& a_list, const std::vector<RCSB::Atom*>& b_list);

void ValidateObj::_checkChirality(RCSB::Residue* residue, RCSB::Molecule* mol)
{
       try {
            const ConnectFormat& drug = _ccDic->find_drug(residue->ResName());
            if (!drug.has_chiral_center()) return;

            std::vector<RCSB::Atom*> atom_list;
            std::map<std::string, std::vector<RCSB::Atom*> > link_mapping;
            std::vector<std::vector<RCSB::Atom*> > atom_lists, pair_lists;
            std::vector<std::vector<double> > coords, found_coords;
            std::vector<std::string> atomnames;
            double ang;
            Value value;

            std::string cs;
            link_mapping.clear();
            for (std::list<_LINK>::const_iterator lpos = _links.begin(); lpos != _links.end(); ++lpos) {
                 // if (lpos->mol_index != mol->index()) continue;
                 if (!lpos->SymOP_1.empty() && lpos->SymOP_1 != "1_555") continue;
                 if (!lpos->SymOP_2.empty() && lpos->SymOP_2 != "1_555") continue;

                 cs = CompositeIndex::getIndex(lpos->fstAtom->pdb_chnid(), lpos->fstAtom->pdb_resnum(),
                                               lpos->fstAtom->pdb_resnam(), lpos->fstAtom->pdb_atmnam());
                 std::map<std::string, std::vector<RCSB::Atom*> >::iterator mpos = link_mapping.find(cs);
                 if (mpos != link_mapping.end())
                      mpos->second.push_back(lpos->sndAtom);
                 else {
                      atom_list.clear();
                      atom_list.push_back(lpos->sndAtom);
                      link_mapping.insert(std::make_pair(cs, atom_list));
                 }

                 cs = CompositeIndex::getIndex(lpos->sndAtom->pdb_chnid(), lpos->sndAtom->pdb_resnum(),
                                               lpos->sndAtom->pdb_resnam(), lpos->sndAtom->pdb_atmnam());
                 mpos = link_mapping.find(cs);
                 if (mpos != link_mapping.end())
                      mpos->second.push_back(lpos->fstAtom);
                 else {
                      atom_list.clear();
                      atom_list.push_back(lpos->fstAtom);
                      link_mapping.insert(std::make_pair(cs, atom_list));
                 }
            }

            const std::vector<AtomFormat>& atoms = drug.atoms();
            for (std::vector<AtomFormat>::const_iterator apos = atoms.begin(); apos != atoms.end(); ++apos) {
                 if (apos->stereo_config().empty() || apos->stereo_config() == "N" ||
                    (apos->atomtype() != "C" && apos->atomtype() != "N")) continue;

                 drug.getChiralCenterAtoms(apos->atomname(), atomnames, coords);
                 if (atomnames.size() < 4) continue;

                 residue->find_atom(atomnames[0], atom_list);
                 if (atom_list.empty()) continue;

                 atom_lists.clear();
                 atom_lists.push_back(atom_list);
                 found_coords.clear();
                 found_coords.push_back(coords[0]);
                 bool used_linked_atom = false;
                 for (unsigned int j = 1; j < atomnames.size(); ++j) {
                      residue->find_atom(atomnames[j], atom_list);
                      if (atom_list.empty() && !used_linked_atom) {
                           cs = CompositeIndex::getIndex(atom_lists[0][0]->pdb_chnid(), atom_lists[0][0]->pdb_resnum(),
                                                         atom_lists[0][0]->pdb_resnam(), atom_lists[0][0]->pdb_atmnam());
                           std::map<std::string, std::vector<RCSB::Atom*> >::iterator mpos = link_mapping.find(cs);
                           if (mpos != link_mapping.end() && is_linked(atom_lists[0], mpos->second)) {
                                atom_list = mpos->second;
                                used_linked_atom = true;
                           }
                      }
                      if (atom_list.empty()) continue;
                      atom_lists.push_back(atom_list);
                      found_coords.push_back(coords[j]);
                      if (atom_lists.size() == 4) break;
                 }
                 if (atom_lists.size() != 4) continue;

                 GetPairList(atom_lists, pair_lists);
                 if (pair_lists.empty()) continue;

                 double value_ref = calc_chiral_volume(found_coords);

                 for (std::vector<std::vector<RCSB::Atom*> >::const_iterator pos = pair_lists.begin(); pos != pair_lists.end(); ++pos) {
                      double volume = calc_chiral_volume(*pos, ang);
                      if (volume * value_ref >= 0.0 && fabs(volume) > 0.7) continue;

                      value.clear();
                      value.set_Mol_ID(mol->Mol_ID());
                      value.set_val(volume);
                      if (fabs(volume) <= 0.7)
                           value.set_sval("PLANAR");
                      else value.set_sval("WRONG HAND");
                      for (unsigned int j = 0; j < pos->size(); ++j) {
                           value.insert_atom((*pos)[j]);
                      }
                      _violated_chiral_center.push_back(value);
                 }
            }
       } catch (const std::exception& exc) {}
}

static bool is_linked(const std::vector<RCSB::Atom*>& a_list, const std::vector<RCSB::Atom*>& b_list)
{
       for (std::vector<RCSB::Atom*>::const_iterator apos = a_list.begin(); apos != a_list.end(); ++apos) {
            for (std::vector<RCSB::Atom*>::const_iterator bpos = b_list.begin(); bpos != b_list.end(); ++bpos) {
                 if (!(*apos)->alt_loc().empty() && !(*bpos)->alt_loc().empty() && (*apos)->alt_loc() != (*bpos)->alt_loc()) continue;
                 if (BondUtil::is_a_bond(*apos, *bpos)) return true;
            }
       }
       return false;
}
