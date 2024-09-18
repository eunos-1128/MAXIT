/*
FILE:     CNA_Reposition_Water.C
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
#include <math.h>

#include "ChainIDNumberAssignment.h"

void ChainIDNumberAssignment::__reposition_waters(xtal& mycrys, RCSB::Molecule* mol, const std::list<RCSB::Residue*>& water_lists)
{
       std::vector<_WATER_RESIDUE> water_residues;
       water_residues.clear();

       if (water_lists.empty()) return;

       mycrys.markAtoms("non_branch_polar_atom");

       _get_water_residues(water_residues, mycrys, water_lists);

       grid *gr = mycrys.crysGrid();
       int elts_x = gr->eltsSide_x();
       int elts_y = gr->eltsSide_y();
       int elts_z = gr->eltsSide_z();

       double cr[3];
       cr[0] = mycrys.aval() / elts_x;
       cr[1] = mycrys.bval() / elts_y;
       cr[2] = mycrys.cval() / elts_z;

       for (std::vector<_WATER_RESIDUE>::iterator pos = water_residues.begin(); pos != water_residues.end(); ++pos) {
            double lo = 0;
            double hi = _water_cutoff_distance;
            int count = 0;
            while (!_found_contact(*pos, cr, lo, hi)) {
                 lo += _water_cutoff_distance;
                 hi += _water_cutoff_distance;
                 count++;
                 if (count > 30) break;
            }
       }

       mycrys.unmarkAtoms();

       if (_multi_model_chain_flag)
            _mol_water_residues.push_back(std::make_pair(mol, water_residues));
       else __reposition_waters(mol, water_residues);
}

void ChainIDNumberAssignment::_get_water_residues(std::vector<_WATER_RESIDUE>& water_residues, xtal& mycrys, const std::list<RCSB::Residue*>& water_lists)
{
       water_residues.clear();
       water_residues.reserve(water_lists.size());

       std::vector<symop> oplist;
       mycrys.getSymopList(oplist);

       cryatom newat;
       newat.setcrys(&mycrys);

       _WATER_RESIDUE water_residue;
       for (std::list<RCSB::Residue*>::const_iterator pos = water_lists.begin(); pos != water_lists.end(); ++pos) {
            RCSB::Atom* atom = (*pos)->find_atom("O");
            if (atom == NULL) continue;

            newat.clear();
            newat.setatom(atom, 100000, (*pos)->position(), (*pos)->index(), WATER_TYPE);

            water_residue.moving_flag = false;
            water_residue.atom = atom;
            water_residue.ref_atom = NULL;
            water_residue.residue = *pos;
            water_residue.cryatoms.clear();
            water_residue.cryatoms.reserve(oplist.size()+1);
            water_residue.pdb_chnid = atom->pdb_chnid();
            water_residue.dist = 1000.0;
            water_residue.sym = 0;
            water_residue.lx = 0;
            water_residue.ly = 0;
            water_residue.lz = 0;

            for (std::vector<symop>::const_iterator spos = oplist.begin(); spos != oplist.end(); ++spos) {
                 cryatom symat = newat.symrel(*spos);
                 symat.moveToCenterCell();
                 symat.setGridCoords();
                 water_residue.cryatoms.push_back(symat);
            }
            newat.moveToCenterCell();
            newat.setGridCoords();
            water_residue.cryatoms.push_back(newat);
            water_residues.push_back(water_residue);
       }
}

bool ChainIDNumberAssignment::_found_contact(_WATER_RESIDUE& water_residue, double *cr, const double& lo, const double &hi)
{
       int cells[3];
       for (int i = 0; i < 3; ++i) {
            cells[i] = int (floor (1 + (hi/(cr[i]))));
       }

       std::multimap<double, CONTACT> asym_contacts, sym_contacts;
       asym_contacts.clear();
       sym_contacts.clear();

       double low = lo - 0.001;
       if (low < 0.5) low = 0.5;
       std::list<CONTACT> found_list;
       for (std::vector<cryatom>::iterator pos = water_residue.cryatoms.begin(); pos != water_residue.cryatoms.end(); ++pos) {
            found_list.clear();
            pos->neighborAtoms(found_list, low, hi, cells, 0);
            if (found_list.empty()) continue;
            for (std::list<CONTACT>::iterator lpos = found_list.begin(); lpos != found_list.end(); ++lpos) {
                 lpos->sym = pos->symnum();
                 if (lpos->sym || lpos->lx || lpos->ly || lpos->lz)
                      sym_contacts.insert(std::make_pair(lpos->dist, *lpos));
                 else asym_contacts.insert(std::make_pair(lpos->dist, *lpos));
            }
       }

       if (asym_contacts.empty() && sym_contacts.empty()) return false;

       if (sym_contacts.empty()) {
            CONTACT contact_a = asym_contacts.begin()->second;
            water_residue.ref_atom = contact_a.b_atm;
            water_residue.pdb_chnid = contact_a.b_atm->pdb_chnid();
            water_residue.dist = contact_a.dist;
            return true;
       }

       CONTACT contact_s = sym_contacts.begin()->second;
       if (contact_s.dist > _water_cutoff_distance) {
            if (asym_contacts.empty()) return false;

            CONTACT contact_a = asym_contacts.begin()->second;
            water_residue.ref_atom = contact_a.b_atm;
            water_residue.pdb_chnid = contact_a.b_atm->pdb_chnid();
            water_residue.dist = contact_a.dist;
            return true;
       }

       bool status = false;
       if (asym_contacts.empty())
            status = _moving_water(water_residue, contact_s);
       else {
            CONTACT contact_a = asym_contacts.begin()->second;
            water_residue.ref_atom = contact_a.b_atm;
            water_residue.pdb_chnid = contact_a.b_atm->pdb_chnid();
            water_residue.dist = contact_a.dist;
            if (contact_s.dist < contact_a.dist && fabs(contact_a.dist - contact_s.dist) > 0.001) {
                 _moving_water(water_residue, contact_s);
            }
            status = true;
       }
       return status;
}

void ChainIDNumberAssignment::__reposition_waters(RCSB::Molecule* mol, std::vector<_WATER_RESIDUE>& water_residues)
{
       std::map<std::string, std::multimap<double, RCSB::Residue*> > ordered_waters;
       ordered_waters.clear();
       for (std::vector<_WATER_RESIDUE>::iterator pos = water_residues.begin(); pos != water_residues.end(); ++pos) {
            pos->residue->set_pdb_chnid(pos->pdb_chnid);
            RCSB::Atom* atom = pos->residue->GetFirstAtom();
            while (atom) {
                 atom->set_chnid(pos->pdb_chnid);
                 atom->set_pdb_chnid(pos->pdb_chnid);
                 atom = pos->residue->GetNextAtom();
            }

            std::map<std::string, std::multimap<double, RCSB::Residue*> >::iterator mpos = ordered_waters.find(pos->pdb_chnid);
            if (mpos != ordered_waters.end())
                 mpos->second.insert(std::make_pair(pos->dist, pos->residue));
            else {
                 std::multimap<double, RCSB::Residue*> t_map;
                 t_map.clear();
                 t_map.insert(std::make_pair(pos->dist, pos->residue));
                 ordered_waters.insert(std::make_pair(pos->pdb_chnid, t_map));
            }
       }

       if (ordered_waters.empty()) return;

       mol->ClearNonAsymAndWaterChains();
       for (std::map<std::string, std::multimap<double, RCSB::Residue*> >::const_iterator mpos = ordered_waters.begin(); mpos != ordered_waters.end(); ++mpos) {
            RCSB::Chain* chain = new RCSB::Chain;
            std::multimap<double, RCSB::Residue*>::const_iterator mmpos = mpos->second.begin();
            RCSB::Atom* atom = mmpos->second->GetFirstAtom();
            chain->set_ChainID(atom->chnid());
            chain->set_PDB_ChainID(atom->pdb_chnid());
            chain->set_PUB_ChainID(atom->pub_chnid());
            chain->set_chain_type("HETAS");
            for (mmpos = mpos->second.begin(); mmpos != mpos->second.end(); ++mmpos) {
                 chain->insert_a_residue(mmpos->second, true);
            }
            mol->insert_a_chain(chain);
       }
}

bool ChainIDNumberAssignment::_moving_water(_WATER_RESIDUE& water_residue, const CONTACT& contact)
{
       std::vector<RCSB::Atom*> atom_list;
       atom_list.clear();
       RCSB::Atom* atom = water_residue.residue->GetFirstAtom();
       while (atom) {
            if (atom->has_anisou()) return false; // skip moving water if existing ANISOU
            atom_list.push_back(atom);
            atom = water_residue.residue->GetNextAtom();
       }
       if (atom_list.empty()) return true;

       water_residue.moving_flag = true;
       water_residue.ref_atom = contact.b_atm;
       water_residue.pdb_chnid = contact.b_atm->pdb_chnid();
       water_residue.dist = contact.dist;
       water_residue.sym = contact.sym;
       water_residue.lx = -contact.lx;
       water_residue.ly = -contact.ly;
       water_residue.lz = -contact.lz;

       _MOVING_ATOM change_atom;
       change_atom.symmetry.clear();
       change_atom.symmetry_as_xyz.clear();
       change_atom.sym = water_residue.sym;
       change_atom.lx = water_residue.lx;
       change_atom.ly = water_residue.ly;
       change_atom.lz = water_residue.lz;

       COORD coord = atom_list[0]->orig();
       _cell->symmetry_operation(coord, contact.sym, -contact.lx, -contact.ly, -contact.lz);
       change_atom.atom = atom_list[0];
       change_atom.old_coord = atom_list[0]->orig();
       _moving_waters->push_back(change_atom);
       atom_list[0]->set_orig(coord);

       for (unsigned int i = 1; i < atom_list.size(); ++i) {
            coord = atom_list[i]->orig();
            _cell->symmetry_operation(coord, contact.sym, -contact.lx, -contact.ly, -contact.lz);
            double dist1 = cal_distance(atom_list[0]->orig(), coord);
            double dist2 = cal_distance(atom_list[0]->orig(), atom_list[i]->orig());
            if (dist2 < dist1) continue;

            change_atom.atom = atom_list[i];
            change_atom.old_coord = atom_list[i]->orig();
            _moving_waters->push_back(change_atom);
            atom_list[i]->set_orig(coord);
       }

       return true;
}
