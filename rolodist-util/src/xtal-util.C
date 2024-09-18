/*
FILE:     xtal-util.C
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
#include <math.h>

#include "CompositeIndex.h"
#include "utillib.h"
#include "xtal.h"

static void set_block(xtal& mycrys, RCSB::Molecule* mol /* , const std::set<std::string>& chain_types, const std::set<std::string>& res_names,
                      const bool& remove_hydrogen_flag */);
static void set_block(xtal& mycrys, const std::vector<RCSB::Residue*>& residues);
static int contact_filter_check(const CONTACT& ptr, const double& hi);

bool set_mycrys(xtal& mycrys, const CrySymmetry& cell, RCSB::Molecule* mol, const std::set<std::string>& chain_types, const std::set<std::string>& res_names,
                const bool& remove_hydrogen_flag, const bool& exculde_water_flag, const bool& asym_only_flag)
{
       if (!cell.is_artifical()) mycrys.set_crystal(cell);

       bool found_atom = false;
       std::vector<RCSB::Residue*> residues;
       RCSB::Chain* chain = mol->GetFirstChain();
       while (chain) {
            if (exculde_water_flag && chain->chain_type() == "HETAS") {
                 RCSB::Residue* res = chain->GetFirstResidue();
                 if (res->ResName() == "HOH" || res->ResName() == "DOD") {
                      chain = mol->GetNextChain();
                      continue;
                 }
            }
            if (!chain_types.empty() && chain_types.find(chain->chain_type()) == chain_types.end()) {
                 chain = mol->GetNextChain();
                 continue;
            }

            chain->GetFirstResidueList(residues);
            while (!residues.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator pos = residues.begin(); pos != residues.end(); ++pos) {
                      if (!res_names.empty() && res_names.find((*pos)->ResName()) == res_names.end()) continue;

                      RCSB::Atom* atom = (*pos)->GetFirstAtom();
                      while (atom) {
                           if (((*pos)->ResName() != "D8U") && ((remove_hydrogen_flag && atom->is_hydrogen()) || atom->is_9999_flag())) {
                                atom = (*pos)->GetNextAtom();
                                continue;
                           }
                           found_atom = true;
                           int chn_type = NONPOLY_TYPE;
                           if ((chain->chain_type() == "ATOMN") || (chain->chain_type() == "ATOMP"))
                                chn_type = POLYMER_TYPE;
                           else if (chain->chain_type() == "ATOMS")
                                chn_type = BRANCH_TYPE;
                           else if ((chain->chain_type() == "HETAS") && (((*pos)->ResName() == "HOH") || ((*pos)->ResName()  == "DOD")))
                                chn_type = WATER_TYPE;
                           mycrys.add_atom(atom, chain->index(), (*pos)->position(), (*pos)->index(), chn_type);
                           atom = (*pos)->GetNextAtom();
                      }
                 }
                 chain->GetNextResidueList(residues);
            }
            chain = mol->GetNextChain();
       }
       if (!found_atom) return found_atom;
       
       if (cell.is_artifical())
            set_block(mycrys, mol /* , chain_types, res_names, remove_hydrogen_flag */ );
       if (asym_only_flag)
            mycrys.genAtoms();
       else mycrys.genSymmAtoms();
       mycrys.moveAllAtomsToCenter();
       mycrys.placeAtomsInGrid();

       return found_atom;
}

void set_mycrys(xtal& mycrys, const CrySymmetry& cell, const std::vector<RCSB::Residue*>& residues, const bool& remove_hydrogen_flag, const bool& asym_only_flag)
{
       if (!cell.is_artifical()) mycrys.set_crystal(cell);

       for (std::vector<RCSB::Residue*>::const_iterator pos = residues.begin(); pos != residues.end(); ++pos) {
            RCSB::Atom* atom = (*pos)->GetFirstAtom();
            while (atom) {
                 if (((*pos)->ResName() != "D8U") && ((remove_hydrogen_flag && atom->is_hydrogen()) || atom->is_9999_flag())) {
                      atom = (*pos)->GetNextAtom();
                      continue;
                 }
                 int chn_type = NONPOLY_TYPE;
                 if (((*pos)->chain_type() == "ATOMN") || ((*pos)->chain_type() == "ATOMP"))
                      chn_type = POLYMER_TYPE;
                 else if ((*pos)->chain_type() == "ATOMS")
                      chn_type = BRANCH_TYPE;
                 else if (((*pos)->chain_type() == "HETAS") && (((*pos)->ResName() == "HOH") || ((*pos)->ResName()  == "DOD")))
                      chn_type = WATER_TYPE;
                 mycrys.add_atom(atom, (*pos)->chain_index(), (*pos)->position(), (*pos)->index(), chn_type);
                 atom = (*pos)->GetNextAtom();
            }
       }
       
       if (cell.is_artifical())
            set_block(mycrys, residues);
       if (asym_only_flag)
            mycrys.genAtoms();
       else mycrys.genSymmAtoms();
       mycrys.moveAllAtomsToCenter();
       mycrys.placeAtomsInGrid();
}

void get_contact(xtal& mycrys, std::vector<CONTACT>& contact_list, const double& lo, const double& hi, const std::string& inter, const std::string& fromtype,
                 const std::string& totype, const int& DeleteMark, const int& filter_flag, const bool& remove_all_flag)
{
       contact_list.clear();

       int remove_flag = SAME_RESIDUE_FLAG;
       if (remove_all_flag) remove_flag = SAME_RESIDUE_FLAG | SAME_CHAIN_FLAG;
       std::list<CONTACT> pl;
       mycrys.getContactList(pl, fromtype, totype, lo, hi, remove_flag);
       if (pl.empty()) return;

       int asym = 0;
       if (inter == "asym") 
            asym = 1;
       else if (inter == "sym")
            asym = 2;

       contact_list.reserve(pl.size());

       std::map<std::string, int> _CIndex;
       _CIndex.clear();

       for (std::list<CONTACT>::iterator ptr = pl.begin(); ptr != pl.end(); ++ptr) {
            int i = 1;
            if (asym == 1) {
                 if (!ptr->sym && !ptr->lx && !ptr->ly && !ptr->lz) i = 0;
            } else if (asym == 2) {
                 if (ptr->sym || ptr->lx || ptr->ly || ptr->lz) i = 0;
            } else i = 0;

            if (i > 0) continue;

            if (filter_flag) i = contact_filter_check(*ptr, hi);

            if (i > 0) continue;

            bool add = true;
            if (DeleteMark) {
                 std::string cs;
                 if (atoi(ptr->a_atm->atnum().c_str()) < atoi(ptr->b_atm->atnum().c_str())) {
                      cs = ptr->a_atm->atnum() + "_" + ptr->b_atm->atnum();
                 } else {
                      cs = ptr->b_atm->atnum() + "_" + ptr->a_atm->atnum();
                 }
                 std::map<std::string, int>::iterator pos = _CIndex.find(cs);
                 if (pos != _CIndex.end()) {
                      if ((!ptr->sym && !ptr->lx && !ptr->ly && !ptr->lz && (contact_list[pos->second].sym > 1 || contact_list[pos->second].lx ||
                          contact_list[pos->second].ly || contact_list[pos->second].lz)) || (contact_list[pos->second].dist - ptr->dist) > 0.001) {
                           contact_list[pos->second] = *ptr;
                           contact_list[pos->second].type = i;
                           contact_list[pos->second].sym = ptr->sym + 1;
                      }
                      add = false;
                 } else _CIndex.insert(std::make_pair(cs, contact_list.size()));
            }
            if (add) {
                 ptr->sym++;
                 ptr->type = i;
                 contact_list.push_back(*ptr);
            }
       }
}

void get_sym_contact_with_special_position_residues(xtal& mycrys, std::vector<CONTACT>& contact_list, std::map<std::string, std::map<std::string,
                 std::set<std::string> > >& residues, const double& lo, const double& hi, const std::string& fromtype, const std::string& totype)
{
       contact_list.clear();

       // first key: residue ID - pdbchnid_resname_resnum_inscode
       // second key: atom ID - atmnam_altloc
       // values: symmetry operators
       residues.clear();

       int remove_flag = SAME_RESIDUE_FLAG | SAME_CHAIN_FLAG;
       std::list<CONTACT> pl;
       mycrys.getContactList(pl, fromtype, totype, lo, hi, remove_flag);
       if (pl.empty()) return;

       contact_list.reserve(pl.size());

       std::map<std::string, int> _CIndex;
       _CIndex.clear();

       std::map<std::string, std::set<std::string> > tmp_map;
       std::set<std::string> t_set;
       std::string cs;

       for (std::list<CONTACT>::iterator ptr = pl.begin(); ptr != pl.end(); ++ptr) {
            int i = 1;
            if (ptr->sym || ptr->lx || ptr->ly || ptr->lz) i = 0;
            if (i > 0) continue;

            if (!ptr->a_atm->is_hydrogen() && (ptr->dist < SPECIAL_WATER_DIST_CUTOFF) && (ptr->a_atm == ptr->b_atm)) {
                 cs = CompositeIndex::getIndex(ptr->a_atm->pdb_chnid(), ptr->a_atm->restype(), ptr->a_atm->pdb_resnum(), ptr->a_atm->ins_code());
                 std::string atomID = ptr->a_atm->pdb_atmnam() + "_" + ptr->a_atm->alt_loc();
                 std::string sym = String::IntToString(ptr->sym) + "_" + String::IntToString(ptr->lx + 5) + String::IntToString(ptr->ly + 5)
                                 + String::IntToString(ptr->lz + 5);
                 std::map<std::string, std::map<std::string, std::set<std::string> > >::iterator mpos = residues.find(cs);
                 if (mpos != residues.end()) {
                      std::map<std::string, std::set<std::string> >::iterator mmpos = mpos->second.find(atomID);
                      if (mmpos != mpos->second.end()) mmpos->second.insert(sym);
                      else {
                           t_set.clear();
                           t_set.insert(sym);
                           mpos->second.insert(std::make_pair(atomID, t_set));
                      }
                 } else {
                      t_set.clear();
                      t_set.insert(sym);
                      tmp_map.clear();
                      tmp_map.insert(std::make_pair(atomID, t_set));
                      residues.insert(std::make_pair(cs, tmp_map));
                 }
            }

            i = contact_filter_check(*ptr, hi);
            if (i > 0) continue;

            bool add = true;
            if (atoi(ptr->a_atm->atnum().c_str()) < atoi(ptr->b_atm->atnum().c_str())) {
                 cs = ptr->a_atm->atnum() + "_" + ptr->b_atm->atnum();
            } else {
                 cs = ptr->b_atm->atnum() + "_" + ptr->a_atm->atnum();
            }
            std::map<std::string, int>::iterator pos = _CIndex.find(cs);
            if (pos != _CIndex.end()) {
                 if ((!ptr->sym && !ptr->lx && !ptr->ly && !ptr->lz && (contact_list[pos->second].sym > 1 || contact_list[pos->second].lx ||
                      contact_list[pos->second].ly || contact_list[pos->second].lz)) || (contact_list[pos->second].dist - ptr->dist) > 0.001) {
                      contact_list[pos->second] = *ptr;
                      contact_list[pos->second].type = i;
                      contact_list[pos->second].sym = ptr->sym + 1;
                 }
                 add = false;
            } else _CIndex.insert(std::make_pair(cs, contact_list.size()));

            if (add) {
                 ptr->sym++;
                 ptr->type = i;
                 contact_list.push_back(*ptr);
            }
       }
}

bool find_closest_contact(cryatom& newat, double *cr, const double& lo, const double &hi, CONTACT& contact)
{
       contact.a_atm = NULL;
       contact.b_atm = NULL;
       contact.a_res_index = -1;
       contact.b_res_index = -1;
       contact.sym = contact.lx = contact.ly = contact.lz = 0;
       contact.type = contact.mol_index = 0;
       contact.chn_type = 0;
       contact.dist = 0;

       int cells[3];
       for (int i = 0; i < 3; ++i) {
            cells[i] = int (floor (1 + (hi/(cr[i]))));
       }

       std::multimap<double, CONTACT> asym_contacts;
       asym_contacts.clear();

       double low = lo - 0.001;
       if (low < 0.5) low = 0.5;

       std::list<CONTACT> found_list;
       found_list.clear();

       newat.neighborAtoms(found_list, low, hi, cells, 0);
       if (found_list.empty()) return false;
       for (std::list<CONTACT>::iterator lpos = found_list.begin(); lpos != found_list.end(); ++lpos) {
            lpos->sym = newat.symnum();
            if (lpos->sym || lpos->lx || lpos->ly || lpos->lz) continue;
            asym_contacts.insert(std::make_pair(lpos->dist, *lpos));
       }

       if (asym_contacts.empty()) return false;

       contact = asym_contacts.begin()->second;
       return true;
}

void print_contact(const std::vector<CONTACT>& contact)
{
       printf("contact=%d\n", (int) contact.size());
       for (std::vector<CONTACT>::const_iterator pos = contact.begin(); pos != contact.end(); ++pos) {
            std::string a = pos->a_atm->alt_loc();
            if (!pos->a_atm->alt_loc().empty()) a = "(" + pos->a_atm->alt_loc() + ")";
            std::string b = pos->b_atm->alt_loc();
            if (!pos->b_atm->alt_loc().empty()) b = "(" + pos->b_atm->alt_loc() + ")";
            printf("%5s%s %3s %5s %s --- %5s%s %3s %5s %s --- %d_%d%d%d %.3f\n",
                   pos->a_atm->atmtype().c_str(), a.c_str(), pos->a_atm->pdb_resnam().c_str(), pos->a_atm->pdb_resnum().c_str(), pos->a_atm->pdb_chnid().c_str(),
                   pos->b_atm->atmtype().c_str(), b.c_str(), pos->b_atm->pdb_resnam().c_str(), pos->b_atm->pdb_resnum().c_str(), pos->b_atm->pdb_chnid().c_str(),
                   pos->sym, pos->lx, pos->ly, pos->lz, pos->dist);
       }
}

void print_contact(const std::list<CONTACT>& contact)
{
       printf("contact=%d\n", (int) contact.size());
       for (std::list<CONTACT>::const_iterator pos = contact.begin(); pos != contact.end(); ++pos) {
            std::string a = pos->a_atm->alt_loc();
            if (!pos->a_atm->alt_loc().empty()) a = "(" + pos->a_atm->alt_loc() + ")";
            std::string b = pos->b_atm->alt_loc();
            if (!pos->b_atm->alt_loc().empty()) b = "(" + pos->b_atm->alt_loc() + ")";
            printf("%5s%s %3s %5s %s --- %5s%s %3s %5s %s --- %d_%d%d%d %.3f\n",
                   pos->a_atm->atmtype().c_str(), a.c_str(), pos->a_atm->pdb_resnam().c_str(), pos->a_atm->pdb_resnum().c_str(), pos->a_atm->pdb_chnid().c_str(),
                   pos->b_atm->atmtype().c_str(), b.c_str(), pos->b_atm->pdb_resnam().c_str(), pos->b_atm->pdb_resnum().c_str(), pos->b_atm->pdb_chnid().c_str(),
                   pos->sym, pos->lx, pos->ly, pos->lz, pos->dist);
       }
}

static void set_block(xtal& mycrys, RCSB::Molecule* mol /* , const std::set<std::string>& chain_types,
                      const std::set<std::string>& res_names, 
                      const bool& remove_hydrogen_flag */ )
{
       double Xmin = 1000.0;
       double Ymin = 1000.0;
       double Zmin = 1000.0;
       double Xmax = -1000.0;
       double Ymax = -1000.0;
       double Zmax = -1000.0;

       std::vector<RCSB::Residue*> residues;
       RCSB::Chain* chain = mol->GetFirstChain();
       while (chain) {
/*
            if (!chain_types.empty() && chain_types.find(chain->chain_type()) !=
                 chain_types.end()) {
                 chain = mol->GetNextChain();
                 continue;
            }
*/
            chain->GetFirstResidueList(residues);
            while (!residues.empty()) {
                 for (std::vector<RCSB::Residue*>::const_iterator pos = residues.begin(); pos != residues.end(); ++pos) {
/*
                      if (!res_names.empty() && res_names.find((*pos)->ResName()) !=
                           res_names.end()) continue;
*/
                      RCSB::Atom* atom = (*pos)->GetFirstAtom();
                      while (atom) {
                           if (atom->is_9999_flag()) {
                           // if (remove_hydrogen_flag && atom->is_hydrogen()) {
                                atom = (*pos)->GetNextAtom();
                                continue;
                           }

                           if (atom->orig().x < Xmin) Xmin = atom->orig().x;
                           if (atom->orig().x > Xmax) Xmax = atom->orig().x;
                           if (atom->orig().y < Ymin) Ymin = atom->orig().y;
                           if (atom->orig().y > Ymax) Ymax = atom->orig().y;
                           if (atom->orig().z < Zmin) Zmin = atom->orig().z;
                           if (atom->orig().z > Zmax) Zmax = atom->orig().z;
                           atom = (*pos)->GetNextAtom();
                      }
                 }
                 chain->GetNextResidueList(residues);
            }
            chain = mol->GetNextChain();
       }

       mycrys.set_block(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
}

static void set_block(xtal& mycrys, const std::vector<RCSB::Residue*>& residues)
{
       double Xmin = 1000.0;
       double Ymin = 1000.0;
       double Zmin = 1000.0;
       double Xmax = -1000.0;
       double Ymax = -1000.0;
       double Zmax = -1000.0;

       for (std::vector<RCSB::Residue*>::const_iterator pos = residues.begin(); pos != residues.end(); ++pos) {
            RCSB::Atom* atom = (*pos)->GetFirstAtom();
            while (atom) {
                 if (atom->is_9999_flag()) {
                      atom = (*pos)->GetNextAtom();
                      continue;
                 }
                 if (atom->orig().x < Xmin) Xmin = atom->orig().x;
                 if (atom->orig().x > Xmax) Xmax = atom->orig().x;
                 if (atom->orig().y < Ymin) Ymin = atom->orig().y;
                 if (atom->orig().y > Ymax) Ymax = atom->orig().y;
                 if (atom->orig().z < Zmin) Zmin = atom->orig().z;
                 if (atom->orig().z > Zmax) Zmax = atom->orig().z;
                 atom = (*pos)->GetNextAtom();
            }
       }

       mycrys.set_block(Xmin, Ymin, Zmin, Xmax, Ymax, Zmax);
}

static int contact_filter_check(const CONTACT& ptr, const double& hi)
{
       int i = 0;
/*
       if (ptr.dist <= 2.2 && !ptr.sym && !ptr.lx && !ptr.ly && !ptr.lz &&
           ptr.a_atm->atmtype() == "SG" && ptr.b_atm->atmtype() == "SG" &&
          (ptr.b_atm->restype() == "CYS" || ptr.b_atm->restype() == "DCY") &&
          (ptr.a_atm->restype() == "CYS" || ptr.a_atm->restype() == "DCY"))
            i = 1;
*/
       if ((ptr.a_atm->type() == "ATOMP" || ptr.a_atm->type() == "ATOMN") && ptr.a_atm->type() == ptr.b_atm->type() &&
            ptr.a_atm->chnid() == ptr.b_atm->chnid() && (atoi(ptr.a_atm->resnum().c_str()) + 1 == atoi(ptr.b_atm->resnum().c_str()) ||
            atoi(ptr.a_atm->resnum().c_str()) == atoi(ptr.b_atm->resnum().c_str()) + 1) && ptr.dist >= 1.8)
            i = 1;
       else if (ptr.a_atm == ptr.b_atm && ptr.dist < 0.3)
            i = -1;
       else if (atof(ptr.a_atm->occ().c_str()) < 1.0 && atof(ptr.b_atm->occ().c_str()) < 1.0) {
            if (ptr.sym && ptr.a_atm->restype() == ptr.b_atm->restype() &&
                ptr.a_atm->chnid() == ptr.b_atm->chnid() &&
                ptr.a_atm->resnum() == ptr.b_atm->resnum()) { 
                 i = -1;
            }
       } else if (ptr.a_atm->atom_type() == "D" && ptr.b_atm->atom_type() == "D") {
            if (ptr.dist >= 1.415) i = 1;
       } else if ((ptr.a_atm->atom_type() == "H" || ptr.a_atm->atom_type() == "D") && (ptr.b_atm->atom_type() == "H" || ptr.b_atm->atom_type() == "D")) {
            if (ptr.dist >= 1.35) i = 1;
       } else if ((ptr.a_atm->atom_type() == "H" || ptr.b_atm->atom_type() == "H" || ptr.a_atm->atom_type() == "D" || ptr.b_atm->atom_type() == "D") &&
                 ptr.dist >= 1.6) {
            i = 1;
       } else if ((ptr.a_atm->atom_type() == "BA" || ptr.a_atm->atom_type() == "MG" || ptr.a_atm->atom_type() == "NA" || ptr.a_atm->atom_type() == "NI" ||
                   ptr.a_atm->atom_type() == "CU" || ptr.a_atm->atom_type() == "CO" || ptr.a_atm->atom_type() == "CA" || ptr.a_atm->atom_type() == "HG" ||
                   ptr.a_atm->atom_type() == "MN" || ptr.a_atm->atom_type() == "PT" || ptr.a_atm->atom_type() == "ZN" || ptr.b_atm->atom_type() == "BA" ||
                   ptr.b_atm->atom_type() == "MG" || ptr.b_atm->atom_type() == "NA" || ptr.b_atm->atom_type() == "NI" || ptr.b_atm->atom_type() == "CU" ||
                   ptr.b_atm->atom_type() == "CO" || ptr.b_atm->atom_type() == "CA" || ptr.b_atm->atom_type() == "HG" || ptr.b_atm->atom_type() == "MN" ||
                   ptr.b_atm->atom_type() == "PT" || ptr.b_atm->atom_type() == "ZN") && ptr.dist >= 1.7)
            i = 1;
       else if (hi <= 2.2 && ptr.dist >= 2.195)
            i = 1;

       return i;
}
