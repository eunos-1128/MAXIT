/*
FILE:     CheckDistantWater.C
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

#include "utillib.h"
#include "ValidateObj.h"
#include "xtal.h"

static void _checkwaters(xtal& mycrys, const int& mol_id, const std::list<RCSB::Residue*>& water_lists,
                       std::list<std::pair<RCSB::Atom*, std::vector<std::string> > >& out_range_waters);
static void _check_residues(xtal& mycrys, const std::string& Mol_ID, const std::list<RCSB::Residue*>& residue_lists,
                            std::list<std::pair<std::string, RCSB::Residue*> >& distant_residue_lists);

void ValidateObj::CheckDistantWaters()
{
       _cal_out_range_water_flag = true;
       _out_range_waters.clear();
       std::vector<RCSB::Chain*> polymer_chains, carbohydrate_chains, nonpolymer_chains;
       std::list<RCSB::Residue*> water_lists;
       for (std::vector<RCSB::Molecule*>::iterator pos = _molecules.begin(); pos != _molecules.end(); ++pos) {
            (*pos)->GetChainsAndWaters(polymer_chains, carbohydrate_chains, nonpolymer_chains, water_lists);
            _check_distant_waters(*pos, water_lists);
       }
}

void ValidateObj::_check_distant_waters(RCSB::Molecule* mol, const std::list<RCSB::Residue*>& water_lists)
{
       
       if (water_lists.empty()) return;

       std::set<std::string> chain_types;
       chain_types.clear();
       std::set<std::string> res_names;
       res_names.clear();

       CrySymmetry cell;
       cell.clear();
       cell.setArtificalFlag(2);

       xtal mycrys;
       if (!set_mycrys(mycrys, cell, mol, chain_types, res_names, true, true, true)) return;

       _checkwaters(mycrys, mol->Mol_ID(), water_lists, _out_range_waters);
}

void ValidateObj::_check_distant_nonpolymer_residues(RCSB::Molecule* mol, const std::string& Mol_ID, const std::list<RCSB::Residue*>& residue_lists,
                                                     std::list<std::pair<std::string, RCSB::Residue*> >& distant_residue_lists)
{
       if (residue_lists.empty()) return;

       std::set<std::string> chain_types;
       chain_types.clear();
       chain_types.insert("ATOMN");
       chain_types.insert("ATOMP");
       chain_types.insert("ATOMS");
/*
       chain_types.insert("HETAD");
       chain_types.insert("HETAI");
       chain_types.insert("HETIC");
       chain_types.insert("HETAC");
       chain_types.insert("HETAIN");
       chain_types.insert("HETAS");
*/
       std::set<std::string> res_names;
       res_names.clear();

       CrySymmetry cell;
       cell.clear();
       cell.setArtificalFlag(2);

       xtal mycrys;
       if (!set_mycrys(mycrys, cell, mol, chain_types, res_names, true, true, true)) return;

       _check_residues(mycrys, Mol_ID, residue_lists, distant_residue_lists);
}

static void _checkwaters(xtal& mycrys, const int& mol_id, const std::list<RCSB::Residue*>& water_lists,
                       std::list<std::pair<RCSB::Atom*, std::vector<std::string> > >& out_range_waters)
{
       if (water_lists.empty()) return;

       mycrys.markAtoms("polar_atom");

       grid *gr = mycrys.crysGrid();
       int elts_x = gr->eltsSide_x();
       int elts_y = gr->eltsSide_y();
       int elts_z = gr->eltsSide_z();

       double cr[3];
       cr[0] = mycrys.aval() / elts_x;
       cr[1] = mycrys.bval() / elts_y;
       cr[2] = mycrys.cval() / elts_z;

       cryatom newat;
       newat.setcrys(&mycrys);

       CONTACT contact;
       std::vector<std::string> data;

       for (std::list<RCSB::Residue*>::const_iterator pos = water_lists.begin(); pos != water_lists.end(); ++pos) {
            RCSB::Atom* atom = (*pos)->find_atom("O");
            if (atom == NULL) continue;

            newat.clear();
            newat.setatom(atom, 100000, (*pos)->position(), (*pos)->index(), WATER_TYPE);
            newat.moveToCenterCell();
            newat.setGridCoords();

            double lo = 0;
            double hi = WATER_CUTOFF;
            int count = 0;
            while (count <= 30) {
                 if (find_closest_contact(newat, cr, lo, hi, contact)) {
                      if (contact.dist > WATER_CUTOFF) {
                           data.clear();
                           data.push_back(String::IntToString(mol_id));
                           data.push_back(FloatToString(contact.dist, 0, 2));
                           if (contact.chn_type == POLYMER_TYPE)
                                data.push_back("polymer");
                           else data.push_back("non-polymer");
                           out_range_waters.push_back(std::make_pair(atom, data));
                      }
                      break;
                 }
                 lo += WATER_CUTOFF;
                 hi += WATER_CUTOFF;
                 count++;
            }
       }

       mycrys.unmarkAtoms();
}

static void _check_residues(xtal& mycrys, const std::string& Mol_ID, const std::list<RCSB::Residue*>& residue_lists,
                            std::list<std::pair<std::string, RCSB::Residue*> >& distant_residue_lists)
{
       mycrys.markAtoms("all");

       grid *gr = mycrys.crysGrid();
       int elts_x = gr->eltsSide_x();
       int elts_y = gr->eltsSide_y();
       int elts_z = gr->eltsSide_z();

       double cr[3];
       cr[0] = mycrys.aval() / elts_x;
       cr[1] = mycrys.bval() / elts_y;
       cr[2] = mycrys.cval() / elts_z;

       cryatom newat;
       newat.setcrys(&mycrys);

       CONTACT contact;
       std::vector<std::string> data;

       for (std::list<RCSB::Residue*>::const_iterator pos = residue_lists.begin(); pos != residue_lists.end(); ++pos) {
            bool found_contact = false;
            RCSB::Atom* atom = (*pos)->GetFirstAtom();
            while (atom) {
                 newat.clear();
                 newat.setatom(atom, 100000, (*pos)->position(), (*pos)->index(), NONPOLY_TYPE);
                 newat.moveToCenterCell();
                 newat.setGridCoords();

                 if (find_closest_contact(newat, cr, 0, 15.0, contact)) {
                      if (contact.dist <= 15.0) {
                           found_contact = true;
                           break;
                      }
                 }
                 atom = (*pos)->GetNextAtom();
            }
            if (!found_contact) distant_residue_lists.push_back(std::make_pair(Mol_ID, *pos));
       }

       mycrys.unmarkAtoms();
}
