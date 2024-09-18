/*
FILE:     xtal.C
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

#include "Element.h"
#include "utillib.h"
#include "xtal.h"

#define BLOCK_SIZE 6.0

xtal::xtal()
{
       _init();
}

xtal::~xtal()
{
       clear();
}

void xtal::_init()
{
       is_artifical = true;
       numatoms = 0;
       numatomsAsym = 0;
       a = b = c = 0.0;
       alpha = beta = gamma = 90.0;
       cosalpha = cosbeta = cosgamma = 1.0; 
       atlist = NULL;
       xtalgrid = NULL;
       atl.clear();
       oplist.clear();
       _grid_size = BLOCK_SIZE;
       _Xmin = -1000.0;
       _Ymin = -1000.0;
       _Zmin = -1000.0;
       _Xmax =  1000.0;
       _Ymax =  1000.0;
       _Zmax =  1000.0;
       _polar_atoms.clear();
       _non_branch_polar_atoms.clear();
       _polymer_atoms.clear();
       _branch_atoms.clear();
       _nonpolymer_atoms.clear();
       _water_atoms.clear();
}

void xtal::clear()
{
       if (atlist) delete [] atlist;
       if (xtalgrid) delete xtalgrid;
       _init();
}

void xtal::set_crystal(const CrySymmetry& cell)
{
       if (cell.is_artifical()) return;

       is_artifical = false;

       a = atof(cell.a().c_str());
       b = atof(cell.b().c_str());
       c = atof(cell.c().c_str());
       alpha = atof(cell.alpha().c_str());
       beta = atof(cell.beta().c_str());
       gamma = atof(cell.gamma().c_str());

       cosalpha = float(cos(double(degtorad * alpha)));
       cosbeta = float(cos(double(degtorad * beta)));
       cosgamma = float(cos(double(degtorad * gamma)));

       for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) xform[i][j] = cell.o_to_f[i][j];
            xform[i][3] = cell.o_to_f_v[i];
       }

       symop newsym;
       float oper[3][4];
       const std::vector<NDBSYMMETRY>& symops = cell.symops();
       for (unsigned int i = 0; i < symops.size(); ++i) {
            newsym.setopid(i);
            for (int j = 0; j < 3; ++j) {
                 for (int k = 0; k < 3; ++k) {
                      oper[j][k] = symops[i].rot[j][k];
                 }
                 oper[j][3] = symops[i].trans[j];
            }
            newsym.setop(oper);
            oplist.push_back(newsym);
       }
}

void xtal::set_grid_size(const double& grid_size)
{
       _grid_size = grid_size;
}

void xtal::set_block(const double& xmin, const double& ymin, const double& zmin,
                     const double& xmax, const double& ymax, const double& zmax)
{
       is_artifical = true;
       _Xmin = xmin - _grid_size * 0.5;
       _Ymin = ymin - _grid_size * 0.5;
       _Zmin = zmin - _grid_size * 0.5;
       _Xmax = xmax + _grid_size * 0.5;
       _Ymax = ymax + _grid_size * 0.5;
       _Zmax = zmax + _grid_size * 0.5;
}

void xtal::add_atom(RCSB::Atom* atom, const int& chn_index, const int& position, const int& res_index, const int& chn_type, const bool& ignoreAltLocFlag)
{
       cryatom newat;
       newat.setcrys(this);
       newat.setatom(atom, chn_index, position, res_index, chn_type, ignoreAltLocFlag);
       atl.push_back(newat);
}

void xtal::genAtoms()
{
       atlist = new cryatom[atl.size()];
       int ats = 0;
       std::list<cryatom>::iterator apos;
       for (apos = atl.begin(); apos != atl.end(); ++apos) {
            atlist[ats] = *apos;
            if (atlist[ats].atomtype() == "N" || atlist[ats].atomtype() == "O" || Element::findMetalFlag(atlist[ats].atomtype()) == 2) {
                 _polar_atoms.insert(ats);
                 if (atlist[ats].chn_type() != BRANCH_TYPE) _non_branch_polar_atoms.insert(ats);
            }
            if (atlist[ats].chn_type() == POLYMER_TYPE)
                 _polymer_atoms.insert(ats);
            else if (atlist[ats].chn_type() == BRANCH_TYPE)
                 _branch_atoms.insert(ats);
            else if (atlist[ats].chn_type() == NONPOLY_TYPE)
                 _nonpolymer_atoms.insert(ats);
            else if (atlist[ats].chn_type() == WATER_TYPE)
                 _water_atoms.insert(ats);
            ats++;
       }
       numatomsAsym = ats;
       numatoms = ats;
       atl.clear();
}

// Use the list of symmetry operations on the list of atoms in the asymmetric
// unit to generate the array of atoms in the entire unit cell.  If the
// atoms are known to be in cartesian, then first convert all of them to
// fractional.  Delete the initial list data structure, to save space.  
// Each copy of the asymmetric unit is contiguous.
void xtal::genSymmAtoms()
{
       int numsym = oplist.size();
       if (numsym == 0) numsym = 1;
       numatoms = atl.size() * numsym;
       atlist = new cryatom[numatoms + 1];
       int ats = 0;
       std::list<cryatom>::iterator apos;
       for (apos = atl.begin(); apos != atl.end(); ++apos) {
            atlist[ats] = *apos;
            if (atlist[ats].atomtype() == "N" || atlist[ats].atomtype() == "O" || Element::findMetalFlag(atlist[ats].atomtype()) == 2) {
                 _polar_atoms.insert(ats);
                 if (atlist[ats].chn_type() != BRANCH_TYPE) _non_branch_polar_atoms.insert(ats);
            }
            if (atlist[ats].chn_type() == POLYMER_TYPE)
                 _polymer_atoms.insert(ats);
            else if (atlist[ats].chn_type() == BRANCH_TYPE)
                 _branch_atoms.insert(ats);
            else if (atlist[ats].chn_type() == NONPOLY_TYPE)
                 _nonpolymer_atoms.insert(ats);
            else if (atlist[ats].chn_type() == WATER_TYPE)
                 _water_atoms.insert(ats);
            ats++;
       }
       numatomsAsym = ats;

       if (oplist.empty()) {
            numatoms = ats;
            atl.clear();
            return;
       }

       std::list<symop>::iterator spos;
       spos = oplist.begin(); ++spos; // assumes first symmetry is identity
       for ( ; spos != oplist.end(); ++spos) {
            for (apos = atl.begin(); apos != atl.end(); ++apos) {
                 atlist[ats] = apos->symrel(*spos);
                 if (atlist[ats].atomtype() == "N" || atlist[ats].atomtype() == "O" || Element::findMetalFlag(atlist[ats].atomtype()) == 2) {
                      _polar_atoms.insert(ats);
                      if (atlist[ats].chn_type() != BRANCH_TYPE) _non_branch_polar_atoms.insert(ats);
                 }
                 if (atlist[ats].chn_type() == POLYMER_TYPE)
                      _polymer_atoms.insert(ats);
                 else if (atlist[ats].chn_type() == BRANCH_TYPE)
                      _branch_atoms.insert(ats);
                 else if (atlist[ats].chn_type() == NONPOLY_TYPE)
                      _nonpolymer_atoms.insert(ats);
                 else if (atlist[ats].chn_type() == WATER_TYPE)
                      _water_atoms.insert(ats);
                 ats++;
            }
       }
       numatoms = ats;
       atl.clear();
}

//////////////////////////////////////////////////////////////////////
//		     Data Member Access Methods
//////////////////////////////////////////////////////////////////////

// Return a pointer to the grid partitioning the crystal
grid* xtal::crysGrid() { return xtalgrid; }

void xtal::getSymopList(std::vector<symop>& voplist)
{
       voplist.clear();
       if (oplist.empty()) return;

       std::list<symop>::iterator spos = oplist.begin();
       ++spos; // assumes first symmetry is identity
       for ( ; spos != oplist.end(); ++spos) voplist.push_back(*spos);
}

// Return the cell dimensions of the crystal
const double& xtal::aval() const { return a; }
const double& xtal::bval() const { return b; }
const double& xtal::cval() const { return c; }
const double& xtal::grid_size() const { return _grid_size; }
const double& xtal::Xmin() const { return _Xmin; }
const double& xtal::Ymin() const { return _Ymin; }
const double& xtal::Zmin() const { return _Zmin; }

//////////////////////////////////////////////////////////////////////
//		   Data member manipulation methods
//////////////////////////////////////////////////////////////////////

// Translate all atoms in the crystal into the center unit cell.  Note
// this should only be called after all of the symmetry related atoms
// have been generated.  Otherwise the program will crash (atom array
// not generated yet) and even if it didn't, would not work correctly,
// since after generating the symmetry related atoms, some of them might
// not be in the center unit cell.
void xtal::moveAllAtomsToCenter()
{
       if (is_artifical) return;
       for (int i = 0; i < numatoms; ++i) {
            atlist[i].moveToCenterCell();
       }
}

// Hash each atom into the appropriate grid element.  Again, make sure
// this is done *after* the symmetry related atoms are generated.
void xtal::placeAtomsInGrid()
{
       if (is_artifical) {
            a = _Xmax - _Xmin;
            b = _Ymax - _Ymin;
            c = _Zmax - _Zmin; 
       }

       int eltsx = (int) (a / _grid_size) + 1;
       int eltsy = (int) (b / _grid_size) + 1;
       int eltsz = (int) (c / _grid_size) + 1;
       if (xtalgrid) delete xtalgrid;
       xtalgrid = new grid(eltsx, eltsy, eltsz, this);
       for (int i = 0; i < numatoms; ++i) {
            atlist[i].placeAtomInGrid();
       }
}

///////////////////////////////////////////////////////////////////////
// Methods for selecting atoms with particular feature from entire list
///////////////////////////////////////////////////////////////////////

void xtal::selectAtoms(std::list<cryatom*>& atm_list, const std::string& selectType)
{
       atm_list.clear();
       if (selectType == "all")
            _selectAtoms(atm_list, 0, numatoms);
       else if (selectType == "all_asym")
            _selectAtoms(atm_list, 0, numatomsAsym);
       else if (selectType == "all_sym")
            _selectAtoms(atm_list, numatomsAsym, numatoms);
       else if (selectType == "polar_atom")
            _selectAtoms(atm_list, _polar_atoms);
       else if (selectType == "non_branch_polar_atom")
            _selectAtoms(atm_list, _non_branch_polar_atoms);
       else if (selectType == "polymer")
            _selectAtoms(atm_list, _polymer_atoms);
       else if (selectType == "branch")
            _selectAtoms(atm_list, _branch_atoms);
       else if (selectType == "nonpolymer")
            _selectAtoms(atm_list, _nonpolymer_atoms);
       else if (selectType == "water")
            _selectAtoms(atm_list, _water_atoms);
}

void xtal::markAtoms(const std::string& selectType)
{
       if (selectType == "all")
            _markAtoms(0, numatoms);
       else if (selectType == "all_asym")
            _markAtoms(0, numatomsAsym);
       else if (selectType == "all_sym")
            _markAtoms(numatomsAsym, numatoms);
       else if (selectType == "polar_atom")
            _markAtoms(_polar_atoms);
       else if (selectType == "non_branch_polar_atom")
            _markAtoms(_non_branch_polar_atoms);
       else if (selectType == "polymer")
            _markAtoms(_polymer_atoms);
      else if (selectType == "branch")
            _markAtoms(_branch_atoms);
       else if (selectType == "nonpolymer")
            _markAtoms(_nonpolymer_atoms);
       else if (selectType == "water")
            _markAtoms(_water_atoms);
}

void xtal::unmarkAtoms()
{
       xtalgrid->reset();
}

void xtal::getContactList(std::list<CONTACT>& pl, const std::list<cryatom*>& fromatoms, const double& lo, const double& hi, const int& remove_flag)
{
       pl.clear();

       if (fromatoms.empty()) return;

       int elts_x = xtalgrid->eltsSide_x();
       int elts_y = xtalgrid->eltsSide_y();
       int elts_z = xtalgrid->eltsSide_z();

       double cr[3];
       cr[0] = a / elts_x;
       cr[1] = b / elts_y;
       cr[2] = c / elts_z;

       int cells[3];
       for (int i = 0; i < 3; ++i) {
            cells[i] = int (floor (1 + (hi/(cr[i]))));
       }

       for (std::list<cryatom*>::const_iterator pos = fromatoms.begin(); pos != fromatoms.end(); ++pos) {
            (*pos)->neighborAtoms(pl, lo, hi, cells, remove_flag);
       }
}

void xtal::getContactList(std::list<CONTACT>& pl, const std::string& fromtype, const std::string& totype, const double& lo, const double& hi,
                          const int& remove_flag)
{
       pl.clear();

       std::list<cryatom*> fromatoms;
       selectAtoms(fromatoms, fromtype);
       if (fromatoms.empty()) return;

       markAtoms(totype);
       getContactList(pl, fromatoms, lo, hi, remove_flag);
       unmarkAtoms();
}

void xtal::checkAtoms()
{
       for (int i = 0; i < numatoms; ++i) atlist[i].checkAtom();
}

void xtal::_selectAtoms(std::list<cryatom*>& atm_list, const int& start, const int& end)
{
       for (int i = start; i < end; ++i) atm_list.push_back(&atlist[i]);
}

void xtal::_selectAtoms(std::list<cryatom*>& atm_list, const std::set<int>& atom_set)
{
       if (atom_set.empty()) return;
       for (std::set<int>::const_iterator pos = atom_set.begin(); pos != atom_set.end(); ++pos) {
            atm_list.push_back(&atlist[*pos]);
       }
}

void xtal::_markAtoms(const int& start, const int& end)
{
       for (int i = start; i < end; ++i) {
            atlist[i].markatom();
            xtalgrid->TheGrid()[atlist[i].gridx()][atlist[i].gridy()][atlist[i].gridz()].markcell();
       }
}

void xtal::_markAtoms(const std::set<int>& atom_set)
{
       if (atom_set.empty()) return;
       for (std::set<int>::const_iterator pos = atom_set.begin(); pos != atom_set.end(); ++pos) {
            atlist[*pos].markatom();
            xtalgrid->TheGrid()[atlist[*pos].gridx()][atlist[*pos].gridy()][atlist[*pos].gridz()].markcell();
       }
}
