/*
FILE:     cryatom.C
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
#include <string.h>

#include "cryatom.h"
#include "gridcell.h"
#include "symop.h"
#include "xtal.h"

static double sqr (const double v) { return v * v; };

// Atom constructor -- default
cryatom::cryatom ()
{
       crystal = NULL;
       clear();
}

void cryatom::clear()
{
       marked = false;
       ignoreAltLoc = false;
       atm = NULL;
       _chn_index = 0;
       _position = 0;
       _res_index = 0;
       _chn_type = 0;
       x = y = z = 0.0;
       sym = 0;
       for (int i = 0; i < 3; ++i) {
            transl[i] = 0;
            gridcoords[i] = 0;
       }
}

void cryatom::operator=(const cryatom& a)
{
       marked = a.marked;
       atm = a.atm;
       _chn_index = a._chn_index;
       _position = a._position;
       _res_index = a._res_index;
       _chn_type = a._chn_type;
       x = a.x;
       y = a.y;
       z = a.z;
       sym = a.sym;
       for (int i = 0; i < 3; ++i) {
            transl[i] = a.transl[i];
            gridcoords[i] = a.gridcoords[i];
       }
       crystal = a.crystal;
}

//////////////////////////////////////////////////////////////////////
//                  Data member access functions
//////////////////////////////////////////////////////////////////////

const bool& cryatom::isMarked() const { return marked; }
const std::string& cryatom::atomtype () const { return atm->atom_type(); }
const short& cryatom::chn_type() const { return _chn_type; }
const int& cryatom::chn_index() const  { return _chn_index; }
const double& cryatom::X() const { return x; }
const double& cryatom::Y() const { return y; }
const double& cryatom::Z() const { return z; }
const short& cryatom::symnum() const { return sym; }
const int& cryatom::gridx () const { return gridcoords[0]; }
const int& cryatom::gridy () const { return gridcoords[1]; }
const int& cryatom::gridz () const { return gridcoords[2]; }
const int& cryatom::translx () const { return transl[0]; }
const int& cryatom::transly () const { return transl[1]; }
const int& cryatom::translz () const { return transl[2]; }
xtal* cryatom::atCrys() { return crystal; }
  
//////////////////////////////////////////////////////////////////////
//            Data member assignment functions
//////////////////////////////////////////////////////////////////////

void cryatom::markatom() { marked = true; }
void cryatom::unmarkatom() { marked = false; }
void cryatom::setcrys (xtal* crys) { crystal = crys; }

void cryatom::setatom(RCSB::Atom* atom, const int& chn_index, const int& position, const int& res_index, const int& chn_type, const bool& ignoreAltLocFlag)
{
       _chn_index = chn_index;
       _position  = position;
       _res_index = res_index;
       _chn_type  = chn_type;
       ignoreAltLoc = ignoreAltLocFlag;
       if (crystal->orthogonal_block())
            setatom_orthogonal(atom);
       else setatom_fractional(atom);
}

// Determine the indices of the grid element containing the atom,
// and record them. Basically coord-in-fractional * number of
// elements per side of grid cell
void cryatom::setGridCoords()
{
       if (crystal->orthogonal_block()) {
            gridcoords[0]=int(floor((x - crystal->Xmin()) / crystal->grid_size()));
            gridcoords[1]=int(floor((y - crystal->Ymin()) / crystal->grid_size()));
            gridcoords[2]=int(floor((z - crystal->Zmin()) / crystal->grid_size()));
       } else {
            int elts_x = crystal->crysGrid()->eltsSide_x();
            int elts_y = crystal->crysGrid()->eltsSide_y();
            int elts_z = crystal->crysGrid()->eltsSide_z();
            gridcoords[0]=int(floor(x*elts_x));
            if (gridcoords[0] == elts_x) gridcoords[0]--;
            gridcoords[1]=int(floor(y*elts_y));
            if (gridcoords[1] == elts_y) gridcoords[1]--;
            gridcoords[2]=int(floor(z*elts_z));
            if (gridcoords[2] == elts_z) gridcoords[2]--;
       }
}

// Generate a symmetry related atom, given the symmetry operation.  
// Again, just a matrix multiplication and addition.
cryatom cryatom::symrel(const symop& currop)
{
       cryatom newat = *this;  // copy the current atom, but will reset its
                               // coordinates and symmetry id.

       newat.x = x * currop.op[0][0] + y * currop.op[0][1] + z * currop.op[0][2] + currop.op[0][3];
       newat.y = x * currop.op[1][0] + y * currop.op[1][1] + z * currop.op[1][2] + currop.op[1][3];
       newat.z = x * currop.op[2][0] + y * currop.op[2][1] + z * currop.op[2][2] + currop.op[2][3];
       newat.sym = currop.id;
       return newat;
}

// Determine the translation necessary to move the atom to the center 
// cell and apply it to the atom.  The translation is just the negative
// of the truncated coordinate (+1 if the coordinate is negative in the
// first place).  Think about it.  For example.  Suppose the coordinates
// are (1.5,-0.3,0.4), the translation is (-1,1,0) and the translated
// atom is at (0.5,0.7,0.4).
void cryatom::moveToCenterCell()
{
       if (crystal->orthogonal_block()) return;

       if (x>=0) transl[0] = int(floor(x) * (-1));
       else transl[0] = int(ceil(fabs(x)));
       if (y>=0) transl[1] = int(floor(y) * (-1));
       else transl[1] = int(ceil(fabs(y)));
       if (z>=0) transl[2] = int(floor(z) * (-1));
       else transl[2] = int(ceil(fabs(z)));
       x += transl[0];
       y += transl[1];
       z += transl[2];
}

// Identify the indices of the grid element containing the atom and
// then update the list of atoms in that grid element to contain this atom.
void cryatom::placeAtomInGrid()
{
       setGridCoords();
       crystal->crysGrid()->TheGrid()[gridcoords[0]][gridcoords[1]][gridcoords[2]].addmember(this);
}

void cryatom::neighborAtoms(std::list<CONTACT> &neighlist, const double& distlo, const double& disthi, int* cells, const int& remove_flag)
{
       if (crystal->orthogonal_block())
            _neighborAtoms_orthogonal(neighlist, distlo, disthi, cells, remove_flag);
       else _neighborAtoms(neighlist, distlo, disthi, cells, remove_flag);
}

void cryatom::neighborAtomsInCell(std::list<CONTACT> &neighlist, const double& dsqlimlo, const double& dsqlimhi, const gridcell& tocell,
                                  int* tr, const int& remove_flag)
{
       if (crystal->orthogonal_block())
            _neighborAtomsInCell(neighlist, dsqlimlo, dsqlimhi, tocell, remove_flag);
       else _neighborAtomsInCell(neighlist, dsqlimlo, dsqlimhi, tocell, tr, remove_flag);
}

//////////////////////////////////////////////////////////////////////
//            private member functions
//////////////////////////////////////////////////////////////////////
void cryatom::setatom_orthogonal(RCSB::Atom* atom)
{
       atm = atom;

       x = atom->orig().x;
       y = atom->orig().y;
       z = atom->orig().z;
}

void cryatom::setatom_fractional(RCSB::Atom* atom)
{
       atm = atom;

       x = atom->orig().x * crystal->xform[0][0] + atom->orig().y * crystal->xform[0][1] + atom->orig().z * crystal->xform[0][2] + crystal->xform[0][3];
       y = atom->orig().x * crystal->xform[1][0] + atom->orig().y * crystal->xform[1][1] + atom->orig().z * crystal->xform[1][2] + crystal->xform[1][3];
       z = atom->orig().x * crystal->xform[2][0] + atom->orig().y * crystal->xform[2][1] + atom->orig().z * crystal->xform[2][2] + crystal->xform[2][3];
}

// Compute the SQUARE distance of this atom to another atom.  
// This is useful when we just want to check whether a distance is in
// the desired range.  Save on computing that extra (expensive) sqrt
// for things out of the range.
double cryatom::dsq (const cryatom& other)
{ 
       double crysa=crystal->a;
       double crysb=crystal->b;
       double crysc=crystal->c;
    
       double delx = x - other.x;
       double dely = y - other.y;
       double delz = z - other.z;

       return (sqr(crysa*delx) + sqr (crysb*dely) + sqr (crysc*delz) + 2*(crysa*crysb*delx*dely*crystal->cosgamma) +
              2*(crysa*crysc*delx*delz*crystal->cosbeta) + 2*(crysb*crysc*dely*delz*crystal->cosalpha));
}

double cryatom::dsq_orthogonal(const cryatom& other)
{
       double delx = x - other.x;
       double dely = y - other.y;
       double delz = z - other.z;

       return (sqr(delx) + sqr(dely) + sqr(delz));
}

// This is the heart of the ROLODIST algorithm
// Given:  a Source atom, an upper and lower distance bound, and the number 
//         of grid elements that fit into the upper limit in each dimensions
// Find:   all atoms in the crystal that are within the given distance range
//         Return a list of pairs, with each std::pair containing the given
//            atom and one of the atoms in the distance range.
// I'll explain as I go along.
void cryatom::_neighborAtoms(std::list<CONTACT> &neighlist, const double& distlo, const double& disthi, int* cells, const int& remove_flag)
{
  double dsqlimlo = sqr(distlo);  // Square of lower bound of distance
  double dsqlimhi = sqr(disthi);  // Square of upper bound of distance
      //  Use these to check if candidate neighbor atom is in desired
      //  distance range, rather than actual distance...then we don't have
      //  to use sqrt unless the std::pair is actually in the range.
                            
  grid* gr = crystal->crysGrid();  // The grid representing the unit cell.

  int elts[3];
  elts[0] = gr->eltsSide_x(); // number of grid elements per side of the cell
  elts[1] = gr->eltsSide_y(); // number of grid elements per side of the cell
  elts[2] = gr->eltsSide_z(); // number of grid elements per side of the cell

  int startcell[3]; // indices of grid element to start checking for neighbors
                    // We start the search here, and scan linearly forward
                    // from here, to fill in a "cube" around the atom.
  int theCell[3];   // indices of grid element currently being checked 
                    //    for neighbors.
  int trhold[3]={0,0,0};  // translation which makes startcell reach the
                          // upper bound of distance.  
  int tr[3]={0,0,0};      // translation for currently scanned grid element

  gridcell* tocell = 0;    // grid element that is currently being scanned
  gridcell* nextcell = 0;  // next grid element to scan

  // Now, figure out where to start searching.  We already have the 
  // number of grid elements we need to search in each direction. So 
  // just start at the element that is that many elements *before* the
  // element containing the atom we are searching around, in each dimension.  
  // Since this may take us into a different unit cell, record the number 
  // of unit cells over this element is.
  for (int i = 0; i < 3; ++i) {
       startcell[i] = gridcoords[i]-cells[i];
       if (startcell[i] >= 0) continue;

       if (startcell[i]%elts[i]) {
            trhold[i] = int (floor ((double)((startcell[i]/elts[i]) - 1)));
            startcell[i] = elts[i] + (startcell[i]%elts[i]);
       } else {
            trhold[i] = int (floor ((double)((startcell[i]/elts[i]))));
            startcell[i] = 0;
       }
  }

  // Now search.  Basically cover the entire "cube" which extends the
  // appropriate number of grid elements to the left and to the right,
  // above and below, forward and behind the atom.  However, since we
  // want to be able to jump over the elements we already know have no
  // target atoms in them, we don't just count the elements we visit.
  // Instead, we maintain the total number of cells we have passed 
  // already in each direction (toti, totj, totk).  Thus, if we 
  // hop over three cells, we add three to the total, even though
  // we only visit one.
  theCell[2]=startcell[2];
  tr[2]=trhold[2];
  for (int totk = 0; totk < cells[2]*2 + 1;)
    {
    theCell[1]=startcell[1];  // Go to the beginning of a plane in the cube
    tr[1]=trhold[1];
    for (int totj = 0; totj < cells[1]*2 + 1;)
      {
      theCell[0]=startcell[0]; // Go to the beginning of a line in the plane
      tr[0]=trhold[0];
      for (int toti = 0; toti < cells[0]*2 + 1;)
      {
              tocell = gr->cellAtPos(theCell[0],theCell[1],theCell[2]);
              if (tocell->isMarked()) // If the current element has no
                  // target atoms in it, then there is nothing to search.
                {
                  // compute the set of Target atoms in this element
                  _neighborAtomsInCell(neighlist,dsqlimlo,dsqlimhi,*tocell,tr, remove_flag);
                }
              // Determine the next element to go to from the current one.
              // This is in the same line of the plane, may be the 
              // immediate neighbor or may require a hop.
              nextcell = gr->cellAtPos(tocell->nextx,theCell[1],theCell[2]);
              // If the next element doesn't have any Target atoms, 
              // then if we ever hit this element in a different search, 
              // we should hop past this next, and go directly to its
              // next.  So update this element's next.  
              if (!(nextcell->isMarked()) &&
                  (tocell->jumpx + nextcell->jumpx < elts[0]))
                {
                  // Make this element's next be the next element's next
                  tocell->nextx = nextcell->nextx;
                  // This requires a total hop of the hop that took
                  // this element to its next plus the hop that took
                  // the next element to its next element.
                  tocell->jumpx += nextcell->jumpx;
                  // And if the next element had to go into another unit cell
                  // cell to get to its next, then this element will also
                  // have to go into another unit cell.
                  tocell->jumpcellx += nextcell->jumpcellx;
                }  
              // Now go to the next element (which may be the updated one)
              theCell[0]=tocell->nextx; // set element x index to next x index
              tr[0] += tocell->jumpcellx; // go into next unit cell, in the 
                   // x direction, if that's where the next element is.
              toti += tocell->jumpx; // total distance updated by the hop size
            }
      theCell[1]=tocell->nexty;   // go to the next line in the plane
      tr[1] += tocell->jumpcelly; // which may be in the next unit cell in y
      totj += tocell->jumpy;  
      }
    tr[2] += tocell->jumpcellz; // go to the next plane in the cube
    theCell[2]=tocell->nextz;   // which may be in the next unit cell in z.
    totk += tocell->jumpz;
    }
}

void cryatom::_neighborAtoms_orthogonal(std::list<CONTACT> &neighlist, const double& distlo, const double& disthi, int* cells, const int& remove_flag)
{
       double dsqlimlo = sqr(distlo);  // Square of lower bound of distance
       double dsqlimhi = sqr(disthi);  // Square of upper bound of distance
                            
       grid* gr = crystal->crysGrid();  // The grid representing the unit cell.

       int elts_x = gr->eltsSide_x(); // number of grid elements per side of the cell

       int startcell[3]; // indices of grid element to start checking for neighbors
                         // We start the search here, and scan linearly forward
                         // from here, to fill in a "cube" around the atom.
       int theCell[3];   // indices of grid element currently being checked 
                         //    for neighbors.

       gridcell* tocell = 0;    // grid element that is currently being scanned
       gridcell* nextcell = 0;  // next grid element to scan

       // Now, figure out where to start searching.  We already have the 
       // number of grid elements we need to search in each direction. So 
       // just start at the element that is that many elements *before* the
       // element containing the atom we are searching around, in each dimension.  
       // Since this may take us into a different unit cell, record the number 
       // of unit cells over this element is.
       startcell[0] = gridcoords[0]-cells[0];
       if (startcell[0] < 0) startcell[0] = 0;
       startcell[1] = gridcoords[1]-cells[1];
       if (startcell[1] < 0) startcell[1] = 0;
       startcell[2] = gridcoords[2]-cells[2];
       if (startcell[2] < 0) startcell[2] = 0;

       // Now search.  Basically cover the entire "cube" which extends the
       // appropriate number of grid elements to the left and to the right,
       // above and below, forward and behind the atom.  However, since we
       // want to be able to jump over the elements we already know have no
       // target atoms in them, we don't just count the elements we visit.
       // Instead, we maintain the total number of cells we have passed 
       // already in each direction (toti, totj, totk).  Thus, if we 
       // hop over three cells, we add three to the total, even though
       // we only visit one.
       theCell[2]=startcell[2];
       for (int totk = 0; totk < cells[2]*2 + 1;) {
            theCell[1]=startcell[1];  // Go to the beginning of a plane in the cube
            for (int totj = 0; totj < cells[1]*2 + 1;) {
                theCell[0]=startcell[0]; // Go to the beginning of a line in the plane
                for (int toti = 0; toti < cells[0]*2 + 1;) {
                     tocell = gr->cellAtPos(theCell[0],theCell[1],theCell[2]);
                     if (tocell->isMarked()) {
                          // If the current element has no
                          // target atoms in it, then there is nothing to search.
                          // compute the set of Target atoms in this element
                          _neighborAtomsInCell(neighlist,dsqlimlo,dsqlimhi,*tocell, remove_flag);
                     }
                     // Determine the next element to go to from the current one.
                     // This is in the same line of the plane, may be the 
                     // immediate neighbor or may require a hop.
                     nextcell = gr->cellAtPos(tocell->nextx,theCell[1],theCell[2]);
                     // If the next element doesn't have any Target atoms, 
                     // then if we ever hit this element in a different search, 
                     // we should hop past this next, and go directly to its
                     // next.  So update this element's next.  
                     if (!(nextcell->isMarked()) &&
                         (tocell->jumpx + nextcell->jumpx < elts_x)) {
                          // Make this element's next be the next element's next
                          tocell->nextx = nextcell->nextx;
                          // This requires a total hop of the hop that took
                          // this element to its next plus the hop that took
                          // the next element to its next element.
                          tocell->jumpx += nextcell->jumpx;
                          // And if the next element had to go into another unit cell
                          // cell to get to its next, then this element will also
                          // have to go into another unit cell.
                          tocell->jumpcellx += nextcell->jumpcellx;
                     }
                     if (tocell->nextx < theCell[0]) break;
                     // Now go to the next element (which may be the updated one)
                     theCell[0]=tocell->nextx; // set element x index to next x index
                     // x direction, if that's where the next element is.
                     toti += tocell->jumpx; // total distance updated by the hop size
                }
                if (tocell->nexty < theCell[1]) break;
                theCell[1]=tocell->nexty;
                totj += tocell->jumpy;  
            }
            if (tocell->nextz < theCell[2]) break;
            theCell[2]=tocell->nextz;
            totk += tocell->jumpz;
       }
}

// This is the heart of the heart of the ROLODIST algorithm.
// Given:  a Source atom, upper and lower square distance limits, a grid
//         element and the translation that was used to reach this
//         grid element from the atom
// Find:   the set of target atoms contained in the grid element that 
//         are within the desired distance bounds.  Return a list of
//         pairs consisting of the source atom with each such target atom.
//
// This routine is a wee bit tricky because we want to maintain 
// enough information to return the correct translations that bring the
// std::pair into the desired distance range.
//
// Let's look at an example.
//

// Suppose that there  is a source  atom in the  input asymmetric unit at
// (1.1,-0.2,0.4).  This would have been stored in the grid as an atom at
// (0.1,0.8,0.4), with  translation  (-1,1,0).  (Clearly,  there  *is* an
// atom at (0.1,0.8,0.4), but with respect to the atom *IN THE ASYMMETRIC
// UNIT*,  it has  translation  (-1,1,0).  Now suppose   that there is  a
// target   atom at (0.95,0.9,0.5)     in   the asymmetric  unit.     Now
// (1.1,-0.2,0.4)  is  close to (0.95,0.9,0.5)  (because  it is  close to
// (0.95,-0.1,0.5), which  *also* must be present  in the crystal.)  We'd
// like  to be able to  report this  Source-Target pair,  with the target
// translation (0,-1,0).  But you would  never think that  (0.95,0.9,0.5)
// is close to  (0.1,0.8,0.4), which is how  the atoms are stored in  the
// grid!  So how do we figure this out???

// When we are exploring the  space around (0.1,0.8,0.4), we  essentially
// want to look at all the grid elements  in a cube with opposite corners
// (-0.1,0.6,0.2)   and (0.3,1.0,0.6).  But of   course  there is no such
// thing as (-0.2,0.6,0.2)  in our grid.  There is  just the grid element
// containing (0.8,0.6,0.2),    and  a translation  of    (-1,0,0).  Aha!
// (0.8,0.6,0.2) is starting to sound  closer  to (0.95,0.9,0.5).  As  we
// scan along,  we come   to  the grid  element  containing  the  atom at
// (0.95,0.9,0.5).  In order to calculate the distance, we must recognize
// that in order to reach this grid element  we essentially had to make a
// translation  of (-1,0,0).  So  the  distance  that should actually  be
// calculated is  between  (0.1,0.8,0.4) and  (0.95,0.9,0.5)+(-1,0,0)  or
// (-0.05,0.9,0.5).   This, indeed, is sufficient   to recognize that the
// two points  are close. HOWEVER.  There is  no atom at (0.1,0.8,0.4) in
// our original list of atoms,  and we want to  report distances of atoms
// in our initial   asymmetric unit to   target atoms in   or out of  the
// asymmetric unit.  So  we just translate the source  atom back to where
// it  was (i.e.  apply  the reverse  of  the translation that brought it
// into the center  cell) and apply the same   translation to the  target
// atom.  So the  distance between  (0.1,0.8,0.4) and (-0.05,0.9,0.5)  is
// the same as the  distance between (1.1,-0.2,0.4) and  (0.95,-0.1,0.4).
// So all we have to do is report the source atom, the target atom, and a
// net  translation  for the target   atom  of (1,-1,0) plus  (-1,0,0) or
// (0,-1,0).  
// Voila!
// (Of course, if some translation had been necessary to bring the target
// atom into the center cell, this would also have had to be added into
// the reported net translation.)

void cryatom::_neighborAtomsInCell(std::list<CONTACT> &neighlist, const double& dsqlimlo, const double& dsqlimhi, const gridcell& tocell,
                                  int* tr, const int& remove_flag)
{
       CONTACT contact;        // data structure to hold currently computed pair.
       int translated[3];
       const std::list<cryatom*>& tolist = tocell.membList();
       // i.e. for each atom in the list of target atoms
       //   for the grid element
       for (std::list<cryatom*>::const_iterator ptr = tolist.begin(); ptr != tolist.end(); ++ptr) {
            if (!((*ptr)->isMarked())) continue;
            for (int i = 0; i < 3; ++i) {
                 translated[i] = (*ptr)->transl[i] + tr[i] - transl[i];
            }
            if (!isValidPair(*(*ptr), translated, remove_flag)) continue;

            // candidate target atom, copy of one of the atoms in
            // the grid element, but translated.
            cryatom ToAtomTranslated = *(*ptr);

            // translate the candidate target atom by the amount the grid element has
            ToAtomTranslated.x += tr[0];  
            ToAtomTranslated.y += tr[1];  
            ToAtomTranslated.z += tr[2];  

            // Calculate the distance between the source atom in the central
            // cell and the translated grid element.  Check whether this distance
            // is in the desired range.
            double dsqval = dsq(ToAtomTranslated);
            if ((dsqval >= dsqlimlo) && (dsqval <= dsqlimhi)) {
                 // Calculate the net translation for the target atom:  
                 // its own translation to bring it into the center cell plus
                 // the translation of the grid element minus the translation to
                 // bring the source atom into the center cell.
                 ToAtomTranslated.transl[0] += (tr[0] - transl[0]);
                 ToAtomTranslated.transl[1] += (tr[1] - transl[1]);
                 ToAtomTranslated.transl[2] += (tr[2] - transl[2]);

                 contact.a_atm = atm;
                 contact.b_atm = ToAtomTranslated.atm;
                 contact.sym = ToAtomTranslated.sym;
                 contact.lx = ToAtomTranslated.transl[0];
                 contact.ly = ToAtomTranslated.transl[1];
                 contact.lz = ToAtomTranslated.transl[2];
                 contact.dist = sqrt(dsqval);
                 contact.a_res_index = _res_index;
                 contact.b_res_index = (*ptr)->_res_index;
                 contact.type = 0;
                 contact.mol_index = 0;
                 contact.chn_type = (*ptr)->_chn_type;
                 neighlist.push_back(contact);
            }
       }
}

// This is the heart of the heart of the ROLODIST algorithm.
// Given:  a Source atom, upper and lower square distance limits, a grid
//         element and the translation that was used to reach this
//         grid element from the atom
// Find:   the set of target atoms contained in the grid element that 
//         are within the desired distance bounds.  Return a list of
//         pairs consisting of the source atom with each such target atom.
//

void cryatom::_neighborAtomsInCell(std::list<CONTACT> &neighlist, const double& dsqlimlo, const double& dsqlimhi, const gridcell& tocell,
                                   const int& remove_flag)
{
       CONTACT contact;

       const std::list<cryatom*>& tolist = tocell.membList();
       for (std::list<cryatom*>::const_iterator ptr = tolist.begin(); ptr != tolist.end(); ++ptr) {
            if (!((*ptr)->isMarked())) continue;
            if (!isValidPair(*(*ptr), (*ptr)->transl, remove_flag)) continue;
            double dsqval = dsq_orthogonal(*(*ptr));
            if ((dsqval>=dsqlimlo) && (dsqval <= dsqlimhi)) {
                 contact.a_atm = atm;
                 contact.b_atm = (*ptr)->atm;
                 contact.sym = 0;
                 contact.lx = 0;
                 contact.ly = 0;
                 contact.lz = 0;
                 contact.dist = sqrt(dsqval);
                 contact.a_res_index = _res_index;
                 contact.b_res_index = (*ptr)->_res_index;
                 contact.type = 0;
                 contact.mol_index = 0;
                 contact.chn_type = (*ptr)->_chn_type;
                 neighlist.push_back(contact);
            }
       }
}

bool cryatom::isValidPair(const cryatom& toatom, int* translated, const int& remove_flag)
{
       // remove atom pairs between different conformation.
       if ((atm->alt_loc() != toatom.atm->alt_loc()) && !atm->alt_loc().empty() &&
           !toatom.atm->alt_loc().empty() && !ignoreAltLoc) return false;

       // has different symmetry, it's a valid pair
       if (translated[0] || translated[1] || translated[2] || sym != toatom.sym) return true;

       // remove atom pairs belong to same residue
       if ((remove_flag & SAME_RESIDUE_FLAG) && isSameResidue(toatom)) return false;

       // remove atom pairs belong to neighbor residues of same chain 
       if ((remove_flag & SAME_CHAIN_FLAG) && isLinkedResidue(toatom)) return false;

       // remove atom pairs belong to neighbor phosphate group contacts
       if (isNaPhosphateContact(toatom)) return false;
       return true;
}

bool cryatom::isSameResidue(const cryatom& toatom)
{
       if (toatom._res_index != _res_index) return false;
       return true;
}

bool cryatom::isLinkedResidue(const cryatom& toatom)
{
       if (!(_chn_type & POLYMER_TYPE) || !(toatom._chn_type & POLYMER_TYPE))
            return false;

       if (_chn_index != toatom._chn_index) return false;

       if (abs(_position - toatom._position) > 1) return false;

       if ((atm->atmtype() == "O3'" && toatom.atm->atmtype() == "P") ||
           (atm->atmtype() == "P" && toatom.atm->atmtype() == "O3'") ||
           (atm->atmtype() == "C" && toatom.atm->atmtype() == "N") ||
           (atm->atmtype() == "N" && toatom.atm->atmtype() == "C"))
            return true;

       return false;
}

bool cryatom::isNaPhosphateContact(const cryatom& toatom)
{
       if (!(_chn_type & POLYMER_TYPE) || !(toatom._chn_type & POLYMER_TYPE))
            return false;

       if (_chn_index != toatom._chn_index) return false;

       if (abs(_position - toatom._position) > 1) return false;

       if ((atm->atmtype() == "C3'" && toatom.atm->atmtype() == "P") || (atm->atmtype() == "P" && toatom.atm->atmtype() == "C3'") ||
           (atm->atmtype() == "O3'" && toatom.atm->atmtype() == "OP1") || (atm->atmtype() == "OP1" && toatom.atm->atmtype() == "O3'") ||
           (atm->atmtype() == "O3'" && toatom.atm->atmtype() == "OP2") || (atm->atmtype() == "OP2" && toatom.atm->atmtype() == "O3'") ||
           (atm->atmtype() == "O3'" && toatom.atm->atmtype() == "O5'") || (atm->atmtype() == "O5'" && toatom.atm->atmtype() == "O3'")) return true;

       return false;
}

void cryatom::checkAtom()
{
       printf("%s %s %s%s\n", atm->pdb_chnid().c_str(), atm->pdb_resnam().c_str(), atm->pdb_resnum().c_str(), atm->ins_code().c_str());
}
