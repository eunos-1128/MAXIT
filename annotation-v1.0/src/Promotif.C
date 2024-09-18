/*
FILE:     Promotif.C
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
#include <string.h>
#include <math.h>

#include "Promotif.h"
#include "VectorUtil.h"
#include "xtal.h"

#define AA_BACKBONE  3
#define NITROGEN     0
#define CARBON       1
#define OXYGEN       2
#define HYDROGEN     3

#define MAXENG    -0.5
#define ENGLIM    -9.9
#define KSQ1      0.42
#define KSQ2       0.2
#define KSF        332
#define MAXCAD     8.0
#define ACURCY  0.0001

static const char *aa_atoms[AA_BACKBONE] = {"N", "C", "O"};

Promotif::Promotif()
{
       _clear();
}

Promotif::~Promotif()
{
       _clear();
}

void Promotif::_clear()
{
       _logIo = NULL;
       _mol = NULL;
       _supportfile.clear();
       _residues.clear();
       _contacts.clear();
       _helices.clear();
       _sheets.clear();
}

void Promotif::setLog(LogUtil *logPt)
{
       _logIo = logPt;
}

void Promotif::setMolecule(RCSB::Molecule* mol)
{
       _mol = mol;
}

void Promotif::setSupportFile(const std::string& filename)
{
       _supportfile = filename;
}

const std::list<_HELIX>& Promotif::getHelices() const
{
       return _helices;
}

const std::list<_SHEET>& Promotif::getSheets() const
{
       return _sheets;
}

void Promotif::calculateSecondaryStructure(std::string& complicated_sheet_warning)
{
       complicated_sheet_warning.clear();

       if (_mol == NULL) return;

       _getResidues();

       if (_residues.empty()) return;

       _calContacts();
       _calEnergy();
       _assignSecondaryStruct();
       _assignHELIXandSHEET();
       _findHelices();
       _fineSheets(complicated_sheet_warning);
}

void Promotif::_getResidues()
{
       int num_residues = 0;
       RCSB::Chain* chain = _mol->GetFirstChain();
       while (chain) {
            if (chain->chain_type() != "ATOMP") {
                 chain = _mol->GetNextChain();
                 continue;
            }       
            RCSB::Residue* residue = chain->GetFirstResidue();
            while (residue) {
                 bool found = false;
                 for (int i = 0; i < AA_BACKBONE; ++i) {
                      if (residue->find_atom(aa_atoms[i])) {
                           found = true; 
                           break;
                      }
                 }       
                 if (found) num_residues++; 
                 residue = chain->GetNextResidue();
            }
            chain = _mol->GetNextChain();
       }
       if (!num_residues) return;

       _residues.clear();
       _residues.reserve(num_residues);

       _NRESIDUE nresidue;
       for (int i = 0; i < 4; ++i) {
            nresidue.energy[i] = 0;
            nresidue.hbond[i] = 0;
       }
       memset(nresidue.ssbond, 0, 12);
       for (int i = 0; i < 11; ++i) nresidue.ssbond[i] = ' ';
       for (int i = 0; i < 2; ++i) nresidue.bridpt[i] = 0;
       nresidue.ssbond1 = 0;
 
       COORD coord;
       chain = _mol->GetFirstChain();
       while (chain) {
            if (chain->chain_type() != "ATOMP") {
                 chain = _mol->GetNextChain();
                 continue;
            }
            RCSB::Residue* residue = chain->GetFirstResidue();
            while (residue) {
                 for (int i = 0; i < 4; ++i) nresidue.exist[i] = false;
                 bool found = false;
                 for (int i = 0; i < AA_BACKBONE; ++i) {
                      RCSB::Atom* atom = residue->find_atom(aa_atoms[i]);
                      if (atom) {
                           found = true;
                           nresidue.exist[i] = true;
                           nresidue.coord[i] = atom->orig();
                      }
                 }
                 if (!found) {
                      residue = chain->GetNextResidue();
                      continue;
                 }

                 nresidue.res = residue;
                 nresidue.chn = chain;
                 nresidue.chnbreak = false;
                 if (!_residues.empty()) {
                      int last = _residues.size() - 1;
                      if (!_residues[last].chnbreak) {
                           _residues[last].chnbreak = !is_connect("C", _residues[last].res, "N", nresidue.res, 2.5);
                           if (!_residues[last].chnbreak && _residues[last].exist[CARBON] && _residues[last].exist[OXYGEN] && nresidue.exist[NITROGEN]) {
                                vector_difference(coord, _residues[last].coord[CARBON], _residues[last].coord[OXYGEN]);
                                vector_normalize(coord);
                                vector_sum(nresidue.coord[HYDROGEN], nresidue.coord[NITROGEN], coord);
                                nresidue.exist[HYDROGEN] = 1;
                           }
                      }
                 }
                 _residues.push_back(nresidue);

                 residue = chain->GetNextResidue();
            }
            if (!_residues.empty()) _residues[_residues.size() - 1].chnbreak = true;

            chain = _mol->GetNextChain();
       }
}

void Promotif::_calContacts()
{
       xtal _xtal;
       _xtal.clear();

       double _Xmin = 10000.0;
       double _Ymin = 10000.0;
       double _Zmin = 10000.0;
       double _Xmax = -10000.0;
       double _Ymax = -10000.0;
       double _Zmax = -10000.0;

       for (unsigned int i = 0; i < _residues.size(); ++i) {
            RCSB::Atom* atom = _residues[i].res->find_atom("CA");
            if (!atom) continue;

            _xtal.add_atom(atom, 0, i, i, NONPOLY_TYPE);
            if (atom->orig().x < _Xmin) _Xmin = atom->orig().x;
            if (atom->orig().x > _Xmax) _Xmax = atom->orig().x;
            if (atom->orig().y < _Ymin) _Ymin = atom->orig().y;
            if (atom->orig().y > _Ymax) _Ymax = atom->orig().y;
            if (atom->orig().z < _Zmin) _Zmin = atom->orig().z;
            if (atom->orig().z > _Zmax) _Zmax = atom->orig().z;
       }

       _xtal.set_grid_size(MAXCAD);
       _xtal.set_block(_Xmin, _Ymin, _Zmin, _Xmax, _Ymax, _Zmax);
       _xtal.genAtoms();
       _xtal.placeAtomsInGrid();

       std::list<CONTACT> pl;
       _xtal.getContactList(pl, "all_asym", "all_asym", 0.85, MAXCAD, 0);

       if (pl.empty()) return;

       std::set<unsigned int> tmp_set;
       for (std::list<CONTACT>::const_iterator pos = pl.begin(); pos != pl.end(); ++pos) {
            std::map<unsigned int, std::set<unsigned int> >::iterator mpos = _contacts.find(pos->a_res_index);
            if (mpos != _contacts.end()) mpos->second.insert(pos->b_res_index);
            else {
                 tmp_set.clear();
                 tmp_set.insert(pos->b_res_index);
                 _contacts.insert(std::make_pair(pos->a_res_index, tmp_set));
            }
            mpos = _contacts.find(pos->b_res_index);
            if (mpos != _contacts.end()) mpos->second.insert(pos->a_res_index);
            else {
                 tmp_set.clear();
                 tmp_set.insert(pos->a_res_index);
                 _contacts.insert(std::make_pair(pos->b_res_index, tmp_set));
            }
       }
}

void Promotif::_calEnergy()
{
       for (unsigned int i = 0; i < _residues.size(); ++i) {
            std::map<unsigned int, std::set<unsigned int> >::const_iterator mpos = _contacts.find(i);
            if (mpos == _contacts.end()) continue;
            if (!_residues[i].exist[NITROGEN] || !_residues[i].exist[HYDROGEN] || _residues[i].res->ResName() == "PRO") continue;
            for (unsigned int j = 0; j < _residues.size(); ++j) {
                 if (mpos->second.find(j) == mpos->second.end()) continue;
                 if (!_residues[j].exist[CARBON] || !_residues[j].exist[OXYGEN]) continue;
                 if ((abs((int) i - (int) j) == 1 && _residues[i].chn != _residues[j].chn) || abs((int) i - (int) j) >= 2) {
                      double don = vector_distance(_residues[j].coord[OXYGEN], _residues[i].coord[NITROGEN]);
                      double doh = vector_distance(_residues[j].coord[OXYGEN], _residues[i].coord[HYDROGEN]);
                      double dch = vector_distance(_residues[j].coord[CARBON], _residues[i].coord[HYDROGEN]);
                      double dcn = vector_distance(_residues[j].coord[CARBON], _residues[i].coord[NITROGEN]);
                      if (!_xeqy(don, 0.0) && !_xeqy(doh, 0.0) && !_xeqy(dch, 0.0) && !_xeqy(dcn, 0.0)) {
                           double energy = KSQ1 * KSQ2 * KSF * (1.0 / don + 1.0 / dch - 1.0 / doh - 1.0 / dcn);
                           if (energy < ENGLIM) energy = ENGLIM;
                           if (energy < MAXENG) _engset(energy, i, j);
                      }
                 }
            }
       }
}

bool Promotif::_xeqy(const double& x, const double& y)
{
       return (fabs(x - y) < ACURCY);
}

void Promotif::_engset(const double& energy, const int& donres, const int& accres)
{
       int donpl = _engpl(energy, &(_residues[donres].energy[2]));
       int acpl  = _engpl(energy, &(_residues[accres].energy[0]));
       if (donpl && acpl) {
            if (_residues[donres].hbond[3]) {
                 int rejptr = _residues[donres].hbond[3];
                 _rejbnd(donres, _residues[rejptr - 1], 2);
            }
            if (_residues[accres].hbond[1]) {
                 int rejptr = _residues[accres].hbond[1];
                 _rejbnd(accres, _residues[rejptr - 1], 0);
            }
            if (acpl == 1) {
                 _residues[accres].energy[1] = _residues[accres].energy[0];
                 _residues[accres].hbond[1]  = _residues[accres].hbond[0];
            }
            if (donpl == 1) {
                 _residues[donres].energy[3] = _residues[donres].energy[2];
                 _residues[donres].hbond[3]  = _residues[donres].hbond[2];
            }
            _residues[accres].energy[acpl -1] = energy;
            _residues[accres].hbond[acpl - 1] = donres + 1;
            _residues[donres].energy[donpl+1] = energy;
            _residues[donres].hbond[donpl+ 1] = accres + 1;
       }
}

int Promotif::_engpl(const double& energy, double *honde)
{
       int index = 0;
       if (energy < honde[0]) {
            index = 1;
       } else {
            if (energy < honde[1]) index = 2;
       }
       return index;
}

void Promotif::_rejbnd(const int& rejptr, _NRESIDUE& residue, const int& start)
{
       int rejpl = 1;
       if (residue.hbond[start] != (rejptr + 1)) rejpl++;
       if (rejpl == 1) {
            residue.hbond[start] = residue.hbond[start + 1];
            residue.energy[start] = residue.energy[start + 1];
       }
       residue.hbond[start + 1] = 0;
       residue.energy[start + 1] = 0;
}
