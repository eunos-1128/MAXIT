/*
FILE:     Pdb2Ndb_Struct_Features.C
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

#include "Maxit.h"
#include "NdbToken.h"
#include "utillib.h"

void Maxit::Read_PDB_StructuralFeatures()
{
       _helix.clear();
       _sheet.clear();
       _turn.clear();
       _sltbrgs.clear();
       _modres.clear();
       _site.clear();

       if (_molecules.empty()) return;

       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = _pdb_records.find("HELIX");
       if (mpos != _pdb_records.end()) _pdb_to_ndb_read_HELIX(mpos->second);

       mpos = _pdb_records.find("SHEET");
       if (mpos != _pdb_records.end()) _pdb_to_ndb_read_SHEET(mpos->second);

       mpos = _pdb_records.find("TURN");
       if (mpos != _pdb_records.end()) _pdb_to_ndb_read_TURN(mpos->second);

       InputPDBSSBONDRecord();
       InputPDBLINKRecord();
       _pdb_to_ndb_read_SLTBRG();
       _pdb_to_ndb_read_MODRES();
       _pdb_to_ndb_read_SITE();
}

void Maxit::_pdb_to_ndb_read_HELIX(const std::list<std::vector<std::string> >& records)
{
       if (records.empty()) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > names;
       std::vector<RCSB::Residue*> residues;

       _HELIX helix;
       helix.mol_index = _molecules[0]->index();
       helix.ID.clear();
       helix.initRes = -1;
       helix.endRes = -1;
       helix.helixClass = 0;
       helix.comment.clear();

       for (std::list<std::vector<std::string> >::const_iterator pos = records.begin(); pos != records.end(); ++pos) {
            names.clear();
            data.clear();
            data.push_back((*pos)[4]);
            data.push_back((*pos)[3]);
            data.push_back((*pos)[5]);
            data.push_back((*pos)[6]);
            names.push_back(data);
            data.clear();
            data.push_back((*pos)[8]);
            data.push_back((*pos)[7]);
            data.push_back((*pos)[9]);
            data.push_back((*pos)[10]);
            names.push_back(data);

            _getResiduesfromNames(0, names, helix.mol_index, residues, data);
            if (residues.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 continue;
            }

            helix.ID = (*pos)[2];
            helix.initRes = residues[0]->index();
            helix.endRes  = residues[1]->index();
            helix.helixClass = atoi((*pos)[11].c_str());
            helix.comment = (*pos)[12];
            _helix.push_back(helix);
       }
}

void Maxit::_pdb_to_ndb_read_SHEET(const std::list<std::vector<std::string> >& records)
{
       if (records.empty() || _molecules.empty()) return;

       std::vector<std::string> sheet_ids;
       sheet_ids.clear();
       std::set<std::string> error_sheet_ids;
       error_sheet_ids.clear();
       std::map<std::string, std::pair<int, std::vector<_SHEET_STRAND> > > strands_map;
       strands_map.clear();

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > names;
       std::vector<RCSB::Residue*> residues;

       _SHEET_STRAND strand;
       std::vector<_SHEET_STRAND> strands;

       for (std::list<std::vector<std::string> >::const_iterator pos = records.begin(); pos != records.end(); ++pos) {
            names.clear();
            data.clear();
            data.push_back((*pos)[5]);
            data.push_back((*pos)[4]);
            data.push_back((*pos)[6]);
            data.push_back((*pos)[7]);
            names.push_back(data);
            data.clear();
            data.push_back((*pos)[9]);
            data.push_back((*pos)[8]);
            data.push_back((*pos)[10]);
            data.push_back((*pos)[11]);
            names.push_back(data);
            if (!(*pos)[14].empty() && !(*pos)[15].empty() && !(*pos)[16].empty() && !(*pos)[19].empty() && !(*pos)[20].empty() && !(*pos)[21].empty()) {
                 data.clear();
                 data.push_back((*pos)[15]);
                 data.push_back((*pos)[14]);
                 data.push_back((*pos)[16]);
                 data.push_back((*pos)[17]);
                 names.push_back(data);
                 data.clear();
                 data.push_back((*pos)[20]);
                 data.push_back((*pos)[19]);
                 data.push_back((*pos)[21]);
                 data.push_back((*pos)[22]);
                 names.push_back(data);
            }

            int mol_index = _molecules[0]->index();
            _getResiduesfromNames(0, names, mol_index, residues, data);
            if (residues.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 error_sheet_ids.insert((*pos)[2]);
                 continue;
            }

            strand.strand_id = atoi((*pos)[1].c_str());
            strand.begin_res_index = residues[0]->index();
            strand.end_res_index = residues[1]->index();
            strand.sense = atoi((*pos)[12].c_str());
            if (residues.size() == 4) {
                 strand.curr_hbond_res_index = residues[2]->index();
                 strand.prev_hbond_res_index = residues[3]->index();
                 strand.curr_hbond_atom_name = (*pos)[13];
                 strand.prev_hbond_atom_name = (*pos)[18];
            } else {
                 strand.curr_hbond_atom_name.clear();
                 strand.curr_hbond_res_index = -1;
                 strand.prev_hbond_atom_name.clear();
                 strand.prev_hbond_res_index = -1;
            }

            std::map<std::string, std::pair<int, std::vector<_SHEET_STRAND> > >::iterator mpos = strands_map.find((*pos)[2]);
            if (mpos != strands_map.end()) mpos->second.second.push_back(strand);
            else {
                 strands.clear();
                 strands.push_back(strand);
                 strands_map.insert(std::make_pair((*pos)[2], std::make_pair(atoi((*pos)[3].c_str()), strands)));
                 sheet_ids.push_back((*pos)[2]);
            }
       }

       _SHEET_TOPOLOGY strand_order;
       _SHEET sheet;
       for (std::vector<std::string>::const_iterator vpos = sheet_ids.begin(); vpos != sheet_ids.end(); ++vpos) {
            if (error_sheet_ids.find(*vpos) != error_sheet_ids.end()) continue;
            std::map<std::string, std::pair<int, std::vector<_SHEET_STRAND> > >::const_iterator mpos = strands_map.find(*vpos);
            if (mpos == strands_map.end()) continue;

            sheet.mol_index = _molecules[0]->index();
            sheet.sheetID = *vpos;
            sheet.numStrands = mpos->second.first;
            sheet.complicateFlag = false;
            sheet._strands = mpos->second.second;
            sheet._strand_orders.clear();

            for (unsigned int i = 1; i < sheet._strands.size(); ++i) {
                 if ((sheet._strands[i].curr_hbond_res_index < 0) || (sheet._strands[i].prev_hbond_res_index < 0)) continue;

                 strand_order.sense_number = sheet._strands[i].sense;
                 strand_order.sense_string.clear();
                 if (strand_order.sense_number == 1)
                      strand_order.sense_string = "parallel";
                 else if (strand_order.sense_number == (-1))
                      strand_order.sense_string = "anti-parallel";
                 strand_order.range_id_1 = String::IntToString(sheet._strands[i-1].strand_id);
                 strand_order.range_id_1_res_index = sheet._strands[i].prev_hbond_res_index;
                 strand_order.range_id_1_atom_name = sheet._strands[i].prev_hbond_atom_name;
                 strand_order.range_id_2 = String::IntToString(sheet._strands[i].strand_id);
                 strand_order.range_id_2_res_index = sheet._strands[i].curr_hbond_res_index;
                 strand_order.range_id_2_atom_name = sheet._strands[i].curr_hbond_atom_name;
                 sheet._strand_orders.push_back(strand_order);
            }
            _sheet.push_back(sheet);
       }
}

void Maxit::_pdb_to_ndb_read_TURN(const std::list<std::vector<std::string> >& records)
{
       if (records.empty()) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > names;
       std::vector<RCSB::Residue*> residues;

       _HELIX turn;
       turn.mol_index = _molecules[0]->index();
       turn.ID.clear();
       turn.initRes = -1;
       turn.endRes = -1;
       turn.helixClass = 0;
       turn.comment.clear();

       for (std::list<std::vector<std::string> >::const_iterator pos = records.begin(); pos != records.end(); ++pos) {
            names.clear();
            data.clear();
            data.push_back((*pos)[4]);
            data.push_back((*pos)[3]);
            data.push_back((*pos)[5]);
            data.push_back((*pos)[6]);
            names.push_back(data);
            data.clear();
            data.push_back((*pos)[8]);
            data.push_back((*pos)[7]);
            data.push_back((*pos)[9]);
            data.push_back((*pos)[10]);
            names.push_back(data);

            _getResiduesfromNames(0, names, turn.mol_index, residues, data);
            if (residues.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 continue;
            }

            turn.ID = (*pos)[2];
            turn.initRes = residues[0]->index();
            turn.endRes  = residues[1]->index();
            turn.comment = (*pos)[11];
            _turn.push_back(turn);
       }
}

void Maxit::_pdb_to_ndb_read_SLTBRG()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator mpos = _pdb_records.find("SLTBRG");
       if (mpos == _pdb_records.end()) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > names;
       std::vector<RCSB::Atom*> atoms;
       _SLTBRG sltbrg;

       for (std::list<std::vector<std::string> >::const_iterator pos = mpos->second.begin(); pos != mpos->second.end(); ++pos) {
            names.clear();
            data.clear();
            data.push_back((*pos)[4]);
            data.push_back((*pos)[3]);
            data.push_back((*pos)[5]);
            data.push_back((*pos)[6]);
            data.push_back((*pos)[1]);
            data.push_back((*pos)[2]);
            names.push_back(data);
            data.clear();
            data.push_back((*pos)[10]);
            data.push_back((*pos)[9]);
            data.push_back((*pos)[11]);
            data.push_back((*pos)[12]);
            data.push_back((*pos)[7]);
            data.push_back((*pos)[8]);
            names.push_back(data);

            _getAtomsfromNames(0, names, sltbrg.mol_index, atoms, data);
            if (atoms.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 continue;
            }

            sltbrg.fstAtom = atoms[0];
            sltbrg.sndAtom = atoms[1];
            sltbrg.SymOP_1 = _reformat_symmetry((*pos)[13]);
            sltbrg.SymOP_2 = _reformat_symmetry((*pos)[14]);
            _sltbrgs.push_back(sltbrg);
       }
}

void Maxit::_pdb_to_ndb_read_MODRES()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("MODRES");
       if (ppos == _pdb_records.end()) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > names;
       std::vector<RCSB::Residue*> residues;

       _MODRES modres;
       for (std::list<std::vector<std::string> >::const_iterator pos = ppos->second.begin(); pos != ppos->second.end(); ++pos) {
            names.clear();
            data.clear();
            data.push_back((*pos)[3]);
            data.push_back((*pos)[2]);
            data.push_back((*pos)[4]);
            data.push_back((*pos)[5]);
            names.push_back(data);

            _getResiduesfromNames(0, names, modres.mol_index, residues, data);
            if (residues.empty()) {
                 for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                      _logIo->message("%s\n", vpos->c_str());
                 }
                 continue;
            }

            modres.res_index = residues[0]->index();
            modres.Standard_Name = (*pos)[6];
            modres.details = (*pos)[7];
            _modres.push_back(modres);
       }
}

void Maxit::_pdb_to_ndb_read_SITE()
{
       std::map<std::string, std::list<std::vector<std::string> > >::const_iterator ppos = _pdb_records.find("SITE");
       if (ppos == _pdb_records.end()) return;

       std::vector<std::string> data;
       std::vector<std::vector<std::string> > names;
       std::vector<RCSB::Residue*> residues;

       _SITE site;

       std::vector<std::string> site_ids;
       site_ids.clear();
       std::map<std::string, _SITE> site_mapping;
       site_mapping.clear();

       for (std::list<std::vector<std::string> >::const_iterator pos = ppos->second.begin(); pos != ppos->second.end(); ++pos) {
            std::string site_id = (*pos)[2];
            std::map<std::string, _SITE>::iterator spos = site_mapping.find(site_id);
            if (spos == site_mapping.end()) { 
                 site.mol_index = 0;
                 site.SiteID = site_id;
                 site.code.clear();
                 site.details.clear();
                 site.self_chain_indices.clear();
                 site.self_res_indices.clear();
                 site.associated_residues.clear();

                 site_mapping.insert(std::make_pair(site_id, site));
                 site_ids.push_back(site_id);
                 spos = site_mapping.find(site_id);
            }
            for (int siteId = 0; siteId < 4; ++siteId) {
                 if ((*pos)[siteId * 4 + 4].empty() || (*pos)[siteId * 4 + 5].empty() ||
                     (*pos)[siteId * 4 + 6].empty()) continue;

                 names.clear();
                 data.clear();
                 data.push_back((*pos)[siteId * 4 + 5]);
                 data.push_back((*pos)[siteId * 4 + 4]);
                 data.push_back((*pos)[siteId * 4 + 6]);
                 data.push_back((*pos)[siteId * 4 + 7]);
                 names.push_back(data);

                 _getResiduesfromNames(0, names, site.mol_index, residues, data);
                 if (residues.empty()) {
                      for (std::vector<std::string>::const_iterator vpos = data.begin(); vpos != data.end(); ++vpos) {
                           _logIo->message("%s\n", vpos->c_str());
                      }
                      continue;
                 }

                 spos->second.associated_residues.push_back(std::make_pair(residues[0]->index(), ""));
            }
       }

       for (std::vector<std::string>::const_iterator vpos = site_ids.begin(); vpos != site_ids.end(); ++vpos) {
            std::map<std::string, _SITE>::const_iterator
                mpos = site_mapping.find(*vpos);
            if (mpos == site_mapping.end()) continue;
            if (mpos->second.associated_residues.empty()) continue;
            _site.push_back(mpos->second);
       }
}
