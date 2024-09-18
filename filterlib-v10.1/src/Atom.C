/*
FILE:     Atom.C
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

#include "Atom.h"
#include "Atom_global.h"
#include "CompositeIndex.h"
#include "utillib.h"

using namespace RCSB;

#define FORMAT  "ATOM  %5s %-4s%c%-3s%2s%4s%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n"

void Atom::_clear_data()
{
       _atom.clear();
       _empty_flag.clear();
       _sigatm.clear();
       _anisou.clear();
       _siguij.clear();
       // _deleted = false;
       _ter_flag = 0;
       _lineno = 0;
       _is_9999_flag = false;
       _empty.clear();
       _atom_item_position_mapping.clear();
       _anisotrop_item_position_mapping.clear();
       _sigatm_item_position_mapping.clear();
       _anisou_item_position_mapping.clear();
       _siguij_item_position_mapping.clear();
       _extra_item_value_mapping.clear();
}

void Atom::clear()
{
       _clear_data();

       for (short int i = 0; i <= FORMAT_LONG; ++i) {
            _atom.push_back("");
            _empty_flag.push_back(true);
            _atom_item_position_mapping.insert(std::make_pair(_format_atom_site_pdbx[i], i));
       }
}

Atom& Atom::operator=(const Atom& atom)
{
       if (this != &atom) {
            _clear_data();
            _orig = atom._orig;
            _atom = atom._atom;
            _empty_flag = atom._empty_flag;
            _sigatm = atom._sigatm;
            _anisou = atom._anisou;
            _siguij = atom._siguij;
            _ter_flag = atom._ter_flag;
            _lineno = atom._lineno;
            _is_9999_flag = atom._is_9999_flag;
            _empty = atom._empty;
            _atom_item_position_mapping = atom._atom_item_position_mapping;
            _anisotrop_item_position_mapping = atom._anisotrop_item_position_mapping;
            _sigatm_item_position_mapping = atom._sigatm_item_position_mapping;
            _anisou_item_position_mapping = atom._anisou_item_position_mapping;
            _siguij_item_position_mapping = atom._siguij_item_position_mapping;
       }
       return (*this);
}

const std::string& Atom::type() const       { return _atom[0];  }
const std::string& Atom::atnum() const      { return _atom[1];  }
const std::string& Atom::atmtype() const    { return _atom[2];  }
const std::string& Atom::alt_loc() const    { return _atom[3];  }
const std::string Atom::alt_loc_char() const
{
       if (!_atom[3].empty())
            return _atom[3].substr(0, 1);
       else return " ";
}
const std::string& Atom::restype() const    { return _atom[4];  }
const std::string& Atom::chnid() const      { return _atom[5];  }
const std::string& Atom::resnum() const     { return _atom[6];  }
const std::string& Atom::ins_code() const   { return _atom[7];  }
const std::string Atom::ins_code_char() const
{
       if (!_atom[7].empty())
            return _atom[7].substr(0, 1);
       else return " ";
}
const std::string& Atom::occ() const        { return _atom[11]; }
const std::string& Atom::t_fct() const      { return _atom[12]; }
const std::string& Atom::atom_type() const  { return _atom[13]; }
const std::string& Atom::charge() const     { return _atom[14]; }
const std::string& Atom::pdb_resnum() const { return _atom[15]; }
const std::string& Atom::pdb_resnam() const { return _atom[16]; }
const std::string& Atom::pdb_chnid() const  { return _atom[17]; }
const std::string Atom::pdb_chnid_char() const
{
       if (!_atom[17].empty())
            return _atom[17].substr(0, 1);
       else return " ";
}
const std::string& Atom::pdb_atmnam() const { return _atom[18]; }
const std::string& Atom::pub_resnum() const { return _atom[19]; }
const std::string& Atom::pub_resnam() const { return _atom[20]; }
const std::string& Atom::pub_chnid() const  { return _atom[21]; }
const std::string& Atom::pub_atmnam() const { return _atom[22]; }
const COORD& Atom::orig() const             { return _orig; }
const double& Atom::x() const               { return _orig.x; }
const double& Atom::y() const               { return _orig.y; }
const double& Atom::z() const               { return _orig.z; }
const bool Atom::is_hydrogen() const        { return ((_atom[13] == "H") || (_atom[13] == "D")); }
const bool Atom::has_sigatm() const         { return (!_sigatm.empty()); }
const bool Atom::has_anisou() const         { return (!_anisou.empty()); }
const bool Atom::has_siguij() const         { return (!_siguij.empty()); }
// const bool& Atom::deleted() const           { return _deleted; }
const int& Atom::ter_flag() const           { return _ter_flag; }
const int& Atom::lineno() const             { return _lineno; }
const bool& Atom::is_9999_flag() const      { return _is_9999_flag; }

const std::string Atom::getAtomIndex() const
{
       return CompositeIndex::getIndex(_atom[4], _atom[2]);
}

const std::string Atom::getAtomAltIndex() const
{
       return CompositeIndex::getIndex(_atom[4], _atom[2], _atom[3]);
}

const std::string Atom::getAtomAllIndex(const bool& alt_flag) const
{
      std::vector<std::string> data;
      data.clear();
      data.push_back(_atom[17]);
      data.push_back(_atom[16]);
      data.push_back(_atom[15]);
      data.push_back(_atom[7]);
      data.push_back(_atom[18]);
      if (alt_flag) data.push_back(_atom[3]);
      return CompositeIndex::getIndex(data);
}

const std::string Atom::getAtomAllIndexWithoutAltId() const
{
      std::vector<std::string> data;
      data.clear();
      data.push_back(_atom[17]);
      data.push_back(_atom[16]);
      data.push_back(_atom[15]);
      data.push_back(_atom[7]);
      data.push_back(_atom[18]);
      return CompositeIndex::getIndex(data);
}

const std::string& Atom::getValue(const int& pos) const
{
       if ((pos >= 0) && (pos < FORMAT_LONG)) return _atom[pos];
       return _empty;
}

const std::string& Atom::getOrigValue(const int& pos) const
{
       if ((pos >= 0) && (pos < FORMAT_LONG) && !_empty_flag[pos]) return _atom[pos];
       return _empty;
}

const std::string& Atom::getValue(const std::string& item)
{
       std::map<std::string, short int>::const_iterator pos = _atom_item_position_mapping.find(item);
       if (pos != _atom_item_position_mapping.end()) return _atom[pos->second]; 

       pos = _sigatm_item_position_mapping.find(item);
       if (pos != _sigatm_item_position_mapping.end()) return _sigatm[pos->second];

       std::map<std::string, std::string>::const_iterator mpos = _extra_item_value_mapping.find(item);
       if (mpos != _extra_item_value_mapping.end()) return mpos->second;

       return _empty;
}

const std::string& Atom::getAnisouValue(const std::string& item)
{
       std::map<std::string, short int>::const_iterator pos = _anisotrop_item_position_mapping.find(item);
       if (pos != _anisotrop_item_position_mapping.end()) return _atom[pos->second]; 

       pos = _anisou_item_position_mapping.find(item);
       if (pos != _anisou_item_position_mapping.end()) return _anisou[pos->second];

       pos = _siguij_item_position_mapping.find(item);
       if (pos != _siguij_item_position_mapping.end()) return _siguij[pos->second];

       return _empty;
}

const std::vector<std::string>& Atom::getValue() const
{
       return _atom;
}

void Atom::getAuxiliaryValue(const std::string& type, std::vector<std::string>& vals)
{
       vals.clear();
       if ((type == "SIGATM") && !_sigatm.empty()) {
            vals = _atom;
            vals[0] = "SIGATM";
            vals[8]  = _sigatm[0];
            vals[9]  = _sigatm[1];
            vals[10] = _sigatm[2];
            vals[11] = _sigatm[3];
            vals[12] = _sigatm[4];
       } else if ((type == "ANISOU") && !_anisou.empty()) {
            vals = _atom;
            vals[0] = "ANISOU";
            vals[14] = "";
            vals[8]  = String::IntToString((int) (rint(10000 * atof(_anisou[0].c_str()))));
            vals[9]  = String::IntToString((int) (rint(10000 * atof(_anisou[1].c_str()))));
            vals[10] = String::IntToString((int) (rint(10000 * atof(_anisou[2].c_str()))));
            vals[11] = String::IntToString((int) (rint(10000 * atof(_anisou[3].c_str()))));
            vals[12] = String::IntToString((int) (rint(10000 * atof(_anisou[4].c_str()))));
            vals[13] = String::IntToString((int) (rint(10000 * atof(_anisou[5].c_str()))));
       } else if ((type == "SIGUIJ") && !_siguij.empty()) {
            vals = _atom;
            vals[0] = "SIGUIJ";
            vals[14] = "";
            vals[8]  = String::IntToString((int) (rint(10000 * atof(_siguij[0].c_str()))));
            vals[9]  = String::IntToString((int) (rint(10000 * atof(_siguij[1].c_str()))));
            vals[10] = String::IntToString((int) (rint(10000 * atof(_siguij[2].c_str()))));
            vals[11] = String::IntToString((int) (rint(10000 * atof(_siguij[3].c_str()))));
            vals[12] = String::IntToString((int) (rint(10000 * atof(_siguij[4].c_str()))));
            vals[13] = String::IntToString((int) (rint(10000 * atof(_siguij[5].c_str()))));
       }
}

void Atom::getAtomIndex(std::vector<std::string>& data, const bool& clear_flag, const bool& exclude_alt_loc_flag)
{
       if (clear_flag) data.clear();
       data.push_back(_atom[17]); // pdb_chnid()
       data.push_back(_atom[16]); // pdb_resnam()
       data.push_back(_atom[15]); // pdb_resnum()
       data.push_back(_atom[7]);  // ins_code()
       data.push_back(_atom[18]); // pdb_atmnam()
       if (exclude_alt_loc_flag) return;
       data.push_back(_atom[3]);  // alt_loc()
}

void Atom::getResidueIndex(std::vector<std::string>& data, const bool& clear_flag)
{
       if (clear_flag) data.clear();
       data.push_back(_atom[17]); // pdb_chnid()
       data.push_back(_atom[16]); // pdb_resnam()
       data.push_back(_atom[15]); // pdb_resnum()
       data.push_back(_atom[7]);  // ins_code()
}

void Atom::set_type(const std::string& a)       { _atom[0] = a;  }
void Atom::set_atnum(const int& a)              { _atom[1] = String::IntToString(a); }
void Atom::set_atmtype(const std::string& a)    { _atom[2] = a;  }
void Atom::set_alt_loc(const std::string& a)    { _atom[3] = a;  }
void Atom::set_restype(const std::string& a)    { _atom[4] = a;  }
void Atom::set_chnid(const std::string& a)      { _atom[5] = a;  }
void Atom::set_resnum(const std::string& a)     { _atom[6] = a; }
void Atom::set_ins_code(const std::string& a)   { _atom[7] = a;  }
void Atom::set_atom_type(const std::string& a)  { _atom[13] = a; }
void Atom::set_charge(const std::string& a)     { _atom[14] = a; }
void Atom::set_pdb_resnum(const std::string& a) { _atom[15] = a; }
void Atom::set_pdb_resnam(const std::string& a) { _atom[16] = a; }
void Atom::set_pdb_chnid(const std::string& a)  { _atom[17] = a; }
void Atom::set_pdb_atmnam(const std::string& a) { _atom[18] = a; }
// void Atom::set_deletion()                       { _deleted  = true; }
void Atom::set_ter_flag(const int& lineno)      { _ter_flag = lineno; }

void Atom::set_orig(const COORD& a)
{
       _orig = a;
       _atom[8] = FloatToString(_orig.x, 0, 3);
       _atom[9] = FloatToString(_orig.y, 0, 3);
       _atom[10]= FloatToString(_orig.z, 0, 3);
}

void Atom::set_x(const double& x)
{
       _orig.x = x;
       _atom[8] = FloatToString(_orig.x, 0, 3);
}

void Atom::set_y(const double& y)
{
       _orig.y = y;
       _atom[9] = FloatToString(_orig.y, 0, 3);
}

void Atom::set_z(const double& z)
{
       _orig.z = z;
       _atom[10]= FloatToString(_orig.z, 0, 3);
}

void Atom::setValue(const std::string& val, const int& pos)
{
       if ((pos >= 0) && (pos <= FORMAT_LONG)) _atom[pos] = val;
}

void Atom::setValue(const std::vector<std::string>& vals, const int& lineno)
{
       _lineno = lineno;
       for (unsigned int i = 0; i < _atom.size(); ++i) {
            if (i == vals.size()) break;
            _atom[i] = vals[i];
            if (!_atom[i].empty() && (_atom[i] != "?") && (_atom[i] != ".")) _empty_flag[i] = false;
       }
       if (_atom[11].empty()) _atom[11] = "1.00";
       if (_atom[11][0] == '.') _atom[11] = "0" + _atom[11];
       if (!_atom[13].empty()) String::UpperCase(_atom[13]);
       if (_atom[15].empty()) _atom[15] = _atom[6];
       if (_atom[16].empty()) _atom[16] = _atom[4];
       if (_atom[17].empty()) _atom[17] = _atom[5];
       if (_atom[17] == "?") _atom[17].clear();
       if (_atom[18].empty()) _atom[18] = _atom[2];

       if (_atom[19].empty()) _atom[19] = _atom[15];
       if (_atom[20].empty()) _atom[20] = _atom[16];
       if (_atom[21].empty()) _atom[21] = _atom[17];
       else if (_atom[21] == "?") _atom[21].clear();
       if (_atom[22].empty()) _atom[22] = _atom[18];

       if (_atom[2] != _atom[18]) {
            if (!_atom[2].empty())
                 _atom[18] = _atom[2];
            else _atom[2] = _atom[18];
       }

       if (((_atom[8] == "9999.000") && (_atom[9] == "9999.000") && (_atom[10] == "9999.000")) ||
           ((_atom[8] == "-999.000") && (_atom[9] == "-999.000") && (_atom[10] == "-999.000")) ||
           ((_atom[8] == "999.000")  && (_atom[9] == "999.000")  && (_atom[10] == "999.000"))) _is_9999_flag = true;

       for (unsigned int i = 8; i <= 10; ++i) {
            if (_atom[i].empty()) continue;
            if (_atom[i][0] == '.') _atom[i].insert(0, "0");
            else if (_atom[i].substr(0, 2) == "-.") _atom[i].replace(0, 2, "-0.");
       }

       _orig.x = atof(_atom[8].c_str());
       _orig.y = atof(_atom[9].c_str());
       _orig.z = atof(_atom[10].c_str());
}

void Atom::setAuxiliaryValue(const std::vector<std::string>& vals, const bool convert_flag)
{
       if (vals.empty()) return;

       if (vals[0] == "SIGATM") {
            _sigatm.clear();
            int start = 8;
            if (vals.size() == FORMAT_SIGATM + 1) start = 1;
            for (short int i = 0; i < FORMAT_SIGATM; ++i) {
                 _sigatm_item_position_mapping.insert(std::make_pair(_format_sigatm_pdbx[i], i));
                 _sigatm.push_back(vals[i + start]);
            }
       } else if (vals[0] == "ANISOU") {
            for (short int i = 1; i <= 7; ++i) {
                  _anisotrop_item_position_mapping.insert(std::make_pair(_format_anisotrop_pdbx[i], i));
            }
            _anisotrop_item_position_mapping.insert(std::make_pair(_format_anisotrop_pdbx[14], 13));
            for (short int i = 15; i <= 22; ++i) {
                  _anisotrop_item_position_mapping.insert(std::make_pair(_format_anisotrop_pdbx[i], i));
            }

            _anisou.clear();
            int start = 8;
            if (vals.size() == FORMAT_SIGUIJ + 1) start = 1;
            for (short int i = 0; i < FORMAT_SIGUIJ; ++i) {
                 _anisou_item_position_mapping.insert(std::make_pair(_format_anisotrop_pdbx[i+8], i));
                 if (convert_flag) {
                      double u = atof(vals[i + start].c_str()) / 10000.0;
                      _anisou.push_back(FloatToString(u, 0, 4));
                 } else _anisou.push_back(vals[i + start]);
            }
       } else if (vals[0] == "SIGUIJ") {
            _siguij.clear();
            int start = 8;
            if (vals.size() == FORMAT_SIGUIJ + 1) start = 1;
            for (short int i = 0; i < FORMAT_SIGUIJ; ++i) {
                 _siguij_item_position_mapping.insert(std::make_pair(_format_siguij_pdbx[i], i));
                 if (convert_flag) {
                      double u = atof(vals[i + start].c_str()) / 10000.0;
                      _siguij.push_back(FloatToString(u, 0, 4));
                 } else _siguij.push_back(vals[i + start]);
            }
       }
}

void Atom::setExtraValue(const std::string& item, const std::string& value)
{
       if (item.empty() || value.empty()) return;

       _extra_item_value_mapping.insert(std::make_pair(item, value));
}

void Atom::UpperCase()
{
       if (!_atom[2].empty()) String::UpperCase(_atom[2]);
       if (!_atom[13].empty()) String::UpperCase(_atom[13]);
       if (!_atom[18].empty()) String::UpperCase(_atom[18]);
}

void Atom::writeAtomRecord(std::string& record)
{
       record.clear();

       std::string chn_id = _atom[17]; if (chn_id.empty()) chn_id = " ";
       std::string inscode = _atom[7]; if (inscode.empty()) inscode = " ";
       std::string alt_loc = _atom[3]; if (alt_loc.empty()) alt_loc = " ";

       char buffer[100];
       memset(buffer, 0, 100);
       sprintf(buffer, FORMAT, _atom[1].c_str(), _atom[18].c_str(), alt_loc[0], _atom[16].c_str(), chn_id.c_str(), _atom[15].c_str(),
               inscode[0], _orig.x, _orig.y, _orig.z, atof(_atom[11].c_str()), atof(_atom[12].c_str()));
       record += buffer;
}

double cal_distance(const COORD &coord1, const COORD &coord2)
{
       double dx = coord1.x - coord2.x;
       double dy = coord1.y - coord2.y;
       double dz = coord1.z - coord2.z;

       return (sqrt(dx * dx + dy * dy + dz * dz));
}

double cal_distance(const RCSB::Atom *atom1, const RCSB::Atom *atom2)
{
       return cal_distance(atom1->orig(), atom2->orig());
}

bool is_same_atom(const RCSB::Atom* atom1, const RCSB::Atom* atom2)
{
       if ((atom1->pdb_atmnam() == atom2->pdb_atmnam()) && (atom1->alt_loc() == atom2->alt_loc()) && (atom1->pdb_resnam() == atom2->pdb_resnam()) &&
           (atom1->pdb_chnid() == atom2->pdb_chnid()) && (atom1->pdb_resnum() == atom2->pdb_resnum()) && (atom1->ins_code() == atom2->ins_code()))
            return true;
       return false;
}
